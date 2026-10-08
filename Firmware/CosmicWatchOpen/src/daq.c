// Core 1 acquisition loop.
//
// Core 1 runs with interrupts disabled and everything it executes lives in RAM
// (__not_in_flash_func plus inline register access only), so flash/XIP activity
// and interrupts on core 0 never add latency to the readout.
#include "daq.h"

#include <string.h>

#include "board.h"
#include "daq.pio.h"
#include "hardware/adc.h"
#include "hardware/clocks.h"
#include "hardware/gpio.h"
#include "hardware/pio.h"
#include "hardware/pwm.h"
#include "hardware/structs/timer.h"
#include "hardware/sync.h"
#include "pico/multicore.h"
#include "pico/stdlib.h"

#define RING_SIZE 2048u  // power of two
#define TRIGGER_LOW_TIMEOUT_US 1000u

static cw_event_t ring[RING_SIZE];

static struct {
    volatile uint32_t head;          // written by core 1
    volatile uint32_t tail;          // written by core 0
    volatile uint32_t dropped;
    volatile uint32_t heartbeat;
    volatile uint32_t stuck;
    volatile uint64_t dead_us;
    volatile bool     run;
    volatile bool     pause_req;
    volatile bool     paused;
    // housekeeping ADC request, served between events
    volatile uint32_t hk_channel;
    volatile uint32_t hk_todo;
    volatile uint32_t hk_sum;
    volatile bool     hk_busy;
    // baseline sample request, served only while the PIO waits for a trigger
    volatile bool     bl_req;
    volatile uint32_t bl_val;
} sh;

static daq_config_t cfg;
static PIO pio = pio0;
static uint sm;
static uint armed_pc;   // program counter of the "wait 1 pin" instruction
static bool launched;
static uint32_t led_slice, led_level;

static inline uint64_t now_us(void) {
    uint32_t hi = timer_hw->timerawh;
    for (;;) {
        uint32_t lo = timer_hw->timerawl;
        uint32_t hi2 = timer_hw->timerawh;
        if (hi == hi2) return ((uint64_t)hi << 32) | lo;
        hi = hi2;
    }
}

static inline void spin_us(uint32_t us) {
    uint32_t t0 = timer_hw->timerawl;
    while (timer_hw->timerawl - t0 < us) {
    }
}

static inline uint32_t adc_convert(void) {
    adc_hw->cs |= ADC_CS_START_ONCE_BITS;
    while (!(adc_hw->cs & ADC_CS_READY_BITS)) {
    }
    return adc_hw->result;
}

static inline void adc_mux(uint32_t ch) {
    hw_write_masked(&adc_hw->cs, ch << ADC_CS_AINSEL_LSB, ADC_CS_AINSEL_BITS);
}

static void __not_in_flash_func(read_event)(void) {
    uint64_t t0 = now_us();
    (void)pio->rxf[sm];                          // trigger marker
    while (pio->fstat & (1u << (PIO_FSTAT_RXEMPTY_LSB + sm))) {
    }
    uint32_t coinc = pio->rxf[sm] & 1u;          // after the ~3 us window

    uint32_t adc = adc_convert();                // peak held on C5

    sio_hw->gpio_set = 1u << PIN_PEAK_RESET;     // discharge C5 through Q2
    spin_us(cfg.reset_us);
    uint32_t tw = timer_hw->timerawl;
    while (sio_hw->gpio_in & (1u << PIN_TRIGGER)) {
        if (timer_hw->timerawl - tw > TRIGGER_LOW_TIMEOUT_US) {
            sh.stuck++;
            break;
        }
    }
    spin_us(1);
    sio_hw->gpio_clr = 1u << PIN_PEAK_RESET;
    spin_us(1);
    uint64_t t1 = now_us();
    pio->txf[sm] = 1;                            // re-arm

    sh.dead_us += t1 - t0;
    if (led_level) {
        // channel A = GPIO12 (event LED); the coincidence LED is channel B
        uint32_t cc = led_level | (coinc && cfg.coinc_mode ? led_level << 16 : 0);
        pwm_hw->slice[led_slice].cc = cc;
    }

    uint32_t head = sh.head;
    if (head - sh.tail >= RING_SIZE) {
        sh.dropped++;
        return;
    }
    cw_event_t *ev = &ring[head & (RING_SIZE - 1)];
    ev->t_us = t0;
    ev->dead_us = sh.dead_us;
    ev->adc = (uint16_t)adc;
    ev->coinc = (uint8_t)(cfg.coinc_mode ? coinc : 0);
    __dmb();
    sh.head = head + 1;
}

static void __not_in_flash_func(housekeeping_sample)(void) {
    adc_mux(sh.hk_channel);
    uint32_t v = adc_convert();
    adc_mux(ADC_CH_SIGNAL);
    sh.hk_sum += v;
    __dmb();
    if (--sh.hk_todo == 0) sh.hk_busy = false;
}

// Peak-detector baseline with the trigger disabled: stop the state machine
// while it waits for a trigger, discharge C5, let it settle 50 us (as in the
// boot measurement) and convert once. About 60 us, counted as dead time.
static void __not_in_flash_func(baseline_sample)(void) {
    const uint32_t en = 1u << (PIO_CTRL_SM_ENABLE_LSB + sm);
    hw_clear_bits(&pio->ctrl, en);
    bool idle = (pio->fstat & (1u << (PIO_FSTAT_RXEMPTY_LSB + sm))) && pio->sm[sm].addr == armed_pc;
    if (!idle) {   // an event is in progress: retry later
        hw_set_bits(&pio->ctrl, en);
        return;
    }
    uint64_t t0 = now_us();
    sio_hw->gpio_set = 1u << PIN_PEAK_RESET;
    spin_us(cfg.reset_us);
    sio_hw->gpio_clr = 1u << PIN_PEAK_RESET;
    spin_us(50);
    uint32_t v = adc_convert();
    // a pulse during the measurement would leave a partial peak on C5
    uint32_t tw = timer_hw->timerawl;
    while ((sio_hw->gpio_in & (1u << PIN_TRIGGER)) && timer_hw->timerawl - tw < TRIGGER_LOW_TIMEOUT_US) {
    }
    sio_hw->gpio_set = 1u << PIN_PEAK_RESET;
    spin_us(cfg.reset_us);
    sio_hw->gpio_clr = 1u << PIN_PEAK_RESET;
    spin_us(1);
    hw_set_bits(&pio->ctrl, en);
    sh.dead_us += now_us() - t0;
    sh.bl_val = v;
    __dmb();
    sh.bl_req = false;
}

static void __not_in_flash_func(core1_main)(void) {
    (void)save_and_disable_interrupts();
    for (;;) {
        sh.heartbeat++;
        if (sh.pause_req) {
            uint64_t tp = now_us();
            sh.paused = true;
            while (sh.pause_req) {
                sh.heartbeat++;
            }
            sh.dead_us += now_us() - tp;
            sh.paused = false;
            continue;
        }
        if (!sh.run) continue;
        if (!(pio->fstat & (1u << (PIO_FSTAT_RXEMPTY_LSB + sm)))) {
            read_event();
        } else if (sh.hk_busy) {
            housekeeping_sample();
        } else if (sh.bl_req) {
            baseline_sample();
        }
    }
}

void daq_init(const daq_config_t *c) {
    cfg = *c;

    gpio_init(PIN_TRIGGER);
    gpio_set_dir(PIN_TRIGGER, GPIO_IN);
    gpio_disable_pulls(PIN_TRIGGER);

    gpio_init(PIN_PEAK_RESET);
    gpio_set_dir(PIN_PEAK_RESET, GPIO_OUT);
    gpio_put(PIN_PEAK_RESET, 0);

    led_slice = pwm_gpio_to_slice_num(PIN_LED_EVENT);
    led_level = (uint32_t)c->led_event_level_pct * 10u;   // LED PWM wrap is 999

    adc_select_input(ADC_CH_SIGNAL);

    uint offset = pio_add_program(pio, &cw_daq_program);
    armed_pc = offset + cw_daq_offset_armed;
    sm = (uint)pio_claim_unused_sm(pio, true);
    cw_daq_program_init(pio, sm, offset, PIN_TRIGGER, c->coinc_mode ? c->pin_out : -1,
                        (uint)c->pin_in);
    pio_sm_set_enabled(pio, sm, true);

    multicore_launch_core1(core1_main);
    launched = true;
}

void daq_start(uint64_t run_start_us) {
    (void)run_start_us;
    sh.dead_us = 0;
    pio_sm_clear_fifos(pio, sm);
    pio_sm_put(pio, sm, 1);   // first arm
    __dmb();
    sh.run = true;
}

bool daq_pop(cw_event_t *ev) {
    uint32_t tail = sh.tail;
    if (tail == sh.head) return false;
    __dmb();
    *ev = ring[tail & (RING_SIZE - 1)];
    __dmb();
    sh.tail = tail + 1;
    return true;
}

uint32_t daq_dropped(void) { return sh.dropped; }
uint32_t daq_heartbeat(void) { return sh.heartbeat; }
uint32_t daq_stuck_triggers(void) { return sh.stuck; }
uint64_t daq_dead_us(void) { return sh.dead_us; }

float daq_adc_average(unsigned channel, unsigned n) {
    if (n == 0) n = 1;
    if (!launched) {
        adc_select_input(channel);
        uint32_t sum = 0;
        for (unsigned i = 0; i < n; i++) sum += adc_read();
        adc_select_input(ADC_CH_SIGNAL);
        return (float)sum / (float)n;
    }
    sh.hk_channel = channel;
    sh.hk_sum = 0;
    sh.hk_todo = n;
    __dmb();
    sh.hk_busy = true;
    uint64_t t0 = time_us_64();
    while (sh.hk_busy) {
        if (time_us_64() - t0 > 100000) return -1.0f;   // core 1 not serving
        tight_loop_contents();
    }
    __dmb();
    return (float)sh.hk_sum / (float)n;
}

float daq_baseline_sample(void) {
    if (!launched) return -1.0f;
    sh.bl_req = true;
    uint64_t t0 = time_us_64();
    while (sh.bl_req) {
        if (time_us_64() - t0 > 100000) {   // trigger busy for 100 ms: give up this time
            sh.bl_req = false;
            return -1.0f;
        }
        tight_loop_contents();
    }
    __dmb();
    return (float)sh.bl_val;
}

void daq_pause(void) {
    if (!launched) return;
    sh.pause_req = true;
    while (!sh.paused) tight_loop_contents();
}

void daq_resume(void) {
    if (!launched) return;
    sh.pause_req = false;
    while (sh.paused) tight_loop_contents();
}

float daq_window_us(void) {
    return (float)CW_DAQ_WINDOW_CYCLES * 1e6f / (float)clock_get_hz(clk_sys);
}
