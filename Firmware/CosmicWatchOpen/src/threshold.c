// The original design sets the threshold with a 1 MHz PWM. At 125 MHz that is
// only 125 duty steps (13.2 mV each at the comparator). Here a DMA channel,
// paced by the PWM wrap, cycles through 256 compare values that alternate
// between two adjacent steps (first-order sigma-delta). The RC filter
// (tau = 0.5 ms) averages them, giving ~0.05 mV resolution with the same
// negligible 1 MHz ripple.
#include "threshold.h"

#include <math.h>

#include "board.h"
#include "daq.h"
#include "hardware/clocks.h"
#include "hardware/dma.h"
#include "hardware/gpio.h"
#include "hardware/irq.h"
#include "hardware/pwm.h"
#include "pico/stdlib.h"

#define TABLE_LEN 256u

static uint32_t table[TABLE_LEN] __attribute__((aligned(TABLE_LEN * sizeof(uint32_t))));
static uint slice, chan;
static int dma_ch = -1;
static uint32_t wrap;            // PWM counts per period
static threshold_fit_t fit;
static float set_mV;
static float set_duty;

static void dma_restart_irq(void) {
    if (dma_channel_get_irq0_status((uint)dma_ch)) {
        dma_channel_acknowledge_irq0((uint)dma_ch);
        dma_channel_set_trans_count((uint)dma_ch, 0xFFFFFFFFu, true);
    }
}

static void fill_table(float duty) {
    if (duty < 0.0f) duty = 0.0f;
    if (duty > 1.0f) duty = 1.0f;
    float level = duty * (float)wrap;
    uint32_t base = (uint32_t)level;
    uint32_t frac = (uint32_t)lroundf((level - (float)base) * TABLE_LEN);
    uint32_t acc = 0;
    for (uint32_t i = 0; i < TABLE_LEN; i++) {
        acc += frac;
        uint32_t lv = base;
        if (acc >= TABLE_LEN) {
            acc -= TABLE_LEN;
            lv++;
        }
        if (lv > wrap) lv = wrap;
        table[i] = chan == PWM_CHAN_A ? lv : lv << 16;
    }
    set_duty = duty;
}

void threshold_init(void) {
    gpio_set_function(PIN_THRESH_PWM, GPIO_FUNC_PWM);
    slice = pwm_gpio_to_slice_num(PIN_THRESH_PWM);
    chan = pwm_gpio_to_channel(PIN_THRESH_PWM);
    wrap = clock_get_hz(clk_sys) / THRESH_PWM_HZ;   // 125 at 125 MHz

    pwm_config pc = pwm_get_default_config();
    pwm_config_set_clkdiv_int(&pc, 1);
    pwm_config_set_wrap(&pc, (uint16_t)(wrap - 1));
    pwm_init(slice, &pc, false);
    fill_table(0.0f);
    pwm_set_chan_level(slice, chan, 0);
    pwm_set_enabled(slice, true);

    dma_ch = dma_claim_unused_channel(true);
    dma_channel_config dc = dma_channel_get_default_config((uint)dma_ch);
    channel_config_set_transfer_data_size(&dc, DMA_SIZE_32);
    channel_config_set_read_increment(&dc, true);
    channel_config_set_write_increment(&dc, false);
    channel_config_set_ring(&dc, false, 10);          // 256 x 4 bytes
    channel_config_set_dreq(&dc, pwm_get_dreq(slice));
    dma_channel_configure((uint)dma_ch, &dc, &pwm_hw->slice[slice].cc, table, 0xFFFFFFFFu,
                          false);
    // ~71 minutes per transfer at 1 MHz: re-trigger from the completion interrupt.
    dma_channel_set_irq0_enabled((uint)dma_ch, true);
    irq_add_shared_handler(DMA_IRQ_0, dma_restart_irq, PICO_SHARED_IRQ_HANDLER_DEFAULT_ORDER_PRIORITY);
    irq_set_enabled(DMA_IRQ_0, true);
    dma_channel_start((uint)dma_ch);
}

static float measure_settled(float duty) {
    fill_table(duty);
    sleep_ms(4);                                   // 8 RC time constants
    return threshold_measure_mV();
}

float threshold_measure_mV(void) {
    float raw = daq_adc_average(ADC_CH_THRESH, 64);
    return raw * ADC_VREF_MV / ADC_FULL_SCALE;
}

threshold_fit_t threshold_calibrate(void) {
    // Integer duty steps from 1 to 1.2 V; above that the comparator is useless.
    double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    int n = 0;
    for (uint32_t lv = 1; lv <= wrap; lv++) {
        float duty = (float)lv / (float)wrap;
        float v = measure_settled(duty);
        if (v > 1200.0f) break;
        if (v >= 2490.0f) continue;                   // ADC saturated
        sx += duty; sy += v; sxx += (double)duty * duty; sxy += (double)duty * v;
        syy += (double)v * v;
        n++;
    }
    fit.npoints = n;
    fit.ok = false;
    if (n >= 5) {
        double den = n * sxx - sx * sx;
        double a = (n * sxy - sx * sy) / den;
        double b = (sy - a * sx) / n;
        double ss_tot = syy - sy * sy / n;
        double ss_res = syy - b * sy - a * sxy;
        fit.slope_mV = (float)a;
        fit.offset_mV = (float)b;
        fit.r2 = ss_tot > 0 ? (float)(1.0 - ss_res / ss_tot) : 0.0f;
        // Expected slope: 3.3 V / 2 = 1650 mV per unit duty.
        fit.ok = a > 1000.0 && a < 2500.0 && fit.r2 > 0.99f;
    }
    fill_table(0.0f);
    return fit;
}

void threshold_set_mV(float mV) {
    set_mV = mV;
    float a = fit.ok ? fit.slope_mV : 1650.0f;
    float b = fit.ok ? fit.offset_mV : 0.0f;
    fill_table((mV - b) / a);
}

float threshold_get_set_mV(void) { return set_mV; }

void threshold_trim(float target_mV, float measured_mV) {
    float a = fit.ok ? fit.slope_mV : 1650.0f;
    float d = (target_mV - measured_mV) / a;
    if (fabsf(d) < 0.01f) fill_table(set_duty + d);
}
