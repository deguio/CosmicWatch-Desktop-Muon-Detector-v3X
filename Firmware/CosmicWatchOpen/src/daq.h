// Event acquisition running alone on core 1.
#pragma once
#include <stdbool.h>
#include <stdint.h>

typedef struct {
    uint64_t t_us;      // trigger time, us since boot
    uint64_t dead_us;   // cumulative dead time since the start of the run
    uint16_t adc;       // peak detector, 12 bit
    uint8_t  coinc;     // 1 if the partner detector triggered within the window
    uint8_t  _pad[5];
} cw_event_t;

typedef struct {
    bool     coinc_mode;     // a partner was found at boot
    int      pin_out;        // our coincidence line (-1 in single mode)
    int      pin_in;         // partner coincidence line
    uint32_t reset_us;       // peak-detector discharge time
    uint8_t  led_event_level_pct;
} daq_config_t;

// Prepare PIO, ADC and GPIOs and launch core 1. Acquisition starts with daq_start().
void daq_init(const daq_config_t *cfg);
void daq_start(uint64_t run_start_us);

// Called by core 0. Returns false when no event is pending.
bool daq_pop(cw_event_t *ev);
uint32_t daq_dropped(void);       // events lost because the queue was full
uint32_t daq_heartbeat(void);     // incremented continuously by core 1
uint32_t daq_stuck_triggers(void);
uint64_t daq_dead_us(void);

// Average n ADC conversions of a channel. While the DAQ runs the conversions are
// interleaved by core 1 between events; before daq_init they are done directly.
float daq_adc_average(unsigned channel, unsigned n);

// One ADC0 sample of the discharged peak detector (baseline), taken while the
// trigger is armed and idle; ~63 us of dead time. Returns -1 on timeout.
float daq_baseline_sample(void);

// Park core 1 in RAM (e.g. while writing flash). Paused time counts as dead time.
void daq_pause(void);
void daq_resume(void);

// Window length of the coincidence sampling, microseconds.
float daq_window_us(void);
