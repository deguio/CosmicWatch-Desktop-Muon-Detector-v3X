// Trigger threshold DAC: 1 MHz PWM on GPIO20, filtered by R21/R19/C11.
#pragma once
#include <stdbool.h>

typedef struct {
    float slope_mV;     // V_thr[mV] = slope_mV * duty + offset_mV
    float offset_mV;
    float r2;
    int   npoints;
    bool  ok;
} threshold_fit_t;

void threshold_init(void);
// Sweep the duty cycle, read back ADC1 and fit a straight line.
threshold_fit_t threshold_calibrate(void);
// Set the threshold in mV at the comparator input (fit must be valid).
void threshold_set_mV(float mV);
float threshold_get_set_mV(void);
// Read back the threshold through ADC1, mV.
float threshold_measure_mV(void);
// Small closed-loop correction after a measurement.
void threshold_trim(float target_mV, float measured_mV);
