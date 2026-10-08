// Hardware-independent logic (conversions, config file, data format).
// Compiled both into the firmware and into the host unit tests.
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define CW_NAME_MAX 20   // characters, excluding the terminator

// ---------------------------------------------------------------- conversions
// Peak-detector ADC value expressed in "linear LSB": identity up to ~3490, then
// the saturation curve observed in the output of the original firmware, ending
// at 19430 (~1370 mV) for fully saturated pulses.
float cw_adc_linearize(uint16_t adc);
// The gain correction scales only the signal above the baseline (both in LSB):
//   SiPM[mV] = k * (baseline + (lin(adc) - baseline) * gain_corr)
// With gain_corr = 1 this is k * lin(adc), as in the original firmware.
float cw_sipm_mV(uint16_t adc, float mV_per_lsb, float gain_corr, float baseline_lsb,
                 bool saturation_lut);

// Gain normalisation to tref_C, SiPM overvoltage model:
//   (HV - Vbr(Tref)) / (HV - Vbr(T)),  Vbr(T) = 24.45 V + 21.5 mV/C * (T - 21 C).
// Returns a negative value if the inputs are not physical.
float cw_gain_correction(float hv_V, float temp_C, float tref_C);
// Empirical model, alpha = d ln(MPV)/dT in %/C (negative):  exp(-alpha/100 * (T - Tref)).
float cw_gain_correction_empirical(float alpha_pct_per_C, float temp_C, float tref_C);

// ---------------------------------------------------------------- T/P filter
// Median of the last 5 accepted BMP280 readings. A reading is rejected if it is
// not physical or if T jumps by more than 2 C/s from the last accepted reading; after
// 5 consecutive rejections the filter restarts from the new reading.
#define CW_ENVF_N 5
typedef struct {
    float t[CW_ENVF_N], p[CW_ENVF_N];
    int   n, idx, rejects;
    bool  valid;
    float T, P;   // filtered values (valid only if valid == true)
} cw_envfilter_t;

void cw_envfilter_init(cw_envfilter_t *f);
bool cw_envfilter_push(cw_envfilter_t *f, float T, float P, float dt_s);

// ---------------------------------------------------------------- BMP280
typedef struct {
    uint16_t T1; int16_t T2, T3;
    uint16_t P1; int16_t P2, P3, P4, P5, P6, P7, P8, P9;
} cw_bmp280_calib_t;

void cw_bmp280_parse_calib(cw_bmp280_calib_t *c, const uint8_t raw[24]);
// Bosch reference integer compensation. T in C, P in Pa.
void cw_bmp280_compensate(const cw_bmp280_calib_t *c, int32_t adc_T, int32_t adc_P,
                          float *T, float *P);

// ---------------------------------------------------------------- config.txt
typedef struct {
    char  detector_name[CW_NAME_MAX + 1];
    float trigger_mV;
    int   led_brightness;          // percent
    bool  use_leds, use_oled, use_serial, use_bmp280, use_mpu6050, use_buzzer;
    bool  use_gain_correction;
    float sipm_mV_per_lsb;
    bool  sipm_saturation_lut;
    bool  buzzer_active;           // active buzzer: drive GPIO6 high instead of a 2.7 kHz tone
    bool  buzzer_all_events;       // beep on every event, not only on coincidences
    bool  gain_model_empirical;    // GAIN_CORRECTION_MODEL = empirical (else ov)
    float gain_temp_coeff_pct;     // GAIN_TEMP_COEFF, %/C, used by the empirical model
    float gain_ref_temp_C;         // GAIN_REF_TEMP_C
    int   hk_period_s;             // HK_PERIOD_S, 0 = no housekeeping lines
} cw_config_t;

typedef void (*cw_log_fn)(const char *fmt, ...);

void cw_config_defaults(cw_config_t *c);
// Parses "key = value" lines; '#' starts a comment. Unknown keys and bad values
// are reported through log (may be NULL). Returns the number of keys applied.
int cw_config_parse(cw_config_t *c, const char *text, cw_log_fn log);
// Writes a complete, commented config.txt. Returns its length.
size_t cw_config_write(const cw_config_t *c, char *buf, size_t n);

// Replace characters that would break file names or tab-separated columns.
// Returns true if the name was changed.
bool cw_sanitize_name(char *name);

// "<name>_<M|C>_<NNN>.txt" -> NNN, or -1 if fname does not belong to name.
int cw_parse_run_number(const char *fname, const char *name);

// ---------------------------------------------------------------- data lines
typedef struct {
    uint32_t event;
    uint64_t t_us;           // since the start of the run
    uint64_t dead_us;
    uint8_t  coinc;
    uint16_t adc;
    float    sipm_mV;
    bool     has_env;
    float    temp_C, press_Pa;
    bool     has_imu;
    float    ax, ay, az, gx, gy, gz;
} cw_line_t;

// One tab-separated event line terminated by eol. A non-NULL name is appended
// as last column (USB stream). Returns the length, 0 if it does not fit.
size_t cw_format_event(char *buf, size_t n, const cw_line_t *l, const char *name,
                       const char *eol);
// Column header block, compatible with the original firmware layout.
size_t cw_format_header(char *buf, size_t n, const char *name, const char *fw,
                        bool has_env, bool has_imu, bool usb, const char *eol);
