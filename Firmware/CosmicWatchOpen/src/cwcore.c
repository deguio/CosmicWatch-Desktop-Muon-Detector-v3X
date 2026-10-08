#include "cwcore.h"

#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "board.h"

// ============================================================== conversions

// Knots (ADC, linear LSB) of the saturation region, medians of the original
// firmware output of seven data sets from different detectors.
static const uint16_t lut_adc[] = {3490, 3515, 3535, 3555, 3585, 3605, 3635, 3655, 3685,
                                   3705, 3725, 3735, 3745, 3755, 3765, 3775, 3785, 3795,
                                   3805, 3818};
static const uint16_t lut_lin[] = {3490, 3555, 3611, 3676, 3765, 3834, 3928, 4003, 4093,
                                   4145, 4200, 4239, 4413, 4832, 5282, 5805, 7206, 9753,
                                   12121, 19430};
#define LUT_N (sizeof lut_adc / sizeof lut_adc[0])

float cw_adc_linearize(uint16_t adc) {
    if (adc <= lut_adc[0]) return (float)adc;
    if (adc >= lut_adc[LUT_N - 1]) return (float)lut_lin[LUT_N - 1];
    size_t i = 1;
    while (adc > lut_adc[i]) i++;
    float f = (float)(adc - lut_adc[i - 1]) / (float)(lut_adc[i] - lut_adc[i - 1]);
    return (float)lut_lin[i - 1] + f * (float)(lut_lin[i] - lut_lin[i - 1]);
}

float cw_sipm_mV(uint16_t adc, float mV_per_lsb, float gain_corr, float baseline_lsb,
                 bool saturation_lut) {
    float lsb = saturation_lut ? cw_adc_linearize(adc) : (float)adc;
    return mV_per_lsb * (baseline_lsb + (lsb - baseline_lsb) * gain_corr);
}

float cw_gain_correction(float hv_V, float temp_C, float tref_C) {
    if (!(temp_C > -40.0f && temp_C < 85.0f)) return -1.0f;
    float vbr_t = SIPM_VBR_21C_V + SIPM_VBR_TC_V * (temp_C - 21.0f);
    float vbr_ref = SIPM_VBR_21C_V + SIPM_VBR_TC_V * (tref_C - 21.0f);
    float ov = hv_V - vbr_t, ov_ref = hv_V - vbr_ref;
    if (ov < 1.0f || ov_ref < 1.0f || ov > 10.0f) return -1.0f;
    return ov_ref / ov;
}

float cw_gain_correction_empirical(float alpha_pct_per_C, float temp_C, float tref_C) {
    if (!(temp_C > -40.0f && temp_C < 85.0f)) return -1.0f;
    return expf(-alpha_pct_per_C / 100.0f * (temp_C - tref_C));
}

// ============================================================== T/P filter

void cw_envfilter_init(cw_envfilter_t *f) { memset(f, 0, sizeof *f); }

static float median(const float *v, int n) {
    float s[CW_ENVF_N] = {0};
    for (int i = 0; i < n; i++) {   // insertion sort, n <= 5
        float x = v[i];
        int j = i;
        while (j > 0 && s[j - 1] > x) {
            s[j] = s[j - 1];
            j--;
        }
        s[j] = x;
    }
    return n % 2 ? s[n / 2] : 0.5f * (s[n / 2 - 1] + s[n / 2]);
}

bool cw_envfilter_push(cw_envfilter_t *f, float T, float P, float dt_s) {
    bool physical = T > -40.0f && T < 85.0f && P > 1000.0f && P < 115000.0f;
    if (!physical) {
        f->rejects++;
        return false;
    }
    // compare with the last accepted reading: the median lags behind on a real ramp
    float last = f->t[(f->idx + CW_ENVF_N - 1) % CW_ENVF_N];
    float max_jump = 2.0f * (dt_s > 1.0f ? dt_s : 1.0f) + 0.5f;
    if (f->valid && fabsf(T - last) > max_jump) {
        if (++f->rejects < CW_ENVF_N) return false;
        f->n = 0;   // persistent change: restart from this reading
        f->idx = 0;
    }
    f->rejects = 0;
    f->t[f->idx] = T;
    f->p[f->idx] = P;
    f->idx = (f->idx + 1) % CW_ENVF_N;
    if (f->n < CW_ENVF_N) f->n++;
    f->T = median(f->t, f->n);
    f->P = median(f->p, f->n);
    f->valid = true;
    return true;
}

// ============================================================== BMP280

void cw_bmp280_parse_calib(cw_bmp280_calib_t *c, const uint8_t r[24]) {
#define U16(i) ((uint16_t)(r[i] | (r[(i) + 1] << 8)))
#define S16(i) ((int16_t)U16(i))
    c->T1 = U16(0); c->T2 = S16(2); c->T3 = S16(4);
    c->P1 = U16(6); c->P2 = S16(8); c->P3 = S16(10); c->P4 = S16(12); c->P5 = S16(14);
    c->P6 = S16(16); c->P7 = S16(18); c->P8 = S16(20); c->P9 = S16(22);
#undef U16
#undef S16
}

void cw_bmp280_compensate(const cw_bmp280_calib_t *c, int32_t adc_T, int32_t adc_P,
                          float *T, float *P) {
    int32_t var1 = ((((adc_T >> 3) - ((int32_t)c->T1 * 2))) * ((int32_t)c->T2)) >> 11;
    int32_t var2 = (((((adc_T >> 4) - ((int32_t)c->T1)) * ((adc_T >> 4) - ((int32_t)c->T1))) >> 12) *
                    ((int32_t)c->T3)) >> 14;
    int32_t t_fine = var1 + var2;
    *T = (float)((t_fine * 5 + 128) >> 8) / 100.0f;

    // Bosch reference algorithm; left shifts of signed values written as products.
    int64_t v1 = ((int64_t)t_fine) - 128000;
    int64_t v2 = v1 * v1 * (int64_t)c->P6;
    v2 = v2 + v1 * (int64_t)c->P5 * 131072;
    v2 = v2 + (int64_t)c->P4 * 34359738368LL;
    v1 = ((v1 * v1 * (int64_t)c->P3) >> 8) + v1 * (int64_t)c->P2 * 4096;
    v1 = ((((int64_t)1) << 47) + v1) * ((int64_t)c->P1) >> 33;
    if (v1 == 0) {
        *P = 0.0f;
        return;
    }
    int64_t p = 1048576 - adc_P;
    p = ((p * 2147483648LL - v2) * 3125) / v1;
    v1 = (((int64_t)c->P9) * (p >> 13) * (p >> 13)) >> 25;
    v2 = (((int64_t)c->P8) * p) >> 19;
    p = ((p + v1 + v2) >> 8) + ((int64_t)c->P7) * 16;
    *P = (float)p / 256.0f;
}

// ============================================================== config.txt

void cw_config_defaults(cw_config_t *c) {
    memset(c, 0, sizeof *c);
    strcpy(c->detector_name, "CosmicWatch");
    c->trigger_mV = 80.0f;
    c->led_brightness = 10;
    c->use_leds = c->use_oled = c->use_serial = true;
    c->use_bmp280 = c->use_mpu6050 = true;
    c->use_buzzer = false;
    c->use_gain_correction = false;
    c->sipm_mV_per_lsb = SIPM_MV_PER_LSB_DEFAULT;
    c->sipm_saturation_lut = true;
    c->buzzer_active = true;
    c->buzzer_all_events = false;
    c->gain_model_empirical = false;
    c->gain_temp_coeff_pct = -0.37f;
    c->gain_ref_temp_C = GAIN_REF_TEMP_C;
    c->hk_period_s = 60;
}

static char *trim(char *s) {
    while (isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = 0;
    return s;
}

static bool parse_bool(const char *v, bool *out) {
    static const char *yes[] = {"true", "1", "yes", "on"};
    static const char *no[] = {"false", "0", "no", "off"};
    for (size_t i = 0; i < 4; i++) {
        if (!strcasecmp(v, yes[i])) { *out = true; return true; }
        if (!strcasecmp(v, no[i])) { *out = false; return true; }
    }
    return false;
}

static bool parse_float(const char *v, float *out) {
    char *end;
    float f = strtof(v, &end);
    while (*end && (isspace((unsigned char)*end) || isalpha((unsigned char)*end) || *end == '%'))
        end++;   // allow units such as "80 mV" or "10%"
    if (end == v || *end || !isfinite(f)) return false;
    *out = f;
    return true;
}

int cw_config_parse(cw_config_t *c, const char *text, cw_log_fn log) {
    static const struct { const char *key; int type; size_t off; } keys[] = {
#define B(k, f) {k, 0, offsetof(cw_config_t, f)}
#define F(k, f) {k, 1, offsetof(cw_config_t, f)}
        B("USE_LEDs", use_leds), B("USE_OLED", use_oled), B("USE_SERIAL", use_serial),
        B("USE_BMP280", use_bmp280), B("USE_MPU6050", use_mpu6050),
        B("USE_BUZZER", use_buzzer), B("USE_REALTIME_GAIN_CORRECTION", use_gain_correction),
        B("SIPM_SATURATION_LUT", sipm_saturation_lut), B("BUZZER_ACTIVE", buzzer_active),
        F("trigger_voltage_mV", trigger_mV), F("SiPM_mV_per_LSB", sipm_mV_per_lsb),
        F("GAIN_TEMP_COEFF", gain_temp_coeff_pct), F("GAIN_REF_TEMP_C", gain_ref_temp_C),
#undef B
#undef F
    };
    int applied = 0;
    char line[160];
    const char *p = text;
    while (*p) {
        size_t len = strcspn(p, "\r\n");
        size_t cp = len < sizeof line - 1 ? len : sizeof line - 1;
        memcpy(line, p, cp);
        line[cp] = 0;
        p += len;
        while (*p == '\r' || *p == '\n') p++;

        char *hash = strchr(line, '#');
        if (hash) *hash = 0;
        char *eq = strchr(line, '=');
        char *s = trim(line);
        if (!*s) continue;
        if (!eq) {
            if (log) log("#    config.txt: ignoring line without '=': %s", s);
            continue;
        }
        *eq = 0;
        char *key = trim(s), *val = trim(eq + 1);

        if (!strcasecmp(key, "detector_name")) {
            if (!*val || strlen(val) > CW_NAME_MAX) {
                if (log) log("#    config.txt: detector_name must have 1-%d characters", CW_NAME_MAX);
                continue;
            }
            strcpy(c->detector_name, val);
            if (cw_sanitize_name(c->detector_name) && log)
                log("#    config.txt: detector_name changed to %s (allowed: A-Z a-z 0-9 -)",
                    c->detector_name);
            applied++;
            continue;
        }
        if (!strcasecmp(key, "detector_version")) continue;   // written by the original firmware
        if (!strcasecmp(key, "GAIN_CORRECTION_MODEL")) {
            if (!strcasecmp(val, "ov")) {
                c->gain_model_empirical = false;
                applied++;
            } else if (!strcasecmp(val, "empirical")) {
                c->gain_model_empirical = true;
                applied++;
            } else if (log) {
                log("#    config.txt: bad GAIN_CORRECTION_MODEL '%s' (ov or empirical)", val);
            }
            continue;
        }
        if (!strcasecmp(key, "HK_PERIOD_S")) {
            float f;
            if (parse_float(val, &f) && (f == 0.0f || (f >= 5.0f && f <= 3600.0f))) {
                c->hk_period_s = (int)lroundf(f);
                applied++;
            } else if (log) {
                log("#    config.txt: bad HK_PERIOD_S '%s' (0 or 5-3600)", val);
            }
            continue;
        }
        if (!strcasecmp(key, "BUZZER_EVENTS")) {
            if (!strcasecmp(val, "coincidence") || !strcasecmp(val, "coincidences")) {
                c->buzzer_all_events = false;
                applied++;
            } else if (!strcasecmp(val, "all")) {
                c->buzzer_all_events = true;
                applied++;
            } else if (log) {
                log("#    config.txt: bad BUZZER_EVENTS '%s' (coincidence or all)", val);
            }
            continue;
        }
        if (!strcasecmp(key, "LED_BRIGHTNESS")) {
            float f;
            if (parse_float(val, &f) && f >= 0 && f <= 100) {
                c->led_brightness = (int)lroundf(f);
                applied++;
            } else if (log) {
                log("#    config.txt: bad LED_BRIGHTNESS '%s' (0-100)", val);
            }
            continue;
        }
        size_t k;
        for (k = 0; k < sizeof keys / sizeof keys[0]; k++)
            if (!strcasecmp(key, keys[k].key)) break;
        if (k == sizeof keys / sizeof keys[0]) {
            if (log) log("#    config.txt: unknown key '%s'", key);
            continue;
        }
        void *dst = (char *)c + keys[k].off;
        bool ok = keys[k].type == 0 ? parse_bool(val, (bool *)dst) : parse_float(val, (float *)dst);
        if (ok) applied++;
        else if (log) log("#    config.txt: bad value '%s' for %s", val, key);
    }
    if (c->trigger_mV < 30.0f || c->trigger_mV > 1200.0f) {
        if (log) log("#    config.txt: trigger_voltage_mV %.1f out of range [30, 1200], using 80",
                     (double)c->trigger_mV);
        c->trigger_mV = 80.0f;
    }
    if (!(c->gain_temp_coeff_pct >= -5.0f && c->gain_temp_coeff_pct <= 5.0f)) {
        if (log) log("#    config.txt: GAIN_TEMP_COEFF out of range [-5, 5] %%/C, using -0.37");
        c->gain_temp_coeff_pct = -0.37f;
    }
    if (!(c->gain_ref_temp_C > -40.0f && c->gain_ref_temp_C < 85.0f)) {
        if (log) log("#    config.txt: GAIN_REF_TEMP_C out of range, using %.1f", (double)GAIN_REF_TEMP_C);
        c->gain_ref_temp_C = GAIN_REF_TEMP_C;
    }
    if (!(c->sipm_mV_per_lsb > 0.0f && c->sipm_mV_per_lsb < 10.0f)) {
        if (log) log("#    config.txt: SiPM_mV_per_LSB out of range, using %.4f",
                     (double)SIPM_MV_PER_LSB_DEFAULT);
        c->sipm_mV_per_lsb = SIPM_MV_PER_LSB_DEFAULT;
    }
    return applied;
}

static const char *tf(bool b) { return b ? "true" : "false"; }

size_t cw_config_write(const cw_config_t *c, char *buf, size_t n) {
    int len = snprintf(buf, n,
        "# CosmicWatch v3X configuration (" FW_NAME " " FW_VERSION ")\n"
        "# Lines are 'key = value'; '#' starts a comment.\n"
        "\n"
        "# Up to %d characters: letters, digits and '-'. Used in file names and in\n"
        "# the last column of the USB data stream.\n"
        "detector_name = %s\n"
        "\n"
        "# Comparator threshold in mV, measured at the TriggerThreshold node.\n"
        "# The amplified signal baseline sits at about 25 mV.\n"
        "trigger_voltage_mV = %.1f\n"
        "\n"
        "LED_BRIGHTNESS = %d\n"
        "USE_LEDs = %s\n"
        "USE_OLED = %s\n"
        "USE_SERIAL = %s\n"
        "USE_BMP280 = %s\n"
        "USE_MPU6050 = %s\n"
        "USE_BUZZER = %s\n"
        "# true: buzzer with internal oscillator (sounds with DC), false: passive buzzer\n"
        "BUZZER_ACTIVE = %s\n"
        "# coincidence: beep only on coincident events; all: beep on every event\n"
        "BUZZER_EVENTS = %s\n"
        "\n"
        "# Normalise SiPM[mV] to the gain at GAIN_REF_TEMP_C (needs the BMP280) and move\n"
        "# the threshold accordingly. Models: ov = SiPM overvoltage from datasheet Vbr(T)\n"
        "# and measured HV; empirical = exp(-GAIN_TEMP_COEFF/100 * (T - GAIN_REF_TEMP_C)),\n"
        "# GAIN_TEMP_COEFF = d ln(MPV)/dT in %%/C from a temperature calibration.\n"
        "USE_REALTIME_GAIN_CORRECTION = %s\n"
        "GAIN_CORRECTION_MODEL = %s\n"
        "GAIN_TEMP_COEFF = %.3f\n"
        "GAIN_REF_TEMP_C = %.1f\n"
        "\n"
        "# Housekeeping comment line (T, P, HV, threshold, baseline, correction) every\n"
        "# HK_PERIOD_S seconds in the data stream; 0 = off.\n"
        "HK_PERIOD_S = %d\n"
        "\n"
        "# SiPM[mV] = SiPM_mV_per_LSB * ADC. With SIPM_SATURATION_LUT the values above\n"
        "# ~3490 ADC follow the saturation curve of the original firmware.\n"
        "SiPM_mV_per_LSB = %.4f\n"
        "SIPM_SATURATION_LUT = %s\n",
        CW_NAME_MAX, c->detector_name, (double)c->trigger_mV, c->led_brightness,
        tf(c->use_leds), tf(c->use_oled), tf(c->use_serial), tf(c->use_bmp280),
        tf(c->use_mpu6050), tf(c->use_buzzer), tf(c->buzzer_active),
        c->buzzer_all_events ? "all" : "coincidence", tf(c->use_gain_correction),
        c->gain_model_empirical ? "empirical" : "ov", (double)c->gain_temp_coeff_pct,
        (double)c->gain_ref_temp_C, c->hk_period_s,
        (double)c->sipm_mV_per_lsb, tf(c->sipm_saturation_lut));
    return len < 0 || (size_t)len >= n ? 0 : (size_t)len;
}

bool cw_sanitize_name(char *name) {
    bool changed = false;
    for (char *p = name; *p; p++) {
        if (!isalnum((unsigned char)*p) && *p != '-') {
            *p = '-';
            changed = true;
        }
    }
    return changed;
}

int cw_parse_run_number(const char *fname, const char *name) {
    size_t nl = strlen(name);
    if (strncasecmp(fname, name, nl) != 0) return -1;
    const char *p = fname + nl;
    if (p[0] != '_' || !p[1] || p[2] != '_') return -1;
    char mode = (char)toupper((unsigned char)p[1]);
    if (mode != 'M' && mode != 'C') return -1;
    p += 3;
    if (!isdigit((unsigned char)*p)) return -1;
    long num = 0;
    while (isdigit((unsigned char)*p)) {
        num = num * 10 + (*p++ - '0');
        if (num > 99999) return -1;
    }
    if (strcasecmp(p, ".txt") != 0) return -1;
    return (int)num;
}

// ============================================================== data lines

static int fmt_us(char *b, size_t n, uint64_t us) {
    return snprintf(b, n, "%" PRIu64 ".%06" PRIu64, us / 1000000u, us % 1000000u);
}

size_t cw_format_event(char *buf, size_t n, const cw_line_t *l, const char *name,
                       const char *eol) {
    size_t o = 0;
    int r;
#define APPEND(...)                                         \
    do {                                                    \
        if (o >= n) return 0;                               \
        r = snprintf(buf + o, n - o, __VA_ARGS__);          \
        if (r < 0 || (size_t)r >= n - o) return 0;          \
        o += (size_t)r;                                     \
    } while (0)
    APPEND("%" PRIu32 "\t", l->event);
    if (o >= n) return 0;
    r = fmt_us(buf + o, n - o, l->t_us);
    if (r < 0 || (size_t)r >= n - o) return 0;
    o += (size_t)r;
    APPEND("\t%u\t%u\t%.1f\t", (unsigned)l->coinc, (unsigned)l->adc, (double)l->sipm_mV);
    r = fmt_us(buf + o, n - o, l->dead_us);
    if (r < 0 || (size_t)r >= n - o) return 0;
    o += (size_t)r;
    if (l->has_env) APPEND("\t%.1f\t%.0f", (double)l->temp_C, (double)l->press_Pa);
    if (l->has_imu)
        APPEND("\t%.3f:%.3f:%.3f\t%.1f:%.1f:%.1f", (double)l->ax, (double)l->ay, (double)l->az,
               (double)l->gx, (double)l->gy, (double)l->gz);
    if (name) APPEND("\t%s", name);
    APPEND("%s", eol);
#undef APPEND
    return o;
}

size_t cw_format_header(char *buf, size_t n, const char *name, const char *fw, bool has_env,
                        bool has_imu, bool usb, const char *eol) {
    char cols[200];
    snprintf(cols, sizeof cols, "# Event  Timestamp[s]  Flag  ADC[12b]  SiPM[mV]  Deadtime[s]%s%s%s",
             has_env ? "  Temp[C]  Press[Pa]" : "",
             has_imu ? "  Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]" : "", usb ? "  Name" : "");
    int w = (int)strlen(cols);
    char bar[200];
    memset(bar, '#', (size_t)w);
    bar[w] = 0;

    char l1[120], l2[120], l3[120];
    snprintf(l1, sizeof l1, "CosmicWatch: The Desktop Muon Detector v3X");
    snprintf(l2, sizeof l2, "Firmware: %s", fw);
    snprintf(l3, sizeof l3, "Detector Name: %s", name);
    int p1 = (w - (int)strlen(l1)) / 2, p2 = (w - (int)strlen(l2)) / 2,
        p3 = (w - (int)strlen(l3)) / 2;
    if (p1 < 1) p1 = 1;
    if (p2 < 1) p2 = 1;
    if (p3 < 1) p3 = 1;
    int len = snprintf(buf, n, "%s%s#%*s%s%s#%*s%s%s#%*s%s%s%s%s%s%s", bar, eol, p1 - 1, "", l1,
                       eol, p2 - 1, "", l2, eol, p3 - 1, "", l3, eol, cols, eol, bar, eol);
    return len < 0 || (size_t)len >= n ? 0 : (size_t)len;
}
