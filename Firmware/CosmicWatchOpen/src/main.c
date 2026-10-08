// CosmicWatch v3X open firmware - core 0: boot diagnostics, data output
// (USB + microSD), OLED, environmental sensors, LEDs, buzzer, USB commands.
// Core 1 runs the acquisition (daq.c).
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "board.h"
#include "coinc.h"
#include "cwcore.h"
#include "daq.h"
#include "ff.h"
#include "diskio.h"
#include "flashcfg.h"
#include "hardware/adc.h"
#include "hardware/gpio.h"
#include "hardware/pwm.h"
#include "hardware/watchdog.h"
#include "logger.h"
#include "oled.h"
#include "pico/bootrom.h"
#include "pico/stdlib.h"
#include "pico/unique_id.h"
#include "sensors.h"
#include "threshold.h"
#include "tusb.h"

#define FW_STRING FW_NAME " " FW_VERSION

static cw_config_t cfg;
static coinc_role_t role;
static bool sd_ok, has_env, has_imu;
static env_t env;
static imu_t imu;
static float hv_V, thr_meas_mV, gain_corr = 1.0f;
// T/P filter and BMP280 health
static cw_envfilter_t envf;
static uint32_t bmp_fail, bmp_reject, env_ok_ms;
// peak-detector baseline: value used for the conversion and per-HK-period statistics
static float base_lsb;
static double base_sum, base_sum2;
static uint32_t base_n;
static threshold_fit_t tfit;
static uint64_t run_start_us;
static uint32_t n_events, n_coinc, usb_dropped, usb_backoff_ms;
static uint64_t last_dead_us;
static bool running, usb_was_connected;

// ======================================================================== output

static char bootlog[8192];
static size_t bootlog_len;

static bool usb_write(const char *s, size_t n, uint32_t max_wait_us) {
    if (!stdio_usb_connected()) return false;
    uint32_t t0 = time_us_32();
    while (n) {
        size_t chunk = n > 256 ? 256 : n;
        while (tud_cdc_write_available() < chunk) {
            if (time_us_32() - t0 > max_wait_us || !stdio_usb_connected()) return false;
            tight_loop_contents();
        }
        stdio_put_string(s, (int)chunk, false, false);
        s += chunk;
        n -= chunk;
    }
    return true;
}

// Comment line: USB (CRLF) and, before the run, the boot log; during the run the
// data file. Format strings do not include the line terminator.
static void logline(const char *fmt, ...) {
    char line[256];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(line, sizeof line - 2, fmt, ap);
    va_end(ap);
    if (n < 0) return;
    if (n > (int)sizeof line - 3) n = (int)sizeof line - 3;
    if (!running && bootlog_len + (size_t)n + 1 < sizeof bootlog) {
        memcpy(bootlog + bootlog_len, line, (size_t)n);
        bootlog_len += (size_t)n;
        bootlog[bootlog_len++] = '\n';
    }
    if (running) {
        line[n] = '\n';
        sd_append(line, (size_t)n + 1);
    }
    line[n] = '\r';
    line[n + 1] = '\n';
    usb_write(line, (size_t)n + 2, 20000);
}

static uint32_t ms_since_boot(void) { return to_ms_since_boot(get_absolute_time()); }

#define LOGT(fmt, ...) logline("# (%lums) " fmt, (unsigned long)ms_since_boot(), ##__VA_ARGS__)

static void usb_dump_bootlog(void) {
    // bootlog uses LF; convert to CRLF on the fly
    const char *p = bootlog, *end = bootlog + bootlog_len;
    while (p < end) {
        const char *nl = memchr(p, '\n', (size_t)(end - p));
        size_t n = nl ? (size_t)(nl - p) : (size_t)(end - p);
        if (!usb_write(p, n, 200000) || !usb_write("\r\n", 2, 200000)) return;
        p += n + 1;
        watchdog_update();
    }
}

static void usb_header(void) {
    char h[1024];
    size_t n = cw_format_header(h, sizeof h, cfg.detector_name, FW_STRING, has_env, has_imu, true,
                                "\r\n");
    usb_write(h, n, 200000);
}

// ======================================================================== LEDs, buzzer

static uint led_slice, buz_slice;
static uint32_t led_off_ms, buz_off_ms;
static bool led_lit, buz_on;

static void leds_init(void) {
    gpio_set_function(PIN_LED_EVENT, GPIO_FUNC_PWM);
    gpio_set_function(PIN_LED_COINC, GPIO_FUNC_PWM);
    led_slice = pwm_gpio_to_slice_num(PIN_LED_EVENT);   // GPIO12 = A, GPIO13 = B
    pwm_config c = pwm_get_default_config();
    pwm_config_set_clkdiv_int(&c, 125);                 // 1 MHz counter at 125 MHz
    pwm_config_set_wrap(&c, 999);                       // 1 kHz
    pwm_init(led_slice, &c, true);
    pwm_set_both_levels(led_slice, 0, 0);

    gpio_init(PIN_BUZZER);   // silent until buzzer_init() knows the buzzer type
    gpio_set_dir(PIN_BUZZER, GPIO_OUT);
    gpio_put(PIN_BUZZER, 0);
}

// Active buzzer (internal oscillator): GPIO6 is simply driven high.
// Passive buzzer: GPIO6 outputs a 2.7 kHz square wave.
static void buzzer_init(void) {
    if (cfg.buzzer_active) return;
    gpio_set_function(PIN_BUZZER, GPIO_FUNC_PWM);
    buz_slice = pwm_gpio_to_slice_num(PIN_BUZZER);
    pwm_config b = pwm_get_default_config();
    pwm_config_set_clkdiv_int(&b, 125);
    pwm_config_set_wrap(&b, 369);                       // 2.7 kHz
    pwm_init(buz_slice, &b, true);
    pwm_set_chan_level(buz_slice, pwm_gpio_to_channel(PIN_BUZZER), 0);
}

static void buzzer_set(bool on) {
    if (cfg.buzzer_active) gpio_put(PIN_BUZZER, on);
    else pwm_set_chan_level(buz_slice, pwm_gpio_to_channel(PIN_BUZZER), on ? 185 : 0);
}

#define BEEP_MS 20   // an active buzzer needs a few ms to start oscillating

static void leds_set(unsigned ev_pct, unsigned co_pct) {
    pwm_set_both_levels(led_slice, (uint16_t)(ev_pct * 10), (uint16_t)(co_pct * 10));
}

static void beep(uint32_t ms) {
    buzzer_set(true);
    buz_off_ms = ms_since_boot() + ms;
    buz_on = true;
}

static void outputs_service(uint32_t now) {
    if (led_lit && (int32_t)(now - led_off_ms) >= 0) {
        leds_set(0, 0);
        led_lit = false;
    }
    if (buz_on && (int32_t)(now - buz_off_ms) >= 0) {
        buzzer_set(false);
        buz_on = false;
    }
}

// ======================================================================== OLED

static void fmt_rate(char *b, size_t n, const char *label, double count, double live_s) {
    if (live_s <= 0) live_s = 1e-9;
    double r = count / live_s, e = sqrt(count) / live_s;
    int d = e >= 1 ? 0 : e >= 0.1 ? 1 : e >= 0.01 ? 2 : 3;
    snprintf(b, n, "%s%.*f+/-%.*fHz", label, d, r, d, e);
}

static void oled_update(uint32_t now_ms) {
    if (!oled_present() || !cfg.use_oled) return;
    char l[32];
    uint64_t now = time_us_64();
    double run_s = (double)(now - run_start_us) * 1e-6;
    double live_s = run_s - (double)last_dead_us * 1e-6;
    int x = (int)((now_ms / 60000u) % 3u);   // shift by one pixel per minute (burn-in)

    oled_clear();
    unsigned long s = (unsigned long)run_s;
    snprintf(l, sizeof l, "%-10.10s %c %02lu:%02lu:%02lu", cfg.detector_name, role.found ? 'C' : 'M',
             s / 3600, (s / 60) % 60, s % 60);
    oled_text(0, x, l);
    oled_text(1, x, sd_active() ? sd_filename() : "No microSD logging");
    snprintf(l, sizeof l, "Events: %lu", (unsigned long)n_events);
    oled_text(2, x, l);
    fmt_rate(l, sizeof l, "", n_events, live_s);
    oled_text(3, x, l);
    int line = 4;
    if (role.found) {
        snprintf(l, sizeof l, "Coinc.: %lu", (unsigned long)n_coinc);
        oled_text(line++, x, l);
        fmt_rate(l, sizeof l, "", n_coinc, live_s);
        oled_text(line++, x, l);
    }
    snprintf(l, sizeof l, "Thr %.1fmV HV %.1fV", (double)thr_meas_mV, (double)hv_V);
    oled_text(line++, x, l);
    if (has_env && env.ok) {
        snprintf(l, sizeof l, "%.1f\x7f" "C %.1fhPa", (double)env.temp_C, (double)env.press_Pa / 100.0);
        oled_text(line++, x, l);
    } else if (line < 8) {
        snprintf(l, sizeof l, "Dead %.3f%%", run_s > 0 ? 100.0 * (run_s - live_s) / run_s : 0.0);
        oled_text(line++, x, l);
    }
    oled_flush();
}

static void oled_message(const char *a, const char *b, const char *c) {
    if (!oled_present() || !cfg.use_oled) return;
    oled_clear();
    oled_text(1, 0, a);
    if (b) oled_text(3, 0, b);
    if (c) oled_text(5, 0, c);
    oled_flush();
}

static const char *sensor_state(bool enabled, bool found) {
    return !enabled ? "OFF" : found ? "OK" : "MISSING";
}

// Boot screen shown only when an enabled sensor does not answer.
static void oled_sensor_check(void) {
    bool bmp_bad = cfg.use_bmp280 && !has_env, mpu_bad = cfg.use_mpu6050 && !has_imu;
    if (!(bmp_bad || mpu_bad) || !oled_present() || !cfg.use_oled) return;
    char l[32];
    oled_clear();
    oled_text(0, 0, "    SENSOR CHECK");
    oled_invert_line(0);
    snprintf(l, sizeof l, "BMP280 T/P:  %s", sensor_state(cfg.use_bmp280, has_env));
    oled_text(2, 0, l);
    snprintf(l, sizeof l, "MPU6050 IMU: %s", sensor_state(cfg.use_mpu6050, has_imu));
    oled_text(3, 0, l);
    oled_text(5, 0, "Missing sensor data");
    oled_text(6, 0, "is not written to");
    oled_text(7, 0, "the data file.");
    oled_flush();
    sleep_ms(4000);
}

// ======================================================================== events

static void process_event(const cw_event_t *ev) {
    n_events++;
    if (ev->coinc) n_coinc++;
    last_dead_us = ev->dead_us;

    cw_line_t l = {0};
    l.event = n_events;
    l.t_us = ev->t_us - run_start_us;
    l.dead_us = ev->dead_us;
    l.coinc = ev->coinc;
    l.adc = ev->adc;
    l.sipm_mV = cw_sipm_mV(ev->adc, cfg.sipm_mV_per_lsb, gain_corr, base_lsb, cfg.sipm_saturation_lut);
    if (has_env) {
        l.has_env = true;
        l.temp_C = env.temp_C;
        l.press_Pa = env.press_Pa;
    }
    if (has_imu) {
        l.has_imu = true;
        l.ax = imu.ax; l.ay = imu.ay; l.az = imu.az;
        l.gx = imu.gx; l.gy = imu.gy; l.gz = imu.gz;
    }

    char line[256];
    size_t n = cw_format_event(line, sizeof line, &l, NULL, "\n");
    sd_append(line, n);
    if (cfg.use_serial) {
        n = cw_format_event(line, sizeof line, &l, cfg.detector_name, "\r\n");
        // A host that keeps the port open without reading must not stall the
        // acquisition: after a failed write, skip USB output for one second.
        uint32_t t = ms_since_boot();
        if (stdio_usb_connected()) {
            if ((int32_t)(t - usb_backoff_ms) < 0) {
                usb_dropped++;
            } else if (!usb_write(line, n, 3000)) {
                usb_dropped++;
                usb_backoff_ms = t + 1000;
            }
        }
    }

    uint32_t now = ms_since_boot();
    if (cfg.use_leds) {
        led_off_ms = now + 20;   // switched on by core 1
        led_lit = true;
    }
    if (cfg.use_buzzer && (ev->coinc || cfg.buzzer_all_events)) beep(BEEP_MS);
}

// ======================================================================== housekeeping

static float measure_hv(void) {
    float raw = daq_adc_average(ADC_CH_HV, 64);
    return raw < 0 ? -1.0f : raw * ADC_VREF_MV / ADC_FULL_SCALE * HV_DIVIDER / 1000.0f;
}

static void apply_threshold(void) {
    float target = cfg.trigger_mV;
    if (cfg.use_gain_correction && gain_corr > 0)
        target = AMP_BIAS_MV + (cfg.trigger_mV - AMP_BIAS_MV) / gain_corr;
    threshold_set_mV(target);
}

// Raw BMP280 reading -> median filter -> env (used in the data columns and OLED).
static void env_update(uint32_t now, float dt_s) {
    env_t raw;
    if (!bmp280_read(&raw)) {
        bmp_fail++;
        return;
    }
    if (!cw_envfilter_push(&envf, raw.temp_C, raw.press_Pa, dt_s)) {
        bmp_reject++;
        return;
    }
    env.ok = true;
    env.temp_C = envf.T;
    env.press_Pa = envf.P;
    env_ok_ms = now;
}

static float gain_model(void) {
    if (cfg.gain_model_empirical)
        return cw_gain_correction_empirical(cfg.gain_temp_coeff_pct, env.temp_C, cfg.gain_ref_temp_C);
    return cw_gain_correction(hv_V, env.temp_C, cfg.gain_ref_temp_C);
}

// "# HK ..." comment line: ignored by the parsers, read by the calibration analysis.
// base/base_sd/base_n summarise the baseline samples since the previous HK line.
static void hk_line(const char *src) {
    double t = (double)(time_us_64() - run_start_us) * 1e-6;
    char T[16] = "NA", P[16] = "NA";
    if (has_env && env.ok) {
        snprintf(T, sizeof T, "%.2f", (double)env.temp_C);
        snprintf(P, sizeof P, "%.0f", (double)env.press_Pa);
    }
    double bm = base_lsb, bs = 0;
    if (base_n) {
        bm = base_sum / base_n;
        bs = sqrt(fmax(base_sum2 / base_n - bm * bm, 0.0));
        if (base_n >= 5) base_lsb = (float)bm;
    }
    logline("# HK t=%.3f T=%s P=%s HV=%.3f Thr_set=%.2f Thr_meas=%.2f base=%.2f base_sd=%.2f "
            "base_n=%lu corr=%.5f events=%lu coinc=%lu dead=%.6f lost=%lu/%lu/%lu bmp_fail=%lu "
            "bmp_rej=%lu src=%s",
            t, T, P, (double)hv_V, (double)threshold_get_set_mV(), (double)thr_meas_mV, bm, bs,
            (unsigned long)base_n, (double)gain_corr, (unsigned long)n_events,
            (unsigned long)n_coinc, (double)daq_dead_us() * 1e-6, (unsigned long)daq_dropped(),
            (unsigned long)usb_dropped, (unsigned long)sd_dropped_bytes(), (unsigned long)bmp_fail,
            (unsigned long)bmp_reject, src);
    base_sum = base_sum2 = 0;
    base_n = 0;
}

static void housekeeping(uint32_t now) {
    static uint32_t t_imu, t_env, t_hk, t_base, t_oled, t_gain, t_hb, t_hkline;
    static uint32_t hb_last;
    static float corr_logged = 1.0f;

    if (has_imu && now - t_imu >= 100) {
        t_imu = now;
        mpu6050_read(&imu);
    }
    if (has_env && now - t_env >= 1000) {
        float dt = (float)(now - t_env) / 1000.0f;
        t_env = now;
        env_update(now, dt);
    }
    if (now - t_hk >= 1000) {
        t_hk = now;
        float v = measure_hv();
        if (v > 0) hv_V = v;
        thr_meas_mV = threshold_measure_mV();
    }
    // ~63 us of dead time per sample: one every 5 s (12 per default HK period)
    if (now - t_base >= 5000) {
        t_base = now;
        float b = daq_baseline_sample();
        if (b >= 0) {
            base_sum += b;
            base_sum2 += (double)b * b;
            base_n++;
        }
    }
    if (now - t_gain >= 10000) {
        t_gain = now;
        bool changed = false;
        // only with a recent, filtered temperature
        if (cfg.use_gain_correction && env.ok && now - env_ok_ms < 30000) {
            float c = gain_model();
            if (c > 0.8f && c < 1.25f && fabsf(c - gain_corr) > 1e-4f) {
                gain_corr = c;
                apply_threshold();
                changed = true;   // RC filter needs time: trim at the next round
            }
        }
        // keep the threshold on target against drifts of the 3.3 V PWM level
        if (!changed) {
            float target = threshold_get_set_mV();
            float m = threshold_measure_mV();
            if (m > 0 && fabsf(m - target) > 0.3f) threshold_trim(target, m);
        }
        if (fabsf(gain_corr - corr_logged) > 0.002f) {
            corr_logged = gain_corr;
            hk_line("corr");
        }
    }
    if (cfg.hk_period_s && now - t_hkline >= (uint32_t)cfg.hk_period_s * 1000u) {
        t_hkline = now;
        hk_line("periodic");
    }
    if (now - t_oled >= 1000) {
        t_oled = now;
        oled_update(now);
    }
    if (now - t_hb >= 500) {
        t_hb = now;
        uint32_t hb = daq_heartbeat();
        if (hb != hb_last) watchdog_update();   // core 1 alive
        hb_last = hb;
    }
}

// ======================================================================== USB commands

static void print_status(void) {
    double run_s = (double)(time_us_64() - run_start_us) * 1e-6;
    double dead_s = (double)last_dead_us * 1e-6;
    logline("# %s | detector %s | mode %s | file %s", FW_STRING, cfg.detector_name,
            role.found ? (role.leader ? "C (leader)" : "C (follower)") : "M",
            sd_active() ? sd_filename() : "-");
    logline("# run %.1f s, live %.1f s, events %lu, coincident %lu", run_s, run_s - dead_s,
            (unsigned long)n_events, (unsigned long)n_coinc);
    logline("# threshold set %.2f mV, measured %.2f mV, HV %.2f V, gain corr %.4f, baseline %.2f LSB",
            (double)threshold_get_set_mV(), (double)thr_meas_mV, (double)hv_V, (double)gain_corr,
            (double)base_lsb);
    if (has_env)
        logline("# BMP280: T %.2f C, P %.0f Pa, failed reads %lu, rejected %lu, last good %lu s ago",
                (double)env.temp_C, (double)env.press_Pa, (unsigned long)bmp_fail,
                (unsigned long)bmp_reject, (unsigned long)((ms_since_boot() - env_ok_ms) / 1000));
    logline("# lost: queue %lu, USB %lu, SD %lu bytes; stuck triggers %lu",
            (unsigned long)daq_dropped(), (unsigned long)usb_dropped,
            (unsigned long)sd_dropped_bytes(), (unsigned long)daq_stuck_triggers());
}

static void handle_command(char *cmd) {
    char *arg = strchr(cmd, ' ');
    if (arg) {
        *arg++ = 0;
        while (*arg == ' ') arg++;
    }
    if (!strcmp(cmd, "help")) {
        logline("# commands: status | name <NAME> | threshold <mV> | reboot | bootsel");
    } else if (!strcmp(cmd, "status")) {
        print_status();
    } else if (!strcmp(cmd, "name") && arg && *arg) {
        char nm[CW_NAME_MAX + 1];
        strncpy(nm, arg, CW_NAME_MAX);
        nm[CW_NAME_MAX] = 0;
        cw_sanitize_name(nm);
        bool ok = flashcfg_save(nm, cfg.trigger_mV);
        logline("# name '%s' %s; it is used from the next boot%s", nm,
                ok ? "saved to flash" : "NOT saved", sd_ok ? " (config.txt on the microSD wins)" : "");
    } else if (!strcmp(cmd, "threshold") && arg) {
        float v = strtof(arg, NULL);
        if (v < 30.0f || v > 1200.0f) {
            logline("# threshold must be within 30-1200 mV");
            return;
        }
        cfg.trigger_mV = v;
        apply_threshold();
        bool ok = flashcfg_save(cfg.detector_name, v);
        logline("# threshold set to %.1f mV%s", (double)v, ok ? " and saved to flash" : "");
    } else if (!strcmp(cmd, "reboot")) {
        logline("# rebooting");
        sleep_ms(50);
        watchdog_reboot(0, 0, 0);
    } else if (!strcmp(cmd, "bootsel")) {
        logline("# rebooting into the USB bootloader");
        sleep_ms(50);
        rom_reset_usb_boot(0, 0);
    } else if (*cmd) {
        logline("# unknown command '%s' (try: help)", cmd);
    }
}

static void usb_service(void) {
    static char cmd[64];
    static size_t len;
    bool conn = stdio_usb_connected();
    if (conn && !usb_was_connected && running) {
        sleep_ms(50);   // some terminals drop data sent right after opening the port
        usb_dump_bootlog();
        usb_header();
    }
    usb_was_connected = conn;
    int c;
    while ((c = getchar_timeout_us(0)) >= 0) {
        if (c == '\r' || c == '\n') {
            cmd[len] = 0;
            if (len) handle_command(cmd);
            len = 0;
        } else if (len < sizeof cmd - 1) {
            cmd[len++] = (char)c;
        }
    }
}

// ======================================================================== boot

static void cfg_log(const char *fmt, ...) {
    char b[200];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    logline("%s", b);
}

static void default_name(char *dst) {
    char id[2 * PICO_UNIQUE_BOARD_ID_SIZE_BYTES + 1];
    pico_get_unique_board_id_string(id, sizeof id);
    snprintf(dst, CW_NAME_MAX + 1, "CW-%s", id + strlen(id) - 4);
}

static void measure_baseline(float *mean, float *sd) {
    double s = 0, s2 = 0;
    const int n = 256;
    adc_select_input(ADC_CH_SIGNAL);
    for (int i = 0; i < n; i++) {
        gpio_put(PIN_PEAK_RESET, 1);
        busy_wait_us(5);
        gpio_put(PIN_PEAK_RESET, 0);
        busy_wait_us(50);
        double v = adc_read();
        s += v;
        s2 += v * v;
    }
    *mean = (float)(s / n);
    *sd = (float)sqrt(fmax(s2 / n - (s / n) * (s / n), 0.0));
}

static void boot(void) {
    bool wd = watchdog_caused_reboot();
    stdio_init_all();

    gpio_init(PIN_PEAK_RESET);
    gpio_set_dir(PIN_PEAK_RESET, GPIO_OUT);
    gpio_put(PIN_PEAK_RESET, 1);   // keep C5 discharged during boot
    gpio_init(PIN_SD_DETECT);
    gpio_set_dir(PIN_SD_DETECT, GPIO_IN);
    gpio_pull_up(PIN_SD_DETECT);
    leds_init();

    // Coincidence discovery first, so that both detectors reach it together.
    role = coinc_handshake();
    if (role.found) leds_set(100, 100);   // both LEDs on: partner found

    cw_config_defaults(&cfg);
    float flash_thr;
    bool have_flash = flashcfg_load(cfg.detector_name, sizeof cfg.detector_name, &flash_thr);
    if (have_flash) cfg.trigger_mV = flash_thr;
    else default_name(cfg.detector_name);

    logline("# Welcome to the CosmicWatch v3X open firmware (Version: %s)", FW_STRING);
    if (wd) LOGT("WARNING: the previous run ended with a watchdog reset");
    LOGT("Coincidence handshake: %s", role.found ? (role.leader ? "partner found (leader, out=GPIO0)"
                                                                : "partner found (follower, out=GPIO1)")
                                                 : "no partner, master mode");

    LOGT("Initializing I2C bus (GPIO14/15)");
    i2c_bus_init();
    i2c_scan(cfg_log);
    LOGT("Initializing OLED screen ...");
    if (oled_init()) {
        LOGT("OLED ready (SSD1306 at 0x3C)");
        oled_message("    CosmicWatch v3X", "   open firmware", "      " FW_VERSION);
    } else {
        LOGT("WARNING: OLED not found. Check GPIO14/15 and that GND is an outer OLED pin.");
    }

    // microSD: config.txt and data file
    LOGT("microSD card-detect pin: %s", gpio_get(PIN_SD_DETECT) ? "high" : "low");
    LOGT("Initializing microSD card communication");
    if (disk_initialize(0) == 0) {
        LOGT("Mounting microSD card");
        sd_ok = sd_mount(cfg_log);
    } else {
        LOGT("No microSD card (or card not responding): data only over USB");
    }
    if (sd_ok) {
        LOGT("Reading config.txt");
        sd_load_config(&cfg, cfg_log);
    }
    logline("#    Detector Name set to: %s", cfg.detector_name);
    logline("#    Trigger voltage [mV] set to: %.2f mV", (double)cfg.trigger_mV);
    logline("#    LED_BRIGHTNESS set to: %d%%", cfg.led_brightness);
    logline("#    Use LEDs set to: %s", cfg.use_leds ? "true" : "false");
    logline("#    Use OLED set to: %s", cfg.use_oled ? "true" : "false");
    logline("#    Use SERIAL set to: %s", cfg.use_serial ? "true" : "false");
    logline("#    Use BMP280 set to: %s", cfg.use_bmp280 ? "true" : "false");
    logline("#    Use MPU6050 set to: %s", cfg.use_mpu6050 ? "true" : "false");
    logline("#    Use Buzzer set to: %s (%s buzzer, beep on %s)", cfg.use_buzzer ? "true" : "false",
            cfg.buzzer_active ? "active" : "passive",
            cfg.buzzer_all_events ? "all events" : "coincidences only");
    logline("#    Use REALTIME GAIN CORRECTION set to: %s", cfg.use_gain_correction ? "true" : "false");
    if (cfg.gain_model_empirical)
        logline("#    Gain correction model: empirical, %.3f %%/C, Tref %.1f C",
                (double)cfg.gain_temp_coeff_pct, (double)cfg.gain_ref_temp_C);
    else
        logline("#    Gain correction model: ov (Vbr 24.45 V at 21 C, 21.5 mV/C), Tref %.1f C",
                (double)cfg.gain_ref_temp_C);
    logline("#    Housekeeping line every: %d s%s", cfg.hk_period_s, cfg.hk_period_s ? "" : " (off)");
    logline("#    SiPM conversion: %.4f mV/LSB, saturation LUT %s", (double)cfg.sipm_mV_per_lsb,
            cfg.sipm_saturation_lut ? "on" : "off");
    if (sd_ok) {
        // remember name and threshold for runs without microSD
        char fn[CW_NAME_MAX + 1];
        float ft;
        if (!flashcfg_load(fn, sizeof fn, &ft) || strcmp(fn, cfg.detector_name) || ft != cfg.trigger_mV)
            flashcfg_save(cfg.detector_name, cfg.trigger_mV);
        if (sd_open_run(cfg.detector_name, role.found ? 'C' : 'M', cfg_log)) {
            LOGT("Saving data to microSD card file name: %s", sd_filename());
            sd_space(cfg_log);
        }
    }
    if (!cfg.use_oled) oled_power(false);
    buzzer_init();

    // environmental sensors
    if (cfg.use_bmp280) {
        LOGT("Initializing BMP280 temperature/pressure sensor");
        has_env = bmp280_init(cfg_log);
        cw_envfilter_init(&envf);
        if (has_env) {
            sleep_ms(60);
            env_update(ms_since_boot(), 1.0f);
            if (env.ok)
                LOGT("BMP280: temperature %.2f C, pressure %.2f Pa", (double)env.temp_C, (double)env.press_Pa);
        } else {
            LOGT("WARNING: BMP280 not detected");
        }
    } else {
        LOGT("WARNING: BMP280 disabled in config file.");
    }
    if (cfg.use_mpu6050) {
        LOGT("Initializing MPU-6050 accelerometer/gyroscope");
        has_imu = mpu6050_init(cfg_log);
        if (has_imu) mpu6050_read(&imu);
        else LOGT("WARNING: MPU-6050 not detected");
    } else {
        LOGT("WARNING: MPU6050 disabled in config file.");
    }
    oled_sensor_check();

    // analog front-end
    LOGT("Initializing ADC (VREF = 2.5 V)");
    adc_init();
    adc_gpio_init(PIN_ADC_SIGNAL);
    adc_gpio_init(PIN_ADC_THRESH);
    adc_gpio_init(PIN_ADC_HV);

    LOGT("Initializing trigger threshold PWM + calibration");
    threshold_init();
    tfit = threshold_calibrate();
    if (tfit.ok)
        LOGT("Threshold fit: V_meas_mV = %.6f * duty + %.6f   (R^2=%.4f, N=%d)", (double)tfit.slope_mV,
             (double)tfit.offset_mV, (double)tfit.r2, tfit.npoints);
    else
        LOGT("WARNING: threshold fit FAILED (N=%d), using nominal 1650 mV/duty", tfit.npoints);
    apply_threshold();
    sleep_ms(5);
    thr_meas_mV = threshold_measure_mV();
    threshold_trim(threshold_get_set_mV(), thr_meas_mV);
    sleep_ms(5);
    thr_meas_mV = threshold_measure_mV();
    LOGT("Trigger Threshold (Measured, Expected) : (%.2f mV, %.2f mV)", (double)thr_meas_mV,
         (double)threshold_get_set_mV());
    if (fabsf(thr_meas_mV - threshold_get_set_mV()) > 0.2f * threshold_get_set_mV()) {
        logline("#  -> ERROR: trigger threshold not within 20%% of the set value.");
        logline("#     Check R19, R21, C11 and the PWM line (GPIO20).");
        oled_message("Threshold Error!", NULL, NULL);
    }

    hv_V = measure_hv();
    LOGT("Checking: HV to SiPM voltage: %.2f V", (double)hv_V);
    if (fabsf(hv_V - HV_NOMINAL_V) > 0.2f * HV_NOMINAL_V) {
        logline("#  -> ERROR: HV to SiPM not within 20%% of %.1f V.", (double)HV_NOMINAL_V);
        logline("#     Check the DC-DC booster: U1 (MAX5026), L1, D3, R1, R2, FB1.");
        oled_message("HV Error!", "Should be 30V", "Check DC-DC booster");
    }

    gpio_put(PIN_PEAK_RESET, 0);
    sleep_us(100);
    float sig = daq_adc_average(ADC_CH_SIGNAL, 64) * ADC_VREF_MV / ADC_FULL_SCALE;
    LOGT("Checking: Signal line voltage: %.2f mV", (double)sig);
    float bl, bsd;
    measure_baseline(&bl, &bsd);
    LOGT("Measuring ADC0 baseline: %.2f LSB +/- %.2f LSB", (double)bl, (double)bsd);
    base_lsb = bl;
    if (bl > 600.0f || bsd > 30.0f) {
        logline("#  -> WARNING: unusual baseline. Check ~25 mV across R22 and 2.5 V across C25,");
        logline("#     the peak detector (U2, D4, C5) and the reset MOSFET Q2.");
    }
    if (gpio_get(PIN_TRIGGER)) {
        sleep_ms(2);
        if (gpio_get(PIN_TRIGGER))
            logline("#  -> WARNING: the trigger line is high: threshold below the signal baseline?");
    }
}

int main(void) {
    boot();

    daq_config_t dc = {
        .coinc_mode = role.found,
        .pin_out = role.pin_out,
        .pin_in = role.pin_in,
        .reset_us = 5,
        .led_event_level_pct = (uint8_t)(cfg.use_leds ? cfg.led_brightness : 0),
    };
    if (role.found) leds_set(0, 0);
    daq_init(&dc);
    LOGT("Coincidence window: %.2f us", (double)daq_window_us());
    LOGT("Launching CosmicWatch %s detector!", FW_STRING);

    char h[1024];
    size_t n = cw_format_header(h, sizeof h, cfg.detector_name, FW_STRING, has_env, has_imu, false, "\n");
    if (sd_active()) {
        sd_append(bootlog, bootlog_len);
        sd_append(h, n);
    }
    usb_header();
    usb_was_connected = stdio_usb_connected();

    oled_clear();
    oled_flush();
    running = true;
    run_start_us = time_us_64();
    daq_start(run_start_us);
    hk_line("start");
    watchdog_enable(8000, true);

    cw_event_t ev;
    for (;;) {
        int k = 0;
        while (k++ < 64 && daq_pop(&ev)) process_event(&ev);
        uint32_t now = ms_since_boot();
        outputs_service(now);
        usb_service();
        housekeeping(now);
        sd_service(now);
    }
}
