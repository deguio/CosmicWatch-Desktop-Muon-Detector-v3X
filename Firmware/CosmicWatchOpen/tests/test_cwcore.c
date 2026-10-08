// Host unit tests for src/cwcore.c:  ./run_tests.sh
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "cwcore.h"

static int fails, checks;
#define CHECK(c)                                                         \
    do {                                                                 \
        checks++;                                                        \
        if (!(c)) {                                                      \
            fails++;                                                     \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c);          \
        }                                                                \
    } while (0)
#define NEAR(a, b, tol) CHECK(fabs((double)(a) - (double)(b)) <= (tol))

static int nlog;
static void count_log(const char *fmt, ...) {
    (void)fmt;
    nlog++;
}

static void test_bmp280(void) {
    // Worked example of the Bosch BMP280 datasheet (section 8.2 / 3.12).
    cw_bmp280_calib_t c = {27504, 26435, -1000, 36477, -10685, 3024, 2855, 140, -7, 15500, -14600, 6000};
    float T, P;
    cw_bmp280_compensate(&c, 519888, 415148, &T, &P);
    NEAR(T, 25.08, 0.005);
    // Exact result of the 64-bit integer algorithm: p = 25767233 / 256 Pa (the
    // datasheet table quotes the rounded 100653.27; the double version gives .258).
    NEAR(P, 100653.254, 0.008);   // float spacing at 1e5 is 0.0078

    uint8_t raw[24] = {0x70, 0x6B, 0x43, 0x67, 0x18, 0xFC};   // T1, T2, T3 little endian
    cw_bmp280_calib_t d;
    cw_bmp280_parse_calib(&d, raw);
    CHECK(d.T1 == 27504 && d.T2 == 26435 && d.T3 == -1000);
}

static void test_conversion(void) {
    NEAR(cw_adc_linearize(100), 100, 0);
    NEAR(cw_adc_linearize(3490), 3490, 0);
    NEAR(cw_adc_linearize(4095), 19430, 0);
    NEAR(cw_adc_linearize(3818), 19430, 0);
    float prev = 0;
    int monotonic = 1;
    for (int a = 0; a < 4096; a++) {
        float v = cw_adc_linearize((uint16_t)a);
        if (v < prev) monotonic = 0;
        prev = v;
    }
    CHECK(monotonic);
    // Values written by the original firmware (Data/TOP_C_002.txt, k = 0.0706)
    NEAR(cw_sipm_mV(130, 0.0706f, 1, 55, true), 9.2, 0.1);
    NEAR(cw_sipm_mV(1106, 0.0706f, 1, 55, true), 78.1, 0.5);
    NEAR(cw_sipm_mV(4090, 0.0706f, 1, 55, true), 1371.8, 5);
    NEAR(cw_sipm_mV(4090, 0.0706f, 1, 55, false), 288.8, 0.1);

    NEAR(cw_gain_correction(30.2f, 25.0f, 25.0f), 1.0, 1e-6);
    CHECK(cw_gain_correction(30.2f, 35.0f, 25.0f) > 1.0f);   // warmer -> lower gain -> scale up
    CHECK(cw_gain_correction(24.0f, 25.0f, 25.0f) < 0);      // below breakdown: invalid
    // OV model slope at the reference: 21.5 mV/C / OV, OV = 30.2 - 24.536 V
    NEAR((cw_gain_correction(30.2f, 26.0f, 25.0f) - 1) * 100, 0.0215 / (30.2 - 24.536) * 100, 0.002);

    // empirical model: alpha = -0.4 %/C -> +4.0 % at 10 C above Tref
    NEAR(cw_gain_correction_empirical(-0.4f, 35.0f, 25.0f), exp(0.04), 1e-6);
    NEAR(cw_gain_correction_empirical(-0.4f, 25.0f, 25.0f), 1.0, 0);
    CHECK(cw_gain_correction_empirical(-0.4f, 99.0f, 25.0f) < 0);

    // the correction scales only the part above the baseline
    NEAR(cw_sipm_mV(355, 0.1f, 1.10f, 55, false), 0.1 * (55 + 300 * 1.10), 1e-4);
    NEAR(cw_sipm_mV(355, 0.1f, 1.0f, 55, false), 35.5, 1e-4);
}

static void test_envfilter(void) {
    cw_envfilter_t f;
    cw_envfilter_init(&f);
    CHECK(!f.valid);
    // sentinel written by the original firmware on a failed BMP280 read
    CHECK(!cw_envfilter_push(&f, -5.0f, 0.0f, 1));
    CHECK(cw_envfilter_push(&f, 23.2f, 99626, 1));
    CHECK(f.valid && f.T == 23.2f && f.P == 99626);
    CHECK(!cw_envfilter_push(&f, -5.0f, 0.0f, 1));
    CHECK(cw_envfilter_push(&f, 23.3f, 99620, 1));
    CHECK(cw_envfilter_push(&f, 23.4f, 99610, 1));
    NEAR(f.T, 23.3, 1e-5);   // median of 23.2, 23.3, 23.4
    // isolated spike (physical range but 10 C/s): rejected
    CHECK(!cw_envfilter_push(&f, 33.0f, 99600, 1));
    NEAR(f.T, 23.3, 1e-5);
    // a fast but real ramp (1.5 C/s) is followed
    float T = 23.4f;
    for (int i = 0; i < 10; i++) CHECK(cw_envfilter_push(&f, T += 1.5f, 99600, 1));
    NEAR(f.T, T - 3.0f, 1e-4);   // median lags by two samples
    // a persistent step (e.g. sensor replaced) is accepted after 5 rejections
    for (int i = 0; i < 4; i++) CHECK(!cw_envfilter_push(&f, 5.0f, 99600, 1));
    CHECK(cw_envfilter_push(&f, 5.0f, 99600, 1));
    NEAR(f.T, 5.0, 0);
    // balloon flight: low pressure is physical
    CHECK(cw_envfilter_push(&f, 5.2f, 7000, 1));
}

static void test_config(void) {
    cw_config_t c;
    cw_config_defaults(&c);
    const char *txt =
        "# comment\r\n"
        "detector_name = TOP\r\n"
        "detector_version = v3X.28\n"
        "trigger_voltage_mV = 92.5 mV\n"
        "LED_BRIGHTNESS = 35%\n"
        "USE_LEDs = false\n"
        "use_oled = OFF   # trailing comment\n"
        "USE_BUZZER = true\n"
        "USE_REALTIME_GAIN_CORRECTION = yes\n"
        "SiPM_mV_per_LSB = 0.0700\n"
        "bogus = 3\n"
        "USE_BMP280 = maybe\n";
    nlog = 0;
    int n = cw_config_parse(&c, txt, count_log);
    CHECK(n == 8);
    CHECK(nlog == 2);   // unknown key + bad boolean
    CHECK(!strcmp(c.detector_name, "TOP"));
    NEAR(c.trigger_mV, 92.5, 1e-4);
    CHECK(c.led_brightness == 35);
    CHECK(!c.use_leds && !c.use_oled && c.use_buzzer && c.use_gain_correction);
    CHECK(c.use_bmp280);   // unchanged after the bad value
    NEAR(c.sipm_mV_per_lsb, 0.07, 1e-6);

    // buzzer: defaults, valid and invalid values
    CHECK(c.buzzer_active && !c.buzzer_all_events);
    nlog = 0;
    CHECK(cw_config_parse(&c, "BUZZER_ACTIVE = false\nBUZZER_EVENTS = All\n", count_log) == 2);
    CHECK(!c.buzzer_active && c.buzzer_all_events && nlog == 0);
    CHECK(cw_config_parse(&c, "BUZZER_EVENTS = singles\n", count_log) == 0);
    CHECK(nlog == 1 && c.buzzer_all_events);
    CHECK(cw_config_parse(&c, "BUZZER_EVENTS = coincidence\n", count_log) == 1);
    CHECK(!c.buzzer_all_events);

    // gain model and housekeeping
    CHECK(!c.gain_model_empirical && c.hk_period_s == 60);
    nlog = 0;
    CHECK(cw_config_parse(&c, "GAIN_CORRECTION_MODEL = Empirical\nGAIN_TEMP_COEFF = -0.452\n"
                              "GAIN_REF_TEMP_C = 20\nHK_PERIOD_S = 30\n", count_log) == 4);
    CHECK(nlog == 0 && c.gain_model_empirical && c.hk_period_s == 30);
    NEAR(c.gain_temp_coeff_pct, -0.452, 1e-6);
    NEAR(c.gain_ref_temp_C, 20, 0);
    nlog = 0;
    cw_config_parse(&c, "GAIN_CORRECTION_MODEL = magic\nHK_PERIOD_S = 2\nGAIN_TEMP_COEFF = 12\n", count_log);
    CHECK(nlog == 3 && c.gain_model_empirical && c.hk_period_s == 30);
    NEAR(c.gain_temp_coeff_pct, -0.37, 1e-6);   // out of range -> default
    CHECK(cw_config_parse(&c, "HK_PERIOD_S = 0\n", count_log) == 1 && c.hk_period_s == 0);

    cw_config_parse(&c, "detector_name = My Lab_1\ntrigger_voltage_mV = 5\n", count_log);
    CHECK(!strcmp(c.detector_name, "My-Lab-1"));
    NEAR(c.trigger_mV, 80, 0);   // out of range -> default

    nlog = 0;
    cw_config_parse(&c, "detector_name = ThisNameIsWayTooLongForTheDetector\n", count_log);
    CHECK(nlog == 1 && !strcmp(c.detector_name, "My-Lab-1"));

    // round trip
    char buf[4096];
    cw_config_t a, b;
    cw_config_defaults(&a);
    strcpy(a.detector_name, "BOTTOM");
    a.trigger_mV = 71.5f;
    a.use_mpu6050 = false;
    a.buzzer_active = false;
    a.buzzer_all_events = true;
    a.gain_model_empirical = true;
    a.gain_temp_coeff_pct = -0.452f;
    a.hk_period_s = 120;
    size_t len = cw_config_write(&a, buf, sizeof buf);
    CHECK(len > 0);
    cw_config_defaults(&b);
    nlog = 0;
    cw_config_parse(&b, buf, count_log);
    CHECK(nlog == 0);
    CHECK(!strcmp(b.detector_name, "BOTTOM") && b.trigger_mV == 71.5f && !b.use_mpu6050);
    CHECK(!b.buzzer_active && b.buzzer_all_events);
    CHECK(b.gain_model_empirical && b.gain_temp_coeff_pct == -0.452f && b.hk_period_s == 120);
    CHECK(cw_config_write(&a, buf, 100) == 0);   // too small
}

static void test_names(void) {
    CHECK(cw_parse_run_number("TOP_C_002.txt", "TOP") == 2);
    CHECK(cw_parse_run_number("TOP_M_117.TXT", "TOP") == 117);
    CHECK(cw_parse_run_number("top_m_1000.txt", "TOP") == 1000);
    CHECK(cw_parse_run_number("TOPX_C_002.txt", "TOP") == -1);
    CHECK(cw_parse_run_number("TOP_X_002.txt", "TOP") == -1);
    CHECK(cw_parse_run_number("TOP_C_.txt", "TOP") == -1);
    CHECK(cw_parse_run_number("TOP_C_002.csv", "TOP") == -1);
    CHECK(cw_parse_run_number("config.txt", "TOP") == -1);
    char n[] = "a b\tc_d";
    CHECK(cw_sanitize_name(n) && !strcmp(n, "a-b-c-d"));
}

static void test_format(void) {
    cw_line_t l = {0};
    l.event = 28;
    l.t_us = 9598471;
    l.dead_us = 10511;
    l.coinc = 1;
    l.adc = 527;
    l.sipm_mV = 37.2f;
    char b[256];
    size_t n = cw_format_event(b, sizeof b, &l, NULL, "\n");
    CHECK(n == strlen(b));
    CHECK(!strcmp(b, "28\t9.598471\t1\t527\t37.2\t0.010511\n"));

    l.has_env = true;
    l.temp_C = 24.81f;
    l.press_Pa = 101383.4f;
    l.has_imu = true;
    l.ax = 0.177f; l.ay = 0.034f; l.az = 1.092f;
    l.gx = 63.7f; l.gy = 57.4f; l.gz = 4.8f;
    cw_format_event(b, sizeof b, &l, "AxLab", "\r\n");
    CHECK(!strcmp(b, "28\t9.598471\t1\t527\t37.2\t0.010511\t24.8\t101383\t0.177:0.034:1.092\t63.7:57.4:4.8\tAxLab\r\n"));
    CHECK(cw_format_event(b, 20, &l, "AxLab", "\n") == 0);

    l.t_us = 86400ull * 365 * 1000000ull;   // one year of run time
    cw_format_event(b, sizeof b, &l, NULL, "\n");
    CHECK(!strncmp(b, "28\t31536000.000000\t", 19));

    char h[1024];
    n = cw_format_header(h, sizeof h, "TOP", "CosmicWatch-open 1.0.0", true, true, false, "\n");
    CHECK(n > 0);
    CHECK(strstr(h, "# Event  Timestamp[s]  Flag  ADC[12b]  SiPM[mV]  Deadtime[s]  Temp[C]  Press[Pa]  "
                    "Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]\n") != NULL);
    CHECK(strstr(h, "Detector Name: TOP") != NULL);
    // every header line is a comment
    for (char *p = h; *p; p = strchr(p, '\n') + 1) CHECK(*p == '#');
    printf("%s", h);
}

int main(void) {
    test_bmp280();
    test_conversion();
    test_envfilter();
    test_config();
    test_names();
    test_format();
    printf("%d checks, %d failures\n", checks, fails);
    return fails != 0;
}
