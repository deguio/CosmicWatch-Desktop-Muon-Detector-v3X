// Prints USB-format event lines (with every sensor combination) for check_format.py.
#include <stdio.h>

#include "cwcore.h"

int main(void) {
    char b[256];
    for (int mask = 0; mask < 4; mask++) {
        for (uint32_t i = 1; i <= 3; i++) {
            cw_line_t l = {0};
            l.event = i;
            l.t_us = 1234567ull * i;
            l.dead_us = 11ull * i;
            l.coinc = (uint8_t)(i == 2);
            l.adc = (uint16_t)(100 * i);
            l.sipm_mV = cw_sipm_mV(l.adc, 0.0706f, 1.0f, 55.0f, true);
            l.has_env = mask & 1;
            l.temp_C = 21.37f;
            l.press_Pa = 100512.6f;
            l.has_imu = mask & 2;
            l.ax = -0.005f; l.ay = 0.0f; l.az = -0.998f;
            l.gx = -0.04f; l.gy = 0.1f; l.gz = 0.3f;
            cw_format_event(b, sizeof b, &l, "TOP", "\r\n");
            fputs(b, stdout);
        }
    }
    return 0;
}
