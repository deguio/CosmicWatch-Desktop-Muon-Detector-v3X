// SSD1306 128x64 OLED on I2C, 8 lines x 21 characters.
#pragma once
#include <stdbool.h>
#include <stdint.h>

bool oled_init(void);
bool oled_present(void);
void oled_clear(void);
// Text line (0-7) at a horizontal pixel offset; '\x7f' prints a degree sign.
void oled_text(int line, int x, const char *s);
void oled_invert_line(int line);
bool oled_flush(void);
// Send one of the 8 pages (rows of 8 pixels): ~3 ms on the I2C bus instead of
// ~25 ms for the whole frame, so the caller can serve events in between.
bool oled_flush_page(int page);
void oled_power(bool on);
