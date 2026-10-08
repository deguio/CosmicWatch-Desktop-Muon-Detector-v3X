#include "oled.h"

#include <string.h>

#include "board.h"
#include "font5x7.h"
#include "hardware/i2c.h"

#define OLED_ADDR 0x3C
#define W 128
#define PAGES 8

static uint8_t fb[1 + W * PAGES];   // fb[0] is the I2C data control byte
static bool present;

static bool cmds(const uint8_t *c, size_t n) {
    uint8_t b[2] = {0x00, 0};
    for (size_t i = 0; i < n; i++) {
        b[1] = c[i];
        if (i2c_write_timeout_us(I2C_PORT, OLED_ADDR, b, 2, false, 2000) != 2) return false;
    }
    return true;
}

bool oled_init(void) {
    static const uint8_t init[] = {
        0xAE,        // display off
        0xD5, 0x80,  // clock divide
        0xA8, 0x3F,  // multiplex 64
        0xD3, 0x00,  // display offset
        0x40,        // start line 0
        0x8D, 0x14,  // charge pump on
        0x20, 0x00,  // horizontal addressing
        0xA1,        // segment remap
        0xC8,        // COM scan descending
        0xDA, 0x12,  // COM pins
        0x81, 0x8F,  // contrast
        0xD9, 0xF1,  // pre-charge
        0xDB, 0x40,  // VCOMH
        0xA4,        // follow RAM
        0xA6,        // normal (not inverted)
        0x2E,        // no scrolling
        0xAF,        // display on
    };
    present = cmds(init, sizeof init);
    if (present) {
        oled_clear();
        present = oled_flush();
    }
    return present;
}

bool oled_present(void) { return present; }

void oled_clear(void) { memset(fb + 1, 0, sizeof fb - 1); }

void oled_text(int line, int x, const char *s) {
    if (line < 0 || line >= PAGES) return;
    uint8_t *row = fb + 1 + line * W;
    for (; *s && x <= W - 5; s++, x += 6) {
        unsigned c = (unsigned char)*s;
        if (c < 0x20 || c > 0x7F) c = '?';
        memcpy(row + x, font5x7[c - 0x20], 5);
    }
}

void oled_invert_line(int line) {
    if (line < 0 || line >= PAGES) return;
    uint8_t *row = fb + 1 + line * W;
    for (int i = 0; i < W; i++) row[i] ^= 0xFF;
}

bool oled_flush(void) {
    if (!present) return false;
    static const uint8_t win[] = {0x21, 0, W - 1, 0x22, 0, PAGES - 1};
    if (!cmds(win, sizeof win)) return false;
    fb[0] = 0x40;
    return i2c_write_timeout_us(I2C_PORT, OLED_ADDR, fb, sizeof fb, false, 100000) == (int)sizeof fb;
}

void oled_power(bool on) {
    uint8_t c = on ? 0xAF : 0xAE;
    if (present) cmds(&c, 1);
}
