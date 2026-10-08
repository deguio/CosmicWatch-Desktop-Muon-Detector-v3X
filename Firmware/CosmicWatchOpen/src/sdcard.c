// microSD card in SPI mode on SPI1, glued to FatFs through the diskio API.
#include <stdbool.h>
#include <string.h>

#include "board.h"
#include "ff.h"
#include "diskio.h"
#include "hardware/gpio.h"
#include "hardware/spi.h"
#include "pico/stdlib.h"

#define CMD0   0
#define CMD8   8
#define CMD9   9
#define CMD12 12
#define CMD16 16
#define CMD17 17
#define CMD18 18
#define CMD24 24
#define CMD25 25
#define CMD55 55
#define CMD58 58
#define ACMD41 (0x80 | 41)

#define SPI_SLOW_HZ   400000u
#define SPI_FAST_HZ 12500000u

static volatile DSTATUS stat = STA_NOINIT;
static bool block_addr;   // SDHC/SDXC: sector addressing

static inline void cs_low(void) { gpio_put(PIN_SD_CS, 0); }
static inline void cs_high(void) { gpio_put(PIN_SD_CS, 1); }

static uint8_t xchg(uint8_t b) {
    uint8_t r;
    spi_write_read_blocking(SD_SPI, &b, &r, 1);
    return r;
}

static bool wait_ready(uint32_t timeout_ms) {
    absolute_time_t t = make_timeout_time_ms(timeout_ms);
    do {
        if (xchg(0xFF) == 0xFF) return true;
    } while (!time_reached(t));
    return false;
}

static void deselect(void) {
    cs_high();
    xchg(0xFF);   // release MISO
}

static bool select_card(void) {
    cs_low();
    xchg(0xFF);
    if (wait_ready(500)) return true;
    deselect();
    return false;
}

static uint8_t send_cmd(uint8_t cmd, uint32_t arg) {
    if (cmd & 0x80) {   // ACMD<n> = CMD55 + CMD<n>
        cmd &= 0x7F;
        uint8_t r = send_cmd(CMD55, 0);
        if (r > 1) return r;
    }
    if (cmd != CMD12) {
        deselect();
        if (!select_card()) return 0xFF;
    }
    uint8_t frame[6] = {(uint8_t)(0x40 | cmd), (uint8_t)(arg >> 24), (uint8_t)(arg >> 16),
                        (uint8_t)(arg >> 8), (uint8_t)arg, 0x01};
    if (cmd == CMD0) frame[5] = 0x95;
    if (cmd == CMD8) frame[5] = 0x87;
    spi_write_blocking(SD_SPI, frame, 6);
    if (cmd == CMD12) xchg(0xFF);   // skip stuff byte
    uint8_t r;
    int n = 10;
    do {
        r = xchg(0xFF);
    } while ((r & 0x80) && --n);
    return r;
}

static bool rcv_block(uint8_t *buf, size_t len) {
    absolute_time_t t = make_timeout_time_ms(200);
    uint8_t token;
    do {
        token = xchg(0xFF);
    } while (token == 0xFF && !time_reached(t));
    if (token != 0xFE) return false;
    spi_read_blocking(SD_SPI, 0xFF, buf, len);
    xchg(0xFF);   // CRC
    xchg(0xFF);
    return true;
}

static bool xmit_block(const uint8_t *buf, uint8_t token) {
    if (!wait_ready(500)) return false;
    xchg(token);
    if (token == 0xFD) return true;   // stop transmission token
    spi_write_blocking(SD_SPI, buf, 512);
    xchg(0xFF);   // dummy CRC
    xchg(0xFF);
    uint8_t resp = xchg(0xFF);
    return (resp & 0x1F) == 0x05;
}

DSTATUS disk_initialize(BYTE pdrv) {
    if (pdrv) return STA_NOINIT;

    spi_init(SD_SPI, SPI_SLOW_HZ);
    gpio_set_function(PIN_SD_MISO, GPIO_FUNC_SPI);
    gpio_set_function(PIN_SD_SCK, GPIO_FUNC_SPI);
    gpio_set_function(PIN_SD_MOSI, GPIO_FUNC_SPI);
    gpio_pull_up(PIN_SD_MISO);
    gpio_init(PIN_SD_CS);
    gpio_set_dir(PIN_SD_CS, GPIO_OUT);
    cs_high();
    sleep_ms(2);
    for (int i = 0; i < 10; i++) xchg(0xFF);   // >= 74 clocks with CS high

    stat = STA_NOINIT;
    block_addr = false;
    bool ok = false;
    if (send_cmd(CMD0, 0) == 1) {
        absolute_time_t t = make_timeout_time_ms(1000);
        if (send_cmd(CMD8, 0x1AA) == 1) {   // SD v2
            uint8_t ocr[4];
            spi_read_blocking(SD_SPI, 0xFF, ocr, 4);
            if (ocr[2] == 0x01 && ocr[3] == 0xAA) {
                while (!time_reached(t) && send_cmd(ACMD41, 1u << 30)) {
                }
                if (!time_reached(t) && send_cmd(CMD58, 0) == 0) {
                    spi_read_blocking(SD_SPI, 0xFF, ocr, 4);
                    block_addr = ocr[0] & 0x40;
                    ok = true;
                }
            }
        } else {   // SD v1 or MMC
            uint8_t cmd = send_cmd(ACMD41, 0) <= 1 ? ACMD41 : 1;
            while (!time_reached(t) && send_cmd(cmd, 0)) {
            }
            ok = !time_reached(t) && send_cmd(CMD16, 512) == 0;
        }
    }
    deselect();
    if (ok) {
        spi_set_baudrate(SD_SPI, SPI_FAST_HZ);
        stat &= (DSTATUS)~STA_NOINIT;
    }
    return stat;
}

DSTATUS disk_status(BYTE pdrv) { return pdrv ? STA_NOINIT : stat; }

DRESULT disk_read(BYTE pdrv, BYTE *buff, LBA_t sector, UINT count) {
    if (pdrv || !count) return RES_PARERR;
    if (stat & STA_NOINIT) return RES_NOTRDY;
    uint32_t addr = block_addr ? (uint32_t)sector : (uint32_t)sector * 512u;
    if (count == 1) {
        if (send_cmd(CMD17, addr) == 0 && rcv_block(buff, 512)) count = 0;
    } else if (send_cmd(CMD18, addr) == 0) {
        do {
            if (!rcv_block(buff, 512)) break;
            buff += 512;
        } while (--count);
        send_cmd(CMD12, 0);
    }
    deselect();
    return count ? RES_ERROR : RES_OK;
}

DRESULT disk_write(BYTE pdrv, const BYTE *buff, LBA_t sector, UINT count) {
    if (pdrv || !count) return RES_PARERR;
    if (stat & STA_NOINIT) return RES_NOTRDY;
    uint32_t addr = block_addr ? (uint32_t)sector : (uint32_t)sector * 512u;
    if (count == 1) {
        if (send_cmd(CMD24, addr) == 0 && xmit_block(buff, 0xFE)) count = 0;
    } else if (send_cmd(CMD25, addr) == 0) {
        do {
            if (!xmit_block(buff, 0xFC)) break;
            buff += 512;
        } while (--count);
        if (!xmit_block(NULL, 0xFD)) count = 1;
    }
    deselect();
    return count ? RES_ERROR : RES_OK;
}

DRESULT disk_ioctl(BYTE pdrv, BYTE cmd, void *buff) {
    if (pdrv) return RES_PARERR;
    if (stat & STA_NOINIT) return RES_NOTRDY;
    DRESULT res = RES_ERROR;
    switch (cmd) {
    case CTRL_SYNC:
        if (select_card()) res = RES_OK;
        break;
    case GET_SECTOR_COUNT: {
        uint8_t csd[16];
        if (send_cmd(CMD9, 0) == 0 && rcv_block(csd, 16)) {
            if ((csd[0] >> 6) == 1) {   // CSD v2
                uint32_t c_size = ((uint32_t)(csd[7] & 0x3F) << 16) | ((uint32_t)csd[8] << 8) | csd[9];
                *(LBA_t *)buff = (LBA_t)(c_size + 1) << 10;
            } else {                    // CSD v1
                uint32_t n = (csd[5] & 15) + ((csd[10] & 128) >> 7) + ((csd[9] & 3) << 1) + 2;
                uint32_t c_size = (csd[8] >> 6) + ((uint32_t)csd[7] << 2) + ((uint32_t)(csd[6] & 3) << 10) + 1;
                *(LBA_t *)buff = (LBA_t)c_size << (n - 9);
            }
            res = RES_OK;
        }
        break;
    }
    case GET_BLOCK_SIZE:
        *(DWORD *)buff = 128;
        res = RES_OK;
        break;
    default:
        res = RES_PARERR;
    }
    deselect();
    return res;
}

// No RTC: FatFs uses FF_NORTC_* when FF_FS_NORTC = 1.
