#include "sensors.h"

#include "board.h"
#include "cwcore.h"
#include "hardware/gpio.h"
#include "hardware/i2c.h"
#include "pico/stdlib.h"

#define I2C_TIMEOUT_US 5000u

static uint8_t bmp_addr, mpu_addr;
static cw_bmp280_calib_t bmp_cal;

void i2c_bus_init(void) {
    i2c_init(I2C_PORT, I2C_BAUD);
    gpio_set_function(PIN_I2C_SDA, GPIO_FUNC_I2C);
    gpio_set_function(PIN_I2C_SCL, GPIO_FUNC_I2C);
    gpio_pull_up(PIN_I2C_SDA);
    gpio_pull_up(PIN_I2C_SCL);
}

bool i2c_probe(uint8_t addr) {
    uint8_t d;
    return i2c_read_timeout_us(I2C_PORT, addr, &d, 1, false, I2C_TIMEOUT_US) == 1;
}

static bool reg_write(uint8_t addr, uint8_t reg, uint8_t val) {
    uint8_t b[2] = {reg, val};
    return i2c_write_timeout_us(I2C_PORT, addr, b, 2, false, I2C_TIMEOUT_US) == 2;
}

static bool reg_read(uint8_t addr, uint8_t reg, uint8_t *dst, size_t n) {
    if (i2c_write_timeout_us(I2C_PORT, addr, &reg, 1, true, I2C_TIMEOUT_US) != 1) return false;
    return i2c_read_timeout_us(I2C_PORT, addr, dst, n, false, I2C_TIMEOUT_US * 4) == (int)n;
}

void i2c_scan(sensors_log_fn log) {
    log("# Scanning I2C ...");
    for (uint8_t a = 0x08; a < 0x78; a++) {
        if (!i2c_probe(a)) continue;
        const char *what = "(unknown device)";
        if (a == 0x3C || a == 0x3D) what = "(OLED display: SSD1306)";
        else if (a == 0x76 || a == 0x77) what = "(Temp/Pressure sensor: BMP280/BME280)";
        else if (a == 0x68 || a == 0x69) what = "(Accel/Gyro: MPU-6050)";
        else if (a >= 0x50 && a <= 0x57) what = "(I2C EEPROM 24xx series)";
        log("#  -- Device ACK at 0x%02X  %s", a, what);
    }
}

// ---------------------------------------------------------------- BMP280

bool bmp280_init(sensors_log_fn log) {
    static const uint8_t addrs[] = {0x76, 0x77};
    for (unsigned i = 0; i < 2; i++) {
        uint8_t id;
        if (!reg_read(addrs[i], 0xD0, &id, 1)) continue;
        if (id != 0x58 && id != 0x60) {
            log("#    BMP280: unexpected chip ID 0x%02X at 0x%02X", id, addrs[i]);
            continue;
        }
        bmp_addr = addrs[i];
        uint8_t raw[24];
        if (!reg_write(bmp_addr, 0xE0, 0xB6)) break;   // soft reset
        sleep_ms(5);
        if (!reg_read(bmp_addr, 0x88, raw, 24)) {
            log("#    BMP280: calibration read failed");
            break;
        }
        cw_bmp280_parse_calib(&bmp_cal, raw);
        // config: standby 500 ms, IIR filter x16; ctrl_meas: T x2, P x16, normal mode
        if (!reg_write(bmp_addr, 0xF5, 0x90) || !reg_write(bmp_addr, 0xF4, 0x57)) break;
        log("#    BMP280 chip ID = 0x%02X at address 0x%02X", id, bmp_addr);
        return true;
    }
    bmp_addr = 0;
    return false;
}

bool bmp280_read(env_t *out) {
    uint8_t d[6];
    if (!bmp_addr || !reg_read(bmp_addr, 0xF7, d, 6)) return false;
    int32_t adc_P = (int32_t)((d[0] << 12) | (d[1] << 4) | (d[2] >> 4));
    int32_t adc_T = (int32_t)((d[3] << 12) | (d[4] << 4) | (d[5] >> 4));
    if (adc_T == 0x80000) return false;   // no conversion yet
    cw_bmp280_compensate(&bmp_cal, adc_T, adc_P, &out->temp_C, &out->press_Pa);
    out->ok = true;
    return true;
}

// ---------------------------------------------------------------- MPU-6050

bool mpu6050_init(sensors_log_fn log) {
    static const uint8_t addrs[] = {0x68, 0x69};
    for (unsigned i = 0; i < 2; i++) {
        uint8_t id;
        if (!reg_read(addrs[i], 0x75, &id, 1)) continue;
        mpu_addr = addrs[i];
        if (id != 0x68) log("#    MPU-6050: WHO_AM_I = 0x%02X (clone?), using it anyway", id);
        bool ok = reg_write(mpu_addr, 0x6B, 0x80);          // device reset
        sleep_ms(100);
        ok = ok && reg_write(mpu_addr, 0x6B, 0x01);         // wake, PLL on X gyro
        ok = ok && reg_write(mpu_addr, 0x19, 9);            // 100 Hz sample rate
        ok = ok && reg_write(mpu_addr, 0x1A, 0x03);         // DLPF 44 Hz
        ok = ok && reg_write(mpu_addr, 0x1B, 0x00);         // +-250 deg/s
        ok = ok && reg_write(mpu_addr, 0x1C, 0x00);         // +-2 g
        if (!ok) break;
        log("#    MPU-6050 ID = 0x%02X at address 0x%02X", id, mpu_addr);
        return true;
    }
    mpu_addr = 0;
    return false;
}

bool mpu6050_read(imu_t *out) {
    uint8_t d[14];
    if (!mpu_addr || !reg_read(mpu_addr, 0x3B, d, 14)) return false;
#define S16(i) ((int16_t)((d[i] << 8) | d[(i) + 1]))
    out->ax = S16(0) / 16384.0f;
    out->ay = S16(2) / 16384.0f;
    out->az = S16(4) / 16384.0f;
    out->gx = S16(8) / 131.0f;
    out->gy = S16(10) / 131.0f;
    out->gz = S16(12) / 131.0f;
#undef S16
    out->ok = true;
    return true;
}
