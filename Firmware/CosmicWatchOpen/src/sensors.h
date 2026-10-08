// I2C peripherals: BMP280 (temperature/pressure) and MPU-6050 (accel/gyro).
#pragma once
#include <stdbool.h>
#include <stdint.h>

typedef struct {
    bool  ok;
    float temp_C, press_Pa;
} env_t;

typedef struct {
    bool  ok;
    float ax, ay, az;   // g
    float gx, gy, gz;   // deg/s
} imu_t;

typedef void (*sensors_log_fn)(const char *fmt, ...);

void i2c_bus_init(void);
bool i2c_probe(uint8_t addr);
void i2c_scan(sensors_log_fn log);

bool bmp280_init(sensors_log_fn log);
bool bmp280_read(env_t *out);

bool mpu6050_init(sensors_log_fn log);
bool mpu6050_read(imu_t *out);
