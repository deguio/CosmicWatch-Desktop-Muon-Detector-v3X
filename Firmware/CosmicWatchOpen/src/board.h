// Hardware description of the CosmicWatch v3X main PCB (schematic rev v2.1,
// "The Instruction Manual" Fig. B.2 and PCB/Circuit_Diagram.pdf).
#pragma once

#define FW_NAME     "CosmicWatch-open"
#define FW_VERSION  "1.0.0"

// ---------------------------------------------------------------- GPIO map
#define PIN_COINC_A       0   // CoincidentPin1, RJ45 pin 4
#define PIN_COINC_B       1   // CoincidentPin2, RJ45 pin 3
// GPIO2..5: SPI0 on the RJ45 (reserved for the expansion module, unused)
#define PIN_BUZZER        6   // PWM Buzzer (PWM slice 3 A)
// GPIO7: "TriggerCommand" on the manual's schematic revision, not connected
#define PIN_SD_MISO       8   // SPI1 RX
#define PIN_SD_CS         9   // SPI1 CSn (driven as GPIO)
#define PIN_SD_SCK       10   // SPI1 SCK
#define PIN_SD_MOSI      11   // SPI1 TX
#define PIN_LED_EVENT    12   // LED2 -> D5, 5 mm LED, flashes on every event
#define PIN_LED_COINC    13   // LED1 -> D2, 3 mm LED, flashes on coincident events
#define PIN_I2C_SDA      14   // I2C1: OLED (0x3C), BMP280 (0x76/77), MPU-6050 (0x68/69)
#define PIN_I2C_SCL      15
#define PIN_SD_DETECT    19   // microSD socket card-detect switch
#define PIN_THRESH_PWM   20   // PWM -> R21/R19 divider + C11 -> TriggerThreshold
#define PIN_PEAK_RESET   21   // TriggerReset: gate of Q2, discharges the peak-hold C5
#define PIN_TRIGGER      22   // comparator output (U5B), 1k/2k divided to 3.3 V
#define PIN_ADC_SIGNAL   26   // ADC0: buffered peak detector, x1.3/2.3 (R24/R27)
#define PIN_ADC_THRESH   27   // ADC1: trigger threshold read-back
#define PIN_ADC_HV       28   // ADC2: SiPM bias, divided by (22.1k+1k)/1k

#define ADC_CH_SIGNAL     0
#define ADC_CH_THRESH     1
#define ADC_CH_HV         2

#define I2C_PORT       i2c1
#define I2C_BAUD       400000
#define SD_SPI         spi1

// ---------------------------------------------------------------- analog
#define ADC_VREF_MV        2500.0f   // LM4040-2.5 on ADC_VREF
#define ADC_FULL_SCALE     4096.0f
#define HV_DIVIDER         23.1f     // (R7 + R12) / R12
#define HV_NOMINAL_V       30.2f     // MAX5026: 1.25 V * (1 + R2/R1)
#define AMP_BIAS_MV        24.75f    // VRef * R22 / (R22 + R10), comparator baseline
#define THRESH_PWM_HZ      1000000u  // threshold PWM frequency (as in the manual)

// SiPM (onsemi MicroFC-60035) breakdown voltage, datasheet typ. range 24.2-24.7 V
// at 21 C, temperature coefficient 21.5 mV/C.
#define SIPM_VBR_21C_V     24.45f
#define SIPM_VBR_TC_V      0.0215f
#define GAIN_REF_TEMP_C    25.0f

// Conversion ADC -> SiPM pulse height. 0.0706 mV/LSB reproduces the SiPM[mV]
// column written by the original firmware on v3X boards (overridable in config.txt).
#define SIPM_MV_PER_LSB_DEFAULT 0.0706f
