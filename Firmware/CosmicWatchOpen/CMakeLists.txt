cmake_minimum_required(VERSION 3.13)

set(PICO_BOARD pico)
include(pico_sdk_import.cmake)

project(cosmicwatch_open C CXX ASM)
set(CMAKE_C_STANDARD 11)
pico_sdk_init()

add_executable(cosmicwatch_open
    src/main.c
    src/daq.c
    src/coinc.c
    src/threshold.c
    src/cwcore.c
    src/sensors.c
    src/oled.c
    src/logger.c
    src/sdcard.c
    src/flashcfg.c
    lib/fatfs/ff.c
    lib/fatfs/ffunicode.c
)

pico_generate_pio_header(cosmicwatch_open ${CMAKE_CURRENT_LIST_DIR}/src/daq.pio)

target_include_directories(cosmicwatch_open PRIVATE src lib/fatfs)

target_compile_definitions(cosmicwatch_open PRIVATE
    SYS_CLK_HZ=125000000            # PIO timing of the coincidence window assumes 125 MHz
    CFG_TUD_CDC_TX_BUFSIZE=1024
    CFG_TUD_CDC_RX_BUFSIZE=256
    PICO_STDIO_USB_STDOUT_TIMEOUT_US=20000
    USBD_PRODUCT="CosmicWatch v3X"
)

target_compile_options(cosmicwatch_open PRIVATE
    $<$<COMPILE_LANGUAGE:C>:-Wall -Wextra -Wno-unused-parameter>
)
set_source_files_properties(lib/fatfs/ff.c lib/fatfs/ffunicode.c PROPERTIES COMPILE_OPTIONS "-w")

target_link_libraries(cosmicwatch_open
    pico_stdlib
    pico_multicore
    pico_unique_id
    pico_rand
    hardware_adc
    hardware_dma
    hardware_flash
    hardware_i2c
    hardware_pio
    hardware_pwm
    hardware_spi
    hardware_watchdog
)

pico_enable_stdio_usb(cosmicwatch_open 1)
pico_enable_stdio_uart(cosmicwatch_open 0)

pico_set_program_name(cosmicwatch_open "CosmicWatch v3X open firmware")
pico_set_program_version(cosmicwatch_open "1.0.0")

pico_add_extra_outputs(cosmicwatch_open)
