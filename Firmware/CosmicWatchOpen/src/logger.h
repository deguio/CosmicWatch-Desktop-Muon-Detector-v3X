// microSD card: config.txt and buffered data files.
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "cwcore.h"

typedef void (*logger_log_fn)(const char *fmt, ...);

bool sd_mount(logger_log_fn log);
// Reads config.txt into cfg (creating it with cfg's values if missing).
bool sd_load_config(cw_config_t *cfg, logger_log_fn log);
// Opens <name>_<mode>_<NNN>.txt with NNN one above the highest existing number.
bool sd_open_run(const char *name, char mode, logger_log_fn log);
void sd_space(logger_log_fn log);

bool sd_active(void);
const char *sd_filename(void);
void sd_append(const char *s, size_t n);
// Writes the buffer every second (or when half full) and syncs every 10 s.
void sd_service(uint32_t now_ms);
uint32_t sd_dropped_bytes(void);
