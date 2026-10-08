#include "logger.h"

#include <stdio.h>
#include <string.h>

#include "ff.h"

#define BUF_SIZE        32768u
#define WRITE_CHUNK     8192u
#define WRITE_PERIOD_MS 1000u
#define SYNC_PERIOD_MS  10000u

static FATFS fs;
static FIL fil;
static bool mounted, active;
static char fname[48];
static char buf[BUF_SIZE];
static size_t used;
static uint32_t last_write_ms, last_sync_ms, dropped;
static logger_log_fn logf_;

bool sd_mount(logger_log_fn log) {
    logf_ = log;
    FRESULT fr = f_mount(&fs, "", 1);
    mounted = fr == FR_OK;
    if (!mounted) log("#  -> WARNING: failed to mount the microSD card (FatFs error %d)", (int)fr);
    else log("#    File system: %s", fs.fs_type == FS_EXFAT ? "exFAT" : fs.fs_type == FS_FAT32 ? "FAT32" : "FAT");
    return mounted;
}

bool sd_load_config(cw_config_t *cfg, logger_log_fn log) {
    if (!mounted) return false;
    FIL f;
    FRESULT fr = f_open(&f, "config.txt", FA_READ);
    if (fr == FR_NO_FILE) {
        static char text[2048];
        size_t n = cw_config_write(cfg, text, sizeof text);
        UINT bw = 0;
        fr = f_open(&f, "config.txt", FA_WRITE | FA_CREATE_NEW);
        if (fr == FR_OK) {
            fr = f_write(&f, text, (UINT)n, &bw);
            f_close(&f);
        }
        if (fr == FR_OK && bw == n) log("#    config.txt not found: created with default values");
        else log("#  -> ERROR: could not create config.txt (FatFs error %d)", (int)fr);
        return fr == FR_OK;
    }
    if (fr != FR_OK) {
        log("#  -> ERROR: could not open config.txt (FatFs error %d)", (int)fr);
        return false;
    }
    static char text[4096];
    UINT br = 0;
    fr = f_read(&f, text, sizeof text - 1, &br);
    f_close(&f);
    if (fr != FR_OK) {
        log("#  -> ERROR: could not read config.txt (FatFs error %d)", (int)fr);
        return false;
    }
    text[br] = 0;
    cw_config_parse(cfg, text, log);
    return true;
}

bool sd_open_run(const char *name, char mode, logger_log_fn log) {
    if (!mounted) return false;
    DIR dir;
    FILINFO fi;
    int maxn = 0;
    if (f_opendir(&dir, "") == FR_OK) {
        while (f_readdir(&dir, &fi) == FR_OK && fi.fname[0]) {
            if (fi.fattrib & AM_DIR) continue;
            int n = cw_parse_run_number(fi.fname, name);
            if (n > maxn) maxn = n;
        }
        f_closedir(&dir);
    }
    snprintf(fname, sizeof fname, "%s_%c_%03d.txt", name, mode, maxn + 1);
    FRESULT fr = f_open(&fil, fname, FA_WRITE | FA_CREATE_NEW);
    if (fr != FR_OK) {
        log("#  -> ERROR: could not create %s (FatFs error %d)", fname, (int)fr);
        fname[0] = 0;
        return false;
    }
    active = true;
    used = 0;
    return true;
}

void sd_space(logger_log_fn log) {
    DWORD nfree;
    FATFS *pfs;
    if (!mounted || f_getfree("", &nfree, &pfs) != FR_OK) return;
    uint64_t csz = (uint64_t)pfs->csize * 512u;
    uint64_t total = (uint64_t)(pfs->n_fatent - 2) * csz / 1024u;
    uint64_t freek = (uint64_t)nfree * csz / 1024u;
    log("#    Total space: %llu kB, Used: %llu kB, Free: %llu kB", (unsigned long long)total,
        (unsigned long long)(total - freek), (unsigned long long)freek);
}

bool sd_active(void) { return active; }
const char *sd_filename(void) { return fname; }
uint32_t sd_dropped_bytes(void) { return dropped; }

void sd_append(const char *s, size_t n) {
    if (!active) return;
    if (used + n > BUF_SIZE) {
        dropped += (uint32_t)n;
        return;
    }
    memcpy(buf + used, s, n);
    used += n;
}

static bool write_out(void) {
    UINT bw = 0;
    FRESULT fr = f_write(&fil, buf, (UINT)used, &bw);
    if (fr != FR_OK || bw != used) {
        // one retry from where it stopped
        size_t done = bw;
        UINT bw2 = 0;
        fr = f_write(&fil, buf + done, (UINT)(used - done), &bw2);
        if (fr != FR_OK || done + bw2 != used) {
            active = false;
            f_close(&fil);
            if (logf_) logf_("# WARNING: microSD write failed (FatFs error %d), logging stopped", (int)fr);
            return false;
        }
    }
    used = 0;
    return true;
}

void sd_service(uint32_t now_ms) {
    if (!active) return;
    if (used && (used >= WRITE_CHUNK || now_ms - last_write_ms >= WRITE_PERIOD_MS)) {
        last_write_ms = now_ms;
        if (!write_out()) return;
    }
    if (now_ms - last_sync_ms >= SYNC_PERIOD_MS) {
        last_sync_ms = now_ms;
        FRESULT fr = f_sync(&fil);
        if (fr != FR_OK && f_sync(&fil) != FR_OK) {
            active = false;
            f_close(&fil);
            if (logf_) logf_("# WARNING: microSD sync failed (FatFs error %d), logging stopped", (int)fr);
        }
    }
}
