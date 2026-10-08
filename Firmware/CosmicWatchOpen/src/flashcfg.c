#include "flashcfg.h"

#include <string.h>

#include "cwcore.h"
#include "daq.h"
#include "hardware/flash.h"
#include "hardware/sync.h"
#include "pico/stdlib.h"

#define MAGIC 0x43574F31u   // "CWO1"
#define OFFSET (PICO_FLASH_SIZE_BYTES - FLASH_SECTOR_SIZE)

typedef struct {
    uint32_t magic;
    char     name[CW_NAME_MAX + 4];
    float    trigger_mV;
    uint32_t check;
} record_t;

static uint32_t checksum(const record_t *r) {
    const uint8_t *p = (const uint8_t *)r;
    uint32_t h = 2166136261u;   // FNV-1a over everything but the checksum
    for (size_t i = 0; i < offsetof(record_t, check); i++) h = (h ^ p[i]) * 16777619u;
    return h;
}

bool flashcfg_load(char *name, unsigned name_size, float *trigger_mV) {
    const record_t *r = (const record_t *)(XIP_BASE + OFFSET);
    if (r->magic != MAGIC || r->check != checksum(r)) return false;
    if (memchr(r->name, 0, sizeof r->name) == NULL) return false;
    strncpy(name, r->name, name_size - 1);
    name[name_size - 1] = 0;
    *trigger_mV = r->trigger_mV;
    return true;
}

bool flashcfg_save(const char *name, float trigger_mV) {
    static uint8_t page[FLASH_PAGE_SIZE];
    record_t r;
    memset(&r, 0, sizeof r);
    r.magic = MAGIC;
    strncpy(r.name, name, sizeof r.name - 1);
    r.trigger_mV = trigger_mV;
    r.check = checksum(&r);

    const record_t *cur = (const record_t *)(XIP_BASE + OFFSET);
    if (memcmp(cur, &r, sizeof r) == 0) return true;   // nothing to do

    memset(page, 0xFF, sizeof page);
    memcpy(page, &r, sizeof r);
    daq_pause();
    uint32_t irq = save_and_disable_interrupts();
    flash_range_erase(OFFSET, FLASH_SECTOR_SIZE);
    flash_range_program(OFFSET, page, FLASH_PAGE_SIZE);
    restore_interrupts(irq);
    daq_resume();
    return memcmp((const void *)(XIP_BASE + OFFSET), &r, sizeof r) == 0;
}
