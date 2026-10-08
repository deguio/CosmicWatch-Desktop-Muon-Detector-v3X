// Detector name and threshold kept in the last flash sector, so a detector
// without microSD card keeps its identity (needed by import_data.py, which
// requires distinct names on each USB port).
#pragma once
#include <stdbool.h>

bool flashcfg_load(char *name, unsigned name_size, float *trigger_mV);
// Pauses the acquisition core while erasing/programming.
bool flashcfg_save(const char *name, float trigger_mV);
