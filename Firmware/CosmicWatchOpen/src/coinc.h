// Boot-time discovery of a partner detector on the RJ45 coincidence cable.
#pragma once
#include <stdbool.h>

typedef struct {
    bool found;
    bool leader;   // we drive line A, the partner drives line B
    int  pin_out;
    int  pin_in;
} coinc_role_t;

// Blocks for at most ~3 s. Both detectors must be (re)started within ~2 s.
coinc_role_t coinc_handshake(void);
