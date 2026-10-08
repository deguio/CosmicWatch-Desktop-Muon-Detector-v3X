// The cable is straight-through, so line A of one board is line A of the other.
// Lines are only ever pulled low or released (internal pull-ups on both ends):
//
//   1. listen for a random 50-300 ms: if line A is held low, the partner is
//      announcing -> answer by pulling line B low and become the follower;
//   2. otherwise announce by pulling line A low for up to 2.5 s and wait for
//      line B to go low -> leader;
//   3. no answer -> single ("M") mode.
//
// Afterwards the leader signals events on line A and listens on line B, the
// follower does the opposite. Both detectors run the same acquisition code.
#include "coinc.h"

#include "board.h"
#include "hardware/gpio.h"
#include "pico/rand.h"
#include "pico/stdlib.h"

#define STABLE_MS        5u
#define ANNOUNCE_MS   2500u
#define RELEASE_MS    3000u

static void line_release(uint pin) { gpio_set_dir(pin, GPIO_IN); }
static void line_pull_low(uint pin) {
    gpio_put(pin, 0);
    gpio_set_dir(pin, GPIO_OUT);
}

// True if the line stays low for STABLE_MS (ignores glitches and event pulses).
static bool line_low_stable(uint pin) {
    if (gpio_get(pin)) return false;
    absolute_time_t until = make_timeout_time_ms(STABLE_MS);
    while (!time_reached(until)) {
        if (gpio_get(pin)) return false;
    }
    return true;
}

static bool wait_line_level(uint pin, bool level, uint32_t timeout_ms) {
    absolute_time_t until = make_timeout_time_ms(timeout_ms);
    while (!time_reached(until)) {
        if (gpio_get(pin) == level) return true;
    }
    return false;
}

coinc_role_t coinc_handshake(void) {
    const uint a = PIN_COINC_A, b = PIN_COINC_B;
    coinc_role_t r = {false, false, -1, (int)b};
    for (uint p = a; p <= b; p++) {
        gpio_init(p);
        gpio_put(p, 0);
        gpio_set_dir(p, GPIO_IN);
        gpio_pull_up(p);
    }
    sleep_ms(2);

    uint32_t listen_ms = 50u + get_rand_32() % 250u;
    absolute_time_t until = make_timeout_time_ms(listen_ms);
    while (!time_reached(until)) {
        if (line_low_stable(a)) {
            // follower
            line_pull_low(b);
            bool released = wait_line_level(a, true, RELEASE_MS);
            sleep_ms(STABLE_MS * 2);
            line_release(b);
            if (released) {
                r.found = true;
                r.leader = false;
                r.pin_out = (int)b;
                r.pin_in = (int)a;
            }
            return r;
        }
    }

    // leader candidate
    line_pull_low(a);
    bool answered = false;
    absolute_time_t ann = make_timeout_time_ms(ANNOUNCE_MS);
    while (!time_reached(ann)) {
        if (line_low_stable(b)) {
            answered = true;
            break;
        }
    }
    line_release(a);
    if (answered && wait_line_level(b, true, RELEASE_MS)) {
        r.found = true;
        r.leader = true;
        r.pin_out = (int)a;
        r.pin_in = (int)b;
    }
    return r;
}
