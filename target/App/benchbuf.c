/*
 * benchbuf.c -- see benchbuf.h.
 */

#include <stdint.h>

#include "benchbuf.h"

static uint64_t slots[BENCH_SLOTS][BENCH_MAX_N];

void *
bench_slot(unsigned i)
{
	return (void *)slots[i % BENCH_SLOTS];
}
