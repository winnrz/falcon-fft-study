/*
 * bench.h
 *
 * Per-operation cycle costs, backend-agnostic.
 */

#ifndef BENCH_H__
#define BENCH_H__

#include <stdint.h>

struct bench_row {
	uint32_t fft;      /* forward transform, in place              */
	uint32_t ifft;     /* inverse transform, in place              */
	uint32_t mul;      /* pointwise product in the FFT domain      */
	uint32_t full;     /* 2x FFT + pointwise + iFFT                */
	uint32_t restore;  /* the memcpy each timed body needs first   */
	uint32_t spread;   /* max - min of `full` over the repeats     */
};

void bench_emu(unsigned logn, struct bench_row *r);
void bench_native(unsigned logn, struct bench_row *r);
void bench_myfft(unsigned logn, struct bench_row *r);

/* Runs the sweep across every build and size, printing a table and CSV. */
void bench_run(void);

#endif
