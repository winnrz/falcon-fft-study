/*
 * refdrv.h
 *
 * Backend-agnostic interface to the reference FFT.
 *
 * The fpr type differs between backends -- a uint64_t of raw bits under
 * FALCON_FPEMU, a struct wrapping a double under FALCON_FPNATIVE -- so a
 * single translation unit cannot drive both.  Each backend gets its own
 * driver TU, and they meet here behind an interface made only of plain
 * integer types.
 */

#ifndef REFDRV_H__
#define REFDRV_H__

#include <stdint.h>

#include "kat.h"

struct ref_result {
	uint32_t cycles;    /* full multiply: 2x FFT, pointwise, iFFT   */
	int      margin;    /* worst precision headroom, in bits        */
	unsigned bad;       /* coefficients differing from the oracle   */
};

void ref_emu_case(const struct kat_case *kc, struct ref_result *r);
void ref_native_case(const struct kat_case *kc, struct ref_result *r);
void myfft_case(const struct kat_case *kc, struct ref_result *r);

/* Builds twiddle tables outside any timed region.  Returns 0 on failure. */
int myfft_prepare(void);

#endif
