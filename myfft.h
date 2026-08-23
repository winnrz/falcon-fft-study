/*
 * myfft.h
 *
 * From-scratch implementation of FFT-based polynomial multiplication in
 * the ring Z[X] / (X^N + 1), written for this study.
 *
 * This is an independent implementation.  It shares no code with the
 * Falcon reference in fft.c: it uses the plain C `double` type rather
 * than the reference's `fpr` abstraction, and it generates its own
 * twiddle factors from cos/sin rather than reading fpr_gm_tab.
 *
 * It does deliberately reproduce the reference's FFT slot layout, so
 * that intermediate values can be compared coefficient by coefficient
 * against the reference implementation.  See the derivation in myfft.c.
 *
 * Array convention, matching the reference: a transformed polynomial of
 * n coefficients occupies n doubles, holding n/2 complex values as
 *
 *     f[j]        real part      of the j-th evaluation
 *     f[j + n/2]  imaginary part of the j-th evaluation
 */

#ifndef MYFFT_H__
#define MYFFT_H__

#include <stddef.h>
#include <stdint.h>

/* Largest transform this implementation supports (n = 2^MYFFT_MAX_LOGN). */
#define MYFFT_MAX_LOGN   12

/*
 * Build the twiddle tables for the given logn.  Idempotent: calling it
 * again for a logn already prepared is a no-op.  Returns 1 on success,
 * 0 on allocation failure or out-of-range logn.
 *
 * The transforms below call this automatically, so explicit use is only
 * needed to keep table construction out of a timed region.
 */
int my_fft_init(unsigned logn);

/* Release every table built so far.  Optional; for clean heap profiling. */
void my_fft_free(void);

/*
 * Total heap bytes currently held in twiddle tables, across all logn
 * prepared so far.  For the memory analysis.
 */
size_t my_fft_table_bytes(void);

/* Forward transform, in place.  Time domain -> FFT domain. */
void my_FFT(double *f, unsigned logn);

/* Inverse transform, in place.  FFT domain -> time domain. */
void my_iFFT(double *f, unsigned logn);

/* Pointwise product in the FFT domain: a <- a * b.  Both operands must
   already be transformed.  `a` and `b` may not overlap. */
void my_poly_mul_fft(double *a, const double *b, unsigned logn);

/* Pointwise sum in the FFT domain (equivalently, in the time domain --
   the transform is linear): a <- a + b. */
void my_poly_add(double *a, const double *b, unsigned logn);

/*
 * The complete operation: multiply two integer polynomials modulo
 * X^n + 1 by way of the FFT.  The result is left unrounded so callers
 * can measure how far it sits from the exact integer answer.
 *
 * `out` must have room for n doubles.  Uses one caller-sized scratch
 * buffer internally (see my_poly_mul_scratch_bytes).
 */
void my_poly_mul(const int64_t *a, const int64_t *b, double *out,
	unsigned logn);

/* Scratch bytes my_poly_mul allocates internally for a given logn. */
size_t my_poly_mul_scratch_bytes(unsigned logn);

#endif
