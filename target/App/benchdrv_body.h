/*
 * benchdrv_body.h
 *
 * Per-operation timing for one fpr backend.  Included once per backend
 * TU, after refdrv_body.h has pulled in inner.h and cycles.h, with
 * BENCH_FN naming the function to define.
 *
 * Method
 * ------
 * Transforms are destructive, and re-transforming an already transformed
 * buffer makes the values grow without bound, so the input has to be
 * restored before each repetition.  bench_fftmul.c puts that memcpy
 * inside the timed body and subtracts its separately measured cost,
 * because on the host it cannot cheaply start and stop the clock around
 * a sub-microsecond region.  Here the restore sits outside the timed
 * region altogether: DWT reads cost one cycle and the region boundaries
 * are exact, so there is nothing to subtract and no correction to get
 * wrong.  The restore cost is still measured and reported, since it is
 * the quantity the host figures had subtracted from them.
 *
 * The host harness calibrates each operation to a minimum batch duration
 * and takes the minimum of several batches, because on a general-purpose
 * OS the noise sources -- scheduler preemption, cache eviction,
 * competing processes -- all add time and none subtract it.  None of
 * those exist here: bare metal, interrupts masked, zero flash wait
 * states.  So this harness keeps the minimum-of-repeats idea but drops
 * the calibration machinery, and reports the spread as evidence that the
 * determinism is real rather than assumed.
 */

#include <string.h>

#include "bench.h"
#include "benchbuf.h"

#define BENCH_REPEATS   3

/*
 * Deterministic operands at Falcon's magnitude, |coeff| <= 25.  Built
 * from a formula rather than a stored table so no further buffer is
 * needed.
 */
static void
bench_fill(fpr *a, fpr *b, size_t n)
{
	uint32_t st;
	size_t u;

	st = 0x12345678u;
	for (u = 0; u < n; u++) {
		st = st * 1664525u + 1013904223u;
		a[u] = fpr_of((int64_t)((st >> 16) % 51u) - 25);
		st = st * 1664525u + 1013904223u;
		b[u] = fpr_of((int64_t)((st >> 16) % 51u) - 25);
	}
}

void
BENCH_FN(unsigned logn, struct bench_row *r)
{
	fpr *wa, *wb, *sa, *sb;
	size_t n, bytes;
	uint32_t k, prim, t0, t1, d, best, worst;

	wa = (fpr *)bench_slot(0);
	wb = (fpr *)bench_slot(1);
	sa = (fpr *)bench_slot(2);
	sb = (fpr *)bench_slot(3);

	n = (size_t)1 << logn;
	bytes = n * sizeof *sa;

	bench_fill(sa, sb, n);

	/* --- cost of the restore itself ------------------------------ */
	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		prim = cyc_lock();
		t0 = cyc_read();
		memcpy(wa, sa, bytes);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->restore = best;

	/* --- forward transform --------------------------------------- */
	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		Zf(FFT)(wa, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->fft = best;

	/* --- complete multiply --------------------------------------- */
	best = 0xFFFFFFFFu;
	worst = 0;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		memcpy(wb, sb, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		Zf(FFT)(wa, logn);
		Zf(FFT)(wb, logn);
		Zf(poly_mul_fft)(wa, wb, logn);
		Zf(iFFT)(wa, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
		if (d > worst) {
			worst = d;
		}
	}
	r->full = best;
	r->spread = worst - best;

	/*
	 * The inverse transform and the pointwise product both need
	 * operands that are already in the FFT domain, so the sources
	 * are transformed in place here.  Everything above needed them
	 * in the time domain and has finished with them.
	 */
	Zf(FFT)(sa, logn);
	Zf(FFT)(sb, logn);

	/* --- inverse transform --------------------------------------- */
	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		Zf(iFFT)(wa, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->ifft = best;

	/* --- pointwise product --------------------------------------- */
	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		Zf(poly_mul_fft)(wa, sb, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->mul = best;
}
