/*
 * myfft_drv.c -- build 3: the from-scratch implementation.
 *
 * myfft.c uses the plain C double type throughout and never consults
 * config.h, so it is unaffected by the FALCON_FPEMU / FALCON_FPNATIVE
 * selection: on this target it always lands on libgcc's soft-double
 * routines.  That is the same backend build 2 uses, which is what makes
 * builds 2 and 3 a fair comparison of FFT structure.
 *
 * my_poly_mul() is deliberately not called.  It allocates scratch on
 * every invocation, and that malloc would sit inside the timed region.
 * The transform steps are driven individually here, exactly as
 * bench_fftmul.c does on the host, so the two measure the same thing.
 */

#include <stdint.h>
#include <stddef.h>

#include "myfft.h"
#include "refdrv.h"
#include "cycles.h"
#include "bench.h"
#include "benchbuf.h"

#include <string.h>

#define MARGIN_CAP    60
#define BENCH_REPEATS 3

/*
 * Round half away from zero.  fpr_rint rounds half to even, so the two
 * differ on exact halves; a coefficient landing exactly on .5 would mean
 * the transform had already lost the value, which the margin column
 * would show long before this mattered.
 */
static int64_t
round_int(double x)
{
	return (int64_t)(x < 0.0 ? x - 0.5 : x + 0.5);
}

static int
margin_bits_d(double d)
{
	int bits;

	if (d < 0.0) {
		d = -d;
	}
	if (!(d < 0.5)) {
		return 0;
	}
	for (bits = 0; bits < MARGIN_CAP && d < 0.5; bits++) {
		d *= 2.0;
	}
	return bits;
}

int
myfft_prepare(void)
{
	unsigned logn;

	/*
	 * Every size the benchmark will touch, built here so table
	 * construction and its malloc stay out of every timed region.
	 */
	for (logn = 1; logn <= BENCH_MAX_LOGN; logn++) {
		if (!my_fft_init(logn)) {
			return 0;
		}
	}
	return 1;
}

void
myfft_case(const struct kat_case *kc, struct ref_result *r)
{
	size_t n, u;
	uint32_t prim, t0, t1;
	double *buf_a, *buf_b;

	buf_a = (double *)bench_slot(0);
	buf_b = (double *)bench_slot(1);

	n = (size_t)1 << kc->logn;
	for (u = 0; u < n; u++) {
		buf_a[u] = (double)kc->a[u];
		buf_b[u] = (double)kc->b[u];
	}

	prim = cyc_lock();
	t0 = cyc_read();
	my_FFT(buf_a, kc->logn);
	my_FFT(buf_b, kc->logn);
	my_poly_mul_fft(buf_a, buf_b, kc->logn);
	my_iFFT(buf_a, kc->logn);
	t1 = cyc_read();
	cyc_unlock(prim);

	r->cycles = t1 - t0;
	r->bad = 0;
	r->margin = MARGIN_CAP;

	for (u = 0; u < n; u++) {
		int64_t got;
		int m;

		got = round_int(buf_a[u]);
		if (got != kc->expect[u]) {
			r->bad++;
		}
		m = margin_bits_d(buf_a[u] - (double)got);
		if (m < r->margin) {
			r->margin = m;
		}
	}
}

/* ==================================================================== */
/* Per-operation timing -- mirrors benchdrv_body.h exactly              */
/* ==================================================================== */

static void
bench_fill_d(double *a, double *b, size_t n)
{
	uint32_t st;
	size_t u;

	st = 0x12345678u;
	for (u = 0; u < n; u++) {
		st = st * 1664525u + 1013904223u;
		a[u] = (double)((int32_t)((st >> 16) % 51u) - 25);
		st = st * 1664525u + 1013904223u;
		b[u] = (double)((int32_t)((st >> 16) % 51u) - 25);
	}
}

void
bench_myfft(unsigned logn, struct bench_row *r)
{
	double *wa, *wb, *sa, *sb;
	size_t n, bytes;
	uint32_t k, prim, t0, t1, d, best, worst;

	wa = (double *)bench_slot(0);
	wb = (double *)bench_slot(1);
	sa = (double *)bench_slot(2);
	sb = (double *)bench_slot(3);

	n = (size_t)1 << logn;
	bytes = n * sizeof *sa;

	bench_fill_d(sa, sb, n);

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

	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		my_FFT(wa, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->fft = best;

	best = 0xFFFFFFFFu;
	worst = 0;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		memcpy(wb, sb, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		my_FFT(wa, logn);
		my_FFT(wb, logn);
		my_poly_mul_fft(wa, wb, logn);
		my_iFFT(wa, logn);
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

	my_FFT(sa, logn);
	my_FFT(sb, logn);

	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		my_iFFT(wa, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->ifft = best;

	best = 0xFFFFFFFFu;
	for (k = 0; k < BENCH_REPEATS; k++) {
		memcpy(wa, sa, bytes);
		prim = cyc_lock();
		t0 = cyc_read();
		my_poly_mul_fft(wa, sb, logn);
		t1 = cyc_read();
		cyc_unlock(prim);
		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	r->mul = best;
}
