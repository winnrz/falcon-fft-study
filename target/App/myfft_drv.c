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

#define MAX_LOGN    9
#define MAX_N       (1u << MAX_LOGN)
#define MARGIN_CAP  60

static double buf_a[MAX_N];
static double buf_b[MAX_N];

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
	unsigned t;

	for (t = 0; t < kat_ncases; t++) {
		if (!my_fft_init(kat_cases[t].logn)) {
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
