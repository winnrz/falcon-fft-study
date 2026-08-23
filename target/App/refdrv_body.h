/*
 * refdrv_body.h
 *
 * Shared implementation for both reference drivers.  Included exactly
 * once per driver TU, after that TU has selected its backend, with
 * REF_FN naming the function to define.  Written as a header rather
 * than duplicated so the two backends cannot drift apart: any
 * difference in the measured numbers is then the backend, not the
 * harness around it.
 */

#include <stdint.h>

#include "inner.h"
#include "refdrv.h"
#include "cycles.h"
#include "benchbuf.h"

#define MARGIN_CAP 60

static fpr
absval(fpr x)
{
	return fpr_lt(x, fpr_of(0)) ? fpr_neg(x) : x;
}

/* Doublings between d and the 0.5 rounding cliff.  See verify.c. */
static int
margin_bits(fpr d)
{
	fpr half;
	int bits;

	half = fpr_half(fpr_of(1));
	if (!fpr_lt(d, half)) {
		return 0;
	}
	for (bits = 0; bits < MARGIN_CAP && fpr_lt(d, half); bits++) {
		d = fpr_double(d);
	}
	return bits;
}

void
REF_FN(const struct kat_case *kc, struct ref_result *r)
{
	size_t n, u;
	uint32_t prim, t0, t1;
	fpr *buf_a, *buf_b;

	buf_a = (fpr *)bench_slot(0);
	buf_b = (fpr *)bench_slot(1);

	n = (size_t)1 << kc->logn;
	for (u = 0; u < n; u++) {
		buf_a[u] = fpr_of(kc->a[u]);
		buf_b[u] = fpr_of(kc->b[u]);
	}

	prim = cyc_lock();
	t0 = cyc_read();
	Zf(FFT)(buf_a, kc->logn);
	Zf(FFT)(buf_b, kc->logn);
	Zf(poly_mul_fft)(buf_a, buf_b, kc->logn);
	Zf(iFFT)(buf_a, kc->logn);
	t1 = cyc_read();
	cyc_unlock(prim);

	r->cycles = t1 - t0;
	r->bad = 0;
	r->margin = MARGIN_CAP;

	for (u = 0; u < n; u++) {
		int64_t got;
		int m;

		got = (int64_t)fpr_rint(buf_a[u]);
		if (got != kc->expect[u]) {
			r->bad++;
		}
		m = margin_bits(absval(fpr_sub(buf_a[u], fpr_of(got))));
		if (m < r->margin) {
			r->margin = m;
		}
	}
}
