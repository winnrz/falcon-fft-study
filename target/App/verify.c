/*
 * verify.c
 *
 * Runs the Falcon reference FFT multiply on the target and compares it
 * against known-answer vectors produced by tools/gen_kat.c.
 *
 * This is the check that has to happen before any on-target timing is
 * quoted.  config.h autodetects FALCON_ASM_CORTEXM4 on a Cortex-M4
 * target, which forces the emulated binary64 backend and its
 * hand-written ARM assembly.  That assembly is not the plain-C emulation
 * the host harnesses verified: same intended semantics, different code.
 * Until it reproduces the exact integer products, a cycle count taken
 * from it measures an unknown quantity.
 *
 * Two things are reported per case:
 *
 *   result   whether every rounded coefficient equals the exact integer
 *            the oracle computed.  This is the pass/fail criterion.
 *
 *   margin   the worst-case precision headroom, in bits.  Rounding
 *            recovers the exact integer only while a coefficient sits
 *            closer than 0.5 to it, so 0.5 is the cliff.  The margin is
 *            how many times the largest observed deviation could be
 *            doubled before reaching that cliff.  A margin of 30 means
 *            the transform had about 30 bits of precision in hand; a
 *            margin of 1 would mean it passed by luck.  This is the same
 *            question the host harness's T9 asks, expressed in bits so
 *            it prints as a small integer -- newlib-nano's printf
 *            supports neither %f nor %lld, so a scaled fixed-point
 *            deviation cannot be printed here at all.
 *
 * The cycle figures are indicative only.  They are single measurements
 * with no calibration and no batching, included because they cost
 * nothing to collect here; the real timings come from the benchmark
 * harness.
 */

#include <stdio.h>
#include <stdint.h>

#include "inner.h"
#include "kat.h"
#include "cycles.h"
#include "verify.h"

#define MAX_LOGN   9
#define MAX_N      (1u << MAX_LOGN)

/*
 * Static rather than automatic: 8 KB of fpr does not belong on a 1 KB
 * stack, and static placement keeps the buffers out of any allocator so
 * the measurement is not charged for one.
 */
static fpr fa[MAX_N];
static fpr fb[MAX_N];

/* |x|, via the only primitives the fpr backend exposes. */
static fpr
fpr_absval(fpr x)
{
	return fpr_lt(x, fpr_of(0)) ? fpr_neg(x) : x;
}

/*
 * Precision headroom in bits: how many doublings separate d from the 0.5
 * rounding cliff.  Capped, since a deviation of exactly zero would
 * otherwise never reach the limit.
 */
#define MARGIN_CAP  60

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

static int
run_case(const struct kat_case *kc, unsigned idx)
{
	size_t n, u;
	uint32_t prim, t0, t1, cycles;
	int worst_margin, m;
	unsigned bad;

	n = (size_t)1 << kc->logn;

	for (u = 0; u < n; u++) {
		fa[u] = fpr_of(kc->a[u]);
		fb[u] = fpr_of(kc->b[u]);
	}

	prim = cyc_lock();
	t0 = cyc_read();
	Zf(FFT)(fa, kc->logn);
	Zf(FFT)(fb, kc->logn);
	Zf(poly_mul_fft)(fa, fb, kc->logn);
	Zf(iFFT)(fa, kc->logn);
	t1 = cyc_read();
	cyc_unlock(prim);
	cycles = t1 - t0;

	bad = 0;
	worst_margin = MARGIN_CAP;
	for (u = 0; u < n; u++) {
		int64_t got;

		got = (int64_t)fpr_rint(fa[u]);
		if (got != kc->expect[u]) {
			bad++;
		}

		/* Worst (smallest) precision headroom over the case. */
		m = margin_bits(fpr_absval(fpr_sub(fa[u], fpr_of(got))));
		if (m < worst_margin) {
			worst_margin = m;
		}
	}

	printf("  %4u %5u %6ld  %7s %8d %12lu\r\n",
		idx, kc->logn, (long)kc->bound,
		bad ? "FAIL" : "pass",
		worst_margin, (unsigned long)cycles);

	if (bad) {
		printf("       %u of %lu coefficients wrong\r\n",
			bad, (unsigned long)n);
	}
	return bad ? 1 : 0;
}

int
verify_reference(void)
{
	unsigned t;
	int failures;

	printf("Reference FFT multiply vs exact integer oracle\r\n");
	printf("  backend: FALCON_FPEMU=%d  FPNATIVE=%d  ASM_M4=%d\r\n",
		FALCON_FPEMU, FALCON_FPNATIVE, FALCON_ASM_CORTEXM4);
	printf("  %4s %5s %6s  %7s %8s %12s\r\n",
		"case", "logn", "bound", "result", "margin", "cycles");
	printf("  ---------------------------------------"
		"----------------------\r\n");

	failures = 0;
	for (t = 0; t < kat_ncases; t++) {
		failures += run_case(&kat_cases[t], t);
	}

	printf("  ---------------------------------------"
		"----------------------\r\n");
	if (failures) {
		printf("  %d case(s) FAILED\r\n", failures);
	} else {
		printf("  all %u cases pass\r\n", kat_ncases);
	}
	return failures;
}
