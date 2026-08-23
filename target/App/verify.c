/*
 * verify.c
 *
 * Runs all three builds of the FFT multiply against known-answer vectors
 * produced by tools/gen_kat.c, in a single firmware image so that every
 * number below is taken in one run at one clock setting.
 *
 *   1  reference on Falcon's emulated binary64, which config.h resolves
 *      to hand-written ARM assembly on a Cortex-M4.  Constant-time.
 *   2  reference on the C double type, which the single-precision FPU
 *      cannot execute, so GCC lowers it to libgcc soft-double calls.
 *   3  myfft.c, which always uses plain double and so lands on the same
 *      libgcc routines as build 2.
 *
 * Builds 2 and 3 share a float backend, so the difference between them
 * is FFT structure alone.  Builds 1 and 2 share an FFT, so the
 * difference between them is the price of constant-time arithmetic.
 * Comparing 1 against 3 conflates the two and is not meaningful.
 *
 * Two things are reported per row:
 *
 *   result   whether every rounded coefficient equals the exact integer
 *            the host oracle computed.  This is the pass/fail criterion.
 *
 *   margin   precision headroom in bits: how many times the worst
 *            observed deviation could be doubled before reaching the 0.5
 *            rounding cliff.  A case can pass on luck with a margin of
 *            1; the margin says whether it passed comfortably.  Bits
 *            rather than a scaled deviation because newlib-nano's printf
 *            supports neither %f nor %lld.
 *
 * The cycle figures are single unbatched measurements, taken here
 * because they cost nothing to collect.  They are indicative; the
 * calibrated timings come from the benchmark harness.
 */

#include <stdio.h>
#include <stdint.h>

#include "kat.h"
#include "refdrv.h"
#include "verify.h"

struct build {
	const char *name;
	void      (*run)(const struct kat_case *, struct ref_result *);
};

static const struct build builds[] = {
	{ "1 ref  emu+asm", ref_emu_case    },
	{ "2 ref  softdbl", ref_native_case },
	{ "3 myfft softdbl", myfft_case     },
};

#define NBUILDS  ((unsigned)(sizeof builds / sizeof builds[0]))

int
verify_run(void)
{
	unsigned b, t;
	int failures;
	uint32_t ref_cycles[NBUILDS];
	uint32_t lo[NBUILDS], hi[NBUILDS];

	for (b = 0; b < NBUILDS; b++) {
		ref_cycles[b] = 0;
		lo[b] = 0xFFFFFFFFu;
		hi[b] = 0;
	}

	printf("Correctness and cost, all builds, one run\r\n");
	printf("  %-16s %4s %5s %6s %7s %7s %11s\r\n",
		"build", "case", "logn", "bound", "result", "margin",
		"cycles");
	printf("  --------------------------------------"
		"---------------------------\r\n");

	failures = 0;
	for (b = 0; b < NBUILDS; b++) {
		for (t = 0; t < kat_ncases; t++) {
			struct ref_result r;

			builds[b].run(&kat_cases[t], &r);
			if (r.bad) {
				failures++;
			}

			/*
			 * Remember the operating-point case for the
			 * summary: Falcon-512 size, Falcon-512 operand
			 * magnitudes.
			 */
			if (kat_cases[t].logn == 9
				&& kat_cases[t].bound == 25
				&& ref_cycles[b] == 0)
			{
				ref_cycles[b] = r.cycles;
			}

			/*
			 * Spread across the n=512 cases.  Those cases
			 * differ only in their operand values, so any
			 * spread is data-dependent execution time --
			 * exactly what a constant-time implementation
			 * must not have.
			 */
			if (kat_cases[t].logn == 9) {
				if (r.cycles < lo[b]) {
					lo[b] = r.cycles;
				}
				if (r.cycles > hi[b]) {
					hi[b] = r.cycles;
				}
			}

			printf("  %-16s %4u %5u %6ld %7s %7d %11lu\r\n",
				t == 0 ? builds[b].name : "",
				t, kat_cases[t].logn,
				(long)kat_cases[t].bound,
				r.bad ? "FAIL" : "pass",
				r.margin, (unsigned long)r.cycles);
		}
	}

	printf("  --------------------------------------"
		"---------------------------\r\n");

	/*
	 * Summary at Falcon-512's operating point.  Ratios are scaled by
	 * 100 and printed as integers; there is no %f here.
	 */
	printf("\r\n  At logn=9, |coeff|<=25 -- n=512 multiply\r\n");
	for (b = 0; b < NBUILDS; b++) {
		unsigned long c = (unsigned long)ref_cycles[b];
		unsigned long pct;

		printf("    %-16s %9lu cycles  %4lu ms",
			builds[b].name, c, c / 24000UL);

		/* Ratio against build 2, the shared-backend reference. */
		if (b != 1 && ref_cycles[1] != 0) {
			pct = c * 100UL / (unsigned long)ref_cycles[1];
			printf("   %lu.%02lux vs build 2",
				pct / 100UL, pct % 100UL);
		}
		printf("\r\n");
	}

	printf("\r\n  Meaningful comparisons\r\n");
	if (ref_cycles[1] != 0) {
		unsigned long ct, st;

		ct = (unsigned long)ref_cycles[0] * 100UL
			/ (unsigned long)ref_cycles[1];
		st = (unsigned long)ref_cycles[2] * 100UL
			/ (unsigned long)ref_cycles[1];
		printf("    1 vs 2  same FFT, different backend"
			"    %lu.%02lux  price of constant-time\r\n",
			ct / 100UL, ct % 100UL);
		printf("    3 vs 2  same backend, different FFT"
			"    %lu.%02lux  cost of my structure\r\n",
			st / 100UL, st % 100UL);
	}

	printf("\r\n  Execution-time spread over the n=512 cases\r\n");
	printf("    (same size, different operand values)\r\n");
	for (b = 0; b < NBUILDS; b++) {
		printf("    %-16s %9lu cycles  %s\r\n",
			builds[b].name,
			(unsigned long)(hi[b] - lo[b]),
			hi[b] == lo[b] ? "data-independent"
				: "VARIES with data");
	}

	if (failures) {
		printf("\r\n  %d case(s) FAILED\r\n", failures);
	} else {
		printf("\r\n  all %u checks pass\r\n",
			NBUILDS * kat_ncases);
	}
	return failures;
}
