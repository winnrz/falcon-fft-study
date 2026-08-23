/*
 * bench_fftmul.c
 *
 * Performance measurement for FFT-based polynomial multiplication in
 * Z[X] / (X^N + 1): the from-scratch implementation in myfft.c against
 * the official Falcon reference in fft.c, with the exact integer
 * schoolbook multiply included to show what the FFT buys.
 *
 * Method
 * ------
 * Each operation is calibrated to run for at least MIN_BATCH_NS, then
 * timed over BATCHES repetitions with the minimum taken.  The minimum
 * rather than the mean: this is a deterministic computation with no
 * input-dependent branching, so run-to-run variation is scheduler and
 * cache noise, all of which adds time and none of which subtracts it.
 * The fastest observed run is the closest estimate of the true cost.
 *
 * Transforms are destructive, and re-transforming an already
 * transformed buffer makes the values grow without bound until they
 * reach infinity, which would both invalidate the measurement and risk
 * timing anomalies.  Every timed body therefore restores its input with
 * a memcpy first, and the cost of that memcpy alone is measured
 * separately and subtracted.  The correction is small -- a few percent
 * at n = 512 -- but it is measured rather than assumed.
 *
 * The two full-multiply paths are written out here rather than calling
 * my_poly_mul(), so that both implementations run against preallocated
 * buffers and neither is charged for a malloc the other does not make.
 *
 * Usage:
 *   ./bench_fftmul           human-readable table
 *   ./bench_fftmul --csv     CSV on stdout, for plotting
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "inner.h"
#include "myfft.h"

#define MAX_LOGN      10
#define MIN_BATCH_NS  2.0e7      /* 20 ms per calibrated batch  */
#define BATCHES       7          /* batches timed, minimum kept */

/* ==================================================================== */
/* Clock                                                                */
/* ==================================================================== */

static double
now_ns(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/* ==================================================================== */
/* Work context                                                         */
/* ==================================================================== */

/*
 * One set of buffers per size.  `src_*` hold pristine inputs; the timed
 * bodies copy from them into the `work_*` buffers before doing anything
 * destructive.  `tr_*` hold pre-transformed data, for timing the
 * inverse and the pointwise product in isolation.
 */
struct ctx {
	unsigned logn;
	size_t n;

	int64_t *ia, *ib;              /* integer operands            */

	double *src_a, *src_b;         /* mine: time domain           */
	double *tr_a, *tr_b;           /* mine: FFT domain            */
	double *wa, *wb;               /* mine: scratch               */

	fpr *fsrc_a, *fsrc_b;          /* reference: time domain      */
	fpr *ftr_a, *ftr_b;            /* reference: FFT domain       */
	fpr *fwa, *fwb;                /* reference: scratch          */

	__int128 *sb;                  /* schoolbook output           */
};

typedef void (*op_fn)(struct ctx *);

/* ==================================================================== */
/* Timed bodies                                                         */
/* ==================================================================== */

/* Baselines: the restore cost alone, subtracted from the ops below. */
static void
op_copy1(struct ctx *c)
{
	memcpy(c->wa, c->src_a, c->n * sizeof *c->wa);
}

static void
op_copy2(struct ctx *c)
{
	memcpy(c->wa, c->src_a, c->n * sizeof *c->wa);
	memcpy(c->wb, c->src_b, c->n * sizeof *c->wb);
}

static void
op_fcopy1(struct ctx *c)
{
	memcpy(c->fwa, c->fsrc_a, c->n * sizeof *c->fwa);
}

static void
op_fcopy2(struct ctx *c)
{
	memcpy(c->fwa, c->fsrc_a, c->n * sizeof *c->fwa);
	memcpy(c->fwb, c->fsrc_b, c->n * sizeof *c->fwb);
}

/* --- mine --- */

static void
op_my_fft(struct ctx *c)
{
	memcpy(c->wa, c->src_a, c->n * sizeof *c->wa);
	my_FFT(c->wa, c->logn);
}

static void
op_my_ifft(struct ctx *c)
{
	memcpy(c->wa, c->tr_a, c->n * sizeof *c->wa);
	my_iFFT(c->wa, c->logn);
}

static void
op_my_mulfft(struct ctx *c)
{
	memcpy(c->wa, c->tr_a, c->n * sizeof *c->wa);
	my_poly_mul_fft(c->wa, c->tr_b, c->logn);
}

static void
op_my_full(struct ctx *c)
{
	memcpy(c->wa, c->src_a, c->n * sizeof *c->wa);
	memcpy(c->wb, c->src_b, c->n * sizeof *c->wb);
	my_FFT(c->wa, c->logn);
	my_FFT(c->wb, c->logn);
	my_poly_mul_fft(c->wa, c->wb, c->logn);
	my_iFFT(c->wa, c->logn);
}

/* --- reference --- */

static void
op_ref_fft(struct ctx *c)
{
	memcpy(c->fwa, c->fsrc_a, c->n * sizeof *c->fwa);
	Zf(FFT)(c->fwa, c->logn);
}

static void
op_ref_ifft(struct ctx *c)
{
	memcpy(c->fwa, c->ftr_a, c->n * sizeof *c->fwa);
	Zf(iFFT)(c->fwa, c->logn);
}

static void
op_ref_mulfft(struct ctx *c)
{
	memcpy(c->fwa, c->ftr_a, c->n * sizeof *c->fwa);
	Zf(poly_mul_fft)(c->fwa, c->ftr_b, c->logn);
}

static void
op_ref_full(struct ctx *c)
{
	memcpy(c->fwa, c->fsrc_a, c->n * sizeof *c->fwa);
	memcpy(c->fwb, c->fsrc_b, c->n * sizeof *c->fwb);
	Zf(FFT)(c->fwa, c->logn);
	Zf(FFT)(c->fwb, c->logn);
	Zf(poly_mul_fft)(c->fwa, c->fwb, c->logn);
	Zf(iFFT)(c->fwa, c->logn);
}

/* --- the O(n^2) contrast --- */

static void
op_schoolbook(struct ctx *c)
{
	size_t n = c->n, i, j;
	__int128 *out = c->sb;

	for (i = 0; i < n; i++) {
		out[i] = 0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			__int128 p = (__int128)c->ia[i] * (__int128)c->ib[j];
			size_t k = i + j;

			if (k < n) {
				out[k] += p;
			} else {
				out[k - n] -= p;
			}
		}
	}
}

/* ==================================================================== */
/* Timing driver                                                        */
/* ==================================================================== */

static double
bench(op_fn op, struct ctx *c)
{
	long reps = 1;
	double best = 1e300;
	int b;

	/* Calibrate: grow reps until one batch is long enough to time. */
	for (;;) {
		double t0, el;
		long i;

		t0 = now_ns();
		for (i = 0; i < reps; i++) {
			op(c);
		}
		el = now_ns() - t0;
		if (el >= MIN_BATCH_NS || reps >= (1L << 30)) {
			break;
		}
		reps <<= 1;
	}

	for (b = 0; b < BATCHES; b++) {
		double t0, el;
		long i;

		t0 = now_ns();
		for (i = 0; i < reps; i++) {
			op(c);
		}
		el = now_ns() - t0;
		if (el < best) {
			best = el;
		}
	}
	return best / (double)reps;
}

/* ==================================================================== */
/* Setup                                                                */
/* ==================================================================== */

static uint64_t prng_state = 0x9E3779B97F4A7C15ULL;

static uint64_t
prng_next(void)
{
	prng_state ^= prng_state >> 12;
	prng_state ^= prng_state << 25;
	prng_state ^= prng_state >> 27;
	return prng_state * 0x2545F4914F6CDD1DULL;
}

static int
ctx_init(struct ctx *c, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t u;

	memset(c, 0, sizeof *c);
	c->logn = logn;
	c->n = n;

	c->ia = malloc(n * sizeof *c->ia);
	c->ib = malloc(n * sizeof *c->ib);
	c->src_a = malloc(n * sizeof *c->src_a);
	c->src_b = malloc(n * sizeof *c->src_b);
	c->tr_a = malloc(n * sizeof *c->tr_a);
	c->tr_b = malloc(n * sizeof *c->tr_b);
	c->wa = malloc(n * sizeof *c->wa);
	c->wb = malloc(n * sizeof *c->wb);
	c->fsrc_a = malloc(n * sizeof *c->fsrc_a);
	c->fsrc_b = malloc(n * sizeof *c->fsrc_b);
	c->ftr_a = malloc(n * sizeof *c->ftr_a);
	c->ftr_b = malloc(n * sizeof *c->ftr_b);
	c->fwa = malloc(n * sizeof *c->fwa);
	c->fwb = malloc(n * sizeof *c->fwb);
	c->sb = malloc(n * sizeof *c->sb);

	if (c->ia == NULL || c->ib == NULL || c->src_a == NULL
		|| c->src_b == NULL || c->tr_a == NULL || c->tr_b == NULL
		|| c->wa == NULL || c->wb == NULL || c->fsrc_a == NULL
		|| c->fsrc_b == NULL || c->ftr_a == NULL || c->ftr_b == NULL
		|| c->fwa == NULL || c->fwb == NULL || c->sb == NULL)
	{
		return 0;
	}

	/* Falcon-512-like magnitudes: small operand times large operand. */
	for (u = 0; u < n; u++) {
		c->ia[u] = (int64_t)(prng_next() % 51) - 25;
		c->ib[u] = (int64_t)(prng_next() % 24579) - 12289;
		c->src_a[u] = (double)c->ia[u];
		c->src_b[u] = (double)c->ib[u];
		c->fsrc_a[u] = fpr_of(c->ia[u]);
		c->fsrc_b[u] = fpr_of(c->ib[u]);
	}

	/* Pre-transformed copies, for timing iFFT and the product alone. */
	memcpy(c->tr_a, c->src_a, n * sizeof *c->tr_a);
	memcpy(c->tr_b, c->src_b, n * sizeof *c->tr_b);
	my_FFT(c->tr_a, logn);
	my_FFT(c->tr_b, logn);

	memcpy(c->ftr_a, c->fsrc_a, n * sizeof *c->ftr_a);
	memcpy(c->ftr_b, c->fsrc_b, n * sizeof *c->ftr_b);
	Zf(FFT)(c->ftr_a, logn);
	Zf(FFT)(c->ftr_b, logn);

	return 1;
}

static void
ctx_free(struct ctx *c)
{
	free(c->ia); free(c->ib);
	free(c->src_a); free(c->src_b);
	free(c->tr_a); free(c->tr_b);
	free(c->wa); free(c->wb);
	free(c->fsrc_a); free(c->fsrc_b);
	free(c->ftr_a); free(c->ftr_b);
	free(c->fwa); free(c->fwb);
	free(c->sb);
}

/* ==================================================================== */
/* Reporting                                                            */
/* ==================================================================== */

struct row {
	unsigned logn;
	size_t n;
	double my_fft, my_ifft, my_mulfft, my_full;
	double ref_fft, ref_ifft, ref_mulfft, ref_full;
	double schoolbook;
};

static double
nonneg(double x)
{
	return x < 0.0 ? 0.0 : x;
}

static void
print_table(const struct row *rows, int nrows)
{
	int i;

	printf("\nPer-operation cost, nanoseconds (minimum of %d batches)\n",
		BATCHES);
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	printf("%5s %6s | %10s %10s %10s | %10s %10s %10s\n",
		"logn", "n",
		"my FFT", "my iFFT", "my mul",
		"ref FFT", "ref iFFT", "ref mul");
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	for (i = 0; i < nrows; i++) {
		const struct row *r = &rows[i];

		printf("%5u %6zu | %10.1f %10.1f %10.1f | %10.1f %10.1f %10.1f\n",
			r->logn, r->n,
			r->my_fft, r->my_ifft, r->my_mulfft,
			r->ref_fft, r->ref_ifft, r->ref_mulfft);
	}

	printf("\nComplete integer multiply, nanoseconds\n");
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	printf("%5s %6s | %12s %12s %8s | %14s %10s\n",
		"logn", "n", "mine", "reference", "ratio",
		"schoolbook", "speedup");
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	for (i = 0; i < nrows; i++) {
		const struct row *r = &rows[i];

		printf("%5u %6zu | %12.1f %12.1f %8.2f | %14.1f %9.1fx\n",
			r->logn, r->n, r->my_full, r->ref_full,
			r->ref_full > 0.0 ? r->my_full / r->ref_full : 0.0,
			r->schoolbook,
			r->my_full > 0.0 ? r->schoolbook / r->my_full : 0.0);
	}

	/*
	 * Normalised cost.  If the implementation really is O(n log n)
	 * then time / (n * log2 n) is a constant, and any trend in this
	 * column is the memory hierarchy showing through rather than the
	 * algorithm.  The schoolbook column is normalised by n^2 for the
	 * same reason.
	 */
	printf("\nNormalised cost: FFT multiply / (n log2 n),"
		" schoolbook / n^2  [ns]\n");
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	printf("%5s %6s | %14s %14s | %16s\n",
		"logn", "n", "mine", "reference", "schoolbook");
	printf("%s\n",
		"-------------------------------------------------------"
		"---------------------------");
	for (i = 0; i < nrows; i++) {
		const struct row *r = &rows[i];
		double nlogn = (double)r->n * (double)r->logn;
		double nsq = (double)r->n * (double)r->n;

		printf("%5u %6zu | %14.4f %14.4f | %16.4f\n",
			r->logn, r->n,
			r->my_full / nlogn, r->ref_full / nlogn,
			r->schoolbook / nsq);
	}
}

static void
print_csv(const struct row *rows, int nrows)
{
	int i;

	printf("logn,n,my_fft,my_ifft,my_mulfft,my_full,"
		"ref_fft,ref_ifft,ref_mulfft,ref_full,schoolbook\n");
	for (i = 0; i < nrows; i++) {
		const struct row *r = &rows[i];

		printf("%u,%zu,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
			r->logn, r->n,
			r->my_fft, r->my_ifft, r->my_mulfft, r->my_full,
			r->ref_fft, r->ref_ifft, r->ref_mulfft, r->ref_full,
			r->schoolbook);
	}
}

/* ==================================================================== */

int
main(int argc, char **argv)
{
	int csv = (argc > 1 && strcmp(argv[1], "--csv") == 0);
	struct row rows[MAX_LOGN];
	int nrows = 0;
	unsigned logn;

	if (!csv) {
		printf("FFT polynomial multiplication -- performance\n");
		printf("  reference backend  %s\n",
			FALCON_FPNATIVE ? "FALCON_FPNATIVE (hardware binary64)"
			: "FALCON_FPEMU (software binary64)");
		printf("  reference AVX2     %s\n",
			FALCON_AVX2 ? "enabled" : "disabled");
		printf("  operands           Falcon-512-like, |a|<=25,"
			" |b|<=12289\n");
		printf("  timing             min of %d batches, each >= %.0f ms\n",
			BATCHES, MIN_BATCH_NS / 1e6);
		printf("  restore cost       measured and subtracted\n");
		fflush(stdout);
	}

	/* Build every twiddle table up front, outside the timed regions. */
	for (logn = 1; logn <= MAX_LOGN; logn++) {
		if (!my_fft_init(logn)) {
			fprintf(stderr, "table init failed at logn=%u\n", logn);
			return 1;
		}
	}

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		struct ctx c;
		struct row *r = &rows[nrows];
		double base1, base2, fbase1, fbase2;

		if (!ctx_init(&c, logn)) {
			fprintf(stderr, "allocation failed at logn=%u\n", logn);
			ctx_free(&c);
			return 1;
		}

		base1 = bench(op_copy1, &c);
		base2 = bench(op_copy2, &c);
		fbase1 = bench(op_fcopy1, &c);
		fbase2 = bench(op_fcopy2, &c);

		r->logn = logn;
		r->n = c.n;
		r->my_fft    = nonneg(bench(op_my_fft, &c)    - base1);
		r->my_ifft   = nonneg(bench(op_my_ifft, &c)   - base1);
		r->my_mulfft = nonneg(bench(op_my_mulfft, &c) - base1);
		r->my_full   = nonneg(bench(op_my_full, &c)   - base2);
		r->ref_fft    = nonneg(bench(op_ref_fft, &c)    - fbase1);
		r->ref_ifft   = nonneg(bench(op_ref_ifft, &c)   - fbase1);
		r->ref_mulfft = nonneg(bench(op_ref_mulfft, &c) - fbase1);
		r->ref_full   = nonneg(bench(op_ref_full, &c)   - fbase2);
		r->schoolbook = bench(op_schoolbook, &c);
		nrows++;

		if (!csv) {
			printf("  measured logn=%u (n=%zu)\n", logn, c.n);
			fflush(stdout);
		}
		ctx_free(&c);
	}

	if (csv) {
		print_csv(rows, nrows);
	} else {
		print_table(rows, nrows);

		/* Memory accounting: see also measure_memory.sh. */
		printf("\nWorking memory at logn=9 (n=512)\n");
		printf("%s\n", "--------------------------------------------");
		printf("  transform buffer, mine       %6zu bytes"
			"  (in place, n doubles)\n", (size_t)512 * sizeof(double));
		printf("  transform buffer, reference  %6zu bytes"
			"  (in place, n fpr)\n", (size_t)512 * sizeof(fpr));
		printf("  my twiddle tables, logn=9    %6zu bytes\n",
			(size_t)(2 * 256 + 2 * 128) * sizeof(double));
		printf("  my twiddle tables, all sizes %6zu bytes\n",
			my_fft_table_bytes());
		printf("  reference tables (static)    %6zu bytes"
			"  (fpr_gm_tab + fpr_p2_tab)\n",
			(size_t)2048 * sizeof(fpr) + (size_t)11 * sizeof(fpr));
		printf("  my_poly_mul scratch, logn=9  %6zu bytes"
			"  (heap, per call)\n",
			my_poly_mul_scratch_bytes(9));
		printf("\n  Neither transform allocates: fft.c and myfft.c both\n");
		printf("  work in place on caller-provided buffers.\n");
	}

	my_fft_free();
	return 0;
}
