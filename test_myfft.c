/*
 * test_myfft.c
 *
 * Verification harness for the from-scratch implementation in myfft.c.
 *
 * It runs the same ten properties as test_fftmul.c, with the same PRNG
 * seeds and the same coefficient ranges, so the two outputs can be set
 * side by side.  It adds one group the reference harness cannot have:
 *
 *   T11 cross-check   my_FFT vs the official reference, slot by slot
 *
 * T11 is what the project objectives call for -- "compared with the
 * official reference implementation".  It is stronger than agreeing on
 * products: it checks that the two implementations agree on every
 * intermediate value in the FFT domain, which only holds if the slot
 * layouts match as well as the arithmetic.
 *
 * The exact integer schoolbook oracle is duplicated from test_fftmul.c
 * rather than shared, deliberately: test_fftmul.c is a verified
 * artefact whose results are quoted in the write-up, and it is left
 * untouched.
 *
 * Exit status is 0 if every check passes, 1 otherwise.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include "inner.h"
#include "myfft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define MAX_LOGN   10
#define MAX_N      (1u << MAX_LOGN)

static int total_checks = 0;
static int total_failures = 0;

/* ==================================================================== */
/* Small utilities                                                      */
/* ==================================================================== */

/* Extract the binary64 value behind an fpr, for the T11 cross-check. */
static double
fpr_to_double(fpr x)
{
#if FALCON_FPEMU
	uint64_t bits = (uint64_t)x;
	double d;

	memcpy(&d, &bits, sizeof d);
	return d;
#else
	return x.v;
#endif
}

/* Deterministic PRNG (xorshift64*), matching test_fftmul.c exactly. */
static uint64_t prng_state = 0x9E3779B97F4A7C15ULL;

static void
prng_seed(uint64_t s)
{
	prng_state = s ? s : 0x9E3779B97F4A7C15ULL;
}

static uint64_t
prng_next(void)
{
	prng_state ^= prng_state >> 12;
	prng_state ^= prng_state << 25;
	prng_state ^= prng_state >> 27;
	return prng_state * 0x2545F4914F6CDD1DULL;
}

static int64_t
rand_coeff(int64_t B)
{
	return (int64_t)(prng_next() % (uint64_t)(2 * B + 1)) - B;
}

static void
rand_poly(int64_t *a, size_t n, int64_t B)
{
	size_t u;

	for (u = 0; u < n; u++) {
		a[u] = rand_coeff(B);
	}
}

static size_t
bitrev(size_t x, unsigned bits)
{
	size_t r = 0;

	while (bits--) {
		r = (r << 1) | (x & 1);
		x >>= 1;
	}
	return r;
}

static double
d_abs(double x)
{
	return x < 0.0 ? -x : x;
}

static __int128
i_abs(__int128 x)
{
	return x < 0 ? -x : x;
}

static double
i128_to_double(__int128 x)
{
	int neg = x < 0;
	double d;

	if (neg) {
		x = -x;
	}
	d = (double)(uint64_t)(x >> 64) * 18446744073709551616.0
		+ (double)(uint64_t)x;
	return neg ? -d : d;
}

/* Round a double to the nearest integer, ties away from zero. */
static __int128
d_rint(double x)
{
	return (__int128)(x < 0.0 ? x - 0.5 : x + 0.5);
}

static void
check(const char *name, int ok, const char *fmt, ...)
{
	va_list ap;

	total_checks++;
	if (ok) {
		printf("  PASS  %-34s", name);
	} else {
		total_failures++;
		printf("  FAIL  %-34s", name);
	}
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	printf("\n");
}

static void
banner(const char *s)
{
	printf("\n%s\n", s);
}

/* ==================================================================== */
/* Reference oracle: exact integer schoolbook multiply mod X^n + 1      */
/* ==================================================================== */

/*
 * c[k] = sum_{i+j=k} a[i]b[j]  -  sum_{i+j=k+n} a[i]b[j]
 *
 * Exact in __int128, no floating point anywhere.
 */
static void
schoolbook_mul(const int64_t *a, const int64_t *b, __int128 *c, size_t n)
{
	size_t i, j;

	for (i = 0; i < n; i++) {
		c[i] = 0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			__int128 p = (__int128)a[i] * (__int128)b[j];
			size_t k = i + j;

			if (k < n) {
				c[k] += p;
			} else {
				c[k - n] -= p;
			}
		}
	}
}

/* The operation under test, left unrounded. */
static void
mine_mul(const int64_t *a, const int64_t *b, double *out, unsigned logn)
{
	my_poly_mul(a, b, out, logn);
}

/* ==================================================================== */
/* T1  round trip: iFFT(FFT(a)) == a                                    */
/* ==================================================================== */

static void
t1_roundtrip(void)
{
	unsigned logn;

	banner("T1  round trip   iFFT(FFT(a)) == a");
	prng_seed(1001);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		double *f = malloc(n * sizeof *f);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 64; iter++) {
			rand_poly(a, n, 512);
			for (u = 0; u < n; u++) {
				f[u] = (double)a[u];
			}
			my_FFT(f, logn);
			my_iFFT(f, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(f[u] - (double)a[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)d_rint(f[u]) != a[u]) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "max|err| = %.3e", maxerr);
		}
		free(a); free(f);
	}
}

/* ==================================================================== */
/* T2  linearity: FFT(a + b) == FFT(a) + FFT(b)                         */
/* ==================================================================== */

static void
t2_linearity(void)
{
	unsigned logn;

	banner("T2  linearity    FFT(a+b) == FFT(a) + FFT(b)");
	prng_seed(2002);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		double *fa = malloc(n * sizeof *fa);
		double *fb = malloc(n * sizeof *fb);
		double *fs = malloc(n * sizeof *fs);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 64; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			for (u = 0; u < n; u++) {
				fa[u] = (double)a[u];
				fb[u] = (double)b[u];
				fs[u] = (double)(a[u] + b[u]);
			}
			my_FFT(fa, logn);
			my_FFT(fb, logn);
			my_FFT(fs, logn);
			my_poly_add(fa, fb, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(fa[u] - fs[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if (e > 1e-9) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "max|err| = %.3e", maxerr);
		}
		free(a); free(b); free(fa); free(fb); free(fs);
	}
}

/* ==================================================================== */
/* T3  identity: a * 1 == a                                             */
/* ==================================================================== */

static void
t3_identity(void)
{
	unsigned logn;

	banner("T3  identity     a * 1 == a");
	prng_seed(3003);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *one = calloc(n, sizeof *one);
		double *r = malloc(n * sizeof *r);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		one[0] = 1;
		for (iter = 0; iter < 64; iter++) {
			rand_poly(a, n, 512);
			mine_mul(a, one, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(r[u] - (double)a[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)d_rint(r[u]) != a[u]) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "max|err| = %.3e", maxerr);
		}
		free(a); free(one); free(r);
	}
}

/* ==================================================================== */
/* T4  negacyclic shift: a * X^k == rot_k(a), sign flip on wrap         */
/* ==================================================================== */

static void
t4_negacyclic_shift(void)
{
	unsigned logn;

	banner("T4  negacyclic   a * X^k == rot_k(a), sign flip on wrap");
	prng_seed(4004);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *xk = calloc(n, sizeof *xk);
		double *r = malloc(n * sizeof *r);
		size_t ks[5], nk = 0, u, ki;
		int bad = 0;
		double maxerr = 0.0;

		ks[nk++] = 0;
		ks[nk++] = 1;
		if (n > 2) {
			ks[nk++] = 2;
			ks[nk++] = n / 2;
			ks[nk++] = n - 1;
		}

		rand_poly(a, n, 512);
		for (ki = 0; ki < nk; ki++) {
			size_t k = ks[ki];

			memset(xk, 0, n * sizeof *xk);
			xk[k] = 1;
			mine_mul(a, xk, r, logn);
			for (u = 0; u < n; u++) {
				/*
				 * (a * X^k)[u] is a[u-k] when u >= k, and
				 * -a[u-k+n] when the shift wraps.
				 */
				int64_t want;

				if (u >= k) {
					want = a[u - k];
				} else {
					want = -a[u - k + n];
				}
				double e = d_abs(r[u] - (double)want);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)d_rint(r[u]) != want) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "%zu shifts, max|err| = %.3e",
				nk, maxerr);
		}
		free(a); free(xk); free(r);
	}
}

/* ==================================================================== */
/* T5  commutativity: a * b == b * a                                    */
/* ==================================================================== */

static void
t5_commutativity(void)
{
	unsigned logn;

	banner("T5  commutative  a * b == b * a");
	prng_seed(5005);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		double *r1 = malloc(n * sizeof *r1);
		double *r2 = malloc(n * sizeof *r2);
		size_t u;
		int iter, bad = 0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			mine_mul(a, b, r1, logn);
			mine_mul(b, a, r2, logn);
			for (u = 0; u < n; u++) {
				if (r1[u] != r2[u]) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "32 random pairs");
		}
		free(a); free(b); free(r1); free(r2);
	}
}

/* ==================================================================== */
/* T6  distributivity: a * (b + c) == a*b + a*c                         */
/* ==================================================================== */

static void
t6_distributivity(void)
{
	unsigned logn;

	banner("T6  distributive a * (b+c) == a*b + a*c");
	prng_seed(6006);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		int64_t *c = malloc(n * sizeof *c);
		int64_t *bc = malloc(n * sizeof *bc);
		double *r1 = malloc(n * sizeof *r1);
		double *r2 = malloc(n * sizeof *r2);
		double *r3 = malloc(n * sizeof *r3);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 128);
			rand_poly(b, n, 128);
			rand_poly(c, n, 128);
			for (u = 0; u < n; u++) {
				bc[u] = b[u] + c[u];
			}
			mine_mul(a, bc, r1, logn);
			mine_mul(a, b, r2, logn);
			mine_mul(a, c, r3, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(r1[u] - (r2[u] + r3[u]));

				if (e > maxerr) {
					maxerr = e;
				}
				if (d_rint(r1[u]) != d_rint(r2[u] + r3[u])) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "max|err| = %.3e", maxerr);
		}
		free(a); free(b); free(c); free(bc);
		free(r1); free(r2); free(r3);
	}
}

/* ==================================================================== */
/* T7  FFT layout known-answer test                                     */
/* ==================================================================== */

/*
 * For a = X^k the evaluation at the root held in slot j is exactly
 * w^k with w = exp(i*(2*rev(j)+1)*pi/n), rev() over logn bits.  This is
 * the same closed form test_fftmul.c applies to the reference, so
 * passing it means the two implementations share a slot layout, not
 * merely a correct multiset of evaluations.
 */
static void
t7_layout_kat(void)
{
	unsigned logn;

	banner("T7  FFT layout   FFT(X^k) vs closed form w_j^k");

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		size_t hn = n >> 1;
		double *f = malloc(n * sizeof *f);
		size_t ks[4], nk = 0, u, j, ki;
		int bad = 0;
		double maxerr = 0.0;

		ks[nk++] = 0;
		ks[nk++] = 1;
		if (n > 2) {
			ks[nk++] = 2;
			ks[nk++] = n - 1;
		}

		for (ki = 0; ki < nk; ki++) {
			size_t k = ks[ki];

			for (u = 0; u < n; u++) {
				f[u] = 0.0;
			}
			f[k] = 1.0;
			my_FFT(f, logn);

			for (j = 0; j < hn; j++) {
				size_t root = bitrev(j, logn);
				double ang = (double)(2 * root + 1) * (double)k
					* M_PI / (double)n;
				double e_re = d_abs(f[j] - cos(ang));
				double e_im = d_abs(f[j + hn] - sin(ang));

				if (e_re > maxerr) {
					maxerr = e_re;
				}
				if (e_im > maxerr) {
					maxerr = e_im;
				}
				if (e_re > 1e-11 || e_im > 1e-11) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "%zu monomials, max|err| = %.3e",
				nk, maxerr);
		}
		free(f);
	}
}

/* ==================================================================== */
/* T8  differential test against the exact integer schoolbook           */
/* ==================================================================== */

static void
t8_vs_schoolbook(int iters)
{
	unsigned logn;

	banner("T8  differential my FFT multiply vs exact integer schoolbook");
	prng_seed(8008);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		__int128 *c = malloc(n * sizeof *c);
		double *r = malloc(n * sizeof *r);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;
		int num = iters * (int)(MAX_N / n);

		for (iter = 0; iter < num; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			schoolbook_mul(a, b, c, n);
			mine_mul(a, b, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(r[u] - i128_to_double(c[u]));

				if (e > maxerr) {
					maxerr = e;
				}
				if (d_rint(r[u]) != c[u]) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "%d iters, max|err| = %.3e",
				num, maxerr);
		}
		free(a); free(b); free(c); free(r);
	}
}

/* ==================================================================== */
/* T9  exactness boundary (measurement, not pass/fail)                  */
/* ==================================================================== */

static void
t9_exactness_boundary(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	int64_t *a = malloc(n * sizeof *a);
	int64_t *b = malloc(n * sizeof *b);
	__int128 *c = malloc(n * sizeof *c);
	double *r = malloc(n * sizeof *r);
	int64_t B;

	banner("T9  exactness boundary (measurement, not pass/fail)");
	printf("    logn=%u  n=%zu\n\n", logn, n);
	printf("    %12s  %14s  %14s  %9s  %8s\n",
		"coeff bound", "max|product|", "max|fft err|",
		"round ok", "repr ok");
	prng_seed(9009);

	for (B = 1; B <= 100000000LL; B *= 10) {
		double maxerr = 0.0;
		__int128 maxprod = 0;
		int round_ok = 1, repr_ok = 1;
		size_t u;

		rand_poly(a, n, B);
		rand_poly(b, n, B);
		schoolbook_mul(a, b, c, n);
		mine_mul(a, b, r, logn);
		for (u = 0; u < n; u++) {
			double e = d_abs(r[u] - i128_to_double(c[u]));

			if (e > maxerr) {
				maxerr = e;
			}
			if (i_abs(c[u]) > maxprod) {
				maxprod = i_abs(c[u]);
			}
			if (d_rint(r[u]) != c[u]) {
				round_ok = 0;
			}
		}
		/* Is the exact product itself representable in binary64? */
		if (i128_to_double(maxprod) >= 9007199254740992.0) {
			repr_ok = 0;
		}
		printf("    %12" PRId64 "  %14.4g  %14.4g  %9s  %8s\n",
			B, i128_to_double(maxprod), maxerr,
			round_ok ? "yes" : "NO", repr_ok ? "yes" : "NO");
	}
	printf("\n    Exact integer recovery requires max|fft err| < 0.5.\n");
	free(a); free(b); free(c); free(r);
}

/* ==================================================================== */
/* T10 Falcon-512 realistic operand magnitudes                          */
/* ==================================================================== */

static void
t10_falcon_operands(void)
{
	static const struct {
		const char *name;
		int64_t ba, bb;
	} cases[] = {
		{ "f * g      (small * small)",      25,    25 },
		{ "f * F      (small * large)",      25, 12289 },
		{ "F * G      (large * large)",   12289, 12289 },
	};
	unsigned logn = 9;
	size_t n = (size_t)1 << logn;
	size_t ci;

	banner("T10 Falcon-512  realistic operand magnitudes (logn=9, n=512)");
	prng_seed(10010);

	for (ci = 0; ci < sizeof cases / sizeof cases[0]; ci++) {
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		__int128 *c = malloc(n * sizeof *c);
		double *r = malloc(n * sizeof *r);
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 64; iter++) {
			size_t u;

			rand_poly(a, n, cases[ci].ba);
			rand_poly(b, n, cases[ci].bb);
			schoolbook_mul(a, b, c, n);
			mine_mul(a, b, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(r[u] - i128_to_double(c[u]));

				if (e > maxerr) {
					maxerr = e;
				}
				if (d_rint(r[u]) != c[u]) {
					bad++;
				}
			}
		}
		check(cases[ci].name, bad == 0,
			"max|err| = %.3e  (margin to 0.5: %.3gx)",
			maxerr, maxerr > 0.0 ? 0.5 / maxerr : INFINITY);
		free(a); free(b); free(c); free(r);
	}
}

/* ==================================================================== */
/* T11 cross-check against the official reference implementation        */
/* ==================================================================== */

/*
 * The project objective that this group exists to satisfy: agreement
 * with the official Falcon reference, not merely with an oracle.
 *
 * Three separate comparisons, each strictly stronger than the last as a
 * statement about interoperability:
 *
 *   forward   every FFT-domain slot agrees after the forward transform
 *   inverse   every coefficient agrees after a full round trip
 *   product   every coefficient of the finished product agrees
 *
 * The forward check is the one that constrains slot ordering.  An
 * implementation with a permuted layout would still pass the product
 * check, so the forward comparison is what makes this meaningful.
 *
 * The tolerance is absolute rather than relative because the quantities
 * being compared are evaluations of a polynomial with coefficients up
 * to 512 in magnitude at points on the unit circle, so they are O(n*B)
 * in size and an absolute bound is the honest measure.
 */
static void
t11_vs_reference(void)
{
	unsigned logn;

	banner("T11 cross-check  my implementation vs official Falcon reference");
	prng_seed(11011);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		double *mine = malloc(n * sizeof *mine);
		double *mine2 = malloc(n * sizeof *mine2);
		fpr *ref = malloc(n * sizeof *ref);
		fpr *ref2 = malloc(n * sizeof *ref2);
		size_t u;
		int iter;
		int bad_fwd = 0, bad_inv = 0, bad_prod = 0;
		double err_fwd = 0.0, err_inv = 0.0, err_prod = 0.0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);

			/* --- forward transform, slot by slot --- */
			for (u = 0; u < n; u++) {
				mine[u] = (double)a[u];
				ref[u] = fpr_of(a[u]);
			}
			my_FFT(mine, logn);
			Zf(FFT)(ref, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(mine[u] - fpr_to_double(ref[u]));

				if (e > err_fwd) {
					err_fwd = e;
				}
				if (e > 1e-9) {
					bad_fwd++;
				}
			}

			/* --- inverse, completing the round trip --- */
			my_iFFT(mine, logn);
			Zf(iFFT)(ref, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(mine[u] - fpr_to_double(ref[u]));

				if (e > err_inv) {
					err_inv = e;
				}
				if (e > 1e-9) {
					bad_inv++;
				}
			}

			/* --- the finished product --- */
			for (u = 0; u < n; u++) {
				mine[u] = (double)a[u];
				mine2[u] = (double)b[u];
				ref[u] = fpr_of(a[u]);
				ref2[u] = fpr_of(b[u]);
			}
			my_FFT(mine, logn);
			my_FFT(mine2, logn);
			my_poly_mul_fft(mine, mine2, logn);
			my_iFFT(mine, logn);

			Zf(FFT)(ref, logn);
			Zf(FFT)(ref2, logn);
			Zf(poly_mul_fft)(ref, ref2, logn);
			Zf(iFFT)(ref, logn);

			for (u = 0; u < n; u++) {
				double e = d_abs(mine[u] - fpr_to_double(ref[u]));

				if (e > err_prod) {
					err_prod = e;
				}
				if (e > 1e-6) {
					bad_prod++;
				}
			}
		}
		{
			char nm[64];

			snprintf(nm, sizeof nm, "logn=%u fwd  (slot layout)",
				logn);
			check(nm, bad_fwd == 0, "max|diff| = %.3e", err_fwd);

			snprintf(nm, sizeof nm, "logn=%u inv  (round trip)",
				logn);
			check(nm, bad_inv == 0, "max|diff| = %.3e", err_inv);

			snprintf(nm, sizeof nm, "logn=%u prod (full multiply)",
				logn);
			check(nm, bad_prod == 0, "max|diff| = %.3e", err_prod);
		}
		free(a); free(b); free(mine); free(mine2); free(ref); free(ref2);
	}
}

/* ==================================================================== */

int
main(int argc, char **argv)
{
	int iters = (argc > 1) ? atoi(argv[1]) : 8;

	printf("From-scratch FFT polynomial multiplication -- verification\n");
	printf("  ring      Z[X]/(X^N + 1)\n");
	printf("  under test  myfft.c (plain double, own twiddle tables)\n");
	printf("  oracle      exact __int128 schoolbook mod X^N+1\n");
	printf("  reference   %s\n",
		FALCON_FPNATIVE ? "Falcon fft.c, FALCON_FPNATIVE"
		: "Falcon fft.c, FALCON_FPEMU");

	t1_roundtrip();
	t2_linearity();
	t3_identity();
	t4_negacyclic_shift();
	t5_commutativity();
	t6_distributivity();
	t7_layout_kat();
	t8_vs_schoolbook(iters);
	t10_falcon_operands();
	t11_vs_reference();
	t9_exactness_boundary(9);

	printf("\n==============================================\n");
	printf("  %d checks, %d failures\n", total_checks, total_failures);
	printf("==============================================\n");

	my_fft_free();
	return total_failures == 0 ? 0 : 1;
}
