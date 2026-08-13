/*
 * test_fftmul.c
 *
 * Verification harness for the Falcon-512 FFT-based polynomial
 * multiplication in the ring Z[X] / (X^N + 1).
 *
 * The tests are layered, weakest-but-cheapest first:
 *
 *   T1  round trip            iFFT(FFT(a)) == a
 *   T2  linearity             FFT(a + b) == FFT(a) + FFT(b)
 *   T3  identity              a * 1 == a
 *   T4  negacyclic shift      a * X^k == rot_k(a) with sign flip on wrap
 *   T5  commutativity         a * b == b * a
 *   T6  distributivity        a * (b + c) == a*b + a*c
 *   T7  FFT layout (KAT)      FFT(X^k) against the closed form w_j^k
 *   T8  exact differential    FFT multiply vs exact integer schoolbook
 *   T9  exactness boundary    where rounding stops recovering the integer
 *   T10 Falcon-512 operands   realistic f, g, F, G coefficient ranges
 *
 * T8 is the primary correctness oracle: an exact __int128 schoolbook
 * multiply mod X^N+1 that shares no code and no floating point with the
 * FFT.  It catches essentially every structural error -- wrong twiddles,
 * wrong final scaling, cyclic-vs-negacyclic sign confusion, misapplied
 * stage elision.
 *
 * T7 exists because T8 has one blind spot.  Pointwise multiplication is
 * invariant under permutation of the FFT slots, so an implementation that
 * emits the correct multiset of evaluations in the wrong order, and whose
 * iFFT undoes that same permutation, passes T1 and T8 while being
 * incompatible with the reference layout.  T7 pins the layout down by
 * checking slot values against a closed form.
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

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define MAX_LOGN   10
#define MAX_N      (1u << MAX_LOGN)

/* ==================================================================== */
/* Small utilities                                                      */
/* ==================================================================== */

/*
 * Extract the IEEE-754 binary64 value behind an fpr.  Both backends store
 * standard binary64 bit patterns: FALCON_FPNATIVE wraps a double in a
 * struct, FALCON_FPEMU keeps the raw bits in a uint64_t.
 */
static double
to_double(fpr x)
{
#if FALCON_FPNATIVE
	return x.v;
#else
	double d;
	uint64_t u = (uint64_t)x;
	memcpy(&d, &u, sizeof d);
	return d;
#endif
}

/* Deterministic PRNG (xorshift64*), so every run is reproducible. */
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

/* Uniform integer in [-B, B]. */
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

/* Bit reversal of x over the low `bits` bits. */
static size_t
bitrev(size_t x, unsigned bits)
{
	size_t r = 0;
	unsigned i;

	for (i = 0; i < bits; i++) {
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

/* __int128 has no printf conversion; render via double for display. */
static double
i128_to_double(__int128 x)
{
	return (double)x;
}

/* ==================================================================== */
/* Reference oracle: exact integer schoolbook multiply mod X^n + 1      */
/* ==================================================================== */

/*
 * c[k] = sum_{i+j=k} a[i]b[j]  -  sum_{i+j=k+n} a[i]b[j]
 *
 * Computed in __int128 so it stays exact well past the point where the
 * true product coefficients cease to be representable in binary64.  No
 * floating point is involved anywhere in this function.
 */
static void
schoolbook_mul(const int64_t *a, const int64_t *b, __int128 *c, size_t n)
{
	size_t i, j, k;

	for (k = 0; k < n; k++) {
		c[k] = 0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			k = i + j;
			if (k < n) {
				c[k] += (__int128)a[i] * (__int128)b[j];
			} else {
				c[k - n] -= (__int128)a[i] * (__int128)b[j];
			}
		}
	}
}

/*
 * The operation under test: multiply two integer polynomials through the
 * FFT domain.  Result is left as fpr (unrounded) so callers can inspect
 * the floating-point error before rounding.
 */
static void
fft_mul(const int64_t *a, const int64_t *b, fpr *out, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	fpr *fb = malloc(n * sizeof *fb);
	size_t u;

	for (u = 0; u < n; u++) {
		out[u] = fpr_of(a[u]);
		fb[u] = fpr_of(b[u]);
	}
	Zf(FFT)(out, logn);
	Zf(FFT)(fb, logn);
	Zf(poly_mul_fft)(out, fb, logn);
	Zf(iFFT)(out, logn);

	free(fb);
}

/* ==================================================================== */
/* Test bookkeeping                                                     */
/* ==================================================================== */

static int total_checks = 0;
static int total_failures = 0;

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
		fpr *f = malloc(n * sizeof *f);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 64; iter++) {
			rand_poly(a, n, 512);
			for (u = 0; u < n; u++) {
				f[u] = fpr_of(a[u]);
			}
			Zf(FFT)(f, logn);
			Zf(iFFT)(f, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(to_double(f[u]) - (double)a[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)fpr_rint(f[u]) != a[u]) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "max|err| = %.3e", maxerr);
		}
		free(a);
		free(f);
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
		fpr *fa = malloc(n * sizeof *fa);
		fpr *fb = malloc(n * sizeof *fb);
		fpr *fs = malloc(n * sizeof *fs);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			for (u = 0; u < n; u++) {
				fa[u] = fpr_of(a[u]);
				fb[u] = fpr_of(b[u]);
				fs[u] = fpr_of(a[u] + b[u]);
			}
			Zf(FFT)(fa, logn);
			Zf(FFT)(fb, logn);
			Zf(FFT)(fs, logn);
			for (u = 0; u < n; u++) {
				double lhs = to_double(fs[u]);
				double rhs = to_double(fa[u]) + to_double(fb[u]);
				double e = d_abs(lhs - rhs);
				double tol = 1e-9 + 1e-12 * d_abs(rhs);

				if (e > maxerr) {
					maxerr = e;
				}
				if (e > tol) {
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
		fpr *r = malloc(n * sizeof *r);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;

		one[0] = 1;
		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 512);
			fft_mul(a, one, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(to_double(r[u]) - (double)a[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)fpr_rint(r[u]) != a[u]) {
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
/* T4  negacyclic shift: a * X^k                                        */
/* ==================================================================== */

/*
 * In Z[X]/(X^N+1), multiplying by X^k rotates coefficients left by k and
 * negates everything that wraps past the top.  This is the sharpest cheap
 * test of the X^N = -1 convention: an implementation that accidentally
 * computes a cyclic (X^N - 1) convolution passes T1, T2, T3 and fails here.
 */
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
		int64_t *want = malloc(n * sizeof *want);
		fpr *r = malloc(n * sizeof *r);
		size_t ks[5], nk = 0, u, i, ki;
		int bad = 0;
		double maxerr = 0.0;

		ks[nk++] = 0;
		ks[nk++] = 1;
		if (n > 2) {
			ks[nk++] = 2;
			ks[nk++] = n >> 1;
		}
		ks[nk++] = n - 1;

		rand_poly(a, n, 512);
		for (ki = 0; ki < nk; ki++) {
			size_t k = ks[ki];

			memset(xk, 0, n * sizeof *xk);
			xk[k] = 1;

			for (u = 0; u < n; u++) {
				want[u] = 0;
			}
			for (i = 0; i < n; i++) {
				size_t j = i + k;

				if (j < n) {
					want[j] += a[i];
				} else {
					want[j - n] -= a[i];
				}
			}

			fft_mul(a, xk, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(to_double(r[u]) - (double)want[u]);

				if (e > maxerr) {
					maxerr = e;
				}
				if ((int64_t)fpr_rint(r[u]) != want[u]) {
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
		free(a); free(xk); free(want); free(r);
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
		fpr *ab = malloc(n * sizeof *ab);
		fpr *ba = malloc(n * sizeof *ba);
		size_t u;
		int iter, bad = 0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			fft_mul(a, b, ab, logn);
			fft_mul(b, a, ba, logn);
			for (u = 0; u < n; u++) {
				if ((int64_t)fpr_rint(ab[u])
					!= (int64_t)fpr_rint(ba[u]))
				{
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "32 random pairs");
		}
		free(a); free(b); free(ab); free(ba);
	}
}

/* ==================================================================== */
/* T6  distributivity: a * (b + c) == a*b + a*c                         */
/* ==================================================================== */

static void
t6_distributivity(void)
{
	unsigned logn;

	banner("T6  distributive a*(b+c) == a*b + a*c");
	prng_seed(6006);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		int64_t *c = malloc(n * sizeof *c);
		int64_t *bc = malloc(n * sizeof *bc);
		fpr *r1 = malloc(n * sizeof *r1);
		fpr *r2 = malloc(n * sizeof *r2);
		fpr *r3 = malloc(n * sizeof *r3);
		size_t u;
		int iter, bad = 0;

		for (iter = 0; iter < 32; iter++) {
			rand_poly(a, n, 128);
			rand_poly(b, n, 128);
			rand_poly(c, n, 128);
			for (u = 0; u < n; u++) {
				bc[u] = b[u] + c[u];
			}
			fft_mul(a, bc, r1, logn);
			fft_mul(a, b, r2, logn);
			fft_mul(a, c, r3, logn);
			for (u = 0; u < n; u++) {
				int64_t lhs = (int64_t)fpr_rint(r1[u]);
				int64_t rhs = (int64_t)fpr_rint(r2[u])
					+ (int64_t)fpr_rint(r3[u]);

				if (lhs != rhs) {
					bad++;
				}
			}
		}
		{
			char nm[64];
			snprintf(nm, sizeof nm, "logn=%u (n=%zu)", logn, n);
			check(nm, bad == 0, "32 random triples");
		}
		free(a); free(b); free(c); free(bc);
		free(r1); free(r2); free(r3);
	}
}

/* ==================================================================== */
/* T7  FFT layout, analytic known-answer test                           */
/* ==================================================================== */

/*
 * FFT representation stores evaluations at the roots of X^N+1,
 * w_j = exp(i*(2j+1)*pi/N), with real and imaginary parts split across
 * the array.  The slot-to-root mapping, determined empirically and
 * checked here for every logn, is
 *
 *     slot s        holds  Re(f(w_j))     where j = rev(s)
 *     slot s + N/2  holds  Im(f(w_j))
 *
 * with rev() the bit reversal over logn bits applied to the SLOT index.
 * Only even j occur, which is how the conjugate-pair halving shows up:
 * the N/2 stored slots cover w_0, w_2, w_4, ... in bit-reversed order.
 *
 * (Note: the prose in fft.c's header comment states this as
 * "Re(f(w_j)) -> slot rev(j)/2", which does not reproduce the observed
 * layout.  The form above is what the code actually implements.)
 *
 * For a = X^k the evaluation is exactly w_j^k = exp(i*(2j+1)*k*pi/N), so
 * every slot has a closed form.  This is the only test here that
 * constrains the ordering and the real/imaginary split; T1 and T8 are
 * both blind to a consistent permutation of the slots.
 */
static void
t7_layout_kat(void)
{
	unsigned logn;

	banner("T7  FFT layout   FFT(X^k) vs closed form w_j^k");

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		size_t hn = n >> 1;
		int64_t *a = calloc(n, sizeof *a);
		fpr *f = malloc(n * sizeof *f);
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

			memset(a, 0, n * sizeof *a);
			a[k] = 1;
			for (u = 0; u < n; u++) {
				f[u] = fpr_of(a[u]);
			}
			Zf(FFT)(f, logn);

			for (j = 0; j < hn; j++) {
				/* slot j holds the evaluation at w_{rev(j)} */
				size_t root = bitrev(j, logn);
				double ang = (double)(2 * root + 1) * (double)k
					* M_PI / (double)n;
				double want_re = cos(ang);
				double want_im = sin(ang);
				double got_re = to_double(f[j]);
				double got_im = to_double(f[j + hn]);
				double e_re = d_abs(got_re - want_re);
				double e_im = d_abs(got_im - want_im);

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
		free(a); free(f);
	}
}

/* ==================================================================== */
/* T8  exact differential test against integer schoolbook               */
/* ==================================================================== */

static void
t8_vs_schoolbook(int iters)
{
	unsigned logn;

	banner("T8  differential FFT multiply vs exact integer schoolbook");
	prng_seed(8008);

	for (logn = 1; logn <= MAX_LOGN; logn++) {
		size_t n = (size_t)1 << logn;
		int64_t *a = malloc(n * sizeof *a);
		int64_t *b = malloc(n * sizeof *b);
		__int128 *c = malloc(n * sizeof *c);
		fpr *r = malloc(n * sizeof *r);
		size_t u;
		int iter, bad = 0;
		double maxerr = 0.0;
		/* Scale iteration count so small n get more coverage. */
		int num = iters * (int)(MAX_N / n);

		for (iter = 0; iter < num; iter++) {
			rand_poly(a, n, 512);
			rand_poly(b, n, 512);
			schoolbook_mul(a, b, c, n);
			fft_mul(a, b, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(to_double(r[u])
					- i128_to_double(c[u]));

				if (e > maxerr) {
					maxerr = e;
				}
				if ((__int128)fpr_rint(r[u]) != c[u]) {
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
/* T9  exactness boundary                                               */
/* ==================================================================== */

/*
 * Not a pass/fail test but a measurement: how far the input coefficients
 * can grow before rounding the FFT result stops recovering the exact
 * integer product.  Two distinct failure modes are reported separately.
 *
 *   "round ok"   the rounded FFT result equals the exact product
 *   "repr ok"    the exact product itself fits in binary64 (< 2^53)
 *
 * Once "repr ok" is no, no binary64 algorithm could succeed, so the
 * failure is no longer attributable to FFT error accumulation.
 */
static void
t9_exactness_boundary(unsigned logn)
{
	size_t n = (size_t)1 << logn;
	int64_t *a = malloc(n * sizeof *a);
	int64_t *b = malloc(n * sizeof *b);
	__int128 *c = malloc(n * sizeof *c);
	fpr *r = malloc(n * sizeof *r);
	int64_t B;

	banner("T9  exactness boundary (measurement, not pass/fail)");
	printf("    logn=%u  n=%zu\n\n", logn, n);
	printf("    %12s  %14s  %14s  %9s  %8s\n",
		"coeff bound", "max|product|", "max|fft err|",
		"round ok", "repr ok");
	prng_seed(9009);

	for (B = 1; B <= 100000000LL; B *= 10) {
		size_t u;
		double maxerr = 0.0;
		__int128 maxc = 0;
		int roundok = 1, reprok = 1;

		rand_poly(a, n, B);
		rand_poly(b, n, B);
		schoolbook_mul(a, b, c, n);
		fft_mul(a, b, r, logn);

		for (u = 0; u < n; u++) {
			__int128 m = i_abs(c[u]);
			double e;

			if (m > maxc) {
				maxc = m;
			}
			if (m > (((__int128)1) << 53)) {
				reprok = 0;
			}
			e = d_abs(to_double(r[u]) - i128_to_double(c[u]));
			if (e > maxerr) {
				maxerr = e;
			}
			if ((__int128)fpr_rint(r[u]) != c[u]) {
				roundok = 0;
			}
		}
		printf("    %12" PRId64 "  %14.4g  %14.4g  %9s  %8s\n",
			B, i128_to_double(maxc), maxerr,
			roundok ? "yes" : "NO", reprok ? "yes" : "NO");
	}
	printf("\n    Exact integer recovery requires max|fft err| < 0.5.\n");

	free(a); free(b); free(c); free(r);
}

/* ==================================================================== */
/* T10  Falcon-512 realistic operands                                   */
/* ==================================================================== */

/*
 * The operand magnitudes Falcon-512 actually produces:
 *
 *   f, g   short key polynomials, discrete Gaussian, sigma ~ 4.05
 *          -> coefficients essentially always within +/- 25
 *   F, G   NTRU completion, substantially larger, bounded on the order
 *          of q = 12289
 *
 * The products that matter in signing are f*F, g*G and similar, so the
 * worst realistic case is roughly the small-times-large and
 * large-times-large regimes below.
 */
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
		fpr *r = malloc(n * sizeof *r);
		int iter, bad = 0;
		double maxerr = 0.0;

		for (iter = 0; iter < 64; iter++) {
			size_t u;

			rand_poly(a, n, cases[ci].ba);
			rand_poly(b, n, cases[ci].bb);
			schoolbook_mul(a, b, c, n);
			fft_mul(a, b, r, logn);
			for (u = 0; u < n; u++) {
				double e = d_abs(to_double(r[u])
					- i128_to_double(c[u]));

				if (e > maxerr) {
					maxerr = e;
				}
				if ((__int128)fpr_rint(r[u]) != c[u]) {
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

int
main(int argc, char **argv)
{
	int iters = (argc > 1) ? atoi(argv[1]) : 8;

	printf("Falcon FFT polynomial multiplication -- verification harness\n");
	printf("  ring      Z[X]/(X^N + 1)\n");
	printf("  backend   %s\n",
		FALCON_FPNATIVE ? "FALCON_FPNATIVE (hardware binary64)"
		: "FALCON_FPEMU (software binary64)");
	printf("  AVX2      %s\n", FALCON_AVX2 ? "enabled" : "disabled");
	printf("  oracle    exact __int128 schoolbook mod X^N+1\n");

	t1_roundtrip();
	t2_linearity();
	t3_identity();
	t4_negacyclic_shift();
	t5_commutativity();
	t6_distributivity();
	t7_layout_kat();
	t8_vs_schoolbook(iters);
	t10_falcon_operands();
	t9_exactness_boundary(9);

	printf("\n==============================================\n");
	printf("  %d checks, %d failures\n", total_checks, total_failures);
	printf("==============================================\n");

	return total_failures == 0 ? 0 : 1;
}
