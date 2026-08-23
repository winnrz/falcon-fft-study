/*
 * myfft.c
 *
 * From-scratch FFT-based polynomial multiplication in Z[X] / (X^N + 1).
 *
 * ====================================================================
 * Derivation
 * ====================================================================
 *
 * We want the product of two integer polynomials modulo X^N + 1.  The
 * quotient ring is not the cyclic one: X^N = -1 rather than +1, so a
 * plain length-N DFT does not apply directly.  The standard route is to
 * evaluate at the N complex roots of X^N + 1,
 *
 *     w_j = exp(i * (2j+1) * pi / N),      j = 0 .. N-1
 *
 * multiply the evaluations pointwise, and interpolate back.  Two
 * reductions make this cheap.
 *
 * (1) Conjugate symmetry.  The input polynomials have real (indeed
 *     integer) coefficients, so A(conj(x)) = conj(A(x)).  The roots come
 *     in conjugate pairs, so half the evaluations are redundant.  We
 *     keep the roots with even j, that is
 *
 *         x_t = exp(i * (4t+1) * pi / N),   t = 0 .. N/2 - 1
 *
 *     which contains exactly one member of each conjugate pair.  (The
 *     conjugate of exp(i(4t+1)pi/N) is exp(i(2N-4t-1)pi/N), and
 *     2N-4t-1 = 2(N-2t-1)+1 with N-2t-1 odd, so the discarded odd-j
 *     roots are precisely the conjugates of the kept ones.)
 *
 * (2) Twisting to a cyclic transform.  Write hn = N/2 and note
 *
 *         x_t = zeta * omega^t,   zeta  = exp(i * pi / N)
 *                                 omega = exp(2 * pi * i / hn)
 *
 *     since zeta * omega^t = exp(i*pi/N) * exp(i*4*pi*t/N)
 *                          = exp(i*(4t+1)*pi/N).
 *
 *     Then, splitting the sum over u into u = v and u = v + hn,
 *
 *         A(x_t) = sum_{u=0}^{N-1} a[u] * zeta^u * omega^(t*u)
 *
 *     omega^(t*u) has period hn in u, and zeta^hn = exp(i*pi/2) = i, so
 *
 *         A(x_t) = sum_{v=0}^{hn-1} zeta^v * (a[v] + i*a[v+hn])
 *                                          * omega^(t*v)
 *
 *     which is an ordinary length-hn cyclic DFT of the twisted sequence
 *
 *         c[v] = zeta^v * (a[v] + i * a[v+hn]).
 *
 * So: twist the N real coefficients into hn complex values, run one
 * complex DFT of length hn, and we have all the evaluations we need.
 * The inverse runs the same steps backwards.
 *
 * ====================================================================
 * Output ordering
 * ====================================================================
 *
 * The DFT is done by a decimation-in-frequency radix-2 Cooley-Tukey
 * loop, and the final bit-reversal permutation is deliberately NOT
 * applied.  DIF with natural-order input leaves the output in
 * bit-reversed order, so
 *
 *     slot j  holds  A(x_t)   with  t = rev(j) over (logn - 1) bits
 *
 * This is exactly the layout the Falcon reference implementation
 * produces -- see T7 in test_fftmul.c, which checks the reference
 * against the closed form and finds slot j holding the root with index
 * rev(j) over logn bits, that index always being even, i.e. 2t with
 * t = rev(j) over logn-1 bits.
 *
 * Matching the layout is not required for the product to come out
 * right: pointwise multiplication is invariant under any permutation of
 * the slots that the inverse transform undoes consistently.  It is
 * required for comparing intermediate values against the reference,
 * which is what the cross-check in test_myfft.c does.
 *
 * ====================================================================
 * Storage
 * ====================================================================
 *
 * A transformed polynomial occupies the same n doubles as the time
 * domain one: real parts in f[0 .. hn-1], imaginary parts in
 * f[hn .. n-1].  Both transforms are fully in place and allocate
 * nothing.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "myfft.h"

#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif

/*
 * Disable contraction of floating-point expressions.
 *
 * Without this the compiler is free to fuse a multiply and an add into
 * a single FMA opcode, which rounds once instead of twice.  That is
 * usually more accurate, but it is not symmetric: the complex product
 *
 *     Im = ar*bi + ai*br
 *
 * contracts to fma(ar, bi, ai*br), keeping ar*bi exact and rounding
 * ai*br.  Exchanging the operands keeps the other product exact
 * instead, so a*b and b*a can differ in the last bit and bit-exact
 * commutativity is lost.  Measured here on arm64: T5 fails at every
 * logn > 1 with contraction left on.
 *
 * The Falcon reference disables contraction for the same reason; see
 * inner.h, which notes that GCC has historically ignored the C99
 * pragma and needs its own.  This mirrors that treatment so the two
 * implementations are compiled under the same floating-point rules.
 */
#if defined __clang__
#pragma STDC FP_CONTRACT OFF
#elif defined __GNUC__
#pragma GCC optimize ("fp-contract=off")
#endif

/* ==================================================================== */
/* Twiddle tables                                                       */
/* ==================================================================== */

/*
 * Per-logn tables, built on demand and kept for the life of the
 * process.  Two tables are needed:
 *
 *   tw[v]  = zeta^v = exp(i * pi * v / n),   v = 0 .. hn-1
 *            the twist that converts negacyclic to cyclic
 *
 *   w[m]   = omega^m = exp(2 * pi * i * m / hn),  m = 0 .. hn/2 - 1
 *            the DFT twiddles, shared across every stage by striding
 *
 * Real and imaginary parts are kept in separate arrays: it keeps the
 * butterfly inner loop free of struct member access, and it mirrors how
 * the transformed data itself is stored.
 */
struct twiddle_tab {
	double *tw_re, *tw_im;   /* hn entries    */
	double *w_re,  *w_im;    /* hn/2 entries  */
	size_t bytes;            /* heap held by the four arrays above */
	int ready;
};

static struct twiddle_tab tabs[MYFFT_MAX_LOGN + 1];

int
my_fft_init(unsigned logn)
{
	size_t n, hn, hw, v, m;
	struct twiddle_tab *t;

	if (logn < 1 || logn > MYFFT_MAX_LOGN) {
		return 0;
	}
	t = &tabs[logn];
	if (t->ready) {
		return 1;
	}

	n = (size_t)1 << logn;
	hn = n >> 1;
	hw = hn >> 1;          /* zero when logn == 1; nothing to stride */

	t->tw_re = malloc(hn * sizeof *t->tw_re);
	t->tw_im = malloc(hn * sizeof *t->tw_im);
	t->w_re = malloc((hw ? hw : 1) * sizeof *t->w_re);
	t->w_im = malloc((hw ? hw : 1) * sizeof *t->w_im);
	if (t->tw_re == NULL || t->tw_im == NULL
		|| t->w_re == NULL || t->w_im == NULL)
	{
		free(t->tw_re); free(t->tw_im);
		free(t->w_re); free(t->w_im);
		memset(t, 0, sizeof *t);
		return 0;
	}
	t->bytes = (2 * hn + 2 * (hw ? hw : 1)) * sizeof(double);

	/*
	 * The twist.  Computing each angle from scratch rather than by
	 * repeated multiplication keeps the table accurate to within one
	 * ulp of the true value; iterating zeta^v = zeta^(v-1) * zeta
	 * would accumulate error across hn steps.
	 */
	for (v = 0; v < hn; v++) {
		double ang = M_PI * (double)v / (double)n;
		t->tw_re[v] = cos(ang);
		t->tw_im[v] = sin(ang);
	}

	/* The DFT twiddles, positive exponent to match omega above. */
	for (m = 0; m < hw; m++) {
		double ang = 2.0 * M_PI * (double)m / (double)hn;
		t->w_re[m] = cos(ang);
		t->w_im[m] = sin(ang);
	}

	t->ready = 1;
	return 1;
}

void
my_fft_free(void)
{
	unsigned logn;

	for (logn = 0; logn <= MYFFT_MAX_LOGN; logn++) {
		struct twiddle_tab *t = &tabs[logn];

		if (t->ready) {
			free(t->tw_re); free(t->tw_im);
			free(t->w_re); free(t->w_im);
			memset(t, 0, sizeof *t);
		}
	}
}

size_t
my_fft_table_bytes(void)
{
	size_t total = 0;
	unsigned logn;

	for (logn = 0; logn <= MYFFT_MAX_LOGN; logn++) {
		if (tabs[logn].ready) {
			total += tabs[logn].bytes;
		}
	}
	return total;
}

/* ==================================================================== */
/* Forward transform                                                    */
/* ==================================================================== */

void
my_FFT(double *f, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	const struct twiddle_tab *t;
	size_t len, v;

	if (!my_fft_init(logn)) {
		return;
	}
	t = &tabs[logn];

	/*
	 * Step 1: twist.  f currently holds n real coefficients.  Fold
	 * them into hn complex values c[v] = zeta^v * (a[v] + i*a[v+hn]),
	 * stored as real part in f[v] and imaginary part in f[v+hn].
	 */
	for (v = 0; v < hn; v++) {
		double re = f[v];
		double im = f[v + hn];
		double zr = t->tw_re[v];
		double zi = t->tw_im[v];

		f[v]      = re * zr - im * zi;
		f[v + hn] = re * zi + im * zr;
	}

	/*
	 * Step 2: decimation-in-frequency radix-2 DFT of length hn.
	 *
	 * Each stage halves the transform length.  For a stage of length
	 * `len` the butterfly on the pair (p, p+half) is
	 *
	 *     A' = A + B
	 *     B' = (A - B) * omega_len^k
	 *
	 * and omega_len^k = w[k * (hn/len)], which is why one table of
	 * hn/2 entries serves every stage.
	 *
	 * No bit-reversal pass follows; the bit-reversed output order is
	 * the layout we want.
	 */
	for (len = hn; len >= 2; len >>= 1) {
		size_t half = len >> 1;
		size_t step = hn / len;
		size_t start;

		for (start = 0; start < hn; start += len) {
			size_t k;

			for (k = 0; k < half; k++) {
				size_t p = start + k;
				size_t q = p + half;
				double ar = f[p],      ai = f[p + hn];
				double br = f[q],      bi = f[q + hn];
				double dr = ar - br,   di = ai - bi;
				double wr = t->w_re[k * step];
				double wi = t->w_im[k * step];

				f[p]      = ar + br;
				f[p + hn] = ai + bi;
				f[q]      = dr * wr - di * wi;
				f[q + hn] = dr * wi + di * wr;
			}
		}
	}
}

/* ==================================================================== */
/* Inverse transform                                                    */
/* ==================================================================== */

void
my_iFFT(double *f, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	const struct twiddle_tab *t;
	size_t len, v;
	double scale;

	if (!my_fft_init(logn)) {
		return;
	}
	t = &tabs[logn];

	/*
	 * Step 1: undo the DFT, decimation in time, taking bit-reversed
	 * input back to natural order.  This is the forward butterfly
	 * run backwards:  from A' = A + B and B' = (A-B)*w we recover
	 * A = (A' + B'/w) / 2 and B = (A' - B'/w) / 2.  The halvings are
	 * deferred and applied once at the end as a single 1/hn.
	 *
	 * 1/w is conj(w) because every twiddle is on the unit circle.
	 */
	for (len = 2; len <= hn; len <<= 1) {
		size_t half = len >> 1;
		size_t step = hn / len;
		size_t start;

		for (start = 0; start < hn; start += len) {
			size_t k;

			for (k = 0; k < half; k++) {
				size_t p = start + k;
				size_t q = p + half;
				double ar = f[p],  ai = f[p + hn];
				double br = f[q],  bi = f[q + hn];
				double wr = t->w_re[k * step];
				double wi = -t->w_im[k * step];   /* conj */
				double cr = br * wr - bi * wi;
				double ci = br * wi + bi * wr;

				f[p]      = ar + cr;
				f[p + hn] = ai + ci;
				f[q]      = ar - cr;
				f[q + hn] = ai - ci;
			}
		}
	}

	/*
	 * Step 2: untwist and scale.  c[v] / zeta^v = a[v] + i*a[v+hn],
	 * so the real part is the low coefficient and the imaginary part
	 * the high one.  The 1/hn from the deferred butterfly halvings is
	 * folded in here rather than in a separate pass.
	 */
	scale = 1.0 / (double)hn;
	for (v = 0; v < hn; v++) {
		double re = f[v];
		double im = f[v + hn];
		double zr = t->tw_re[v];
		double zi = -t->tw_im[v];        /* conj(zeta^v) */

		f[v]      = (re * zr - im * zi) * scale;
		f[v + hn] = (re * zi + im * zr) * scale;
	}
}

/* ==================================================================== */
/* Operations in the FFT domain                                         */
/* ==================================================================== */

void
my_poly_mul_fft(double *a, const double *b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t j;

	for (j = 0; j < hn; j++) {
		double ar = a[j], ai = a[j + hn];
		double br = b[j], bi = b[j + hn];

		a[j]      = ar * br - ai * bi;
		a[j + hn] = ar * bi + ai * br;
	}
}

void
my_poly_add(double *a, const double *b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t u;

	for (u = 0; u < n; u++) {
		a[u] += b[u];
	}
}

/* ==================================================================== */
/* The complete integer multiply                                        */
/* ==================================================================== */

size_t
my_poly_mul_scratch_bytes(unsigned logn)
{
	return ((size_t)1 << logn) * sizeof(double);
}

void
my_poly_mul(const int64_t *a, const int64_t *b, double *out, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	double *tmp;
	size_t u;

	tmp = malloc(n * sizeof *tmp);
	if (tmp == NULL) {
		return;
	}
	for (u = 0; u < n; u++) {
		out[u] = (double)a[u];
		tmp[u] = (double)b[u];
	}
	my_FFT(out, logn);
	my_FFT(tmp, logn);
	my_poly_mul_fft(out, tmp, logn);
	my_iFFT(out, logn);
	free(tmp);
}
