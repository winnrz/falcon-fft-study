# Falcon-512 FFT Polynomial Multiplication — Study

A minimal extraction of the Falcon reference implementation, reduced to just
the code needed to implement, verify, and analyse FFT-based polynomial
multiplication in the ring **Z[X] / (X^N + 1)** as used by Falcon-512.

## Scope

Falcon-512 signing and key generation multiply polynomials in `Z[X]/(X^N+1)`
with `N = 512`. Rather than an `O(N²)` schoolbook multiply, the reference
implementation evaluates each polynomial at the `N` complex roots of `X^N+1`
via an FFT, multiplies pointwise, and transforms back — `O(N log N)`.

This repository isolates that operation for study along three axes:
**correctness**, **numerical stability**, and **computational performance**.

## Files

Extracted verbatim from the Falcon reference implementation (round 3,
2021-11-01), MIT licensed, © 2017-2019 Falcon Project, author Thomas Pornin:

| File | Role |
|---|---|
| `fft.c` | The FFT, inverse FFT, and the polynomial ops built on them |
| `fpr.c` | Floating-point backend: constant tables, software binary64 |
| `fpr.h` | `fpr` type and arithmetic, both native and emulated |
| `inner.h` | Internal declarations |
| `config.h` | Build-time configuration macros |

Written for this study:

| File | Role |
|---|---|
| `test_fftmul.c` | Verification harness (this document's subject) |
| `Makefile` | Builds and runs the harness under both FP backends |

`fft.c` depends on `fpr.c` only for two constant tables
(`fpr_gm_tab`, `fpr_p2_tab`); `fpr.c` is fully self-contained. Nothing else
from the reference implementation is required.

## Build and run

```sh
make test         # hardware binary64 backend
make test-both    # run under both hardware and software binary64
make clean
```

Exit status is 0 if every check passes, 1 otherwise.

## What the harness verifies

Ten test groups, ordered weakest-but-cheapest first. Every group runs for
`logn = 1..10`, which exercises the degenerate small-`n` cases and the
non-vector fallback paths as well as Falcon-512's `logn = 9`.

| | Test | Property |
|---|---|---|
| T1 | round trip | `iFFT(FFT(a)) == a` |
| T2 | linearity | `FFT(a+b) == FFT(a) + FFT(b)` |
| T3 | identity | `a * 1 == a` |
| T4 | negacyclic shift | `a * X^k == rot_k(a)`, sign flip on wrap |
| T5 | commutativity | `a * b == b * a` |
| T6 | distributivity | `a * (b+c) == a*b + a*c` |
| T7 | FFT layout (KAT) | `FFT(X^k)` against closed form `w_j^k` |
| T8 | exact differential | FFT multiply vs exact integer schoolbook |
| T9 | exactness boundary | where rounding stops recovering the integer |
| T10 | Falcon-512 operands | realistic `f`, `g`, `F`, `G` magnitudes |

### The primary oracle (T8)

An exact `__int128` schoolbook multiply mod `X^N+1`. It shares no code and no
floating point with the FFT, so it is a genuinely independent reference. It
catches essentially every structural error: wrong twiddle factors, wrong final
scaling, cyclic-versus-negacyclic sign confusion, misapplied stage elision.

### The oracle's blind spot (T7)

**Pointwise multiplication is invariant under permutation of the FFT slots.**
An implementation that emits the correct *multiset* of evaluations in the
wrong order, and whose inverse undoes that same permutation, produces exactly
correct products — and therefore passes T1 and T8 while being incompatible
with the reference layout.

This matters for any re-implementation (in particular a hardware one) that
needs to interoperate with reference intermediate values, or that wants to
feed `poly_split_fft` / `poly_merge_fft`, both of which index the twiddle
table directly and so depend on the ordering.

T7 pins the layout down by checking slot values against a closed form. For
`a = X^k` the evaluation at root `w_j` is exactly `w_j^k = exp(i(2j+1)kπ/N)`.

**The layout, as the code actually implements it:**

```
slot s        holds  Re(f(w_j))    where j = rev(s)
slot s + N/2  holds  Im(f(w_j))
```

with `rev()` the bit reversal over `logn` bits applied to the **slot** index.
Only even `j` occur — that is how the conjugate-pair halving manifests: the
`N/2` stored slots cover `w_0, w_2, w_4, …` in bit-reversed order.

Note that `fft.c`'s own header comment states this as
`Re(f(w_j)) -> slot rev(j)/2`, which does **not** reproduce the observed
layout. Verified empirically for `logn = 3`, slots 0–3 hold roots with
exponents 1, 9, 5, 13, i.e. `j = 0, 4, 2, 6 = rev(s)`. The form above is
what the code does, and T7 checks it for every `logn`.

## Results

83 checks, 0 failures, under both backends. `FALCON_FPNATIVE` and
`FALCON_FPEMU` agree on every reported value — the software binary64
implementation is bit-exact against the hardware FPU here.

### Numerical stability (T9, `n = 512`)

| coeff bound | max product | max FFT error | rounds exactly | representable |
|---|---|---|---|---|
| 1 | 44 | 2.13e-14 | yes | yes |
| 10 | 2 697 | 1.36e-12 | yes | yes |
| 100 | 2.86e+05 | 1.02e-10 | yes | yes |
| 1 000 | 2.54e+07 | 7.45e-09 | yes | yes |
| 10 000 | 2.82e+09 | 9.54e-07 | yes | yes |
| 100 000 | 2.17e+11 | 1.07e-04 | yes | yes |
| 1 000 000 | 2.23e+13 | 9.77e-03 | yes | yes |
| 10 000 000 | 2.93e+15 | 1.25 | **no** | yes |
| 100 000 000 | 2.16e+17 | 96 | **no** | **no** |

Error scales as the square of the coefficient bound — each 10× in input gives
100× in error — consistent with `err ≈ C · B² · N · log₂N · ε`.

Two failure modes must be distinguished. At `B = 10⁷` the FFT accumulation
error alone exceeds the 0.5 needed for exact integer recovery. At `B = 10⁸`
the true product coefficient exceeds `2⁵³`, so no binary64 algorithm could
succeed regardless of accumulation error.

### Falcon-512's actual operating point (T10)

| operand pair | max error | margin to 0.5 |
|---|---|---|
| `f * g` (small × small) | 8.19e-12 | 6.1e+10 × |
| `f * F` (small × large) | 3.73e-09 | 1.3e+08 × |
| `F * G` (large × large) | 1.91e-06 | 2.6e+05 × |

Falcon-512 operates with at least five orders of magnitude of headroom in the
worst case, so double-precision FFT multiplication is comfortably exact for
this parameter set.

## Platform note

Built and verified on arm64 (Apple silicon), where `FALCON_AVX2` is
necessarily disabled — the scalar code paths in `fft.c` are what execute.
On x86-64 the AVX2 branches activate and should be re-verified separately.

## License

The extracted Falcon sources retain their original MIT license; see the
header of each file. Material written for this study is offered under the
same terms.
