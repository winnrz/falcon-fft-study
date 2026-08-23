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

It contains two implementations. The Falcon reference, extracted verbatim, is
the object of study and the comparison target. Alongside it sits an
independent implementation written from scratch for this project, verified
against the same oracle and cross-checked against the reference value by
value. Performance and memory are measured for both.

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
| `test_fftmul.c` | Verification harness for the reference |
| `myfft.c`, `myfft.h` | The from-scratch implementation |
| `test_myfft.c` | Verification harness for it, plus the reference cross-check |
| `bench_fftmul.c` | Performance measurement, both implementations |
| `measure_memory.sh` | Peak RSS, heap profile, leak check |
| `Makefile` | Builds and runs everything |

`fft.c` depends on `fpr.c` only for two constant tables
(`fpr_gm_tab`, `fpr_p2_tab`); `fpr.c` is fully self-contained. Nothing else
from the reference implementation is required.

## Build and run

```sh
make test         # verify the reference
make test-mine    # verify the from-scratch code, cross-check vs reference
make test-all     # both verification harnesses
make test-both    # the reference under both floating-point backends
make bench        # performance
make bench-csv    # same, as bench.csv for plotting
make memory       # peak RSS, heap profile, leak check
make clean
```

Exit status of the verification harnesses is 0 if every check passes,
1 otherwise.

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

## The from-scratch implementation

`myfft.c` implements the same operation independently. It shares no code with
`fft.c`: it uses the plain C `double` type rather than the reference's `fpr`
abstraction, and it generates its own twiddle factors from `cos`/`sin` rather
than reading `fpr_gm_tab`.

### Derivation

Evaluating at the `N` roots of `X^N+1` is not a cyclic DFT, so the transform
is obtained by twisting it into one. Writing `hn = N/2`,

```
x_t   = exp(i(4t+1)π/N) = ζ · ω^t      ζ = exp(iπ/N), ω = exp(2πi/hn)
```

and splitting the coefficient sum over `u = v` and `u = v + hn`, with
`ζ^hn = i`,

```
A(x_t) = Σ_v  ζ^v · (a[v] + i·a[v+hn]) · ω^(tv)
```

which is an ordinary length-`hn` cyclic DFT of the twisted sequence
`c[v] = ζ^v · (a[v] + i·a[v+hn])`. So: fold `N` real coefficients into `hn`
complex ones, run one complex DFT of half the length, and every evaluation
needed is in hand. Conjugate symmetry is what makes the halving legitimate —
the roots with even index contain exactly one member of each conjugate pair.

### Matching the reference layout

The DFT is a decimation-in-frequency radix-2 loop with the final bit-reversal
permutation deliberately omitted. DIF on natural-order input leaves the output
in bit-reversed order, which is precisely the reference's layout:

```
slot j  holds  A(x_t)   with  t = rev(j) over logn-1 bits
```

This falls out of the construction rather than being imposed on it, which is
the useful part: the reference's ordering is not arbitrary, it is what a
standard DIF transform produces.

### Verification

`make test-mine` runs the same ten property groups as the reference harness,
with the same seeds and coefficient ranges, and adds **T11**, a direct
cross-check against the official reference at every size:

| logn = 9 comparison | max abs difference |
|---|---|
| forward transform, slot by slot | 1.18e-11 |
| inverse, completing the round trip | 5.12e-13 |
| finished product, coefficient by coefficient | 5.59e-09 |

113 checks, 0 failures. The forward comparison is the one that constrains slot
ordering — a permuted layout would still pass the product comparison.

Accuracy is close to the reference but not identical, being consistently
around 1.5× worse:

| operand pair | reference | from scratch |
|---|---|---|
| `f * g` | 8.19e-12 | 1.43e-11 |
| `f * F` | 3.73e-09 | 5.82e-09 |
| `F * G` | 1.91e-06 | 2.86e-06 |

Both break at the same place in the T9 sweep (`B = 10⁷`), so the exactness
boundary is a property of the approach rather than of either implementation.

### A finding: contraction breaks commutativity

T5 checks `a*b == b*a` bit-exactly. The from-scratch code failed it at every
`logn > 1` until floating-point contraction was disabled.

The complex product's imaginary part is `ar·bi + ai·br`. With contraction
enabled the compiler fuses this into `fma(ar, bi, ai·br)`, which keeps
`ar·bi` exact and rounds `ai·br`. Exchanging the operands keeps the *other*
product exact, so the two orderings can differ in the last bit.

The reference does not have this problem because `inner.h` disables
contraction explicitly, noting that reproducibility requires it. `myfft.c`
does not include `inner.h` and so had to repeat the treatment. The lesson
generalises: an FFT that must produce reproducible values cannot be compiled
under default floating-point rules.

## Performance

`make bench`. Each operation is calibrated to run at least 20 ms per batch and
timed over 7 batches with the minimum taken. Buffers are restored by `memcpy`
before each destructive operation, and the cost of that restore is measured
separately and subtracted.

### Per-operation cost at `n = 512` (nanoseconds)

| | forward | inverse | pointwise |
|---|---|---|---|
| from scratch | 1478 | 1644 | 95.3 |
| reference | 1501 | 1002 | 94.0 |

The forward transforms and the pointwise product are level. The entire
difference between the two implementations is in the inverse transform, where
the reference is about 1.6× faster — it elides work in the final stages that
the straightforward DIT loop does not.

### Complete integer multiply

| n | from scratch | reference | ratio | schoolbook | speedup |
|---|---|---|---|---|---|
| 8 | 81 | 43 | 1.91 | 53 | 0.7× |
| 16 | 152 | 99 | 1.54 | 199 | 1.3× |
| 64 | 503 | 374 | 1.34 | 3 985 | 7.9× |
| 512 | 4 714 | 4 095 | 1.15 | 224 188 | 47.6× |
| 1024 | 9 955 | 9 109 | 1.09 | 880 188 | 88.4× |

Two things worth reading off this. The from-scratch implementation closes on
the reference as `n` grows, from 1.9× slower at `n = 8` to 1.09× at
`n = 1024`, because the fixed per-call overhead amortises. And the FFT does
not pay for itself until `n = 16` — below that the exact `O(N²)` schoolbook is
simply faster, which is why the asymptotic argument needs the constant factors
attached before it means anything.

### Scaling

Normalised cost, ns per `n·log₂n` for the FFT and per `n²` for the schoolbook:

| n | from scratch | reference | schoolbook |
|---|---|---|---|
| 32 | 1.619 | 1.102 | 1.025 |
| 128 | 1.154 | 0.923 | 0.951 |
| 512 | 1.023 | 0.889 | 0.855 |
| 1024 | 0.972 | 0.890 | 0.839 |

Both columns flatten, which is the empirical confirmation that the
implementations really are `O(N log N)` and `O(N²)` respectively. The
residual downward drift is fixed overhead amortising, not a change in
complexity class.

## Memory

`make memory` measures peak RSS, and on Linux also runs a `massif` heap
profile and a `memcheck` leak check. `make bench` prints the static
accounting, which is the figure that actually characterises the algorithm:

| | bytes |
|---|---|
| transform working buffer, either implementation (`n = 512`) | 4 096 |
| from-scratch twiddle tables, `logn = 9` | 6 144 |
| from-scratch twiddle tables, all sizes 1–10 | 24 560 |
| reference static tables (`fpr_gm_tab` + `fpr_p2_tab`) | 16 472 |

**Neither transform allocates.** Both work in place on caller-provided
buffers; `fft.c` contains no `malloc`, no `calloc`, and no stack arrays. The
working set for a Falcon-512 multiply is therefore two 4 KB operand buffers
plus a constant table, which is what makes the scheme viable on embedded
targets. Process peak RSS (1.4–1.7 MB across the three harnesses) is
dominated by the C runtime and the harnesses' own test buffers, and says
nothing useful about the transform — it is reported only to confirm the
absence of unexpected growth.

## Platform note

Built and verified on arm64 (Apple silicon) with clang, where `FALCON_AVX2` is
necessarily disabled — the scalar code paths in `fft.c` are what execute.
On x86-64 the AVX2 branches activate and should be re-verified separately;
the reference would be expected to gain more from that than the from-scratch
code, so the performance ratios above are the arm64 figures and not
transferable.

Three things need re-running on x86-64 Ubuntu with GCC before the numbers
here can be quoted as final:

- the verification harnesses, to exercise the AVX2 paths in `fft.c`
- the benchmark, since every timing above is arm64 scalar
- `make memory`, whose `massif` and `memcheck` stages need a working
  valgrind and so are skipped on Apple silicon

## License

The extracted Falcon sources retain their original MIT license; see the
header of each file. Material written for this study is offered under the
same terms.
