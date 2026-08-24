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

**The platform of record is an STM32F411 (Cortex-M4F), a Nucleo-F411RE
clocked at 24 MHz.** Every performance figure quoted as a result is measured
there, in cycles. The development host is used for fast iteration during
development and its timings are reported separately and clearly labelled;
they are not the result. The choice matters more than it might appear: the
M4F's FPU is single precision only, so the binary64 arithmetic Falcon's FFT
requires has no hardware behind it and is emulated. That turns out to change
not just the magnitude of the costs but their ordering.

## Headline results

At Falcon-512's operating point — `n = 512`, coefficients within ±25 — on a
Nucleo-F411RE at 24 MHz, measured in cycles with the DWT counter:

| build | cycles | vs build 2 | spread over operand values |
|---|---|---|---|
| 1 reference, Falcon's emulated binary64 in M4 assembly | 3 920 995 | 1.36× | **0** |
| 2 reference, libgcc soft-double | 2 890 054 | — | 8 727 |
| 3 from scratch, libgcc soft-double | 3 101 545 | 1.07× | 9 443 |

- **Constant-time arithmetic costs 1.36×.** Build 1 against build 2 isolates
  it: same FFT, different backend.
- **The from-scratch FFT is 7% behind the reference** on identical
  floating-point routines. Builds 2 and 3 isolate that: same backend,
  different FFT.
- **The constant-time property is measured, not asserted.** The spread column
  is execution-time variation across cases that differ only in operand
  values. Build 1 does not vary at all; the other two do.
- **Precision headroom is 37 bits** at the operating point — the worst
  observed deviation could double 37 times before exact integer recovery
  would fail.
- All three builds are verified on hardware against exact integer vectors
  produced by a `__int128` schoolbook oracle.

A finding worth flagging early: measured on a desktop, the reference's
inverse transform is markedly cheaper than its forward transform, and that is
where its whole advantage lies. Under emulation that advantage
[disappears](#the-references-inverse-transform-advantage-does-not-survive-emulation).
A study run only on a desktop would have drawn the opposite conclusion about
the code that Falcon actually ships to embedded targets.

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
| `bench_fftmul.c` | Performance measurement on the host, both implementations |
| `measure_memory.sh` | Peak RSS, heap profile, leak check |
| `Makefile` | Builds and runs everything on the host |
| `target/` | The Cortex-M4 firmware: CubeMX project plus `target/App` |
| `tools/gen_kat.c` | Generates the target's known-answer vectors |
| `docs/BUILDING.md` | Toolchain, target project, and reproduction notes |
| `bench.csv`, `bench-m4.csv` | Host and target timings |
| `verify-m4.csv` | Target correctness and precision margins, all three builds |

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
1 otherwise. [`docs/BUILDING.md`](docs/BUILDING.md) covers the toolchain, the
target project layout, regenerating the known-answer vectors, and the edits
that must be re-applied after regenerating from CubeMX.

On the target, one firmware image contains all three builds and runs the
verification and the benchmark in sequence, reporting over the ST-Link's
virtual COM port at 115200:

```sh
cd target && make
st-flash --connect-under-reset --reset write build/target.bin 0x8000000
```

`--connect-under-reset` is required; without it the write fails to attach
even though `st-info --probe` reads the chip. The exact integer oracle cannot
cross-compile, since `__int128` does not exist on 32-bit ARM, so
`tools/gen_kat.c` runs it on the host and freezes the results into
`target/App/kat.c`. The oracle is unchanged in strength; only where it runs
has moved.

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
| forward transform, slot by slot | 7.28e-12 |
| inverse, completing the round trip | 5.68e-13 |
| finished product, coefficient by coefficient | 5.12e-09 |

113 checks, 0 failures. The forward comparison is the one that constrains slot
ordering — a permuted layout would still pass the product comparison. The
differences are not expected to reach zero: the two implementations apply the
same operations in a different order, so they round differently even when both
are equally accurate.

Accuracy against the exact integer oracle is what settles which is better, and
on the host the two are now identical at Falcon-512's operating point:

| operand pair | reference | from scratch |
|---|---|---|
| `f * g` | 8.185e-12 | 8.185e-12 |
| `f * F` | 3.725e-09 | 3.725e-09 |
| `F * G` | 1.907e-06 | 1.907e-06 |

This was not true of an earlier version of the implementation, which was
consistently about 1.5× less accurate. The cause was a rounding error in the
twiddle **arguments** rather than in `cos` and `sin` themselves; the diagnosis
and the fix are described under [Precision headroom](#precision-headroom).

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

All figures in this section are Cortex-M4 cycles at 24 MHz with zero flash
wait states and interrupts masked, minimum of three repeats, from
`bench-m4.csv`. The restore `memcpy` each destructive operation needs sits
outside the timed region rather than being subtracted from it.

### Why there are three builds

Two things vary independently on this target: whose FFT runs, and which
software routine performs each `double` operation underneath it. The M4F
cannot execute binary64 in hardware, so there is always a software routine.

| build | FFT | binary64 backend | constant-time |
|---|---|---|---|
| 1 | reference | Falcon's own, in hand-written M4 assembly | yes |
| 2 | reference | libgcc soft-double | no |
| 3 | from scratch | libgcc soft-double | no |

`config.h` selects build 1 automatically: it autodetects
`FALCON_ASM_CORTEXM4` on any Cortex-M4 target and that forces the emulated
backend on. Build 2 requires `FALCON_FPNATIVE` to be set explicitly. Build 3
is unaffected by either, because `myfft.c` uses the plain C `double` type and
never consults `config.h`.

Only two comparisons are meaningful. **Builds 2 and 3 share a backend**, so
the difference between them is FFT structure alone. **Builds 1 and 2 share an
FFT**, so the difference between them is the price of constant-time
arithmetic. Comparing 1 against 3 conflates the two and says nothing.

### Per-operation cost at `n = 512` (cycles)

| build | forward | inverse | pointwise | complete |
|---|---|---|---|---|
| 1 reference, emulated + asm | 1 228 583 | 1 294 571 | 169 258 | 3 920 995 |
| 2 reference, soft-double | 911 274 | 940 414 | 127 465 | 2 890 054 |
| 3 from scratch, soft-double | 977 517 | 1 019 980 | 127 359 | 3 101 545 |

Constant-time arithmetic costs **1.357×**. The from-scratch implementation is
**1.073×** the reference on identical floating-point routines.

The pointwise products agree to within 0.1% — 127 359 against 127 465 — which
is what should happen, since that operation is the same arithmetic on the same
soft-double routines in both. The entire structural gap is therefore in the
transforms: forward +7.3%, inverse +8.5%.

### The reference's inverse-transform advantage does not survive emulation

Measured on the development host, the reference's inverse transform is
markedly *cheaper* than its own forward transform, and that is where its whole
advantage over the from-scratch code lay. On the target that advantage is
gone:

| inverse ÷ forward | host | target |
|---|---|---|
| reference | 0.675 | 1.032 |
| from scratch | 1.114 | 1.043 |

Same source, opposite conclusion. The reference elides work in the final
stages of the inverse transform, and on hardware doubles, where a multiply and
an addition cost about the same, eliding operations is simply a saving. Under
emulation the two are no longer comparable: software binary64 addition has to
align exponents and renormalise, and is not obviously cheaper than software
multiplication. An optimisation that trades one for the other therefore stops
paying.

This is the clearest single reason the platform of record matters. A study
that measured only on a desktop would have concluded that the from-scratch
implementation's inverse transform was its one real weakness, and would have
been wrong about the hardware Falcon actually targets.

### Scaling

Normalised cost, cycles per `n·log₂n` for the complete multiply:

| n | build 1 | build 2 | build 3 |
|---|---|---|---|
| 32 | 818.6 | 604.0 | 676.3 |
| 128 | 839.0 | 617.7 | 673.1 |
| 512 | 850.9 | 627.2 | 673.1 |
| 1024 | 855.2 | 630.5 | 673.4 |

Flat across a 32× range of `n`, which is the empirical confirmation that all
three are `O(N log N)`. Taken as a ratio between the two largest sizes,
against an ideal of 2.222:

| build 1 | build 2 | build 3 |
|---|---|---|
| 2.2334 | 2.2340 | 2.2234 |

### Falcon-1024 (`logn = 10`)

Every measurement above was repeated at `n = 1024`, the parameter set for
Falcon-1024. The complete multiply:

| build | cycles | at 24 MHz | vs build 2 |
|---|---|---|---|
| 1 reference, emulated + asm | 8 757 193 | 365 ms | 1.356× |
| 2 reference, soft-double | 6 456 246 | 269 ms | — |
| 3 from scratch, soft-double | 6 895 876 | 287 ms | 1.068× |

Nothing qualitative changes at the larger size, which is the point of
measuring it. The constant-time premium is 1.356× against 1.357× at `n = 512`
— unchanged to three figures — so the cost of Falcon's emulated arithmetic is
a property of the arithmetic and not of the transform length. The structural
gap between the two FFTs narrows slightly, 1.073× to 1.068×, with the same
shape as before: the pointwise products stay within 0.02% of each other
(254 084 against 254 120) while the transforms account for the whole
difference, forward +6.7% and inverse +7.9%.

Precision behaves as the scaling argument predicts. One extra stage costs one
more bit of headroom at the operating point, 37 at `logn = 9` down to 36 at
`logn = 10` for both implementations, and the wide-bound cases sit at 23 bits
in both. Falcon-1024 therefore retains a very large margin — around 2³⁶ times
the deviation needed to break exact recovery.

The practical reading is that doubling `n` costs slightly more than the ideal
2.222× predicted by `N log N`, at 2.233× for the reference and 2.223× for the
from-scratch code, and that a Falcon-1024 signature's polynomial multiply
costs roughly a quarter of a second on this part at 24 MHz.

### Constant time, measured rather than asserted

The three builds were run against known-answer cases that differ only in
their operand values, at the same size. Execution-time spread across those
cases:

| build | spread |
|---|---|
| 1 reference, emulated + asm | **0 cycles** |
| 2 reference, soft-double | 8 727 cycles |
| 3 from scratch, soft-double | 9 443 cycles |

Build 1 does not vary at all with the data. Builds 2 and 3 do. This is the
security property observed directly on the hardware rather than inferred from
the source.

A second, independent check says the same thing. Measuring the four
operations separately and adding them up should reproduce the separately
measured complete multiply — but the isolated measurements run over different
intermediate values than the combined one, so only a data-independent
implementation can be exactly additive:

| build | 2·forward + pointwise + inverse | complete | difference |
|---|---|---|---|
| 1 | 3 920 995 | 3 920 995 | **0** |
| 2 | 2 890 427 | 2 890 054 | 373 |
| 3 | 3 102 373 | 3 101 545 | 828 |

Build 1 agrees to the cycle. Neither of the others does.

### Precision headroom

Rounding recovers the exact integer product only while a coefficient sits
closer than 0.5 to it. The margin below is how many times the worst observed
deviation could be doubled before reaching that limit — the headroom in bits,
measured on target across ten known-answer cases spanning `logn` 5 to 10. The
full table is `verify-m4.csv`:

| logn | operand bound | reference | from scratch |
|---|---|---|---|
| 5 | ±25 | 40 | 39 |
| 6 | ±25 | 39 | 38 |
| 7 | ±25 | 38 | 38 |
| 8 | ±25 | 37 | 37 |
| 9 | ±25 | 37 | 36 and 37 |
| 10 | ±25 | 36 | 36 |
| 5 | ±2048 | 27 | 27 |
| 9 | ±2048 | 23 | 24 |
| 10 | ±2048 | 23 | 23 |

Two cases are listed at `logn = 9, ±25` because two independent operand draws
were measured there; they disagree by a bit, which is itself informative about
how finely this metric resolves.

#### Where the from-scratch code lost a bit, and why

An earlier version of this implementation trailed the reference by a bit at
most sizes. The cause turned out not to be the one that first suggests itself.
`myfft.c` builds its twiddles by calling `cos` and `sin` at run time where the
reference reads the precomputed `fpr_gm_tab`, so the natural suspicion is that
the library functions are the weak link. Measured against a 60-digit
reference, they are not: on the development host they are very nearly
correctly rounded, contributing about 2⁻⁵³·⁵ of the error.

The error was in the **argument**, before either function was called.
`M_PI * (double)v / (double)n` rounds twice — `M_PI` is already a rounded π,
and the product rounds again — so the angle handed to `cos` was itself several
ulps from the true one, and no accuracy inside the function can recover an
argument that arrived wrong. The twiddle table carried about 2⁻⁵¹·⁴ of error
against a correctly rounded table's 2⁻⁵⁴·⁰, and every butterfly that read it
inherited the difference.

Carrying the angle as an unevaluated sum `hi + lo` removes the second
rounding: π splits into a head and a tail, `fma` recovers the exact residual
of the product, the division by `n` is exact because `n` is a power of two,
and a first-order expansion transfers the residual onto the result. On the
host this closes the gap **exactly** — the from-scratch code now matches the
reference to the last digit at Falcon-512's operating point, where it was
previously 1.5× worse.

On target the improvement is real but smaller: summed across the ten cases
above, the deficit against the reference falls from four bits to one. Two
things account for the difference. The margin is reported as a whole number of
bits while the improvement is worth roughly 0.6 of one, so individual cases
flip either way depending on which side of a rounding boundary they fall —
which is why `logn = 6` moved the wrong way while `logn = 9, ±2048` moved a
bit ahead of the reference. And the target builds its table with newlib's
`cos`, `sin` and `fma` rather than the host's, so the two platforms do not
construct the same table; the residual on target is bounded by newlib's
library quality, which this study has not measured independently.

The honest summary is that the host measurement identifies the mechanism and
the target measurement is consistent with it, but the on-target metric is too
coarsely quantised to confirm it case by case. Eliminating the remaining
difference would mean precomputing the table off-line as `static const
double[]`, which would also make host and target bit-identical and remove the
run-time allocation; it was judged not worth reopening verified code for a
sub-bit effect the measurement cannot resolve.

Thirty-seven bits is a large margin, and it bears directly on whether the
idle single-precision FPU could take over part of the transform: `float`
carries 24 mantissa bits against `double`'s 53, so a naive substitution costs
about 29. That leaves roughly 8 bits in hand. Suggestive, not conclusive —
error growth is not exactly linear in mantissa bits — but enough to say the
question is worth measuring.

### Host figures

Retained for comparison, not as results. Measured on arm64 (Apple silicon)
with clang, in nanoseconds, from `bench.csv`. `FALCON_AVX2` is necessarily
disabled there, so the scalar paths in `fft.c` are what execute.

Complete integer multiply, including the exact schoolbook for reference:

| n | from scratch | reference | ratio | schoolbook | speedup |
|---|---|---|---|---|---|
| 8 | 81 | 43 | 1.89 | 53 | 0.7× |
| 16 | 154 | 98 | 1.57 | 200 | 1.3× |
| 64 | 502 | 366 | 1.37 | 3 944 | 7.9× |
| 512 | 4 728 | 3 783 | 1.25 | 224 688 | 47.5× |
| 1024 | 9 902 | 8 340 | 1.19 | 879 062 | 88.8× |

These were re-measured after the twiddle change described above; the
from-scratch timings moved by less than 0.5%, since table construction happens
outside every timed region. The reference's figures moved more, which is
session-to-session drift on a shared laptop rather than anything about the
code — one more reason the host is not the platform of record here.

The point that does transfer: **the FFT does not pay for itself until
`n = 16`.** Below that the exact `O(N²)` schoolbook is simply faster, which is
why the asymptotic argument needs its constant factors attached before it
means anything. The schoolbook was not re-measured on target.

## Comparison with published work

Two sources report Falcon on Cortex-M4 in enough detail to set these figures
against: Thomas Pornin's *Falcon on ARM Cortex-M4: an Update*
([eprint 2025/123](https://eprint.iacr.org/2025/123.pdf)), by the algorithm's
own author, and the
[pqm4](https://github.com/mupq/pqm4) benchmarking project, which carries
Falcon under its standardised name FN-DSA.

Full citations, the figures relied on, and a snapshot of the pqm4 benchmark
table as it stood on 23 August 2026 are in [`refs/`](refs/README.md).

### The methodology matches

| | this study | Pornin 2025 | pqm4 |
|---|---|---|---|
| board | Nucleo-F411RE | STM32F407G-DISC1 | Nucleo-L4R5ZI |
| clock | 24 MHz | 24 MHz | 24 MHz |
| unit | cycles | cycles | cycles |
| compiler | GCC 15.3, `-O3` | GCC 13.2, `-O2` | GCC |

The 24 MHz choice made here to eliminate flash wait states turns out to be the
established convention. Pornin measures "at 24 MHz speed, as is customary in
the literature… at that relatively low frequency (the board can be used at up
to 168 MHz), caches can be disabled because RAM and ROM (Flash) accesses
normally complete with minimal latency". pqm4 gives the same reason: "All
cycle counts were obtained at 24MHz to avoid wait cycles due to the speed of
the memory controller." The boards differ, so the memory subsystems differ,
but the figures are comparable in the way this literature treats them.

### What the published figures measure is not what this study measures

This is the caveat that governs everything below. Both sources report **whole
scheme operations** — key generation, signing, verification. This study
measures **one primitive**, a single `n = 512` polynomial multiply. The
numbers cannot be equated, only placed in context.

Falcon-512, cycles at 24 MHz:

| | keygen | sign | verify |
|---|---|---|---|
| Pornin 2025 | 71 943 764 | 22 008 433 | 255 306 (orig) / 358 517 (BUFF) |
| pqm4 `m4f`, min | 57 825 106 | 22 280 159 | 390 368 |
| pqm4 `ref`, min | 64 822 516 | 49 325 465 | 714 301 |

The two optimised signing figures agree closely — 22.0 M against 22.3 M —
across different boards, which is a useful check on the comparability of the
convention. The gap between pqm4's `ref` and `m4f` signing, 49.3 M against
22.3 M, is the 2.2× that hand optimisation buys, consistent with Pornin's
report that the current code is "about twice faster than the 2019 code".

Against that, the complete `n = 512` multiply measured here at
**3 920 995 cycles** is roughly **18% of an optimised Falcon-512 signature**,
or 8% of the reference implementation's. That is the honest way to place it:
one polynomial multiply is a substantial minority of a signature, and signing
performs several along with Gaussian sampling and the LDL tree.

One further scoping point: Falcon **verification does not use this FFT at
all.** Pornin's verify path uses an integer NTT modulo *q*, with 4 096 bytes
of NTT tables. The floating-point FFT studied here is used in key generation
and signing.

### Independent support for the addition-versus-multiplication finding

The claim above — that eliding operations stops paying once binary64 is
emulated, because software addition must align exponents and renormalise — is
corroborated by Pornin's own profile of where signing spends its time
(Figure 2, Falcon-512):

| operation | share of signature generation |
|---|---|
| `fpr_add` | 22.78% |
| `fpr_mul` | 21.81% |
| `fpr_add_sub` | 14.89% |
| `fpr_div` | 4.79% |
| `fpr_scaled` | 1.85% |
| `fpr_sqrt` | 1.42% |
| Keccak-f (SHAKE) | 10.67% |

Floating-point arithmetic is **about 67% of a Falcon-512 signature** on this
platform, which is the strongest available justification for studying it in
isolation. And the addition family, `fpr_add` plus `fpr_add_sub` at 37.7%
combined, costs considerably more in aggregate than `fpr_mul` at 21.8%. These
are aggregate shares rather than per-call costs, so they do not prove the
mechanism, but they are consistent with it and inconsistent with the
assumption that a saved multiplication is worth more than an added addition.

### On the single-precision question

The precision headroom measured above — 37 bits at Falcon-512's operating
point, against roughly 29 that `float` would cost — invites the question of
whether the idle single-precision FPU could take over.

For Falcon as a whole, the algorithm's author has already answered it. Pornin
states plainly that the FPU "is not directly usable for Falcon signature
generation, since it supports only 32-bit precision (IEEE 754 type
binary32)", and uses the floating-point registers only as fast-access
storage.

That does not contradict the margin measured here, because the two concern
different operations. The 37 bits is headroom in the **polynomial multiply**.
Falcon's precision-critical work is elsewhere — the Gaussian sampling and the
LDL tree in signing — and nothing in this study measures those. The open
question is therefore the narrow one: whether single precision could serve for
the polynomial multiply alone, inside an implementation that keeps binary64
where the algorithm actually demands it. The broad version is settled, and
should not be presented as open.

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

That embedded claim is now measured rather than argued. Cross-compiled for
Cortex-M4 at `-O3`, the transform code occupies 4 468 bytes of flash for
`fft.c`, 21 704 for `fpr.c` — almost all of it the constant tables — and
2 672 for `myfft.c`, against the F411's 512 KB. A Falcon-512 multiply needs
two 4 KB operand buffers of its 128 KB SRAM. The verification and benchmark
firmware is far larger than any of this, because it carries all three builds,
the frozen known-answer vectors, and a 32 KB heap for `myfft.c`'s twiddle
tables at every size; none of that is the algorithm.

## Platform note

The target is a Nucleo-F411RE: STM32F411RE, Cortex-M4F, 512 KB flash,
128 KB SRAM, driven at 24 MHz from the ST-Link's 8 MHz clock through the PLL.

**24 MHz rather than the part's 100 MHz maximum, deliberately.** Above 90 MHz
the F411 needs three flash wait states, at which point a large and variable
share of every measured cycle is the core stalled on flash and dependent on
whether the ART accelerator happened to hold the loop. Below 30 MHz there are
no wait states at all, so the counts reflect the transform rather than the
memory system. The same reasoning is why comparable published work on these
parts also clocks them down.

Timing is the DWT cycle counter, validated before any FFT code was flashed:
measurement overhead one cycle, and a busy loop of `2n` iterations costing
1.999× one of `n`. The core clock was cross-checked against SysTick, which
derives from HCLK by an independent path.

Two caveats when setting these numbers against published work. The usual
reference board for this kind of measurement is the STM32F407 — 168 MHz,
192 KB SRAM — not the F411 used here, so cycle counts are comparable but the
parts are not identical. And the reference's AVX2 paths never execute on
either platform, so nothing here says anything about x86-64 performance.

## License

Material written for this study is MIT licensed. The vendored Falcon sources
and the STMicroelectronics HAL retain their own terms. See
[`LICENSE`](LICENSE) for the details of each.
