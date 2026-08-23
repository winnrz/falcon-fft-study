# References

Sources for the "Comparison with published work" section of the main README.
Figures are reproduced here so the claims can be checked without opening the
sources, and so they stay fixed even if the upstream numbers move.

Accessed 23 August 2026.

---

## Pornin, *Falcon on ARM Cortex-M4: an Update* (2025)

Thomas Pornin, NCC Group. 29 January 2025. IACR ePrint 2025/123.
<https://eprint.iacr.org/2025/123>

Local copy: `pornin-2025-falcon-cortex-m4-update.pdf` (not committed — see
`.gitignore`; re-fetch from the ePrint URL, which is stable).

Written by Falcon's own designer, which makes it the strongest available
authority on the questions this study touches.

**Setup.** STM32F407G-DISC1 board, STM32F407VGT6 MCU, measured at 24 MHz
"as is customary in the literature"; at that frequency caches can be disabled
because flash and RAM complete with minimal latency. GCC 13.2.1 with
`-O2 -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16`.
Constant-time by design. Averages over 15 000 key pairs and signatures.

**Table 1 — Falcon-512, cycles (RAM in bytes, transient + stack):**

| operation | speed | RAM |
|---|---|---|
| keygen | 71 943 764 | 13 343 + 928 |
| sign | 22 008 433 | 39 967 + 2 048 |
| verify (orig) | 255 306 | 2 079 + 872 |
| verify (BUFF) | 358 517 | 2 079 + 872 |

**Figure 2 — where signature generation spends its time (Falcon-512):**

| | share |
|---|---|
| `fpr_add` | 22.78% |
| `fpr_mul` | 21.81% |
| `fpr_add_sub` | 14.89% |
| Keccak-f (SHAKE) | 10.67% |
| `fpr_div` | 4.79% |
| `fpr_scaled` | 1.85% |
| `fpr_sqrt` | 1.42% |

Floating point therefore accounts for roughly 67% of a Falcon-512 signature.
The addition family at 37.7% combined exceeds `fpr_mul` at 21.8%. These are
aggregate shares across all calls, not per-call costs.

**Two quotations this study relies on directly.**

On the FPU, which settles the broad single-precision question:

> the used M4 CPU includes the optional DSP and floating-point units; the
> latter is not directly usable for Falcon signature generation, since it
> supports only 32-bit precision (IEEE 754 type "binary32"); however, the
> floating-point registers are a convenient fast-access storage area, which
> the C compiler uses, and which we also leverage in our assembly routines.

On verification, which does **not** use the floating-point FFT:

> if using only signature verification, then the code footprint totals 11840
> bytes. This includes 4096 bytes for the NTT tables

Implementation: <https://github.com/pornin/c-fn-dsa>

---

## pqm4

<https://github.com/mupq/pqm4>

Snapshot: `pqm4-benchmarks-2026-08-23.md`, taken from
<https://raw.githubusercontent.com/mupq/pqm4/master/benchmarks.md>. Committed
because upstream is a moving target and a dissertation needs a fixed
reference point.

**Setup.** Default board Nucleo-L4R5ZI. "All cycle counts were obtained at
24MHz to avoid wait cycles due to the speed of the memory controller."

**Naming.** Falcon appears as **FN-DSA**, its standardised name —
`fndsa_provisional-512` and `-1024`. Searching pqm4 for "falcon" finds
nothing.

**fndsa_provisional-512, cycles (10 executions):**

| impl | keygen | sign | verify |
|---|---|---|---|
| m4f | AVG 67 693 338 / MIN 57 825 106 | AVG 22 469 685 / MIN 22 280 159 | AVG 396 949 / MIN 390 368 |
| ref | AVG 85 699 591 / MIN 64 822 516 | AVG 49 522 949 / MIN 49 325 465 | AVG 731 387 / MIN 714 301 |

**Memory (bytes):** m4f keygen 27 772, sign 41 952, verify 2 976;
ref keygen 27 676, sign 82 276, verify 5 308.

---

## Lead not followed up

*Optimized Falcon Verify on Cortex-M4 for Post-Quantum secure UAV
communications*, ScienceDirect S2405959524001401. Reports Falcon-512
verification gains over pqm4 and the reference using NTT-based polynomial
multiplication with Plantard modular arithmetic. Not consulted for this
study, since verification uses the integer NTT rather than the floating-point
FFT under study here. Relevant if the scope ever widens to verification.

---

## BibTeX

```bibtex
@misc{pornin2025falconm4,
  author       = {Thomas Pornin},
  title        = {Falcon on {ARM} {Cortex-M4}: an Update},
  howpublished = {Cryptology {ePrint} Archive, Paper 2025/123},
  year         = {2025},
  month        = jan,
  url          = {https://eprint.iacr.org/2025/123},
  note         = {Accessed 2026-08-23}
}

@misc{pqm4,
  author       = {Matthias J. Kannwischer and Richard Petri and
                  Joost Rijneveld and Peter Schwabe and Ko Stoffelen},
  title        = {{pqm4}: Testing and Benchmarking {NIST} {PQC} on
                  {ARM} {Cortex-M4}},
  howpublished = {\url{https://github.com/mupq/pqm4}},
  note         = {Benchmark snapshot accessed 2026-08-23}
}
```
