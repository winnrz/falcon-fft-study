# CLAUDE.md

Working notes for this repository. `README.md` documents the study itself —
the maths, the test design, and the results. This file covers only what is
needed to build, run and extend the code correctly, and the things that are
easy to get wrong.

## What this project is

An MSc study (7V0007, MMU) of FFT-based polynomial multiplication in
`Z[X]/(X^N + 1)`, the operation at the heart of Falcon-512. It compares a
from-scratch implementation (`myfft.c`) against the official Falcon
reference (`fft.c`, `fpr.c`), for both correctness and cost.

**Target platform: STM32F411 (Cortex-M4F).** The x86-64 / Ubuntu / GCC
baseline named in the signed proposal was dropped deliberately in August 2026
after a supervisor meeting. Do not propose reinstating it. Note the signed
proposal still says Ubuntu/GCC, so a proposal amendment may be outstanding.

Host builds remain useful for fast verification during development, but the
Cortex-M4F is the platform of record for every benchmark number.

## Layout

| file | role |
|---|---|
| `fft.c`, `fpr.c`, `inner.h`, `fpr.h`, `config.h` | official Falcon reference, unmodified |
| `myfft.c`, `myfft.h` | from-scratch implementation, shares no code with the reference |
| `test_fftmul.c` | verifies the reference (T1–T10) |
| `test_myfft.c` | verifies `myfft.c` (T1–T10) plus T11, a slot-by-slot cross-check against the reference |
| `bench_fftmul.c` | timing, host only so far |
| `measure_memory.sh` | peak RSS, massif, leak check |

## Conventions

- C99, tabs for indent (8-wide), BSD style: return type on its own line
  above the function name. Comments are block `/* */`, wrapped near column 72.
  Match the surrounding density — these files are heavily commented on purpose.
- The Makefile is `.POSIX`. Keep it portable; no GNU-make-only constructs.
- **`test_fftmul.c` is frozen.** It is a verified artefact whose results are
  quoted in the write-up. The exact `__int128` schoolbook oracle is duplicated
  into `test_myfft.c` rather than shared, deliberately. Do not refactor the
  two harnesses together.
- `myfft.c` deliberately reproduces the reference's FFT slot layout so that
  intermediate values can be compared coefficient by coefficient. Do not
  "simplify" the layout — T7 and T11 both depend on it.

## Building for the host

```
make            build every harness
make test-all   run both verification harnesses
make bench      timings
make test-both  reference under both float backends
```

## Building for Cortex-M4

Toolchain (macOS):

```
brew install --cask gcc-arm-embedded    # arm-none-eabi-gcc
brew install stlink                     # st-info, st-flash
```

Core flags:

```
-mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16
```

All three FFT sources cross-compile clean with zero warnings under
`-Wall -Wextra -Wshadow -Wundef`. Sizes at `-O3`: `fft.o` 4468 B text,
`fpr.o` 21704 B, `myfft.o` 2400 B text + 312 B bss. Comfortable against the
F411's 512 KB flash / 128 KB SRAM.

### The float backend gotcha

`config.h` autodetects `FALCON_ASM_CORTEXM4 = 1` whenever the target is
`-mcpu=cortex-m4`, and that **forces `FALCON_FPEMU = 1` by default**. The host
default is the opposite. Consequences:

- `-DFALCON_FPEMU=1` is a **no-op on target**. Objects built with and without
  it are byte-identical. The `emu` Makefile target's model inverts here.
- The on-target reference path is emulated binary64 **in hand-written M4
  assembly**, which is *not* the plain-C emulation verified bit-exact on the
  host. It must pass T1–T10 on target in its own right before any timing
  taken from it is quoted.

### The three-build matrix

The M4F FPU is single-precision only (`fpv4-sp-d16`), so `double` has no
hardware support and every double operation is emulated. Two things vary
independently — whose FFT runs, and which software-double routine it calls:

| # | FFT code | double backend | flags | constant-time |
|---|---|---|---|---|
| 1 | reference | Falcon's own, M4 assembly | *(none — the default)* | yes |
| 2 | reference | libgcc soft-double | `-DFALCON_FPNATIVE=1` | no |
| 3 | `myfft.c` | libgcc soft-double | *(always — uses plain `double`)* | no |

- **2 vs 3** is the fair algorithm comparison: identical float backend, so any
  gap is FFT structure alone.
- **1 vs 2** prices the constant-time emulation.

Quoting 1 against 3 conflates the two effects. Don't.

## Benchmarking on target

- Measure **cycles** via the DWT counter (`CoreDebug->DEMCR |= TRCENA`, then
  `DWT->CYCCNT`), not wall time. Cycles are what the comparable literature
  reports.
- **Clock the board at 24 MHz, not 100 MHz.** The F411 needs 3 flash wait
  states above 90 MHz, at which point a large share of measured cycles is
  flash stall and ART cache behaviour rather than FFT cost. 24 MHz gives zero
  wait states. This also matches the pqm4 convention — verify the exact figure
  against pqm4 before quoting comparability.
- The STM32 HAL enables **SysTick at 1 kHz** by default. Disable its interrupt
  around timed regions or it injects noise into a measurement whose method is
  "take the minimum".
- Carry over the method in `bench_fftmul.c`: calibrate each operation to a
  minimum batch duration, time several batches, keep the **minimum**, and
  measure-and-subtract the input-restore `memcpy` rather than assuming it is
  negligible.

### Porting obstacles

- `__int128` does not exist on 32-bit ARM. The exact schoolbook oracle (T8)
  cannot cross-compile as written — either implement it on 64×64→128 helpers
  or keep T8 host-side and flash only the properties that fit.
- `myfft.c` calls `malloc` for twiddle tables (`my_fft_init`) and for
  per-call scratch inside `my_poly_mul` (`myfft.c:430`). On target, prefer
  precomputing the twiddles on the host and emitting them as
  `static const double[]` so they live in flash rather than SRAM — this also
  drops the `cos`/`sin` calls and the double-precision libm they pull in.
  The per-call scratch `malloc` sits inside the timed path and must go.

## Hardware

Nucleo-F411RE, onboard ST-LINK/V2-1.

```
st-info --probe
```

should report `chipid: 0x431`, `dev-type: STM32F411xC_xE`, 512 KB flash,
128 KB SRAM. UART output arrives over the ST-Link virtual COM port
(`/dev/cu.usbmodem*` on macOS, USART2 on the board side).

No mass-storage volume appears for this board, so flash over SWD with
`st-flash` rather than drag-and-drop.

The project scaffold is generated by **standalone STM32CubeMX** with
Toolchain/IDE set to **Makefile** — not STM32CubeIDE, which locks that
dropdown to itself and cannot emit a Makefile project.

## Stale documentation

`README.md`'s "Platform note" section still describes arm64/clang host results
and lists x86-64 Ubuntu re-runs as outstanding work. That section predates the
Cortex-M4 pivot and needs rewriting once on-target numbers exist.
