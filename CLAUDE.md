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
after a supervisor meeting. Do not propose reinstating it. The signed
proposal still names Ubuntu/GCC, but the supervisor is aware of the pivot and
no amendment is being pursued — do not raise it as an open item.

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
| `bench-m4.csv` | target timings, `logn` 1-10, all three builds |
| `verify-m4.csv` | target correctness and precision margins, 10 KAT cases |

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
`fpr.o` 21704 B, `myfft.o` 2672 B text + 312 B bss. Comfortable against the
F411's 512 KB flash / 128 KB SRAM.

### The float backend gotcha

`config.h` autodetects `FALCON_ASM_CORTEXM4 = 1` whenever the target is
`-mcpu=cortex-m4`, and that **forces `FALCON_FPEMU = 1` by default**. The host
default is the opposite. Consequences:

- `-DFALCON_FPEMU=1` is a **no-op on target**. Objects built with and without
  it are byte-identical. The `emu` Makefile target's model inverts here.
- The on-target reference path is emulated binary64 **in hand-written M4
  assembly**, which is *not* the plain-C emulation verified bit-exact on the
  host. It has since been checked in its own right against the exact integer
  oracle on target and passes; see `target/App/verify.c`.

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

### Porting obstacles, and how they were handled

- **`__int128` does not exist on 32-bit ARM**, so the exact schoolbook oracle
  cannot cross-compile. Reimplementing 128-bit arithmetic on target would put
  new, unverified code exactly where the trust is meant to come from, so the
  oracle stays on the host: `tools/gen_kat.c` runs it and freezes the results
  into `target/App/kat.c`. The oracle is unchanged in strength; only where it
  runs moved.
- **`myfft.c` allocates.** `my_fft_init` mallocs twiddle tables and
  `my_poly_mul` mallocs scratch on every call (`myfft.c:430`). Handled by
  raising `_Min_Heap_Size` and calling `myfft_prepare()` at startup so table
  construction sits outside every timed region, and by driving the transform
  steps directly instead of calling `my_poly_mul`, exactly as
  `bench_fftmul.c` does on the host. Still worth doing eventually:
  precomputing the twiddles on the host as `static const double[]` would move
  them from SRAM to flash and drop the `cos`/`sin` calls with their
  double-precision libm.
- **newlib-nano's `printf` supports neither `%f` nor `%lld`.** Anything the
  firmware reports must be expressible as a 32-bit integer. This is why the
  precision margin is reported in bits rather than as a scaled deviation.

## Current state

The port is working and all three builds are verified on hardware. At
Falcon-512's operating point (`logn = 9`, `|coeff| <= 25`, n = 512 multiply):

| build | cycles | vs build 2 | time spread over cases |
|---|---|---|---|
| 1 reference, FPEMU + M4 asm | 3 920 995 | 1.36x | **0 cycles** |
| 2 reference, libgcc soft-double | 2 890 054 | — | 8 727 cycles |
| 3 `myfft.c`, libgcc soft-double | 3 101 545 | 1.07x | 9 443 cycles |

Constant-time arithmetic costs 1.36x; the from-scratch FFT is 7% behind the
reference on identical float routines. The spread column measures the
constant-time property directly — those cases differ only in operand values.

Run-to-run noise on these figures is +/-1 cycle, so differences of a few
hundred are real. Full tables: `bench-m4.csv` (timings) and `verify-m4.csv`
(correctness and precision margins), both captured from the firmware that is
currently on the board.

Precision headroom is 37 bits for the reference at the operating point and
36-37 for `myfft.c`, over ten KAT cases spanning `logn` 5-10. See the twiddle
note below for why the two are no longer cleanly separated.

### What is done

- Port working, all three builds verified on hardware against exact integer
  vectors generated by `tools/gen_kat.c`.
- Per-operation benchmark harness sweeping `fft`, `ifft`, pointwise and the
  full multiply across `logn` 1–10 for all three builds, in one firmware
  image at one clock setting. Results in `bench-m4.csv`.
- `README.md` rewritten around the target results, with the host figures
  retained but explicitly labelled as not the result.
- Comparison against Pornin (ePrint 2025/123) and pqm4 written up, with
  sources and a pinned pqm4 snapshot in `refs/`.
- **The `myfft.c` precision gap is diagnosed and fixed.** The cause was not
  `cos`/`sin` — measured against a 60-digit reference those are very nearly
  correctly rounded. It was the *argument*: `M_PI * v / n` rounds twice, so
  the angle reaching libm was already several ulps out, and the table
  carried 2^-51.4 of error against a correctly rounded 2^-54.0.
  `cos_sin_pi` in `myfft.c` now carries the angle as an unevaluated
  `hi + lo` pair (pi split into head and tail, `fma` for the exact product
  residual, exact division by the power-of-two `n`, first-order expansion
  onto the result). On the host this closes the gap exactly — `myfft.c` now
  matches the reference to the last digit at the operating point. On target
  the summed deficit over ten cases falls from 4 bits to 1; it does not
  vanish because the margin is quantised to whole bits while the improvement
  is worth ~0.6, and because the target builds its table with newlib's
  `cos`/`sin`/`fma` rather than the host's, so the two platforms do not
  construct the same table.
- KAT coverage extended from 5 cases at two sizes to 10 spanning `logn` 5-10,
  so the precision trend can be read rather than inferred from two points.
- `verify-m4.csv` captured, so the precision claim has an artefact behind it
  the way the timings always did.

### Outstanding

1. **The single-precision question, in its narrow form only.** Pornin states
   the M4F's FPU "is not directly usable for Falcon signature generation,
   since it supports only 32-bit precision", so the broad version is settled
   by the algorithm's author — do not present it as open. What the 37-bit
   margin measured here actually speaks to is whether single precision could
   serve for the **polynomial multiply alone**, inside an implementation
   keeping binary64 for the Gaussian sampling and LDL tree, which this study
   does not measure. See `refs/README.md` for the quotation.
2. **Optional: precompute the twiddles as `static const double[]`.**
   Deliberately not done. It would make host and target construct
   bit-identical tables, remove the last dependence on newlib's libm, and
   remove `my_fft_init`'s `malloc` along with the `_Min_Heap_Size` linker
   edit that has to be re-applied after every CubeMX regeneration. The cost
   is ~96 KB of flash at `MYFFT_MAX_LOGN = 12` and reopening code that is
   currently verified and working, to chase a sub-bit effect the on-target
   metric cannot resolve. Judged not worth it; recorded so the reasoning is
   not lost.

The proposal amendment is **not** outstanding: the supervisor is aware of the
Cortex-M4 pivot and no amendment is being pursued. Do not re-raise it.

### Where things physically stand

The board currently holds the verification-plus-benchmark firmware built from
the current tree — ten KAT cases, the fixed twiddles: it runs the KAT check and
the full sweep on a loop, reporting over USART2. Reading it back needs no
rebuild, only a reader attached to `/dev/cu.usbmodem*` at 115200 — or
`st-flash reset` to restart the cycle.

macOS resets the port's termios when the last descriptor closes, so a plain
`stty -f ... 115200` followed by `cat` reads garbage: the settings are gone
before `cat` opens the device. Set the line discipline on the same descriptor
you then read from (`termios.tcsetattr` in a short Python script works).

Everything the write-up needs is already captured in `bench-m4.csv` and
`verify-m4.csv`, and the per-build object sizes come from
`arm-none-eabi-size`, so the board is not needed again unless a new figure is
wanted.

## The `target/` project

Generated by **standalone STM32CubeMX** (`target/target.ioc`) with
Toolchain/IDE set to **Makefile** — not STM32CubeIDE, which locks that
dropdown to itself and cannot emit a Makefile project. Regenerate by opening
the `.ioc`; "Keep User Code when re-generating" is on, so anything inside the
`/* USER CODE BEGIN */` … `/* USER CODE END */` markers survives.

Configuration: HSE bypass 8 MHz (from the ST-Link MCO) → PLL M=4, N=72, P=6 →
**24 MHz SYSCLK**, AHB /1, and `FLASH_LATENCY_0`. Only SYS (Serial Wire) and
USART2 are enabled. Drivers are vendored under `target/Drivers` so the tree
builds without CubeMX installed.

```
cd target && make          # -> build/target.elf, .hex, .bin
```

The study's own code lives in `target/App/`, outside the directories CubeMX
regenerates:

| file | role |
|---|---|
| `cycles.[ch]` | DWT `CYCCNT` timing, plus `cyc_lock`/`cyc_unlock` to mask interrupts across a measured region |
| `console.c` | `__io_putchar` retarget so `printf` reaches USART2 |
| `app.[ch]` | `app_main()`, called from `main()`'s USER CODE 2 marker |
| `kat.[ch]` | **generated** by `tools/gen_kat.c` — do not edit |
| `ref_{emu,native}_{fft,fpr}.c` | four-line wrappers instantiating the unmodified reference under distinct `FALCON_PREFIX` values |
| `refdrv_body.h` | shared driver implementation, included once per backend |
| `ref_{emu,native}_drv.c` | backend selection + `REF_FN` name |
| `myfft_drv.c` | build 3 driver |
| `verify.[ch]` | runs every build against every KAT case |

**All three builds live in one firmware image.** `inner.h` documents
`FALCON_PREFIX` as existing so several versions can coexist in one
application, and there are no non-`Zf()` globals to collide. Each backend
needs its own driver TU because `fpr` is a different type under each
(`uint64_t` vs a struct), so one TU cannot call both. The backend is chosen
by `#define` inside the wrapper `.c` files rather than by compiler flag,
because the CubeMX Makefile compiles every source with identical flags.

`tools/gen_kat.c` regenerates the vectors:

```
clang -O2 -o gen_kat tools/gen_kat.c
./gen_kat target/App/kat.c target/App/kat.h
```

The linker script's `_Min_Heap_Size` is raised from the generated `0x200` to
`0x8000` for `myfft.c`'s twiddle tables, which need 24 560 bytes across
`logn` 1-10 and so will not fit in `0x4000`. Like the Makefile, the linker script
is CubeMX-generated and this change must be re-applied after regeneration.

`target/Makefile` carries a block marked **"study additions -- RE-APPLY AFTER
ANY CubeMX REGENERATION"**. CubeMX rewrites the Makefile wholesale, so that
block and the `OPT` setting must be restored by hand after any regeneration.
`OPT` is already switched from the generated `-Og` to **`-O3`**, to match the
host benchmarks.

Still outstanding: the link uses `-specs=nano.specs`, and newlib-nano's
`printf` drops floating-point conversions unless `-u _printf_float` is added
to `LDFLAGS`. Irrelevant for integer cycle counts, needed as soon as a harness
prints error magnitudes.

## Hardware

Nucleo-F411RE, onboard ST-LINK/V2-1.

```
st-info --probe
```

should report `chipid: 0x431`, `dev-type: STM32F411xC_xE`, 512 KB flash,
128 KB SRAM. UART output arrives over the ST-Link virtual COM port
(`/dev/cu.usbmodem*` on macOS, USART2 on the board side).

No mass-storage volume appears for this board, so flash over SWD:

```
st-flash --connect-under-reset --reset write build/target.bin 0x8000000
```

**`--connect-under-reset` is required.** Without it the write fails with
`Failed to read core_id` even though `st-info --probe` reads the chip fine.
Note also that a running STM32CubeIDE can hold the probe; close it if
connection problems persist.

Reading the board's output does not need a terminal emulator; set the line
discipline and read the character device directly (115200 8N1). The bring-up
firmware repeats its report every 3 s so the port can be attached at any time
without having to catch the boot.

### Verified on target

The DWT counter was validated before any FFT code was flashed: measurement
overhead 1 cycle, `busy(2n)/busy(n)` = 1.999, and `HAL_Delay(100)` measured at
~2.405 Mcycles. That last figure reads ~0.3% above 24 MHz because `HAL_Delay`
adds a tick of margin and starts mid-tick, so it waits 100-101 ms; the result
is consistent with the clock being exactly 24.000 MHz.

