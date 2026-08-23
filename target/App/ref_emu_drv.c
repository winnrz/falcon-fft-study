/*
 * ref_emu_drv.c -- build 1: reference FFT on Falcon's own emulated
 * binary64, which on a Cortex-M4 target is hand-written ARM assembly.
 * Constant-time.  This is the backend a deployed Falcon would use here.
 */

#define FALCON_PREFIX  falcon_emu
#define REF_FN         ref_emu_case
#define BENCH_FN       bench_emu

#include "refdrv_body.h"
#include "benchdrv_body.h"
