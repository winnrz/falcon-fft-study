/*
 * ref_native_drv.c -- build 2: the same reference FFT on the C double
 * type, which the M4F's single-precision FPU cannot execute, so GCC
 * lowers every operation to a libgcc soft-double call.  Faster than the
 * emulated backend but not constant-time.
 *
 * Build 2 exists to be compared against build 3: myfft.c also runs on
 * libgcc soft-double, so the two share a float backend and any
 * difference between them is FFT structure alone.
 */

#define FALCON_PREFIX     falcon_native
#define FALCON_FPNATIVE   1
#define REF_FN            ref_native_case
#define BENCH_FN          bench_native

#include "refdrv_body.h"
#include "benchdrv_body.h"
