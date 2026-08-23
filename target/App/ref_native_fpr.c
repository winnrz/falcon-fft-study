/*
 * ref_native_fpr.c -- generated-by-hand wrapper, not a copy.
 *
 * Instantiates the unmodified reference fpr.c under its own symbol
 * prefix so both float backends can coexist in one firmware and be
 * measured in the same run.  inner.h documents FALCON_PREFIX as existing
 * for precisely this purpose.
 *
 * The backend is selected here rather than on the compiler command line
 * because the CubeMX Makefile compiles every source with identical
 * flags.
 */

#define FALCON_PREFIX     falcon_native
#define FALCON_FPNATIVE   1

#include "../../fpr.c"
