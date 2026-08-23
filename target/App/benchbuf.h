/*
 * benchbuf.h
 *
 * Shared working buffers.
 *
 * Each backend needs its own buffers in principle -- fpr is a uint64_t
 * of raw bits under FALCON_FPEMU, a struct wrapping a double under
 * FALCON_FPNATIVE, and myfft.c uses plain double -- but all three are
 * eight bytes wide and the builds run one after another, never at the
 * same time.  Giving each its own set at n = 1024 would cost 96 KB of
 * the F411's 128 KB; sharing costs 32 KB.
 *
 * The storage is handed out through a function returning void *, defined
 * in its own translation unit, so the compiler cannot see the underlying
 * declared type at the point of use and the cast is not an aliasing
 * violation it can act on.  There is no link-time optimisation in this
 * build to defeat that.
 */

#ifndef BENCHBUF_H__
#define BENCHBUF_H__

#define BENCH_MAX_LOGN   10
#define BENCH_MAX_N      (1u << BENCH_MAX_LOGN)
#define BENCH_SLOTS      4

/* Slot i, 8-byte aligned, BENCH_MAX_N elements of 8 bytes. */
void *bench_slot(unsigned i);

#endif
