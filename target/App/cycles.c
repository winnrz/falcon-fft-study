/*
 * cycles.c -- see cycles.h.
 */

#include "cycles.h"

void
cyc_init(void)
{
	/*
	 * TRCENA gates the whole trace block, DWT included.  It can be
	 * set by software with no debugger attached, so the counter
	 * works in a standalone binary.
	 */
	CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;

	DWT->CYCCNT = 0;
	DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}

uint32_t
cyc_overhead(void)
{
	uint32_t best, i;

	best = 0xFFFFFFFFu;
	for (i = 0; i < 64; i++) {
		uint32_t prim, t0, t1, d;

		prim = cyc_lock();
		t0 = cyc_read();
		t1 = cyc_read();
		cyc_unlock(prim);

		d = t1 - t0;
		if (d < best) {
			best = d;
		}
	}
	return best;
}
