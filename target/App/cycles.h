/*
 * cycles.h
 *
 * Cycle-accurate timing for the Cortex-M4, via the DWT unit's CYCCNT
 * register.  CYCCNT increments once per core clock, so it measures the
 * quantity the comparable literature reports and is independent of what
 * the PLL happens to be set to.
 *
 * The counter is 32 bits and wraps every 2^32 cycles -- about three
 * minutes at 24 MHz.  Unsigned subtraction of two reads is correct
 * across a single wrap, so no interval shorter than that needs care.
 *
 * Timed regions should run with interrupts masked.  The HAL enables
 * SysTick at 1 kHz, and a tick landing inside a measured region adds
 * its handler's cycles to the result; the effect is small but it is
 * pure noise in a measurement whose method is to take the minimum.
 * Use cyc_lock()/cyc_unlock() around anything being measured, except
 * where the code under test needs interrupts itself (HAL_Delay does).
 */

#ifndef CYCLES_H__
#define CYCLES_H__

#include <stdint.h>

#include "stm32f4xx.h"

/* Enable the trace unit and start CYCCNT.  Call once at startup. */
void cyc_init(void);

/*
 * Cost of the measurement itself: the minimum observed delta between
 * two back-to-back reads.  Subtract it from measured intervals.
 */
uint32_t cyc_overhead(void);

static inline uint32_t
cyc_read(void)
{
	return DWT->CYCCNT;
}

static inline uint32_t
cyc_lock(void)
{
	uint32_t prim;

	prim = __get_PRIMASK();
	__disable_irq();
	return prim;
}

static inline void
cyc_unlock(uint32_t prim)
{
	__set_PRIMASK(prim);
}

#endif
