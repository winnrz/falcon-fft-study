/*
 * app.c
 *
 * Bring-up firmware for the Cortex-M4 port.  It verifies the two pieces
 * every later measurement depends on, before any FFT code is involved:
 *
 *   1. stdout reaches the host over USART2.
 *   2. The DWT cycle counter runs, and counts what it should.
 *
 * Check 2 is the interesting one.  A cycle counter that is merely
 * incrementing proves very little -- it could be counting at the wrong
 * rate, and every timing taken afterwards would be silently wrong by a
 * constant factor.  Two independent checks pin it down:
 *
 *   Linearity   twice the work should cost twice the cycles.  This
 *               catches a counter that is running but not tracking
 *               work, and it needs no knowledge of the clock rate.
 *
 *   Clock rate  HAL_Delay(100) is timed by SysTick, which is derived
 *               from HCLK independently of DWT.  Measuring it with
 *               CYCCNT and dividing by the elapsed 0.1 s recovers the
 *               core clock.  It should come out at 24 MHz -- the value
 *               configured in CubeMX.  If the PLL were misconfigured,
 *               or HSE had failed over to the 16 MHz HSI, this line
 *               would say so immediately.
 */

#include <stdio.h>
#include <stdint.h>

#include "main.h"
#include "app.h"
#include "cycles.h"
#include "verify.h"

/*
 * Deliberately unoptimisable work.  The volatile forces the loop to
 * survive -O3, so the cycle count scales with n.
 */
static void
busy(uint32_t n)
{
	volatile uint32_t k = 0;

	while (n--) {
		k++;
	}
}

static uint32_t
time_busy(uint32_t n, uint32_t overhead)
{
	uint32_t prim, t0, t1, d;

	prim = cyc_lock();
	t0 = cyc_read();
	busy(n);
	t1 = cyc_read();
	cyc_unlock(prim);

	d = t1 - t0;
	return d > overhead ? d - overhead : 0;
}

static void
report(void)
{
	uint32_t overhead, c1, c2, t0, t1, delay_cycles, implied_hz;

	printf("\r\n");
	printf("Falcon FFT study -- Cortex-M4 bring-up\r\n");
	printf("======================================\r\n");

	printf("  SystemCoreClock   %lu Hz\r\n",
		(unsigned long)SystemCoreClock);

	overhead = cyc_overhead();
	printf("  measurement cost  %lu cycles\r\n",
		(unsigned long)overhead);

	/* Linearity: 2x the iterations should be ~2x the cycles. */
	c1 = time_busy(1000, overhead);
	c2 = time_busy(2000, overhead);
	printf("  busy(1000)        %lu cycles\r\n", (unsigned long)c1);
	printf("  busy(2000)        %lu cycles\r\n", (unsigned long)c2);
	printf("  ratio x1000       %lu  (expect ~2000)\r\n",
		(unsigned long)(c1 ? (c2 * 1000u) / c1 : 0));

	/*
	 * Clock rate.  Interrupts stay enabled here: HAL_Delay is driven
	 * by the SysTick interrupt and would hang without them.
	 */
	t0 = cyc_read();
	HAL_Delay(100);
	t1 = cyc_read();
	delay_cycles = t1 - t0;
	implied_hz = delay_cycles * 10u;    /* 100 ms -> 1 s */

	printf("  100 ms measured   %lu cycles\r\n",
		(unsigned long)delay_cycles);
	printf("  implied clock     %lu Hz  (expect 24000000)\r\n",
		(unsigned long)implied_hz);

	printf("\r\n");
}

void
app_main(void)
{
	console_init();
	cyc_init();

	/*
	 * Repeated rather than printed once, so the report can be read
	 * by attaching to the port at any time without having to catch
	 * the boot.
	 */
	for (;;) {
		report();
		verify_reference();
		printf("\r\n");
		HAL_Delay(5000);
	}
}
