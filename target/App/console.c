/*
 * console.c
 *
 * Retargets stdout to USART2, which on the Nucleo-F411RE is wired to
 * the ST-Link's virtual COM port and so appears on the host as
 * /dev/cu.usbmodem* at 115200 8N1.
 *
 * newlib's _write() in Core/Src/syscalls.c is declared weak and calls
 * __io_putchar() per byte, so defining __io_putchar here is enough to
 * capture everything printf emits.
 */

#include <stdio.h>

#include "main.h"
#include "usart.h"

int
__io_putchar(int ch)
{
	uint8_t c;

	c = (uint8_t)ch;
	HAL_UART_Transmit(&huart2, &c, 1, HAL_MAX_DELAY);
	return ch;
}

void
console_init(void)
{
	/*
	 * Unbuffered: output appears as it is produced rather than when
	 * a buffer fills, which matters if the board is reset or halted
	 * partway through a run.
	 */
	setvbuf(stdout, NULL, _IONBF, 0);
}
