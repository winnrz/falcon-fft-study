/*
 * bench.c
 *
 * Sweeps every build across every size and reports per-operation cycle
 * costs.
 *
 * Output is a human-readable table followed by CSV lines, each prefixed
 * with "CSV," so they can be lifted straight out of a captured serial
 * log with grep.  The CSV columns match bench.csv on the host, with
 * cycles in place of nanoseconds, so the two can be set side by side.
 *
 * The reference caps out at logn = 10, which is Falcon-1024; logn = 9 is
 * Falcon-512 and is the row that matters most.
 */

#include <stdio.h>
#include <stdint.h>

#include "bench.h"
#include "benchbuf.h"

struct bench_build {
	const char *name;
	const char *tag;
	void      (*run)(unsigned, struct bench_row *);
};

static const struct bench_build blds[] = {
	{ "1 ref  emu+asm",  "ref_emu",    bench_emu    },
	{ "2 ref  softdbl",  "ref_native", bench_native },
	{ "3 myfft softdbl", "myfft",      bench_myfft  },
};

#define NB  ((unsigned)(sizeof blds / sizeof blds[0]))

/*
 * Results are collected once and printed twice.  Running the sweep again
 * for the CSV would double the runtime and, worse, let the two renderings
 * disagree if anything varied between them.
 */
static struct bench_row rows[NB][BENCH_MAX_LOGN + 1];

void
bench_run(void)
{
	unsigned b, logn;

	for (b = 0; b < NB; b++) {
		for (logn = 1; logn <= BENCH_MAX_LOGN; logn++) {
			blds[b].run(logn, &rows[b][logn]);
		}
	}

	printf("\r\nPer-operation cost in cycles, minimum of 3 repeats\r\n");
	printf("  24 MHz, zero flash wait states, interrupts masked\r\n\r\n");
	printf("  %-16s %4s %6s %9s %9s %8s %10s %8s\r\n",
		"build", "logn", "n", "fft", "ifft", "mul", "full",
		"spread");
	printf("  ------------------------------------------"
		"---------------------------------\r\n");

	for (b = 0; b < NB; b++) {
		for (logn = 1; logn <= BENCH_MAX_LOGN; logn++) {
			const struct bench_row *r = &rows[b][logn];

			printf("  %-16s %4u %6lu %9lu %9lu %8lu %10lu %8lu\r\n",
				logn == 1 ? blds[b].name : "",
				logn, (unsigned long)1UL << logn,
				(unsigned long)r->fft,
				(unsigned long)r->ifft,
				(unsigned long)r->mul,
				(unsigned long)r->full,
				(unsigned long)r->spread);
		}
	}

	printf("  ------------------------------------------"
		"---------------------------------\r\n");

	printf("\r\nCSV,build,logn,n,fft,ifft,mul,full,restore,spread\r\n");
	for (b = 0; b < NB; b++) {
		for (logn = 1; logn <= BENCH_MAX_LOGN; logn++) {
			const struct bench_row *r = &rows[b][logn];

			printf("CSV,%s,%u,%lu,%lu,%lu,%lu,%lu,%lu,%lu\r\n",
				blds[b].tag, logn,
				(unsigned long)1UL << logn,
				(unsigned long)r->fft,
				(unsigned long)r->ifft,
				(unsigned long)r->mul,
				(unsigned long)r->full,
				(unsigned long)r->restore,
				(unsigned long)r->spread);
		}
	}
	printf("\r\n");
}
