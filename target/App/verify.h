/*
 * verify.h -- on-target correctness check for the Falcon reference FFT.
 */

#ifndef VERIFY_H__
#define VERIFY_H__

/* Runs every KAT case, prints a table.  Returns 0 if all passed. */
int verify_reference(void);

#endif
