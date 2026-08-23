/*
 * verify.h -- on-target correctness and cost for all three builds.
 */

#ifndef VERIFY_H__
#define VERIFY_H__

/* Runs every build against every KAT case.  Returns 0 if all passed. */
int verify_run(void);

#endif
