/*
 * app.h
 *
 * Entry point for the study's own code.  main() is CubeMX-generated and
 * is left alone apart from a call to app_main() inside its USER CODE
 * markers; everything below that lives here, outside the directories
 * CubeMX regenerates.
 */

#ifndef APP_H__
#define APP_H__

/* Runs the bring-up report; does not return. */
void app_main(void);

/* From console.c. */
void console_init(void);

#endif
