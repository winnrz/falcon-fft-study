# Falcon FFT polynomial multiplication study.
#
# Targets:
#   make            build every harness
#   make test       verify the Falcon reference implementation
#   make test-mine  verify the from-scratch implementation, and
#                   cross-check it against the reference
#   make test-all   run both verification harnesses
#   make test-both  run the reference harness under both FP backends
#   make bench      measure performance
#   make bench-csv  same, as CSV in bench.csv, for plotting
#   make memory     measure memory (peak RSS, massif, leak check)
#   make emu        build with the software binary64 backend
#   make clean      remove build products

.POSIX:

CC      = clang
CFLAGS  = -Wall -Wextra -Wshadow -Wundef -O3
LD      = clang
LDFLAGS =
LIBS    = -lm

REF_OBJ  = fft.o fpr.o
MINE_OBJ = myfft.o

all: test_fftmul test_myfft bench_fftmul

# --- the three harnesses ---------------------------------------------

test_fftmul: test_fftmul.o $(REF_OBJ)
	$(LD) $(LDFLAGS) -o test_fftmul test_fftmul.o $(REF_OBJ) $(LIBS)

test_myfft: test_myfft.o $(MINE_OBJ) $(REF_OBJ)
	$(LD) $(LDFLAGS) -o test_myfft test_myfft.o $(MINE_OBJ) \
		$(REF_OBJ) $(LIBS)

bench_fftmul: bench_fftmul.o $(MINE_OBJ) $(REF_OBJ)
	$(LD) $(LDFLAGS) -o bench_fftmul bench_fftmul.o $(MINE_OBJ) \
		$(REF_OBJ) $(LIBS)

# --- software binary64 build -----------------------------------------
#
# Built in a separate directory so the two backends never share object
# files.  Note that only the reference is affected: myfft.c always uses
# the hardware double type, so test_myfft_emu compares a hardware-double
# implementation against a software-binary64 reference.

emu:
	mkdir -p build-emu
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c -o build-emu/fft.o fft.c
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c -o build-emu/fpr.o fpr.c
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c \
		-o build-emu/test_fftmul.o test_fftmul.c
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c \
		-o build-emu/test_myfft.o test_myfft.c
	$(CC) $(CFLAGS) -I. -c -o build-emu/myfft.o myfft.c
	$(LD) $(LDFLAGS) -o test_fftmul_emu build-emu/test_fftmul.o \
		build-emu/fft.o build-emu/fpr.o $(LIBS)
	$(LD) $(LDFLAGS) -o test_myfft_emu build-emu/test_myfft.o \
		build-emu/myfft.o build-emu/fft.o build-emu/fpr.o $(LIBS)

# --- running ---------------------------------------------------------

test: test_fftmul
	./test_fftmul

test-mine: test_myfft
	./test_myfft

test-all: test_fftmul test_myfft
	./test_fftmul
	@echo
	./test_myfft

test-both: test_fftmul emu
	./test_fftmul
	@echo
	./test_fftmul_emu

bench: bench_fftmul
	./bench_fftmul

bench-csv: bench_fftmul
	./bench_fftmul --csv > bench.csv
	@echo "wrote bench.csv"

memory: test_fftmul test_myfft bench_fftmul
	./measure_memory.sh

clean:
	-rm -rf $(REF_OBJ) $(MINE_OBJ) \
		test_fftmul test_fftmul.o test_fftmul_emu \
		test_myfft test_myfft.o test_myfft_emu \
		bench_fftmul bench_fftmul.o \
		memory-results build-emu

# --- objects ---------------------------------------------------------

fft.o: fft.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o fft.o fft.c

fpr.o: fpr.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o fpr.o fpr.c

myfft.o: myfft.c myfft.h
	$(CC) $(CFLAGS) -c -o myfft.o myfft.c

test_fftmul.o: test_fftmul.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o test_fftmul.o test_fftmul.c

test_myfft.o: test_myfft.c config.h inner.h fpr.h myfft.h
	$(CC) $(CFLAGS) -c -o test_myfft.o test_myfft.c

bench_fftmul.o: bench_fftmul.c config.h inner.h fpr.h myfft.h
	$(CC) $(CFLAGS) -c -o bench_fftmul.o bench_fftmul.c
