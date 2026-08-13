# Falcon FFT polynomial multiplication study.
#
# Targets:
#   make            build the verification harness (default backend)
#   make test       build and run it
#   make test-both  run under both floating-point backends
#   make emu        build with the software (constant-time) binary64 backend
#   make clean      remove build products

.POSIX:

CC      = clang
CFLAGS  = -Wall -Wextra -Wshadow -Wundef -O3
LD      = clang
LDFLAGS =
LIBS    = -lm

OBJ = fft.o fpr.o

all: test_fftmul

test_fftmul: test_fftmul.o $(OBJ)
	$(LD) $(LDFLAGS) -o test_fftmul test_fftmul.o $(OBJ) $(LIBS)

# Software binary64 build, in a separate directory so the two backends
# never share object files.
emu:
	mkdir -p build-emu
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c -o build-emu/fft.o fft.c
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c -o build-emu/fpr.o fpr.c
	$(CC) $(CFLAGS) -DFALCON_FPEMU=1 -I. -c -o build-emu/test_fftmul.o test_fftmul.c
	$(LD) $(LDFLAGS) -o test_fftmul_emu build-emu/test_fftmul.o \
		build-emu/fft.o build-emu/fpr.o $(LIBS)

test: test_fftmul
	./test_fftmul

test-both: test_fftmul emu
	./test_fftmul
	@echo
	./test_fftmul_emu

clean:
	-rm -rf $(OBJ) test_fftmul test_fftmul.o test_fftmul_emu build-emu

fft.o: fft.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o fft.o fft.c

fpr.o: fpr.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o fpr.o fpr.c

test_fftmul.o: test_fftmul.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) -c -o test_fftmul.o test_fftmul.c
