#!/bin/sh
#
# measure_memory.sh
#
# Memory measurement for the FFT polynomial multiplication study.
#
# Three separate questions, measured by three different means:
#
#   1. Peak resident set size of the whole process.  Cheap, but coarse:
#      it includes the C runtime, the test harness's own buffers, and
#      the static twiddle tables, so it overstates what the transform
#      itself needs.
#
#   2. Heap profile over time (valgrind/massif).  Attributes allocations
#      to call sites, which separates the harness's buffers from the
#      implementation's.  Linux only in practice.
#
#   3. Leak check (valgrind/memcheck).  Confirms my_fft_free() releases
#      every twiddle table and the harnesses release their buffers.
#
# The figure that actually characterises the algorithm -- bytes of
# working memory for one transform -- is not any of these; it is the
# static accounting printed by ./bench_fftmul, because both transforms
# operate in place on caller-provided buffers and allocate nothing.
# These measurements exist to confirm that claim from the outside.
#
# Usage:  ./measure_memory.sh [binary ...]
#         defaults to the three harnesses.

set -u

BINS="${*:-./test_fftmul ./test_myfft ./bench_fftmul}"
OUT=memory-results
mkdir -p "$OUT"

uname_s=$(uname -s)

echo "Memory measurement"
echo "  platform  $uname_s"
echo "  binaries  $BINS"
echo "  output    $OUT/"
echo

# --------------------------------------------------------------------
# 1. Peak resident set size
# --------------------------------------------------------------------

echo "=== Peak resident set size ==="
for b in $BINS; do
	if [ ! -x "$b" ]; then
		echo "  skip $b (not built)"
		continue
	fi
	case "$uname_s" in
	Linux)
		# GNU time reports peak RSS in kilobytes.
		if command -v /usr/bin/time >/dev/null 2>&1; then
			kb=$(/usr/bin/time -f '%M' "$b" 2>&1 >/dev/null \
				| tail -1)
			echo "  $b: ${kb} KB peak RSS"
		else
			echo "  $b: /usr/bin/time not available"
		fi
		;;
	Darwin)
		# BSD time reports maximum resident set size in bytes.
		bytes=$(/usr/bin/time -l "$b" 2>&1 >/dev/null \
			| awk '/maximum resident set size/ {print $1}')
		if [ -n "$bytes" ]; then
			echo "  $b: $((bytes / 1024)) KB peak RSS"
		else
			echo "  $b: could not parse /usr/bin/time -l"
		fi
		;;
	*)
		echo "  $b: unsupported platform for RSS measurement"
		;;
	esac
done
echo

# --------------------------------------------------------------------
# 2. Heap profile
# --------------------------------------------------------------------

echo "=== Heap profile (massif) ==="
if command -v valgrind >/dev/null 2>&1; then
	for b in $BINS; do
		if [ ! -x "$b" ]; then
			continue
		fi
		name=$(basename "$b")
		valgrind --tool=massif --massif-out-file="$OUT/massif.$name" \
			--time-unit=B "$b" >/dev/null 2>&1
		if [ -f "$OUT/massif.$name" ]; then
			peak=$(awk -F= '/^mem_heap_B=/ {if ($2+0 > m) m = $2+0}
				END {print m}' "$OUT/massif.$name")
			echo "  $name: $peak bytes peak heap" \
				"-> $OUT/massif.$name"
			if command -v ms_print >/dev/null 2>&1; then
				ms_print "$OUT/massif.$name" \
					> "$OUT/massif.$name.txt" 2>/dev/null
				echo "      readable report:" \
					"$OUT/massif.$name.txt"
			fi
		fi
	done
else
	echo "  valgrind not installed."
	case "$uname_s" in
	Linux)  echo "  install with: sudo apt install valgrind" ;;
	Darwin) echo "  valgrind is unreliable on Apple silicon;" \
			"run this step on the Ubuntu machine." ;;
	esac
fi
echo

# --------------------------------------------------------------------
# 3. Leak check
# --------------------------------------------------------------------

echo "=== Leak check (memcheck) ==="
if command -v valgrind >/dev/null 2>&1; then
	for b in $BINS; do
		if [ ! -x "$b" ]; then
			continue
		fi
		name=$(basename "$b")
		valgrind --leak-check=full --errors-for-leak-kinds=definite \
			--error-exitcode=9 "$b" \
			> "$OUT/memcheck.$name.txt" 2>&1
		rc=$?
		if [ "$rc" -eq 9 ]; then
			echo "  $name: LEAKS FOUND -> $OUT/memcheck.$name.txt"
		else
			echo "  $name: clean"
		fi
	done
else
	echo "  valgrind not installed; skipped."
fi
echo

echo "Static accounting of the transform's own working set:"
echo "  ./bench_fftmul  (last section of its output)"
