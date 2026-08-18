#!/bin/sh
# coverage.sh -- what the smoke tests actually reach in this library.
#
# Coverage is a map of the UNTESTED, not a score to raise. A line the suite
# never executes cannot be protected by it, however green the run looks.
#
# This fork carries capability work that upstream does not have -- calibration,
# cine I/O, uncertainty, the status-flag taxonomy, auto-ROI, crack residual --
# and each arrived with a smoke test written alongside it. What no one had
# checked until now is how much of that code those tests actually touch.
#
# Builds into its own directory, because --coverage changes the object code and
# sharing a build tree means never being sure which library you have.
#
# Usage: tools/coverage.sh [--html]

set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$here/build-coverage"
out="$build/report"

want_html=0
want_full=0
for arg in "$@"; do
    case "$arg" in
        --html) want_html=1 ;;
        --full) want_full=1 ;;
    esac
done

# Two of the twenty tests are multi-minute numerical runs, and instrumenting
# them is disproportionately expensive: gcov counters are shared across threads,
# so an instrumented parallel run contends on cache lines rather than computing.
# Measured: uncertainty_smoke_test takes 225 s optimised and had not finished
# after 47 MINUTES instrumented; sequence_tracker_smoke_test takes 186 s and had
# not finished after 26. Serialising does not rescue them -- OpenCV keeps its own
# threads -- so they are skipped by default and the report says which source
# files that leaves unmeasured.
#
# --full includes them. Expect it to take well over an hour.
SLOW_TESTS='uncertainty_smoke_test|sequence_tracker_smoke_test'
SLOW_SOURCES='src/oc_uncertainty.cpp and src/oc_sequence_tracker.cpp'


if ! command -v gcovr >/dev/null 2>&1; then
    echo "coverage.sh: gcovr not found. Install it with: sudo apt install gcovr" >&2
    exit 1
fi

echo "=== Building instrumented ==="
cmake -S "$here" -B "$build" -G Ninja \
    -DOPENCORR_BUILD_SMOKE_TEST=ON -DOPENCORR_COVERAGE=ON >/dev/null
cmake --build "$build"

echo ""
echo "=== Running the smoke tests, single-threaded ==="
# OMP_NUM_THREADS=1 is not a shortcut here, it is what makes this finish.
#
# gcov counters are shared, and every OpenMP thread increments the counters for
# the same lines, so an instrumented parallel run spends its time contending on
# cache lines rather than computing. Measured: uncertainty_smoke_test takes
# 225 s optimised, and had not finished after 47 MINUTES instrumented at ~1000%
# CPU, on a test that is a few minutes of actual work.
#
# Serialising removes the contention. Coverage does not care how many threads
# executed a line, only that one did, so nothing is lost from the measurement.
#
# Run from the repository root: several tests read example images by a path
# relative to it.
if [ "$want_full" -eq 1 ]; then
    ( cd "$here" && OMP_NUM_THREADS=1 ctest --test-dir "$build" --output-on-failure )
else
    echo "Skipping the two heaviest tests (pass --full to include them):"
    echo "  uncertainty_smoke_test, sequence_tracker_smoke_test"
    echo "So $SLOW_SOURCES will read as uncovered here even though they"
    echo "are tested. That is a limit of this run, not a gap in the suite."
    echo ""
    ( cd "$here" && OMP_NUM_THREADS=1 ctest --test-dir "$build" \
        --output-on-failure -E "$SLOW_TESTS" )
fi

mkdir -p "$out"

echo ""
echo "=== Coverage: src/ ==="
# The smoke tests themselves are excluded -- measuring how thoroughly a test
# runs itself flatters the total and says nothing.
gcovr \
    --root "$here" \
    --filter "$here/src/" \
    --exclude '.*_smoke_test\.cpp' \
    --print-summary \
    --sort=uncovered-percent \
    --txt "$out/coverage.txt" \
    --json-summary "$out/coverage.json" \
    "$build"

if [ "$want_html" -eq 1 ]; then
    gcovr --root "$here" --filter "$here/src/" --exclude '.*_smoke_test\.cpp' \
        --html-details "$out/index.html" "$build" >/dev/null
    echo ""
    echo "HTML report: $out/index.html"
fi

echo ""
if [ "$want_full" -eq 0 ]; then
    echo "NOTE: $SLOW_SOURCES are not measured in this run."
    echo "      Their own tests were skipped for time. Use --full to include them."
fi
echo "Per-file detail: $out/coverage.txt"
