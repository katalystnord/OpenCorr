#!/usr/bin/env bash
# Run the smoke tests under the address and undefined-behaviour sanitizers.
#
# ⚑ WHAT THIS CATCHES THAT NOTHING ELSE DOES. A read or write one element past
# the end of a vector returns whatever sits there and changes nothing
# downstream: no assertion can see it, and a mutation sweep reports it as a
# survivor with nothing to be done about it. That is the largest single class of
# survivor left in SurView and PatternFab after every case worth writing had
# been written, and it is the reason this exists.
#
# It matters more here than in either of those. This library is index and
# pointer arithmetic in OpenMP hot loops, and a defect in it is a WRONG NUMBER
# rather than a wrong message.
#
# Its own build directory: the instrumentation changes the object code, so
# producing the library anyone links and checking it are different jobs.
set -euo pipefail

cd "$(dirname "$0")/.."

cmake -S . -B build-asan -G Ninja \
    -DOPENCORR_BUILD_SMOKE_TEST=ON \
    -DOPENCORR_SANITIZE=ON
cmake --build build-asan

# -fno-sanitize-recover is compiled in, so a diagnostic fails the test rather
# than printing and carrying on.
ctest --test-dir build-asan --output-on-failure "$@"
