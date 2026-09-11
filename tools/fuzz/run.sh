#!/usr/bin/env bash
# Fuzz the .cine header parser: arbitrary bytes in, and the sanitizers decide.
#
# ⚑ WHY THIS FORMAT. A .cine is a high-speed camera's own binary format, and it
# is the one place in this library where a stranger's bytes decide how much
# memory is allocated and how far a pointer walks. Everything else is read
# through OpenCV or written by us.
#
# Needs clang, which ships libFuzzer; nothing else to install. The corpus is the
# example files that ship, so the fuzzer starts from valid structure and mutates
# outwards rather than guessing a 44-byte header from noise.
#
# Findings land in tools/fuzz/findings/ as the exact bytes that caused them.
# Re-run one with:  ./build/cine_header_fuzz tools/fuzz/findings/<file>
set -euo pipefail

cd "$(dirname "$0")/../.."

SECONDS_TO_RUN="${1:-300}"
BUILD=/tmp/cine_header_fuzz
CORPUS=$(mktemp -d)
FINDINGS=tools/fuzz/findings

mkdir -p "$FINDINGS"
cp examples/cine/*.cine "$CORPUS"/

clang++ -std=c++17 -g -O1 -DUSE_FLOAT_STORAGE=1 \
    -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=all \
    -I deps/hypercine \
    -o "$BUILD" tools/fuzz/cine_header_fuzz.cpp deps/hypercine/hypercine.cpp

# A refusal is the correct outcome for garbage and the parser says so on stdout,
# loudly and often; only the sanitizers' verdict matters here.
"$BUILD" "$CORPUS" \
    -artifact_prefix="$FINDINGS/" \
    -max_total_time="$SECONDS_TO_RUN" \
    -print_final_stats=1
