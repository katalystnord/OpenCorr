// Feed arbitrary bytes to the .cine header parser and see whether it survives.
//
// ⚑ WHY THIS FILE AND NOT ANOTHER. A .cine is a high-speed camera's own binary
// format, read field by field out of a file whose contents nobody here controls
// -- the one place in this library where a stranger's bytes decide how much
// memory is allocated and how far a pointer walks. Fork issue #18 found three
// crash bugs in read_header() by READING it: a division by zero, an
// out-of-bounds read, and a truncation check that overflows 32-bit. All three
// are the shape a fuzzer finds in an afternoon without anyone suspecting them
// first, which is the argument for having one.
//
// A refusal is the CORRECT outcome for garbage, so exceptions are caught and
// discarded here. What must not happen is a crash, an out-of-bounds access, or
// an allocation sized from a number the file simply asserted. The sanitizers
// are what decide that, not this file.
//
// The parser takes a path rather than a buffer, so each input is written to a
// file first. Keep it on tmpfs (the default here) or the run is bounded by the
// disk rather than by the code under test.
//
// Build and run with tools/fuzz/run.sh.

#include "hypercine.h"

#include <cstdint>
#include <cstdio>
#include <fstream>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    static const char *path = "/tmp/hypercine_fuzz_input.cine";

    {
        std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::trunc);
        if (!file.is_open()) {
            return 0;
        }
        file.write(reinterpret_cast<const char *>(data), static_cast<std::streamsize>(size));
    }

    try {
        hypercine::HyperCine cine(path);
    } catch (...) {
        // Refusing a malformed file is the right answer, and the overwhelming
        // majority of inputs are malformed. Only a fault reaches the sanitizer.
    }
    return 0;
}
