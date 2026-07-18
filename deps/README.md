# Vendored dependencies

- `nanoflann/nanoflann.hpp` — [jlblancoc/nanoflann](https://github.com/jlblancoc/nanoflann) v1.7.0, BSD license. Header-only, vendored directly since OpenCorr has no upstream CMake dependency-fetch mechanism. Matches the version pinned in SurView's `CLAUDE.md`.
- `hypercine/` — [dicengine/hypercine](https://github.com/dicengine/hypercine), BSD-3-Clause (NTESS). Reads `.cine` files (Phantom/Vision Research high-speed camera format). Vendored as source (not header-only); `storage_t` fixed to `float` via `USE_FLOAT_STORAGE` (set in `CMakeLists.txt`) to match `Image2D`'s `Eigen::MatrixXf` representation. Wrapped by `src/oc_cine.h/.cpp` (`Cine2D`), see [issue #4](https://github.com/katalystnord/OpenCorr/issues/4).
