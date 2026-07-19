# Third-party licenses

This fork (`katalystnord/OpenCorr`) ports algorithms and data-model designs
from two other BSD-3-Clause-licensed open-source projects, in addition to
OpenCorr's own MPL-2.0 license (see `LICENSE`). Per BSD-3-Clause, this file
carries the required copyright notice, list of conditions, and disclaimer
for each. Individual source files with ported content carry a short pointer
comment back to this file rather than repeating the full text in every
location; see each such file's own header comment for exactly what was
ported from where.

A third project, `dicengine/hypercine`, is vendored verbatim (not ported)
under `deps/hypercine/` -- its BSD-3-Clause notice is already carried
in full in every file there and in `deps/hypercine/LICENSE`, so it is not
repeated below.

## DICe (dicengine/dice)

Ported into: `oc_shape.h`/`.cpp`, `oc_simplex.h`/`.cpp`,
`oc_phase_correlation.h`/`.cpp`, `oc_camera_calibrator.h`/`.cpp`,
`oc_subset.h`/`.cpp`, `oc_uncertainty.h`/`.cpp`, `oc_sequence_tracker.h`/`.cpp`.

```
              Digital Image Correlation Engine (DICe)
                Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## ncorr_2D_cpp (justinblaber/ncorr_2D_cpp)

Ported into: `oc_reliability_guided.h`/`.cpp`, `oc_sequence_tracker.h`/`.cpp`.

```
Copyright (c) 2014, Justin Blaber
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution
    * Neither the name of the Georgia Institute of Technology nor the names
      of its contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```
