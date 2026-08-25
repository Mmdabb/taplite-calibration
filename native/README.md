# Bundled TAPLite kernel

This directory is the minimum native source set required by
`taplite-calibration`. It was migrated from the TAPLite4MPO
`dev/auto-calibration` working tree on 2026-08-24, including the latest refined
volume-envelope implementation. This repository is now the authoritative home
for the packaged calibration implementation. The packaged default fit mode is
`refined_fixed_point`; legacy and further-development fit modes require an
explicit native settings override.

Included here:

- `TAPLite.cpp` and `TAPLite.h`: assignment kernel and dedicated
  `AutoCalibrationAPI` integration point;
- `AutoCalibration*`: refined calibration engine, inverse/oracle logic, and the
  isolated further-development helper still referenced by the engine;
- `binding.cpp`: Python extension surface for ordinary assignment and dedicated
  auto calibration;
- focused native tests and the original TAPLite MIT license.

Not migrated: native binaries, old build trees, CBI packages, NVTA inputs,
historical outputs, unrelated TAPLite examples, or executable launchers.

Normal users build through `pip`; `setup.py` compiles this source into
`taplite_calibration._native`. The CMake file is provided for focused native
development and C++ unit tests.
