# NVTA-style project wiring

Copy `calibration.toml` to a data project created with
`python -m taplite_calibration init`. Keep the network and observations outside
the Python package repository. If this data project sits beside a TAPLite4MPO
development checkout, the example's `../TAPLite4MPO` path remains portable.

Remove `fallback_qvdf_dictionary` when all directed node pairs are present in
AM, MD, and PM and no extra base-dictionary coverage needs to be retained.
Remove `pytaplite_path` after installing a compatible TAPLite4MPO wheel.

The `[prepare.*]` paths are fallback sources. A run first validates the normal
`inputs/` paths. When they pass, those raw paths are not opened. When they fail,
the package warns, creates a run-local prepared copy, validates it again, and
only then starts ODME or auto calibration. Preparation is always active and is
not a selectable mode.
