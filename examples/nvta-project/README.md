# NVTA-style project wiring

Copy `calibration.toml` to a data project created with
`python -m taplite_calibration init`. Keep the network and observations outside
the Python package repository. The installed platform wheel supplies the native
TAPLite kernel; no neighboring TAPLite4MPO checkout is required.

Remove `fallback_qvdf_dictionary` when all directed node pairs are present in
AM, MD, and PM and no extra base-dictionary coverage needs to be retained.

The `[prepare.*]` paths are fallback sources. A run first validates the normal
`inputs/` paths. When they pass, those raw paths are not opened. When they fail,
the package warns, creates a run-local prepared copy, validates it again, and
only then starts ODME or auto calibration. Preparation is always active and is
not a selectable mode.
