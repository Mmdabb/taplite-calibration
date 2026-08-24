# Changelog

## 0.2.0 - 2026-08-24

- Added an automatic prepared-input QA gate before every ODME or
  auto-calibration run. Preparation is always active and is not a public mode.
- Added daily screen-count normalization with configurable source columns.
- Added portable OD policy assembly and optional streaming policy generation
  from DTAC-v2 route pools with integer path-screen link multiplicity.
- Added CBI average-weekday/candidate-episode preprocessing, period slicing,
  canonical TMC-link selection, managed-lane exclusion, GP link targets, and
  calibration-ready scenario enrichment.
- Added explicit warning/repair/re-QA behavior and fail-fast errors when raw
  sources cannot repair an invalid prepared input.
- Added a checked-in, synthetic three-zone quickstart dataset and executable
  end-to-end example that requires no proprietary data or TAPLite executable.
- Hardened OD demand and settings updates for pandas 3 read-only arrays and
  strict mixed-type assignment behavior.

## 0.1.0 - 2026-08-24

- Added joint NT/AM/MD/PM screen-count ODME with cell bounds, KL regularization,
  optional per-period total-variation budgets, and exact period/mode marginals.
- Added TAPLite refined QVDF auto-calibration orchestration through the native
  `pytaplite.auto_calibrate` API.
- Added `odme`, `auto-calibration`, and ordered `both` modes.
- Added portable OD-factor and calibrated QVDF dictionaries, run manifests,
  logs, validation, and a deterministic smoke backend.
