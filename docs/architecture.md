# Architecture

The distribution uses a `src/` layout and keeps algorithms separate from the
CLI. `taplite_calibration.pipeline` is the public orchestrator. It loads and
validates relative paths, creates an immutable run directory, dispatches
stages, and records status and portable artifacts.

```text
calibration.toml
       │
       ▼
 prepared-input QA gate
       │
       ├── PASS ───────────────► reuse prepared inputs
       ├── FAIL + raw sources ─► preprocess ─► re-QA
       └── FAIL + no sources ──► record error and exit
                                      │
                                      ▼
       mode=odme ──────────────────► ODME
       mode=auto-calibration ──────► auto calibration
       mode=both ──────────────────► ODME ─► adjusted scenarios ─► auto calibration
                                      │
                                      ▼
                              run manifest + artifacts
```

## ODME stage

The ODME module loads sparse OD-to-screen policy matrices and applies bounded
projected mirror descent. Alternating bounded multiplicative projections return
each candidate to the original productions and attractions. An optional convex
retraction to the original OD enforces a per-period total-variation budget.

The stage writes adjusted demand scenarios, per-mode factor shards, a portable
factor index dictionary, conservation audits, optimization history, and daily
screen comparisons. Network files are hard-linked when the platform permits;
mutable demand files are copied before modification.

## Prepare stage

The prepare stage is a run-local gate with structured before/after QA. Daily
screen data can be normalized from arbitrary source headings. Sparse policies
can be staged from a validated source or built by streaming DTAC-v2. The policy
builder counts every screen-member link on a path, preserving integer
path-screen incidence before path shares are aggregated to OD coefficients.

The auto-target builder consumes CBI producer outputs without importing or
bundling CBI. It constructs full-period mean speed and S3 targets for canonical
GP links, adds candidate-episode P/vt2 when available, slices long episodes at
period boundaries, records mapping choices, marks managed/unmapped links U, and
sets the TAPLite QVDF/profile execution contract.

## Auto-calibration stage

The distribution compiles the minimum native source set under `native/` into
`taplite_calibration._native`. Its binding exposes ordinary `AssignmentAPI`
and the separate `AutoCalibrationAPI`; no external `pytaplite` package or
TAPLite executable participates in a packaged run.

The auto-calibration module stages one period at a time and invokes the
dedicated API in an isolated Python worker. The worker reads network and demand
once, retains every inner equilibrium's DTAC-v2 route/arrival state in memory,
writes only accepted final outputs, and exits to release process-scoped regional
memory. Large immutable inputs are hard-linked where possible. Each solve gets
the selected processor budget, while AM/MD/PM remain sequential.

After all three periods finish, a one-time finalizer merges accepted native
audit parameters into calibrated link files. It then rebuilds consistent speed
anchors and writes a portable node-pair QVDF dictionary. No CBI, NVTA network,
old build tree, or calibration run is package data.

## Failure behavior

Run directories are never overwritten. A failure leaves `run_manifest.json`
with `status=failed` and a traceback file for debugging. Completed runs include
the configuration SHA-256, package/runtime versions, stage results, `run.log`,
and an artifact index.
