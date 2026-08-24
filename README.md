# taplite-calibration

`taplite-calibration` is a small, installable workflow package for regional
model input preparation and two calibration tasks:

1. package-specific target preparation and QA;
2. joint four-period screen-count OD matrix estimation (ODME); and
3. refined TAPLite QVDF auto calibration.

It can run either calibration task independently or `both`. Every execution
first applies the same prepared-input QA gate automatically. In `both` mode,
the order is prepare/QA, ODME, and auto calibration; ODME's adjusted NT/AM/MD/PM
scenarios become the demand input to auto calibration. The package does not
contain TAPLite4MPO source, CBI packages, NVTA networks, observations, compiled
binaries, or historical runs.

## Install

Install the pure-Python wheel:

```powershell
python -m pip install taplite_calibration-0.2.0-py3-none-any.whl
```

For production auto calibration, also install a compatible TAPLite4MPO wheel
whose `pytaplite` module exposes `auto_calibrate`. A local development checkout
can instead be referenced with the relative `pytaplite_path` configuration.
ODME alone needs only the package's normal NumPy, pandas, and SciPy dependencies.

On managed Windows systems, use the module form below. It calls Python and the
native `.pyd` API without depending on a separate TAPLite executable or the
console-script `.exe` launcher.

```powershell
python -m taplite_calibration --help
```

## Start a project

```powershell
python -m taplite_calibration init my-calibration-project
```

The command creates configuration and empty input/output directories only. It
does not copy networks or other large datasets.

```text
my-calibration-project/
├── calibration.toml
├── inputs/
│   ├── scenarios/{nt,am,md,pm}/
│   ├── odme/policies/{nt,am,md,pm}/
│   ├── observations/
│   ├── auto-calibration/
│   └── raw/                         # optional producer/source inputs
└── outputs/runs/<run-id>/
    ├── run_manifest.json
    ├── logs/run.log
    ├── 00-prepare/
    │   ├── manifest.json
    │   ├── odme/                     # normalized screens/policies when needed
    │   ├── auto-target-audits/       # TMC/link selection audit when needed
    │   └── prepared-scenarios/       # calibration-ready scenarios when needed
    ├── 01-odme/                       # when selected
    │   ├── results/
    │   └── adjusted-scenarios/
    ├── 01- or 02-auto-calibration/    # number reflects stage order
    │   ├── period-runs/{am,md,pm}/
    │   └── finalized/assignment/{am,md,pm}/
    └── artifacts/
        ├── artifact_index.json
        └── calibrated_qvdf_node_pair_dict.npy
```

All paths stored in `calibration.toml`, manifests, and dictionary indexes are
relative. A configuration may use `..` to reference a separately managed data
directory without copying it into this repository.

## Run

For a complete public example that needs no NVTA data or TAPLite executable,
see [`examples/quickstart`](examples/quickstart). After installing the package
from the clone, run:

```powershell
python examples/quickstart/run_example.py
```

Validate before starting a regional solve:

```powershell
python -m taplite_calibration validate --project my-calibration-project
```

Run the configured mode:

```powershell
python -m taplite_calibration run `
  --project my-calibration-project `
  --run-id nvta-production-01
```

Override the mode or processor count without editing TOML:

```powershell
python -m taplite_calibration run `
  --project my-calibration-project `
  --mode odme `
  --processors 16
```

Valid modes are `odme`, `auto-calibration`, and `both`. Preparation is not a
mode or an option: the package always checks and, when possible, prepares the
inputs required by the selected mode. Processor values are restricted to
1--20. AM, MD, and PM auto-calibration solves execute sequentially to avoid
regional memory contention; each native solve uses the full selected processor
budget internally.

The equivalent Python API is:

```python
from pathlib import Path
from taplite_calibration import run_project

result = run_project(
    Path("my-calibration-project"),
    mode="both",
    processors=16,
    run_id="nvta-production-01",
)
print(result.manifest_path)
```

## Method contracts

### Preparation and QA

Before a downstream stage starts, the package performs schema, coverage,
positivity, consistency, and QVDF-contract checks on its prepared inputs:

```text
prepared-input QA
├── PASS ───────────────► reuse inputs and continue
└── FAIL
    ├── raw sources sufficient ─► warn, preprocess, re-QA, continue
    └── raw sources insufficient ► error, record failure, exit
```

ODME preparation normalizes the observed daily screen table and either stages
existing sparse policies or streams DTAC-v2 route pools to build them. The
streaming builder uses binary link membership and integer path-screen
multiplicity, so a route contributes once for every screen-member link it
contains.

Auto-calibration preparation reads average-weekday TMC profiles and candidate
episodes, optionally slices congestion episodes at AM/MD/PM boundaries,
computes period-average speed and S3, maps the frozen canonical TMC selection to
links, excludes managed facilities as class U, and enriches clean assignment
scenarios. Actual and virtual CBI roots are supported, including the expanded
coverage snapshot layout used by the latest NVTA run.

Preparation is run-local and never overwrites raw or previously prepared
inputs. Both the failed initial QA and passing post-preparation QA are retained
in `00-prepare/manifest.json`.

### ODME

ODME compares the sum of NT, AM, MD, and PM fixed-route screen contributions
directly with observed daily screen totals. It does not create period screen
targets and does not use CUBE or Vs3 as targets or priors. Each period/mode
keeps its original positive support, productions, attractions, and total
demand; no demand moves between periods. Cell factor bounds, KL regularization,
and an optional common per-period total-variation budget are configurable.

For the screen operator, `M[s,a]` is binary for counted link `a`, while path
incidence `A[s,r] = sum(a in r) M[s,a]` is integer. A vehicle contributes once
for each counted member link it traverses.

### Auto calibration

Auto calibration calls `pytaplite.auto_calibrate` directly. Every inner
assignment uses `column_output=2` in memory, `column_file_output=0`, and no
route/vehicle output. Intermediate outer iterations are not serialized and
reread. Input links must use `vdf_type=2` and `qvdf_profile_mode=1`; the package
fails before solving if either contract is violated. Final AM/MD/PM modeled
average speeds are used to rebuild the shared QVDF speed anchors, and the final
parameter dictionary uses the existing node-pair schema.

The `smoke` backend exists only for tests of orchestration and output contracts.
Never use it for a production calibration.

## Configuration and publication

See [configuration](docs/configuration.md), [architecture](docs/architecture.md),
and [publishing](docs/publishing.md). A data-free NVTA-style example is under
[`examples/nvta-project`](examples/nvta-project).

## License

MIT. See [LICENSE](LICENSE).
