# taplite-calibration

`taplite-calibration` is a small, installable workflow package for regional
model input preparation and two calibration tasks:

1. package-specific target preparation and QA;
2. joint four-period screen-count OD matrix estimation (ODME); and
3. refined TAPLite QVDF auto calibration.

It can run either calibration task independently or `both`. Every execution
first applies the same prepared-input QA gate automatically. In `both` mode,
the order is prepare/QA, ODME, and auto calibration; ODME's adjusted NT/AM/MD/PM
scenarios become the demand input to auto calibration. The package contains the
minimum TAPLite assignment/refined-calibration C++ source needed to build its
native module. It does not contain CBI packages, NVTA networks, observations,
historical runs, prebuilt executables, or old TAPLite build trees.

## Install

Install the platform wheel that matches your Python and operating system:

```powershell
python -m pip install taplite_calibration-0.3.2-cp39-cp39-win_amd64.whl
```

The wheel already contains the TAPLite assignment kernel and refined
auto-calibration API. No `pytaplite` dependency, adjacent TAPLite4MPO checkout,
or TAPLite executable is required. Source installs compile the same bundled C++
code and therefore need a C++17 toolchain; managed Windows users should prefer
the prebuilt wheel. ODME uses the package's NumPy, pandas, and SciPy components.

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
│   ├── auto-calibration/refined_auto_calibration_settings.csv
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
    │   │   ├── od_factor_dictionary.npy
    │   │   └── factor_shards/*_od_factors.npz
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

`od_factor_dictionary.npy` is a NumPy dictionary indexing each period/mode's
compact per-OD factor shard; every shard includes origin, destination, factor,
original volume, and adjusted volume arrays. Auto calibration writes the
directly loadable `calibrated_qvdf_node_pair_dict.npy`, keyed by
`(from_node_id, to_node_id)` with all AM/MD/PM QVDF parameters.

## Run

For a complete public workflow example that needs no NVTA data or TAPLite executable,
see [`examples/quickstart`](examples/quickstart). After installing the package
from the clone, run:

```powershell
python examples/quickstart/run_example.py
```

The quickstart uses the deterministic smoke backend so it can focus on
prepare/QA, ODME, handoff, and artifact structure. To exercise a real native
TAPLite solve on the bundled four-node network, run:

```powershell
python examples/native-kernel/run_example.py --output native-example-output
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

Auto calibration calls the bundled
`taplite_calibration._native.run_auto_calibration_in_dir` entry point. Each
period runs in a short-lived Python worker for reliable cleanup of TAPLite's
regional route store; this is still an in-process `.pyd`/`.so` call and never
launches a TAPLite executable. Every inner assignment uses `column_output=2` in
memory, `column_file_output=0`, and no route/vehicle output. Intermediate outer
iterations are not serialized and reread. Input links must use `vdf_type=2` and
`qvdf_profile_mode=1`; the package fails before solving if either contract is
violated. Final AM/MD/PM modeled average speeds rebuild the shared QVDF speed
anchors, and the final parameter dictionary uses the existing node-pair schema.

The same wheel exposes ordinary TAPLite assignment independently of the
calibration workflow:

```python
from taplite_calibration import assign, native_status

print(native_status(16))
result = assign("path/to/scenario", settings_overrides={"number_of_processors": 16})
print(result.summary())
```

The `smoke` backend exists only for tests of orchestration and output contracts.
Never use it for a production calibration.

## Configuration and publication

See [configuration](docs/configuration.md), [architecture](docs/architecture.md),
and [publishing](docs/publishing.md). A data-free NVTA-style example is under
[`examples/nvta-project`](examples/nvta-project).

## License

MIT. See [LICENSE](LICENSE).
