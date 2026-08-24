# Configuration reference

Configuration uses TOML schema version 1. Every path value must be relative to
the project directory supplied to `--project`.

## Project and pipeline

```toml
schema_version = 1

[project]
input_dir = "inputs"
output_dir = "outputs"

[pipeline]
mode = "both"
processors = 20
```

`mode` is `odme`, `auto-calibration`, or `both`. The command-line `--mode` and
`--processors` options override these values for one run. Preparation is an
always-active input gate, not a mode. The package enforces a maximum of 20
processors.

## Prepare and QA

The QA gate runs automatically before every selected calibration. The prepare
subtables provide raw-source wiring used only if prepared QA fails.

```toml
[prepare.odme]
source_scenario_root = "inputs/raw/scenarios"
daily_screen_source_csv = "inputs/raw/screen-counts.csv"
screen_id_column = "Screen Code"
observed_column = "Obs.2025 (Apply Growth to 2018 COG Cnts)"

# Choose one policy preparation source:
policy_source_root = "inputs/raw/prebuilt-od-screen-policies"
# route_run_root = "inputs/raw/dtac-route-runs"

[prepare.auto_calibration]
source_scenario_root = "inputs/raw/scenarios"
coverage_root = "inputs/raw/tmc-observation-coverage"
departure_profile_csv = "inputs/raw/departure_profiles.csv"
calibration_settings_csv = "inputs/raw/refined_volume_floor_settings.csv"
episode_period_policy = "split_intersection"
```

`daily_screen_source_csv` is normalized to the two-column ODME target schema.
Whitespace and line breaks in a configured source-column name are normalized,
which supports the multiline VDOT/COG observed-count heading.

For OD policies, `policy_source_root` stages an existing validated four-period
policy set. Alternatively, `route_run_root/{nt,am,md,pm}` can contain
`route_columns.bin`, `node.csv`, `link.csv`, `mode_type.csv`, and demand inputs.
The package then streams DTAC-v2 without loading the regional route pool into
memory.

For auto-calibration targets, `coverage_root` supports the expanded NVTA
snapshot convention (`cbi/actual`, optional `cbi/virtual`, frozen canonical
mapping, and GP/managed decisions). A generic layout can instead set
`cbi_actual_root`, optional `cbi_virtual_root`, `canonical_mapping_csv`, and
optional `facility_mapping_csv` explicitly. Raw data is read but never modified.

If prepared QA passes, all raw-source settings are ignored. If it fails and the
required source wiring is incomplete, the run records an error and exits.

## ODME

```toml
[odme]
policy_root = "inputs/odme/policies"
scenario_root = "inputs/scenarios"
daily_target_csv = "inputs/observations/daily_screens.csv"
daily_target_column = "observed_daily_vehicle_volume"
factor_cap = 0.20
kl_weight = 0.02
# tv_budget = 0.03
max_iterations = 12
initial_step = 0.15
maximum_step = 0.35
minimum_step = 0.0001
minimum_improvement = 1e-11
max_line_search = 8
projection_iterations = 10000
projection_tolerance = 0.00002
```

`factor_cap=0.20` gives cell factors `[0.80, 1.20]`. `tv_budget`, when present,
is a common normalized redistribution budget applied separately to every
period. For example, `0.03` permits a 3% period-level total-variation share.

Each period scenario contains `mode_type.csv`. Every listed `demand_file` must
be a CSV or DTAB-v1 binary. For each mode, the corresponding policy is:

```text
<policy_root>/<period>/od_screen_policy_<mode_type>.npz
```

Policy stores contain `origin`, `destination`, `q0`, `screen_ids`,
`zone_external`, and CSR arrays `data`, `indices`, and `indptr`. All periods and
modes must share their screen and zone axes.

The daily-target CSV contains a unique `screen_id` and the configured positive
target column. CUBE and Vs3 proportions are intentionally absent from the
objective.

## Auto calibration

```toml
[auto_calibration]
scenario_root = "inputs/scenarios"
backend = "native"
periods = ["am", "md", "pm"]
timeout_seconds = 86400
# calibration_settings_csv = "inputs/auto-calibration/refined_settings.csv"
# fallback_qvdf_dictionary = "inputs/auto-calibration/base_qvdf.npy"
# pytaplite_path = "../TAPLite4MPO"
```

AM, MD, and PM are calibrated as one consistent set. If
`calibration_settings_csv` is omitted, each scenario must contain its own
`auto_calibration_settings.csv`. `fallback_qvdf_dictionary` is needed only when
some node pairs do not exist in all three period networks or when additional
base-dictionary coverage must be preserved.

`pytaplite_path` supports an adjacent development checkout. Installed wheels do
not need it. The path must point to the directory from which `import pytaplite`
works, not to an executable.

The native preflight requires every link row to have `vdf_type=2` and
`qvdf_profile_mode=1`. The speed-anchor dictionary and network conversion are
upstream responsibilities; this package consumes the prepared, consistent
scenario and then rebuilds final AM/MD/PM boundary speeds after calibration.
