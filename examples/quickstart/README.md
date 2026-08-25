# Working quickstart

This is a complete, synthetic three-zone example. It exercises the automatic
input QA gate, joint NT/AM/MD/PM ODME, adjusted-scenario handoff, and the
auto-calibration output contract. It contains no NVTA or proprietary data.

The auto-calibration backend is deliberately `smoke`: it verifies package
orchestration and artifacts without spending time on a native TAPLite solve.
Never use that backend for a production calibration; use the separate
`examples/native-kernel` network to test the real bundled C++ API.

From the repository root:

```powershell
python -m pip install -e .
python -m taplite_calibration validate --project examples/quickstart
python examples/quickstart/run_example.py
```

Validation should report `PASS`. The run script creates a unique directory
under `outputs/runs/`, verifies that ODME moved the daily screen volume closer
to the observation, verifies the calibrated QVDF dictionary, and prints the
key result paths. Re-running it creates another run and does not overwrite an
earlier result.

The example starts with 120 vehicles at its single daily screen and targets
132. Each period/mode keeps 60 total trips and the same productions and
attractions. One synthetic link has an episode target (`E`) and the other has a
no-episode target (`N`).

To recreate all CSV and NPZ input files deterministically:

```powershell
python examples/quickstart/generate_inputs.py
```
