# Contributing

Create a focused branch, keep datasets and generated outputs out of Git, add or
update tests, and run `python -m pytest` before opening a pull request. Public
interfaces and configuration keys require documentation and a changelog entry.

Use the `smoke` backend only for orchestration tests. Changes to the assignment
or calibration path must also pass `tests/test_native_kernel.py` against the
bundled extension. Do not commit local `.pyd`/`.so` files or native build trees;
release wheels are built from `native/` by GitHub Actions.
