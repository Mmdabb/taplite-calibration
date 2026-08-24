# Contributing

Create a focused branch, keep datasets and generated outputs out of Git, add or
update tests, and run `python -m pytest` before opening a pull request. Public
interfaces and configuration keys require documentation and a changelog entry.

Use the `smoke` backend only in tests. Production behavior must be verified
against a compatible `pytaplite.auto_calibrate` native wheel.
