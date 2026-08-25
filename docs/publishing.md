# Build and publish

## Local checks

From the repository root:

```powershell
python -m pytest
python -m build --sdist
python -m twine check dist/*.tar.gz
```

Building a wheel locally compiles the bundled C++17 extension. On managed
Windows machines where compiler executables are blocked, use the wheel built by
GitHub Actions instead of attempting a local wheel build. A source distribution
can still be assembled with:

```powershell
python setup.py sdist
```

Release artifacts are a source distribution plus platform- and Python-specific
wheels, for example:

```text
taplite_calibration-0.3.1-cp39-cp39-win_amd64.whl
taplite_calibration-0.3.1-cp312-cp312-manylinux_2_17_x86_64.whl
taplite_calibration-0.3.1-cp312-cp312-macosx_11_0_arm64.whl
taplite_calibration-0.3.1.tar.gz
```

Inspect and install the wheel into a temporary environment before publication.
The package test suite covers the QA pass, repair, and unrecoverable-error
branches, integer path-screen multiplicity, a real tiny ODME problem, the smoke
workflow backend, and a real tiny native TAPLite auto calibration. The source
distribution ships the minimum kernel source and a tiny public test network,
but no prebuilt executable or proprietary network.

`examples/quickstart` is a complete synthetic dataset retained in the source
distribution so users can validate a clone without proprietary data. Its
`smoke` backend is for workflow verification only.

## Repository setup

The repository is `Mmdabb/taplite-calibration`. CI compiles and tests the native
module on Windows and Linux. The release workflow uses `cibuildwheel` to create
Windows, Linux, and macOS wheels and publishes them with the source distribution.

Protect the default branch, require the CI check, and create a `pypi`
environment. For PyPI Trusted Publishing, configure the PyPI project to trust
the repository's `.github/workflows/publish.yml` workflow and `pypi`
environment. The workflow publishes only when a GitHub release is published.

## First release

1. Confirm the distribution name is still available on PyPI.
2. Update `CHANGELOG.md` and the version in `pyproject.toml` and
   `src/taplite_calibration/__init__.py` together.
3. Run tests and build both artifacts.
4. Create and push tag `v0.3.1`.
5. Publish a GitHub release from that tag.
6. Confirm the Trusted Publishing workflow and then install from PyPI in a
   clean environment.

Do not commit NVTA inputs, outputs, fallback dictionaries, locally compiled
`.pyd`/`.so` files, credentials, or PyPI tokens. Native wheels belong in GitHub
release/PyPI artifacts, not source control. Trusted Publishing does not require
a stored API token.
