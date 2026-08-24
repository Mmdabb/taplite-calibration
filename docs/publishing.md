# Build and publish

## Local checks

From the repository root:

```powershell
python -m pytest
python -m build
python -m twine check dist/*
```

If the `build` package is unavailable on a managed offline machine, the source
tree also includes a small compatibility `setup.py`:

```powershell
python setup.py sdist bdist_wheel
```

The expected release files are a pure-Python wheel and source distribution:

```text
dist/taplite_calibration-0.2.0-py3-none-any.whl
dist/taplite-calibration-0.2.0.tar.gz
```

Inspect and install the wheel into a temporary environment before publication.
The package test suite covers the QA pass, repair, and unrecoverable-error
branches, integer path-screen multiplicity, a real tiny ODME problem, and a
deterministic auto-calibration contract backend. It does not ship a network or
compiled TAPLite binary.

`examples/quickstart` is a complete synthetic dataset retained in the source
distribution so users can validate a clone without proprietary data. Its
`smoke` backend is for workflow verification only.

## Repository setup

Create `Mmdabb/taplite-calibration` on GitHub, then initialize this prepared
directory as that repository and push it. The project URLs in `pyproject.toml`
already point there.

Protect the default branch, require the CI check, and create a `pypi`
environment. For PyPI Trusted Publishing, configure the PyPI project to trust
the repository's `.github/workflows/publish.yml` workflow and `pypi`
environment. The workflow publishes only when a GitHub release is published.

## First release

1. Confirm the distribution name is still available on PyPI.
2. Update `CHANGELOG.md` and the version in `pyproject.toml` and
   `src/taplite_calibration/__init__.py` together.
3. Run tests and build both artifacts.
4. Create and push tag `v0.2.0`.
5. Publish a GitHub release from that tag.
6. Confirm the Trusted Publishing workflow and then install from PyPI in a
   clean environment.

Do not commit NVTA inputs, outputs, fallback dictionaries, native `.pyd` files,
credentials, or PyPI tokens. Trusted Publishing does not require a stored API
token.
