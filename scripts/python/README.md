# Python examples

These scripts demonstrate the `splinefit` Python API.

## Setup

Create a virtual environment and install `splinefit` from PyPI:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install splinefit numpy
```

Alternatively, to build from source (requires the Rust toolchain):

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install maturin numpy
maturin develop --features python
```

## Running

With the virtual environment activated:

```sh
python3 scripts/python/smoothing_spline.py
python3 scripts/python/integration_and_roots.py
python3 scripts/python/cardinal_spline.py
```

## VSCode integration

If VSCode shows "splinefit could not be resolved", select the project's
virtual environment as the Python interpreter:

1. Open the Command Palette (**Cmd+Shift+P** on macOS, **Ctrl+Shift+P** on Linux/Windows)
2. Type **Python: Select Interpreter**
3. Choose the `.venv/bin/python` entry from this project

This tells both the editor (for autocompletion / type checking) and the
integrated terminal (for running scripts) to use the environment where
`splinefit` is installed.
