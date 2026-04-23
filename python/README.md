# splinefit

B-spline curve fitting for Python, powered by a native Rust engine. Fit smoothing, interpolating, or cardinal cubic splines to data and evaluate, integrate, or find roots -- with the same numerical accuracy as SciPy's FITPACK, at native speed.

Built from a pure-Rust translation of Paul Dierckx' classic [FITPACK](http://www.netlib.org/dierckx/) library (the same Fortran engine behind SciPy's `splrep`/`splev`).

## Installation

```sh
pip install splinefit
```

Requires Python >= 3.9 and NumPy.

## Quick start

```python
import numpy as np
from splinefit import CubicSpline

# Sample data: sin(x) on [0, 2pi]
x = np.linspace(0, 2 * np.pi, 50)
y = np.sin(x)

# Fit a smoothing spline (rms = 0.05)
spline = CubicSpline.smoothing(x, y, rms=0.05)
print(spline)  # CubicSpline(num_knots=11, domain=[0.000000, 6.283185])

# Evaluate at 200 points
x_new = np.linspace(0, 2 * np.pi, 200)
y_fit = spline.evaluate(x_new)
```

## API

### Constructors

All constructors take NumPy arrays (or array-like inputs) for `x` and `y`. `x` must be strictly increasing with at least 4 points.

#### `CubicSpline.smoothing(x, y, rms)`

Fit a smoothing spline. `rms` is the target root-mean-square residual in the same units as `y`. Smaller values produce more knots and a tighter fit.

```python
spline = CubicSpline.smoothing(x, y, rms=0.05)
```

#### `CubicSpline.interpolating(x, y)`

Fit an interpolating spline that passes exactly through every data point.

```python
spline = CubicSpline.interpolating(x, y)
```

#### `CubicSpline.cardinal(x, y, dt)`

Fit a spline on a fixed equidistant knot grid with spacing `dt`.

```python
spline = CubicSpline.cardinal(x, y, dt=0.5)
```

### Methods

#### `spline.evaluate(x) -> numpy.ndarray`

Evaluate the spline at each point in `x`.

```python
y_fit = spline.evaluate(np.array([0.5, 1.0, 1.5]))
```

#### `spline.integral(a, b) -> float`

Compute the definite integral of the spline over `[a, b]`.

```python
area = spline.integral(0, np.pi)  # integral of sin(x) from 0 to pi ~ 2.0
```

#### `spline.roots() -> numpy.ndarray`

Find all interior zeros of the spline, returned in ascending order. Zeros at the domain boundaries may not be included.

```python
zeros = spline.roots()  # e.g. array([3.14159...])
```

#### `spline.knots() -> numpy.ndarray`

Return the knot vector.

#### `spline.coefficients() -> numpy.ndarray`

Return the B-spline coefficients.

#### `spline.num_knots -> int`

Number of knots (property).

### String representation

```python
>>> spline
CubicSpline(num_knots=11, domain=[0.000000, 6.283185])
```

## Comparison with SciPy

`splinefit` uses the same FITPACK algorithms as SciPy but compiled as a native Rust extension rather than wrapping Fortran via `f2py`. The API is simpler:

| splinefit | SciPy equivalent |
|---|---|
| `CubicSpline.smoothing(x, y, rms)` | `splrep(x, y, s=m*rms**2)` |
| `CubicSpline.interpolating(x, y)` | `splrep(x, y, s=0)` |
| `spline.evaluate(x)` | `splev(x, tck)` |
| `spline.integral(a, b)` | `splint(a, b, tck)` |
| `spline.roots()` | `sproot(tck)` |

## How it works

This package is compiled from the Rust crate [`splinefit`](https://crates.io/crates/splinefit) using [PyO3](https://pyo3.rs) and [maturin](https://www.maturin.rs). It runs the same numerical algorithms as SciPy's FITPACK, translated line-by-line from Paul Dierckx' original Fortran into Rust.

## License

Apache-2.0 OR MIT
