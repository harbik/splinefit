<!-- cargo-rdme start -->

# splinefit: B-Spline Curve and Surface Fitting

A pure-Rust B-spline fitting library with an ergonomic high-level API, backed by a
line-by-line Rust translation of Paul Dierckx' classic Fortran fitting library.

No Fortran compiler, no C compiler, no build-time dependencies — just Rust.

The original Fortran library was written by Paul Dierckx in the mid-1980s and remains
one of the most complete and mathematically rigorous B-spline fitting libraries available.
It is described in detail in his book:
[Curve and Surface Fitting with Splines](https://www.google.com/books/edition/Curve_and_Surface_Fitting_with_Splines/-RIQ3SR0sZMC?hl=en),
Oxford University Press, 1993.

Numerical output is verified against [scipy](https://docs.scipy.org/doc/scipy/reference/interpolate.html),
which wraps the original Fortran library.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
splinefit = "0.3.0"
```

## Example

Fit a smoothing cubic spline to noisy data, then evaluate it:

```rust
use splinefit::{CubicSplineFit, evaluate};
use std::f64::consts::PI;

let m = 20usize;
let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();

let spline = CubicSplineFit::new(x.clone(), y)
    .smoothing_spline(0.05)
    .expect("fit failed");

let y_fit = evaluate::evaluate(&spline, &x).expect("evaluate failed");
assert!((y_fit[0] - x[0].sin()).abs() < 0.1);
```

## Fit types

| Type alias | Description |
|---|---|
| `CubicSplineFit` | Smoothing / interpolating 1-D cubic spline |
| `LinearSplineFit` / `QuinticSplineFit` | Same, degree 1 or 5 |
| `CubicSplineFit2D` | Parametric cubic curve in 2-D |
| `CubicSplineFit3D` | Parametric cubic curve in 3-D |
| `ClosedCubicSplineFit2D` | Closed (periodic) parametric curve in 2-D |
| `ClosedCubicSplineFit3D` | Closed (periodic) parametric curve in 3-D |

Each fit type exposes three methods:

| Method | Description |
|---|---|
| `.interpolating_spline()` | Exact interpolation through all data points |
| `.smoothing_spline(rms)` | Fewest knots with RMS error ≤ `rms` |
| `.cardinal_spline(dt)` | Weighted least-squares on a fixed equidistant knot grid |

## Relationship to `dierckx-sys` and `splinify`

`dierckx-sys` wrapped the original Fortran sources via `gfortran`, requiring a Fortran
compiler at build time. `splinify` used `dierckx-sys` as its backend.

This crate (`splinefit`) combines the high-level API from `splinify` with a pure-Rust
translation of the Fortran engine, so it builds anywhere Rust builds with no external
toolchain dependencies.

## Performance

The pure-Rust engine matches or outperforms SciPy's Fortran-compiled FITPACK.
Benchmarks on Apple M1 (macOS), fitting and evaluating `sin(x)` on `[0, 2π]`:

**1-D smoothing fit**

| Data points | Rust | SciPy/Fortran |
|---:|---:|---:|
| 50 | 17 µs | 21 µs |
| 200 | 57 µs | 62 µs |
| 1 000 | 279 µs | 282 µs |
| 5 000 | 1.35 ms | 1.4 ms |

**1-D evaluation (10× data points)**

| Data points | Rust | SciPy/Fortran |
|---:|---:|---:|
| 50 | 12 µs | 18 µs |
| 200 | 50 µs | 61 µs |
| 1 000 | 254 µs | 307 µs |
| 5 000 | 1.27 ms | 1.5 ms |

**2-D parametric evaluation (10× data points)**

| Data points | Rust | SciPy/Fortran |
|---:|---:|---:|
| 50 | 12 µs | 34 µs |
| 200 | 50 µs | 121 µs |
| 1 000 | 248 µs | 611 µs |

Reproduce with `cargo bench` (Rust) and `python3 benches/scipy_bench.py` (SciPy).

## WebAssembly / JavaScript

This crate compiles to WebAssembly via [`wasm-pack`](https://rustwasm.github.io/wasm-pack/).
Enable the `wasm` feature to expose a `CubicSpline` class to JavaScript.

```html
<script type="module">
  import init, { CubicSpline } from "https://esm.sh/splinefit";
  await init();

  // Fit a smoothing spline to sin(x)
  const x = Float64Array.from({ length: 50 }, (_, i) => i * 2 * Math.PI / 49);
  const y = x.map(Math.sin);
  const spline = CubicSpline.smoothing(x, y, 0.05);

  // Evaluate at 200 points
  const xNew = Float64Array.from({ length: 200 }, (_, i) => i * 2 * Math.PI / 199);
  const yFit = spline.evaluate(xNew);

  // Integration and root-finding
  const area = spline.integral(0, Math.PI);  // ∫₀^π sin(x) dx ≈ 2
  const zeros = spline.roots();               // interior zeros
</script>
```

## License

The original Dierckx Fortran algorithms (on which this translation is based) were
downloaded from [netlib.org](http://www.netlib.org/dierckx/) and carry no license
restrictions. Please acknowledge Paul Dierckx' work in your projects.

All Rust source in this repository is &copy; 2021–2026 Harbers Bik LLC, dual-licensed under:

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
* MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

<!-- cargo-rdme end -->
