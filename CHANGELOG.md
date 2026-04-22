# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-04-22

### Added

- Pure-Rust translation of Paul Dierckx' Fortran FITPACK library (`dierckx.rs`),
  eliminating all build-time dependencies on a Fortran or C compiler.
- `SplineCurveFit<K>` — 1-D smoothing / interpolating spline fits via `curfit_`.
- `ParameterSplineCurveFit<K, N>` — open parametric curve fits in N dimensions via `concur_`.
- `ClosedParameterSplineCurveFit<K, N>` — closed/periodic parametric curve fits via `clocur_`.
- Type aliases for common degrees and dimensions: `CubicSplineFit`, `CubicSplineFit2D`,
  `ClosedCubicSplineFit2D`, and their Linear / Quintic variants.
- Three fit methods on every builder type: `interpolating_spline`, `smoothing_spline`,
  and `cardinal_spline`.
- `smoothing_spline_optimize` on parametric builders for iterative RMS-guided refinement.
- `evaluate::evaluate` — evaluates a `SplineCurve` via `splev_` (1-D) or `curev_` (N-D).
- CSV helpers: `read_csv_xy`, `read_csv_uxy`, `write_csv_xy`.
- `SplineCurveData` — serialisable snapshot of knot vector and B-spline coefficients.
- `xtask` with `cargo xtask test` and `cargo xtask doc` commands.
- Numerical correctness verified against SciPy reference values.

[Unreleased]: https://github.com/harbik/splinefit/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/harbik/splinefit/releases/tag/v0.1.0
