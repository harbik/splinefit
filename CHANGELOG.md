# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-04-22

### Fixed

- **`fpcons` (open parametric curve engine)** — the x- and y-coefficient blocks were
  written to the `c` array with stride `nest` (buffer capacity) but `curev_` and the
  high-level `ParameterSplineCurveFit` wrapper read them with stride `n` (actual knot
  count).  When `n < nest` — the common case — y-coordinates were read from uninitialised
  memory and came out as zero.  Fixed all six affected loops: QR update, back-substitution
  (knot phase), two fp-computation loops, and the regularisation Givens rotation and
  back-substitution.
- **`fpclos` (closed curve engine)** — two index expressions of the form
  `c[j1v-1+(j2-1)*nest]` double-counted the dimension offset that the `l0 += nest`
  accumulator already advances.  Reduced to `c[j1v-1]` at both sites.
- **`fpintb` (B-spline integral helper)** — `ia` and `ib` were `usize`, causing a
  wrapping-subtract panic in debug builds when the lower integration endpoint falls in
  the first knot span (`lk2 == 1`, so `lk2 - 2` underflows).  Changed both to `i32`.
- **`sproot_`** — boundary-validity loop used `t[jv - n_v]` which underflowed on the
  second iteration; rewritten as a straightforward `for i in 0..3` check on the first
  and last three knots.
- **`fpadpo`** — `fpinst` output coefficients were written into a temporary `Vec` that
  was immediately dropped; the slice now lives long enough and is copied back into `c`.
- **`cardinal_spline`** (`SplineCurveFit` and `ParameterSplineCurveFit`) — the `scan`
  accumulator started at `tb` instead of `tb + dt`, causing the left-boundary knot `tb`
  to appear `K+2` times instead of the required `K+1`.
- **`curfit`** — the `nest` value was recomputed from `m * k + 1` after a knot
  replacement instead of being derived from the actual buffer length, which could produce
  an inconsistent system size.
- **`concur` / `clocur`** — after a user-supplied knot vector was installed, the internal
  `t` buffer was not resized to `nest`; subsequent writes could exceed the buffer.
- **`ParameterSplineCurveFit::u()` / `xn()`** — returned slices sized by the knot count
  `n` instead of the data-point count `m`, producing truncated or out-of-bounds results
  before fitting.
- **`FitError::new`** — visibility tightened to `pub(crate)` and usage made consistent
  across all modules.

### Changed

- SPDX licence headers (`Apache-2.0 OR MIT`) added to every source file.
- Copyright year range updated to 2021–2026 in all source files and `README.md`.
- Unresolved doc links to type aliases (e.g. `[CubicSplineFit2D]`) fixed by qualifying
  them as `crate::TypeName` in submodule doc comments.
- Test bound in `curfit_sin_smoothing` tightened from `max_err < 0.2` to `< 0.1`
  (actual observed error ≈ 0.076).

### Tests

- `splev_constant_spline` — hand-built constant cubic B-spline; verifies `splev_` returns
  the constant value at every evaluation point.
- `splint_constant_and_sin` — verifies `splint_` on a constant spline (exact to 1 e-13)
  and on an interpolating sin spline (integral ≈ 2 on [0, π]).
- `sproot_sin_zeros` — fits an interpolating cubic spline to sin on [0, 2π] and verifies
  `sproot_` finds a zero near π.
- `insert_preserves_curve` — inserts a knot via `insert_` and confirms evaluation is
  unchanged to 1 e-13.
- `concur_unit_circle_interpolation` — fits a 2-D cubic parametric spline to a 3/4-circle
  arc with `concur_`, evaluates with `curev_`, and verifies the radius stays within 0.05
  of 1.0 at 50 dense parameter values.

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

[Unreleased]: https://github.com/harbik/splinefit/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/harbik/splinefit/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/harbik/splinefit/releases/tag/v0.1.0
