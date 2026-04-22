# CLAUDE.md — splinefit

## Project identity

**Crate name:** `splinefit`
**Directory:** `/Users/gerardharbers/Projects/splinefit`
**Repository:** <https://github.com/harbik/splinefit>
**Owner:** Gerard Harbers / Harbers Bik LLC

This is a pure-Rust B-spline fitting library. It merges two previously separate crates:

- `dierckx-rs` — a line-by-line Rust translation of Paul Dierckx' Fortran FITPACK library
- `splinify` — a high-level Rust wrapper that previously used `dierckx-sys` (a Fortran FFI crate) as its backend

The merger means there are now no build-time dependencies on a Fortran or C compiler.

## Architecture

```
src/
├── lib.rs          — public API surface: re-exports all high-level types and type aliases,
│                     defines FitError, Result. The dierckx module is private.
├── dierckx.rs      — private backend: ~3 200-line translation of the Fortran library.
│                     All internal helper fns are private; the 11 C-ABI entry points are
│                     pub unsafe extern "C" (visible within the crate only).
├── curfit.rs       — SplineCurveFit<K>: 1-D smoothing / interpolating spline fits
├── concur.rs       — ParameterSplineCurveFit<K,N>: parametric curve fits (open)
├── clocur.rs       — ClosedParameterSplineCurveFit<K,N>: closed/periodic curve fits
├── evaluate.rs  — evaluate(): calls splev_ (N=1) or curev_ (N>1) to eval a SplineCurve
└── util.rs         — CSV helpers (read_csv_xy, read_csv_uxy, write_csv_xy), SplineCurveData
```

## Key design decisions

**The `dierckx` module is intentionally private.** The public API is the high-level
`SplineCurveFit` / `ParameterSplineCurveFit` / `ClosedParameterSplineCurveFit` types.
Do not re-export the low-level `curfit_`, `splev_`, etc. from `lib.rs`.

**`SplineCurve<K, N>` comes from the `spliny` crate** (also owned by Gerard Harbers).
It is the output type of all fit methods. Evaluation is done via `evaluate::evaluate()`,
not via `spliny`'s own `eval()` method (whose return type is not a plain `f64`).

**The C-ABI functions use Fortran-style pointer arguments** (`*const i32`, `*mut f64`, etc.).
When calling them from the high-level wrappers, use the raw pointer helpers on `Vec`:

- `vec.as_ptr()` for read-only inputs
- `vec.as_mut_ptr()` for in/out parameters
- Note: `clocur_` takes `u: *mut f64` (it writes computed parameter values back).

**Workspace sizing formulas** (must match the Fortran source exactly):

| Function | `lwrk` formula |
|---|---|
| `curfit_` | `m*(k+1) + nest*(7+3*k)` |
| `clocur_` | `m*(k+1) + nest*(7+idim+5*k)` |
| `concur_` | `m*(k+1) + nest*(6+idim+3*k)` |

`nest` is the maximum number of knots allocated. For `curfit_`: `nest = m + k + 1` (s=0)
or larger. For `clocur_`: `nest = m + 2*k`.

## The `dierckx.rs` backend

The file is ~3 200 lines and is structured in seven sections:

| Section | Content |
|---|---|
| 1 | Core primitives: `fpbspl`, `fpgivs`, `fprota`, `fpback`, `fpbacp` |
| 2 | Knot utilities: `fpdisc`, `fpknot`, `fprati`, `fpchec`, `fpchep`, `fpched` |
| 3 | Basis function utilities: `fpader`, `fpintb`, `fpcuro`, `fpinst` |
| 4 | Fourier utilities: `fppara`, `fpcsin`, `fpbfou` (dead_code — kept for completeness) |
| 5 | Core fitting algorithms: `fpcurf` (main 1-D engine), `fpclos`, `fpcons` |
| 6 | Endpoint-constrained / closed-curve helpers: `fppocu`, `fpadpo` |
| 7 | Public C-ABI wrappers: `curfit_`, `splev_`, `curev_`, `spalde_`, `cualde_`, `spalder_`, `clocur_`, `concur_`, `splint_`, `sproot_`, `insert_` |

The whole file carries `#![allow(clippy::all, non_snake_case, unused_variables, unused_assignments, unused_mut)]`
to suppress noise from the mechanical Fortran translation. Do not remove this.

Functions `fppara`, `fpcsin`, and `fpbfou` (Section 4) are marked `#[allow(dead_code)]`
because they implement the Fourier-spline routines that have not yet been wired up to
a high-level wrapper. Keep them; they will be needed when surface/Fourier fitting is added.

## Translation conventions (Fortran → Rust)

- Fortran `real` compiled with `-fdefault-real-8` → `f64`
- Fortran `integer` → `i32` (or `usize` for indexing)
- Fortran 1-based arrays → 0-based in Rust (`t[i-1]` in Rust = `t(i)` in Fortran)
- Fortran 2D column-major `a(nest, k)` → `a[(i-1) + (j-1)*nest]`
- Numerical correctness is verified against SciPy (which wraps the original Fortran)

## Tests

**Internal unit tests** (in `src/dierckx.rs`, `#[cfg(test)] mod tests`):

- `fpcurf_sin_smoothing_direct` — calls internal `fpcurf` directly
- `curfit_sin_interpolation` — scipy reference, ier=-1, n=13
- `curfit_sin_smoothing` — scipy reference, s=0.05, n=11
- `curfit_quadratic_exact_fit` — degree-2 fits x² exactly
- `curfit_large_s_gives_polynomial` — large s → ier=-2, n=2*(k+1)
- `clocur_unit_circle_interpolation` — closed curve on unit circle
- `curfit_iopt1_restart_consistent` — warm-restart call

**External integration tests** (`tests/integration.rs`):

- `cubic_smoothing_spline_sin` — exercises the high-level `CubicSplineFit` API end-to-end

Run all tests: `cargo test`

## Adding a new high-level fit type

1. Identify which low-level function in `dierckx.rs` implements it (e.g. `percur_` for
   periodic 1-D curves).
2. If the Fortran routine has not been translated yet, add it to `dierckx.rs` Section 7
   following the existing `pub unsafe extern "C"` pattern.
3. Create a new `src/<name>.rs` file modelled on `curfit.rs` or `clocur.rs`.
   Import the low-level function as `use crate::dierckx::<fn_name>;`.
4. Declare the module and re-export it in `lib.rs`.
5. Add type aliases in `lib.rs` if appropriate (`pub type CubicFoo = Foo<3>`).
6. Write an internal `#[cfg(test)]` test against a scipy reference value.
7. Write an external `tests/integration.rs` test against the high-level API.

## Dependencies

| Crate | Purpose |
|---|---|
| `spliny` | `SplineCurve<K,N>` output type (also by Gerard Harbers) |
| `csv` | CSV helpers in `util.rs` |
| `serde` / `serde_json` | `SplineCurveData` serialisation |
| `bitflags` | (available for future flag-style options) |
| `plotters` | optional (`plot` feature) for example visualisations |

## What NOT to do

- Do not add `dierckx-sys` as a dependency — the whole point of this crate is to avoid it.
- Do not expose the `dierckx` module or its functions publicly.
- Do not change the numeric algorithms in `dierckx.rs` without a scipy reference test
  confirming the output is still correct to ≥ 1e-10.
- Do not remove `fppara`, `fpcsin`, or `fpbfou` — they are intentionally kept.
