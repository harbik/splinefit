# dierckx-rs: Pure-Rust B-Spline Curve and Surface Fitting

A pure-Rust implementation of Paul Dierckx' B-spline fitting library.
No Fortran compiler, no C compiler, no build-time dependencies — just Rust.

The original Fortran library was written by Paul Dierckx in the mid-1980s and remains
one of the most complete and mathematically rigorous B-spline fitting libraries available.
It is described in detail in his book:
[Curve and Surface Fitting with Splines](https://www.google.com/books/edition/Curve_and_Surface_Fitting_with_Splines/-RIQ3SR0sZMC?hl=en),
Oxford University Press, 1993.

This crate is a line-by-line translation of the Fortran sources into safe (and where
necessary unsafe) Rust, preserving the original algorithms exactly.  Numerical output
is verified against [scipy](https://docs.scipy.org/doc/scipy/reference/interpolate.html),
which wraps the original Fortran library.

## Functions

| Function | Description |
|---|---|
| `curfit_` | Least-squares / smoothing spline for scattered 1-D data |
| `splev_` | Evaluate a B-spline (or its derivatives) |
| `spalde_` | Evaluate all derivatives of a B-spline at a point |
| `spalder_` | Evaluate derivative of a B-spline |
| `splint_` | Definite integral of a B-spline |
| `sproot_` | Real zeros of a cubic B-spline |
| `insert_` | Insert a knot into a B-spline |
| `curev_` | Evaluate a parametric spline curve |
| `cualde_` | All derivatives of a parametric spline curve |
| `clocur_` | Least-squares / smoothing closed parametric curve |
| `concur_` | Parametric spline curve with derivative endpoint constraints |

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
dierckx-rs = "0.1.0"
```

All functions use the same calling convention as the original Fortran routines
(Fortran-compatible C ABI, pointer arguments), so existing code targeting the
`dierckx-sys` crate can be ported by changing the crate name and import.

## Example

```rust
use std::f64::consts::PI;
use dierckx_rs::{curfit_, splev_};

let m = 20usize;
let k = 3i32;
let nest = m + k as usize + 1;
let lwrk = m * (k as usize + 1) + nest * (7 + 3 * k as usize);

let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
let w = vec![1.0f64; m];
let s = 0.05f64;

let (iopt, m_i, k_i, nest_i, lwrk_i) = (0i32, m as i32, k, nest as i32, lwrk as i32);
let (xb, xe) = (x[0], x[m - 1]);

let mut n = 0i32;
let mut t   = vec![0.0f64; nest];
let mut c   = vec![0.0f64; nest];
let mut fp  = 0.0f64;
let mut wrk = vec![0.0f64; lwrk];
let mut iwrk = vec![0i32; nest];
let mut ier  = 0i32;

unsafe {
    curfit_(
        &iopt, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
        &xb, &xe, &k_i, &s, &nest_i,
        &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
        wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier,
    );
}
assert!(ier <= 0, "curfit failed: ier={ier}");
```

## Relationship to `dierckx-sys`

The `dierckx-sys` crate wraps the original Fortran sources using `gfortran` and the
`cc` build crate, requiring a Fortran compiler at build time.  This crate (`dierckx-rs`)
replaces that with a pure-Rust translation, so it builds anywhere Rust builds.

## License

The original Dierckx Fortran algorithms (on which this translation is based) were
downloaded from [netlib.org](http://www.netlib.org/dierckx/) and carry no license
restrictions. Please acknowledge Paul Dierckx' work in your projects.

All Rust source in this repository is &copy; 2021 Harbers Bik LLC, dual-licensed under:

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
* MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
