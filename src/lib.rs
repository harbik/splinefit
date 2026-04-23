// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

//! # splinefit: B-Spline Curve and Surface Fitting
//!
//! A pure-Rust B-spline fitting library with an ergonomic high-level API, backed by a
//! line-by-line Rust translation of Paul Dierckx' classic Fortran fitting library.
//!
//! No Fortran compiler, no C compiler, no build-time dependencies — just Rust.
//!
//! The original Fortran library was written by Paul Dierckx in the mid-1980s and remains
//! one of the most complete and mathematically rigorous B-spline fitting libraries available.
//! It is described in detail in his book:
//! [Curve and Surface Fitting with Splines](https://www.google.com/books/edition/Curve_and_Surface_Fitting_with_Splines/-RIQ3SR0sZMC?hl=en),
//! Oxford University Press, 1993.
//!
//! Numerical output is verified against [scipy](https://docs.scipy.org/doc/scipy/reference/interpolate.html),
//! which wraps the original Fortran library.
//!
//! ## Usage
//!
//! Add this to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! splinefit = "0.2.0"
//! ```
//!
//! ## Example
//!
//! Fit a smoothing cubic spline to noisy data, then evaluate it:
//!
//! ```rust
//! use splinefit::{CubicSplineFit, evaluate};
//! use std::f64::consts::PI;
//!
//! let m = 20usize;
//! let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
//! let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
//!
//! let spline = CubicSplineFit::new(x.clone(), y)
//!     .smoothing_spline(0.05)
//!     .expect("fit failed");
//!
//! let y_fit = evaluate::evaluate(&spline, &x).expect("evaluate failed");
//! assert!((y_fit[0] - x[0].sin()).abs() < 0.1);
//! ```
//!
//! ## Fit types
//!
//! | Type alias | Description |
//! |---|---|
//! | `CubicSplineFit` | Smoothing / interpolating 1-D cubic spline |
//! | `LinearSplineFit` / `QuinticSplineFit` | Same, degree 1 or 5 |
//! | `CubicSplineFit2D` | Parametric cubic curve in 2-D |
//! | `CubicSplineFit3D` | Parametric cubic curve in 3-D |
//! | `ClosedCubicSplineFit2D` | Closed (periodic) parametric curve in 2-D |
//! | `ClosedCubicSplineFit3D` | Closed (periodic) parametric curve in 3-D |
//!
//! Each fit type exposes three methods:
//!
//! | Method | Description |
//! |---|---|
//! | `.interpolating_spline()` | Exact interpolation through all data points |
//! | `.smoothing_spline(rms)` | Fewest knots with RMS error ≤ `rms` |
//! | `.cardinal_spline(dt)` | Weighted least-squares on a fixed equidistant knot grid |
//!
//! ## Relationship to `dierckx-sys` and `splinify`
//!
//! `dierckx-sys` wrapped the original Fortran sources via `gfortran`, requiring a Fortran
//! compiler at build time. `splinify` used `dierckx-sys` as its backend.
//!
//! This crate (`splinefit`) combines the high-level API from `splinify` with a pure-Rust
//! translation of the Fortran engine, so it builds anywhere Rust builds with no external
//! toolchain dependencies.
//!
//! ## Performance
//!
//! The pure-Rust engine matches or outperforms SciPy's Fortran-compiled FITPACK.
//! Benchmarks on Apple M1 (macOS), fitting and evaluating `sin(x)` on `[0, 2π]`:
//!
//! **1-D smoothing fit**
//!
//! | Data points | Rust | SciPy/Fortran |
//! |---:|---:|---:|
//! | 50 | 17 µs | 21 µs |
//! | 200 | 57 µs | 62 µs |
//! | 1 000 | 279 µs | 282 µs |
//! | 5 000 | 1.35 ms | 1.4 ms |
//!
//! **1-D evaluation (10× data points)**
//!
//! | Data points | Rust | SciPy/Fortran |
//! |---:|---:|---:|
//! | 50 | 12 µs | 18 µs |
//! | 200 | 50 µs | 61 µs |
//! | 1 000 | 254 µs | 307 µs |
//! | 5 000 | 1.27 ms | 1.5 ms |
//!
//! **2-D parametric evaluation (10× data points)**
//!
//! | Data points | Rust | SciPy/Fortran |
//! |---:|---:|---:|
//! | 50 | 12 µs | 34 µs |
//! | 200 | 50 µs | 121 µs |
//! | 1 000 | 248 µs | 611 µs |
//!
//! Reproduce with `cargo bench` (Rust) and `python3 benches/scipy_bench.py` (SciPy).
//!
//! ## License
//!
//! The original Dierckx Fortran algorithms (on which this translation is based) were
//! downloaded from [netlib.org](http://www.netlib.org/dierckx/) and carry no license
//! restrictions. Please acknowledge Paul Dierckx' work in your projects.
//!
//! All Rust source in this repository is &copy; 2021–2026 Harbers Bik LLC, dual-licensed under:
//!
//! * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
//! * MIT license ([LICENSE-MIT](LICENSE-MIT))
//!
//! at your option.
//!
//! ## Contribution
//!
//! Unless you explicitly state otherwise, any contribution intentionally submitted
//! for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
//! dual licensed as above, without any additional terms or conditions.

// Private backend: pure-Rust translation of the Dierckx Fortran library.
mod dierckx;

// Public high-level API.
pub mod curfit;
pub use curfit::*;

pub mod concur;
pub use concur::*;

pub mod clocur;
pub use clocur::*;

pub mod evaluate;

pub mod ops;

pub mod util;
pub use util::*;

use std::error;
use std::fmt;

pub type Result<T> = std::result::Result<T, Box<dyn error::Error>>;

// Single Output Spline Fit
pub type LinearSplineFit = SplineCurveFit<1>;
pub type CubicSplineFit = SplineCurveFit<3>;
pub type QuinticSplineFit = SplineCurveFit<5>;

// Multi-Output Parametrized Curve Fits
pub type LinearSplineFit1D = ParameterSplineCurveFit<1,1>;
pub type CubicSplineFit1D = ParameterSplineCurveFit<3,1>;
pub type QuinticSplineFit1D = ParameterSplineCurveFit<5,1>;
pub type LinearSplineFit2D = ParameterSplineCurveFit<1,2>;
pub type CubicSplineFit2D = ParameterSplineCurveFit<3,2>;
pub type QuinticSplineFit2D = ParameterSplineCurveFit<5,2>;
pub type LinearSplineFit3D = ParameterSplineCurveFit<1,3>;
pub type CubicSplineFit3D = ParameterSplineCurveFit<3,3>;
pub type QuinticSplineFit3D = ParameterSplineCurveFit<5,3>;

// Closed Periodic Curve Fits
pub type ClosedLinearSplineFit2D = ClosedParameterSplineCurveFit<1,2>;
pub type ClosedCubicSplineFit2D = ClosedParameterSplineCurveFit<3,2>;
pub type ClosedQuinticSplineFit2D = ClosedParameterSplineCurveFit<5,2>;
pub type ClosedLinearSplineFit3D = ClosedParameterSplineCurveFit<1,3>;
pub type ClosedCubicSplineFit3D = ClosedParameterSplineCurveFit<3,3>;
pub type ClosedQuinticSplineFit3D = ClosedParameterSplineCurveFit<5,3>;

#[derive(Debug, Clone)]
pub struct FitError(i32);

impl FitError {
    pub(crate) fn new(ierr: i32) -> Self { Self(ierr) }
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.0 {
            // Dierckx
           -2 => write!(f, "normal return for weighted least squares spline, fp upper bound for smoothing factor"),
           -1 => write!(f, "normal return for interpolating spline"),
            0 => write!(f, "normal return"),
            1 => write!(f, "out of storage space; nest too small (m/2); or s too small"),
            2 => write!(f, "smoothing spline error, s too small"),
            3 => write!(f, "reached iteration limit (20) for finding smoothing spline; s too small"),
           10 => write!(f, "invalid input data; check if -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)"),
           // this library
           200 => write!(f, "N should be between 1 and 10"),
           201 => write!(f, "need at least 2 parameter values"),
           202 => write!(f, "incorrect size of coordinate array xn"),
           203 => write!(f, "wrong size for weights array"),
           204 => write!(f, "too many derivative constraints supplied"),
           205 => write!(f, "cardinal spline spacing too large: select smaller interval"),
           206 => write!(f, "smoothing_spline not converged"),
           207 => write!(f, "failed to initialize smoothing_spline"),
           208 => write!(f, "K should be 1, 3 or 5"),
            _ => write!(f, "unknown error"),
        }
    }
}

impl error::Error for FitError {}
