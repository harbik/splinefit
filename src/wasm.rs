// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

//! WebAssembly bindings for splinefit.
//!
//! Enabled by the `wasm` feature flag.  Build with:
//!
//! ```sh
//! wasm-pack build --features wasm
//! ```
//!
//! Exposes [`CubicSpline`] — a thin wrapper around the Rust `SplineCurve<3, 1>` —
//! to JavaScript.

use wasm_bindgen::prelude::*;
use spliny::SplineCurve;

use crate::{CubicSplineFit, evaluate, ops};

/// A cubic B-spline curve fitted to 1-D data.
///
/// Construct via `CubicSpline.smoothing(x, y, rms)`, `CubicSpline.interpolating(x, y)`,
/// or `CubicSpline.cardinal(x, y, dt)`.  Then call `.evaluate()`, `.integral()`,
/// `.roots()`, or inspect `.knots()` and `.coefficients()`.
#[wasm_bindgen]
pub struct CubicSpline {
    inner: SplineCurve<3, 1>,
}

#[wasm_bindgen]
impl CubicSpline {
    /// Fit a smoothing cubic spline.
    ///
    /// `x` and `y` must have the same length (≥ 4) with `x` strictly increasing.
    /// `rms` controls the trade-off between smoothness and closeness of fit:
    /// smaller values produce more knots and a closer fit.
    pub fn smoothing(x: &[f64], y: &[f64], rms: f64) -> Result<CubicSpline, JsError> {
        let s = CubicSplineFit::new(x.to_vec(), y.to_vec())
            .smoothing_spline(rms)
            .map_err(|e| JsError::new(&e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Fit an interpolating cubic spline that passes through every data point.
    ///
    /// `x` and `y` must have the same length (≥ 4) with `x` strictly increasing.
    pub fn interpolating(x: &[f64], y: &[f64]) -> Result<CubicSpline, JsError> {
        let s = CubicSplineFit::new(x.to_vec(), y.to_vec())
            .interpolating_spline()
            .map_err(|e| JsError::new(&e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Fit a cardinal cubic spline on an equidistant knot grid with spacing `dt`.
    ///
    /// `x` and `y` must have the same length (≥ 4) with `x` strictly increasing.
    pub fn cardinal(x: &[f64], y: &[f64], dt: f64) -> Result<CubicSpline, JsError> {
        let s = CubicSplineFit::new(x.to_vec(), y.to_vec())
            .cardinal_spline(dt)
            .map_err(|e| JsError::new(&e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Evaluate the spline at each value in `x`, returning the fitted values.
    pub fn evaluate(&self, x: &[f64]) -> Result<Vec<f64>, JsError> {
        evaluate::evaluate(&self.inner, x)
            .map_err(|e| JsError::new(&e.to_string()))
    }

    /// Compute the definite integral of the spline over `[a, b]`.
    pub fn integral(&self, a: f64, b: f64) -> f64 {
        ops::integral(&self.inner, a, b)
    }

    /// Find all interior zeros of the spline, returned in ascending order.
    pub fn roots(&self) -> Result<Vec<f64>, JsError> {
        ops::roots(&self.inner)
            .map_err(|e| JsError::new(&e.to_string()))
    }

    /// Return the knot vector.
    pub fn knots(&self) -> Vec<f64> {
        self.inner.t.clone()
    }

    /// Return the B-spline coefficients.
    pub fn coefficients(&self) -> Vec<f64> {
        self.inner.c.clone()
    }

    /// Return the number of knots.
    #[wasm_bindgen(getter)]
    pub fn num_knots(&self) -> usize {
        self.inner.t.len()
    }
}
