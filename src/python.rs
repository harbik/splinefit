// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

//! Python bindings for splinefit via PyO3.
//!
//! Enabled by the `python` feature flag.  Build with:
//!
//! ```sh
//! maturin develop --features python
//! ```
//!
//! Exposes [`CubicSpline`] to Python as a native extension module.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
use spliny::SplineCurve;

use crate::{CubicSplineFit, evaluate, ops};

/// Convert a PyReadonlyArray1 to Vec<f64>, returning a Python error on failure.
fn to_vec(arr: &PyReadonlyArray1<f64>) -> PyResult<Vec<f64>> {
    arr.to_vec()
        .map_err(|e| PyValueError::new_err(format!("failed to convert array: {e}")))
}

/// Validate and convert x/y arrays for fitting.
fn extract_xy(x: &PyReadonlyArray1<f64>, y: &PyReadonlyArray1<f64>) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let xv = to_vec(x)?;
    let yv = to_vec(y)?;
    if xv.len() != yv.len() {
        return Err(PyValueError::new_err(
            format!("x and y must have the same length, got {} and {}", xv.len(), yv.len())
        ));
    }
    if xv.len() < 4 {
        return Err(PyValueError::new_err(
            format!("need at least 4 data points, got {}", xv.len())
        ));
    }
    Ok((xv, yv))
}

/// A cubic B-spline curve fitted to 1-D data.
///
/// Construct via class methods:
///
///     s = CubicSpline.smoothing(x, y, rms=0.05)
///     s = CubicSpline.interpolating(x, y)
///     s = CubicSpline.cardinal(x, y, dt=0.5)
///
/// Then evaluate, integrate, or find roots:
///
///     y_fit = s.evaluate(x_new)
///     area  = s.integral(0.0, 3.0)
///     zeros = s.roots()
#[pyclass]
pub struct CubicSpline {
    inner: SplineCurve<3, 1>,
}

#[pymethods]
impl CubicSpline {
    /// Fit a smoothing cubic spline.
    ///
    /// Parameters
    /// ----------
    /// x : numpy.ndarray
    ///     Strictly increasing abscissae (at least 4 points).
    /// y : numpy.ndarray
    ///     Ordinate values, same length as `x`.
    /// rms : float
    ///     Target root-mean-square residual. Smaller values produce
    ///     more knots and a closer fit.
    ///
    /// Returns
    /// -------
    /// CubicSpline
    #[staticmethod]
    fn smoothing(x: PyReadonlyArray1<f64>, y: PyReadonlyArray1<f64>, rms: f64) -> PyResult<Self> {
        let (xv, yv) = extract_xy(&x, &y)?;
        let s = CubicSplineFit::new(xv, yv)
            .smoothing_spline(rms)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Fit an interpolating cubic spline through every data point.
    ///
    /// Parameters
    /// ----------
    /// x : numpy.ndarray
    ///     Strictly increasing abscissae (at least 4 points).
    /// y : numpy.ndarray
    ///     Ordinate values, same length as `x`.
    ///
    /// Returns
    /// -------
    /// CubicSpline
    #[staticmethod]
    fn interpolating(x: PyReadonlyArray1<f64>, y: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let (xv, yv) = extract_xy(&x, &y)?;
        let s = CubicSplineFit::new(xv, yv)
            .interpolating_spline()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Fit a cardinal cubic spline on an equidistant knot grid.
    ///
    /// Parameters
    /// ----------
    /// x : numpy.ndarray
    ///     Strictly increasing abscissae (at least 4 points).
    /// y : numpy.ndarray
    ///     Ordinate values, same length as `x`.
    /// dt : float
    ///     Knot spacing.
    ///
    /// Returns
    /// -------
    /// CubicSpline
    #[staticmethod]
    fn cardinal(x: PyReadonlyArray1<f64>, y: PyReadonlyArray1<f64>, dt: f64) -> PyResult<Self> {
        let (xv, yv) = extract_xy(&x, &y)?;
        let s = CubicSplineFit::new(xv, yv)
            .cardinal_spline(dt)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(CubicSpline { inner: s })
    }

    /// Evaluate the spline at the given points.
    ///
    /// Parameters
    /// ----------
    /// x : numpy.ndarray
    ///     Points at which to evaluate, within the spline's domain.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_slice = x.as_slice()
            .map_err(|e| PyValueError::new_err(format!("failed to read array: {e}")))?;
        let y = evaluate::evaluate(&self.inner, x_slice)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PyArray1::from_vec(py, y))
    }

    /// Compute the definite integral over [a, b].
    ///
    /// Parameters
    /// ----------
    /// a : float
    ///     Lower integration bound.
    /// b : float
    ///     Upper integration bound.
    ///
    /// Returns
    /// -------
    /// float
    fn integral(&self, a: f64, b: f64) -> f64 {
        ops::integral(&self.inner, a, b)
    }

    /// Find all interior zeros of the spline.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Zeros in ascending order. Boundary zeros may not be included.
    fn roots<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let z = ops::roots(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PyArray1::from_vec(py, z))
    }

    /// The knot vector.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    fn knots<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.t.clone())
    }

    /// The B-spline coefficients.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    fn coefficients<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.c.clone())
    }

    /// Number of knots.
    #[getter]
    fn num_knots(&self) -> usize {
        self.inner.t.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "CubicSpline(num_knots={}, domain=[{:.6}, {:.6}])",
            self.inner.t.len(),
            self.inner.t.first().unwrap_or(&0.0),
            self.inner.t.last().unwrap_or(&0.0),
        )
    }
}

/// splinefit — B-spline curve fitting
///
/// A pure-Rust B-spline fitting library exposed to Python. Backed by a
/// line-by-line translation of Paul Dierckx' Fortran FITPACK library
/// (the same engine behind SciPy's splrep/splev).
///
/// Example
/// -------
/// >>> import numpy as np
/// >>> from splinefit import CubicSpline
/// >>> x = np.linspace(0, 2 * np.pi, 50)
/// >>> y = np.sin(x)
/// >>> s = CubicSpline.smoothing(x, y, rms=0.05)
/// >>> y_fit = s.evaluate(x)
#[pymodule]
pub fn splinefit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CubicSpline>()?;
    Ok(())
}
