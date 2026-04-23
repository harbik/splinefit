// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

//! Post-fit operations on [`SplineCurve`](spliny::SplineCurve) values.
//!
//! These free functions mirror SciPy's `splint` / `sproot` / `insert` functional
//! interface and complement the [`evaluate`](crate::evaluate) module.
//!
//! | Function | SciPy equivalent | Description |
//! |---|---|---|
//! | [`integral`] | `splint(a, b, tck)` | Definite integral over `[a, b]` |
//! | [`roots`] | `sproot(tck)` | All zeros of a cubic 1-D spline |
//! | [`insert_knot`] | `insert(x, tck)` | Insert a knot, refining the representation |

use super::{FitError, Result};
use crate::dierckx::{insert_, splint_, sproot_};
use spliny::SplineCurve;

/// Compute the definite integral of a 1-D spline from `a` to `b`.
///
/// The integration range is clamped to the spline's domain if `a` or `b` lie outside it.
///
/// # Example
/// ```
/// use splinefit::{CubicSplineFit, ops};
///
/// let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
/// let y: Vec<f64> = x.iter().map(|&xi| xi * xi).collect(); // y = x²
/// let s = CubicSplineFit::new(x, y).interpolating_spline().unwrap();
///
/// let val = ops::integral(&s, 0.0, 3.0); // ∫₀³ x² dx = 9.0
/// assert!((val - 9.0).abs() < 1e-10);
/// ```
pub fn integral<const K: usize>(s: &SplineCurve<K, 1>, a: f64, b: f64) -> f64 {
    let n = s.t.len() as i32;
    let k = K as i32;
    // SplineCurve stores only nk1 = n - k - 1 coefficients; splint_ indexes into
    // c[0..n-1], so pad the stored coefficients with k+1 trailing zeros.
    let mut c_full = vec![0.0f64; s.t.len()];
    c_full[..s.c.len()].copy_from_slice(&s.c);
    let mut wrk = vec![0.0f64; s.t.len()];
    unsafe { splint_(s.t.as_ptr(), &n, c_full.as_ptr(), &k, &a, &b, wrk.as_mut_ptr()) }
}

/// Find all zeros of a cubic 1-D spline within its domain, returned in ascending order.
///
/// Only cubic (`K = 3`) 1-D splines are supported — this restriction is enforced at the
/// type level.  Returns an error if the knot vector is malformed (fewer than 8 knots or
/// non-strictly-increasing interior knots).
///
/// # Example
/// ```
/// use splinefit::{CubicSplineFit, ops};
/// use std::f64::consts::PI;
///
/// let m = 17usize;
/// let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
/// let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
/// let s = CubicSplineFit::new(x, y).interpolating_spline().unwrap();
///
/// let zeros = ops::roots(&s).unwrap();
/// assert!(zeros.iter().any(|&z| (z - PI).abs() < 1e-6));
/// ```
pub fn roots(s: &SplineCurve<3, 1>) -> Result<Vec<f64>> {
    let n = s.t.len() as i32;
    let mest = n; // safe upper bound on number of roots
    let mut zeros = vec![0.0f64; n as usize];
    let mut m = 0i32;
    let mut ier = 0i32;
    // SplineCurve<3,1> stores only nk1 = n - 4 coefficients; sproot_ indexes into
    // c[0..n-1], so pad the stored coefficients with k+1 = 4 trailing zeros.
    let mut c_full = vec![0.0f64; n as usize];
    c_full[..s.c.len()].copy_from_slice(&s.c);
    unsafe {
        sproot_(
            s.t.as_ptr(),
            &n,
            c_full.as_ptr(),
            zeros.as_mut_ptr(),
            &mest,
            &mut m,
            &mut ier,
        );
    }
    if ier != 0 {
        return Err(FitError::new(ier).into());
    }
    zeros.truncate(m as usize);
    Ok(zeros)
}

/// Insert a single knot at position `x` into a spline, returning a new spline that
/// represents the same function over a finer knot grid.
///
/// The new spline has one additional knot and (per dimension) one additional B-spline
/// coefficient, but evaluates identically to the original.  This is useful for
/// increasing the local resolution of a spline representation without refitting.
///
/// Works for any degree `K` and any output dimension `N`.
///
/// Returns an error if `x` lies outside the spline's interior domain or if the knot
/// cannot be inserted (e.g. the knot already has full multiplicity).
///
/// # Example
/// ```
/// use splinefit::{CubicSplineFit, evaluate, ops};
/// use std::f64::consts::PI;
///
/// let m = 20usize;
/// let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
/// let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
/// let s = CubicSplineFit::new(x.clone(), y).smoothing_spline(0.05).unwrap();
///
/// let s2 = ops::insert_knot(&s, PI).unwrap();
/// let y1 = evaluate::evaluate(&s,  &x).unwrap();
/// let y2 = evaluate::evaluate(&s2, &x).unwrap();
/// for (a, b) in y1.iter().zip(y2.iter()) {
///     assert!((a - b).abs() < 1e-13);
/// }
/// ```
pub fn insert_knot<const K: usize, const N: usize>(
    s: &SplineCurve<K, N>,
    x: f64,
) -> Result<SplineCurve<K, N>> {
    let n = s.t.len();
    let k1 = K + 1;
    let nk1 = n - k1;      // coefficients per dimension in the input
    let n_new = n + 1;
    let nk1_new = nk1 + 1; // coefficients per dimension in the output
    let nest = n_new as i32;
    let n_i32 = n as i32;
    let k_i32 = K as i32;
    let iopt = 0i32;

    let mut tt = vec![0.0f64; n_new];
    let mut nn = 0i32;
    let mut new_c = vec![0.0f64; N * nk1_new];

    for dim in 0..N {
        // Reconstruct the full coefficient array for this dimension:
        // [c[0..nk1], 0, 0, ..., 0]  (k+1 trailing zeros, as in the clamped representation).
        let mut c_full = vec![0.0f64; n];
        c_full[..nk1].copy_from_slice(&s.c[dim * nk1..(dim + 1) * nk1]);

        let mut cc = vec![0.0f64; n_new];
        let mut ier = 0i32;
        unsafe {
            insert_(
                &iopt,
                s.t.as_ptr(),
                &n_i32,
                c_full.as_ptr(),
                &k_i32,
                &x,
                tt.as_mut_ptr(),
                &mut nn,
                cc.as_mut_ptr(),
                &nest,
                &mut ier,
            );
        }
        if ier != 0 {
            return Err(FitError::new(ier).into());
        }
        new_c[dim * nk1_new..(dim + 1) * nk1_new].copy_from_slice(&cc[..nk1_new]);
    }

    tt.truncate(nn as usize);
    Ok(SplineCurve::new(tt, new_c))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CubicSplineFit, evaluate};
    use std::f64::consts::PI;

    #[test]
    fn integral_x_squared() {
        // ∫₀³ x² dx = 9.0  — interpolating spline through x² should be exact
        let x: Vec<f64> = (0..=9).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi * xi).collect();
        let s = CubicSplineFit::new(x, y).interpolating_spline().unwrap();
        let val = integral(&s, 0.0, 3.0);
        assert!((val - 9.0).abs() < 1e-10, "∫₀³ x² = {val:.15}");
    }

    #[test]
    fn integral_sin_full_period() {
        // ∫₀^2π sin(x) dx = 0  (one full period)
        let m = 20usize;
        let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let s = CubicSplineFit::new(x, y).smoothing_spline(0.01).unwrap();
        let val = integral(&s, 0.0, 2.0 * PI);
        assert!(val.abs() < 0.01, "∫₀^2π sin ≈ 0, got {val:.6}");
    }

    #[test]
    fn roots_sin_on_two_pi() {
        let m = 17usize;
        let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let s = CubicSplineFit::new(x, y).interpolating_spline().unwrap();
        let zeros = roots(&s).unwrap();
        // sproot_ searches the open interior (t[k+1], t[n-k-1]), so boundary zeros
        // at x=0 and x=2π are not returned.  The interior zero near π must be present.
        assert!(zeros.iter().any(|&z| (z - PI).abs() < 0.01), "no zero near π");
    }

    #[test]
    fn insert_knot_preserves_evaluation_1d() {
        let m = 20usize;
        let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let s = CubicSplineFit::new(x.clone(), y).smoothing_spline(0.05).unwrap();
        let s2 = insert_knot(&s, PI).unwrap();
        assert_eq!(s2.t.len(), s.t.len() + 1, "knot count should increase by 1");
        let y1 = evaluate::evaluate(&s,  &x).unwrap();
        let y2 = evaluate::evaluate(&s2, &x).unwrap();
        for (i, (a, b)) in y1.iter().zip(y2.iter()).enumerate() {
            assert!((a - b).abs() < 1e-13, "point {i}: {a} vs {b}");
        }
    }

    #[test]
    fn insert_knot_preserves_evaluation_2d() {
        use crate::CubicSplineFit2D;
        let m = 12usize;
        let u: Vec<f64> = (0..m).map(|i| i as f64 * 1.5 * PI / (m - 1) as f64).collect();
        let xn: Vec<f64> = u.iter().flat_map(|&ui| [ui.cos(), ui.sin()]).collect();
        let s = CubicSplineFit2D::new(u.clone(), xn).unwrap()
            .smoothing_spline(0.01).unwrap();
        let s2 = insert_knot(&s, u[m / 2]).unwrap();
        assert_eq!(s2.t.len(), s.t.len() + 1);
        let y1 = evaluate::evaluate(&s,  &u).unwrap();
        let y2 = evaluate::evaluate(&s2, &u).unwrap();
        for (i, (a, b)) in y1.iter().zip(y2.iter()).enumerate() {
            assert!((a - b).abs() < 1e-12, "coord {i}: {a} vs {b}");
        }
    }
}
