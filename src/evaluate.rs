// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

use super::{FitError, Result};
use spliny::SplineCurve;
use crate::dierckx::{splev_, curev_};


/// Evaluate a [`SplineCurve`] at each value in `x` (or `u` for parametric curves).
///
/// For 1-D splines (`N = 1`) this calls the Dierckx `splev` routine.
/// For parametric curves (`N > 1`) it calls `curev`, and the output is the
/// interleaved coordinate vector `[x₀, y₀, …, x₁, y₁, …]` of length `N * x.len()`.
///
/// Returns an error if the Dierckx routine signals an invalid argument.
pub fn evaluate<const K: usize, const N: usize>(s: &SplineCurve<K,N>, x: &[f64]) -> Result<Vec<f64>> {
    let (ierr, y)  =
        match N {
            1  => splev(s, x),
            _ => curev(s, x),
        };
    if ierr<=0 {
        Ok(y)
    } else {
        Err(FitError::new(ierr).into())
    }
}

fn splev<const K: usize, const N: usize>(s: &SplineCurve<K,N>, x: &[f64]) -> (i32, Vec<f64>) {
    let k = K as i32;
    let m = x.len() as i32;
    let mut y_v = vec![0.0; m as usize];
    let n = s.t.len() as i32;
    let mut ierr = 0;
    unsafe {
        splev_(
            s.t.as_ptr(),
            &n,
            s.c.as_ptr(),
            &k,
            x.as_ptr(),
            y_v.as_mut_ptr(),
            &m,
            &mut ierr
        );
    }
    (ierr, y_v)
}

fn curev<const K: usize, const N: usize>(s: &SplineCurve<K,N>, u: &[f64]) -> (i32, Vec<f64>) {
    let k = K as i32;
    let idim = N as i32;
    let m = u.len() as i32;
    let mxy = m * idim;
    let mut xy = vec![0.0; mxy as usize];
    let n = s.t.len();
    let n_i32 = n as i32;
    let nk1 = n - K - 1;
    // curev_ uses stride n per dimension; SplineCurve stores only nk1 coefficients
    // per dimension (the k+1 trailing zeros are stripped by the From conversion).
    // Reconstruct the full n*N array with k+1 trailing zeros per dimension.
    let nc = n * N;
    let mut c_padded = vec![0.0f64; nc];
    for dim in 0..N {
        c_padded[dim * n .. dim * n + nk1]
            .copy_from_slice(&s.c[dim * nk1 .. (dim + 1) * nk1]);
    }
    let nc_i32 = nc as i32;
    let mut ierr = 0;
    unsafe {
        curev_(
            &idim,
            s.t.as_ptr(),
            &n_i32,
            c_padded.as_ptr(),
            &nc_i32,
            &k,
            u.as_ptr(),
            &m,
            xy.as_mut_ptr(),
            &mxy,
            &mut ierr
        );
    }
    (ierr, xy)
}
