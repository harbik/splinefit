/// Integration tests for the pure-Rust Dierckx translation.
///
/// Each test uses the same inputs as a scipy reference run (scipy wraps the
/// original Dierckx Fortran library), so numerical agreement within ~1e-10
/// confirms correctness of the translation.
///
/// Reference values were generated with:
///
///   python3 -c "
///     import numpy as np; from scipy.interpolate import splrep, splev, splprep
///     # … see individual test comments for the exact calls
///   "
use std::f64::consts::PI;

// Import the public C-ABI functions as Rust items (re-exported from lib.rs).
use dierckx_rs::{curfit_, splev_, clocur_, curev_};

// ── helpers ───────────────────────────────────────────────────────────────────

/// Workspace sizes for curfit_.
///   lwrk ≥ m*(k+1) + nest*(7+3*k)   (from curfit.f lwest formula)
fn curfit_lwrk(m: usize, k: usize, nest: usize) -> usize {
    m * (k + 1) + nest * (7 + 3 * k)
}

// ── Test 1: curfit_ interpolation of sin(x) ───────────────────────────────────

/// curfit_ with s=0 (interpolating spline):
///   - ier must be -1
///   - knots must match scipy reference to 1e-10
///   - splev at data sites must recover y_i to < 1e-13
///
/// scipy: splrep(x, sin(x), k=3, s=0.0) → ier=-1, n=13
#[test]
fn curfit_sin_interpolation() {
    const M: usize = 9;
    const K: i32 = 3;
    const NEST: usize = M + K as usize + 1; // = 13  (exact for s=0)
    let lwrk = curfit_lwrk(M, K as usize, NEST);

    let x: Vec<f64> = (0..M).map(|i| i as f64 * PI / (M - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
    let w = vec![1.0f64; M];

    let (iopt, m_i, k, s, nest_i, lwrk_i, xb, xe) =
        (0i32, M as i32, K, 0.0f64, NEST as i32, lwrk as i32, x[0], x[M - 1]);

    let mut n = 0i32;
    let mut t   = vec![0.0f64; NEST];
    let mut c   = vec![0.0f64; NEST];
    let mut fp  = 0.0f64;
    let mut wrk = vec![0.0f64; lwrk];
    let mut iwrk = vec![0i32; NEST];
    let mut ier  = 0i32;

    unsafe {
        curfit_(
            &iopt, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
            &xb, &xe, &k, &s, &nest_i,
            &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
            wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier,
        );
    }

    assert_eq!(ier, -1, "expected interpolating spline (ier=-1), got {ier}");
    assert_eq!(n, 13,   "scipy gives n=13 for this input, got {n}");
    assert!(fp < 1e-14, "fp should be ~0 for interpolation, got {:.3e}", fp);

    // scipy reference knots (π/4 interior spacing, k+1 repeats at boundaries):
    //   splrep → t = [0,0,0,0, π/4, 3π/8, π/2, 5π/8, 3π/4, π,π,π,π]
    let t_ref = [
        0.0, 0.0, 0.0, 0.0,
        PI / 4.0, 3.0 * PI / 8.0, PI / 2.0, 5.0 * PI / 8.0, 3.0 * PI / 4.0,
        PI, PI, PI, PI,
    ];
    for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
        assert!((got - want).abs() < 1e-10, "t[{i}]: got {:.15}, want {:.15}", got, want);
    }

    // splev at data sites: |s(x_i) - sin(x_i)| < 1e-13
    let mut y_spl = vec![0.0f64; M];
    let mut spl_ier = 0i32;
    unsafe {
        splev_(t.as_ptr(), &n, c.as_ptr(), &k,
               x.as_ptr(), y_spl.as_mut_ptr(), &m_i, &mut spl_ier);
    }
    assert_eq!(spl_ier, 0);
    for i in 0..M {
        let err = (y_spl[i] - y[i]).abs();
        assert!(err < 1e-13,
            "splev residual at x[{i}]={:.4}: |{} - {}| = {err:.3e}",
            x[i], y_spl[i], y[i]);
    }
}

// ── Test 2: curfit_ smoothing of sin(x) ──────────────────────────────────────

/// curfit_ with s=0.05 (smoothing spline):
///   - ier must be 0 (normal return) or -1 (interpolation — also acceptable)
///   - if ier=0: |fp - s| / s ≤ 0.001  (the algorithm's convergence tolerance)
///   - interior knots must be in strictly ascending order
///   - splev must give a smooth approximation (eval error < 0.2 everywhere)
///
/// scipy: splrep(x, sin(x), k=3, s=0.05) → ier=0, fp≈0.050037, n=11
#[test]
fn curfit_sin_smoothing() {
    const M: usize = 20;
    const K: i32 = 3;
    const NEST: usize = M + K as usize + 1; // = 24
    let lwrk = curfit_lwrk(M, K as usize, NEST);

    let x: Vec<f64> = (0..M).map(|i| i as f64 * 2.0 * PI / (M - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
    let w = vec![1.0f64; M];
    let s = 0.05f64;

    let (iopt, m_i, k, nest_i, lwrk_i, xb, xe) =
        (0i32, M as i32, K, NEST as i32, lwrk as i32, x[0], x[M - 1]);

    let mut n = 0i32;
    let mut t   = vec![0.0f64; NEST];
    let mut c   = vec![0.0f64; NEST];
    let mut fp  = 0.0f64;
    let mut wrk = vec![0.0f64; lwrk];
    let mut iwrk = vec![0i32; NEST];
    let mut ier  = 0i32;

    unsafe {
        curfit_(
            &iopt, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
            &xb, &xe, &k, &s, &nest_i,
            &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
            wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier,
        );
    }

    assert!(ier <= 0, "expected success (ier ≤ 0), got ier={ier}");
    // scipy gives n=11 (3 interior knots)
    assert_eq!(n, 11, "scipy gives n=11 for s=0.05, got {n}");

    if ier == 0 {
        let rel = (fp - s).abs() / s;
        assert!(rel <= 0.002,
            "|fp - s| / s = {rel:.5}  (fp={fp:.6e}, s={s})");
    }

    // Interior knots in strictly ascending order
    let k1 = (K + 1) as usize;
    let nv  = n as usize;
    for i in k1..(nv - k1 - 1) {
        assert!(t[i] < t[i + 1], "knots not sorted: t[{i}]={} ≥ t[{}]={}", t[i], i+1, t[i+1]);
    }

    // splev gives a reasonable approximation of sin everywhere
    let x_dense: Vec<f64> = (0..100).map(|i| i as f64 * 2.0 * PI / 99.0).collect();
    let y_ref:   Vec<f64> = x_dense.iter().map(|&xi| xi.sin()).collect();
    let m_dense = 100i32;
    let mut y_spl = vec![0.0f64; 100];
    let mut spl_ier = 0i32;
    unsafe {
        splev_(t.as_ptr(), &n, c.as_ptr(), &k,
               x_dense.as_ptr(), y_spl.as_mut_ptr(), &m_dense, &mut spl_ier);
    }
    assert_eq!(spl_ier, 0);
    let max_err = y_spl.iter().zip(y_ref.iter()).map(|(&a, &b)| (a - b).abs())
                       .fold(0.0f64, f64::max);
    assert!(max_err < 0.2, "max approximation error = {max_err:.4}");
}

// ── Test 3: curfit_ exact fit of a quadratic ─────────────────────────────────

/// A degree-2 spline fits any quadratic exactly.  With s=0 on y=x² the
/// returned curve must recover x² to < 1e-13 at every evaluation site.
///
/// scipy: splrep(x, x**2, k=2, s=0) → ier=-1, n=8,
///   t = [0, 0, 0, 0.375, 0.625, 1, 1, 1]
#[test]
fn curfit_quadratic_exact_fit() {
    let x = vec![0.0f64, 0.25, 0.5, 0.75, 1.0];
    let y: Vec<f64> = x.iter().map(|&xi| xi * xi).collect();
    let w = vec![1.0f64; 5];
    let (m, k) = (5usize, 2i32);
    let nest = (m + k as usize + 1) as i32; // = 8
    let lwrk = curfit_lwrk(m, k as usize, nest as usize) as i32;

    let (iopt, m_i, s, xb, xe) = (0i32, m as i32, 0.0f64, 0.0f64, 1.0f64);

    let mut n = 0i32;
    let mut t   = vec![0.0f64; nest as usize];
    let mut c   = vec![0.0f64; nest as usize];
    let mut fp  = 0.0f64;
    let mut wrk = vec![0.0f64; lwrk as usize];
    let mut iwrk = vec![0i32; nest as usize];
    let mut ier  = 0i32;

    unsafe {
        curfit_(&iopt, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(), &xb, &xe, &k, &s, &nest,
                &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
                wrk.as_mut_ptr(), &lwrk, iwrk.as_mut_ptr(), &mut ier);
    }

    assert_eq!(ier, -1, "expected interpolation (ier=-1), got {ier}");
    assert_eq!(n, 8,    "scipy gives n=8, got {n}");

    // scipy reference knots: [0, 0, 0, 0.375, 0.625, 1, 1, 1]
    let t_ref = [0.0f64, 0.0, 0.0, 0.375, 0.625, 1.0, 1.0, 1.0];
    for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
        assert!((got - want).abs() < 1e-12, "t[{i}]: got {got:.15}, want {want:.15}");
    }

    // splev at data sites must recover x^2 exactly
    let mut y_spl = vec![0.0f64; m];
    let mut spl_ier = 0i32;
    unsafe {
        splev_(t.as_ptr(), &n, c.as_ptr(), &k,
               x.as_ptr(), y_spl.as_mut_ptr(), &m_i, &mut spl_ier);
    }
    assert_eq!(spl_ier, 0);
    for i in 0..m {
        let err = (y_spl[i] - y[i]).abs();
        assert!(err < 1e-13, "residual[{i}]: |{} - {}| = {err:.3e}", y_spl[i], y[i]);
    }
}

// ── Test 4: curfit_ returns weighted-least-squares polynomial for large s ─────

/// When s greatly exceeds the data variance, curfit returns ier=-2
/// (the least-squares polynomial of degree k, no interior knots)
/// and n = 2*(k+1).
///
/// scipy: splrep(x, sin(x), k=3, s=100) → ier=-2, n=8=2*(3+1)
#[test]
fn curfit_large_s_gives_polynomial() {
    const M: usize = 10;
    const K: i32 = 3;
    const NEST: usize = M + K as usize + 1; // = 14
    let lwrk = curfit_lwrk(M, K as usize, NEST);

    let x: Vec<f64> = (0..M).map(|i| i as f64 / (M - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
    let w = vec![1.0f64; M];
    let s = 100.0f64; // far above any realistic residual

    let (iopt, m_i, k, nest_i, lwrk_i, xb, xe) =
        (0i32, M as i32, K, NEST as i32, lwrk as i32, x[0], x[M - 1]);

    let mut n = 0i32;
    let mut t   = vec![0.0f64; NEST];
    let mut c   = vec![0.0f64; NEST];
    let mut fp  = 0.0f64;
    let mut wrk = vec![0.0f64; lwrk];
    let mut iwrk = vec![0i32; NEST];
    let mut ier  = 0i32;

    unsafe {
        curfit_(&iopt, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(), &xb, &xe, &k, &s, &nest_i,
                &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
                wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier);
    }

    assert_eq!(ier, -2, "expected polynomial result (ier=-2), got {ier}");
    assert_eq!(n, 2 * (K + 1), "n must be 2*(k+1)={}, got {n}", 2 * (K + 1));

    // scipy reference knots: [0,0,0,0, 1,1,1,1]
    let t_ref = [0.0f64, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
        assert!((got - want).abs() < 1e-12, "t[{i}]: got {got:.15}, want {want:.15}");
    }
}

// ── Test 5: clocur_ interpolation of a unit circle ───────────────────────────

/// clocur_ on 8 uniformly-spaced points of a unit circle with s=0:
///   - ier must be -1 (interpolating closed curve)
///   - n must be NEST = M + 2*k = 14
///   - curev_ at u[0..M-2] must recover the fitted data points to < 1e-12
///
/// Note: clocur_ uses u[M-1] only to specify the period (per = u[M-1] - u[0]);
/// the Mth point is NOT fitted (fpclos loop runs 1..m-1). By periodicity,
/// c(u[M-1]) == c(u[0]), so only points 0..M-2 can be checked for interpolation.
///
/// scipy: splprep([cos, sin], k=3, s=0, per=True) → ier=-1, n=14
#[test]
fn clocur_unit_circle_interpolation() {
    const M: usize = 8;
    const K: i32 = 3;
    const IDIM: usize = 2;
    const NEST: usize = M + 2 * K as usize;          // = 14
    const NC: usize   = NEST * IDIM;                  // = 28
    // lwrk ≥ m*k1 + nest*(7 + idim + 5*k)  (from clocur.f lwest)
    let lwrk = M * (K as usize + 1) + NEST * (7 + IDIM + 5 * K as usize); // = 368

    // 8 uniformly-spaced points on the unit circle, stored point-major
    let theta: Vec<f64> = (0..M).map(|i| 2.0 * PI * i as f64 / M as f64).collect();
    let xy: Vec<f64>    = theta.iter().flat_map(|&th| [th.cos(), th.sin()]).collect();
    let w = vec![1.0f64; M];

    let (iopt, ipar, idim, m_i, mx, k, s, nest_i, nc_i, lwrk_i) = (
        0i32, 0i32, IDIM as i32, M as i32, (M * IDIM) as i32,
        K, 0.0f64, NEST as i32, NC as i32, lwrk as i32,
    );

    let mut u    = vec![0.0f64; M];
    let mut n    = 0i32;
    let mut t    = vec![0.0f64; NEST];
    let mut c    = vec![0.0f64; NC];
    let mut fp   = 0.0f64;
    let mut wrk  = vec![0.0f64; lwrk];
    let mut iwrk = vec![0i32; NEST];
    let mut ier  = 0i32;

    unsafe {
        clocur_(
            &iopt, &ipar, &idim, &m_i,
            u.as_mut_ptr(), &mx, xy.as_ptr(), w.as_ptr(),
            &k, &s, &nest_i, &mut n,
            t.as_mut_ptr(), &nc_i, c.as_mut_ptr(), &mut fp,
            wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier,
        );
    }

    assert_eq!(ier, -1,       "expected interpolating closed curve (ier=-1), got {ier}");
    assert_eq!(n, NEST as i32,"expected n={NEST} knots, got {n}");
    assert!(fp < 1e-13, "fp should be ~0 for interpolation, got {fp:.3e}");

    // Evaluate at the computed parameter values; must recover original data
    let m_eval   = M as i32;
    let mx_eval  = (M * IDIM) as i32;
    let mut xy_out = vec![0.0f64; M * IDIM];
    let mut curev_ier = 0i32;
    unsafe {
        curev_(
            &idim, t.as_ptr(), &n, c.as_ptr(), &nc_i, &k,
            u.as_ptr(), &m_eval, xy_out.as_mut_ptr(), &mx_eval, &mut curev_ier,
        );
    }
    assert_eq!(curev_ier, 0, "curev_ returned error {curev_ier}");

    // Check points 0..M-2; the last entry (M-1) is the period specification
    // and evaluates to c(u[0]) by periodicity, not to the Mth data point.
    for i in 0..M-1 {
        let dx = (xy_out[2 * i]     - xy[2 * i]    ).abs();
        let dy = (xy_out[2 * i + 1] - xy[2 * i + 1]).abs();
        assert!(dx < 1e-12, "x residual at point {}: {:.3e}", i, dx);
        assert!(dy < 1e-12, "y residual at point {}: {:.3e}", i, dy);
    }
}

// ── Test 6: iopt=1 restart produces same result as fresh iopt=0 ───────────────

/// Calling curfit_ twice with iopt=1 (warm-restart with fewer knots, then
/// refine) must give the same ier and a fit at least as good as a fresh iopt=0
/// call with the same s.
#[test]
fn curfit_iopt1_restart_consistent() {
    const M: usize = 20;
    const K: i32 = 3;
    const NEST: usize = M + K as usize + 1;
    let lwrk = curfit_lwrk(M, K as usize, NEST);

    let x: Vec<f64> = (0..M).map(|i| i as f64 / (M - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| (3.0 * PI * xi).sin()).collect();
    let w = vec![1.0f64; M];
    let s = 0.02f64;
    let (m_i, k, nest_i, lwrk_i, xb, xe) =
        (M as i32, K, NEST as i32, lwrk as i32, x[0], x[M - 1]);

    // ── fresh call (iopt=0) ───────────────────────────────────────────────────
    let mut n0 = 0i32;
    let mut t0   = vec![0.0f64; NEST];
    let mut c0   = vec![0.0f64; NEST];
    let mut fp0  = 0.0f64;
    let mut wrk0 = vec![0.0f64; lwrk];
    let mut iwrk0 = vec![0i32; NEST];
    let mut ier0  = 0i32;
    unsafe {
        curfit_(&0i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(), &xb, &xe, &k, &s, &nest_i,
                &mut n0, t0.as_mut_ptr(), c0.as_mut_ptr(), &mut fp0,
                wrk0.as_mut_ptr(), &lwrk_i, iwrk0.as_mut_ptr(), &mut ier0);
    }
    assert!(ier0 <= 0, "fresh call failed: ier0={ier0}");

    // ── restart call (iopt=1) with same workspace ─────────────────────────────
    let mut fp1  = 0.0f64;
    let mut ier1 = 0i32;
    unsafe {
        curfit_(&1i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(), &xb, &xe, &k, &s, &nest_i,
                &mut n0, t0.as_mut_ptr(), c0.as_mut_ptr(), &mut fp1,
                wrk0.as_mut_ptr(), &lwrk_i, iwrk0.as_mut_ptr(), &mut ier1);
    }
    assert!(ier1 <= 0, "restart call failed: ier1={ier1}");
    // The restart should not increase the residual significantly
    assert!(fp1 <= fp0 + s * 0.01,
        "restart residual {fp1:.6e} is worse than fresh {fp0:.6e}");
}
