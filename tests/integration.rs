// Low-level integration tests live in src/dierckx.rs (#[cfg(test)] mod tests).
// High-level API tests live here.
use splinefit::{CubicSplineFit, evaluate};

#[test]
fn cubic_smoothing_spline_sin() {
    use std::f64::consts::PI;
    let m = 20usize;
    let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();

    let spline = CubicSplineFit::new(x.clone(), y)
        .smoothing_spline(0.05)
        .expect("smoothing spline failed");

    let y_fit = evaluate::evaluate(&spline, &x).expect("evaluate failed");
    let max_err = y_fit.iter().zip(x.iter())
        .map(|(&yf, &xi)| (yf - xi.sin()).abs())
        .fold(0.0f64, f64::max);
    assert!(max_err < 0.15, "max error = {max_err:.4}");
}
