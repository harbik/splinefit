// Low-level integration tests live in src/dierckx.rs (#[cfg(test)] mod tests).
// High-level API tests live here.
use splinefit::{CubicSplineFit, evaluate};

#[test]
fn cardinal_spline_sin() {
    let m = 50usize;
    let x: Vec<f64> = (0..m).map(|i| i as f64 * 10.0 / (m - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();

    // Finer spacing → more knots → better fit
    let spline_fine = CubicSplineFit::new(x.clone(), y.clone()).cardinal_spline(0.5).unwrap();
    let y_fit = evaluate::evaluate(&spline_fine, &x).unwrap();
    let max_err = y.iter().zip(y_fit.iter()).map(|(a,b)| (a-b).abs()).fold(0.0f64, f64::max);
    assert!(max_err < 1e-3, "dt=0.5 max error = {max_err:.2e}, expected < 1e-3");

    // Coarser spacing → fewer knots → larger error
    let spline_coarse = CubicSplineFit::new(x.clone(), y.clone()).cardinal_spline(2.0).unwrap();
    let y_fit2 = evaluate::evaluate(&spline_coarse, &x).unwrap();
    let max_err2 = y.iter().zip(y_fit2.iter()).map(|(a,b)| (a-b).abs()).fold(0.0f64, f64::max);
    assert!(max_err2 < 0.1, "dt=2.0 max error = {max_err2:.2e}, expected < 0.1");
    assert!(max_err2 > max_err, "coarser knots should give larger error");
}

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
