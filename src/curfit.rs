use std::iter::repeat;
use crate::dierckx::curfit_;
use spliny::SplineCurve;
use super::FitError;
use crate::Result;

/// Builder for fitting a 1-D B-spline to a set of `(x, y)` data points.
///
/// The degree of the spline is controlled by the const generic `K` (must be odd: 1, 3, or 5).
/// Use the type aliases [`LinearSplineFit`](crate::LinearSplineFit),
/// [`CubicSplineFit`](crate::CubicSplineFit), or [`QuinticSplineFit`](crate::QuinticSplineFit)
/// rather than naming the generic directly.
///
/// Construct with [`SplineCurveFit::new`], optionally call [`set_weights`](Self::set_weights),
/// then call one of the three fit methods to obtain a [`SplineCurve`].
pub struct SplineCurveFit<const K: usize> {
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,

    t: Vec<f64>,
    c: Vec<f64>,
    n: i32,

    e_rms: Option<f64>,

    wrk: Vec<f64>,
    iwrk: Vec<i32>,
}

impl<const K: usize> SplineCurveFit<K> {

    /// Create a new fit builder from `x` and `y` data vectors.
    ///
    /// Both vectors must have the same length and at least `K + 1` elements.
    /// `x` values must be strictly increasing.
    /// Default weights of `1.0` are used for all points; call
    /// [`set_weights`](Self::set_weights) to override.
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let m = x.len();
        let w_vec = vec![1.0; m];
        assert!(y.len() == m);

        let nest = m * K + 1;
        let t_vec = vec![0.0; nest];
        let c_vec = vec![0.0; nest];
        let n = nest as i32;

        let iwrk_vec = vec![0i32; nest];

        let lwrk = m * (K + 1) + nest * (7 + 3 * K);
        let wrk_vec = vec![0f64; lwrk];

        Self { x, y, w: w_vec, t: t_vec, c: c_vec, n, wrk: wrk_vec, iwrk: iwrk_vec, e_rms: None }
    }

    /// Override the per-point weights used in the least-squares fit.
    ///
    /// `weights` must have the same length as the `x`/`y` vectors supplied to [`new`](Self::new).
    /// All weights must be positive.
    pub fn set_weights(mut self, weights: Vec<f64>) -> Result<Self> {
        if weights.len() == self.x.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(FitError::new(203).into())
        }
    }

    fn curfit(&mut self, iopt: i32, e_rms: Option<f64>, knots: Option<Vec<f64>>) -> i32 {
        let k = K as i32;
        let m = self.x.len() as i32;
        if let Some(knots) = knots {
            self.n = knots.len() as i32;
            self.t = knots;
        }
        let nest = self.t.len() as i32;
        let lwrk = self.wrk.len() as i32;
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms { m as f64 * e.powi(2) } else { 0.0 };
        let mut ierr = 0;

        unsafe {
            curfit_(
                &iopt, &m,
                self.x.as_ptr(), self.y.as_ptr(), self.w.as_ptr(),
                &self.x[0], &self.x[m as usize - 1],
                &k, &s, &nest, &mut self.n,
                self.t.as_mut_ptr(), self.c.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(), &lwrk, self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        // RMS residual: only meaningful when all weights are 1.0.
        // With non-unit weights fp is the weighted sum of squared residuals,
        // so this value is an approximation when weights differ.
        self.e_rms = Some((fp / m as f64).sqrt());
        ierr
    }

    /// Fit a weighted least-squares spline on a fixed equidistant knot grid.
    ///
    /// Knots are spaced `dt` apart and aligned to integer multiples of `dt`,
    /// covering the range of `x`.  This is useful when you want direct control
    /// over knot placement rather than letting the algorithm choose.
    pub fn cardinal_spline(mut self, dt: f64) -> Result<SplineCurve<K, 1>> {
        let m = self.x.len();
        let tb = (self.x[0] / dt).ceil() * dt;
        let te = (self.x[m - 1] / dt).floor() * dt;
        let n = ((te - tb) / dt).round() as usize;
        if n == 0 { return Err(FitError::new(205).into()) }

        let mut t: Vec<f64> = Vec::with_capacity(n + 2 * (K + 1) + 1);
        t.extend(
            repeat(tb).take(K + 1)
                .chain(repeat(dt).scan(tb + dt, |s, dx| {
                    let t = *s;
                    *s += dx;
                    if t < te { Some(t) } else { None }
                }))
                .chain(repeat(te).take(K + 1)),
        );

        let ierr = self.curfit(-1, Some(0.0), Some(t));
        if ierr <= 0 { Ok(self.into()) } else { Err(FitError::new(ierr).into()) }
    }

    /// Fit a spline that passes exactly through all data points.
    ///
    /// Equivalent to setting the smoothing factor `s = 0`.
    pub fn interpolating_spline(mut self) -> Result<SplineCurve<K, 1>> {
        let ierr = self.curfit(0, Some(0.0), None);
        if ierr <= 0 { Ok(self.into()) } else { Err(FitError::new(ierr).into()) }
    }

    /// Fit a smoothing spline with the fewest knots whose RMS error is ≤ `rms`.
    ///
    /// `rms` is the root-mean-square residual target in the same units as `y`.
    /// A smaller value produces a closer fit with more knots; a larger value
    /// produces a smoother curve with fewer knots.
    pub fn smoothing_spline(mut self, rms: f64) -> Result<SplineCurve<K, 1>> {
        let ierr = self.curfit(0, Some(rms), None);
        if ierr <= 0 { Ok(self.into()) } else { Err(FitError::new(ierr).into()) }
    }
}

impl<const K: usize> From<SplineCurveFit<K>> for SplineCurve<K, 1> {
    fn from(mut sp: SplineCurveFit<K>) -> Self {
        sp.t.truncate(sp.n as usize);
        sp.t.shrink_to_fit();
        sp.c.truncate(sp.n as usize - (K + 1));
        sp.c.shrink_to_fit();
        Self::new(sp.t, sp.c)
    }
}
