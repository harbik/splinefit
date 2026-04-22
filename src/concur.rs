//! Builder for fitting a parametric B-spline curve through multi-dimensional data.

use std::iter::repeat;
use crate::dierckx::concur_;
use super::{FitError};
use crate::Result;
use spliny::SplineCurve;


/// Builder for fitting a parametric B-spline curve to multi-dimensional `(u, x₁, …, xₙ)` data.
///
/// The curve is parameterized by `u` and lives in an `N`-dimensional space.  The degree `K`
/// must be odd (1, 3, or 5) and `N` must be between 1 and 10.
/// Use the type aliases such as [`CubicSplineFit2D`] rather than naming the generics directly.
///
/// Coordinates are interleaved in a flat `xn` vector: `[x₀, y₀, x₁, y₁, …]` for `N = 2`.
///
/// Optionally call [`weights`](Self::weights), [`begin_constraints`](Self::begin_constraints),
/// or [`end_constraints`](Self::end_constraints) before calling a fit method.
#[derive(Clone)]
pub struct ParameterSplineCurveFit<const K:usize, const N:usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>,    // weight factors,
    xb: Vec<f64>,
    xe: Vec<f64>,

    t: Vec<f64>,
    c: Vec<f64>,
    e_rms: Option<f64>,
    n: i32,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
    xx: Vec<f64>,
    cp: Vec<f64>,
    ib: i32,
    ie: i32,
    m: i32,
    mx: i32,
    nest: i32,
    k: i32,
    idim: i32,
}


impl<const K:usize, const N:usize> ParameterSplineCurveFit<K, N> {

    /// Create a new fit builder from a parameter vector `u` and interleaved coordinate vector `xn`.
    ///
    /// `u` must be strictly increasing with at least 2 elements.
    /// `xn` must have length `N * u.len()`, with coordinates laid out as
    /// `[x₀, y₀, …, x₁, y₁, …]` (all dimensions of point 0 first, then point 1, etc.).
    pub fn new(
        u: Vec<f64>,
        xn: Vec<f64>,
    ) -> Result<Self> {

        let k = K as i32;
        if ![1,3,5].contains(&(k as i32)) { return Err(FitError::new(208).into()) };
        let idim = if (1..=10).contains(&N) { N as i32 } else {
                return Err(FitError::new(200).into())
            };
        let m = u.len() as i32; // number of coordinates
        if m<2 {return Err(FitError::new(201).into())};
        let mx = m * idim;
        if xn.len() as i32 != mx { return Err(FitError::new(202).into())}
        let w_vec = vec![1.0;m as usize];

        let xb = Vec::new();
        let ib = 0;

        let xe = Vec::new();
        let ie = 0;

        let nest = m+k+1 + 2*(k-1);
        let n = nest;  // length of tc
        let t_vec = vec![0.0; nest as usize];
        let c_vec = vec![0.0; (nest * idim) as usize];

        let iwrk_vec = vec![0i32; nest as usize];

        let wrk_vec = vec![0f64; (m*(k+1)+nest*(6+idim+3*k)) as usize];
        let xx_vec = vec![0.0; (idim*m) as usize];
        let cp_vec = vec![0.0; (2 * (k+1) * idim) as usize];

        Ok(Self { u, xn, w: w_vec, xb, xe, t: t_vec, c: c_vec, wrk: wrk_vec, iwrk: iwrk_vec, xx: xx_vec, cp: cp_vec, ib, ie, m, mx, nest, k, idim, n, e_rms: None})

    }

    /// Return the parameter values `u` used in the fit.
    pub fn u(&self) -> Vec<f64> {
        self.u[..self.m as usize].to_vec()
    }

    /// Return the interleaved coordinate data `xn` used in the fit.
    pub fn xn(&self) -> Vec<f64> {
        self.xn[..(self.m * self.idim) as usize].to_vec()
    }

    /// Override the per-point weights. `weights` must have the same length as `u`.
    pub fn weights(mut self, weights: Vec<f64>) -> Result<Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(FitError::new(203).into())
        }
    }

    /// Set derivative constraints at the start of the curve.
    ///
    /// `ub` is an array of `D` derivative vectors, one per order starting from the first
    /// derivative.  `D` may be at most `(K + 1) / 2`.
    pub fn begin_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xb = ub.iter().flatten().cloned().collect();
            self.ib = D as i32 -1;
            Ok(self)
        } else {
            Err(FitError::new(204).into())
        }
    }

    /// Set derivative constraints at the end of the curve.
    ///
    /// Mirror of [`begin_constraints`](Self::begin_constraints) applied at the terminal knot.
    pub fn end_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xe = ub.iter().flatten().cloned().collect();
            self.ie = D as i32 -1;
            Ok(self)
        } else {
            Err(FitError::new(204).into())
        }
    }

    fn concur(&mut self, iopt:i32, e_rms:Option<f64>, knots: Option<Vec<f64>>) ->  i32 {
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            self.m as f64 * e.powi(2)
        } else {
            0.0
        };

        let nb = self.xb.len() as i32;
        let ne = self.xe.len() as i32;
        let np = self.cp.len() as i32;
        let nc = self.c.len() as i32;
        let lwrk = self.wrk.len() as i32;

        if let Some(knots) = knots {
            self.n = knots.len() as i32;
            self.t = knots;
            self.t.resize(self.nest as usize, 0.0);
        }
        let mut ierr = 0;
        unsafe {
            concur_(
                &iopt,
                &self.idim,
                &self.m,
                self.u.as_ptr(),
                &self.mx,
                self.xn.as_ptr(),
                self.xx.as_mut_ptr(),
                self.w.as_ptr(),
                &self.ib,
                self.xb.as_ptr(),
                &nb,
                &self.ie,
                self.xe.as_ptr(),
                &ne,
                &self.k,
                &s,
                &self.nest,
                &mut self.n,
                self.t.as_mut_ptr(),
                &nc,
                self.c.as_mut_ptr(),
                &np,
                self.cp.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        self.e_rms = Some((fp/self.m as f64).sqrt());
        ierr
    }


    /// Fit a weighted least-squares spline on a fixed equidistant knot grid.
    ///
    /// Knots are spaced `dt` apart along `u`, aligned to integer multiples of `dt`.
    pub fn cardinal_spline(mut self, dt:f64) -> Result<SplineCurve<K,N>>{
        let m = self.u.len();
        let tb = ((self.u[0]*(1.0+f64::EPSILON))/dt).ceil() * dt;
        let te = ((self.u[m-1]*(1.0-f64::EPSILON))/dt).floor() * dt;
        let n = ((te - tb)/dt).round() as usize;
        if n == 0 { return Err(FitError::new(205).into()) };
        let mut t: Vec<f64> = Vec::with_capacity(n + 2 * (K + 1) + 1);
        t.extend( repeat(tb).take(K+1) // begin padding, needed for spline evaluation
            .chain(
                repeat(dt).scan(tb + dt,
                    |s, dx|{
                        let t=*s;
                        *s+=dx;
                        if t<te {
                            Some(t)
                        } else {
                            None
                        }
                    })
            )
            .chain(
                repeat(te).take(K+1) // end padding
            )
        );

        let ierr = self.concur(-1, Some(0.0),Some(t));
        if ierr <= 0 {
            Ok(self.into())
        } else {
            Err(FitError::new(ierr).into())
        }
    }

    /// Fit a spline that passes exactly through all data points.
    pub fn interpolating_spline(mut self) -> Result<SplineCurve<K,N>> {
        let ierr = self.concur(0, Some(0.0),None);
        if ierr<=0  {
            Ok(self.into())
        } else {
            Err(FitError::new(ierr).into())
        }
    }

    /// Fit a smoothing spline with the fewest knots whose RMS error is ≤ `rms`.
    pub fn smoothing_spline(mut self, rms: f64) -> Result<SplineCurve<K,N>>{
        let ierr= self.concur(0, Some(rms), None);
        if ierr>0 {
            Err(FitError::new(ierr).into())
        } else {
            Ok(self.into())
        }
    }

    /// Iteratively refine a smoothing spline by progressively tightening the RMS target.
    ///
    /// Starts at `rms_start`, multiplies the target by `rms_scale_ratio` each step (use a value
    /// < 1.0 to tighten), and stops when `converged` returns `true` or `n_iter` steps are reached.
    ///
    /// The `converged` callback receives `(n_knots, added_knots, rms, rms_improvement)`.
    /// `n_iter` defaults to 40 if `None`.
    pub fn smoothing_spline_optimize(mut self,
            rms_start: f64,
            rms_scale_ratio: f64,
            converged: impl Fn(i32, i32, f64, f64) -> bool,
            n_iter: Option<usize>,
        ) -> Result<SplineCurve<K,N>>{
        let n_iter = n_iter.unwrap_or(40);
        let ierr= self.concur(0, Some(rms_start), None);
        if ierr>0 {
            return Err(FitError::new(ierr).into())
        }
        let mut rms = self.e_rms.unwrap();
        let mut n_prev;
        let mut rms_prev;
        for _ in 0..n_iter {
            n_prev = self.n;
            rms_prev = rms;
            let ierr= self.concur(1, Some(rms * rms_scale_ratio), None);
            rms = self.e_rms.unwrap();
            if ierr>0 {
                return Err(FitError::new(ierr).into())
            }
            if converged(self.n, self.n-n_prev, rms, rms_prev - rms) {
                let ierr = self.concur(0, Some(rms_prev), None);
                if ierr>0 {
                    return Err(FitError::new(ierr).into())
                } else {
                    return Ok(self.into())
                }
            };
        }
        Err(FitError::new(206).into())
    }

}


impl<const K:usize, const N:usize> From<ParameterSplineCurveFit<K,N>> for SplineCurve<K,N> {
    fn from(mut sp: ParameterSplineCurveFit<K,N>) -> Self {
        sp.t.truncate(sp.n as usize);
        sp.t.shrink_to_fit();

        sp.c.truncate((sp.n*sp.idim) as usize);

        for dim in 0..sp.idim as usize {
            let ib = (dim+1) * (sp.n-sp.k-1) as usize;
            let ie = ib + sp.k as usize + 1;
            sp.c.drain(ib..ie);
        }
        sp.c.shrink_to_fit();

        Self::new(
            sp.t,
            sp.c,
        )
    }
}
