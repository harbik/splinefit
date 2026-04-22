//! Builder for fitting a closed (periodic) parametric B-spline curve through multi-dimensional data.

use super::FitError;
use crate::Result;
use crate::dierckx::clocur_;
use spliny::SplineCurve;

/// Builder for fitting a closed (periodic) parametric B-spline curve to multi-dimensional data.
///
/// Like [`ParameterSplineCurveFit`](crate::concur::ParameterSplineCurveFit) but for curves that
/// form a closed loop.  The first and last data points **must coincide** to within 1e-10.
///
/// Degree `K` must be odd (1, 3, or 5) and `N` must be between 1 and 10.
/// Use the type aliases such as [`ClosedCubicSplineFit2D`] rather than naming the generics directly.
#[derive(Clone)]
pub struct ClosedParameterSplineCurveFit<const K: usize, const N: usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>, // weight factors
    ipar: i32,

    t: Vec<f64>,
    c: Vec<f64>,
    e_rms: Option<f64>,
    n: i32,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries,
    m: i32,
    mx: i32,
    nest: i32,
    k: i32,
    idim: i32,
}

impl<const K: usize, const N: usize> ClosedParameterSplineCurveFit<K, N> {
    /// Create a new fit builder from a parameter vector `u` and interleaved coordinate vector `xn`.
    ///
    /// `u` must be strictly increasing with at least 2 elements.
    /// `xn` must have length `N * u.len()`, laid out as `[x₀, y₀, …, x₁, y₁, …]`.
    /// The first and last points must coincide (to within 1e-10) to form a closed loop.
    pub fn new(u: Vec<f64>, xn: Vec<f64>) -> Result<Self> {
        let k = K as i32;
        if ![1, 3, 5].contains(&k) {
            return Err(FitError::new(208).into());
        };
        let idim = if (1..=10).contains(&N) {
            N as i32
        } else {
            return Err(FitError::new(200).into());
        };
        let m = u.len() as i32;
        if m < 2 {
            return Err(FitError::new(201).into());
        };
        let mx = m * idim;
        if xn.len() as i32 != mx {
            return Err(FitError::new(202).into());
        }

        // Validate that first and last points coincide
        let n_dim = N;
        let m_usize = m as usize;
        for d in 0..n_dim {
            let first = xn[d];
            let last = xn[(m_usize - 1) * n_dim + d];
            if (first - last).abs() > 1e-10 {
                return Err(
                    format!("clocur requires first and last points to coincide, but dimension {d} differs: {first} vs {last}").into()
                );
            }
        }

        let ipar = 1; // user-supplied parameter values
        let w_vec = vec![1.0; m as usize];

        let nest = m + 2 * k;
        let n = 0;
        let t_vec = vec![0.0; nest as usize];
        let c_vec = vec![0.0; (nest * idim) as usize];

        let iwrk_vec = vec![0i32; nest as usize];

        let lwrk = m * (k + 1) + nest * (7 + idim + 5 * k);
        let wrk_vec = vec![0f64; lwrk as usize];

        Ok(Self {
            u,
            xn,
            w: w_vec,
            ipar,
            t: t_vec,
            c: c_vec,
            wrk: wrk_vec,
            iwrk: iwrk_vec,
            m,
            mx,
            nest,
            k,
            idim,
            n,
            e_rms: None,
        })
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

    fn clocur(&mut self, iopt: i32, e_rms: Option<f64>, knots: Option<Vec<f64>>) -> i32 {
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            self.m as f64 * e.powi(2)
        } else {
            0.0
        };

        let nc = self.c.len() as i32;
        let lwrk = self.wrk.len() as i32;

        if let Some(knots) = knots {
            self.n = knots.len() as i32;
            self.t = knots;
            self.t.resize(self.nest as usize, 0.0);
        }
        let mut ierr = 0;
        unsafe {
            clocur_(
                &iopt,
                &self.ipar,
                &self.idim,
                &self.m,
                self.u.as_mut_ptr(),
                &self.mx,
                self.xn.as_ptr(),
                self.w.as_ptr(),
                &self.k,
                &s,
                &self.nest,
                &mut self.n,
                self.t.as_mut_ptr(),
                &nc,
                self.c.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        self.e_rms = Some((fp / self.m as f64).sqrt());
        ierr
    }

    /// Fit a closed spline that passes exactly through all data points.
    pub fn interpolating_spline(mut self) -> Result<SplineCurve<K, N>> {
        let ierr = self.clocur(0, Some(0.0), None);
        if ierr <= 0 {
            Ok(self.into())
        } else {
            Err(FitError::new(ierr).into())
        }
    }

    /// Fit a closed smoothing spline with the fewest knots whose RMS error is ≤ `rms`.
    pub fn smoothing_spline(mut self, rms: f64) -> Result<SplineCurve<K, N>> {
        let ierr = self.clocur(0, Some(rms), None);
        if ierr > 0 {
            Err(FitError::new(ierr).into())
        } else {
            Ok(self.into())
        }
    }

    /// Iteratively refine a closed smoothing spline by progressively tightening the RMS target.
    ///
    /// Starts at `rms_start`, multiplies the target by `rms_scale_ratio` each step (use a value
    /// < 1.0 to tighten), and stops when `converged` returns `true` or `n_iter` steps are reached.
    ///
    /// The `converged` callback receives `(n_knots, added_knots, rms, rms_improvement)`.
    /// `n_iter` defaults to 40 if `None`.
    pub fn smoothing_spline_optimize(
        mut self,
        rms_start: f64,
        rms_scale_ratio: f64,
        converged: impl Fn(i32, i32, f64, f64) -> bool,
        n_iter: Option<usize>,
    ) -> Result<SplineCurve<K, N>> {
        let n_iter = n_iter.unwrap_or(40);
        let ierr = self.clocur(0, Some(rms_start), None);
        if ierr > 0 {
            return Err(FitError::new(ierr).into());
        }
        let mut rms = self.e_rms.unwrap();
        let mut n_prev;
        let mut rms_prev;
        for _ in 0..n_iter {
            n_prev = self.n;
            rms_prev = rms;
            let ierr = self.clocur(1, Some(rms * rms_scale_ratio), None);
            rms = self.e_rms.unwrap();
            if ierr > 0 {
                return Err(FitError::new(ierr).into());
            }
            if converged(self.n, self.n - n_prev, rms, rms_prev - rms) {
                let ierr = self.clocur(0, Some(rms_prev), None);
                if ierr > 0 {
                    return Err(FitError::new(ierr).into());
                } else {
                    return Ok(self.into());
                }
            };
        }
        Err(FitError::new(206).into())
    }
}

impl<const K: usize, const N: usize> From<ClosedParameterSplineCurveFit<K, N>>
    for SplineCurve<K, N>
{
    fn from(mut sp: ClosedParameterSplineCurveFit<K, N>) -> Self {
        sp.t.truncate(sp.n as usize);
        sp.t.shrink_to_fit();

        sp.c.truncate((sp.n * sp.idim) as usize);

        for dim in 0..sp.idim as usize {
            let ib = (dim + 1) * (sp.n - sp.k - 1) as usize;
            let ie = ib + sp.k as usize + 1;
            sp.c.drain(ib..ie);
        }
        sp.c.shrink_to_fit();

        Self::new(sp.t, sp.c)
    }
}
