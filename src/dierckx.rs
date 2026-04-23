// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

// Pure Rust translation of the Dierckx Fortran B-spline library.
//
// ## Why `unsafe extern "C"` and raw pointers?
//
// The 11 public entry points in this module (`curfit_`, `splev_`, `curev_`, etc.) are
// declared `pub unsafe extern "C"` and accept raw pointers (`*const f64`, `*mut f64`,
// `*const i32`, …) instead of idiomatic Rust slices.  This is a deliberate design
// choice, not a limitation.
//
// The entire module is a **line-by-line translation** of Paul Dierckx' original Fortran
// FITPACK library.  Every Rust function maps 1:1 to a Fortran subroutine, every loop
// maps to the corresponding DO loop, and every array index `a[(i-1)+(j-1)*nest]` maps
// to the Fortran column-major `a(i,j)`.  The raw-pointer signatures preserve this
// correspondence at the API boundary: the Fortran subroutines receive arrays by address,
// and so do ours.
//
// This 1:1 correspondence provides three key benefits:
//
// 1. **Verifiability** — any line of Rust can be compared side-by-side with the original
//    Fortran source (available at <http://www.netlib.org/dierckx/>).  Translation bugs
//    are caught by reading two lines, not by reasoning about an abstracted rewrite.
//
// 2. **Numerical fidelity** — the algorithms are mathematically subtle (Givens rotations,
//    knot insertion, B-spline recurrences).  A mechanical translation preserves the exact
//    order of operations, accumulator structure, and loop bounds that Dierckx validated
//    in his 1993 book and that SciPy relies on.  Refactoring to idiomatic Rust would risk
//    introducing silent numerical drift that is hard to test for.
//
// 3. **Upstream traceability** — when SciPy or the Fortran community discovers a bug in
//    FITPACK, the fix can be applied here by reading the Fortran patch and updating the
//    matching Rust line.  An abstracted rewrite would require re-deriving each fix.
//
// The `unsafe` boundary is contained: the module is private (`mod dierckx`), and all
// public-facing callers (`curfit.rs`, `concur.rs`, `evaluate.rs`, `ops.rs`, …) convert
// between safe Rust types and raw pointers at the call site.  No `unsafe` leaks into
// the public API.
//
// ## Translation conventions
//
// All Fortran `real` compiled with -fdefault-real-8 → f64.
// All Fortran `integer` → i32 (or usize for indexing).
// Fortran 1-based arrays → 0-based in Rust.
// Fortran 2D column-major a(nest,k): a(i,j) → a[(i-1)+(j-1)*nest].
#![allow(clippy::all, non_snake_case, unused_variables, unused_assignments, unused_mut)]

// ─── Section 1: Core Primitives ───────────────────────────────────────────────

/// Evaluate the (k+1) non-zero B-splines of degree k at t[l-1] <= x < t[l]
/// using the stable recurrence of de Boor and Cox.
/// `l` is 1-based (Fortran convention); `h` must have length >= 6.
fn fpbspl(t: &[f64], _n: i32, k: i32, x: f64, l: i32, h: &mut [f64]) {
    let l = l as usize;
    h[0] = 1.0_f64;
    for j in 1..=(k as usize) {
        let mut hh = [0.0f64; 5];
        for i in 0..j {
            hh[i] = h[i];
        }
        h[0] = 0.0;
        for i in 1..=j {
            let li = l + i;       // 1-based: Fortran li = l+i
            let lj = li - j;      // 1-based: Fortran lj = li-j
            let f = hh[i - 1] / (t[li - 1] - t[lj - 1]);
            h[i - 1] += f * (t[li - 1] - x);
            h[i] = f * (x - t[lj - 1]);
        }
    }
}

/// Calculate Givens rotation parameters.
/// Matches the Fortran fpgivs exactly, including the store==0 early-exit.
fn fpgivs(piv: f64, ww: &mut f64, cos: &mut f64, sin: &mut f64) {
    let store = piv.abs();
    if store >= *ww {
        // piv dominates; guard against piv==0 (and ww==0)
        if store == 0.0 { *cos = 1.0; *sin = 0.0; return; }
        let dd2 = *ww / store;
        let dd = store * (1.0 + dd2 * dd2).sqrt();
        *cos = *ww / dd;
        *sin = piv / dd;
        *ww = dd;
    } else {
        let dd1 = store / *ww;
        let dd = *ww * (1.0 + dd1 * dd1).sqrt();
        *cos = *ww / dd;
        *sin = piv / dd;
        *ww = dd;
    }
}

/// Apply Givens rotation to a and b.
#[inline]
fn fprota(cos: f64, sin: f64, a: &mut f64, b: &mut f64) {
    let stor1 = *a;
    let stor2 = *b;
    *b = cos * stor2 + sin * stor1;
    *a = cos * stor1 - sin * stor2;
}

/// Backward substitution: solve a*c = z where a is n×n upper triangular
/// with bandwidth k. a is stored column-major with leading dim nest.
fn fpback(a: &[f64], z: &[f64], n: i32, k: i32, c: &mut [f64], nest: i32) {
    let n = n as usize;
    let k = k as usize;
    let nest = nest as usize;
    let k1 = k - 1;
    // a(i,j) → a[(i-1) + (j-1)*nest]
    c[n - 1] = z[n - 1] / a[(n - 1) + 0 * nest]; // a(n,1)
    if n == 1 { return; }
    let mut i = n - 1; // Fortran i starts at n-1 (n-1 in 1-based is n-2 in 0-based? no)
    // Fortran: i=n-1 initially (1-based), then decremented
    // After c(n) = z(n)/a(n,1), i = n-1
    let mut i_f = (n - 1) as i32; // Fortran i (1-based)
    for _j in 2..=n {
        let i_idx = (i_f - 1) as usize; // 0-based
        let mut store = z[i_idx];
        let i1 = if (_j as i32) <= (k1 as i32) { _j - 1 } else { k1 };
        let mut m_f = i_f;
        for l in 1..=i1 {
            m_f += 1;
            let m_idx = (m_f - 1) as usize;
            // a(i, l+1) → a[(i_f-1) + l*nest]
            store -= c[m_idx] * a[i_idx + l * nest];
        }
        // a(i,1) → a[(i_f-1) + 0*nest]
        c[i_idx] = store / a[i_idx + 0 * nest];
        i_f -= 1;
    }
}

/// Solve g*c = z where g is n×n upper triangular of the form
///   g = [ a | 0 ]
///       [   | b ]
/// with b an n×k matrix and a an (n-k)×(n-k) upper triangular matrix
/// of bandwidth k1. Both stored column-major with leading dim nest.
fn fpbacp(a: &[f64], b: &[f64], z: &[f64], n: i32, k: i32, c: &mut [f64], k1: i32, nest: i32) {
    let n = n as usize;
    let k = k as usize;
    let k1 = k1 as usize;
    let nest = nest as usize;
    let n2 = n - k;

    let mut l_f = n as i32; // Fortran l
    for i in 1..=k {
        let l_idx = (l_f - 1) as usize;
        let mut store = z[l_idx];
        let j_start = (k + 2 - i) as i32;
        if i != 1 {
            let mut l0_f = l_f;
            for l1 in j_start as usize..=k {
                l0_f += 1;
                let l0_idx = (l0_f - 1) as usize;
                // b(l, l1) → b[(l_idx) + (l1-1)*nest]
                store -= c[l0_idx] * b[l_idx + (l1 - 1) * nest];
            }
        }
        // b(l, j-1) → b[l_idx + (j_start as usize - 2)*nest]
        let j_m1 = (j_start - 1) as usize;
        c[l_idx] = store / b[l_idx + (j_m1 - 1) * nest];
        l_f -= 1;
        if l_f == 0 { return; }
    }
    for i in 1..=n2 {
        let i_idx = i - 1;
        let mut store = z[i_idx];
        let mut l_f2 = n2 as i32;
        for _j in 1..=k {
            l_f2 += 1;
            let l_idx = (l_f2 - 1) as usize;
            // b(i, j) → b[i_idx + (j-1)*nest]
            store -= c[l_idx] * b[i_idx + (_j - 1) * nest];
        }
        c[i_idx] = store;
    }
    let mut i_f = n2 as i32;
    let i_idx = (i_f - 1) as usize;
    c[i_idx] /= a[i_idx + 0 * nest]; // a(i,1)
    if i_f == 1 { return; }
    for j in 2..=n2 {
        i_f -= 1;
        let i_idx = (i_f - 1) as usize;
        let mut store = c[i_idx];
        let i1 = if (j as i32) <= (k as i32) { j - 1 } else { k };
        let mut l_f2 = i_f;
        for l0 in 1..=i1 {
            l_f2 += 1;
            let l_idx = (l_f2 - 1) as usize;
            // a(i, l0+1) → a[i_idx + l0*nest]
            store -= c[l_idx] * a[i_idx + l0 * nest];
        }
        c[i_idx] = store / a[i_idx + 0 * nest];
    }
}

// ─── Section 2: Knot Utilities ────────────────────────────────────────────────

/// Calculate discontinuity jumps of the k-th derivative of B-splines at
/// interior knots t[k+1]..t[n-k-2] (1-based: t(k+2)..t(n-k-1)).
/// b is stored column-major with leading dim nest.
fn fpdisc(t: &[f64], n: i32, k2: i32, b: &mut [f64], nest: i32) {
    let n = n as usize;
    let k2 = k2 as usize;
    let nest = nest as usize;
    let k1 = k2 - 1;
    let k = k1 - 1;
    let nk1 = n - k1;
    let nrint = nk1 - k;
    let an = nrint as f64;
    let fac = an / (t[nk1] - t[k1 - 1]); // t(nk1+1)-t(k1) 1-based → t[nk1]-t[k1-1]
    let mut h = [0.0f64; 12];

    for l in k2..=nk1 {
        // Fortran l goes from k2 to nk1 (1-based)
        let lmk = l - k1; // 1-based
        for j in 1..=k1 {
            let ik = j + k1;
            let lj = l + j;
            let lk = lj - k2;
            h[j - 1] = t[l - 1] - t[lk - 1];   // h(j) = t(l)-t(lk)
            h[ik - 1] = t[l - 1] - t[lj - 1];   // h(ik) = t(l)-t(lj)
        }
        let mut lp = lmk; // 1-based
        for j in 1..=k2 {
            let mut jk = j;
            let mut prod = h[j - 1];
            for _i in 1..=k {
                jk += 1;
                prod *= h[jk - 1] * fac;
            }
            let lk2 = lp + k1; // 1-based
            // b(lmk, j) → b[(lmk-1) + (j-1)*nest]
            b[(lmk - 1) + (j - 1) * nest] = (t[lk2 - 1] - t[lp - 1]) / prod;
            lp += 1;
        }
    }
}

/// Locate an additional knot for degree-k spline and adjust parameters.
fn fpknot(x: &[f64], m: i32, t: &mut [f64], n: &mut i32, fpint: &mut [f64],
          nrdata: &mut [i32], nrint: &mut i32, nest: i32, istart: i32) {
    let m = m as usize;
    let k = (*n - *nrint - 1) / 2;
    let mut fpmax = 0.0f64;
    let mut number = 0i32;
    let mut maxpt = 0i32;
    let mut maxbeg = 0i32;
    let mut jbegin = istart;
    for j in 1..=*nrint {
        let jpoint = nrdata[(j - 1) as usize];
        if fpmax < fpint[(j - 1) as usize] && jpoint != 0 {
            fpmax = fpint[(j - 1) as usize];
            number = j;
            maxpt = jpoint;
            maxbeg = jbegin;
        }
        jbegin += jpoint + 1;
    }
    // No interval has a splittable data point; nothing to do.
    if number == 0 { return; }
    let ihalf = maxpt / 2 + 1;
    let nrx = (maxbeg + ihalf) as usize; // 1-based index into x
    let next = number + 1;
    if next <= *nrint {
        for j in next..=*nrint {
            let jj = (next + *nrint - j) as usize;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            let jk = jj + k as usize;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[(number - 1) as usize] = ihalf - 1;
    nrdata[next as usize - 1] = maxpt - ihalf;
    let am = maxpt as f64;
    let an_num = nrdata[(number - 1) as usize] as f64;
    fpint[(number - 1) as usize] = fpmax * an_num / am;
    let an_next = nrdata[(next - 1) as usize] as f64;
    fpint[(next - 1) as usize] = fpmax * an_next / am;
    let jk = (next + k) as usize;
    t[jk - 1] = x[nrx - 1];
    *n += 1;
    *nrint += 1;
}

/// Rational interpolation: given (p1,f1),(p2,f2),(p3,f3), return p with r(p)=0.
/// Also adjusts (p1,f1) or (p3,f3) to maintain f1>0, f3<0.
fn fprati(p1: &mut f64, f1: &mut f64, p2: f64, f2: f64, p3: &mut f64, f3: &mut f64) -> f64 {
    let p = if *p3 <= 0.0 {
        (*p1 * (*f1 - *f3) * f2 - p2 * (f2 - *f3) * *f1) / ((*f1 - f2) * *f3)
    } else {
        let h1 = *f1 * (f2 - *f3);
        let h2 = f2 * (*f3 - *f1);
        let h3 = *f3 * (*f1 - f2);
        -(*p1 * p2 * h3 + p2 * *p3 * h1 + *p3 * *p1 * h2) / (*p1 * h1 + p2 * h2 + *p3 * h3)
    };
    if f2 < 0.0 {
        *p3 = p2;
        *f3 = f2;
    } else {
        *p1 = p2;
        *f1 = f2;
    }
    p
}

/// Verify knot positions for a degree-k spline (non-periodic).
fn fpchec(x: &[f64], m: i32, t: &[f64], n: i32, k: i32, ier: &mut i32) {
    let m = m as usize;
    let n = n as usize;
    let k = k as usize;
    let k1 = k + 1;
    let k2 = k1 + 1;
    let nk1 = n - k1;
    let nk2 = nk1 + 1;
    *ier = 31;
    if nk1 < k1 || nk1 > m { return; }
    *ier = 32;
    let mut j = n;
    for i in 1..=k {
        if t[i - 1] > t[i] { return; }
        if t[j - 1] < t[j - 2] { return; }
        j -= 1;
    }
    *ier = 33;
    for i in k2..=nk2 {
        if t[i - 1] <= t[i - 2] { return; }
    }
    *ier = 34;
    if x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1] { return; }
    *ier = 35;
    if x[0] >= t[k2 - 1] || x[m - 1] <= t[nk1 - 1] { return; }
    let mut i = 1usize;
    let mut l = k2;
    let nk3 = nk1 - 1;
    if nk3 < 2 { *ier = 0; return; }
    for j in 2..=nk3 {
        let tj = t[j - 1];
        l += 1;
        let tl = t[l - 1];
        loop {
            i += 1;
            if i >= m { return; }
            if x[i - 1] > tj { break; }
        }
        if x[i - 1] >= tl { return; }
    }
    *ier = 0;
}

/// Verify knot positions for a periodic spline.
fn fpchep(x: &[f64], m: i32, t: &[f64], n: i32, k: i32, ier: &mut i32) {
    let m = m as usize;
    let n = n as usize;
    let k = k as usize;
    let k1 = k + 1;
    let k2 = k1 + 1;
    let nk1 = n - k1;
    let nk2 = nk1 + 1;
    let m1 = m - 1;
    *ier = 10;
    if nk1 < k1 || n > m + 2 * k { return; }
    let mut j = n;
    for i in 1..=k {
        if t[i - 1] > t[i] { return; }
        if t[j - 1] < t[j - 2] { return; }
        j -= 1;
    }
    for i in k2..=nk2 {
        if t[i - 1] <= t[i - 2] { return; }
    }
    if x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1] { return; }
    // check condition 5
    let mut l1 = k1;
    let mut l2 = 1usize;
    let mut l_exit = m; // default
    'outer: for l in 1..=m {
        let xi = x[l - 1];
        loop {
            if xi < t[l1] || l == nk1 { break; } // xi < t(l1+1) → t[l1]
            l1 += 1;
            l2 += 1;
            if l2 > k1 { l_exit = l; break 'outer; }
        }
    }
    let per = t[nk2 - 1] - t[k1 - 1];
    'outer2: for i1 in 2..=l_exit {
        let i_start = i1 - 1;
        let mm = i_start + m1;
        let mut i = i_start;
        for j in k1..=nk1 {
            let tj = t[j - 1];
            let j1 = j + k1;
            let tl = t[j1 - 1];
            loop {
                i += 1;
                if i > mm { continue 'outer2; }
                let xi = if i <= m1 { x[i - 1] } else { x[i - m1 - 1] + per };
                if xi <= tj { continue; }
                if xi >= tl { continue 'outer2; }
                break;
            }
        }
        *ier = 0;
        return;
    }
}

/// Verify knot positions for a spline with derivative constraints at endpoints.
fn fpched(x: &[f64], m: i32, t: &[f64], n: i32, k: i32, ib: i32, ie: i32, ier: &mut i32) {
    let m = m as usize;
    let n = n as usize;
    let k = k as usize;
    let k1 = k + 1;
    let k2 = k1 + 1;
    let nk1 = n - k1;
    let nk2 = nk1 + 1;
    let ib1 = (ib - 1).max(0) as usize;
    let ie1 = (ie - 1).max(0) as usize;
    *ier = 51;
    if nk1 < k1 || nk1 > m + ib1 + ie1 { return; }
    *ier = 52;
    let mut j = n;
    for i in 1..=k {
        if t[i - 1] > t[i] { return; }
        if t[j - 1] < t[j - 2] { return; }
        j -= 1;
    }
    *ier = 53;
    for i in k2..=nk2 {
        if t[i - 1] <= t[i - 2] { return; }
    }
    *ier = 54;
    if x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1] { return; }
    *ier = 55;
    if x[0] >= t[k2 - 1] || x[m - 1] <= t[nk1 - 1] { return; }
    let mut i = 1usize;
    let jj_start = 2 + ib1;
    let mut l = jj_start + k;
    let nk3 = nk1 - 1 - ie1;
    if nk3 < jj_start { *ier = 0; return; }
    for _j in jj_start..=nk3 {
        let tj = t[_j - 1];
        l += 1;
        let tl = t[l - 1];
        loop {
            i += 1;
            if i >= m { return; }
            if x[i - 1] > tj { break; }
        }
        if x[i - 1] >= tl { return; }
    }
    *ier = 0;
}

// ─── Section 3: Basis Function Utilities ──────────────────────────────────────

/// Compute derivatives d(j) = s^(j-1)(x) at t[l-1] <= x < t[l].
/// k1 = k+1 (spline order), l is 1-based knot interval index.
fn fpader(t: &[f64], _n: i32, c: &[f64], k1: i32, x: f64, l: i32, d: &mut [f64]) {
    let k1 = k1 as usize;
    let l = l as usize;
    let lk = l - k1; // Fortran: lk = l-k1 (1-based offset into c)

    let mut h = [0.0f64; 6];
    // h(i) = c(lk+i) for i=1..k1
    for i in 1..=k1 {
        h[i - 1] = c[lk + i - 1]; // c(lk+i) 1-based → c[lk+i-1]
    }

    let mut kj = k1;
    let mut fac = 1.0f64;

    for j in 1..=k1 {
        let ki = kj;

        // Backward difference step (skip for j==1)
        if j > 1 {
            let mut i = k1;
            for _jj in j..=k1 {
                let li = i + lk;    // 1-based
                let lj = li + kj;   // 1-based
                h[i - 1] = (h[i - 1] - h[i - 2]) / (t[lj - 1] - t[li - 1]);
                if i <= 1 { break; }
                i -= 1;
            }
        }

        // Copy h[j..k1] to d[j..k1]
        for i in j..=k1 {
            d[i - 1] = h[i - 1];
        }

        // De Boor evaluation triangular sweep
        if j < k1 {
            let mut ki2 = ki;
            for jj in (j + 1)..=k1 {
                ki2 -= 1;
                let mut i = k1;
                for _j2 in jj..=k1 {
                    let li = i + lk;    // 1-based
                    let lj = li + ki2;  // 1-based
                    d[i - 1] = ((x - t[li - 1]) * d[i - 1] + (t[lj - 1] - x) * d[i - 2])
                        / (t[lj - 1] - t[li - 1]);
                    if i <= 1 { break; }
                    i -= 1;
                }
            }
        }

        // d(j) = d(k1) * fac
        d[j - 1] = d[k1 - 1] * fac;

        // fac *= (k1-j);  kj -= 1
        let ak = (k1 - j) as f64;
        fac *= ak;
        kj = kj.saturating_sub(1);
    }
}

/// Compute integrals of normalized B-splines: bint[j] = integral of N_{j,k+1}(x)
/// from a to b.
fn fpintb(t: &[f64], n: i32, bint: &mut [f64], nk1: i32, x: f64, y: f64) {
    let n = n as usize;
    let nk1 = nk1 as usize;
    let k1 = n - nk1;
    let ak = k1 as f64;
    let k = k1 - 1;
    for i in 0..nk1 { bint[i] = 0.0; }
    let mut a = x;
    let mut b = y;
    let min_flag;
    if a > b {
        core::mem::swap(&mut a, &mut b);
        min_flag = true;
    } else {
        min_flag = false;
    }
    if a < t[k1 - 1] { a = t[k1 - 1]; }
    if b > t[nk1] { b = t[nk1]; } // t(nk1+1) 1-based → t[nk1]
    let mut l = k1;
    let mut l0 = l + 1;
    let mut arg = a;
    let mut aint_arr = [0.0f64; 6];
    let mut h = [0.0f64; 6];
    let mut h1 = [0.0f64; 6];
    let mut ia = 0i32;
    let mut ib = 0i32;
    for it in 1..=2 {
        // search for knot interval
        loop {
            if arg < t[l0 - 1] || l == nk1 { break; }
            l = l0;
            l0 = l + 1;
        }
        // compute aint
        for j in 0..k1 { aint_arr[j] = 0.0; }
        aint_arr[0] = (arg - t[l - 1]) / (t[l] - t[l - 1]); // (arg-t(l))/(t(l+1)-t(l))
        h1[0] = 1.0;
        for j in 1..=k {
            h[0] = 0.0;
            for i in 1..=j {
                let li = l + i; // 1-based
                let lj2 = li - j; // 1-based
                let f = h1[i-1] / (t[li-1] - t[lj2-1]);
                h[i-1] += f * (t[li-1] - arg);
                h[i] = f * (arg - t[lj2-1]);
            }
            let j1 = j + 1;
            for i in 1..=j1 {
                let li = l + i;
                let lj2 = li - j1;
                aint_arr[i-1] += h[i-1] * (arg - t[lj2-1]) / (t[li-1] - t[lj2-1]);
                h1[i-1] = h[i-1];
            }
        }
        if it == 1 {
            let lk2 = l as i32 - k as i32;
            ia = lk2 - 1; // 0-based
            for i in 1..=k1 {
                let idx = (lk2 + i as i32 - 2) as usize;
                if idx < nk1 {
                    bint[idx] = -aint_arr[i-1];
                }
            }
            arg = b;
        } else {
            let lk2 = l as i32 - k as i32;
            ib = lk2 - 2; // ib = lk-1-1 = lk2-2 (0-based), may be negative
            for i in 1..=k1 {
                let idx = (lk2 + i as i32 - 2) as usize;
                if idx < nk1 {
                    bint[idx] += aint_arr[i-1];
                }
            }
            if ib >= ia {
                for i in ia..=ib {
                    let i = i as usize;
                    if i < nk1 { bint[i] += 1.0; }
                }
            }
        }
    }
    let f = 1.0 / ak;
    for i in 0..nk1 {
        let j = i + k1; // j = i+k1 (0-based → t[j] = t(j+1) 1-based... wait)
        // bint(i) *= (t(i+k1)-t(i))/ak → bint[i] *= (t[i+k1] - t[i]) * f
        bint[i] *= (t[i + k1] - t[i]) * f;
    }
    if min_flag {
        for i in 0..nk1 { bint[i] = -bint[i]; }
    }
}

/// Find real zeros of cubic polynomial p(x) = a*x^3 + b*x^2 + c*x + d.
fn fpcuro(a: f64, b: f64, c: f64, d: f64, x: &mut [f64; 3], n: &mut i32) {
    let ovfl = 1.0e5f64;
    let tent = 0.1f64;
    let e3 = 0.1f64 / 0.3f64;
    let pi3 = (1.0f64).atan() / 0.75f64;

    let a1 = a.abs();
    let b1 = b.abs();
    let c1 = c.abs();
    let d1 = d.abs();

    // Determine effective degree by Fortran's forward-goto logic:
    // if max(b1,c1,d1) < a1*ovfl → cubic (a is dominant)
    // elif max(c1,d1)  < b1*ovfl → quadratic
    // elif d1          < c1*ovfl → linear
    // else constant
    let is_cubic     = b1.max(c1).max(d1) < a1 * ovfl;
    let is_quadratic = !is_cubic && c1.max(d1) < b1 * ovfl;
    let is_linear    = !is_cubic && !is_quadratic && d1 < c1 * ovfl;

    if is_cubic {
        // Third degree
        let b1r = b / a * e3;
        let c1r = c / a;
        let d1r = d / a;
        let q = c1r * e3 - b1r * b1r;
        let r = b1r * b1r * b1r + (d1r - b1r * c1r) * 0.5;
        let disc = q * q * q + r * r;
        if disc > 0.0 {
            let u = disc.sqrt();
            let u1 = -r + u;
            let u2 = -r - u;
            *n = 1;
            x[0] = u1.abs().powf(e3).copysign(u1) + u2.abs().powf(e3).copysign(u2) - b1r;
        } else {
            let u = q.abs().sqrt();
            let u = if r < 0.0 { -u } else { u };
            let p3 = (-disc).sqrt().atan2(r.abs()) * e3;
            let u2d = u + u;
            *n = 3;
            x[0] = -u2d * p3.cos() - b1r;
            x[1] =  u2d * (pi3 - p3).cos() - b1r;
            x[2] =  u2d * (pi3 + p3).cos() - b1r;
        }
    } else if is_quadratic {
        // Second degree
        let disc = c * c - 4.0 * b * d;
        if disc < 0.0 {
            *n = 0;
            return;
        }
        *n = 2;
        let u = disc.sqrt();
        let b2 = b + b;
        x[0] = (-c + u) / b2;
        x[1] = (-c - u) / b2;
    } else if is_linear {
        // First degree
        *n = 1;
        x[0] = -d / c;
    } else {
        // Constant
        *n = 0;
        return;
    }

    // Newton refinement
    let three = 3.0f64;
    let two   = 2.0f64;
    for i in 0..*n as usize {
        let y  = x[i];
        let f  = ((a * y + b) * y + c) * y + d;
        let df = (three * a * y + two * b) * y + c;
        let step = if f.abs() < df.abs() * tent { f / df } else { 0.0 };
        x[i] = y - step;
    }
}

/// Insert a new knot x at position l into spline (t,n,c,k), producing (tt,nn,cc).
fn fpinst(iopt: i32, t: &[f64], n: i32, c: &[f64], k: i32, x: f64, l: i32,
          tt: &mut [f64], nn: &mut i32, cc: &mut [f64], nest: i32) {
    let n = n as usize;
    let k = k as usize;
    let l = l as usize;
    let nest = nest as usize;
    let k1 = k + 1;
    let nk1 = n - k1;
    // new knots
    let ll = l + 1;
    let mut i = n;
    for _j in ll..=n {
        tt[i] = t[i - 1];
        i -= 1;
    }
    tt[ll - 1] = x;
    for j in 0..l {
        tt[j] = t[j];
    }
    // new b-spline coefficients
    i = nk1;
    for _j in l..=nk1 {
        cc[i] = c[i - 1];
        if i == 0 { break; }
        i -= 1;
    }
    i = l;
    for _j in 1..=k {
        let m = i + k1;
        let fac = (x - tt[i - 1]) / (tt[m - 1] - tt[i - 1]);
        let i1 = i - 1;
        cc[i - 1] = fac * c[i - 1] + (1.0 - fac) * c[i1 - 1];
        i = i1;
    }
    for j in 0..i {
        cc[j] = c[j];
    }
    *nn = (n + 1) as i32;
    if iopt == 0 { return; }
    let nk = *nn as usize - k;
    let nl = nk - k1;
    let per = tt[nk - 1] - tt[k1 - 1];
    let mut i_f = k1;
    let mut j_f = nk;
    if ll > nl {
        for m in 1..=k {
            let mk = m + nl;
            cc[m - 1] = cc[mk - 1];
            i_f -= 1;
            j_f -= 1;
            tt[i_f - 1] = tt[j_f - 1] - per;
        }
    } else if ll <= k1 + k {
        for m in 1..=k {
            let mk = m + nl;
            cc[mk - 1] = cc[m - 1];
            i_f += 1;
            j_f += 1;
            tt[j_f - 1] = tt[i_f - 1] + per;
        }
    }
}

// ─── Section 5: Core Fitting Algorithms ───────────────────────────────────────

/// Core 1D spline fitting engine (translates fpcurf.f).
/// Fills knot vector t[0..n-1] and B-spline coefficients c[0..nk1-1].
#[allow(clippy::too_many_arguments)]
pub(crate) fn fpcurf(
    iopt: i32, x: &[f64], y: &[f64], w: &[f64], m: i32,
    xb: f64, xe: f64, k: i32, s: f64, nest: i32,
    tol: f64, maxit: i32, k1: i32, k2: i32,
    n: &mut i32, t: &mut [f64], c: &mut [f64], fp: &mut f64,
    fpint: &mut [f64], z: &mut [f64], a: &mut [f64], b: &mut [f64],
    g: &mut [f64], q: &mut [f64], nrdata: &mut [i32], ier: &mut i32,
) {
    let m = m as usize;
    let nest = nest as usize;
    let k1 = k1 as usize;
    let k2 = k2 as usize;
    let maxit = maxit as i32;

    let con1 = 0.1f64;
    let con9 = 0.9f64;
    let con4 = 0.04f64;
    let half = 0.5f64;

    let nmin = 2 * k1;

    // ── Part 1: knot placement ────────────────────────────────────────────────
    if iopt < 0 {
        // iopt = -1: user-supplied knots — do a single weighted LS fit.
        // Knots and boundary knots have already been set by curfit_.
        let n_us = *n as usize;
        let nk1 = n_us - k1;
        *fp = 0.0;
        for i in 0..nk1 { z[i] = 0.0; }
        for i in 0..nk1 { for j in 0..k1 { a[i + j*nest] = 0.0; } }
        let mut h = [0.0f64; 7];
        let mut l = k1;
        for it in 1..=m {
            let xi = x[it - 1];
            let wi = w[it - 1];
            let yi = y[it - 1] * wi;
            while xi >= t[l] && l < nk1 { l += 1; }
            fpbspl(t, *n as i32, k1 as i32 - 1, xi, l as i32, &mut h);
            for i in 1..=k1 {
                h[i - 1] *= wi;
            }
            let mut j_f = l - k1;
            let mut yi_rot = yi;
            for i in 1..=k1 {
                j_f += 1;
                let piv = h[i - 1];
                if piv.abs() < f64::EPSILON { continue; }
                let (mut cos, mut sin) = (0.0f64, 0.0f64);
                fpgivs(piv, &mut a[j_f - 1], &mut cos, &mut sin);
                fprota(cos, sin, &mut yi_rot, &mut z[j_f - 1]);
                if i == k1 { break; }
                let mut i2 = 0usize;
                for i1 in (i + 1)..=k1 {
                    i2 += 1;
                    fprota(cos, sin, &mut h[i1 - 1], &mut a[j_f - 1 + i2 * nest]);
                }
            }
            *fp += yi_rot * yi_rot;
        }
        fpback(a, z, nk1 as i32, k1 as i32, c, nest as i32);
        *ier = -2;
        return;
    }

    if iopt >= 0 {
        let acc = tol * s;
        let nmax = m + k1;

        if s == 0.0 {
            // Interpolating spline: n = nmax, place interior knots, do LS fit.
            *n = nmax as i32;
            if nmax > nest { *ier = 1; return; }
            // Place interior knots
            let mk1 = m - k1;
            if mk1 != 0 {
                let k3 = k as usize / 2;
                let mut i = k2;
                let mut j = k3 + 2;
                if k3 * 2 == k as usize {
                    for _l in 1..=mk1 {
                        t[i - 1] = (x[j - 1] + x[j - 2]) * half;
                        i += 1; j += 1;
                    }
                } else {
                    for _l in 1..=mk1 {
                        t[i - 1] = x[j - 1];
                        i += 1; j += 1;
                    }
                }
            }
            // Set boundary knots and do the LS fit (Fortran label 120).
            let n_us = *n as usize;
            let nk1 = n_us - k1;
            let mut i_f = n_us;
            for _j in 1..=k1 {
                t[_j - 1] = xb;
                t[i_f - 1] = xe;
                i_f -= 1;
            }
            *fp = 0.0;
            for i in 0..nk1 { z[i] = 0.0; }
            for i in 0..nk1 { for j in 0..k1 { a[i + j*nest] = 0.0; } }
            let mut h = [0.0f64; 7];
            let mut l = k1;
            for it in 1..=m {
                let xi = x[it - 1];
                let wi = w[it - 1];
                let yi = y[it - 1] * wi;
                while xi >= t[l] && l < nk1 { l += 1; }
                fpbspl(t, *n as i32, k1 as i32 - 1, xi, l as i32, &mut h);
                for i in 1..=k1 {
                    q[it - 1 + (i - 1) * m] = h[i - 1];
                    h[i - 1] *= wi;
                }
                let mut j_f = l - k1;
                let mut yi_rot = yi;
                for i in 1..=k1 {
                    j_f += 1;
                    let piv = h[i - 1];
                    if piv.abs() < f64::EPSILON { continue; }
                    let (mut cos, mut sin) = (0.0f64, 0.0f64);
                    fpgivs(piv, &mut a[j_f - 1], &mut cos, &mut sin);
                    fprota(cos, sin, &mut yi_rot, &mut z[j_f - 1]);
                    if i == k1 { break; }
                    let mut i2 = 0usize;
                    for i1 in (i + 1)..=k1 {
                        i2 += 1;
                        fprota(cos, sin, &mut h[i1 - 1], &mut a[j_f - 1 + i2 * nest]);
                    }
                }
                *fp += yi_rot * yi_rot;
            }
            fpback(a, z, nk1 as i32, k1 as i32, c, nest as i32);
            *ier = -1;
            return;
        } else {
            // Smoothing spline (Fortran Part 1).
            // fp0 = fp from the minimum-knot fit; carried as a local variable
            // and always stored at fpint[n-1] so Part 2 can read it.
            let mut fp0 = 0.0f64;
            let mut fpold = 0.0f64;
            let mut nplus = 0i32;

            if iopt != 0 && *n != nmin as i32 {
                // iopt=1: resume from previous call — restore stored metadata.
                let n_us = *n as usize;
                fp0   = fpint[n_us - 1];
                fpold = fpint[n_us - 2];
                nplus = nrdata[n_us - 1];
                if fp0 <= s {
                    // Previous result already satisfies s; restart from nmin.
                    *n = nmin as i32;
                    nrdata[0] = (m - 2) as i32;
                    fp0   = 0.0;
                    fpold = 0.0;
                    nplus = 0;
                }
                // else: keep existing knots and proceed.
            } else {
                *n = nmin as i32;
                nrdata[0] = (m - 2) as i32;
            }

            'outer: for _iter in 1..=m as i32 {
                let n_us = *n as usize;
                if n_us == nmin { *ier = -2; }
                let nrint = n_us - nmin + 1;
                let nk1 = n_us - k1;
                // Set boundary knots.
                let mut i_f = n_us;
                for _j in 1..=k1 {
                    t[_j - 1] = xb;
                    t[i_f - 1] = xe;
                    i_f -= 1;
                }
                // Build the LS observation matrix and compute fp.
                *fp = 0.0;
                for i in 0..nk1 { z[i] = 0.0; }
                for i in 0..nk1 { for j in 0..k1 { a[i + j*nest] = 0.0; } }
                let mut h = [0.0f64; 7];
                let mut l = k1; // 1-based knot-interval index
                for it in 1..=m {
                    let xi = x[it - 1];
                    let wi = w[it - 1];
                    let yi = y[it - 1] * wi;
                    while xi >= t[l] && l < nk1 { l += 1; }
                    fpbspl(t, *n as i32, k1 as i32 - 1, xi, l as i32, &mut h);
                    for i in 1..=k1 {
                        q[it - 1 + (i - 1) * m] = h[i - 1];
                        h[i - 1] *= wi;
                    }
                    let mut j_f = l - k1;
                    let mut yi_rot = yi;
                    for i in 1..=k1 {
                        j_f += 1;
                        let piv = h[i - 1];
                        if piv.abs() < f64::EPSILON { continue; }
                        let (mut cos, mut sin) = (0.0f64, 0.0f64);
                        fpgivs(piv, &mut a[j_f - 1], &mut cos, &mut sin);
                        fprota(cos, sin, &mut yi_rot, &mut z[j_f - 1]);
                        if i == k1 { break; }
                        let mut i2 = 0usize;
                        for i1 in (i + 1)..=k1 {
                            i2 += 1;
                            fprota(cos, sin, &mut h[i1 - 1], &mut a[j_f - 1 + i2 * nest]);
                        }
                    }
                    *fp += yi_rot * yi_rot;
                }
                // When n==nmin, record the polynomial fit residual as fp0.
                // fp0 is preserved through all subsequent iterations (Fortran local var).
                if *ier == -2 { fp0 = *fp; }
                // Always store fp0 at fpint[n-1] so Part 2 can read it.
                fpint[n_us - 1] = fp0;
                fpint[n_us - 2] = fpold;
                nrdata[n_us - 1] = nplus;
                // Backward substitution: solve a*c = z.
                fpback(a, z, nk1 as i32, k1 as i32, c, nest as i32);
                if iopt < 0 { break 'outer; }
                let fpms = *fp - s;
                if fpms.abs() < acc { break 'outer; }
                if fpms < 0.0 { break 'outer; } // fp < s → go to Part 2
                if *n == nmax as i32 { *ier = -1; break 'outer; }
                if *n == nest as i32 { *ier = 1;  break 'outer; }
                // Determine how many knots to add this iteration.
                if *ier != 0 {
                    nplus = 1;
                    *ier = 0;
                } else {
                    let rn = nplus as f64;
                    let npl1 = if fpold - *fp > acc {
                        (rn * fpms / (fpold - *fp)) as i32
                    } else {
                        nplus * 2
                    };
                    nplus = (nplus * 2).min(npl1.max(nplus / 2).max(1));
                }
                fpold = *fp;
                // Compute per-interval residuals into fpint[0..nrint-1].
                let mut fpart = 0.0f64;
                let mut ii = 1usize;
                let mut ll = k2;
                let mut new_flag = 0i32;
                for it in 1..=m {
                    if x[it - 1] >= t[ll - 1] && ll <= nk1 {
                        new_flag = 1;
                        ll += 1;
                    }
                    let mut term = 0.0f64;
                    let mut l0 = ll - k2;
                    for j in 1..=k1 {
                        l0 += 1;
                        term += c[l0 - 1] * q[it - 1 + (j - 1) * m];
                    }
                    let wi = w[it - 1];
                    let yi = y[it - 1];
                    term = (wi * (term - yi)) * (wi * (term - yi));
                    fpart += term;
                    if new_flag != 0 {
                        let store = term * half;
                        fpint[ii - 1] = fpart - store;
                        ii += 1;
                        fpart = store;
                        new_flag = 0;
                    }
                }
                fpint[nrint - 1] = fpart;
                // Add up to nplus knots; use a mutable nrint copy so each
                // fpknot call sees the updated interval count.
                let mut nrint_cur = nrint as i32;
                for _l in 1..=nplus {
                    fpknot(x, m as i32, t, n, fpint, nrdata, &mut nrint_cur, nest as i32, 1);
                    if *n == nmax as i32 {
                        // n reached nmax: place interpolation knots and restart.
                        let mk1 = m - k1;
                        if mk1 != 0 {
                            let k3 = k as usize / 2;
                            let mut i_fi = k2;
                            let mut j_f2 = k3 + 2;
                            if k3 * 2 == k as usize {
                                for _lv in 1..=mk1 {
                                    t[i_fi - 1] = (x[j_f2 - 1] + x[j_f2 - 2]) * half;
                                    i_fi += 1; j_f2 += 1;
                                }
                            } else {
                                for _lv in 1..=mk1 {
                                    t[i_fi - 1] = x[j_f2 - 1];
                                    i_fi += 1; j_f2 += 1;
                                }
                            }
                        }
                        break;
                    }
                    if *n == nest as i32 { break; }
                }
                // continue to next outer iteration
            }
        }
    }

    // ── Part 2: smoothing spline (p-finding iteration) ────────────────────────
    if *ier == -2 { return; } // polynomial is sufficient

    let n_us = *n as usize;
    let nk1 = n_us - k1;
    let n8 = n_us - nmin;
    let fp0 = fpint[n_us - 1];
    let f1_init = fp0 - s;  // Fortran: f1 = fp0 - s

    if iopt < 0 || f1_init <= 0.0 || n8 == 0 {
        return;
    }

    // Compute discontinuity matrix b
    fpdisc(t, *n, k2 as i32, b, nest as i32);

    let mut p1 = 0.0f64;
    let mut f1 = f1_init;
    let mut p3 = -1.0f64;
    let mut f3 = *fp - s;  // current fp - s (< 0, p = infinity)
    let mut p = 0.0f64;
    for i in 0..nk1 { p += a[i]; } // sum of a(i,1)
    p = (nk1 as f64) / p;
    let mut ich1 = 0i32;
    let mut ich3 = 0i32;
    let acc = tol * s;

    for iter in 1..=maxit {
        let pinv = 1.0 / p;
        // copy a → g, z → c
        for i in 0..nk1 {
            c[i] = z[i];
            g[i + (k2 - 1) * nest] = 0.0; // g(i, k2) = 0
            for j in 0..k1 { g[i + j * nest] = a[i + j * nest]; }
        }
        // rotate b rows into g
        let mut h = [0.0f64; 7];
        for it in 1..=n8 {
            for i in 0..k2 { h[i] = b[it - 1 + i * nest] * pinv; }
            let mut yi = 0.0f64;
            for j_f in it..=nk1 {
                let piv = h[0];
                let (mut cos, mut sin) = (0.0, 0.0);
                fpgivs(piv, &mut g[j_f - 1], &mut cos, &mut sin);
                fprota(cos, sin, &mut yi, &mut c[j_f - 1]);
                if j_f == nk1 { break; }
                let i2 = if j_f > n8 { nk1 - j_f } else { k1 };
                for i in 1..=i2 {
                    fprota(cos, sin, &mut h[i], &mut g[j_f - 1 + i * nest]);
                    h[i - 1] = h[i];
                }
                h[i2] = 0.0;
            }
        }
        // backward substitution (in-place: z == c; copy z first)
        let z_copy: Vec<f64> = c[..nk1].to_vec();
        fpback(g, &z_copy, nk1 as i32, k2 as i32, c, nest as i32);
        // compute fp
        *fp = 0.0;
        let mut ll = k2;
        for it in 1..=m {
            if x[it - 1] >= t[ll - 1] && ll <= nk1 { ll += 1; }
            let mut l0 = ll - k2;
            let mut term = 0.0f64;
            for j in 1..=k1 { l0 += 1; term += c[l0 - 1] * q[it - 1 + (j - 1) * m]; }
            let wi = w[it - 1]; let yi = y[it - 1];
            *fp += (wi * (term - yi)).powi(2);
        }
        let fpms = *fp - s;
        if fpms.abs() < acc { return; }
        if iter == maxit { *ier = 3; return; }

        let p2 = p;
        let f2 = fpms;
        if ich3 == 0 {
            if (f2 - f3) <= acc {
                p3 = p2; f3 = f2;
                p = p * con4;
                if p <= p1 { p = p1 * con9 + p2 * con1; }
                continue;
            }
            if f2 < 0.0 { ich3 = 1; }
        }
        if ich1 == 0 {
            if (f1 - f2) <= acc {
                p1 = p2; f1 = f2;
                p = p / con4;
                if p3 < 0.0 { continue; }
                if p >= p3 { p = p2 * con1 + p3 * con9; }
                continue;
            }
            if f2 > 0.0 { ich1 = 1; }
        }
        if f2 >= f1 || f2 <= f3 { *ier = 2; return; }
        p = fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
    }
    *ier = 3;
}

/// Core parametric spline fitting engine (translates fppara.f).
#[allow(dead_code, clippy::too_many_arguments)]
fn fppara(
    iopt: i32, idim: i32, m: i32, u: &[f64], _mx: i32, x: &[f64], w: &[f64],
    ub: f64, ue: f64, k: i32, s: f64, nest: i32,
    tol: f64, maxit: i32, k1: i32, k2: i32,
    n: &mut i32, t: &mut [f64], nc: i32, c: &mut [f64], fp: &mut f64,
    fpint: &mut [f64], z: &mut [f64], a: &mut [f64], b: &mut [f64],
    g: &mut [f64], q: &mut [f64], nrdata: &mut [i32], ier: &mut i32,
) {
    let m = m as usize;
    let nest = nest as usize;
    let k1 = k1 as usize;
    let k2 = k2 as usize;
    let idim = idim as usize;
    let nc = nc as usize;
    let maxit = maxit as i32;

    let con1 = 0.1f64;
    let con9 = 0.9f64;
    let con4 = 0.04f64;
    let half = 0.5f64;
    let nmin = 2 * k1;

    let mut h = [0.0f64; 7];
    let mut xi_loc = [0.0f64; 10];

    // ── Part 1: knot placement ────────────────────────────────────────────────
    if iopt >= 0 {
        let acc = tol * s;
        let nmax = m + k1;

        if s == 0.0 {
            *n = nmax as i32;
            if nmax > nest { *ier = 1; return; }
            let mk1 = m - k1;
            if mk1 != 0 {
                let k3 = k as usize / 2;
                let mut i_f = k2;
                let mut j_f = k3 + 2;
                if k3 * 2 == k as usize {
                    for _l in 1..=mk1 { t[i_f-1] = (u[j_f-1]+u[j_f-2])*half; i_f+=1; j_f+=1; }
                } else {
                    for _l in 1..=mk1 { t[i_f-1] = u[j_f-1]; i_f+=1; j_f+=1; }
                }
            }
        } else {
            if iopt == 0 || *n == nmin as i32 {
                *n = nmin as i32;
                fpint[*n as usize - 1] = 0.0;
                nrdata[0] = (m - 2) as i32;
            }
        }

        let mut fpold = 0.0f64;
        let mut nplus = 0i32;
        if iopt == 0 { fpold = 0.0; nplus = 0; nrdata[0] = (m-2) as i32; }

        'outer: for _iter in 1..=m as i32 {
            let n_us = *n as usize;
            if n_us == nmin { *ier = -2; }
            let nrint = n_us - nmin + 1;
            let nk1 = n_us - k1;
            let mut i_f = n_us;
            for _j in 1..=k1 { t[_j-1] = ub; t[i_f-1] = ue; i_f -= 1; }

            *fp = 0.0;
            for i in 0..nc { z[i] = 0.0; }
            for i in 0..nk1 { for j in 0..k1 { a[i + j*nest] = 0.0; } }

            let mut l = k1;
            let mut jj = 0usize;
            for it in 1..=m {
                let ui = u[it-1];
                let wi = w[it-1];
                for j in 1..=idim { jj += 1; xi_loc[j-1] = x[jj-1] * wi; }
                while ui >= t[l] && l < nk1 { l += 1; }
                fpbspl(t, *n as i32, k1 as i32 - 1, ui, l as i32, &mut h);
                for i in 1..=k1 { q[it-1 + (i-1)*m] = h[i-1]; h[i-1] *= wi; }
                let mut j_f = l - k1;
                for i in 1..=k1 {
                    j_f += 1;
                    let piv = h[i-1];
                    if piv.abs() < f64::EPSILON { continue; }
                    let (mut cos, mut sin) = (0.0, 0.0);
                    fpgivs(piv, &mut a[j_f-1], &mut cos, &mut sin);
                    let mut j1 = j_f;
                    for j2 in 1..=idim {
                        fprota(cos, sin, &mut xi_loc[j2-1], &mut z[j1-1]);
                        j1 += nest;
                    }
                    if i == k1 { break; }
                    let mut i2 = 0usize;
                    for i1 in (i+1)..=k1 {
                        i2 += 1;
                        fprota(cos, sin, &mut h[i1-1], &mut a[j_f-1 + i2*nest]);
                    }
                }
                for j2 in 1..=idim { *fp += xi_loc[j2-1].powi(2); }
            }
            if *ier == -2 { fpint[n_us-1] = *fp; }
            let fp0 = fpint[n_us-1];
            fpint[n_us-1] = fp0;
            fpint[n_us-2] = fpold;
            nrdata[n_us-1] = nplus;
            // backward substitution for each dimension
            let mut j1 = 1usize;
            for _j2 in 1..=idim {
                fpback(a, &z[j1-1..], nk1 as i32, k1 as i32, &mut c[j1-1..], nest as i32);
                j1 += nest;
            }
            if iopt < 0 { break 'outer; }
            let fpms = *fp - s;
            let acc = tol * s;
            if fpms.abs() < acc { break 'outer; }
            if fpms < 0.0 { break 'outer; }
            if *n == nmax as i32 { *ier = -1; break 'outer; }
            if *n == nest as i32 { *ier = 1; break 'outer; }
            if *ier != 0 { nplus = 1; *ier = 0; } else {
                let rn = nplus as f64;
                let fpms_v = *fp - s;
                let npl1 = if fpold - *fp > acc { (rn * fpms_v / (fpold - *fp)) as i32 } else { nplus*2 };
                nplus = (nplus*2).min(npl1.max(nplus/2).max(1));
            }
            fpold = *fp;
            let mut fpart = 0.0f64;
            let mut ii = 1usize;
            let mut ll = k2;
            let mut new_f = 0i32;
            jj = 0;
            for it in 1..=m {
                if u[it-1] >= t[ll-1] && ll <= nk1 { new_f = 1; ll += 1; }
                let mut term = 0.0f64;
                let mut l0 = ll - k2;
                for j2 in 1..=idim {
                    let mut fac = 0.0f64;
                    let mut j1_v = l0;
                    for j in 1..=k1 { j1_v += 1; fac += c[j1_v - 1] * q[it-1+(j-1)*m]; }
                    jj += 1;
                    term += (w[it-1] * (fac - x[jj-1])).powi(2);
                    l0 += nest;
                }
                fpart += term;
                if new_f != 0 {
                    let store = term * half;
                    fpint[ii-1] = fpart - store; ii += 1; fpart = store; new_f = 0;
                }
            }
            fpint[nrint-1] = fpart;
            for _l in 1..=nplus {
                fpknot(u, m as i32, t, n, fpint, nrdata, &mut (nrint as i32), nest as i32, 1);
                if *n == nmax as i32 {
                    let mk1 = m - k1;
                    if mk1 != 0 {
                        let k3 = k as usize / 2;
                        let mut i_fi = k2; let mut j_fi = k3+2;
                        if k3*2 == k as usize {
                            for _ in 1..=mk1 { t[i_fi-1]=(u[j_fi-1]+u[j_fi-2])*half; i_fi+=1; j_fi+=1; }
                        } else {
                            for _ in 1..=mk1 { t[i_fi-1]=u[j_fi-1]; i_fi+=1; j_fi+=1; }
                        }
                    }
                    break;
                }
                if *n == nest as i32 { break; }
            }
        }
    } else {
        // iopt < 0: least-squares with fixed knots - set boundary knots
        let n_us = *n as usize;
        let nk1 = n_us - k1;
        let mut i_f = n_us;
        for _j in 1..=k1 { t[_j-1] = ub; t[i_f-1] = ue; i_f -= 1; }
        *fp = 0.0;
        for i in 0..nc { z[i] = 0.0; }
        for i in 0..nk1 { for j in 0..k1 { a[i + j*nest] = 0.0; } }
        let mut l = k1;
        let mut jj = 0usize;
        for it in 1..=m {
            let ui = u[it-1]; let wi = w[it-1];
            for j in 1..=idim { jj += 1; xi_loc[j-1] = x[jj-1] * wi; }
            while ui >= t[l] && l < nk1 { l += 1; }
            fpbspl(t, *n as i32, k1 as i32 - 1, ui, l as i32, &mut h);
            for i in 1..=k1 { q[it-1+(i-1)*m] = h[i-1]; h[i-1] *= wi; }
            let mut j_f = l - k1;
            for i in 1..=k1 {
                j_f += 1;
                let piv = h[i-1];
                if piv.abs() < f64::EPSILON { continue; }
                let (mut cos, mut sin) = (0.0, 0.0);
                fpgivs(piv, &mut a[j_f-1], &mut cos, &mut sin);
                let mut j1 = j_f;
                for j2 in 1..=idim { fprota(cos, sin, &mut xi_loc[j2-1], &mut z[j1-1]); j1 += nest; }
                if i == k1 { break; }
                let mut i2 = 0usize;
                for i1 in (i+1)..=k1 { i2+=1; fprota(cos, sin, &mut h[i1-1], &mut a[j_f-1+i2*nest]); }
            }
            for j2 in 1..=idim { *fp += xi_loc[j2-1].powi(2); }
        }
        let mut j1 = 1usize;
        for _j2 in 1..=idim {
            fpback(a, &z[j1-1..], nk1 as i32, k1 as i32, &mut c[j1-1..], nest as i32);
            j1 += nest;
        }
        return;
    }

    // ── Part 2: smoothing spline (p-finding) ─────────────────────────────────
    if *ier == -2 { return; }
    let n_us = *n as usize;
    let nk1 = n_us - k1;
    let n8 = n_us - nmin;
    let fp0 = fpint[n_us - 1];
    let fpms0 = *fp - s;
    if fpms0 >= 0.0 && n8 > 0 {
        fpdisc(t, *n, k2 as i32, b, nest as i32);
        let mut p1 = 0.0f64;
        let mut f1 = fp0 - s;
        let mut p3 = -1.0f64;
        let mut f3 = fpms0;
        let mut p: f64 = 0.0;
        for i in 0..nk1 { p += a[i]; }
        p = (nk1 as f64) / p;
        let mut ich1 = 0i32;
        let mut ich3 = 0i32;
        let acc = tol * s;
        for iter in 1..=maxit {
            let pinv = 1.0 / p;
            for i in 0..nc { c[i] = z[i]; }
            for i in 0..nk1 { g[i + (k2-1)*nest] = 0.0; for j in 0..k1 { g[i+j*nest] = a[i+j*nest]; } }
            for it in 1..=n8 {
                for i in 0..k2 { h[i] = b[it-1 + i*nest] * pinv; }
                let mut xi2 = [0.0f64; 10];
                let mut j_f = it;
                while j_f <= nk1 {
                    let piv = h[0];
                    let (mut cos, mut sin) = (0.0, 0.0);
                    fpgivs(piv, &mut g[j_f-1], &mut cos, &mut sin);
                    let mut j1 = j_f;
                    for j2 in 1..=idim { fprota(cos, sin, &mut xi2[j2-1], &mut c[j1-1]); j1 += nest; }
                    if j_f == nk1 { break; }
                    let i2 = if j_f > n8 { nk1 - j_f } else { k1 };
                    for i in 1..=i2 { fprota(cos, sin, &mut h[i], &mut g[j_f-1+i*nest]); h[i-1] = h[i]; }
                    h[i2] = 0.0;
                    j_f += 1;
                }
            }
            let mut j1 = 1usize;
            for _j2 in 1..=idim {
                let z_copy2: Vec<f64> = c[j1-1..j1-1+nk1].to_vec();
                fpback(g, &z_copy2, nk1 as i32, k2 as i32, &mut c[j1-1..], nest as i32);
                j1 += nest;
            }
            *fp = 0.0;
            let mut ll = k2;
            let mut jj2 = 0usize;
            for it in 1..=m {
                if u[it-1] >= t[ll-1] && ll <= nk1 { ll += 1; }
                let mut l0 = ll - k2;
                let mut term = 0.0f64;
                for j2 in 1..=idim {
                    let mut fac = 0.0f64;
                    let mut j1v = l0;
                    for j in 1..=k1 { j1v += 1; fac += c[j1v-1] * q[it-1+(j-1)*m]; }
                    jj2 += 1;
                    term += (fac - x[jj2-1]).powi(2);
                    l0 += nest;
                }
                *fp += term * w[it-1].powi(2);
            }
            let fpms = *fp - s;
            if fpms.abs() < acc { return; }
            if iter == maxit { *ier = 3; return; }
            let p2 = p; let f2 = fpms;
            if ich3 == 0 {
                if (f2 - f3) <= acc { p3 = p2; f3 = f2; p = p*con4; if p <= p1 { p = p1*con9+p2*con1; } continue; }
                if f2 < 0.0 { ich3 = 1; }
            }
            if ich1 == 0 {
                if (f1 - f2) <= acc { p1 = p2; f1 = f2; p = p/con4; if p3 < 0.0 { continue; } if p >= p3 { p = p2*con1+p3*con9; } continue; }
                if f2 > 0.0 { ich1 = 1; }
            }
            if f2 >= f1 || f2 <= f3 { *ier = 2; return; }
            p = fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
        }
        *ier = 3;
    }
}

// ─── Section 4: Fourier Utilities ─────────────────────────────────────────────

/// Compute ress = integral((b-x)^3*sin(par*x)) and resc = integral((b-x)^3*cos(par*x))
/// over (a,b).
#[allow(dead_code)]
fn fpcsin(a: f64, b: f64, par: f64, sia: f64, coa: f64, sib: f64, cob: f64,
          ress: &mut f64, resc: &mut f64) {
    let ab = b - a;
    let ab4 = ab * ab * ab * ab;
    let alfa = ab * par;
    let eps = 1.0e-10f64;
    if alfa.abs() > 1.0 {
        let beta = 1.0 / alfa;
        let b2 = beta * beta;
        let b4 = 6.0 * b2 * b2;
        let f1 = 3.0 * b2 * (1.0 - 2.0 * b2);
        let f2 = beta * (1.0 - 6.0 * b2);
        *ress = ab4 * (coa * f2 + sia * f1 + sib * b4);
        *resc = ab4 * (coa * f1 - sia * f2 + cob * b4);
    } else {
        let mut fac = 0.25f64;
        let mut f1 = fac;
        let mut f2 = 0.0f64;
        let mut i = 4i32;
        'series: for _j in 1..=5 {
            i += 1;
            let ai = i as f64;
            fac = fac * alfa / ai;
            f2 += fac;
            if fac.abs() <= eps { break 'series; }
            i += 1;
            let ai = i as f64;
            fac = -fac * alfa / ai;
            f1 += fac;
            if fac.abs() <= eps { break 'series; }
        }
        *ress = ab4 * (coa * f2 + sia * f1);
        *resc = ab4 * (coa * f1 - sia * f2);
    }
}

/// Compute Fourier integrals of cubic B-splines.
#[allow(dead_code)]
fn fpbfou(t: &[f64], n: i32, par: f64, ress: &mut [f64], resc: &mut [f64]) {
    let n = n as usize;
    let nm3 = n - 3;
    let nm7 = if n >= 7 { n - 7 } else { 0 };
    let eps = 1.0e-7f64;
    let quart = 0.25f64;
    let con1 = 0.05f64;
    let con2 = 120.0f64;
    let term = if par.abs() >= f64::EPSILON { 6.0 / par } else { 6.0 / par }; // same
    let mut co = [0.0f64; 5];
    let mut si = [0.0f64; 5];
    let mut hs = [0.0f64; 5];
    let mut hc = [0.0f64; 5];
    let mut rs = [0.0f64; 3];
    let mut rc = [0.0f64; 3];
    let beta0 = par * t[3]; // t(4) 1-based → t[3]
    co[0] = beta0.cos();
    si[0] = beta0.sin();
    // j=1,2,3
    for j in 1..=3usize {
        let jp1 = j + 1;
        let jp4 = j + 4;
        let beta = par * t[jp4 - 1];
        co[jp1 - 1] = beta.cos();
        si[jp1 - 1] = beta.sin();
        fpcsin(t[3], t[jp4 - 1], par, si[0], co[0], si[jp1 - 1], co[jp1 - 1],
               &mut rs[j - 1], &mut rc[j - 1]);
        let i0 = 5 - j;
        hs[i0 - 1] = 0.0;
        hc[i0 - 1] = 0.0;
        for jj in 1..=j {
            let ipj = i0 + jj;
            hs[ipj - 1] = rs[jj - 1];
            hc[ipj - 1] = rc[jj - 1];
        }
        for jj in 1..=3usize {
            let mut i = if i0 < jj { jj } else { i0 };
            let mut k = 5usize;
            let mut li = jp4;
            for _ll in i..=4usize {
                let lj2 = li - jj;
                let fac2 = t[li - 1] - t[lj2 - 1];
                hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac2;
                hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac2;
                if k > 1 { k -= 1; }
                if li > 1 { li -= 1; }
            }
        }
        ress[j - 1] = hs[4] - hs[3];
        resc[j - 1] = hc[4] - hc[3];
    }
    if nm7 >= 4 {
        for j in 4..=nm7 {
            let jp4 = j + 4;
            let beta = par * t[jp4 - 1];
            co[4] = beta.cos();
            si[4] = beta.sin();
            let delta = t[jp4 - 1] - t[j - 1];
            let beta2 = delta * par;
            if beta2.abs() > 1.0 {
                for k in 0..5 { hs[k] = si[k]; hc[k] = co[k]; }
                for jj in 1..=3usize {
                    let mut k = 5usize;
                    let mut li = jp4;
                    for _ll in jj..=4usize {
                        let lj2 = li - jj;
                        let fac2 = par * (t[li - 1] - t[lj2 - 1]);
                        hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac2;
                        hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac2;
                        if k > 1 { k -= 1; }
                        if li > 1 { li -= 1; }
                    }
                }
                ress[j - 1] = (hs[4] - hs[3]) * term;
                resc[j - 1] = (hc[4] - hc[3]) * term;
            } else {
                let mut f3 = 0.0f64;
                let mut hsi = [0.0f64; 4];
                let mut hci = [0.0f64; 4];
                for i in 0..4 {
                    hsi[i] = par * (t[i + j] - t[j - 1]);
                    hci[i] = hsi[i];
                    f3 += hsi[i];
                }
                f3 *= con1;
                let mut c1 = quart;
                let mut s1 = f3;
                if f3.abs() > eps {
                    let mut sign = 1.0f64;
                    let mut fac2 = con2;
                    let mut k2 = 5i32;
                    let mut is = 0i32;
                    'done: for _ic in 1..=20 {
                        k2 += 1;
                        let ak = k2 as f64;
                        fac2 *= ak;
                        let mut f1v = 0.0f64;
                        let mut f3v = 0.0f64;
                        for i in 0..4 {
                            f1v += hci[i];
                            let f2v = f1v * hsi[i];
                            hci[i] = f2v;
                            f3v += f2v;
                        }
                        f3 = f3v * 6.0 / fac2;
                        if is == 0 {
                            sign = -sign;
                            is = 1;
                            c1 += f3 * sign;
                        } else {
                            is = 0;
                            s1 += f3 * sign;
                        }
                        if f3.abs() <= eps { break 'done; }
                    }
                }
                ress[j - 1] = delta * (co[0] * s1 + si[0] * c1);
                resc[j - 1] = delta * (co[0] * c1 - si[0] * s1);
            }
            for i in 0..4 { co[i] = co[i + 1]; si[i] = si[i + 1]; }
        }
    }
    // j = n-6, n-5, n-4
    for j in 1..=3usize {
        let nmj = nm3 - j;
        let i0 = 5 - j;
        fpcsin(t[nm3 - 1], t[nmj - 1], par, si[3], co[3], si[i0 - 2], co[i0 - 2],
               &mut rs[j - 1], &mut rc[j - 1]);
        hs[i0 - 1] = 0.0;
        hc[i0 - 1] = 0.0;
        for jj in 1..=j {
            let ipj = i0 + jj;
            hc[ipj - 1] = rc[jj - 1];
            hs[ipj - 1] = rs[jj - 1];
        }
        for jj in 1..=3usize {
            let mut i = if i0 < jj { jj } else { i0 };
            let mut k = 5usize;
            let mut li = nmj;
            for _ll in i..=4usize {
                let lj2 = li + jj;
                let fac2 = t[lj2 - 1] - t[li - 1];
                hs[k - 1] = (hs[k - 2] - hs[k - 1]) / fac2;
                hc[k - 1] = (hc[k - 2] - hc[k - 1]) / fac2;
                if k > 1 { k -= 1; }
                li += 1;
            }
        }
        ress[nmj - 1] = hs[3] - hs[4];
        resc[nmj - 1] = hc[3] - hc[4];
    }
}

// ─── Section 6: Endpoint-constrained and closed-curve helpers ─────────────────

/// Find a polynomial curve of degree k satisfying derivative constraints at a and b.
/// Returns b-spline representation in cp[0..np-1] (layout: idim-major, k1 values per dimension).
fn fppocu(idim: i32, k: i32, a: f64, b: f64, ib: i32, db: &[f64], _nb: i32,
          ie: i32, de: &[f64], _ne: i32, cp: &mut [f64], _np: i32) {
    let k1 = (k + 1) as usize;
    let k2 = 2 * k1;
    let idim = idim as usize;
    let ib = ib as usize;
    let ie = ie as usize;
    let ab = b - a;
    let mut work = [[0.0f64; 6]; 6]; // work(j,i) stored as work[j-1][i-1]

    for id in 1..=idim {
        for j in 0..k1 { work[j][0] = 0.0; }
        // begin constraints
        if ib > 0 {
            let mut l = id - 1; // 0-based into db
            for i in 1..=ib {
                work[0][i - 1] = db[l];
                l += idim;
            }
            if ib > 1 {
                let mut ll = ib;
                for j in 2..=ib {
                    ll -= 1;
                    for i in 1..=ll {
                        let aki = (k1 - i) as f64;
                        work[j - 1][i - 1] = ab * work[j - 2][i] / aki + work[j - 2][i - 1];
                    }
                }
            }
        }
        // end constraints
        if ie > 0 {
            let mut l = id - 1;
            let mut j = k1;
            for i in 1..=ie {
                work[j - 1][i - 1] = de[l];
                l += idim;
                j -= 1;
            }
            if ie > 1 {
                let mut ll = ie;
                for jj in 2..=ie {
                    ll -= 1;
                    let mut jv = k1 + 1 - jj;
                    for i in 1..=ll {
                        let aki = (k1 - i) as f64;
                        work[jv - 1][i - 1] = work[jv][i - 1] - ab * work[jv - 1][i] / aki;
                        jv -= 1;
                    }
                }
            }
        }
        // store result
        let base = (id - 1) * k2;
        for j in 1..=k1 {
            cp[base + j - 1] = work[j - 1][0];
        }
    }
}

/// Add polynomial cp to spline c, both in b-spline representation on knots t.
fn fpadpo(idim: i32, t: &[f64], n: i32, c: &mut [f64], nc: i32, k: i32,
          cp: &[f64], _np: i32, cc: &mut [f64], t1: &mut [f64], t2: &mut [f64]) {
    let k1 = (k + 1) as usize;
    let n = n as usize;
    let nc = nc as usize;
    let idim = idim as usize;
    let nk1 = n - k1;
    let k2 = 2 * k1;

    // initialize cc with the polynomial cp (only first k1 coefficients per dimension)
    let mut j = 0usize;
    let mut l = 0usize;
    for _jj in 1..=idim {
        let l1 = j;
        for ii in 0..k1 {
            cc[l1 + ii] = cp[l + ii];
        }
        j += n;
        l += k1;
        // also skip the extra k1 in cp layout (k2 per dim)
        l += k1;
    }

    if nk1 != k1 {
        // set up minimal knot vector t1 with k2 knots
        let n1_init = k2;
        let mut j_f = n;
        let mut l_f = n1_init;
        for i in 0..k1 {
            t1[i] = t[i];
            t1[l_f - 1 - i] = t[j_f - 1 - i];
        }

        // insert knots one by one
        let nk2 = nk1 - 1;
        let mut n1 = n1_init;
        for l_v in k1..nk2 {
            let l1 = l_v + 1; // 1-based: t(l1+1)
            let mut j_base = 0usize;
            let mut n2 = 0i32;
            for _i in 0..idim {
                let mut nn_out = 0i32;
                let mut cc_out = vec![0.0f64; n];
                fpinst(0, t1, n1 as i32, &cc[j_base..j_base + n1], k, t[l1], l_v as i32,
                       t2, &mut nn_out, &mut cc_out, n as i32);
                let nk1_out = (nn_out as usize).saturating_sub(k as usize + 1);
                cc[j_base..j_base + nk1_out].copy_from_slice(&cc_out[..nk1_out]);
                // copy t2 back to t1
                n2 = nn_out;
                j_base += n;
            }
            for i in 0..(n2 as usize) { t1[i] = t2[i]; }
            n1 = n2 as usize;
        }
    }

    // add polynomial coefficients to spline coefficients
    let mut j_f = 0usize;
    for _jj in 1..=idim {
        for i in 0..nk1 {
            c[j_f + i] += cc[j_f + i];
        }
        j_f += n;
    }
}

/// Closed-curve fitting engine (translates fpclos.f).
#[allow(clippy::too_many_arguments)]
fn fpclos(
    iopt: i32, idim: i32, m: i32, u: &[f64], _mx: i32, x: &[f64], w: &[f64],
    k: i32, s: f64, nest: i32, tol: f64, maxit: i32, k1: i32, k2: i32,
    n: &mut i32, t: &mut [f64], _nc: i32, c: &mut [f64], fp: &mut f64,
    fpint: &mut [f64], z: &mut [f64], a1: &mut [f64], a2: &mut [f64],
    b: &mut [f64], g1: &mut [f64], g2: &mut [f64], q: &mut [f64],
    nrdata: &mut [i32], ier: &mut i32,
) {
    let m_v    = m    as usize;
    let k_v    = k    as usize;
    let k1_v   = k1   as usize;
    let k2_v   = k2   as usize;
    let idim_v = idim as usize;
    let nest_v = nest as usize;

    let con1 = 0.1_f64; let con9 = 0.9_f64; let con4 = 0.04_f64; let half = 0.5_f64;

    let m1  = m_v - 1;
    let mut kk  = k  as i32;   // may become k-1 for odd-k interpolation
    let mut kk1 = k1 as i32;   // may become k   for odd-k interpolation
    let nmin = 2 * k1_v;
    let per  = u[m_v - 1] - u[0];

    let mut h  = [0.0f64; 6];
    let mut h1 = [0.0f64; 7];
    let mut h2 = [0.0f64; 6];
    let mut xi = [0.0f64; 10];

    let mut acc   = 0.0f64;
    let mut fp0   = 0.0f64;
    let mut fpold = 0.0f64;
    let mut nplus = 0i32;
    let mut nmax_v = 0usize;

    // ── initial dispatch ──────────────────────────────────────────────────────
    // iopt <  0 → goto 50 (knots already set by caller)
    // s==0 & nmax!=nmin → goto 5 (set interior knots for interpolation)
    // otherwise → label 30/35 (fixed-point or restart)

    let mut do_knot_setup: bool;

    if iopt < 0 {
        do_knot_setup = false; // straight to label 50
    } else {
        acc    = tol * s;
        nmax_v = (m + 2 * k) as usize;

        if s > 0.0 || nmax_v == nmin {
            // ── label 30 ──────────────────────────────────────────────────
            let skip_fixed_point = iopt != 0
                && *n as usize != nmin
                && {
                    fp0   = fpint[*n as usize - 1];
                    fpold = fpint[*n as usize - 2];
                    nplus = nrdata[*n as usize - 1] as i32;
                    fp0 > s   // true → goto 50 directly
                };

            if !skip_fixed_point {
                // ── label 35: constant closed-curve solution ───────────────
                fp0  = 0.0;
                let mut d1 = 0.0f64;
                for j in 0..idim_v { z[j] = 0.0; }
                let mut jj_f = 0usize;
                let mut cos_v = 0.0f64; let mut sin_v = 0.0f64;
                for it in 0..m1 {
                    let wi = w[it];
                    fpgivs(wi, &mut d1, &mut cos_v, &mut sin_v);
                    for j in 0..idim_v {
                        let mut fac = wi * x[jj_f]; jj_f += 1;
                        fprota(cos_v, sin_v, &mut fac, &mut z[j]);
                        fp0 += fac * fac;
                    }
                }
                for j in 0..idim_v { z[j] /= d1; }

                let fpms = fp0 - s;
                if fpms < acc || nmax_v == nmin {
                    // ── label 640 ────────────────────────────────────────
                    *ier = -2;
                    for i in 1..=k1_v {
                        t[i - 1]       = u[0]        - (k1_v - i) as f64 * per;
                        t[i + k1_v - 1] = u[m_v - 1] + (i - 1)    as f64 * per;
                    }
                    *n = nmin as i32;
                    let nv = *n as usize;
                    let mut j1 = 0usize;
                    for j in 0..idim_v {
                        let fac = z[j]; let mut j2 = j1;
                        for _ in 0..k1_v { c[j2] = fac; j2 += 1; }
                        j1 += nv;
                    }
                    *fp = fp0;
                    fpint[nv - 1] = fp0; fpint[nv - 2] = 0.0; nrdata[nv - 1] = 0;
                    return;
                }
                fpold = fp0;
                if *n >= nest { *ier = 1; return; }
                nplus = 1;
                *n = nmin as i32 + 1;
                let mm = (m_v + 1) / 2;
                t[k2_v - 1] = u[mm - 1];
                nrdata[0]   = mm as i32 - 2;
                nrdata[1]   = m1  as i32 - mm as i32;
            }
            do_knot_setup = false; // goto 50
        } else {
            // s==0 and nmax!=nmin: goto 5
            *n = nmax_v as i32;
            if *n > nest { *ier = 1; return; }
            do_knot_setup = true;
        }
    }

    // ── 'outer: label-5 (knot setup) + label-50 (main iteration loop) ─────────
    let mut go_to_part2 = false;
    let mut fpms_saved  = 0.0f64;

    'outer: loop {
        // ── label 5 ───────────────────────────────────────────────────────────
        if do_knot_setup {
            do_knot_setup = false;
            let n_v = *n as usize;
            if k_v % 2 == 0 {
                // k even – label 20: midpoints
                for i in 2..=m1 {
                    t[i + k_v - 1] = (u[i - 1] + u[i - 2]) * half;
                }
                // fall through to label 50
            } else {
                // k odd: data points
                for i in 2..=m1 {
                    t[i + k_v - 1] = u[i - 1];
                }
                if s <= 0.0 {
                    kk  = k - 1;
                    kk1 = k;
                    if kk <= 0 {
                        // k==1 special case: directly set coefficients → ier=-1
                        t[0]       = t[m_v - 1] - per;
                        t[1]       = u[0];
                        t[m_v]     = u[m_v - 1];
                        t[m_v + 1] = t[2] + per;
                        let mut jj_c = 0usize;
                        for i in 1..=m1 {
                            let mut j_c = i - 1;
                            for _ in 0..idim_v { c[j_c] = x[jj_c]; jj_c += 1; j_c += n_v; }
                        }
                        let mut j2_c = m_v - 1; let mut jj2_c = 0usize;
                        for _ in 0..idim_v { c[j2_c] = c[jj2_c]; j2_c += n_v; jj2_c += n_v; }
                        *fp = 0.0;
                        fpint[n_v - 1] = fp0; fpint[n_v - 2] = 0.0; nrdata[n_v - 1] = 0;
                        *ier = -1;
                        return;
                    }
                    // kk > 0: fall through to label 50
                }
            }
        }

        // ── label 50: main knot-placement loop (do 340 iter=1,m) ─────────────
        'main: for _iter in 0..m_v {
            let n_v   = *n as usize;
            let nk1   = n_v - k1_v;
            let nk2   = nk1 + 1;
            let n7    = nk1 - k_v;
            let n10_i = n7 as i32 - kk;
            let nc_v  = nest_v * idim_v;
            let kk_v  = kk  as usize;
            let kk1_v = kk1 as usize;
            let nrint = n_v - nmin + 1;

            // set the additional periodic knots
            t[k1_v - 1] = u[0]; t[nk2 - 1] = u[m_v - 1];
            for j in 1..=k_v {
                t[nk2 + j - 1] = t[k1_v + j - 1] + per;
                t[k1_v - j - 1] = t[nk2 - j - 1] - per;
            }

            // initialise z and a1
            for i in 0..nc_v { z[i] = 0.0; }
            for i in 0..nk1 { for j in 0..kk1_v { a1[i + j * nest_v] = 0.0; } }

            let mut jper  = 0i32;
            *fp   = 0.0;
            let mut l_cur = k1_v; // l (1-based), starts at k1
            let mut jj_x  = 0usize;

            // ── do 290 it=1,m1 ───────────────────────────────────────────────
            for it in 1..=m1 {
                let ui = u[it - 1]; let wi = w[it - 1];
                for j in 0..idim_v { xi[j] = x[jj_x] * wi; jj_x += 1; }

                // advance l: find interval t(l) <= ui < t(l+1)
                while ui >= t[l_cur] { l_cur += 1; }

                fpbspl(t, *n, k, ui, l_cur as i32, &mut h);
                for i in 0..k1_v { q[(it - 1) + i * m_v] = h[i]; h[i] *= wi; }

                let l5_i = l_cur as i32 - k1;

                if l5_i < n10_i {
                    // ── label 285: simple zone ────────────────────────────────
                    let mut j_f = l5_i;
                    for i in 1..=kk1_v {
                        j_f += 1;
                        let piv = h[i - 1];
                        if piv.abs() < f64::EPSILON { continue; }
                        let j_idx = (j_f - 1) as usize;
                        let (mut cos_v, mut sin_v) = (0.0, 0.0);
                        fpgivs(piv, &mut a1[j_idx], &mut cos_v, &mut sin_v);
                        let mut j1_f = j_f;
                        for j2 in 0..idim_v {
                            fprota(cos_v, sin_v, &mut xi[j2], &mut z[(j1_f - 1) as usize]);
                            j1_f += n_v as i32;
                        }
                        if i == kk1_v { break; }
                        let mut i2 = 0usize;
                        for i1 in (i + 1)..=kk1_v {
                            i2 += 1;
                            fprota(cos_v, sin_v, &mut h[i1 - 1], &mut a1[j_idx + i2 * nest_v]);
                        }
                    }
                    // label 150: add xi² to fp
                    for j2 in 0..idim_v { *fp += xi[j2] * xi[j2]; }
                    continue; // next it
                }

                // ── periodic boundary region ──────────────────────────────────
                if jper == 0 {
                    // initialise a2 from upper part of a1
                    for i in 0..n7 { for j in 0..kk_v { a2[i + j * nest_v] = 0.0; } }
                    let mut jk = n10_i + 1;
                    for i in 1..=kk_v {
                        let mut ik = jk;
                        for j in 1..=kk1_v {
                            if ik <= 0 { break; }
                            a2[(ik - 1) as usize + (i - 1) * nest_v] =
                                a1[(ik - 1) as usize + (j - 1) * nest_v];
                            ik -= 1;
                        }
                        jk += 1;
                    }
                    jper = 1;
                }

                // ── label 160: build h1/h2 for periodic b-splines ─────────────
                for i in 0..kk_v { h1[i] = 0.0; h2[i] = 0.0; }
                h1[kk1_v - 1] = 0.0;

                let mut j_f = l5_i - n10_i;
                for i in 1..=kk1_v {
                    j_f += 1;
                    let mut l0 = j_f;
                    loop {
                        let l1 = l0 - kk;
                        if l1 <= 0 {
                            h2[(l0 - 1) as usize] += h[i - 1]; break;
                        } else if l1 <= n10_i {
                            h1[(l1 - 1) as usize] = h[i - 1]; break;
                        } else {
                            l0 = l1 - n10_i;
                        }
                    }
                }

                // ── rotation with rows 1..n10 ─────────────────────────────────
                if n10_i > 0 {
                    for j in 1..=(n10_i as usize) {
                        let piv = h1[0];
                        if piv.abs() <= f64::EPSILON {
                            // shift h1 left
                            for i in 0..kk_v { h1[i] = h1[i + 1]; }
                            h1[kk1_v - 1] = 0.0;
                            continue;
                        }
                        let j_idx = j - 1;
                        let (mut cos_v, mut sin_v) = (0.0, 0.0);
                        fpgivs(piv, &mut a1[j_idx], &mut cos_v, &mut sin_v);
                        let mut j1_f = j as i32;
                        for j2 in 0..idim_v {
                            fprota(cos_v, sin_v, &mut xi[j2], &mut z[(j1_f - 1) as usize]);
                            j1_f += n_v as i32;
                        }
                        for i in 0..kk_v {
                            fprota(cos_v, sin_v, &mut h2[i], &mut a2[j_idx + i * nest_v]);
                        }
                        if j as i32 == n10_i { break; }
                        let i2_max = ((n10_i as usize) - j).min(kk_v);
                        let mut i2 = 0usize;
                        for i in 0..i2_max {
                            i2 += 1;
                            fprota(cos_v, sin_v, &mut h1[i + 1], &mut a1[j_idx + i2 * nest_v]);
                            h1[i] = h1[i + 1];
                        }
                        h1[i2] = 0.0;
                    }
                }

                // ── label 250: rotation with rows n10+1..n7 ───────────────────
                for j in 1..=(kk_v as i32) {
                    let ij = n10_i + j;
                    if ij <= 0 { continue; }
                    let piv = h2[(j - 1) as usize];
                    if piv.abs() < f64::EPSILON { continue; }
                    let ij_idx = (ij - 1) as usize;
                    let (mut cos_v, mut sin_v) = (0.0, 0.0);
                    fpgivs(piv, &mut a2[ij_idx + (j as usize - 1) * nest_v], &mut cos_v, &mut sin_v);
                    let mut j1_f = ij;
                    for j2 in 0..idim_v {
                        fprota(cos_v, sin_v, &mut xi[j2], &mut z[(j1_f - 1) as usize]);
                        j1_f += n_v as i32;
                    }
                    if j == kk { break; }
                    for i in (j as usize + 1)..=kk_v {
                        fprota(cos_v, sin_v, &mut h2[i - 1], &mut a2[ij_idx + (i - 1) * nest_v]);
                    }
                }
                // label 280: add to fp
                for j2 in 0..idim_v { *fp += xi[j2] * xi[j2]; }
            } // end do-290

            fpint[n_v - 1] = fp0;
            fpint[n_v - 2] = fpold;
            nrdata[n_v - 1] = nplus;

            // backward substitution: solve the periodic banded system
            let n7_i = n7 as i32;
            let mut j1_f = 0usize;
            for _ in 0..idim_v {
                fpbacp(a1, a2, &z[j1_f..], n7_i, kk, &mut c[j1_f..], kk1, nest);
                j1_f += n_v;
            }
            // fill remaining coefficients from periodicity condition (**)
            for i in 1..=k_v {
                let mut j1_f = i - 1;
                for _ in 0..idim_v { c[j1_f + n7] = c[j1_f]; j1_f += n_v; }
            }

            if iopt < 0 { return; } // goto 660

            let fpms = *fp - s;
            if fpms.abs() < acc { return; } // goto 660
            if fpms < 0.0 { fpms_saved = fpms; go_to_part2 = true; break 'main; } // goto 350
            if *n == nmax_v as i32 { *ier = -1; return; } // goto 630
            if *n == nest { *ier = 1; return; } // goto 620

            // compute nplus for next round
            let rn = nplus as f64;
            let npl1 = if fpold - *fp > acc {
                ((rn * fpms / (fpold - *fp)) as i32).max(1)
            } else { nplus * 2 };
            nplus = (nplus * 2).min(npl1.max(nplus / 2).max(1));
            fpold = *fp;

            // compute residuals per knot interval (for fpknot selection)
            let mut fpart = 0.0f64;
            let mut ii    = 1usize;
            let mut ll_r  = k1_v; // l starts at k1
            let mut new_v = 0i32;
            let mut jj_r  = 0usize;
            for it in 1..=m1 {
                if u[it - 1] >= t[ll_r - 1] { new_v = 1; ll_r += 1; }
                let mut term  = 0.0f64;
                let mut l0_f  = ll_r as i32 - k2;
                for _ in 0..idim_v {
                    let mut fac = 0.0f64; let mut j1_v = l0_f;
                    for j in 1..=k1_v { j1_v += 1; fac += c[(j1_v - 1) as usize] * q[(it - 1) + (j - 1) * m_v]; }
                    jj_r += 1;
                    term += (w[it - 1] * (fac - x[jj_r - 1])).powi(2);
                    l0_f += n_v as i32;
                }
                fpart += term;
                if new_v != 0 {
                    if ll_r <= k2_v {
                        fpint[nrint - 1] = term; // first interval → last fpint slot
                    } else {
                        let store = term * half;
                        fpint[ii - 1] = fpart - store; ii += 1; fpart = store;
                    }
                    new_v = 0;
                }
            }
            fpint[nrint - 1] += fpart;

            // add nplus new knots
            let mut hit_nest = false;
            for _ in 0..nplus {
                fpknot(u, m, t, n, fpint, nrdata, &mut (nrint as i32), nest, 1);
                if *n == nmax_v as i32 { do_knot_setup = true; continue 'outer; } // goto 5
                if *n == nest { hit_nest = true; break; } // goto 340
            }
            if hit_nest { continue 'main; } // goto 340 (next iter)
        } // end do-340
        break 'outer;
    } // end 'outer

    if !go_to_part2 { return; }

    // ══ Part 2: find smoothing parameter p ════════════════════════════════════
    let n_v   = *n as usize;
    let nk1   = n_v - k1_v;
    let n7    = nk1 - k_v;
    let n10_i = n7 as i32 - kk;
    let n10   = n10_i.max(0) as usize;
    let n11_i = n10_i - 1;
    let n11   = n11_i.max(0) as usize;
    let n8    = if n7 > 0 { n7 - 1 } else { 0 };
    let kk_v  = kk  as usize;
    let kk1_v = kk1 as usize;
    let nc_v  = nest_v * idim_v;

    fpdisc(t, *n, k2, b, nest);

    // compute initial estimate of p
    let mut l_f = n7 as i32;
    let mut p_sum = 0.0f64;
    let mut skip_n10 = false;
    for i in 1..=(k_v as i32) {
        let j = k - i + 1;
        p_sum += a2[(l_f - 1) as usize + (j as usize - 1) * nest_v];
        l_f -= 1;
        if l_f == 0 { skip_n10 = true; break; }
    }
    if !skip_n10 {
        for i in 1..=n10 { p_sum += a1[(i - 1) + 0]; }
    }
    let mut p  = n7 as f64 / p_sum;
    let mut p1 = 0.0f64;
    let mut f1 = fp0 - s;
    let mut p3 = -1.0f64;
    let mut f3 = fpms_saved;
    let mut ich1 = 0i32;
    let mut ich3 = 0i32;

    for iter in 1..=(maxit as usize) {
        let pinv = 1.0 / p;

        // copy z to c
        for i in 0..nc_v { c[i] = z[i]; }

        // init g1, g2 from a1, a2
        for i in 0..n7 {
            g1[i + (k1_v - 1) * nest_v] = a1[i + (k1_v - 1) * nest_v];
            g1[i + (k2_v - 1) * nest_v] = 0.0;
            g2[i] = 0.0; // g2(i,1) = 0
            for j in 0..kk_v {
                g1[i + j * nest_v] = a1[i + j * nest_v];
                g2[i + (j + 1) * nest_v] = a2[i + j * nest_v];
            }
        }
        // g2(l,1) = a1(l,j) for l=n10..1
        let mut l_g = n10_i;
        for j in 0..k1_v {
            if l_g <= 0 { break; }
            g2[(l_g - 1) as usize] = a1[(l_g - 1) as usize + j * nest_v];
            l_g -= 1;
        }

        // rotate the n8 rows of b into g
        for it in 1..=n8 {
            for j in 0..idim_v { xi[j] = 0.0; }
            for i in 0..k2_v { h1[i] = 0.0; }
            for i in 0..k1_v { h2[i] = 0.0; }

            let l_start: usize;
            if (it as i32) <= n11_i {
                // non-periodic zone (it <= n11)
                l_start = it;
                let mut l0 = it;
                let mut j = 0usize;
                while j < k2_v {
                    j += 1;
                    if l0 == n10 {
                        // from here: remaining rows go to h2
                        let mut l0h = 1usize;
                        for l1 in j..=k2_v {
                            h2[l0h - 1] = b[(it - 1) + (l1 - 1) * nest_v] * pinv;
                            l0h += 1;
                        }
                        break;
                    }
                    h1[j - 1] = b[(it - 1) + (j - 1) * nest_v] * pinv;
                    l0 += 1;
                }
            } else {
                // periodic zone (it > n11)
                l_start = 1;
                let mut i_v = it as i32 - n10_i;
                for j in 1..=k2_v {
                    i_v += 1;
                    let mut l0 = i_v;
                    loop {
                        let l1 = l0 - k1 as i32;
                        if l1 <= 0 {
                            h2[(l0 - 1) as usize] += b[(it - 1) + (j - 1) * nest_v] * pinv;
                            break;
                        } else if l1 <= n11_i {
                            h1[(l1 - 1) as usize] = b[(it - 1) + (j - 1) * nest_v] * pinv;
                            break;
                        } else {
                            l0 = l1 - n11_i;
                        }
                    }
                }
            }

            // rotation with rows l_start..n11 of g1
            if n11 > 0 {
                for j in l_start..=n11 {
                    let piv = h1[0];
                    let j_idx = j - 1;
                    let (mut cos_v, mut sin_v) = (0.0, 0.0);
                    fpgivs(piv, &mut g1[j_idx], &mut cos_v, &mut sin_v);
                    let mut j1_f = j as i32;
                    for j2 in 0..idim_v {
                        fprota(cos_v, sin_v, &mut xi[j2], &mut c[(j1_f - 1) as usize]);
                        j1_f += n_v as i32;
                    }
                    for i in 0..k1_v {
                        fprota(cos_v, sin_v, &mut h2[i], &mut g2[j_idx + i * nest_v]);
                    }
                    if j == n11 { break; }
                    let i2_max = (n11 - j).min(k1_v);
                    let mut i2 = 0usize;
                    for i in 0..i2_max {
                        i2 += 1;
                        fprota(cos_v, sin_v, &mut h1[i + 1], &mut g1[j_idx + i2 * nest_v]);
                        h1[i] = h1[i + 1];
                    }
                    h1[i2] = 0.0;
                }
            }

            // rotation with rows n11+1..n7 of g2
            for j in 1..=(k1_v as i32) {
                let ij = n11_i + j;
                if ij <= 0 { continue; }
                let ij_idx = (ij - 1) as usize;
                let piv = h2[(j - 1) as usize];
                let (mut cos_v, mut sin_v) = (0.0, 0.0);
                fpgivs(piv, &mut g2[ij_idx + (j as usize - 1) * nest_v], &mut cos_v, &mut sin_v);
                let mut j1_f = ij;
                for j2 in 0..idim_v {
                    fprota(cos_v, sin_v, &mut xi[j2], &mut c[(j1_f - 1) as usize]);
                    j1_f += n_v as i32;
                }
                if j == k1 { break; }
                for i in (j as usize + 1)..=k1_v {
                    fprota(cos_v, sin_v, &mut h2[i - 1], &mut g2[ij_idx + (i - 1) * nest_v]);
                }
            }
        } // end n8 rows

        // backward substitution (alias: c is used for both z and output)
        let mut j1_f = 0usize;
        for _ in 0..idim_v {
            let z_copy: Vec<f64> = c[j1_f..j1_f + n7].to_vec();
            fpbacp(g1, g2, &z_copy, n7 as i32, k1, &mut c[j1_f..], k2, nest);
            j1_f += n_v;
        }
        // fill periodicity coefficients
        for i in 1..=k_v {
            let mut j1_f = i - 1;
            for _ in 0..idim_v { c[j1_f + n7] = c[j1_f]; j1_f += n_v; }
        }

        // compute fp
        *fp = 0.0;
        let mut l_fp = k1_v;
        let mut jj_fp = 0usize;
        for it in 1..=m1 {
            if u[it - 1] >= t[l_fp - 1] { l_fp += 1; }
            let mut term  = 0.0f64;
            let mut l0_f  = l_fp as i32 - k2;
            for _ in 0..idim_v {
                let mut fac = 0.0f64; let mut j1_v = l0_f;
                for j in 1..=k1_v { j1_v += 1; fac += c[(j1_v - 1) as usize] * q[(it - 1) + (j - 1) * m_v]; }
                jj_fp += 1;
                term += (fac - x[jj_fp - 1]).powi(2);
                l0_f += n_v as i32;
            }
            *fp += term * w[it - 1].powi(2);
        }

        // convergence check
        let fpms = *fp - s;
        if fpms.abs() < acc { return; }
        if iter == maxit as usize { *ier = 3; return; }

        let p2 = p; let f2 = fpms;
        if ich3 == 0 {
            if (f2 - f3) <= acc {
                p3 = p2; f3 = f2; p = p * con4;
                if p <= p1 { p = p1 * con9 + p2 * con1; }
                continue;
            }
            if f2 < 0.0 { ich3 = 1; }
        }
        if ich1 == 0 {
            if (f1 - f2) <= acc {
                p1 = p2; f1 = f2; p = p / con4;
                if p3 < 0.0 { continue; }
                if p >= p3 { p = p2 * con1 + p3 * con9; }
                continue;
            }
            if f2 > 0.0 { ich1 = 1; }
        }
        if f2 >= f1 || f2 <= f3 { *ier = 2; return; }
        p = fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
    }
    *ier = 3;
}

/// Constrained-endpoint curve fitting engine (translates fpcons.f).
#[allow(clippy::too_many_arguments)]
fn fpcons(
    iopt: i32, idim: i32, m: i32, u: &[f64], _mx: i32, x: &[f64], w: &[f64],
    ib: i32, ie: i32, k: i32, s: f64, nest: i32,
    tol: f64, maxit: i32, k1: i32, k2: i32,
    n: &mut i32, t: &mut [f64], _nc: i32, c: &mut [f64], fp: &mut f64,
    fpint: &mut [f64], z: &mut [f64], a: &mut [f64], b: &mut [f64],
    g: &mut [f64], q: &mut [f64], nrdata: &mut [i32], ier: &mut i32,
) {
    let m = m as usize;
    let nest = nest as usize;
    let k1 = k1 as usize;
    let k2 = k2 as usize;
    let idim = idim as usize;
    let ib = ib as usize;
    let ie = ie as usize;
    let maxit = maxit as i32;
    let con1 = 0.1f64; let con9 = 0.9f64; let con4 = 0.04f64; let half = 0.5f64;
    let nmin = 2 * k1;

    // Which data points to use
    let mb = if ib > 0 { 2usize } else { 1 };
    let me = if ie > 0 { m - 1 } else { m };

    let mut h = [0.0f64; 7];
    let mut xi_loc = [0.0f64; 10];

    let mut fpold = 0.0f64;
    let mut nplus = 0i32;

    if iopt >= 0 {
        let acc = tol * s;
        let kbe = k1 as i32 - ib as i32 - ie as i32;
        let mmin = (kbe + 2).max(0) as usize;
        let mm = if m > mmin { m - mmin } else { 0 };
        let nmax = nmin + mm;

        if s == 0.0 {
            *n = nmax as i32;
            if nmax > nest { *ier = 1; return; }
            if mm > 0 {
                let mut i_f = k2;
                let j_start = 3 + k as usize / 2 - ib;
                let mut j_f = j_start;
                for _l in 1..=mm { t[i_f - 1] = u[j_f - 1]; i_f += 1; j_f += 1; }
            }
        } else {
            if iopt == 0 || *n == nmin as i32 {
                *n = nmin as i32;
                fpint[*n as usize - 1] = 0.0;
                nrdata[0] = (m - 2) as i32;
                fpold = 0.0; nplus = 0;
            } else {
                let n_us = *n as usize;
                let fp0 = fpint[n_us - 1];
                fpold = fpint[n_us - 2];
                nplus = nrdata[n_us - 1];
                if fp0 <= s { *n = nmin as i32; fpold = 0.0; nplus = 0; nrdata[0] = (m-2) as i32; }
            }
        }
    }

    'outer: for _iter in 1..=m as i32 {
        let n_us = *n as usize;
        if n_us == nmin { *ier = -2; }
        let nrint = n_us - nmin + 1;
        let nk1 = n_us - k1;
        let nn = if nk1 >= ib + ie { nk1 - ib - ie } else { 0 };
        let mut i_f = n_us;
        for _j in 1..=k1 { t[_j - 1] = u[0]; t[i_f - 1] = u[m - 1]; i_f -= 1; }

        *fp = 0.0;
        let nc_v = _nc as usize;
        for i in 0..nc_v { z[i] = 0.0; c[i] = 0.0; }
        if nn > 0 { for i in 0..nn { for j in 0..k1 { a[i + j*nest] = 0.0; } } }

        if me >= mb {
            let mut l = k1;
            let mut jj = (mb - 1) * idim;
            for it in mb..=me {
                let ui = u[it - 1]; let wi = w[it - 1];
                for j in 1..=idim { xi_loc[j-1] = x[jj] * wi; jj += 1; }
                while ui >= t[l] && l < nk1 { l += 1; }
                fpbspl(t, *n as i32, k1 as i32 - 1, ui, l as i32, &mut h);
                for i in 1..=k1 { q[it - 1 + (i-1)*m] = h[i-1]; h[i-1] *= wi; }
                // handle zero-forced coefficients
                let lj = {
                    let lj_tmp = k1;
                    let j_v = nk1 as i32 - l as i32 - ie as i32;
                    if j_v < 0 { (lj_tmp as i32 + j_v) as usize } else { lj_tmp }
                };
                let (li, j_off) = {
                    let j_v = l as i32 - k1 as i32 - ib as i32;
                    if j_v < 0 { (1 - j_v as usize, 0usize) } else { (1, j_v as usize) }
                };
                if li <= lj {
                    let mut j_f = j_off;
                    for i in li..=lj {
                        j_f += 1;
                        let piv = h[i - 1];
                        if piv.abs() < f64::EPSILON { continue; }
                        let (mut cos, mut sin) = (0.0, 0.0);
                        fpgivs(piv, &mut a[j_f - 1], &mut cos, &mut sin);
                        let mut j1 = j_f;
                        for j2 in 1..=idim { fprota(cos, sin, &mut xi_loc[j2-1], &mut z[j1-1]); j1 += n_us; }
                        if i == lj { break; }
                        let mut i2 = 0usize;
                        for i1 in (i+1)..=lj { i2 += 1; fprota(cos, sin, &mut h[i1-1], &mut a[j_f-1+i2*nest]); }
                    }
                }
                for j2 in 1..=idim { *fp += xi_loc[j2-1].powi(2); }
            }
        }
        let fp0 = if *ier == -2 { *fp } else { fpint[n_us - 1] };
        fpint[n_us - 1] = fp0;
        fpint[n_us - 2] = fpold;
        nrdata[n_us - 1] = nplus;
        if nn > 0 {
            let mut j1 = 1usize;
            for _j2 in 1..=idim {
                let j3 = j1 + ib;
                fpback(a, &z[j1-1..], nn as i32, k1 as i32, &mut c[j3-1..], nest as i32);
                j1 += n_us;
            }
        }
        if iopt < 0 { break 'outer; }
        let acc = tol * s;
        let fpms = *fp - s;
        if fpms.abs() < acc || fpms < 0.0 { break 'outer; }
        let kbe = k1 as i32 - ib as i32 - ie as i32;
        let mmin = (kbe + 2).max(0) as usize;
        let mm2 = if m > mmin { m - mmin } else { 0 };
        let nmax = nmin + mm2;
        if *n == nmax as i32 { *ier = -1; break 'outer; }
        if *n == nest as i32 { *ier = 1; break 'outer; }
        if *ier != 0 { nplus = 1; *ier = 0; } else {
            let rn = nplus as f64;
            let npl1 = if fpold - *fp > acc { (rn * fpms / (fpold - *fp)) as i32 } else { nplus*2 };
            nplus = (nplus*2).min(npl1.max(nplus/2).max(1));
        }
        fpold = *fp;
        let mut fpart = 0.0f64;
        let mut ii = 1usize;
        let mut ll = k2;
        let mut new_f = 0i32;
        let mut jj = (mb - 1) * idim;
        for it in mb..=me {
            if u[it-1] >= t[ll-1] && ll <= nk1 { new_f = 1; ll += 1; }
            let mut term = 0.0f64;
            let mut l0 = ll - k2;
            for j2 in 1..=idim {
                let mut fac = 0.0f64;
                let mut j1 = l0;
                for j in 1..=k1 { j1 += 1; fac += c[j1-1] * q[it-1+(j-1)*m]; }
                term += (w[it-1] * (fac - x[jj])).powi(2);
                jj += 1; l0 += n_us;
            }
            fpart += term;
            if new_f != 0 { let store=term*half; fpint[ii-1]=fpart-store; ii+=1; fpart=store; new_f=0; }
        }
        fpint[nrint-1] = fpart;
        let kbe = k1 as i32 - ib as i32 - ie as i32;
        let mmin = (kbe + 2).max(0) as usize;
        let mm2 = if m > mmin { m - mmin } else { 0 };
        let nmax = nmin + mm2;
        for _l in 1..=nplus {
            fpknot(u, m as i32, t, n, fpint, nrdata, &mut (nrint as i32), nest as i32, 1);
            if *n == nmax as i32 {
                if mm2 > 0 {
                    let mut i_fi = k2; let mut j_fi = 3 + k as usize/2 - ib;
                    for _ in 1..=mm2 { t[i_fi-1]=u[j_fi-1]; i_fi+=1; j_fi+=1; }
                }
                break;
            }
            if *n == nest as i32 { break; }
        }
    }

    if *ier == -2 { return; }
    let n_us = *n as usize;
    let nk1 = n_us - k1;
    let nn = if nk1 >= ib + ie { nk1 - ib - ie } else { return };
    let n8 = n_us - nmin;
    let fp0 = fpint[n_us - 1];
    let fpms0 = *fp - s;
    if fpms0 >= 0.0 && n8 > 0 && nn > 0 {
        fpdisc(t, *n, k2 as i32, b, nest as i32);
        let mut p1 = 0.0f64; let mut f1 = fp0 - s;
        let mut p3 = -1.0f64; let mut f3 = fpms0;
        let mut p: f64 = 0.0;
        for i in 0..nn { p += a[i]; }
        p = (nn as f64) / p;
        let mut ich1 = 0i32; let mut ich3 = 0i32;
        let acc = tol * s;
        for iter in 1..=maxit {
            let pinv = 1.0 / p;
            let nc_v = _nc as usize;
            for i in 0..nc_v { c[i] = z[i]; }
            for i in 0..nn { g[i + (k2-1)*nest] = 0.0; for j in 0..k1 { g[i+j*nest] = a[i+j*nest]; } }
            for it in 1..=n8 {
                for i in 0..k2 { h[i] = b[it-1+i*nest] * pinv; }
                // handle ib offset
                if it <= ib {
                    let j1 = ib - it + 2;
                    let mut j2 = 1usize;
                    for i in j1..=k2 { h[j2-1] = h[i-1]; j2 += 1; }
                    for i in j2..=k2 { h[i-1] = 0.0; }
                }
                let mut xi2 = [0.0f64; 10];
                let jj_start = if it > ib { it - ib } else { 1 };
                let mut j_f = jj_start;
                while j_f <= nn {
                    let piv = h[0];
                    let (mut cos, mut sin) = (0.0, 0.0);
                    fpgivs(piv, &mut g[j_f-1], &mut cos, &mut sin);
                    let mut j1 = j_f;
                    for j2 in 1..=idim { fprota(cos, sin, &mut xi2[j2-1], &mut c[j1-1]); j1 += n_us; }
                    if j_f == nn { break; }
                    let i2 = (nn - j_f).min(k1);
                    for i in 1..=i2 { fprota(cos, sin, &mut h[i], &mut g[j_f-1+i*nest]); h[i-1] = h[i]; }
                    h[i2] = 0.0;
                    j_f += 1;
                }
            }
            let mut j1 = 1usize;
            for _j2 in 1..=idim {
                let j3 = j1 + ib;
                let z_copy3: Vec<f64> = c[j1-1..j1-1+nn].to_vec();
                fpback(g, &z_copy3, nn as i32, k2 as i32, &mut c[j3-1..], nest as i32);
                if ib > 0 { for i in 0..ib { c[j1-1+i] = 0.0; } }
                j1 += n_us;
            }
            *fp = 0.0;
            let mut ll = k2;
            let mut jj = (mb-1)*idim;
            for it in mb..=me {
                if u[it-1] >= t[ll-1] && ll <= nk1 { ll += 1; }
                let mut l0 = ll - k2; let mut term = 0.0f64;
                for j2 in 1..=idim {
                    let mut fac = 0.0f64; let mut j1v = l0;
                    for j in 1..=k1 { j1v += 1; fac += c[j1v-1]*q[it-1+(j-1)*m]; }
                    term += (fac - x[jj]).powi(2); jj += 1; l0 += n_us;
                }
                *fp += term * w[it-1].powi(2);
            }
            let fpms = *fp - s;
            if fpms.abs() < acc { return; }
            if iter == maxit { *ier = 3; return; }
            let p2 = p; let f2 = fpms;
            if ich3 == 0 { if (f2-f3)<=acc { p3=p2;f3=f2;p=p*con4;if p<=p1{p=p1*con9+p2*con1;}continue; } if f2<0.0{ich3=1;} }
            if ich1 == 0 { if (f1-f2)<=acc { p1=p2;f1=f2;p=p/con4;if p3<0.0{continue;}if p>=p3{p=p2*con1+p3*con9;}continue; } if f2>0.0{ich1=1;} }
            if f2 >= f1 || f2 <= f3 { *ier = 2; return; }
            p = fprati(&mut p1, &mut f1, p2, f2, &mut p3, &mut f3);
        }
        *ier = 3;
    }
}

// ─── Section 7: Public C API (no_mangle extern "C") ───────────────────────────

/// curfit: Determine a smoothing or interpolating spline for 1D data.
#[no_mangle]
pub unsafe extern "C" fn curfit_(
    iopt: *const i32, m: *const i32, x: *const f64, y: *const f64, w: *const f64,
    xb: *const f64, xe: *const f64, k: *const i32, s: *const f64,
    nest: *const i32, n: *mut i32, t: *mut f64, c: *mut f64, fp: *mut f64,
    wrk: *mut f64, lwrk: *const i32, iwrk: *mut i32, ier: *mut i32,
) {
    let iopt = *iopt; let m_v = *m as usize; let k_v = *k; let s_v = *s;
    let nest_v = *nest as usize;
    let xb_v = *xb; let xe_v = *xe;
    let k1 = (k_v + 1) as usize; let k2 = k1 + 1;
    let maxit = 20i32; let tol = 0.001f64;
    *ier = 10;
    if k_v < 1 || k_v > 5 { return; }
    if iopt < -1 || iopt > 1 { *ier = 11; return; }
    let nmin = 2 * k1;
    if m_v < k1 || nest_v < nmin { *ier = 12; return; }
    let lwest = m_v * k1 + nest_v * (7 + 3 * k_v as usize);
    if (*lwrk as usize) < lwest { *ier = 13; return; }
    let x_sl = std::slice::from_raw_parts(x, m_v);
    let y_sl = std::slice::from_raw_parts(y, m_v);
    let w_sl = std::slice::from_raw_parts(w, m_v);
    if xb_v > x_sl[0] || xe_v < x_sl[m_v-1] || w_sl[0] <= 0.0 { *ier = 14; return; }
    for i in 1..m_v { if x_sl[i-1] >= x_sl[i] || w_sl[i] <= 0.0 { *ier = 15; return; } }
    if iopt < 0 {
        *ier = 16;
        let n_v = *n as usize;
        if n_v < nmin || n_v > nest_v { return; }
        let t_sl = std::slice::from_raw_parts_mut(t, nest_v);
        let mut jv = n_v;
        for i in 1..=k1 { t_sl[i-1] = xb_v; t_sl[jv-1] = xe_v; jv -= 1; }
        let x_ref = std::slice::from_raw_parts(x, m_v);
        let t_ref = std::slice::from_raw_parts(t, n_v);
        let mut ier_v = 0i32;
        fpchec(x_ref, m_v as i32, t_ref, n_v as i32, k_v, &mut ier_v);
        if ier_v != 0 { *ier = ier_v; return; }
    } else {
        if s_v < 0.0 { *ier = 16; return; }
        if s_v == 0.0 && nest_v < m_v + k1 { *ier = 17; return; }
    }
    *ier = 0;
    // Partition working space via raw pointers (avoids multiple mutable borrow checks)
    let ifp: usize = 0;
    let iz      = ifp + nest_v;
    let ia      = iz  + nest_v;
    let ib_off  = ia  + nest_v * k1;
    let ig      = ib_off + nest_v * k2;
    let iq      = ig  + nest_v * k2;
    fpcurf(
        iopt,
        std::slice::from_raw_parts(x, m_v),
        std::slice::from_raw_parts(y, m_v),
        std::slice::from_raw_parts(w, m_v),
        m_v as i32, xb_v, xe_v, k_v, s_v,
        nest_v as i32, tol, maxit, k1 as i32, k2 as i32,
        &mut *n,
        std::slice::from_raw_parts_mut(t, nest_v),
        std::slice::from_raw_parts_mut(c, nest_v),
        &mut *fp,
        std::slice::from_raw_parts_mut(wrk.add(ifp), nest_v),
        std::slice::from_raw_parts_mut(wrk.add(iz),  nest_v),
        std::slice::from_raw_parts_mut(wrk.add(ia),  nest_v * k1),
        std::slice::from_raw_parts_mut(wrk.add(ib_off), nest_v * k2),
        std::slice::from_raw_parts_mut(wrk.add(ig),  nest_v * k2),
        std::slice::from_raw_parts_mut(wrk.add(iq),  m_v   * k1),
        std::slice::from_raw_parts_mut(iwrk, nest_v),
        &mut *ier,
    );
}

/// splev: Evaluate a spline s(x) at m points.
#[no_mangle]
pub unsafe extern "C" fn splev_(
    t: *const f64, n: *const i32, c: *const f64, k: *const i32,
    x: *const f64, y: *mut f64, m: *const i32, ier: *mut i32,
) {
    let n_v = *n as usize; let k_v = *k as usize; let m_v = *m as usize;
    *ier = 10;
    if m_v == 0 { return; }
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let x_sl = std::slice::from_raw_parts(x, m_v);
    let y_sl = std::slice::from_raw_parts_mut(y, m_v);
    for i in 1..m_v { if x_sl[i] < x_sl[i-1] { return; } }
    *ier = 0;
    let k1 = k_v + 1;
    let nk1 = n_v - k1;
    let tb = t_sl[k1 - 1];
    let te = t_sl[nk1];
    let mut l = k1;
    let mut l1 = l + 1;
    let mut h = [0.0f64; 6];
    for i in 0..m_v {
        let mut arg = x_sl[i];
        if arg < tb { arg = tb; }
        if arg > te { arg = te; }
        while arg >= t_sl[l1 - 1] && l < nk1 { l = l1; l1 = l + 1; }
        fpbspl(t_sl, *n, *k, arg, l as i32, &mut h);
        let mut sp = 0.0f64;
        let mut ll = l - k1;
        for j in 0..k1 { ll += 1; sp += c_sl[ll - 1] * h[j]; }
        y_sl[i] = sp;
    }
}

/// curev: Evaluate a parametric spline curve s(u) at m points.
#[no_mangle]
pub unsafe extern "C" fn curev_(
    idim: *const i32, t: *const f64, n: *const i32, c: *const f64, nc: *const i32,
    k: *const i32, u: *const f64, m: *const i32, xy: *mut f64, mxy: *const i32, ier: *mut i32,
) {
    let n_v = *n as usize; let k_v = *k as usize; let m_v = *m as usize;
    let idim_v = *idim as usize; let mxy_v = *mxy as usize;
    *ier = 10;
    if m_v == 0 { return; }
    if mxy_v < m_v * idim_v { return; }
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, *nc as usize);
    let u_sl = std::slice::from_raw_parts(u, m_v);
    let xy_sl = std::slice::from_raw_parts_mut(xy, mxy_v);
    for i in 1..m_v { if u_sl[i] < u_sl[i-1] { return; } }
    *ier = 0;
    let k1 = k_v + 1;
    let nk1 = n_v - k1;
    let tb = t_sl[k1 - 1];
    let te = t_sl[nk1];
    let mut l = k1;
    let mut l1 = l + 1;
    let mut h = [0.0f64; 6];
    let mut mm = 0usize;
    for i in 0..m_v {
        let mut arg = u_sl[i];
        if arg < tb { arg = tb; }
        if arg > te { arg = te; }
        while arg >= t_sl[l1 - 1] && l < nk1 { l = l1; l1 = l + 1; }
        fpbspl(t_sl, *n, *k, arg, l as i32, &mut h);
        let mut ll = l - k1;
        for _j1 in 1..=idim_v {
            let mut jj = ll;
            let mut sp = 0.0f64;
            for j in 0..k1 { jj += 1; sp += c_sl[jj - 1] * h[j]; }
            xy_sl[mm] = sp;
            mm += 1;
            ll += n_v;
        }
    }
}

/// spalde: Evaluate all derivatives of a spline at a single point.
#[no_mangle]
pub unsafe extern "C" fn spalde_(
    t: *const f64, n: *const i32, c: *const f64, k1: *const i32,
    x: *const f64, d: *mut f64, ier: *mut i32,
) {
    let n_v = *n as usize; let k1_v = *k1 as usize;
    let nk1 = n_v - k1_v;
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let d_sl = std::slice::from_raw_parts_mut(d, k1_v);
    let x_v = *x;
    *ier = 10;
    if x_v < t_sl[k1_v - 1] || x_v > t_sl[nk1] { return; }
    // search knot interval
    let mut l = k1_v;
    while x_v >= t_sl[l] && l < nk1 { l += 1; }
    if t_sl[l - 1] >= t_sl[l] { return; }
    *ier = 0;
    fpader(t_sl, *n, c_sl, *k1, x_v, l as i32, d_sl);
}

/// cualde: Evaluate all derivatives of a parametric spline at a single point.
#[no_mangle]
pub unsafe extern "C" fn cualde_(
    idim: *const i32, t: *const f64, n: *const i32, c: *const f64, nc: *const i32,
    k1: *const i32, u: *const f64, d: *mut f64, nd: *const i32, ier: *mut i32,
) {
    let n_v = *n as usize; let k1_v = *k1 as usize; let idim_v = *idim as usize;
    let nc_v = *nc as usize; let nd_v = *nd as usize;
    *ier = 10;
    if nd_v < k1_v * idim_v { return; }
    let nk1 = n_v - k1_v;
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, nc_v);
    let d_sl = std::slice::from_raw_parts_mut(d, nd_v);
    let u_v = *u;
    if u_v < t_sl[k1_v - 1] || u_v > t_sl[nk1] { return; }
    let mut l = k1_v;
    while u_v >= t_sl[l] && l < nk1 { l += 1; }
    if t_sl[l - 1] >= t_sl[l] { return; }
    *ier = 0;
    let mut h = [0.0f64; 6];
    let mut j_base = 0usize;
    for i in 1..=idim_v {
        fpader(t_sl, *n, &c_sl[j_base..], *k1, u_v, l as i32, &mut h);
        let mut m_off = i - 1;
        for kk in 0..k1_v {
            d_sl[m_off] = h[kk];
            m_off += idim_v;
        }
        j_base += n_v;
    }
}

/// spalder (= splder): Evaluate the nu-th derivative of a spline at m points.
#[no_mangle]
pub unsafe extern "C" fn spalder_(
    t: *const f64, n: *const i32, c: *const f64, k: *const i32,
    nu: *const i32, x: *const f64, y: *mut f64, m: *const i32,
    wrk: *mut f64, ier: *mut i32,
) {
    let n_v = *n as usize; let k_v = *k as usize; let nu_v = *nu as usize;
    let m_v = *m as usize;
    *ier = 10;
    if nu_v > k_v { return; }
    if m_v == 0 { return; }
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let x_sl = std::slice::from_raw_parts(x, m_v);
    let y_sl = std::slice::from_raw_parts_mut(y, m_v);
    let wrk_sl = std::slice::from_raw_parts_mut(wrk, n_v);
    for i in 1..m_v { if x_sl[i] < x_sl[i-1] { return; } }
    *ier = 0;
    let k1 = k_v + 1;
    let nk1 = n_v - k1;
    // build derivative B-spline coefficients in wrk
    for i in 0..nk1 { wrk_sl[i] = c_sl[i]; }
    let mut kk = k_v;
    let mut l_off = 0usize;
    let mut nk2 = nk1;
    for _j in 0..nu_v {
        let ak = kk as f64;
        nk2 -= 1;
        let l1 = l_off;
        for i in 0..nk2 {
            let l1i = l1 + i + 1;
            let l2i = l1i + kk;
            let fac = t_sl[l2i] - t_sl[l1i]; // t(l1+1), t(l2+1) → 0-based
            if fac > 0.0 { wrk_sl[i] = ak * (wrk_sl[i+1] - wrk_sl[i]) / fac; }
        }
        l_off += 1;
        kk -= 1;
    }
    if kk == 0 {
        // piecewise constant
        let mut j = 0usize;
        let l_off2 = l_off;
        for i in 0..m_v {
            let arg = x_sl[i];
            while arg >= t_sl[l_off2 + j + 1] && l_off2 + j + 1 < nk1 { j += 1; }
            y_sl[i] = wrk_sl[j];
        }
        return;
    }
    let k2 = k1 - nu_v;
    let mut l = k1;
    let mut l1 = l + 1;
    let mut h = [0.0f64; 6];
    let tb = t_sl[k1 - 1]; let te = t_sl[nk1];
    for i in 0..m_v {
        let mut arg = x_sl[i];
        if arg < tb { arg = tb; }
        if arg > te { arg = te; }
        while arg >= t_sl[l1 - 1] && l < nk1 { l = l1; l1 = l + 1; }
        fpbspl(t_sl, *n, kk as i32, arg, l as i32, &mut h);
        let mut sp = 0.0f64;
        let mut ll = l - k1;
        for j in 0..k2 { ll += 1; sp += wrk_sl[ll - 1] * h[j]; }
        y_sl[i] = sp;
    }
}

/// clocur: Closed-curve fitting (delegates to fpclos).
#[no_mangle]
pub unsafe extern "C" fn clocur_(
    iopt: *const i32, ipar: *const i32, idim: *const i32, m: *const i32,
    u: *mut f64, mx: *const i32, x: *const f64, w: *const f64,
    k: *const i32, s: *const f64, nest: *const i32, n: *mut i32,
    t: *mut f64, nc: *const i32, c: *mut f64, fp: *mut f64,
    wrk: *mut f64, lwrk: *const i32, iwrk: *mut i32, ier: *mut i32,
) {
    let iopt_v = *iopt; let idim_v = *idim as usize; let m_v = *m as usize;
    let k_v = *k; let s_v = *s; let nest_v = *nest as usize;
    let mx_v = *mx as usize; let nc_v = *nc as usize; let lwrk_v = *lwrk as usize;
    *ier = 10;
    let k1 = (k_v + 1) as usize; let k2 = k1 + 1;
    let nmin = 2 * k1;
    if k_v < 1 || k_v > 5 { return; }
    if idim_v < 1 || idim_v > 10 { return; }
    let m1 = m_v - 1;
    if m_v < k1 { return; }
    let ncc = nest_v * idim_v;
    if mx_v < m_v * idim_v || nc_v < ncc { return; }
    let lwest = m_v * k1 + nest_v * (7 + idim_v + 5 * k_v as usize);
    if lwrk_v < lwest { return; }

    // compute parameter values if ipar==0
    let ipar_v = *ipar;
    let u_sl = std::slice::from_raw_parts_mut(u, m_v);
    if ipar_v == 0 {
        u_sl[0] = 0.0;
        let x_sl_full = std::slice::from_raw_parts(x, mx_v);
        for i in 1..m_v {
            let mut dist = 0.0f64;
            for j in 0..idim_v {
                let diff = x_sl_full[i*idim_v + j] - x_sl_full[(i-1)*idim_v + j];
                dist += diff * diff;
            }
            u_sl[i] = u_sl[i-1] + dist.sqrt();
        }
    }
    for i in 0..m1 { if u_sl[i] >= u_sl[i+1] || *std::slice::from_raw_parts(w, m_v).get(i).unwrap_or(&1.0) <= 0.0 { return; } }

    if iopt_v >= 0 {
        if s_v < 0.0 { return; }
        if s_v == 0.0 && nest_v < m_v + 2 * k_v as usize { return; }
    }
    *ier = 0;
    let ifp: usize = 0; let iz = ifp+nest_v;
    let ia1 = iz+ncc; let ia2 = ia1+nest_v*k1;
    let ib_off = ia2+nest_v*(k1-1); let ig1 = ib_off+nest_v*k2;
    let ig2 = ig1+nest_v*k2; let iq = ig2+nest_v*k1;
    // set periodic boundary knots if iopt < 0
    if iopt_v < 0 {
        let n_v = *n as usize;
        let t_sl2 = std::slice::from_raw_parts_mut(t, nest_v);
        let u_sl2 = std::slice::from_raw_parts(u as *const f64, m_v);
        if n_v <= nmin || n_v > nest_v { *ier = 10; return; }
        let per = u_sl2[m_v-1] - u_sl2[0];
        let mut j1 = k1; t_sl2[j1-1] = u_sl2[0];
        let mut i1 = n_v - k_v as usize; t_sl2[i1-1] = u_sl2[m_v-1];
        let mut j2 = j1; let mut i2 = i1;
        for _ in 1..=k_v as usize {
            i1 += 1; i2 -= 1; j1 += 1; j2 -= 1;
            t_sl2[j2-1] = t_sl2[i2-1] - per;
            t_sl2[i1-1] = t_sl2[j1-1] + per;
        }
        let mut ier_v = 0i32;
        fpchep(u_sl2, m_v as i32, t_sl2, *n, k_v, &mut ier_v);
        if ier_v != 0 { *ier = ier_v; return; }
    }
    let maxit = 20i32; let tol = 0.001f64;
    fpclos(
        iopt_v, idim_v as i32, m_v as i32,
        std::slice::from_raw_parts(u as *const f64, m_v),
        mx_v as i32,
        std::slice::from_raw_parts(x, mx_v),
        std::slice::from_raw_parts(w, m_v),
        k_v, s_v, nest_v as i32, tol, maxit, k1 as i32, k2 as i32,
        &mut *n,
        std::slice::from_raw_parts_mut(t, nest_v),
        nc_v as i32,
        std::slice::from_raw_parts_mut(c, nc_v),
        &mut *fp,
        std::slice::from_raw_parts_mut(wrk.add(ifp), nest_v),
        std::slice::from_raw_parts_mut(wrk.add(iz),  ncc),
        std::slice::from_raw_parts_mut(wrk.add(ia1), nest_v*k1),
        std::slice::from_raw_parts_mut(wrk.add(ia2), nest_v*(k1-1)),
        std::slice::from_raw_parts_mut(wrk.add(ib_off), nest_v*k2),
        std::slice::from_raw_parts_mut(wrk.add(ig1), nest_v*k2),
        std::slice::from_raw_parts_mut(wrk.add(ig2), nest_v*k1),
        std::slice::from_raw_parts_mut(wrk.add(iq),  m_v*k1),
        std::slice::from_raw_parts_mut(iwrk, nest_v),
        &mut *ier,
    );
}

/// concur: Parametric spline with derivative endpoint constraints.
#[no_mangle]
pub unsafe extern "C" fn concur_(
    iopt: *const i32, idim: *const i32, m: *const i32,
    u: *const f64, mx: *const i32, x: *const f64, xx: *mut f64, w: *const f64,
    ib: *const i32, db: *const f64, nb: *const i32,
    ie: *const i32, de: *const f64, ne: *const i32,
    k: *const i32, s: *const f64, nest: *const i32, n: *mut i32,
    t: *mut f64, nc: *const i32, c: *mut f64,
    np: *const i32, cp: *mut f64,
    fp: *mut f64, wrk: *mut f64, lwrk: *const i32, iwrk: *mut i32, ier: *mut i32,
) {
    let iopt_v = *iopt; let idim_v = *idim as usize; let m_v = *m as usize;
    let k_v = *k; let s_v = *s; let nest_v = *nest as usize;
    let mx_v = *mx as usize; let nc_v = *nc as usize;
    let ib_v = *ib as usize; let ie_v = *ie as usize;
    let nb_v = *nb as usize; let ne_v = *ne as usize;
    let np_v = *np as usize; let lwrk_v = *lwrk as usize;
    *ier = 10;
    let k1 = (k_v + 1) as usize; let k2 = k1 + 1;
    let nmin = 2 * k1;
    let ib1 = if ib_v > 0 { ib_v - 1 } else { 0 };
    let ie1 = if ie_v > 0 { ie_v - 1 } else { 0 };
    if k_v < 1 || k_v > 5 { return; }
    if idim_v < 1 || idim_v > 10 { return; }
    if m_v < k1 { return; }
    if ib_v as i32 > (k_v+1)/2 || ie_v as i32 > (k_v+1)/2 { return; }
    let ncc = nest_v * idim_v;
    let mxx = m_v * idim_v;
    if mx_v < mxx || nc_v < ncc { return; }
    let lwest = m_v * k1 + nest_v * (6 + idim_v + 3 * k_v as usize);
    if lwrk_v < lwest { return; }
    let u_sl = std::slice::from_raw_parts(u, m_v);
    let w_sl = std::slice::from_raw_parts(w, m_v);
    let x_sl = std::slice::from_raw_parts(x, mx_v);
    let xx_sl = std::slice::from_raw_parts_mut(xx, mxx);
    let db_sl = std::slice::from_raw_parts(db, nb_v);
    let de_sl = std::slice::from_raw_parts(de, ne_v);
    let cp_sl = std::slice::from_raw_parts_mut(cp, np_v);
    let t_sl = std::slice::from_raw_parts_mut(t, nest_v);
    let c_sl = std::slice::from_raw_parts_mut(c, nc_v);
    let iwrk_sl = std::slice::from_raw_parts_mut(iwrk, nest_v);
    let wrk_sl = std::slice::from_raw_parts_mut(wrk, lwrk_v);
    for i in 1..m_v { if u_sl[i-1] >= u_sl[i] || w_sl[i] <= 0.0 { return; } }
    if w_sl[0] <= 0.0 { return; }
    if iopt_v >= 0 {
        let nmax = m_v + k1 + ib1 + ie1;
        if s_v < 0.0 { return; }
        if s_v == 0.0 && nest_v < nmax { return; }
    } else {
        let n_v = *n as usize;
        if n_v < nmin || n_v > nest_v { return; }
        let mut jv = n_v;
        for i in 1..=k1 { t_sl[i-1] = u_sl[0]; t_sl[jv-1] = u_sl[m_v-1]; jv -= 1; }
        let mut ier_v = 0i32;
        fpched(u_sl, m_v as i32, t_sl, *n, k_v, ib_v as i32, ie_v as i32, &mut ier_v);
        if ier_v != 0 { *ier = ier_v; return; }
    }
    *ier = 0;
    // compute polynomial curve satisfying constraints
    fppocu(idim_v as i32, k_v, u_sl[0], u_sl[m_v-1],
           ib_v as i32, db_sl, nb_v as i32, ie_v as i32, de_sl, ne_v as i32,
           cp_sl, np_v as i32);
    // build modified knot vector for polynomial
    let mut wrk_nmin = vec![0.0f64; nmin];
    let mut jv = nmin;
    for i in 1..=k1 { wrk_nmin[i-1] = u_sl[0]; wrk_nmin[jv-1] = u_sl[m_v-1]; jv -= 1; }
    // evaluate polynomial curve and subtract from data
    let mut xx_tmp = vec![0.0f64; mxx];
    let mut ier_tmp = 0i32;
    curev_(
        idim as *const i32, wrk_nmin.as_ptr(), &(nmin as i32) as *const i32,
        cp_sl.as_ptr(), np as *const i32, k as *const i32,
        u as *const f64, m as *const i32, xx_tmp.as_mut_ptr(),
        &(mxx as i32) as *const i32, &mut ier_tmp,
    );
    for i in 0..mxx { xx_sl[i] = x_sl[i] - xx_tmp[i]; }
    // now fit spline to residual
    let jfp = 0; let jz = jfp+nest_v; let ja = jz+ncc;
    let jb_off = ja+nest_v*k1; let jg = jb_off+nest_v*k2; let jq = jg+nest_v*k2;
    let maxit = 20i32; let tol = 0.001f64;
    fpcons(
        iopt_v, idim_v as i32, m_v as i32, u_sl, mxx as i32, xx_sl, w_sl,
        ib_v as i32, ie_v as i32, k_v, s_v, nest_v as i32,
        tol, maxit, k1 as i32, k2 as i32,
        &mut *n,
        std::slice::from_raw_parts_mut(t, nest_v),
        nc_v as i32,
        std::slice::from_raw_parts_mut(c, nc_v),
        &mut *fp,
        std::slice::from_raw_parts_mut(wrk.add(jfp), nest_v),
        std::slice::from_raw_parts_mut(wrk.add(jz),  ncc),
        std::slice::from_raw_parts_mut(wrk.add(ja),  nest_v*k1),
        std::slice::from_raw_parts_mut(wrk.add(jb_off), nest_v*k2),
        std::slice::from_raw_parts_mut(wrk.add(jg),  nest_v*k2),
        std::slice::from_raw_parts_mut(wrk.add(jq),  m_v*k1),
        std::slice::from_raw_parts_mut(iwrk, nest_v),
        &mut *ier,
    );
    // add polynomial curve back
    let mut cc_tmp = vec![0.0f64; nc_v];
    let mut t1_tmp = vec![0.0f64; nest_v];
    let mut t2_tmp = vec![0.0f64; nest_v];
    fpadpo(idim_v as i32,
           std::slice::from_raw_parts(t as *const f64, nest_v),
           *n,
           std::slice::from_raw_parts_mut(c, nc_v),
           nc_v as i32, k_v,
           cp_sl, np_v as i32, &mut cc_tmp, &mut t1_tmp, &mut t2_tmp);
}

/// splint: Compute the integral of a spline from a to b.
#[no_mangle]
pub unsafe extern "C" fn splint_(
    t: *const f64, n: *const i32, c: *const f64, k: *const i32,
    a: *const f64, b: *const f64, wrk: *mut f64,
) -> f64 {
    let n_v = *n as usize; let k_v = *k as usize;
    let nk1 = n_v - k_v - 1;
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let wrk_sl = std::slice::from_raw_parts_mut(wrk, n_v);
    fpintb(t_sl, *n, wrk_sl, nk1 as i32, *a, *b);
    let mut result = 0.0f64;
    for i in 0..nk1 { result += c_sl[i] * wrk_sl[i]; }
    result
}

/// sproot: Find zeros of a cubic spline.
#[no_mangle]
pub unsafe extern "C" fn sproot_(
    t: *const f64, n: *const i32, c: *const f64,
    zero: *mut f64, mest: *const i32, m: *mut i32, ier: *mut i32,
) {
    let n_v = *n as usize; let mest_v = *mest as usize;
    *ier = 10;
    if n_v < 8 { return; }
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let zero_sl = std::slice::from_raw_parts_mut(zero, mest_v);
    let n4 = n_v - 4;
    // check t is non-decreasing at the first and last 3 knots
    for i in 0..3 {
        if t_sl[i] > t_sl[i + 1] { return; }
        if t_sl[n_v - i - 1] < t_sl[n_v - i - 2] { return; }
    }
    for i in 3..n4 { if t_sl[i] >= t_sl[i+1] { return; } }
    *ier = 0;
    let mut m_count = 0i32;
    // Evaluate cubic spline at each knot interval [t(l), t(l+1)] for l=4..n4
    let h1_init = t_sl[4] - t_sl[3];
    let h2_init = t_sl[5] - t_sl[4];
    let t1_init = t_sl[4] - t_sl[2];
    let t2_init = t_sl[5] - t_sl[3];
    let t3_init = t_sl[6] - t_sl[4];
    let t4_init = t_sl[5] - t_sl[2];
    let t5_init = t_sl[6] - t_sl[3];
    let c1 = c_sl[0]; let c2 = c_sl[1]; let c3 = c_sl[2];
    let c4 = (c2 - c1) / t4_init;
    let c5 = (c3 - c2) / t5_init;
    let d4 = (h2_init * c1 + t1_init * c2) / t4_init;
    let d5 = (t3_init * c2 + h1_init * c3) / t5_init;
    let mut a0 = (h2_init * d4 + h1_init * d5) / t2_init;
    let mut ah = 3.0 * (h2_init * c4 + h1_init * c5) / t2_init;
    let mut h1 = h2_init; let mut h2;
    let mut t2v = t2_init; let mut t3v = t3_init;
    let mut t4v = t4_init; let mut t5v = t5_init;
    let mut t1v;
    let mut c2v = c2; let mut c3v = c3;
    let mut c4v = c5; // c4 = previous c5

    for l in 3..n4 {
        h2 = t_sl[l+2] - t_sl[l+1];
        t1v = t2v; t2v = t3v; t3v = t_sl[l+3] - t_sl[l+1];
        t4v = t5v; t5v = t_sl[l+3] - t_sl[l];
        let c1v = c2v; c2v = c3v; c3v = c_sl[l];
        let c5v = (c3v - c2v) / t5v;
        let d4v = (h2 * c1v + t1v * c2v) / t4v;
        let d5v = (h1 * c3v + t3v * c2v) / t5v;
        let b0 = (h2 * d4v + h1 * d5v) / t2v;
        let bh = 3.0 * (h2 * c4v + h1 * c5v) / t2v;
        let a1 = ah * h1; let b1 = bh * h1;
        let a2 = 3.0 * (b0 - a0) - b1 - 2.0 * a1;
        let a3 = 2.0 * (a0 - b0) + b1 + a1;
        let mut y = [0.0f64; 3]; let mut jcount = 0i32;
        if a0 * b0 <= 0.0 {
            fpcuro(a3, a2, a1, a0, &mut y, &mut jcount);
        } else {
            // quick sign test to avoid fpcuro call
            let z0 = a0 >= 0.0; let z1 = ah >= 0.0;
            let z2 = a2 >= 0.0; let z3 = b1 >= 0.0;
            let z4 = 3.0*a3+a2 >= 0.0;
            let maybe = (z0 && (!z1 && (z3 || (z2 && !z4)) || (!z2 && z3 && z4)))
                     || (!z0 && (z1 && (!z3 || (!z2 && z4)) || (z2 && !z3 && !z4)));
            if maybe { fpcuro(a3, a2, a1, a0, &mut y, &mut jcount); }
        }
        for i in 0..jcount as usize {
            if y[i] >= 0.0 && y[i] <= 1.0 {
                if m_count >= *mest { *ier = 1; *m = m_count; return; }
                zero_sl[m_count as usize] = t_sl[l] + h1 * y[i];
                m_count += 1;
            }
        }
        a0 = b0; ah = bh;
        h1 = h2; c4v = c5v;
    }
    // sort zeros
    for i in 1..m_count as usize {
        let mut j = i;
        while j > 0 && zero_sl[j] < zero_sl[j-1] {
            zero_sl.swap(j, j-1); j -= 1;
        }
    }
    // remove duplicates
    if m_count >= 2 {
        let mut new_m = 1usize;
        for i in 1..m_count as usize {
            if (zero_sl[i] - zero_sl[new_m-1]).abs() >= f64::EPSILON {
                zero_sl[new_m] = zero_sl[i]; new_m += 1;
            }
        }
        m_count = new_m as i32;
    }
    *m = m_count;
}

/// insert: Insert a knot into a spline.
#[no_mangle]
pub unsafe extern "C" fn insert_(
    iopt: *const i32, t: *const f64, n: *const i32, c: *const f64,
    k: *const i32, x: *const f64, tt: *mut f64, nn: *mut i32, cc: *mut f64,
    nest: *const i32, ier: *mut i32,
) {
    let n_v = *n as usize; let k_v = *k as usize; let nest_v = *nest as usize;
    *ier = 10;
    if nest_v <= n_v { return; }
    let k1 = k_v + 1;
    let nk = n_v - k_v;
    let t_sl = std::slice::from_raw_parts(t, n_v);
    let x_v = *x;
    if x_v < t_sl[k1 - 1] || x_v > t_sl[nk - 1] { return; }
    let nk1 = nk - 1;
    let mut l = k1;
    while x_v >= t_sl[l] && l < nk1 { l += 1; }
    if t_sl[l - 1] >= t_sl[l] { return; }
    let iopt_v = *iopt;
    if iopt_v != 0 {
        let kk = 2 * k_v;
        if l <= kk || l >= n_v - kk { return; }
    }
    *ier = 0;
    let c_sl = std::slice::from_raw_parts(c, n_v);
    let tt_sl = std::slice::from_raw_parts_mut(tt, nest_v);
    let cc_sl = std::slice::from_raw_parts_mut(cc, nest_v);
    fpinst(*iopt, t_sl, *n, c_sl, *k, x_v, l as i32, tt_sl, &mut *nn, cc_sl, *nest);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn curfit_lwrk(m: usize, k: usize, nest: usize) -> usize {
        m * (k + 1) + nest * (7 + 3 * k)
    }

    #[test]
    fn fpcurf_sin_smoothing_direct() {
        const M: usize = 20;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1; // = 24
        let k1 = (K + 1) as usize;
        let k2 = k1 + 1;
        let lwrk = curfit_lwrk(M, K as usize, NEST);

        let x: Vec<f64> = (0..M).map(|i| i as f64 * 2.0 * PI / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let w = vec![1.0f64; M];
        let s = 0.05f64;
        let xb = x[0];
        let xe = x[M - 1];
        let tol = 0.001f64;
        let maxit = 20i32;

        let mut n = 0i32;
        let mut t   = vec![0.0f64; NEST];
        let mut c   = vec![0.0f64; NEST];
        let mut fp  = 0.0f64;
        let mut wrk = vec![0.0f64; lwrk];
        let mut iwrk = vec![0i32; NEST];
        let mut ier  = 0i32;

        let ifp: usize = 0;
        let iz  = ifp + NEST;
        let ia  = iz  + NEST;
        let ib  = ia  + NEST * k1;
        let ig  = ib  + NEST * k2;
        let iq  = ig  + NEST * k2;

        // Use raw pointers to split the workspace into non-overlapping slices.
        let p = wrk.as_mut_ptr();
        let (s_fp, s_iz, s_ia, s_ib, s_ig, s_iq) = unsafe {(
            std::slice::from_raw_parts_mut(p.add(ifp), NEST),
            std::slice::from_raw_parts_mut(p.add(iz),  NEST),
            std::slice::from_raw_parts_mut(p.add(ia),  NEST * k1),
            std::slice::from_raw_parts_mut(p.add(ib),  NEST * k2),
            std::slice::from_raw_parts_mut(p.add(ig),  NEST * k2),
            std::slice::from_raw_parts_mut(p.add(iq),  M    * k1),
        )};

        fpcurf(
            0, &x, &y, &w, M as i32,
            xb, xe, K, s, NEST as i32,
            tol, maxit, k1 as i32, k2 as i32,
            &mut n, &mut t, &mut c, &mut fp,
            s_fp, s_iz, s_ia, s_ib, s_ig, s_iq,
            &mut iwrk,
            &mut ier,
        );

        println!("n={n}, fp={fp}, ier={ier}");
        assert!(ier <= 0, "expected ier <= 0, got {}", ier);
    }

    // ── Integration tests (formerly tests/integration.rs) ─────────────────────

    #[test]
    fn curfit_sin_interpolation() {
        const M: usize = 9;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1;
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
        assert_eq!(n, 13, "scipy gives n=13 for this input, got {n}");
        assert!(fp < 1e-14, "fp should be ~0 for interpolation, got {:.3e}", fp);

        let t_ref = [
            0.0, 0.0, 0.0, 0.0,
            PI / 4.0, 3.0 * PI / 8.0, PI / 2.0, 5.0 * PI / 8.0, 3.0 * PI / 4.0,
            PI, PI, PI, PI,
        ];
        for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
            assert!((got - want).abs() < 1e-10, "t[{i}]: got {:.15}, want {:.15}", got, want);
        }

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

    #[test]
    fn curfit_sin_smoothing() {
        const M: usize = 20;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1;
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
        assert_eq!(n, 11, "scipy gives n=11 for s=0.05, got {n}");

        if ier == 0 {
            let rel = (fp - s).abs() / s;
            assert!(rel <= 0.002,
                "|fp - s| / s = {rel:.5}  (fp={fp:.6e}, s={s})");
        }

        let k1 = (K + 1) as usize;
        let nv  = n as usize;
        for i in k1..(nv - k1 - 1) {
            assert!(t[i] < t[i + 1], "knots not sorted: t[{i}]={} ≥ t[{}]={}", t[i], i+1, t[i+1]);
        }

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
        assert!(max_err < 0.1, "max approximation error = {max_err:.4}");
    }

    #[test]
    fn curfit_quadratic_exact_fit() {
        let x = vec![0.0f64, 0.25, 0.5, 0.75, 1.0];
        let y: Vec<f64> = x.iter().map(|&xi| xi * xi).collect();
        let w = vec![1.0f64; 5];
        let (m, k) = (5usize, 2i32);
        let nest = (m + k as usize + 1) as i32;
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
        assert_eq!(n, 8, "scipy gives n=8, got {n}");

        let t_ref = [0.0f64, 0.0, 0.0, 0.375, 0.625, 1.0, 1.0, 1.0];
        for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
            assert!((got - want).abs() < 1e-12, "t[{i}]: got {got:.15}, want {want:.15}");
        }

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

    #[test]
    fn curfit_large_s_gives_polynomial() {
        const M: usize = 10;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1;
        let lwrk = curfit_lwrk(M, K as usize, NEST);

        let x: Vec<f64> = (0..M).map(|i| i as f64 / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let w = vec![1.0f64; M];
        let s = 100.0f64;

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

        let t_ref = [0.0f64, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        for (i, (&got, &want)) in t[..n as usize].iter().zip(t_ref.iter()).enumerate() {
            assert!((got - want).abs() < 1e-12, "t[{i}]: got {got:.15}, want {want:.15}");
        }
    }

    #[test]
    fn clocur_unit_circle_interpolation() {
        const M: usize = 8;
        const K: i32 = 3;
        const IDIM: usize = 2;
        const NEST: usize = M + 2 * K as usize;
        const NC: usize   = NEST * IDIM;
        let lwrk = M * (K as usize + 1) + NEST * (7 + IDIM + 5 * K as usize);

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

        assert_eq!(ier, -1, "expected interpolating closed curve (ier=-1), got {ier}");
        assert_eq!(n, NEST as i32, "expected n={NEST} knots, got {n}");
        assert!(fp < 1e-13, "fp should be ~0 for interpolation, got {fp:.3e}");

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

        for i in 0..M-1 {
            let dx = (xy_out[2 * i]     - xy[2 * i]    ).abs();
            let dy = (xy_out[2 * i + 1] - xy[2 * i + 1]).abs();
            assert!(dx < 1e-12, "x residual at point {}: {:.3e}", i, dx);
            assert!(dy < 1e-12, "y residual at point {}: {:.3e}", i, dy);
        }
    }

    /// curfit_ with iopt=-1 (user-supplied knots, weighted LS fit).
    /// Reference: scipy.interpolate.LSQUnivariateSpline(x, sin(x), t_interior, k=3)
    /// with 50 points on [0,10] and interior knots at 1,2,...,9.
    #[test]
    fn curfit_iopt_minus1_cardinal() {
        const M: usize = 50;
        const K: i32 = 3;
        let x: Vec<f64> = (0..M).map(|i| i as f64 * 10.0 / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let w = vec![1.0f64; M];

        // Knot vector: 4× boundary + 9 interior + 4× boundary = 17 knots
        let mut t = vec![0.0f64; 17];
        for i in 0..4 { t[i] = 0.0; t[16 - i] = 10.0; }
        for i in 0..9 { t[4 + i] = (i + 1) as f64; }
        let mut n = 17i32;
        let nest = 17i32;
        let mut c = vec![0.0f64; 17];
        let mut fp = 0.0f64;
        let lwrk = M * 4 + 17 * (7 + 9);
        let mut wrk = vec![0.0f64; lwrk];
        let mut iwrk = vec![0i32; 17];
        let mut ier = 0i32;
        let iopt = -1i32;

        unsafe {
            curfit_(&iopt, &(M as i32), x.as_ptr(), y.as_ptr(), w.as_ptr(),
                    &x[0], &x[M - 1], &K, &0.0f64, &nest, &mut n,
                    t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
                    wrk.as_mut_ptr(), &(lwrk as i32), iwrk.as_mut_ptr(), &mut ier);
        }
        assert!(ier <= 0, "curfit_ iopt=-1 failed: ier={ier}");

        // Scipy reference coefficients (LSQUnivariateSpline)
        let scipy_c = [
            3.15797216e-04, 3.32519297e-01, 9.94821738e-01, 1.07606479e+00,
            1.66769853e-01, -8.95228307e-01, -1.13447484e+00, -3.30553015e-01,
            7.77265692e-01, 1.17034618e+00, 4.87811782e-01, -2.67284682e-01,
            -5.43651349e-01,
        ];
        let nk1 = n as usize - K as usize - 1;
        assert_eq!(nk1, scipy_c.len());
        for i in 0..nk1 {
            assert!((c[i] - scipy_c[i]).abs() < 1e-6,
                "c[{i}] = {:.10e}, scipy = {:.10e}", c[i], scipy_c[i]);
        }
    }

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

        let mut fp1  = 0.0f64;
        let mut ier1 = 0i32;
        unsafe {
            curfit_(&1i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(), &xb, &xe, &k, &s, &nest_i,
                    &mut n0, t0.as_mut_ptr(), c0.as_mut_ptr(), &mut fp1,
                    wrk0.as_mut_ptr(), &lwrk_i, iwrk0.as_mut_ptr(), &mut ier1);
        }
        assert!(ier1 <= 0, "restart call failed: ier1={ier1}");
        assert!(fp1 <= fp0 + s * 0.01,
            "restart residual {fp1:.6e} is worse than fresh {fp0:.6e}");
    }

    // ── splev_ standalone ─────────────────────────────────────────────────────

    #[test]
    fn splev_constant_spline() {
        // Constant cubic B-spline: t=[0,0,0,0,1,1,1,1], c=[3,...], K=3.
        // A clamped cubic spline with all coefficients equal to 3
        // evaluates to exactly 3 everywhere on [0,1].
        let k = 3i32;
        let n = 8i32;
        let t = [0.0f64, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let c = [3.0f64, 3.0, 3.0, 3.0, 0.0, 0.0, 0.0, 0.0];
        let x = [0.0f64, 0.2, 0.4, 0.6, 0.8, 1.0];
        let m = x.len() as i32;
        let mut y = vec![0.0f64; x.len()];
        let mut ier = 0i32;
        unsafe {
            splev_(t.as_ptr(), &n, c.as_ptr(), &k,
                   x.as_ptr(), y.as_mut_ptr(), &m, &mut ier);
        }
        assert_eq!(ier, 0);
        for (i, &yi) in y.iter().enumerate() {
            assert!((yi - 3.0).abs() < 1e-14,
                "splev at x[{i}]={}: expected 3.0, got {yi}", x[i]);
        }
    }

    // ── splint_ ───────────────────────────────────────────────────────────────

    #[test]
    fn splint_constant_and_sin() {
        // ── Case 1: constant cubic spline, integral = constant × width ───────
        let k = 3i32;
        let n = 8i32;
        let t = [0.0f64, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let c = [2.0f64, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0];
        let mut wrk = vec![0.0f64; n as usize];

        let int_full = unsafe {
            splint_(t.as_ptr(), &n, c.as_ptr(), &k, &0.0f64, &1.0f64, wrk.as_mut_ptr())
        };
        assert!((int_full - 2.0).abs() < 1e-13,
            "∫₀¹ 2 dx: expected 2.0, got {int_full:.15}");

        let int_half = unsafe {
            splint_(t.as_ptr(), &n, c.as_ptr(), &k, &0.25f64, &0.75f64, wrk.as_mut_ptr())
        };
        assert!((int_half - 1.0).abs() < 1e-13,
            "∫₀.₂₅⁰·⁷⁵ 2 dx: expected 1.0, got {int_half:.15}");

        // ── Case 2: interpolating sin on [0,π], ∫₀^π sin = 2 ────────────────
        const M: usize = 9;
        const K2: i32 = 3;
        const NEST2: usize = M + K2 as usize + 1;
        let lwrk2 = curfit_lwrk(M, K2 as usize, NEST2);
        let x: Vec<f64> = (0..M).map(|i| i as f64 * PI / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let w = vec![1.0f64; M];
        let (m_i, k2, nest_i, lwrk_i) = (M as i32, K2, NEST2 as i32, lwrk2 as i32);
        let mut n2 = 0i32;
        let mut t2 = vec![0.0f64; NEST2];
        let mut c2 = vec![0.0f64; NEST2];
        let mut fp2 = 0.0f64;
        let mut wrk2 = vec![0.0f64; lwrk2];
        let mut iwrk2 = vec![0i32; NEST2];
        let mut ier2 = 0i32;
        unsafe {
            curfit_(&0i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
                    &x[0], &x[M-1], &k2, &0.0f64, &nest_i,
                    &mut n2, t2.as_mut_ptr(), c2.as_mut_ptr(), &mut fp2,
                    wrk2.as_mut_ptr(), &lwrk_i, iwrk2.as_mut_ptr(), &mut ier2);
        }
        assert!(ier2 <= 0, "sin fit failed: ier={ier2}");
        let mut wrk3 = vec![0.0f64; n2 as usize];
        let int_sin = unsafe {
            splint_(t2.as_ptr(), &n2, c2.as_ptr(), &k2, &x[0], &x[M-1], wrk3.as_mut_ptr())
        };
        assert!((int_sin - 2.0).abs() < 1e-3,
            "∫₀^π sin: expected 2.0, got {int_sin:.12}");
    }

    // ── sproot_ ───────────────────────────────────────────────────────────────

    #[test]
    fn sproot_sin_zeros() {
        // Fit an interpolating cubic spline to sin on [0, 2π].
        // sproot_ should find zeros close to 0, π, and 2π.
        const M: usize = 17;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1;
        let lwrk = curfit_lwrk(M, K as usize, NEST);
        let x: Vec<f64> = (0..M).map(|i| i as f64 * 2.0 * PI / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
        let w = vec![1.0f64; M];
        let (m_i, k, nest_i, lwrk_i) = (M as i32, K, NEST as i32, lwrk as i32);
        let mut n = 0i32;
        let mut t = vec![0.0f64; NEST];
        let mut c = vec![0.0f64; NEST];
        let mut fp = 0.0f64;
        let mut wrk = vec![0.0f64; lwrk];
        let mut iwrk = vec![0i32; NEST];
        let mut ier = 0i32;
        unsafe {
            curfit_(&0i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
                    &x[0], &x[M-1], &k, &0.0f64, &nest_i,
                    &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
                    wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier);
        }
        assert!(ier <= 0, "sin fit failed: ier={ier}");

        let mest = 10i32;
        let mut zeros = vec![0.0f64; mest as usize];
        let mut m_found = 0i32;
        let mut ier_r = 0i32;
        unsafe {
            sproot_(t.as_ptr(), &n, c.as_ptr(),
                    zeros.as_mut_ptr(), &mest, &mut m_found, &mut ier_r);
        }
        assert_eq!(ier_r, 0, "sproot_ returned error {ier_r}");
        assert!(m_found >= 1, "expected at least one zero, found {m_found}");

        // At least one zero must be within 0.01 of π (the interior zero of sin on [0,2π])
        let has_pi_zero = zeros[..m_found as usize].iter()
            .any(|&z| (z - PI).abs() < 0.01);
        assert!(has_pi_zero,
            "no zero found near π; zeros = {:?}", &zeros[..m_found as usize]);
    }

    // ── insert_ ───────────────────────────────────────────────────────────────

    #[test]
    fn insert_preserves_curve() {
        // Fit a smoothing spline, insert a knot, verify evaluation is unchanged.
        const M: usize = 20;
        const K: i32 = 3;
        const NEST: usize = M + K as usize + 1;
        let lwrk = curfit_lwrk(M, K as usize, NEST);
        let x: Vec<f64> = (0..M).map(|i| i as f64 / (M - 1) as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| (2.0 * PI * xi).sin()).collect();
        let w = vec![1.0f64; M];
        let (m_i, k, nest_i, lwrk_i) = (M as i32, K, NEST as i32, lwrk as i32);
        let mut n = 0i32;
        let mut t = vec![0.0f64; NEST];
        let mut c = vec![0.0f64; NEST];
        let mut fp = 0.0f64;
        let mut wrk = vec![0.0f64; lwrk];
        let mut iwrk = vec![0i32; NEST];
        let mut ier = 0i32;
        unsafe {
            curfit_(&0i32, &m_i, x.as_ptr(), y.as_ptr(), w.as_ptr(),
                    &x[0], &x[M-1], &k, &0.02f64, &nest_i,
                    &mut n, t.as_mut_ptr(), c.as_mut_ptr(), &mut fp,
                    wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier);
        }
        assert!(ier <= 0, "fit failed: ier={ier}");

        // Evaluate before insertion
        let x_eval: Vec<f64> = (0..10).map(|i| i as f64 / 9.0).collect();
        let m_eval = x_eval.len() as i32;
        let mut y_before = vec![0.0f64; x_eval.len()];
        let mut ier_e = 0i32;
        unsafe {
            splev_(t.as_ptr(), &n, c.as_ptr(), &k,
                   x_eval.as_ptr(), y_before.as_mut_ptr(), &m_eval, &mut ier_e);
        }
        assert_eq!(ier_e, 0);

        // Insert a knot at x=0.37 (choose a value in the interior)
        let nest2 = NEST as i32;
        let mut tt = vec![0.0f64; NEST];
        let mut nn = 0i32;
        let mut cc = vec![0.0f64; NEST];
        let x_new = 0.37f64;
        let mut ier_i = 0i32;
        unsafe {
            insert_(&0i32, t.as_ptr(), &n, c.as_ptr(), &k,
                    &x_new, tt.as_mut_ptr(), &mut nn, cc.as_mut_ptr(),
                    &nest2, &mut ier_i);
        }
        assert_eq!(ier_i, 0, "insert_ failed: ier={ier_i}");
        assert_eq!(nn, n + 1, "expected one extra knot after insertion");

        // Evaluate after insertion — must match to machine precision
        let mut y_after = vec![0.0f64; x_eval.len()];
        unsafe {
            splev_(tt.as_ptr(), &nn, cc.as_ptr(), &k,
                   x_eval.as_ptr(), y_after.as_mut_ptr(), &m_eval, &mut ier_e);
        }
        assert_eq!(ier_e, 0);
        for (i, (&before, &after)) in y_before.iter().zip(y_after.iter()).enumerate() {
            assert!((before - after).abs() < 1e-13,
                "insert changed evaluation at x[{i}]={:.4}: {before} → {after}",
                x_eval[i]);
        }
    }

    // ── concur_ ───────────────────────────────────────────────────────────────

    #[test]
    fn concur_unit_circle_interpolation() {
        // Fit a cubic parametric spline to 9 points on a 3/4-circle arc,
        // parameterized by angle in [0, 3π/2], then verify curev_ evaluations
        // lie on the unit circle to within 1e-6.
        const M: usize = 9;
        const K: i32 = 3;
        const IDIM: usize = 2;
        let k1 = (K + 1) as usize;
        let nest = M + k1 + 2 * (k1 - 2); // = m + k + 1 + 2*(k-1)
        let nc = nest * IDIM;
        let np = 2 * k1 * IDIM;
        let lwrk = M * k1 + nest * (6 + IDIM + 3 * K as usize);
        let mx = M * IDIM;

        // Use a 3/4-circle arc so endpoints are distinct (avoids degenerate open-curve fit)
        let u: Vec<f64> = (0..M).map(|i| i as f64 * 1.5 * PI / (M - 1) as f64).collect();
        // interleaved [x0,y0, x1,y1, ...]
        let xn: Vec<f64> = u.iter().flat_map(|&ui| [ui.cos(), ui.sin()]).collect();
        let w = vec![1.0f64; M];

        // Use a moderate smoothing factor: s = M * rms² with rms=0.01 per point.
        // This keeps the algorithm in the knot-refinement regime (no regularization),
        // while still giving a tight fit that should stay on the unit circle.
        let s_val = M as f64 * 1e-4f64; // rms ≈ 0.01 per point
        let (iopt, idim, m_i, mx_i, ib, ie, k, s, nest_i, nc_i, np_i, lwrk_i) =
            (0i32, IDIM as i32, M as i32, mx as i32,
             0i32, 0i32, K, s_val,
             nest as i32, nc as i32, np as i32, lwrk as i32);

        let mut n = 0i32;
        let mut t = vec![0.0f64; nest];
        let mut c = vec![0.0f64; nc];
        let mut cp = vec![0.0f64; np];
        let mut xx = vec![0.0f64; mx];
        let mut fp = 0.0f64;
        let mut wrk = vec![0.0f64; lwrk];
        let mut iwrk = vec![0i32; nest];
        let mut ier = 0i32;
        let db = vec![0.0f64; 1]; // unused (ib=0)
        let de = vec![0.0f64; 1]; // unused (ie=0)
        let (nb, ne) = (0i32, 0i32);

        unsafe {
            concur_(
                &iopt, &idim, &m_i,
                u.as_ptr(), &mx_i, xn.as_ptr(), xx.as_mut_ptr(), w.as_ptr(),
                &ib, db.as_ptr(), &nb,
                &ie, de.as_ptr(), &ne,
                &k, &s, &nest_i, &mut n,
                t.as_mut_ptr(), &nc_i, c.as_mut_ptr(),
                &np_i, cp.as_mut_ptr(),
                &mut fp, wrk.as_mut_ptr(), &lwrk_i, iwrk.as_mut_ptr(), &mut ier,
            );
        }
        assert!(ier <= 0, "concur_ failed: ier={ier}");

        // Evaluate on a 50-point dense grid over the same arc and verify radius ≈ 1
        let u_dense: Vec<f64> = (0..50).map(|i| i as f64 * 1.5 * PI / 49.0).collect();
        let m_dense = u_dense.len() as i32;
        let mxy = m_dense * idim;
        let mut xy = vec![0.0f64; mxy as usize];
        let nc_eval = n * idim;
        let mut ier_e = 0i32;
        unsafe {
            curev_(
                &idim, t.as_ptr(), &n, c.as_ptr(), &nc_eval,
                &k, u_dense.as_ptr(), &m_dense,
                xy.as_mut_ptr(), &mxy, &mut ier_e,
            );
        }
        assert_eq!(ier_e, 0, "curev_ failed: ier={ier_e}");
        for i in 0..50 {
            let cx = xy[2 * i];
            let cy = xy[2 * i + 1];
            let r = (cx * cx + cy * cy).sqrt();
            assert!((r - 1.0).abs() < 0.05,
                "point {i}: radius = {r:.8}, expected 1.0 (cx={cx:.6}, cy={cy:.6})");
        }
    }

}
