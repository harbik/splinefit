// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2021-2026, Harbers Bik LLC

//! Benchmarks for spline fitting and evaluation.
//!
//! Run with:  cargo bench
//!
//! Compare against the SciPy/Fortran reference by running:
//!   python3 benches/scipy_bench.py

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use splinefit::{CubicSplineFit, CubicSplineFit2D, evaluate};
use std::f64::consts::PI;

/// Generate `m` evenly-spaced points of sin(x) on [0, 2π].
fn sin_data(m: usize) -> (Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = (0..m).map(|i| i as f64 * 2.0 * PI / (m - 1) as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| xi.sin()).collect();
    (x, y)
}

/// Generate `m` evenly-spaced points on a unit circle arc [0, 1.5π].
fn circle_data(m: usize) -> (Vec<f64>, Vec<f64>) {
    let u: Vec<f64> = (0..m).map(|i| i as f64 * 1.5 * PI / (m - 1) as f64).collect();
    let xn: Vec<f64> = u.iter().flat_map(|&ui| [ui.cos(), ui.sin()]).collect();
    (u, xn)
}

fn bench_1d_fit(c: &mut Criterion) {
    let mut group = c.benchmark_group("1d_fit");
    for &m in &[50, 200, 1000, 5000] {
        let (x, y) = sin_data(m);
        group.bench_with_input(BenchmarkId::new("smoothing", m), &m, |b, _| {
            b.iter(|| {
                CubicSplineFit::new(x.clone(), y.clone())
                    .smoothing_spline(black_box(0.05))
                    .unwrap()
            });
        });
        group.bench_with_input(BenchmarkId::new("interpolating", m), &m, |b, _| {
            b.iter(|| {
                CubicSplineFit::new(x.clone(), y.clone())
                    .interpolating_spline()
                    .unwrap()
            });
        });
    }
    group.finish();
}

fn bench_1d_evaluate(c: &mut Criterion) {
    let mut group = c.benchmark_group("1d_evaluate");
    for &m in &[50, 200, 1000, 5000] {
        let (x, y) = sin_data(m);
        let spline = CubicSplineFit::new(x.clone(), y)
            .smoothing_spline(0.05)
            .unwrap();
        // Evaluate at 10× the data points
        let x_eval: Vec<f64> = (0..m * 10)
            .map(|i| i as f64 * 2.0 * PI / (m * 10 - 1) as f64)
            .collect();
        group.bench_with_input(BenchmarkId::new("splev", m), &m, |b, _| {
            b.iter(|| evaluate::evaluate(black_box(&spline), black_box(&x_eval)).unwrap());
        });
    }
    group.finish();
}

fn bench_2d_fit(c: &mut Criterion) {
    let mut group = c.benchmark_group("2d_fit");
    for &m in &[50, 200, 1000] {
        let (u, xn) = circle_data(m);
        group.bench_with_input(BenchmarkId::new("smoothing", m), &m, |b, _| {
            b.iter(|| {
                CubicSplineFit2D::new(u.clone(), xn.clone())
                    .unwrap()
                    .smoothing_spline(black_box(0.01))
                    .unwrap()
            });
        });
    }
    group.finish();
}

fn bench_2d_evaluate(c: &mut Criterion) {
    let mut group = c.benchmark_group("2d_evaluate");
    for &m in &[50, 200, 1000] {
        let (u, xn) = circle_data(m);
        let spline = CubicSplineFit2D::new(u.clone(), xn)
            .unwrap()
            .smoothing_spline(0.01)
            .unwrap();
        let u_eval: Vec<f64> = (0..m * 10)
            .map(|i| i as f64 * 1.5 * PI / (m * 10 - 1) as f64)
            .collect();
        group.bench_with_input(BenchmarkId::new("curev", m), &m, |b, _| {
            b.iter(|| evaluate::evaluate(black_box(&spline), black_box(&u_eval)).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench_1d_fit, bench_1d_evaluate, bench_2d_fit, bench_2d_evaluate);
criterion_main!(benches);
