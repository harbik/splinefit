#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0 OR MIT
# Copyright (c) 2021-2026, Harbers Bik LLC

"""
SciPy/Fortran FITPACK benchmark — counterpart to benches/spline_bench.rs.

Run:  python3 benches/scipy_bench.py

Requires: pip install scipy numpy
"""

import time
import numpy as np
from scipy.interpolate import splrep, splev, splprep

REPEATS = 200  # number of iterations per timing


def time_it(fn, repeats=REPEATS):
    """Run fn() `repeats` times and return (mean_seconds, result_of_last_call)."""
    # warm-up
    result = fn()
    t0 = time.perf_counter()
    for _ in range(repeats):
        result = fn()
    elapsed = time.perf_counter() - t0
    return elapsed / repeats, result


def bench_1d(m):
    x = np.linspace(0, 2 * np.pi, m)
    y = np.sin(x)

    # --- smoothing fit (s = 0.05 * m, since SciPy's s is sum-of-squares) ---
    s_val = 0.05**2 * m
    t_fit, tck = time_it(lambda: splrep(x, y, s=s_val))

    # --- interpolating fit (s = 0) ---
    t_interp, _ = time_it(lambda: splrep(x, y, s=0))

    # --- evaluation at 10x points ---
    x_eval = np.linspace(0, 2 * np.pi, m * 10)
    t_eval, _ = time_it(lambda: splev(x_eval, tck))

    return t_fit, t_interp, t_eval


def bench_2d(m):
    u = np.linspace(0, 1.5 * np.pi, m)
    x = np.cos(u)
    y = np.sin(u)

    # --- smoothing fit ---
    s_val = 0.01**2 * m
    t_fit, (tck, _u) = time_it(lambda: splprep([x, y], u=u, s=s_val))

    # --- evaluation at 10x points ---
    u_eval = np.linspace(0, 1.5 * np.pi, m * 10)
    t_eval, _ = time_it(lambda: splev(u_eval, tck))

    return t_fit, t_eval


def fmt(seconds):
    if seconds < 1e-6:
        return f"{seconds * 1e9:8.1f} ns"
    elif seconds < 1e-3:
        return f"{seconds * 1e6:8.1f} µs"
    else:
        return f"{seconds * 1e3:8.1f} ms"


def main():
    print(f"SciPy FITPACK benchmark ({REPEATS} iterations each)")
    print("=" * 65)

    print("\n1-D spline (sin on [0, 2π])")
    print(f"{'m':>6}  {'smooth fit':>12}  {'interp fit':>12}  {'eval 10×m':>12}")
    print("-" * 50)
    for m in [50, 200, 1000, 5000]:
        t_fit, t_interp, t_eval = bench_1d(m)
        print(f"{m:6d}  {fmt(t_fit)}  {fmt(t_interp)}  {fmt(t_eval)}")

    print("\n2-D parametric spline (circle arc)")
    print(f"{'m':>6}  {'smooth fit':>12}  {'eval 10×m':>12}")
    print("-" * 36)
    for m in [50, 200, 1000]:
        t_fit, t_eval = bench_2d(m)
        print(f"{m:6d}  {fmt(t_fit)}  {fmt(t_eval)}")


if __name__ == "__main__":
    main()
