#!/usr/bin/env python3
"""Fit a smoothing spline to noisy sin(x) data and compare with the true function.

Setup (run once):
    python3 -m venv .venv && source .venv/bin/activate
    pip install splinefit numpy
"""

import numpy as np
from splinefit import CubicSpline

# Generate noisy sin(x) data
rng = np.random.default_rng(42)
x = np.linspace(0, 2 * np.pi, 50)
y_true = np.sin(x)
y_noisy = y_true + rng.normal(0, 0.1, len(x))

# Fit with different smoothing levels
for rms in [0.05, 0.1, 0.2]:
    spline = CubicSpline.smoothing(x, y_noisy, rms=rms)
    y_fit = spline.evaluate(x)
    max_err = np.max(np.abs(y_true - y_fit))
    print(f"rms={rms:.2f}  knots={spline.num_knots:3d}  max error vs true: {max_err:.4f}")

# Interpolating spline (passes through every noisy point)
spline_interp = CubicSpline.interpolating(x, y_noisy)
y_interp = spline_interp.evaluate(x)
residual = np.max(np.abs(y_noisy - y_interp))
print(f"\nInterpolating: knots={spline_interp.num_knots}, max residual: {residual:.2e}")
