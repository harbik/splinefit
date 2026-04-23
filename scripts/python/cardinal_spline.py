#!/usr/bin/env python3
"""Compare cardinal splines at different knot spacings on data with a sharp feature.

Setup (run once):
    python3 -m venv .venv && source .venv/bin/activate
    pip install splinefit numpy
"""

import numpy as np
from splinefit import CubicSpline

# A function with a sharp feature: Gaussian bump on a sine wave
x = np.linspace(0, 10, 200)
y = np.sin(x) + 2 * np.exp(-((x - 5) ** 2) / 0.2)

print("Cardinal spline fits with varying knot spacing:")
print(f"{'dt':>8}  {'knots':>6}  {'max error':>10}  {'rms error':>10}")
print("-" * 42)

for dt in [2.0, 1.0, 0.5, 0.25, 0.1]:
    spline = CubicSpline.cardinal(x, y, dt=dt)
    y_fit = spline.evaluate(x)
    errors = np.abs(y - y_fit)
    print(f"{dt:8.2f}  {spline.num_knots:6d}  {errors.max():10.6f}  {np.sqrt(np.mean(errors**2)):10.6f}")

# Interpolating spline captures the peak exactly
spline_exact = CubicSpline.interpolating(x, y)
y_fit = spline_exact.evaluate(x)
peak_idx = np.argmax(y)
print(f"\nInterpolating ({spline_exact.num_knots} knots):")
print(f"  Peak at x={x[peak_idx]:.1f}: true={y[peak_idx]:.6f}, fit={y_fit[peak_idx]:.6f}")
print(f"  Max residual: {np.max(np.abs(y - y_fit)):.2e}")
