#!/usr/bin/env python3
"""Demonstrate integration and root-finding on a spline."""

import numpy as np
from splinefit import CubicSpline

# Fit an interpolating spline to sin(x) on [0, 2pi]
x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)
spline = CubicSpline.interpolating(x, y)

# Integration
integral_0_pi = spline.integral(0, np.pi)       # should be ~2.0
integral_full = spline.integral(0, 2 * np.pi)   # should be ~0.0
print(f"integral [0, pi]:   {integral_0_pi:.10f}  (exact: 2.0)")
print(f"integral [0, 2pi]:  {integral_full:.10f}  (exact: 0.0)")

# Root-finding
zeros = spline.roots()
print(f"\nZeros found: {len(zeros)}")
for z in zeros:
    print(f"  x = {z:.8f}  (sin(x) = {np.sin(z):.2e})")

# Compare with known zeros: pi
pi_zero = zeros[np.argmin(np.abs(zeros - np.pi))]
print(f"\nClosest zero to pi: {pi_zero:.10f}  (error: {abs(pi_zero - np.pi):.2e})")
