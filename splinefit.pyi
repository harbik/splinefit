"""Type stubs for the splinefit native Python module."""

import numpy as np
import numpy.typing as npt

class CubicSpline:
    """A cubic B-spline curve fitted to 1-D data.

    Construct via class methods: ``smoothing``, ``interpolating``, or ``cardinal``.
    Then evaluate, integrate, or find roots.

    Example
    -------
    >>> import numpy as np
    >>> from splinefit import CubicSpline
    >>> x = np.linspace(0, 2 * np.pi, 50)
    >>> y = np.sin(x)
    >>> s = CubicSpline.smoothing(x, y, rms=0.05)
    >>> y_fit = s.evaluate(x)
    """

    @staticmethod
    def smoothing(
        x: npt.ArrayLike, y: npt.ArrayLike, rms: float
    ) -> CubicSpline:
        """Fit a smoothing cubic spline.

        Parameters
        ----------
        x : array_like
            Strictly increasing abscissae (at least 4 points).
        y : array_like
            Ordinate values, same length as ``x``.
        rms : float
            Target root-mean-square residual. Smaller values produce
            more knots and a closer fit.

        Returns
        -------
        CubicSpline
        """
        ...

    @staticmethod
    def interpolating(x: npt.ArrayLike, y: npt.ArrayLike) -> CubicSpline:
        """Fit an interpolating cubic spline through every data point.

        Parameters
        ----------
        x : array_like
            Strictly increasing abscissae (at least 4 points).
        y : array_like
            Ordinate values, same length as ``x``.

        Returns
        -------
        CubicSpline
        """
        ...

    @staticmethod
    def cardinal(
        x: npt.ArrayLike, y: npt.ArrayLike, dt: float
    ) -> CubicSpline:
        """Fit a cardinal cubic spline on an equidistant knot grid.

        Parameters
        ----------
        x : array_like
            Strictly increasing abscissae (at least 4 points).
        y : array_like
            Ordinate values, same length as ``x``.
        dt : float
            Knot spacing.

        Returns
        -------
        CubicSpline
        """
        ...

    def evaluate(self, x: npt.ArrayLike) -> np.ndarray:
        """Evaluate the spline at the given points.

        Parameters
        ----------
        x : array_like
            Points at which to evaluate, within the spline's domain.

        Returns
        -------
        numpy.ndarray
            Fitted values, same length as ``x``.
        """
        ...

    def integral(self, a: float, b: float) -> float:
        """Compute the definite integral over [a, b].

        Parameters
        ----------
        a : float
            Lower integration bound.
        b : float
            Upper integration bound.

        Returns
        -------
        float
        """
        ...

    def roots(self) -> np.ndarray:
        """Find all interior zeros of the spline.

        Returns
        -------
        numpy.ndarray
            Zeros in ascending order. Boundary zeros may not be included.
        """
        ...

    def knots(self) -> np.ndarray:
        """Return the knot vector.

        Returns
        -------
        numpy.ndarray
        """
        ...

    def coefficients(self) -> np.ndarray:
        """Return the B-spline coefficients.

        Returns
        -------
        numpy.ndarray
        """
        ...

    @property
    def num_knots(self) -> int:
        """Number of knots."""
        ...
