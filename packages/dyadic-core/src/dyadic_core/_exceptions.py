__all__ = [
    "ConvergenceError",
    "DomainError",
    "DyadicError",
    "NonInvertibleError",
]


class DyadicError(ValueError):
    """Base for all dyadic-core errors."""


class NonInvertibleError(DyadicError):
    """Raised when an inversion is attempted on a non-invertible element."""


class ConvergenceError(DyadicError):
    """Raised when an iterative method fails to converge."""


class DomainError(DyadicError):
    """Raised when a function argument is outside its domain of convergence."""
