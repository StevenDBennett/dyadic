"""
core.py (backward-compatibility re-exports)

All code has been moved to focused submodules.  This file re-exports
everything for code that does ``from dyadic_core.core import ...``.
"""

from dyadic_core._dlog import (  # noqa: F401
    DLogNewtonStep,
    dlog_bootstrap,
    dlog_residual_tracking,
    two_adic_dlog,
    two_adic_log5,
)
from dyadic_core._dual import (  # noqa: F401
    DualNumber,
    TwoAdicProcessor,
    run_all_tests,
)
from dyadic_core._exceptions import (  # noqa: F401
    ConvergenceError,
    DomainError,
    DyadicError,
    NonInvertibleError,
)
from dyadic_core._modular import (  # noqa: F401
    dual_add,
    modinv_newton,
)
from dyadic_core._series import (  # noqa: F401
    g0,
    padic_exp,
    padic_log,
)
from dyadic_core._util import (  # noqa: F401
    bitmask,
    mat_det,
    mat_mul,
    valuation,
)
