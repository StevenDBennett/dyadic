"""
popcount_compression.py
-----------------------
Exploratory: correlation between popcount gap and delta popcount.

Moved from dyadic_math.padic_roots — no theoretical grounding,
no test coverage.  Kept for reference.
"""

import random

import numpy as np


def popcount_compression(k: int, n_trials: int = 100) -> dict[str, float]:
    """Measure correlation between popcount gap and sum of delta popcounts."""
    popcounts_num = []
    popcounts_delta = []

    for _ in range(n_trials):
        e_boot = random.randrange(0, 1 << (k - 2))
        e_true = random.randrange(0, 1 << (k - 2))
        e0 = (e_boot - e_true) & ((1 << (k - 2)) - 1)

        pc0 = e_boot.bit_count()
        pc_true = e_true.bit_count()
        pc_delta = e0.bit_count()

        popcounts_num.append(pc0 - pc_true)
        popcounts_delta.append(pc_delta)

    r = np.corrcoef(popcounts_num, popcounts_delta)[0, 1]
    return {
        "pearson_r": float(r),
        "n_trials": n_trials,
    }
