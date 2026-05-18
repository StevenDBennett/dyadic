"""Property-based tests for dyadic_core invariants."""
from dyadic_core import (
    DualNumber,
    TwoAdicProcessor,
    bitmask,
    modinv_newton,
    two_adic_dlog,
    valuation,
)
from hypothesis import given
from hypothesis.strategies import integers


@given(integers(min_value=1, max_value=20))
def test_bitmask_positive(k):
    assert bitmask(k) == (1 << k) - 1


@given(integers(min_value=0, max_value=30))
def test_valuation_powers_of_two(n):
    assert valuation(1 << n) == n


@given(integers(min_value=1, max_value=1000).filter(lambda x: x % 2 == 1))
def test_valuation_odd(n):
    assert valuation(n) == 0


@given(
    integers(min_value=3, max_value=16),
    integers(min_value=1, max_value=1000).filter(lambda x: x % 2 == 1),
)
def test_modinv_newton_correctness(k, a):
    a = a % (1 << k)
    if a % 2 == 0:
        return
    inv = modinv_newton(a, k)
    assert (a * inv) & bitmask(k) == 1


@given(integers(min_value=3, max_value=14))
def test_dual_number_roundtrip_from_int(k):
    mask = bitmask(k)
    for n in range(0, min(1 << k, 256)):
        d = DualNumber(n, k)
        assert d.value == n & mask
        if not d.is_zero:
            assert d.verify()


@given(
    integers(min_value=3, max_value=12),
    integers(min_value=0, max_value=10),
    integers(min_value=0, max_value=1),
    integers(min_value=0, max_value=100),
)
def test_dual_number_from_coords_roundtrip(k, v, alpha, e):
    if v >= k:
        return
    d = DualNumber.from_coords(v, alpha, e, k)
    assert d.verify()
    assert d.v == v
    assert d.alpha == alpha & 1


@given(
    integers(min_value=4, max_value=14),
    integers(min_value=1, max_value=100).filter(lambda x: x % 2 == 1),
)
def test_dlog_roundtrip(k, e_true):
    a = pow(5, e_true, 1 << k)
    if a & 3 != 1:
        return
    result = two_adic_dlog(a, k)
    assert result is not None
    alpha, e = result
    recomputed = pow(5, e, 1 << k)
    if alpha:
        recomputed = (-recomputed) & bitmask(k)
    assert recomputed == a & bitmask(k)


@given(
    integers(min_value=4, max_value=12),
    integers(min_value=1, max_value=30),
)
def test_exponential_isometry(k, e):
    order = 1 << (k - 2)
    e = e % order
    if e == 0:
        return
    diff = (pow(5, e, 1 << k) - 1) & bitmask(k)
    v2_diff = valuation(diff) if diff != 0 else k
    v2_e = valuation(e)
    assert v2_diff == v2_e + 2


@given(
    integers(min_value=4, max_value=12),
    integers(min_value=1, max_value=100).filter(lambda x: x % 2 == 1),
    integers(min_value=1, max_value=100).filter(lambda x: x % 2 == 1),
)
def test_processor_mul_commutative(k, a, b):
    proc = TwoAdicProcessor(k)
    da = DualNumber(a, k)
    db = DualNumber(b, k)
    ab = proc.mul(da, db)
    ba = proc.mul(db, da)
    assert ab.value == ba.value


@given(
    integers(min_value=4, max_value=12),
    integers(min_value=1, max_value=100).filter(lambda x: x % 2 == 1),
)
def test_processor_inv_identity(k, a):
    proc = TwoAdicProcessor(k)
    da = DualNumber(a, k)
    inv = proc.inv(da)
    prod = proc.mul(da, inv)
    assert prod.value == 1


@given(
    integers(min_value=4, max_value=10),
    integers(min_value=1, max_value=100).filter(lambda x: x % 2 == 1),
    integers(min_value=0, max_value=10),
)
def test_processor_pow_correctness(k, a, n):
    proc = TwoAdicProcessor(k)
    da = DualNumber(a, k)
    result = proc.pow(da, n)
    expected = pow(a, n, 1 << k)
    assert result.value == expected
