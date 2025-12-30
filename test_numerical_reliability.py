"""
Tests for numerical reliability framework.
"""

from sage.numerical.root_reliability import residual_check
from sage.all import sqrt


def test_residual_check_reliable():
    """
    Test that a correct numerical root is marked reliable.
    """
    f = lambda x: x^2 - 2
    r = sqrt(2).n()
    result = residual_check(f, r)
    assert result.reliable is True
    assert result.residual < 1e-10


def test_residual_check_unreliable():
    """
    Test that an incorrect root is marked unreliable.
    """
    f = lambda x: x^2 - 2
    r = 1.0
    result = residual_check(f, r)
    assert result.reliable is False

