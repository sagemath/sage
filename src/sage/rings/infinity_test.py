import pytest

from sage.rings.infinity import Infinity, InfinityRing, infinity, minus_infinity
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.qqbar import AA
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.real_lazy import RLF
from sage.rings.real_mpfi import RIF
from sage.rings.real_mpfr import RR, RealField
from sage.symbolic.ring import SR


def generate_elements(ring):
    elements = [-1e3, 99.9999, 0, 1, 100000]
    elements += [-SR(2).sqrt(), SR.pi(), 3 ** (-QQ.one() / 3)]
    elements.append(ring.an_element())
    elements.extend(ring.some_elements())

    valid_elements = []
    for element in elements:
        try:
            valid_elements.append(ring(element))
        except (ValueError, TypeError):
            continue
    return valid_elements


@pytest.mark.parametrize(
    "element",
    [
        element
        for ring in [ZZ, QQ, RDF, RR, RealField(200), RLF, RIF, AA]
        for element in generate_elements(ring)
    ],
)
def test_comparison(element):
    """
    Check comparison with infinity.

    Various attempts are made to generate elements of ``ring``. An
    assertion is triggered if one of these elements does not compare
    correctly with plus/minus infinity.
    """
    assert minus_infinity < element
    assert element > minus_infinity
    assert element < Infinity
    assert infinity > element
    assert minus_infinity <= element
    assert element >= minus_infinity
    assert element <= infinity
    assert infinity >= element


def test_comparison_number_field():
    """
    Comparison with number fields does not work.
    This is a known bug.
    """
    x = polygen(ZZ, "x")
    sqrt3 = NumberField(x**2 - 3, "sqrt3").gen()
    element = 1 + sqrt3
    with pytest.raises(TypeError):
        assert not (minus_infinity < element < infinity)


def test_comparison_symbolic_ring():
    """
    The symbolic ring handles its own infinities, but answers
    ``False`` (meaning: cannot decide) already for some very
    elementary comparisons. This is a known bug.
    """
    element = SR.an_element()
    assert not (minus_infinity < element)
    assert not (element > minus_infinity)
    assert not (element < infinity)
    assert not (infinity > element)
    assert not (minus_infinity <= element)
    assert not (element >= minus_infinity)
    assert not (element <= infinity)
    assert not (infinity >= element)


@pytest.mark.parametrize(
    "pos_inf", [infinity, float("+inf"), RLF(infinity), RIF(infinity), SR(infinity)]
)
def test_signed_infinity(pos_inf):
    """
    Test consistency of infinity representations.

    There are different possible representations of infinity in
    Sage. These are all consistent with the infinity ring, that is,
    compare with infinity in the expected way. See also :issue:`14045`
    """
    assert InfinityRing(pos_inf) is infinity
    assert InfinityRing(-pos_inf) is minus_infinity
    assert infinity == pos_inf
    assert not (infinity > pos_inf)
    assert not (infinity < pos_inf)
    assert minus_infinity == -pos_inf
    assert not (minus_infinity > -pos_inf)
    assert not (minus_infinity < -pos_inf)
    assert pos_inf > -pos_inf
    assert infinity > -pos_inf
    assert pos_inf > minus_infinity


def test_issue_14045():
    assert InfinityRing(float("+inf")) is infinity
    assert InfinityRing(float("-inf")) is minus_infinity
    assert not (infinity > float("+inf"))
    assert infinity == float("+inf")
