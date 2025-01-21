"""
Witt rings: general constructor
"""
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.rings.padics.witt_ring import (
    WittRing_finite_field,
    WittRing_non_p_typical,
    WittRing_p_invertible,
    WittRing_p_typical,
)
from sage.sets.primes import Primes

_CommutativeRings = CommutativeRings()
_Primes = Primes()


def WittRing(base_ring, prec=1, p=None, algorithm='auto'):
    """
    Return the appropriate Witt ring, depending on the input.

    EXAMPLES::

        sage: WittRing(QQ,p=5)
        Ring of 5-Witt Vectors of length 1 over Rational Field
        sage: WittRing(GF(3))
        Ring of Witt Vectors of length 1 over Finite Field of size 3
        sage: WittRing(GF(3)['t'])
        Ring of Witt Vectors of length 1 over Univariate Polynomial Ring in t
        over Finite Field of size 3
        sage: WittRing(Qp(7), prec=30, p=5)
        Ring of 5-Witt Vectors of length 30 over 7-adic Field
        with capped relative precision 20

    TESTS::

        sage: A = SymmetricGroup(3).algebra(QQ)
        sage: WittRing(A)
        Traceback (most recent call last):
        ...
        TypeError: Symmetric group algebra of order 3 over Rational Field
        is not a commutative ring
        sage: WittRing(QQ)
        Traceback (most recent call last):
        ...
        ValueError: Rational Field has non-prime characteristic
        and no prime was supplied

        sage: WittRing(QQ, p=5, algorithm='moon')
        Traceback (most recent call last):
        ...
        ValueError: algorithm must be one of 'none', 'auto',
        'standard', 'finotti'
    """
    if base_ring not in _CommutativeRings:
        raise TypeError(f'{base_ring} is not a commutative ring')

    char = base_ring.characteristic()
    if p is None:
        if char not in _Primes:
            raise ValueError(f'{base_ring} has non-prime characteristic '
                             'and no prime was supplied')
        else:
            prime = char
    else:
        prime = p

    if algorithm is None:
        algorithm = 'none'
    elif algorithm not in ['none', 'auto', 'standard', 'finotti']:
        raise ValueError("algorithm must be one of 'none', 'auto', "
                         "'standard', 'finotti'")

    if prime == char:  # p-typical
        if base_ring in Fields().Finite():
            # TODO: document that this ignores the choice of algorithm
            return WittRing_finite_field(base_ring.field(), prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'finotti'
            return WittRing_p_typical(base_ring, prec, prime,
                                      algorithm=algorithm,
                                      category=_CommutativeRings)
    else:  # non-p-typical
        if algorithm == 'finotti':
            raise ValueError("The 'finotti' algorithm only works "
                             "for p-typical Witt Rings.")
        if base_ring(prime).is_unit():
            # TODO: document that this ignores the choice of algorithm
            return WittRing_p_invertible(base_ring, prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'standard'
            return WittRing_non_p_typical(base_ring, prec, prime,
                                          algorithm=algorithm,
                                          category=_CommutativeRings)
