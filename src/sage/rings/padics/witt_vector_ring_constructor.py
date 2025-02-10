r"""
Constructor for Witt vector rings

This module provides the function :func:`WittVectorRing`, which constructs
rings of Witt vectors with coefficients in any commuatitve ring.
"""
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.rings.padics.witt_vector_ring import (
    WittVectorRing_finite_field,
    WittVectorRing_non_p_typical,
    WittVectorRing_p_invertible,
    WittVectorRing_p_typical,
)
from sage.sets.primes import Primes

_CommutativeRings = CommutativeRings()
_Primes = Primes()


def WittVectorRing(base_ring, prec=1, p=None, algorithm='auto'):
    """
    Return the appropriate Witt vector ring, depending on the input.

    EXAMPLES::

        sage: WittVectorRing(QQ,p=5)
        Ring of 5-Witt Vectors of length 1 over Rational Field
        sage: WittVectorRing(GF(3))
        Ring of Witt Vectors of length 1 over Finite Field of size 3
        sage: WittVectorRing(GF(3)['t'])
        Ring of Witt Vectors of length 1 over Univariate Polynomial Ring in t
        over Finite Field of size 3
        sage: WittVectorRing(Qp(7), prec=30, p=5)
        Ring of 5-Witt Vectors of length 30 over 7-adic Field
        with capped relative precision 20

    TESTS::

        sage: A = SymmetricGroup(3).algebra(QQ)
        sage: WittVectorRing(A)
        Traceback (most recent call last):
        ...
        TypeError: Symmetric group algebra of order 3 over Rational Field
        is not a commutative ring
        sage: WittVectorRing(QQ)
        Traceback (most recent call last):
        ...
        ValueError: Rational Field has non-prime characteristic
        and no prime was supplied

        sage: WittVectorRing(QQ, p=5, algorithm='moon')
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
            return WittVectorRing_finite_field(base_ring.field(), prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'finotti'
            return WittVectorRing_p_typical(base_ring, prec, prime,
                                      algorithm=algorithm,
                                      category=_CommutativeRings)
    else:  # non-p-typical
        if algorithm == 'finotti':
            raise ValueError("The 'finotti' algorithm only works "
                             "for p-typical Witt vector rings.")
        if base_ring(prime).is_unit():
            # TODO: document that this ignores the choice of algorithm
            return WittVectorRing_p_invertible(base_ring, prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'standard'
            return WittVectorRing_non_p_typical(base_ring, prec, prime,
                                          algorithm=algorithm,
                                          category=_CommutativeRings)
