r"""
Constructor for Witt vector rings

This module provides the function :func:`WittVectorRing`, which constructs
rings of Witt vectors with coefficients in any commuatitve ring.
"""
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.rings.padics.witt_vector_ring import (
    WittVectorRing_base,
    WittVectorRing_char_p,
    WittVectorRing_finite_field,
)
from sage.sets.primes import Primes

_CommutativeRings = CommutativeRings()
_Primes = Primes()


def WittVectorRing(base_ring, prec=1, p=None, algorithm=None):
    """
    Return the appropriate Witt vector ring, depending on the input.

    EXAMPLES::

        sage: WittVectorRing(QQ,p=5)
        Ring of truncated 5-typical Witt vectors of length 1 over Rational Field
        sage: WittVectorRing(GF(3))
        Ring of truncated 3-typical Witt vectors of length 1 over Finite Field of size 3
        sage: WittVectorRing(GF(3)['t'])
        Ring of truncated 3-typical Witt vectors of length 1 over Univariate Polynomial Ring in t
        over Finite Field of size 3
        sage: WittVectorRing(Qp(7), prec=30, p=5)
        Ring of truncated 5-typical Witt vectors of length 30 over 7-adic Field
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
        ValueError: algorithm must be one of None, 'standard',
        'p_invertible', 'finotti', 'Zq_isomorphism'
    """
    if base_ring not in _CommutativeRings:
        raise TypeError(f'{base_ring} is not a commutative ring')

    if not is_Integer (prec):
        raise TypeError(f'{prec} is not an integer')
    elif prec <= 0:
        raise ValueError(f'{prec} must be positive')

    char = base_ring.characteristic()
    if p is None:
        if char not in _Primes:
            raise ValueError(f'{base_ring} has non-prime characteristic '
                             'and no prime was supplied')
        else:
            prime = char
    else:
        prime = p

    if algorithm not in [None, 'standard', 'p_invertible', 'finotti',
                         'Zq_isomorphism']:
        raise ValueError(
            "algorithm must be one of None, 'standard', "
            "'p_invertible', 'finotti', 'Zq_isomorphism'")

    if prime == char:
        if algorithm == 'p_invertible':
            raise ValueError (
                "The 'p_invertible' algorithm only works when p is a unit "
                "in the ring of coefficients.")
        elif base_ring in Fields().Finite():
            if not algorithm:
                algorithm = 'Zq_isomorphism'
            return WittVectorRing_finite_field(
                base_ring.field(), prec, prime, algorithm=algorithm,
                category=_CommutativeRings)
        else:
            if not algorithm:
                algorithm = 'finotti'
            elif algorithm == 'Zq_isomorphism':
                raise ValueError(
                    "The 'Zq_isomorphism' algorithm only works when the "
                    "coefficient ring is a finite field of characteristic p.")
            return WittVectorRing_char_p(base_ring, prec, prime,
                                      algorithm=algorithm,
                                      category=_CommutativeRings)
    else:
        if algorithm == 'finotti':
            raise ValueError("The 'finotti' algorithm only works for "
                             "coefficients rings of characteristic p.")
        elif algorithm == 'Zq_isomorphism':
            raise ValueError(
                "The 'Zq_isomorphism' algorithm only works when the "
                "coefficient ring is a finite field of characteristic p.")
        if base_ring(prime).is_unit():
            if not algorithm:
                algorithm = 'p_invertible'
            return WittVectorRing_base(
                base_ring, prec, prime, algorithm='p_invertible',
                category=_CommutativeRings)
        else:
            if not algorithm:
                algorithm = 'standard'
            elif algorithm == 'p_invertible':
                raise ValueError (
                    "The 'p_invertible' algorithm only works when p is a "
                    "unit in the ring of coefficients.")
            return WittVectorRing_base(base_ring, prec, prime,
                                          algorithm=algorithm,
                                          category=_CommutativeRings)
