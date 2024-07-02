# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Linear code constructors that do not preserve the structural information

This file contains a variety of constructions which builds the generator matrix
of special (or random) linear codes and wraps them in a
:class:`sage.coding.linear_code.LinearCode` object. These constructions are
therefore not rich objects such as
:class:`sage.coding.grs_code.GeneralizedReedSolomonCode`.

All codes available here can be accessed through the ``codes`` object::

    sage: codes.random_linear_code(GF(2), 5, 2)
    [5, 2]  linear code over GF(2)

REFERENCES:

- [HP2003]_

AUTHORS:

- David Joyner (2007-05): initial version

- David Joyner (2008-02): added cyclic codes, Hamming codes

- David Joyner (2008-03): added BCH code, LinearCodeFromCheckmatrix, ReedSolomonCode, WalshCode,
  DuadicCodeEvenPair, DuadicCodeOddPair, QR codes (even and odd)

- David Joyner (2008-09) fix for bug in BCHCode reported by F. Voloch

- David Joyner (2008-10) small docstring changes to WalshCode and walsh_matrix
"""
# ****************************************************************************
#       Copyright (C) 2007 David Joyner <wdjoyner@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import gcd, quadratic_residues
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import random_matrix
from sage.misc.misc_c import prod
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sequence import Sequence, Sequence_generic

from .linear_code import LinearCode

############### utility functions ################


def _is_a_splitting(S1, S2, n, return_automorphism=False):
    r"""
    Check whether ``(S1,S2)`` is a splitting of `\ZZ/n\ZZ`.

    A splitting of `R = \ZZ/n\ZZ` is a pair of subsets of `R` which is a
    partition of `R \\backslash \{0\}` and such that there exists an element `r`
    of `R` such that `r S_1 = S_2` and `r S_2 = S_1` (where `r S` is the
    point-wise multiplication of the elements of `S` by `r`).

    Splittings are useful for computing idempotents in the quotient
    ring `Q = GF(q)[x]/(x^n-1)`.

    INPUT:

    - ``S1``, ``S2`` -- disjoint sublists partitioning ``[1, 2, ..., n-1]``

    - ``n`` -- integer

    - ``return_automorphism`` -- boolean (default: ``False``); whether to
      return the automorphism exchanging `S_1` and `S_2`

    OUTPUT:

    If ``return_automorphism is False`` (default) the function returns boolean values.

    Otherwise, it returns a pair ``(b, r)`` where ``b`` is a boolean indicating
    whether `S1`, `S2` is a splitting of `n`, and `r` is such that `r S_1 = S_2`
    and `r S_2 = S_1` (if `b` is ``False``, `r` is equal to ``None``).

    EXAMPLES::

        sage: from sage.coding.code_constructions import _is_a_splitting
        sage: _is_a_splitting([1,2],[3,4],5)
        True
        sage: _is_a_splitting([1,2],[3,4],5,return_automorphism=True)
        (True, 4)

        sage: _is_a_splitting([1,3],[2,4,5,6],7)
        False
        sage: _is_a_splitting([1,3,4],[2,5,6],7)
        False

        sage: for P in SetPartitions(6,[3,3]):
        ....:     res,aut= _is_a_splitting(P[0],P[1],7,return_automorphism=True)
        ....:     if res:
        ....:         print((aut, P))
        (3, {{1, 2, 4}, {3, 5, 6}})
        (6, {{1, 2, 3}, {4, 5, 6}})
        (6, {{1, 3, 5}, {2, 4, 6}})
        (6, {{1, 4, 5}, {2, 3, 6}})

    We illustrate now how to find idempotents in quotient rings::

        sage: n = 11; q = 3
        sage: C = Zmod(n).cyclotomic_cosets(q); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: _is_a_splitting(S1,S2,11)
        True
        sage: F = GF(q)
        sage: P.<x> = PolynomialRing(F,"x")
        sage: I = Ideal(P,[x^n-1])
        sage: Q.<x> = QuotientRing(P,I)
        sage: i1 = -sum([x^i for i in S1]); i1
        2*x^9 + 2*x^5 + 2*x^4 + 2*x^3 + 2*x
        sage: i2 = -sum([x^i for i in S2]); i2
        2*x^10 + 2*x^8 + 2*x^7 + 2*x^6 + 2*x^2
        sage: i1^2 == i1
        True
        sage: i2^2 == i2
        True
        sage: (1-i1)^2 == 1-i1
        True
        sage: (1-i2)^2 == 1-i2
        True

    We return to dealing with polynomials (rather than elements of
    quotient rings), so we can construct cyclic codes::

        sage: P.<x> = PolynomialRing(F,"x")
        sage: i1 = -sum([x^i for i in S1])
        sage: i2 = -sum([x^i for i in S2])
        sage: i1_sqrd = (i1^2).quo_rem(x^n-1)[1]
        sage: i1_sqrd  == i1
        True
        sage: i2_sqrd = (i2^2).quo_rem(x^n-1)[1]
        sage: i2_sqrd  == i2
        True
        sage: C1 = codes.CyclicCode(length = n, generator_pol = gcd(i1, x^n - 1))
        sage: C2 = codes.CyclicCode(length = n, generator_pol = gcd(1-i2, x^n - 1))
        sage: C1.dual_code().systematic_generator_matrix() == C2.systematic_generator_matrix()
        True

    This is a special case of Theorem 6.4.3 in [HP2003]_.
    """
    R = IntegerModRing(n)
    S1 = {R(x) for x in S1}
    S2 = {R(x) for x in S2}

    # we first check whether (S1,S2) is a partition of R - {0}
    if (len(S1) + len(S2) != n-1 or len(S1) != len(S2) or
        R.zero() in S1 or R.zero() in S2 or not S1.isdisjoint(S2)):
        if return_automorphism:
            return False, None
        else:
            return False

    # now that we know that (S1,S2) is a partition, we look for an invertible
    # element b that maps S1 to S2 by multiplication
    for b in Integer(n).coprime_integers(n):
        if b >= 2 and all(b * x in S2 for x in S1):
            if return_automorphism:
                return True, b
            else:
                return True
    if return_automorphism:
        return False, None
    else:
        return False

def _lift2smallest_field(a):
    """
    INPUT:

    - ``a`` -- an element of a finite field GF(q)

    OUTPUT: the element b of the smallest subfield F of GF(q) for
    which F(b)=a.

    EXAMPLES::

        sage: from sage.coding.code_constructions import _lift2smallest_field
        sage: FF.<z> = GF(3^4,"z")
        sage: a = z^10
        sage: _lift2smallest_field(a)
        (2*z + 1, Finite Field in z of size 3^2)
        sage: a = z^40
        sage: _lift2smallest_field(a)
        (2, Finite Field of size 3)

    AUTHORS:

    - John Cremona
    """
    FF = a.parent()
    k = FF.degree()
    if k == 1:
        return a, FF
    pol = a.minimal_polynomial()
    d = pol.degree()
    if d == k:
        return a, FF
    p = FF.characteristic()
    F = GF((p, d), "z")
    b = pol.roots(F, multiplicities=False)[0]
    return b, F


def permutation_action(g, v):
    r"""
    Return permutation of rows `g * v`.

    Works on lists, matrices,
    sequences and vectors (by permuting coordinates). The code requires
    switching from `i` to `i+1` (and back again) since the :class:`SymmetricGroup`
    is, by convention, the symmetric group on the "letters" `1`, `2`, ...,
    `n` (not `0`, `1`, ..., `n-1`).

    EXAMPLES::

        sage: # needs sage.groups
        sage: V = VectorSpace(GF(3),5)
        sage: v = V([0,1,2,0,1])
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: permutation_action(g,v)
        (1, 2, 0, 0, 1)
        sage: g = G([()])
        sage: permutation_action(g,v)
        (0, 1, 2, 0, 1)
        sage: g = G([(1,2,3,4,5)])
        sage: permutation_action(g,v)
        (1, 2, 0, 1, 0)
        sage: L = Sequence([1,2,3,4,5])
        sage: permutation_action(g,L)
        [2, 3, 4, 5, 1]
        sage: MS = MatrixSpace(GF(3),3,7)
        sage: A = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
        sage: S5 = SymmetricGroup(5)
        sage: g = S5([(1,2,3)])
        sage: A
        [1 0 0 0 1 1 0]
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        sage: permutation_action(g,A)
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        [1 0 0 0 1 1 0]

    It also works on lists and is a "left action"::

        sage: # needs sage.groups
        sage: v = [0,1,2,0,1]
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: gv = permutation_action(g,v); gv
        [1, 2, 0, 0, 1]
        sage: permutation_action(g,v) == g(v)
        True
        sage: h = G([(3,4)])
        sage: gv = permutation_action(g,v)
        sage: hgv = permutation_action(h,gv)
        sage: hgv == permutation_action(h*g,v)
        True

    AUTHORS:

    - David Joyner, licensed under the GPL v2 or greater.
    """
    v_type_list = False
    if isinstance(v, list):
        v_type_list = True
        v = Sequence(v)
    if isinstance(v, Sequence_generic):
        V = v.universe()
    else:
        V = v.parent()
    n = len(list(v))
    gv = [v[g(i + 1) - 1] for i in range(n)]
    if v_type_list:
        return gv
    return V(gv)


def walsh_matrix(m0):
    """
    This is the generator matrix of a Walsh code. The matrix of
    codewords correspond to a Hadamard matrix.

    EXAMPLES::

        sage: walsh_matrix(2)
        [0 0 1 1]
        [0 1 0 1]
        sage: walsh_matrix(3)
        [0 0 0 0 1 1 1 1]
        [0 0 1 1 0 0 1 1]
        [0 1 0 1 0 1 0 1]
        sage: C = LinearCode(walsh_matrix(4)); C
        [16, 4] linear code over GF(2)
        sage: C.spectrum()
        [1, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0]

    This last code has minimum distance 8.

    REFERENCES:

    - :wikipedia:`Hadamard_matrix`
    """
    m = int(m0)
    if m == 1:
        return matrix(GF(2), 1, 2, [ 0, 1])
    if m > 1:
        row2 = [x.list() for x in walsh_matrix(m-1).augment(walsh_matrix(m-1)).rows()]
        return matrix(GF(2), m, 2**m, [[0]*2**(m-1) + [1]*2**(m-1)] + row2)
    raise ValueError("%s must be an integer > 0." % m0)

##################### main constructions #####################

def DuadicCodeEvenPair(F,S1,S2):
    r"""
    Construct the "even pair" of duadic codes associated to the
    "splitting" (see the docstring for ``_is_a_splitting``
    for the definition) S1, S2 of n.

    .. warning::

       Maybe the splitting should be associated to a sum of
       q-cyclotomic cosets mod n, where q is a *prime*.

    EXAMPLES::

        sage: from sage.coding.code_constructions import _is_a_splitting
        sage: n = 11; q = 3
        sage: C = Zmod(n).cyclotomic_cosets(q); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: _is_a_splitting(S1,S2,11)
        True
        sage: codes.DuadicCodeEvenPair(GF(q),S1,S2)
        ([11, 5] Cyclic Code over GF(3),
         [11, 5] Cyclic Code over GF(3))
    """
    from sage.misc.stopgap import stopgap
    stopgap("The function DuadicCodeEvenPair has several issues which may cause wrong results", 25896)

    from .cyclic_code import CyclicCode
    n = len(S1) + len(S2) + 1
    if not _is_a_splitting(S1,S2,n):
        raise TypeError("%s, %s must be a splitting of %s." % (S1,S2,n))
    q = F.order()
    k = Mod(q,n).multiplicative_order()
    FF = GF(q**k,"z")
    z = FF.gen()
    zeta = z**((q**k-1)/n)
    P1 = PolynomialRing(FF,"x")
    x = P1.gen()
    g1 = prod([x-zeta**i for i in S1+[0]])
    g2 = prod([x-zeta**i for i in S2+[0]])
    P2 = PolynomialRing(F,"x")
    x = P2.gen()
    gg1 = P2([_lift2smallest_field(c)[0] for c in g1.coefficients(sparse=False)])
    gg2 = P2([_lift2smallest_field(c)[0] for c in g2.coefficients(sparse=False)])
    C1 = CyclicCode(length=n, generator_pol=gg1)
    C2 = CyclicCode(length=n, generator_pol=gg2)
    return C1,C2

def DuadicCodeOddPair(F,S1,S2):
    """
    Construct the "odd pair" of duadic codes associated to the
    "splitting" S1, S2 of n.

    .. warning::

       Maybe the splitting should be associated to a sum of
       q-cyclotomic cosets mod n, where q is a *prime*.

    EXAMPLES::

        sage: from sage.coding.code_constructions import _is_a_splitting
        sage: n = 11; q = 3
        sage: C = Zmod(n).cyclotomic_cosets(q); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: _is_a_splitting(S1,S2,11)
        True
        sage: codes.DuadicCodeOddPair(GF(q),S1,S2)
        ([11, 6] Cyclic Code over GF(3),
         [11, 6] Cyclic Code over GF(3))

    This is consistent with Theorem 6.1.3 in [HP2003]_.
    """
    from sage.misc.stopgap import stopgap
    stopgap("The function DuadicCodeOddPair has several issues which may cause wrong results", 25896)

    from .cyclic_code import CyclicCode
    n = len(S1) + len(S2) + 1
    if not _is_a_splitting(S1,S2,n):
        raise TypeError("%s, %s must be a splitting of %s." % (S1,S2,n))
    q = F.order()
    k = Mod(q,n).multiplicative_order()
    FF = GF(q**k,"z")
    z = FF.gen()
    zeta = z**((q**k-1)/n)
    P1 = PolynomialRing(FF,"x")
    x = P1.gen()
    g1 = prod([x-zeta**i for i in S1+[0]])
    g2 = prod([x-zeta**i for i in S2+[0]])
    j = sum([x**i/n for i in range(n)])
    P2 = PolynomialRing(F,"x")
    x = P2.gen()
    coeffs1 = [_lift2smallest_field(c)[0] for c in (g1+j).coefficients(sparse=False)]
    coeffs2 = [_lift2smallest_field(c)[0] for c in (g2+j).coefficients(sparse=False)]
    gg1 = P2(coeffs1)
    gg2 = P2(coeffs2)
    gg1 = gcd(gg1, x**n - 1)
    gg2 = gcd(gg2, x**n - 1)
    C1 = CyclicCode(length=n, generator_pol=gg1)
    C2 = CyclicCode(length=n, generator_pol=gg2)
    return C1,C2

def ExtendedQuadraticResidueCode(n,F):
    r"""
    The extended quadratic residue code (or XQR code) is obtained from
    a QR code by adding a check bit to the last coordinate. (These
    codes have very remarkable properties such as large automorphism
    groups and duality properties - see [HP2003]_, Section 6.6.3-6.6.4.)

    INPUT:

    - ``n`` -- an odd prime

    - ``F`` -- a finite prime field whose order must be a
      quadratic residue modulo `n`

    OUTPUT: an extended quadratic residue code

    EXAMPLES::

        sage: C1 = codes.QuadraticResidueCode(7, GF(2))
        sage: C2 = C1.extended_code()
        sage: C3 = codes.ExtendedQuadraticResidueCode(7, GF(2)); C3
        Extension of [7, 4] Cyclic Code over GF(2)
        sage: C2 == C3
        True
        sage: C = codes.ExtendedQuadraticResidueCode(17, GF(2))
        sage: C
        Extension of [17, 9] Cyclic Code over GF(2)
        sage: C3 = codes.QuadraticResidueCodeOddPair(7, GF(2))[0]
        sage: C3x = C3.extended_code()
        sage: C4 = codes.ExtendedQuadraticResidueCode(7, GF(2))
        sage: C3x == C4
        True

    AUTHORS:

    - David Joyner (07-2006)
    """
    C = QuadraticResidueCodeOddPair(n,F)[0]
    return C.extended_code()

def from_parity_check_matrix(H):
    r"""
    Return the linear code that has ``H`` as a parity check matrix.

    If ``H`` has dimensions `h \times n` then the linear code will have
    dimension `n-h` and length `n`.

    EXAMPLES::

        sage: C = codes.HammingCode(GF(2), 3); C
        [7, 4] Hamming Code over GF(2)
        sage: H = C.parity_check_matrix(); H
        [1 0 1 0 1 0 1]
        [0 1 1 0 0 1 1]
        [0 0 0 1 1 1 1]
        sage: C2 = codes.from_parity_check_matrix(H); C2
        [7, 4] linear code over GF(2)
        sage: C2.systematic_generator_matrix() == C.systematic_generator_matrix()
        True
    """
    Cd = LinearCode(H)
    return Cd.dual_code()

def QuadraticResidueCode(n,F):
    r"""
    A quadratic residue code (or QR code) is a cyclic code whose
    generator polynomial is the product of the polynomials
    `x-\alpha^i` (`\alpha` is a primitive
    `n`-th root of unity; `i` ranges over the set of
    quadratic residues modulo `n`).

    See :class:`QuadraticResidueCodeEvenPair` and
    :class:`QuadraticResidueCodeOddPair` for a more general construction.

    INPUT:

    - ``n`` -- an odd prime

    - ``F`` -- a finite prime field whose order must be a
      quadratic residue modulo `n`

    OUTPUT: a quadratic residue code

    EXAMPLES::

        sage: C = codes.QuadraticResidueCode(7, GF(2))
        sage: C
        [7, 4] Cyclic Code over GF(2)
        sage: C = codes.QuadraticResidueCode(17, GF(2))
        sage: C
        [17, 9] Cyclic Code over GF(2)
        sage: C1 = codes.QuadraticResidueCodeOddPair(7, GF(2))[0]
        sage: C2 = codes.QuadraticResidueCode(7, GF(2))
        sage: C1 == C2
        True
        sage: C1 = codes.QuadraticResidueCodeOddPair(17, GF(2))[0]
        sage: C2 = codes.QuadraticResidueCode(17, GF(2))
        sage: C1 == C2
        True

    AUTHORS:

    - David Joyner (11-2005)
    """
    return QuadraticResidueCodeOddPair(n,F)[0]

def QuadraticResidueCodeEvenPair(n,F):
    r"""
    Quadratic residue codes of a given odd prime length and base ring
    either don't exist at all or occur as 4-tuples - a pair of
    "odd-like" codes and a pair of "even-like" codes. If `n > 2` is prime
    then (Theorem 6.6.2 in [HP2003]_) a QR code exists over `\GF{q}` iff q is a
    quadratic residue mod `n`.

    They are constructed as "even-like" duadic codes associated the
    splitting `(Q,N)` mod `n`, where `Q` is the set of nonzero quadratic
    residues and `N` is the non-residues.

    EXAMPLES::

        sage: codes.QuadraticResidueCodeEvenPair(17, GF(13))  # known bug (#25896)
        ([17, 8] Cyclic Code over GF(13),
         [17, 8] Cyclic Code over GF(13))
        sage: codes.QuadraticResidueCodeEvenPair(17, GF(2))
        ([17, 8] Cyclic Code over GF(2),
         [17, 8] Cyclic Code over GF(2))
        sage: codes.QuadraticResidueCodeEvenPair(13, GF(9,"z"))  # known bug (#25896)
        ([13, 6] Cyclic Code over GF(9),
         [13, 6] Cyclic Code over GF(9))
        sage: C1,C2 = codes.QuadraticResidueCodeEvenPair(7, GF(2))
        sage: C1.is_self_orthogonal()
        True
        sage: C2.is_self_orthogonal()
        True
        sage: C3 = codes.QuadraticResidueCodeOddPair(17, GF(2))[0]
        sage: C4 = codes.QuadraticResidueCodeEvenPair(17, GF(2))[1]
        sage: C3.systematic_generator_matrix() == C4.dual_code().systematic_generator_matrix()
        True

    This is consistent with Theorem 6.6.9 and Exercise 365 in [HP2003]_.

    TESTS::

        sage: codes.QuadraticResidueCodeEvenPair(14,Zmod(4))
        Traceback (most recent call last):
        ...
        ValueError: the argument F must be a finite field
        sage: codes.QuadraticResidueCodeEvenPair(14, GF(2))
        Traceback (most recent call last):
        ...
        ValueError: the argument n must be an odd prime
        sage: codes.QuadraticResidueCodeEvenPair(5, GF(2))
        Traceback (most recent call last):
        ...
        ValueError: the order of the finite field must be a quadratic residue modulo n
    """
    from sage.arith.srange import srange
    from sage.categories.finite_fields import FiniteFields
    if F not in FiniteFields():
        raise ValueError("the argument F must be a finite field")
    q = F.order()
    n = Integer(n)
    if n <= 2 or not n.is_prime():
        raise ValueError("the argument n must be an odd prime")
    Q = quadratic_residues(n)
    Q.remove(0)       # nonzero quad residues
    N = [x for x in srange(1, n) if x not in Q]   # nonzero quad non-residues
    if q not in Q:
        raise ValueError("the order of the finite field must be a quadratic residue modulo n")
    return DuadicCodeEvenPair(F,Q,N)


def QuadraticResidueCodeOddPair(n,F):
    r"""
    Quadratic residue codes of a given odd prime length and base ring
    either don't exist at all or occur as 4-tuples - a pair of
    "odd-like" codes and a pair of "even-like" codes. If n 2 is prime
    then (Theorem 6.6.2 in [HP2003]_) a QR code exists over `\GF{q} iff `q` is a
    quadratic residue mod `n`.

    They are constructed as "odd-like" duadic codes associated the
    splitting `(Q,N)` mod `n`, where `Q` is the set of nonzero quadratic
    residues and `N` is the non-residues.

    EXAMPLES::

        sage: codes.QuadraticResidueCodeOddPair(17, GF(13))  # known bug (#25896)
        ([17, 9] Cyclic Code over GF(13),
         [17, 9] Cyclic Code over GF(13))
        sage: codes.QuadraticResidueCodeOddPair(17, GF(2))
        ([17, 9] Cyclic Code over GF(2),
         [17, 9] Cyclic Code over GF(2))
        sage: codes.QuadraticResidueCodeOddPair(13, GF(9,"z"))  # known bug (#25896)
        ([13, 7] Cyclic Code over GF(9),
         [13, 7] Cyclic Code over GF(9))
        sage: C1 = codes.QuadraticResidueCodeOddPair(17, GF(2))[1]
        sage: C1x = C1.extended_code()
        sage: C2 = codes.QuadraticResidueCodeOddPair(17, GF(2))[0]
        sage: C2x = C2.extended_code()
        sage: C2x.spectrum(); C1x.spectrum()
        [1, 0, 0, 0, 0, 0, 102, 0, 153, 0, 153, 0, 102, 0, 0, 0, 0, 0, 1]
        [1, 0, 0, 0, 0, 0, 102, 0, 153, 0, 153, 0, 102, 0, 0, 0, 0, 0, 1]
        sage: C3 = codes.QuadraticResidueCodeOddPair(7, GF(2))[0]
        sage: C3x = C3.extended_code()
        sage: C3x.spectrum()
        [1, 0, 0, 0, 14, 0, 0, 0, 1]

    This is consistent with Theorem 6.6.14 in [HP2003]_.

    TESTS::

        sage: codes.QuadraticResidueCodeOddPair(9, GF(2))
        Traceback (most recent call last):
        ...
        ValueError: the argument n must be an odd prime
    """
    from sage.arith.srange import srange
    from sage.categories.finite_fields import FiniteFields
    if F not in FiniteFields():
        raise ValueError("the argument F must be a finite field")
    q = F.order()
    n = Integer(n)
    if n <= 2 or not n.is_prime():
        raise ValueError("the argument n must be an odd prime")
    Q = quadratic_residues(n)
    Q.remove(0)       # nonzero quad residues
    N = [x for x in srange(1, n) if x not in Q]   # nonzero quad non-residues
    if q not in Q:
        raise ValueError("the order of the finite field must be a quadratic residue modulo n")
    return DuadicCodeOddPair(F,Q,N)


def random_linear_code(F, length, dimension):
    r"""
    Generate a random linear code of length ``length``, dimension ``dimension``
    and over the field ``F``.

    This function is Las Vegas probabilistic: always correct, usually fast.
    Random matrices over the ``F`` are drawn until one with full rank is hit.

    If ``F`` is infinite, the distribution of the elements in the random
    generator matrix will be random according to the distribution of
    ``F.random_element()``.

    EXAMPLES::

        sage: C = codes.random_linear_code(GF(2), 10, 3)
        sage: C
        [10, 3] linear code over GF(2)
        sage: C.generator_matrix().rank()
        3
    """
    while True:
        G = random_matrix(F, dimension, length)
        if G.rank() == dimension:
            return LinearCode(G)

def ToricCode(P,F):
    r"""
    Let `P` denote a list of lattice points in
    `\ZZ^d` and let `T` denote the set of all
    points in `(F^x)^d` (ordered in some fixed way). Put
    `n=|T|` and let `k` denote the dimension of the
    vector space of functions `V = \mathrm{Span}\{x^e \ |\ e \in P\}`.
    The associated toric code `C` is the evaluation code which
    is the image of the evaluation map

    .. MATH::

        \operatorname{eval}_T : V \rightarrow F^n,


    where `x^e` is the multi-index notation
    (`x=(x_1,...,x_d)`, `e=(e_1,...,e_d)`, and
    `x^e = x_1^{e_1}...x_d^{e_d}`), where
    `\operatorname{eval}_T (f(x)) = (f(t_1),...,f(t_n))`, and where
    `T=\{t_1,...,t_n\}`. This function returns the toric
    codes discussed in [Joy2004]_.

    INPUT:

    - ``P`` -- all the integer lattice points in a polytope
      defining the toric variety

    - ``F`` -- a finite field

    OUTPUT: toric code with length `n`, dimension `k` over field `F`

    EXAMPLES::

         sage: C = codes.ToricCode([[0,0],[1,0],[2,0],[0,1],[1,1]], GF(7))
         sage: C
         [36, 5] linear code over GF(7)
         sage: C.minimum_distance()                                                     # needs sage.groups
         24
         sage: C.minimum_distance(algorithm='guava')  # optional - gap_package_guava
         ...24
         sage: C = codes.ToricCode([[-2,-2],[-1,-2],[-1,-1],[-1,0],
         ....:                      [0,-1],[0,0],[0,1],[1,-1],[1,0]], GF(5))
         sage: C
         [16, 9] linear code over GF(5)
         sage: C.minimum_distance()                                                     # needs sage.groups
         6
         sage: C.minimum_distance(algorithm='guava')  # optional - gap_package_guava
         6
         sage: C = codes.ToricCode([[0,0],[1,1],[1,2],[1,3],[1,4],[2,1],
         ....:                      [2,2],[2,3],[3,1],[3,2],[4,1]], GF(8,"a"))
         sage: C
         [49, 11] linear code over GF(8)

    This is in fact a [49,11,28] code over `\GF{8}`. If you type next
    ``C.minimum_distance()`` and wait overnight (!), you
    should get 28.

    AUTHOR:

    - David Joyner (07-2006)
    """
    from sage.combinat.tuple import Tuples
    mset = [x for x in F if x != 0]
    d = len(P[0])
    pts = Tuples(mset, d).list()
    n = len(pts)  # (q-1)^d
    k = len(P)
    e = P[0]
    B = []
    for e in P:
        tmpvar = [prod([t[i]**e[i] for i in range(d)]) for t in pts]
        B.append(tmpvar)
    # now B0 *should* be a full rank matrix
    MS = MatrixSpace(F, k, n)
    return LinearCode(MS(B))


def WalshCode(m):
    r"""
    Return the binary Walsh code of length `2^m`.

    The matrix
    of codewords correspond to a Hadamard matrix. This is a (constant
    rate) binary linear `[2^m,m,2^{m-1}]` code.

    EXAMPLES::

        sage: C = codes.WalshCode(4); C
        [16, 4] linear code over GF(2)
        sage: C = codes.WalshCode(3); C
        [8, 3] linear code over GF(2)
        sage: C.spectrum()
        [1, 0, 0, 0, 7, 0, 0, 0, 0]
        sage: C.minimum_distance()                                                      # needs sage.libs.gap
        4
        sage: C.minimum_distance(algorithm='gap')  # check d=2^(m-1)                    # needs sage.libs.gap
        4

    REFERENCES:

    - :wikipedia:`Hadamard_matrix`

    - :wikipedia:`Walsh_code`
    """
    return LinearCode(walsh_matrix(m), d=2**(m - 1))
