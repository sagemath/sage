# sage_setup: distribution = sagemath-modules
"""
Tornaria methods for computing with quadratic forms
"""
# ****************************************************************************
#       Copyright (C) 2007 Gonzalo Tornaria
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import (CRT_vectors,
                             factor,
                             gcd,
                             hilbert_symbol,
                             kronecker as kronecker_symbol,
                             prime_to_m_part)
from sage.misc.functional import is_odd
from sage.misc.misc_c import prod
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ


# TO DO -- Add second argument
#  def __call__(self,v,w=None):
#    if w is None:
#        return half(v * self._matrix_() * v)
#    else:
#      return v * self._matrix_() * w


def disc(self):
    r"""
    Return the discriminant of the quadratic form, defined as

    - `(-1)^n {\rm det}(B)` for even dimension `2n`
    - `{\rm det}(B)/2` for odd dimension

    where `2Q(x) = x^t B x`.

    This agrees with the usual discriminant for binary and ternary quadratic forms.

    EXAMPLES::

        sage: DiagonalQuadraticForm(ZZ, [1]).disc()
        1
        sage: DiagonalQuadraticForm(ZZ, [1,1]).disc()
        -4
        sage: DiagonalQuadraticForm(ZZ, [1,1,1]).disc()
        4
        sage: DiagonalQuadraticForm(ZZ, [1,1,1,1]).disc()
        16
    """
    if is_odd(self.dim()):
        # This is not so good for characteristic 2.
        return self.base_ring()(self.det() / 2)
    return (-1)**(self.dim() // 2) * self.det()


def content(self):
    r"""
    Return the GCD of the coefficients of the quadratic form.

    .. warning::

        Only works over Euclidean domains (probably just `\ZZ`).

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1, 1])
        sage: Q.matrix().gcd()
        2
        sage: Q.content()
        1
        sage: DiagonalQuadraticForm(ZZ, [1, 1]).is_primitive()
        True
        sage: DiagonalQuadraticForm(ZZ, [2, 4]).is_primitive()
        False
        sage: DiagonalQuadraticForm(ZZ, [2, 4]).primitive()
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 1 0 ]
        [ * 2 ]
    """
    return self.gcd()


# in quadratic_form.py
# def is_primitive(self):
#     """
#     Checks if the form is a multiple of another form... only over ZZ for now.
#     """
#     return self.content() == 1


# in quadratic_form.py
# def primitive(self):
#     """
#     Return a primitive quadratic forms in the similarity class of the given form.
#
#     This only works when we have GCDs... so over ZZ.
#     """
#     c = self.content()
#     new_coeffs = [self.base_ring()(a/c)  for a in self.__coeffs]
#     return QuadraticForm(self.base_ring(), self.dim(), new_coeffs)


def adjoint(self):
    """
    This gives the adjoint (integral) quadratic form associated to the
    given form, essentially defined by taking the adjoint of the matrix.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1,2,5])
        sage: Q.adjoint()
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 5 -2 ]
        [ * 1 ]

    ::

        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q.adjoint()
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 39 2 8 ]
        [ * 19 4 ]
        [ * * 8 ]
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm as QuadraticForm

    if is_odd(self.dim()):
        return QuadraticForm(self.matrix().adjoint_classical() * 2)
    return QuadraticForm(self.matrix().adjoint_classical())


def antiadjoint(self):
    """
    This gives an (integral) form such that its adjoint is the given form.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q.adjoint().antiadjoint()
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 0 -1 ]
        [ * 2 -1 ]
        [ * * 5 ]
        sage: Q.antiadjoint()                                                           # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: not an adjoint
    """
    n = self.dim()
    R = self.base_ring()
    try:
        d = R(self.disc()**(ZZ.one() / (n - 1)))
        if is_odd(n):
            return self.adjoint().scale_by_factor(R.one() / 4 / d**(n - 2))
        else:
            return self.adjoint().scale_by_factor(R.one() / d**(n - 2))
    except TypeError:
        raise ValueError("not an adjoint")


def is_adjoint(self) -> bool:
    """
    Determine if the given form is the adjoint of another form.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q.is_adjoint()                                                            # needs sage.symbolic
        False
        sage: Q.adjoint().is_adjoint()
        True
    """
    try:
        self.antiadjoint()
    except ValueError:
        return False
    return True


def reciprocal(self):
    r"""
    This gives the reciprocal quadratic form associated to the given form.

    This is defined as the multiple of the primitive adjoint with the same
    content as the given form.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,37])
        sage: Q.reciprocal()
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 37 0 0 ]
        [ * 37 0 ]
        [ * * 1 ]
        sage: Q.reciprocal().reciprocal()
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 0 0 ]
        [ * 1 0 ]
        [ * * 37 ]
        sage: Q.reciprocal().reciprocal() == Q
        True
    """
    return self.adjoint().primitive() . scale_by_factor(self.content())


def omega(self):
    r"""
    Return the content of the adjoint of the primitive associated quadratic form.

    Ref: See Dickson's "Studies in Number Theory".

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,37])
        sage: Q.omega()
        4
    """
    return self.primitive().adjoint().content()


def delta(self):
    r"""
    Return the omega of the adjoint form.

    This is the same as the omega of the reciprocal form.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,37])
        sage: Q.delta()
        148
    """
    return self.adjoint().omega()


def level__Tornaria(self):
    r"""
    Return the level of the quadratic form.

    This is defined as

    - level(`B`)    for even dimension,
    - level(`B`)/4  for odd dimension,

    where `2Q(x) = x^t\cdot B\cdot x`.

    This agrees with the usual level for even dimension.

    EXAMPLES::

        sage: DiagonalQuadraticForm(ZZ, [1]).level__Tornaria()
        1
        sage: DiagonalQuadraticForm(ZZ, [1,1]).level__Tornaria()
        4
        sage: DiagonalQuadraticForm(ZZ, [1,1,1]).level__Tornaria()
        1
        sage: DiagonalQuadraticForm(ZZ, [1,1,1,1]).level__Tornaria()
        4
    """
    return self.base_ring()(abs(self.disc()) / self.omega() / self.content()**self.dim())


def discrec(self):
    r"""
    Return the discriminant of the reciprocal form.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,37])
        sage: Q.disc()
        148
        sage: Q.discrec()
        5476
        sage: [4 * 37, 4 * 37^2]
        [148, 5476]
    """
    return self.reciprocal().disc()


# Rational equivalence

def hasse_conductor(self):
    """
    Return the Hasse conductor.

    This is the product of all primes where the Hasse invariant equals `-1`.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q.hasse_invariant(2)
        -1
        sage: Q.hasse_invariant(37)
        -1
        sage: Q.hasse_conductor()
        74

        sage: DiagonalQuadraticForm(ZZ, [1, 1, 1]).hasse_conductor()                    # needs sage.libs.pari
        1
        sage: QuadraticForm(ZZ, 3, [2, -2, 0, 2, 0, 5]).hasse_conductor()               # needs sage.libs.pari
        10
    """
    return prod([x for x, _ in factor(2 * self.level())
                 if self.hasse_invariant(x) == -1])


def clifford_invariant(self, p):
    """
    Return the Clifford invariant.

    This is the class in the Brauer group of the Clifford algebra for
    even dimension, of the even Clifford Algebra for odd dimension.

    See Lam (AMS GSM 67) p. 117 for the definition, and p. 119 for the formula
    relating it to the Hasse invariant.

    EXAMPLES:

    For hyperbolic spaces, the Clifford invariant is +1::

        sage: # needs sage.libs.pari
        sage: H = QuadraticForm(ZZ, 2, [0, 1, 0])
        sage: H.clifford_invariant(2)
        1
        sage: (H + H).clifford_invariant(2)
        1
        sage: (H + H + H).clifford_invariant(2)
        1
        sage: (H + H + H + H).clifford_invariant(2)
        1
    """
    n = self.dim() % 8
    if n == 1 or n == 2:
        s = 1
    elif n == 3 or n == 4:
        s = hilbert_symbol(-1, -self.disc(), p)
    elif n == 5 or n == 6:
        s = hilbert_symbol(-1, -1, p)
    elif n == 7 or n == 0:
        s = hilbert_symbol(-1, self.disc(), p)
    return s * self.hasse_invariant(p)


def clifford_conductor(self):
    """
    Return the product of all primes where the Clifford invariant is `-1`.

    .. NOTE::

        For ternary forms, this is the discriminant of the
        quaternion algebra associated to the quadratic space
        (i.e. the even Clifford algebra).

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q.clifford_invariant(2)
        1
        sage: Q.clifford_invariant(37)
        -1
        sage: Q.clifford_conductor()
        37

        sage: DiagonalQuadraticForm(ZZ, [1, 1, 1]).clifford_conductor()                 # needs sage.libs.pari
        2
        sage: QuadraticForm(ZZ, 3, [2, -2, 0, 2, 0, 5]).clifford_conductor()            # needs sage.libs.pari
        30

    For hyperbolic spaces, the Clifford conductor is 1::

        sage: # needs sage.libs.pari
        sage: H = QuadraticForm(ZZ, 2, [0, 1, 0])
        sage: H.clifford_conductor()
        1
        sage: (H + H).clifford_conductor()
        1
        sage: (H + H + H).clifford_conductor()
        1
        sage: (H + H + H + H).clifford_conductor()
        1
    """
    return prod([x for x, _ in factor(2 * self.level())
                 if self.clifford_invariant(x) == -1])


# Genus theory

def basiclemma(self, M):
    """
    Find a number represented by ``self`` and coprime to `M`.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [2, 1, 3])
        sage: Q.basiclemma(6)
        71
    """
    a = self(self.basiclemmavec(M))
    assert gcd(a, M) == 1
    return a


def basiclemmavec(self, M):
    """
    Find a vector where the value of the quadratic form is coprime to `M`.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [2, 1, 5])
        sage: Q.basiclemmavec(10)
        (6, 5)
        sage: Q(_)
        227
    """
    V = FreeModule(self.base_ring(), self.dim())
    vec = []
    mod = []
    M0 = abs(M)
    if M0 == 1:
        return V.zero()

    for i in range(self.dim()):
        M1 = prime_to_m_part(M0, self[i, i])
        if M1 != 1:
            vec.append(V.gen(i))
            mod.append(M1)
        M0 = M0 / M1
        if M0 == 1:
            return tuple(CRT_vectors(vec, mod))

    for i in range(self.dim()):
        for j in range(i):
            M1 = prime_to_m_part(M0, self[i, j])
            if M1 != 1:
                vec.append(V.i + V.j)
                mod.append(M1)
            M0 = M0 / M1
            if M0 == 1:
                return tuple(CRT_vectors(vec, mod))

    raise ValueError("not primitive form")


# FIXME: get the rules for validity of characters straight...
# p=2 might be bad!!!
def xi(self, p):
    """
    Return the value of the genus characters Xi_p... which may be missing one character.
    We allow -1 as a prime.

    REFERENCES:

    Dickson's "Studies in the Theory of Numbers"

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 1, 1, 14, 3, 14])
        sage: Q2 = QuadraticForm(ZZ, 3, [2, -1, 0, 2, 0, 50])
        sage: [Q1.omega(), Q2.omega()]
        [5, 5]
        sage: [Q1.hasse_invariant(5),                   # equivalent over Q_5           # needs sage.libs.pari
        ....:  Q2.hasse_invariant(5)]
        [1, 1]
        sage: [Q1.xi(5), Q2.xi(5)]                      # not equivalent over Z_5       # needs sage.libs.pari
        [1, -1]
    """
    if self.dim() == 2 and self.disc() % p:
        raise ValueError("not a valid character")
    if self.dim() >= 3 and self.omega() % p:
        raise ValueError("not a valid character")
    if (p == -1) or (p == 2):
        return kronecker_symbol(p, self.basiclemma(2))
    return kronecker_symbol(self.basiclemma(p), p)


def xi_rec(self, p):
    """
    Return Xi(`p`) for the reciprocal form.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: Q1 = QuadraticForm(ZZ, 3, [1, 1, 1, 14, 3, 14])
        sage: Q2 = QuadraticForm(ZZ, 3, [2, -1, 0, 2, 0, 50])
        sage: [Q1.clifford_conductor(),                 # equivalent over Q
        ....:  Q2.clifford_conductor()]
        [3, 3]
        sage: Q1.is_locally_equivalent_to(Q2)           # not in the same genus
        False
        sage: [Q1.delta(), Q2.delta()]
        [480, 480]
        sage: factor(480)
        2^5 * 3 * 5
        sage: list(map(Q1.xi_rec, [-1,2,3,5]))
        [-1, -1, -1, 1]
        sage: list(map(Q2.xi_rec, [-1,2,3,5]))
        [-1, -1, -1, -1]
    """
    return self.reciprocal().xi(p)


def lll(self):
    """
    Return an LLL-reduced form of `Q` (using PARI).

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, range(1,11))
        sage: Q.is_definite()
        True
        sage: Q.lll()                                                                   # needs sage.libs.pari
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 1 0 ]
        [ * 4 3 3 ]
        [ * * 6 3 ]
        [ * * * 6 ]
    """
    return self(self.matrix().LLL_gram())


def representation_number_list(self, B):
    """
    Return the vector of representation numbers < B.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1,1,1,1])
        sage: Q.representation_number_list(10)                                          # needs sage.libs.pari
        [1, 16, 112, 448, 1136, 2016, 3136, 5504, 9328, 12112]
    """
    from sage.libs.pari.all import pari

    ans = pari(1).concat(self.__pari__().qfrep(B - 1, 1) * 2)
    return ans.sage()


def representation_vector_list(self, B, maxvectors=10**8):
    """
    Find all vectors `v` where `Q(v) < B`.

    This only works for positive definite quadratic forms.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: Q = DiagonalQuadraticForm(ZZ, [1, 1])
        sage: Q.representation_vector_list(10)
        [[(0, 0)],
         [(0, 1), (0, -1), (1, 0), (-1, 0)],
         [(1, 1), (-1, -1), (1, -1), (-1, 1)],
         [],
         [(0, 2), (0, -2), (2, 0), (-2, 0)],
         [(1, 2), (-1, -2), (1, -2), (-1, 2), (2, 1), (-2, -1), (2, -1), (-2, 1)],
         [],
         [],
         [(2, 2), (-2, -2), (2, -2), (-2, 2)],
         [(0, 3), (0, -3), (3, 0), (-3, 0)]]
        sage: list(map(len, _))
        [1, 4, 4, 0, 4, 8, 0, 0, 4, 4]
        sage: Q.representation_number_list(10)
        [1, 4, 4, 0, 4, 8, 0, 0, 4, 4]

    TESTS::

        sage: R = QuadraticForm(ZZ, 2, [-4,-3,0])
        sage: R.representation_vector_list(10)                                          # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        PariError: domain error in minim0: form is not positive definite
    """
    n, _, vs = self.__pari__().qfminim(2 * (B - 1), maxvectors)

    if n != 2 * len(vs):
        raise RuntimeError("insufficient number of vectors")
    ms = [[] for _ in range(B)]
    ms[0] = [vector([0] * self.dim())]
    for v in vs.sage().columns():
        ms[int(self(v))] += [v, -v]
    return ms


# zeros

def is_zero(self, v, p=0) -> bool:
    r"""
    Determine if the vector `v` is on the conic `Q(x) = 0` (mod `p`).

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q1.is_zero([0,1,0], 2)
        True
        sage: Q1.is_zero([1,1,1], 2)
        True
        sage: Q1.is_zero([1,1,0], 2)
        False
    """
    norm = self(v)
    if p != 0:
        norm = norm % p
    return norm == 0


def is_zero_nonsingular(self, v, p=0) -> bool:
    """
    Determine if the vector `v` is on the conic `Q(x) = 0` (mod `p`),
    and that this point is non-singular point of the conic.

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q1.is_zero_nonsingular([1,1,1], 2)
        True
        sage: Q1.is_zero([1, 19, 2], 37)
        True
        sage: Q1.is_zero_nonsingular([1, 19, 2], 37)
        False
    """
    if not self.is_zero(v, p):
        return False
    vm = vector(self.base_ring(), v) * self.matrix()
    if p != 0:
        vm = vm % p
    return (vm != 0)


def is_zero_singular(self, v, p=0) -> bool:
    """
    Determine if the vector `v` is on the conic `Q(x) = 0` (mod `p`),
    and that this point is singular point of the conic.

    EXAMPLES::

        sage: Q1 = QuadraticForm(ZZ, 3, [1, 0, -1, 2, -1, 5])
        sage: Q1.is_zero([1,1,1], 2)
        True
        sage: Q1.is_zero_singular([1,1,1], 2)
        False
        sage: Q1.is_zero_singular([1, 19, 2], 37)
        True
    """
    if not self.is_zero(v, p):
        return False
    vm = vector(self.base_ring(), v) * self.matrix()
    if p != 0:
        vm = vm % p
    return bool(vm == 0)
