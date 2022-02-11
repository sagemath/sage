r"""
Brandt Modules

Introduction
------------

The construction of Brandt modules provides us with a method to
compute modular forms, as outlined in Pizer's paper [Piz1980]_.

Given a prime number `p` and a positive integer `M` with `p\nmid M`,
the *Brandt module* `B(p, M)` is the free abelian group on right ideal
classes of a quaternion order of level `pM` in the quaternion algebra
ramified precisely at the places `p` and `\infty`. This Brandt module
carries a natural Hecke action given by Brandt matrices. There exists
a non-canonical Hecke algebra isomorphism between `B(p, M)` and a
certain subspace of `S_{2}(\Gamma_0(pM))` containing the newforms.

Quaternion Algebras
-------------------

A quaternion algebra over `\QQ` is a central simple algebra of
dimension 4 over `\QQ`. Such an algebra `A` is said to be ramified at
a place `v` of `\QQ` if and only if `A \otimes \QQ_v` is a division
algebra. Otherwise `A` is said to be split at `v`.

``A = QuaternionAlgebra(p)`` returns the quaternion algebra `A` over
`\QQ` ramified precisely at the places `p` and `\infty`.

``A = QuaternionAlgebra(a, b)`` returns the quaternion algebra `A`
over `\QQ` with basis `\{1, i, j, k\}` such that `i^2 = a`, `j^2 = b`
and `ij = -ji = k.`

An order `R` in a quaternion algebra `A` over `\QQ` is a 4-dimensional
lattice in `A` which is also a subring containing the identity.  A
maximal order is one that is not properly contained in another order.

A particularly important kind of orders are those that have a level;
see Definition 1.2 in [Piz1980]_.  This is a positive integer `N` such
that every prime that ramifies in `A` divides `N` to an odd power.
The maximal orders are those that have level equal to the discriminant
of `A`.

``R = A.maximal_order()`` returns a maximal order `R` in the quaternion
algebra `A.`

A right `\mathcal{O}`-ideal `I` is a lattice in `A` such that for
every prime `p` there exists `a_p\in A_p^*` with `I_p =
a_p\mathcal{O}_p`. Two right `\mathcal{O}`-ideals `I` and `J` are said
to belong to the same class if `I=aJ` for some `a \in A^*`. Left
`\mathcal{O}`-ideals are defined in a similar fashion.

The right order of `I` is the subring of `A` consisting of elements
`a` with `Ia \subseteq I`.

Brandt Modules
--------------

``B = BrandtModule(p, M=1)`` returns the Brandt module associated to
the prime number `p` and the integer `M`, with `p` not dividing `M`.

``A = B.quaternion_algebra()`` returns the quaternion algebra attached
to `B`; this is the quaternion algebra over `\QQ` ramified exactly at
`p` and `\infty`.

``O = B.order_of_level_N()`` returns an order `\mathcal{O}` of level
`N = pM` in `A`.

``B.right_ideals()`` returns a tuple of representatives for all right
ideal classes of `\mathcal{O}`.

The implementation of this method is especially interesting. It
depends on the construction of a Hecke module defined as a free
abelian group on right ideal classes of a quaternion algebra with the
following action:

.. MATH::

    T_n[I] = \sum_{\phi} [J]

where `(n,pM)=1` and the sum is over cyclic `\mathcal{O}`-module
homomorphisms `\phi\colon I\rightarrow J` of degree `n` up to
isomorphism of `J`. Equivalently one can sum over the inclusions of
the submodules `J \rightarrow n^{-1}I`. The rough idea is to start
with the trivial ideal class containing the order `\mathcal{O}`
itself. Using the method ``cyclic_submodules(self, I, q)`` one then
repeatedly computes `T_q([\mathcal{O}])` for some prime `q` not
dividing the level of `\mathcal{O}` and tests for equivalence among
the resulting ideals. A theorem of Serre asserts that one gets a
complete set of ideal class representatives after a finite number of
repetitions.

One can prove that two ideals `I` and `J` are equivalent if and only
if there exists an element `\alpha \in I \overline{J}` such
`N(\alpha)=N(I)N(J)`.

``is_equivalent(I,J)`` returns true if `I` and `J` are equivalent. This
method first compares the theta series of `I` and `J`. If they are the
same, it computes the theta series of the lattice `I\overline(J)`. It
returns true if the `n^{th}` coefficient of this series is nonzero
where `n=N(J)N(I)`.

The theta series of a lattice `L` over the quaternion algebra `A` is
defined as

.. MATH::

    \theta_L(q)=\sum_{x \in L} q^{\frac{N(x)}{N(L)}}

``L.theta_series(T,q)`` returns a power series representing `\theta_L(q)`
up to a precision of `\mathcal{O}(q^{T+1})`.


Hecke Structure
---------------

The Hecke structure defined on the Brandt module is given by the
Brandt matrices which can be computed using the definition of the
Hecke operators given earlier.

``hecke_matrix_from_defn(self,n)`` returns the matrix of the n-th Hecke
operator `B_{0}(n)` acting on self, computed directly from the
definition.

However, one can efficiently compute Brandt matrices using theta
series. In fact, let `\{I_{1},.....,I_{h}\}` be a set of right
`\mathcal{O}`-ideal class representatives. The (i,j) entry in the
Brandt matrix `B_{0}(n)` is the product of the `n^{th}` coefficient in
the theta series of the lattice `I_{i}\overline{I_{j}}` and the first
coefficient in the theta series of the lattice
`I_{i}\overline{I_{i}}`.

``compute_hecke_matrix_brandt(self,n)`` returns the n-th Hecke matrix,
computed using theta series.

EXAMPLES::

    sage: B = BrandtModule(23)

    sage: B.maximal_order()
    Order of Quaternion Algebra (-1, -23) with base ring Rational Field with basis (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)

    sage: B.right_ideals()
    (Fractional ideal (2 + 2*j, 2*i + 2*k, 4*j, 4*k), Fractional ideal (2 + 2*j, 2*i + 6*k, 8*j, 8*k), Fractional ideal (2 + 10*j + 8*k, 2*i + 8*j + 6*k, 16*j, 16*k))

    sage: B.hecke_matrix(2)
    [1 2 0]
    [1 1 1]
    [0 3 0]

    sage: B.brandt_series(3)
    [1/4 + q + q^2 + O(q^3)     1/4 + q^2 + O(q^3)           1/4 + O(q^3)]
    [  1/2 + 2*q^2 + O(q^3) 1/2 + q + q^2 + O(q^3)   1/2 + 3*q^2 + O(q^3)]
    [          1/6 + O(q^3)     1/6 + q^2 + O(q^3)       1/6 + q + O(q^3)]


REFERENCES:

- [Piz1980]_
- [Koh2000]_

Further Examples
----------------

We decompose a Brandt module over both `\ZZ` and `\QQ`. ::

    sage: B = BrandtModule(43, base_ring=ZZ); B
    Brandt module of dimension 4 of level 43 of weight 2 over Integer Ring
    sage: D = B.decomposition()
    sage: D
    [
    Subspace of dimension 1 of Brandt module of dimension 4 of level 43 of weight 2 over Integer Ring,
    Subspace of dimension 1 of Brandt module of dimension 4 of level 43 of weight 2 over Integer Ring,
    Subspace of dimension 2 of Brandt module of dimension 4 of level 43 of weight 2 over Integer Ring
    ]
    sage: D[0].basis()
    ((0, 0, 1, -1),)
    sage: D[1].basis()
    ((1, 2, 2, 2),)
    sage: D[2].basis()
    ((1, 1, -1, -1), (0, 2, -1, -1))
    sage: B = BrandtModule(43, base_ring=QQ); B
    Brandt module of dimension 4 of level 43 of weight 2 over Rational Field
    sage: B.decomposition()[2].basis()
    ((1, 0, -1/2, -1/2), (0, 1, -1/2, -1/2))

AUTHORS:

- Jon Bober
- Alia Hamieh
- Victoria de Quehen
- William Stein
- Gonzalo Tornaria
"""

# ****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# imports
from sage.misc.misc_c import prod
from sage.misc.verbose import verbose
from sage.rings.all import Integer, ZZ, QQ, PolynomialRing, GF, CommutativeRing

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra, basis_for_quaternion_lattice
from sage.algebras.quatalg.quaternion_algebra_cython import rational_matrix_from_rational_quaternions

from sage.arith.all import gcd, factor, prime_divisors, kronecker, next_prime
from sage.modular.hecke.all import (AmbientHeckeModule, HeckeSubmodule,
                                    HeckeModuleElement)
from sage.modular.dirichlet import TrivialCharacter
from sage.matrix.all import MatrixSpace, matrix
from sage.structure.richcmp import richcmp, richcmp_method
from sage.misc.cachefunc import cached_method


cache = {}


def BrandtModule(N, M=1, weight=2, base_ring=QQ, use_cache=True):
    """
    Return the Brandt module of given weight associated to the prime
    power `p^r` and integer `M`, where `p` and `M` are coprime.

    INPUT:

    - `N` -- a product of primes with odd exponents
    - `M` -- an integer coprime to `q` (default: 1)
    - ``weight`` -- an integer that is at least 2 (default: 2)
    - ``base_ring`` -- the base ring (default: ``QQ``)
    - ``use_cache`` -- whether to use the cache (default: ``True``)

    OUTPUT:

    a Brandt module

    EXAMPLES::

        sage: BrandtModule(17)
        Brandt module of dimension 2 of level 17 of weight 2 over Rational Field
        sage: BrandtModule(17,15)
        Brandt module of dimension 32 of level 17*15 of weight 2 over Rational Field
        sage: BrandtModule(3,7)
        Brandt module of dimension 2 of level 3*7 of weight 2 over Rational Field
        sage: BrandtModule(3,weight=2)
        Brandt module of dimension 1 of level 3 of weight 2 over Rational Field
        sage: BrandtModule(11, base_ring=ZZ)
        Brandt module of dimension 2 of level 11 of weight 2 over Integer Ring
        sage: BrandtModule(11, base_ring=QQbar)
        Brandt module of dimension 2 of level 11 of weight 2 over Algebraic Field

    The ``use_cache`` option determines whether the Brandt module returned
    by this function is cached::

        sage: BrandtModule(37) is BrandtModule(37)
        True
        sage: BrandtModule(37,use_cache=False) is BrandtModule(37,use_cache=False)
        False

    TESTS:

    Note that `N` and `M` must be coprime::

        sage: BrandtModule(3,15)
        Traceback (most recent call last):
        ...
        ValueError: M must be coprime to N

    Only weight 2 is currently implemented::

        sage: BrandtModule(3,weight=4)
        Traceback (most recent call last):
        ...
        NotImplementedError: weight != 2 not yet implemented

    Brandt modules are cached::

        sage: B = BrandtModule(3,5,2,ZZ)
        sage: B is BrandtModule(3,5,2,ZZ)
        True
    """
    N, M, weight = Integer(N), Integer(M), Integer(weight)
    if not N.is_prime():
        raise NotImplementedError("Brandt modules currently only implemented when N is a prime")
    if M < 1:
        raise ValueError("M must be positive")
    if M.gcd(N) != 1:
        raise ValueError("M must be coprime to N")
    if weight < 2:
        raise ValueError("weight must be at least 2")
    if not isinstance(base_ring, CommutativeRing):
        raise TypeError("base_ring must be a commutative ring")
    key = (N, M, weight, base_ring)
    if use_cache:
        if key in cache:  # TODO: re-enable caching!
            return cache[key]
    if weight != 2:
        raise NotImplementedError("weight != 2 not yet implemented")
    B = BrandtModule_class(*key)
    if use_cache:
        cache[key] = B
    return B


def class_number(p, r, M):
    r"""
    Return the class number of an order of level `N = p^r M` in the
    quaternion algebra over `\QQ` ramified precisely at `p` and infinity.

    This is an implementation of Theorem 1.12 of [Piz1980]_.

    INPUT:

    - `p` -- a prime
    - `r` -- an odd positive integer (default: 1)
    - `M` -- an integer coprime to `q` (default: 1)

    OUTPUT:

    Integer

    EXAMPLES::

        sage: sage.modular.quatalg.brandt.class_number(389,1,1)
        33
        sage: sage.modular.quatalg.brandt.class_number(389,1,2)  # TODO -- right?
        97
        sage: sage.modular.quatalg.brandt.class_number(389,3,1)  # TODO -- right?
        4892713
    """
    N = M * p**r
    D = prime_divisors(M)
    s = 0
    t = 0
    if N % 4:
        s = (1 - kronecker(-4, p)) / 4 * prod(1 + kronecker(-4, q) for q in D)
    if N % 9:
        t = (1 - kronecker(-3, p)) / 3 * prod(1 + kronecker(-3, q) for q in D)
    h = (N / Integer(12)) * (1 - 1 / p) * prod(1 + 1 / q for q in D) + s + t
    return Integer(h)


def maximal_order(A):
    """
    Return a maximal order in the quaternion algebra ramified
    at `p` and infinity.

    This is an implementation of Proposition 5.2 of [Piz1980]_.

    INPUT:

    - `A` -- quaternion algebra ramified precisely at `p` and infinity

    OUTPUT:

    a maximal order in `A`

    EXAMPLES::

        sage: A = BrandtModule(17).quaternion_algebra()
        sage: sage.modular.quatalg.brandt.maximal_order(A)
        Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis (1/2 + 1/2*i, 1/2*j - 1/2*k, -1/3*i + 1/3*k, -k)

        sage: A = QuaternionAlgebra(17,names='i,j,k')
        sage: A.maximal_order()
        Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis (1/2 + 1/2*i, 1/2*j - 1/2*k, -1/3*i + 1/3*k, -k)
    """
    return A.maximal_order()


def basis_for_left_ideal(R, gens):
    """
    Return a basis for the left ideal of `R` with given generators.

    INPUT:

    - `R` -- quaternion order
    - ``gens`` -- list of elements of `R`

    OUTPUT:

    list of four elements of `R`

    EXAMPLES::

        sage: B = BrandtModule(17); A = B.quaternion_algebra(); i,j,k = A.gens()
        sage: sage.modular.quatalg.brandt.basis_for_left_ideal(B.maximal_order(), [i+j,i-j,2*k,A(3)])
        [1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k]
        sage: sage.modular.quatalg.brandt.basis_for_left_ideal(B.maximal_order(), [3*(i+j),3*(i-j),6*k,A(3)])
        [3/2 + 1/2*i + k, i + 2*k, 3/2*j + 3/2*k, 3*k]
    """
    return basis_for_quaternion_lattice([b * g for b in R.basis() for g in gens])


def right_order(R, basis):
    """
    Given a basis for a left ideal `I`, return the right order in the
    quaternion order `R` of elements `x` such that `I x` is contained in `I`.

    INPUT:

    - `R` -- order in quaternion algebra
    - ``basis`` -- basis for an ideal `I`

    OUTPUT:

    order in quaternion algebra

    EXAMPLES:

    We do a consistency check with the ideal equal to a maximal order::

        sage: B = BrandtModule(17); basis = sage.modular.quatalg.brandt.basis_for_left_ideal(B.maximal_order(), B.maximal_order().basis())
        sage: sage.modular.quatalg.brandt.right_order(B.maximal_order(), basis)
        Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis (1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k)
        sage: basis
        [1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k]

        sage: B = BrandtModule(17); A = B.quaternion_algebra(); i,j,k = A.gens()
        sage: basis = sage.modular.quatalg.brandt.basis_for_left_ideal(B.maximal_order(), [i*j-j])
        sage: sage.modular.quatalg.brandt.right_order(B.maximal_order(), basis)
        Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis (1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k)
    """
    # Compute matrix of multiplication by each element of the basis.
    B = R.basis()
    Z = R.quaternion_algebra()
    M = MatrixSpace(QQ, 4)

    # I = matrix with rows the given basis for I
    I = M([list(f) for f in basis])

    # psi = matrix of right multiplication on each basis element
    psi = [M([list(f * x) for x in Z.basis()]) for f in basis]

    # invert them
    psi_inv = [x**(-1) for x in psi]

    # apply the four inverses to I
    W = [I * x for x in psi_inv]

    # The right order is the intersection of the row span of the W with the row span of B.
    X = M([list(b) for b in B]).row_module(ZZ)
    for A in W:
        X = X.intersection(A.row_module(ZZ))
    C = [Z(list(b)) for b in X.basis()]
    return Z.quaternion_order(C)


def quaternion_order_with_given_level(A, level):
    """
    Return an order in the quaternion algebra A with given level.

    This is implemented only when the base field is the rational numbers.

    INPUT:

    - ``level`` -- The level of the order to be returned. Currently this
      is only implemented when the level is divisible by at
      most one power of a prime that ramifies in this quaternion algebra.

    EXAMPLES::

        sage: from sage.modular.quatalg.brandt import quaternion_order_with_given_level, maximal_order
        sage: A.<i,j,k> = QuaternionAlgebra(5)
        sage: level = 2 * 5 * 17
        sage: O = quaternion_order_with_given_level(A, level)
        sage: M = maximal_order(A)
        sage: L = O.free_module()
        sage: N = M.free_module()
        sage: L.index_in(N) == level/5  #check that the order has the right index in the maximal order
        True
    """
    if A.base_ring() is not QQ:
        raise NotImplementedError("base field must be rational numbers")

    if len(A.ramified_primes()) > 1:
        raise NotImplementedError("Currently this algorithm only works when the quaternion algebra is only ramified at one finite prime.")

    # (The algorithm we use is similar to that in Magma (by David Kohel).)
    # in the following magma code, M denotes the level
    level = abs(level)
    N = A.discriminant()
    N1 = gcd(level, N)
    M1 = level / N1

    O = maximal_order(A)
    # if N1 != 1:
    #     # we do not know why magma does the following, so we do not do it.
    #     for p in A.ramified_primes():
    #         if not (level % p**2):
    #             raise NotImplementedError("Currently sage can only compute orders whose level is divisible by at most one power of any prime that ramifies in the quaternion algebra")

    #     P = basis_for_left_ideal(O, [N1] + [x * y - y * x
    #                                         for x in A.basis()
    #                                         for y in A.basis()])
    #     O = A.quaternion_order(P)

    fact = factor(M1)
    B = O.basis()

    for (p, r) in fact:
        a = int(-p) // 2
        for v in GF(p)**4:
            x = sum([int(v[i] + a) * B[i] for i in range(4)])
            D = x.reduced_trace()**2 - 4 * x.reduced_norm()
            # x = O.random_element((-p/2).floor(), (p/2).ceil())
            if kronecker(D, p) == 1:
                break
        X = PolynomialRing(GF(p), 'x').gen()
        a = ZZ((X**2 - ZZ(x.reduced_trace()) * X + ZZ(x.reduced_norm())).roots()[0][0])
        I = basis_for_left_ideal(O, [p**r, (x - a)**r])
        O = right_order(O, I)
        # right_order returns the RightOrder of I inside O, so we
        # do not need to do another intersection

    return O


class BrandtSubmodule(HeckeSubmodule):
    def _repr_(self):
        """
        Return string representation of this Brandt submodule.

        EXAMPLES::

            sage: BrandtModule(11)[0]._repr_()
            'Subspace of dimension 1 of Brandt module of dimension 2 of level 11 of weight 2 over Rational Field'
        """
        return "Subspace of dimension %s of %s" % (self.dimension(), self.ambient_module())


class BrandtModuleElement(HeckeModuleElement):
    def __init__(self, parent, x):
        """
        EXAMPLES::

            sage: B = BrandtModule(37)
            sage: x = B([1,2,3]); x
            (1, 2, 3)
            sage: parent(x)
            Brandt module of dimension 3 of level 37 of weight 2 over Rational Field
        """
        if isinstance(x, HeckeModuleElement):
            x = x.element()
        HeckeModuleElement.__init__(self, parent, parent.free_module()(x))

    def _richcmp_(self, other, op):
        """
        EXAMPLES::

            sage: B = BrandtModule(13,5)
            sage: B.0
            (1, 0, 0, 0, 0, 0)
            sage: B.0 == B.1
            False
            sage: B.0 == 0
            False
            sage: B(0) == 0
            True
            sage: B.0 + 2*B.1 == 2*B.1 + B.0
            True
            sage: loads(dumps(B.0)) == B.0
            True
        """
        return richcmp(self.element(), other.element(), op)

    def monodromy_pairing(self, x):
        """
        Return the monodromy pairing of ``self`` and ``x``.

        EXAMPLES::

            sage: B = BrandtModule(5,13)
            sage: B.monodromy_weights()
            (1, 3, 1, 1, 1, 3)
            sage: (B.0 + B.1).monodromy_pairing(B.0 + B.1)
            4

        TESTS:

        One check for :trac:`12866`::

            sage: Br = BrandtModule(2,7)
            sage: g1, g2 = Br.basis()
            sage: g = g1 - g2
            sage: g.monodromy_pairing(g)
            6
        """
        B = self.parent()
        w = B.monodromy_weights()
        x = B(x).element()
        v = self.element()
        return sum(x[i] * v[i] * w[i] for i in range(len(v)))

    def __mul__(self, right):
        """
        Return the monodromy pairing of ``self`` and ``right``.

        EXAMPLES::

            sage: B = BrandtModule(7,10)
            sage: B.monodromy_weights()
            (1, 1, 1, 2, 1, 1, 2, 1, 1, 1)
            sage: B.0 * B.0
            1
            sage: B.3 * B.3
            2
            sage: (B.0+B.3) * (B.0 + B.1 + 2*B.3)
            5
        """
        return self.monodromy_pairing(right)

    def _add_(self, right):
        """
        Return the sum of ``self`` and ``right``.

        EXAMPLES::

            sage: B = BrandtModule(11)
            sage: B.0 + B.1 # indirect doctest
            (1, 1)
        """
        return BrandtModuleElement(self.parent(), self.element() + right.element())

    def _sub_(self, right):
        """
        Return the difference of ``self`` and ``right``.

        EXAMPLES::

            sage: B = BrandtModule(11)
            sage: B.0 - B.1 # indirect doctest
            (1, -1)
        """
        return BrandtModuleElement(self.parent(), self.element() - right.element())

    def _neg_(self):
        """
        Return the opposite of ``self``.

        EXAMPLES::

            sage: B = BrandtModule(11)
            sage: -B.0 # indirect doctest
            (-1, 0)
        """
        return BrandtModuleElement(self.parent(), -self.element())


@richcmp_method
class BrandtModule_class(AmbientHeckeModule):
    """
    A Brandt module.

    EXAMPLES::

        sage: BrandtModule(3, 10)
        Brandt module of dimension 4 of level 3*10 of weight 2 over Rational Field
    """
    def __init__(self, N, M, weight, base_ring):
        """
        INPUT:

        - N -- ramification number (coprime to M)
        - M -- auxiliary level
        - weight -- integer 2
        - base_ring -- the base ring

        EXAMPLES::

            sage: BrandtModule(3, 5, weight=2, base_ring=ZZ)
            Brandt module of dimension 2 of level 3*5 of weight 2 over Integer Ring
        """
        assert weight == 2
        self.__N = N
        self.__M = M
        if not N.is_prime():
            raise NotImplementedError("right now N must be prime")
        rank = class_number(N, 1, M)
        self.__key = (N, M, weight, base_ring)
        AmbientHeckeModule.__init__(self, base_ring, rank, N * M, weight=2)
        self._populate_coercion_lists_(coerce_list=[self.free_module()])

    Element = BrandtModuleElement

    def _submodule_class(self):
        """
        Return the Python class of submodules of this ambient Brandt module.

        EXAMPLES::

            sage: BrandtModule(37)._submodule_class()
            <class 'sage.modular.quatalg.brandt.BrandtSubmodule'>
        """
        return BrandtSubmodule

    @cached_method
    def free_module(self):
        """
        Return the underlying free module of the Brandt module.

        EXAMPLES::

            sage: B = BrandtModule(10007,389)
            sage: B.free_module()
            Vector space of dimension 325196 over Rational Field
        """
        return self.base_ring() ** self.dimension()

    def N(self):
        """
        Return ramification level `N`.

        EXAMPLES::

            sage: BrandtModule(7,5,2,ZZ).N()
            7
        """
        return self.__N

    def M(self):
        """
        Return the auxiliary level (prime to `p` part) of the quaternion
        order used to compute this Brandt module.

        EXAMPLES::

            sage: BrandtModule(7,5,2,ZZ).M()
            5
        """
        return self.__M

    def character(self):
        r"""
        The character of this space.

        Always trivial.

        EXAMPLES::

            sage: BrandtModule(11,5).character()
            Dirichlet character modulo 55 of conductor 1 mapping 12 |--> 1, 46 |--> 1
        """
        return TrivialCharacter(self.__N * self.__M)

    def _repr_(self):
        """
        Return string representation of this Brandt module.

        EXAMPLES::

            sage: BrandtModule(7,5,2,ZZ)._repr_()
            'Brandt module of dimension 4 of level 7*5 of weight 2 over Integer Ring'
        """
        aux = '' if self.__M == 1 else '*%s' % self.__M
        txt = "Brandt module of dimension %s of level %s%s of weight %s over %s"
        return txt % (self.rank(), self.__N, aux, self.weight(), self.base_ring())

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: BrandtModule(37, 5, 2, ZZ) == BrandtModule(37, 5, 2, QQ)
            False
            sage: BrandtModule(37, 5, 2, ZZ) == BrandtModule(37, 5, 2, ZZ)
            True
            sage: BrandtModule(37, 5, 2, ZZ) == loads(dumps(BrandtModule(37, 5, 2, ZZ)))
            True
        """
        if not isinstance(other, BrandtModule_class):
            return NotImplemented

        return richcmp((self.__M, self.__N, self.weight(), self.base_ring()),
                       (other.__M, other.__N, other.weight(), other.base_ring()),
                       op)

    @cached_method
    def quaternion_algebra(self):
        r"""
        Return the quaternion algebra `A` over `\QQ` ramified precisely at
        `p` and infinity used to compute this Brandt module.

        EXAMPLES::

            sage: BrandtModule(997).quaternion_algebra()
            Quaternion Algebra (-2, -997) with base ring Rational Field
            sage: BrandtModule(2).quaternion_algebra()
            Quaternion Algebra (-1, -1) with base ring Rational Field
            sage: BrandtModule(3).quaternion_algebra()
            Quaternion Algebra (-1, -3) with base ring Rational Field
            sage: BrandtModule(5).quaternion_algebra()
            Quaternion Algebra (-2, -5) with base ring Rational Field
            sage: BrandtModule(17).quaternion_algebra()
            Quaternion Algebra (-3, -17) with base ring Rational Field
        """
        return QuaternionAlgebra(self.N())

    @cached_method
    def maximal_order(self):
        """
        Return a maximal order in the quaternion algebra associated to this Brandt module.

        EXAMPLES::

            sage: BrandtModule(17).maximal_order()
            Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis (1/2 + 1/2*i, 1/2*j - 1/2*k, -1/3*i + 1/3*k, -k)
            sage: BrandtModule(17).maximal_order() is BrandtModule(17).maximal_order()
            True
        """
        return maximal_order(self.quaternion_algebra())

    @cached_method
    def order_of_level_N(self):
        """
        Return an order of level `N = p^{2 r + 1} M` in the
        quaternion algebra.

        EXAMPLES::

            sage: BrandtModule(7).order_of_level_N()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)
            sage: BrandtModule(7,13).order_of_level_N()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j + 12*k, 1/2*i + 9/2*k, j + 11*k, 13*k)
            sage: BrandtModule(7,3*17).order_of_level_N()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j + 35*k, 1/2*i + 65/2*k, j + 19*k, 51*k)
        """
        return quaternion_order_with_given_level(self.quaternion_algebra(), self.level())

    def cyclic_submodules(self, I, p):
        """
        Return a list of rescaled versions of the fractional right
        ideals `J` such that `J` contains `I` and the quotient has
        group structure the product of two cyclic groups of order `p`.

        We emphasize again that `J` is rescaled to be integral.

        INPUT:

        - `I` -- ideal I in R = self.order_of_level_N()
        - `p` -- prime `p` coprime to self.level()

        OUTPUT:

        list of the `p+1` fractional right R-ideals that contain I
        such that J/I is GF(p) x GF(p).

        EXAMPLES::

            sage: B = BrandtModule(11)
            sage: I = B.order_of_level_N().unit_ideal()
            sage: B.cyclic_submodules(I, 2)
            [Fractional ideal (1/2 + 3/2*j + k, 1/2*i + j + 1/2*k, 2*j, 2*k),
             Fractional ideal (1/2 + 1/2*i + 1/2*j + 1/2*k, i + k, j + k, 2*k),
             Fractional ideal (1/2 + 1/2*j + k, 1/2*i + j + 3/2*k, 2*j, 2*k)]
            sage: B.cyclic_submodules(I, 3)
            [Fractional ideal (1/2 + 1/2*j, 1/2*i + 5/2*k, 3*j, 3*k),
             Fractional ideal (1/2 + 3/2*j + 2*k, 1/2*i + 2*j + 3/2*k, 3*j, 3*k),
             Fractional ideal (1/2 + 3/2*j + k, 1/2*i + j + 3/2*k, 3*j, 3*k),
             Fractional ideal (1/2 + 5/2*j, 1/2*i + 1/2*k, 3*j, 3*k)]
            sage: B.cyclic_submodules(I, 11)
            Traceback (most recent call last):
            ...
            ValueError: p must be coprime to the level
        """
        if not Integer(p).is_prime():
            raise ValueError("p must be a prime")
        if not(self.level() % p):
            raise ValueError("p must be coprime to the level")

        R = self.order_of_level_N()
        A = R.quaternion_algebra()
        B = R.basis()
        V = GF(p)**4

        # step 1: Compute alpha, beta, and the matrix of their action on I/pI.
        # NOTE: Move this code to orders once we have it all working...
        try:
            alpha, beta = self.__cyclic_submodules[p]
            compute = False
        except AttributeError:
            self.__cyclic_submodules = {}
            compute = True
        except KeyError:
            compute = True

        if compute:
            d = R.free_module().basis_matrix().determinant()
            S = None
            for v in V:
                if not v:
                    continue
                alpha = sum(Integer(v[i]) * B[i] for i in range(4))
                # If the quadratic polynomial over GF(p) given by
                #      X^2  -  alpha.reduced_trace() * X  +  alpha.reduced_norm()
                # is not irreducible, we try again with a new element.
                if p == 2:
                    # special case p == 2, since there is a unique quadratic irreducible poly.
                    if alpha.reduced_trace() % 2 == 0 or alpha.reduced_norm() % 2 == 0:
                        continue
                else:
                    # check if the discriminant is a square -- if so, poly is reducible
                    b = alpha.reduced_trace()
                    c = alpha.reduced_norm()
                    if kronecker(b * b - 4 * c, p) != -1:
                        continue
                for w in V:
                    if not w:
                        continue
                    beta = sum(Integer(w[i]) * B[i] for i in range(4))
                    v = [A(1), alpha, beta, alpha * beta]
                    M = rational_matrix_from_rational_quaternions(v)
                    e = M.determinant()
                    if e and not((d / e).valuation(p)):
                        S = A.quaternion_order(v)
                        break
                if S is not None:
                    break
            self.__cyclic_submodules[p] = (alpha, beta)

        # right multiplication by X changes something to be written
        # in terms of the basis for I.
        Y = I.basis_matrix()
        X = Y**(-1)

        # Compute the matrix of right multiplication by alpha acting on
        # our fixed choice of basis for this ideal.

        M_alpha = (matrix([(i * alpha).coefficient_tuple()
                           for i in I.basis()]) * X).change_ring(GF(p))
        M_beta = (matrix([(i * beta).coefficient_tuple()
                          for i in I.basis()]) * X).change_ring(GF(p))

        # step 2: Find j such that if f=I[j], then mod 2 we have span(I[0],alpha*I[i])
        #         has trivial intersection with span(I[j],alpha*I[j]).
        #
        # In terms of our matrices alpha, beta, we can now think of I/p*I
        # as being the GF(p)^4 that M_alpha and M_beta naturally act on,
        # and I[0], I[1], I[2], I[3] correspond to the standard basis.
        #
        # We try each of the standard basis vectors.
        W0 = V.span([V.gen(0), V.gen(0) * M_alpha])
        assert W0.dimension() == 2
        j = None
        for i in range(1, 4):
            Wi = V.span([V.gen(i), V.gen(i) * M_alpha])
            if Wi.dimension() == 2 and W0.intersection(Wi).dimension() == 0:
                j = i
                break
        assert j is not None, "bug -- couldn't find basis"

        # step 3: Enumerate the elements of P^1(GF(p^2)), recording each
        #         cyclic submodule of degree p.
        answer = []
        f = V.gen(0)
        g = V.gen(j)
        M2_4 = MatrixSpace(GF(p), 4)
        M2_2 = MatrixSpace(QQ, 2, 4)
        Yp = p * Y
        from sage.algebras.quatalg.quaternion_algebra_cython import \
            rational_quaternions_from_integral_matrix_and_denom
        for v in [f + g * (a + b * M_alpha)
                  for a in GF(p) for b in GF(p)] + [g]:
            v0 = v
            v1 = v * M_alpha
            v2 = v * M_beta
            v3 = v1 * M_beta
            W = M2_4([v0, v1, v2, v3], coerce=False)
            if W.rank() == 2:
                gen_mat = Yp.stack(M2_2([v0.lift() * Y, v1.lift() * Y],
                                        coerce=False))
                gen_mat, d = gen_mat._clear_denom()
                H = gen_mat._hnf_pari(0, include_zero_rows=False)
                gens = tuple(rational_quaternions_from_integral_matrix_and_denom(A, H, d))
                answer.append(R.right_ideal(gens, check=False))
                if len(answer) == p + 1:
                    break
        return answer

    def hecke_matrix(self, n, algorithm='default', sparse=False, B=None):
        """
        Return the matrix of the `n`-th Hecke operator.

        INPUT:

        - `n` -- integer

        - ``algorithm`` -- string (default: 'default')

           - 'default' -- let Sage guess which algorithm is best

           - 'direct' -- use cyclic subideals (generally much
             better when you want few Hecke operators and the
             dimension is very large); uses 'theta' if n divides
             the level.

           - 'brandt' -- use Brandt matrices (generally much
             better when you want many Hecke operators and the
             dimension is very small; bad when the dimension
             is large)

        - ``sparse`` -- bool (default: ``False``)

        - `B` -- integer or ``None`` (default: ``None``); in direct
          algorithm, use theta series to this precision as an initial
          check for equality of ideal classes.

        EXAMPLES::

            sage: B = BrandtModule(3,7); B.hecke_matrix(2)
            [0 3]
            [1 2]
            sage: B.hecke_matrix(5, algorithm='brandt')
            [0 6]
            [2 4]
            sage: t = B.hecke_matrix(11, algorithm='brandt', sparse=True); t
            [ 6  6]
            [ 2 10]
            sage: type(t)
            <class 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
            sage: B.hecke_matrix(19, algorithm='direct', B=2)
            [ 8 12]
            [ 4 16]
        """
        n = ZZ(n)
        if n <= 0:
            raise IndexError("n must be positive.")
        if n not in self._hecke_matrices:
            if algorithm == 'default':
                try:
                    pr = len(self.__brandt_series_vectors[0][0])
                except (AttributeError, IndexError):
                    pr = 0
                if n <= pr:
                    # already trivially know the Hecke operator in this case
                    algorithm = 'brandt'
                if algorithm == 'default':  # still don't know
                    algorithm = 'direct'

            if self.level().gcd(n) != 1:
                algorithm = 'brandt'

            if algorithm == 'direct':
                T = self._compute_hecke_matrix(n, sparse=sparse, B=B)
            elif algorithm == 'brandt':
                T = self._compute_hecke_matrix_brandt(n, sparse=sparse)
            else:
                raise ValueError("unknown algorithm '%s'" % algorithm)
            T.set_immutable()
            self._hecke_matrices[n] = T
        return self._hecke_matrices[n]

    def _compute_hecke_matrix_prime(self, p, sparse=False, B=None):
        """
        Return matrix of the `p`-th Hecke operator on self.  The matrix
        is always computed using the direct algorithm.

        INPUT:

        - `p` -- prime number

        - `B` -- integer or None (default: None); in direct algorithm,
          use theta series to this precision as an initial check for
          equality of ideal classes.

        - ``sparse`` -- bool (default: False); whether matrix should be sparse

        EXAMPLES::

            sage: B = BrandtModule(37)
            sage: t = B._compute_hecke_matrix_prime(2); t
            [1 1 1]
            [1 0 2]
            [1 2 0]
            sage: type(t)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: type(B._compute_hecke_matrix_prime(2,sparse=True))
            <class 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
        """
        return self._compute_hecke_matrix_directly(n=p, B=B, sparse=sparse)

    def _compute_hecke_matrix_directly(self, n, B=None, sparse=False):
        """
        Given an integer `n` coprime to the level, return the matrix of
        the n-th Hecke operator on self, computed on our fixed basis
        by directly using the definition of the Hecke action in terms
        of fractional ideals.

        INPUT:

        - `n` -- integer, coprime to level

        - ``sparse`` -- bool (default: False); whether matrix should be sparse

        EXAMPLES::

            sage: B = BrandtModule(37)
            sage: t = B._compute_hecke_matrix_directly(2); t
            [1 1 1]
            [1 0 2]
            [1 2 0]
            sage: type(t)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: type(B._compute_hecke_matrix_directly(2,sparse=True))
            <class 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>

        You can't compute the Hecke operator for n not coprime to the level using this function::

            sage: B._compute_hecke_matrix_directly(37)
            Traceback (most recent call last):
            ...
            ValueError: n must be coprime to the level

        The generic function (which uses theta series) does work, though::

            sage: B.hecke_matrix(37)
            [1 0 0]
            [0 0 1]
            [0 1 0]

        An example where the Hecke operator isn't symmetric::

            sage: B = BrandtModule(43)
            sage: B._compute_hecke_matrix_directly(2)
            [1 2 0 0]
            [1 0 1 1]
            [0 1 0 2]
            [0 1 2 0]
            sage: B._compute_hecke_matrix_brandt(2)
            [1 2 0 0]
            [1 0 1 1]
            [0 1 0 2]
            [0 1 2 0]
        """
        level = self.level()
        if gcd(n, level) != 1:
            raise ValueError("n must be coprime to the level")

        # For rigor it does not matter at all what bound we chose.
        # This B is used only for the first phase of checking equality
        # of ideals modulo equivalence -- we always provably check
        # equivalence if the theta series are the same up to this
        # bound.
        if B is None:
            B = self.dimension() // 2 + 5

        T = matrix(self.base_ring(), self.dimension(), sparse=sparse)
        C = self.right_ideals()
        theta_dict = self._theta_dict(B)
        # I think the runtime of this algorithm is now dominated by
        # computing theta series of ideals.  The computation of
        # cyclic submodules is a lower order term.

        # TODO: temporary!! -- it's not sufficiently *optimized* to be
        # sure this is best in these cases.
        # d = lcm([a.denominator() for a in self.order_of_level_N().basis()])
        # q = self._smallest_good_prime()
        # if gcd(2*d*q,n) == 1:
        #     use_fast_alg = True
        # else:
        #     use_fast_alg = False

        use_fast_alg = False

        last_percent = 0
        for r in range(len(C)):
            percent_done = 100 * r // len(C)
            if percent_done != last_percent:
                if not(percent_done % 5):
                    verbose("percent done: %s" % percent_done)
                last_percent = percent_done
            if use_fast_alg:
                v = C[r].cyclic_right_subideals(n)
            else:
                v = self.cyclic_submodules(C[r], n)
            for J in v:
                J_theta = tuple(J.theta_series_vector(B))
                v = theta_dict[J_theta]
                if len(v) == 1:
                    T[r, v[0]] += 1
                else:
                    for i in v:
                        if C[i].is_equivalent(J, 0):
                            T[r, i] += 1
                            break
        return T

    @cached_method
    def _theta_dict(self, B):
        """
        Return a dictionary from theta series vectors of degree `B` to
        list of integers `i`, where the key is the vector of
        coefficients of the normalized theta series of the `i`th right
        ideal, as indexed by ``self.right_ideals()``.

        INPUT:

        - `B` -- positive integer, precision of theta series vectors

        OUTPUT:

        dictionary

        EXAMPLES:

        In this example the theta series determine the ideal classes::

            sage: B = BrandtModule(5,11); B
            Brandt module of dimension 4 of level 5*11 of weight 2 over Rational Field
            sage: sorted(list(B._theta_dict(5).items()))
            [((1, 0, 0, 4, 0), [3]),
             ((1, 0, 0, 4, 2), [2]),
             ((1, 0, 2, 0, 6), [1]),
             ((1, 2, 4, 0, 6), [0])]

        In this example, the theta series does not determine the ideal class::

             sage: sorted(list(BrandtModule(37)._theta_dict(6).items()))
             [((1, 0, 2, 2, 6, 4), [1, 2]), ((1, 2, 2, 4, 2, 4), [0])]
        """
        C = self.right_ideals()
        theta_dict = {}
        for i in range(len(C)):
            I_theta = tuple(C[i].theta_series_vector(B))
            if I_theta in theta_dict:
                theta_dict[I_theta].append(i)
            else:
                theta_dict[I_theta] = [i]
        return theta_dict

    def _compute_hecke_matrix_brandt(self, n, sparse=False):
        """
        Return the n-th Hecke matrix, computed using Brandt matrices
        (theta series).

        When the n-th Hecke operator is requested, we computed theta
        series to precision `2n+20`, since it only takes slightly
        longer, and this means that any Hecke operator `T_m` can
        quickly be computed, for `m<2n+20`.

        INPUT:

        - n -- integer, coprime to level
        - sparse -- bool (default: ``False``); whether matrix should be sparse

        EXAMPLES::

            sage: B = BrandtModule(3,17)
            sage: B._compute_hecke_matrix_brandt(3)
            [0 1 0 0]
            [1 0 0 0]
            [0 0 0 1]
            [0 0 1 0]
            sage: B._compute_hecke_matrix_brandt(5)
            [4 1 1 0]
            [1 4 0 1]
            [2 0 2 2]
            [0 2 2 2]
            sage: B._compute_hecke_matrix_brandt(5).fcp()
            (x - 6) * (x - 3) * (x^2 - 3*x - 2)

        """
        # we go out to 2*n+20 for efficiency, since it takes only a
        # little longer, but saves a lot of time if one computes
        # successive Hecke operators, which is a very common thing to
        # do.
        B = self._brandt_series_vectors()
        if len(B[0][0]) <= n:
            B = self._brandt_series_vectors(2 * n + 10)
        m = len(B)
        K = self.base_ring()
        return matrix(K, m, m, {(i, j): K(B[j][i][n])
                                for i in range(m)
                                for j in range(m)}, sparse=sparse)

    @cached_method
    def _smallest_good_prime(self):
        """
        Return the smallest prime number that does not divide the level.

        EXAMPLES::

            sage: BrandtModule(17,6)._smallest_good_prime()
            5
        """
        level = self.level()
        p = ZZ(2)
        while not(level % p):
            p = next_prime(p)
        return p

    @cached_method
    def right_ideals(self, B=None):
        """
        Return sorted tuple of representatives for the equivalence
        classes of right ideals in ``self``.

        OUTPUT:

        sorted tuple of fractional ideals

        EXAMPLES::

            sage: B = BrandtModule(23)
            sage: B.right_ideals()
            (Fractional ideal (2 + 2*j, 2*i + 2*k, 4*j, 4*k),
             Fractional ideal (2 + 2*j, 2*i + 6*k, 8*j, 8*k),
             Fractional ideal (2 + 10*j + 8*k, 2*i + 8*j + 6*k, 16*j, 16*k))

        TESTS::

            sage: B = BrandtModule(1009)
            sage: Is = B.right_ideals()
            sage: n = len(Is)
            sage: prod(not Is[i].is_equivalent(Is[j]) for i in range(n) for j in range(i))
            1
        """
        p = self._smallest_good_prime()
        R = self.order_of_level_N()
        I = R.unit_ideal()
        I = R.right_ideal([4 * x for x in I.basis()])

        if B is None:
            B = self.dimension() // 2 + 5

        ideals = [I]
        ideals_theta = {tuple(I.theta_series_vector(B)): [I]}
        new_ideals = [I]

        newly_computed_ideals = []
        got_something_new = True

        while got_something_new:
            got_something_new = False
            newly_computed_ideals = []
            for I in new_ideals:
                L = self.cyclic_submodules(I, p)
                for J in L:
                    is_new = True
                    J_theta = tuple(J.theta_series_vector(B))
                    if J_theta in ideals_theta:
                        for K in ideals_theta[J_theta]:
                            if J.is_equivalent(K, 0):
                                is_new = False
                                break
                    if is_new:
                        newly_computed_ideals.append(J)
                        ideals.append(J)
                        if J_theta in ideals_theta:
                            ideals_theta[J_theta].append(J)
                        else:
                            ideals_theta[J_theta] = [J]
                        verbose("found %s of %s ideals" % (len(ideals), self.dimension()), level=2)
                        if len(ideals) >= self.dimension():
                            ideals = tuple(sorted(ideals))
                            self.__right_ideals = ideals
                            return ideals
                        got_something_new = True
            new_ideals = list(newly_computed_ideals)

        return tuple(sorted(ideals))

    @cached_method
    def _ideal_products(self, diagonal_only=False):
        """
        Return all products of right ideals, which are used in computing
        the Brandt matrices.

        This function is used internally by the Brandt matrices
        algorithms.

        INPUT:

        - ``diagonal_only`` -- bool (default: ``False``) if ``True`` returns
          only the diagonal ideal products

        OUTPUT:

        list of ideals

        EXAMPLES::

            sage: B = BrandtModule(37)
            sage: B._ideal_products()
            [[Fractional ideal (8 + 8*j + 8*k, 4*i + 8*j + 4*k, 16*j, 16*k)],
             [Fractional ideal (8 + 24*j + 8*k, 4*i + 8*j + 4*k, 32*j, 32*k),
              Fractional ideal (16 + 16*j + 48*k, 4*i + 8*j + 36*k, 32*j + 32*k, 64*k)],
             [Fractional ideal (8 + 24*j + 24*k, 4*i + 24*j + 4*k, 32*j, 32*k),
              Fractional ideal (8 + 4*i + 16*j + 28*k, 8*i + 16*j + 8*k, 32*j, 64*k),
              Fractional ideal (16 + 16*j + 16*k, 4*i + 24*j + 4*k, 32*j + 32*k, 64*k)]]
            sage: B._ideal_products(diagonal_only=True)
            [Fractional ideal (8 + 8*j + 8*k, 4*i + 8*j + 4*k, 16*j, 16*k),
             Fractional ideal (16 + 16*j + 48*k, 4*i + 8*j + 36*k, 32*j + 32*k, 64*k),
             Fractional ideal (16 + 16*j + 16*k, 4*i + 24*j + 4*k, 32*j + 32*k, 64*k)]
        """
        L = self.right_ideals()
        n = len(L)
        if not n:
            return matrix(self.base_ring()[['q']], 0)

        # 1. Compute the diagonal
        D = [I.multiply_by_conjugate(I) for I in L]

        if diagonal_only:
            return D

        # 2. Compute the rest of the products
        P = []
        for i in range(n):
            v = [L[i].multiply_by_conjugate(L[j]) for j in range(i)]
            v.append(D[i])
            P.append(v)
        return P

    def _brandt_series_vectors(self, prec=None):
        """
        Return Brandt series coefficient vectors out to precision *at least* prec.

        EXAMPLES::

            sage: B = BrandtModule(37, use_cache=False)
            sage: B._brandt_series_vectors(5)
            [[(1/2, 1, 1, 2, 1), (1/2, 0, 1, 1, 3), (1/2, 0, 1, 1, 3)],
             [(1/2, 0, 1, 1, 3), (1/2, 1, 0, 0, 3), (1/2, 0, 2, 3, 1)],
             [(1/2, 0, 1, 1, 3), (1/2, 0, 2, 3, 1), (1/2, 1, 0, 0, 3)]]

        If you have computed to higher precision and ask for a lower
        precision, the higher precision is still returned::

            sage: B._brandt_series_vectors(2)
            [[(1/2, 1, 1, 2, 1), (1/2, 0, 1, 1, 3), (1/2, 0, 1, 1, 3)],
             [(1/2, 0, 1, 1, 3), (1/2, 1, 0, 0, 3), (1/2, 0, 2, 3, 1)],
             [(1/2, 0, 1, 1, 3), (1/2, 0, 2, 3, 1), (1/2, 1, 0, 0, 3)]]
        """
        if prec is None:
            try:
                return self.__brandt_series_vectors
            except AttributeError:
                prec = 2
        elif prec < 2:
            raise ValueError("prec must be at least 2")
        L = self.right_ideals()
        if not L:
            return [[]]
        try:
            if len(self.__brandt_series_vectors[0][0]) >= prec:
                return self.__brandt_series_vectors
        except AttributeError:
            pass

        n = len(L)
        # 1. Compute the theta series
        theta = [[I.theta_series_vector(prec) for I in x]
                 for x in self._ideal_products()]

        # 2. Compute the number e_j
        e = [theta[j][j][1] for j in range(n)]

        B = [[0 for _ in range(n)] for _ in range(n)]

        # 3. Make the Brandt matrix series
        for i in range(n):
            B[i][i] = theta[i][i] / e[i]
            for j in range(i):
                B[j][i] = theta[i][j] / e[j]
                B[i][j] = theta[i][j] / e[i]

        self.__brandt_series_vectors = B
        return B

    def brandt_series(self, prec, var='q'):
        r"""
        Return matrix of power series `\sum T_n q^n` to the given
        precision.

        Note that the Hecke operators in this series are
        always over `\QQ`, even if the base ring of this Brandt module
        is not `\QQ`.

        INPUT:

        - ``prec`` -- positive integer
        - ``var`` -- string (default: `q`)

        OUTPUT:

        matrix of power series with coefficients in `\QQ`

        EXAMPLES::

            sage: B = BrandtModule(11)
            sage: B.brandt_series(2)
            [1/4 + q + O(q^2)     1/4 + O(q^2)]
            [    1/6 + O(q^2) 1/6 + q + O(q^2)]
            sage: B.brandt_series(5)
            [1/4 + q + q^2 + 2*q^3 + 5*q^4 + O(q^5)   1/4 + 3*q^2 + 3*q^3 + 3*q^4 + O(q^5)]
            [  1/6 + 2*q^2 + 2*q^3 + 2*q^4 + O(q^5)         1/6 + q + q^3 + 4*q^4 + O(q^5)]


        Asking for a smaller precision works::

            sage: B.brandt_series(3)
            [1/4 + q + q^2 + O(q^3)   1/4 + 3*q^2 + O(q^3)]
            [  1/6 + 2*q^2 + O(q^3)       1/6 + q + O(q^3)]
            sage: B.brandt_series(3,var='t')
            [1/4 + t + t^2 + O(t^3)   1/4 + 3*t^2 + O(t^3)]
            [  1/6 + 2*t^2 + O(t^3)       1/6 + t + O(t^3)]
        """
        A = self._brandt_series_vectors(prec)
        R = QQ[[var]]
        n = len(A[0])
        return matrix(R, n, n,
                      [[R(x.list()[:prec], prec) for x in Y] for Y in A])

    @cached_method
    def eisenstein_subspace(self):
        """
        Return the 1-dimensional subspace of ``self`` on which the Hecke
        operators `T_p` act as `p+1` for `p` coprime to the level.

        .. NOTE::

            This function assumes that the base field has
            characteristic 0.

        EXAMPLES::

            sage: B = BrandtModule(11); B.eisenstein_subspace()
            Subspace of dimension 1 of Brandt module of dimension 2 of level 11 of weight 2 over Rational Field
            sage: B.eisenstein_subspace() is B.eisenstein_subspace()
            True
            sage: BrandtModule(3,11).eisenstein_subspace().basis()
            ((1, 1),)
            sage: BrandtModule(7,10).eisenstein_subspace().basis()
            ((1, 1, 1, 1/2, 1, 1, 1/2, 1, 1, 1),)
            sage: BrandtModule(7,10,base_ring=ZZ).eisenstein_subspace().basis()
            ((2, 2, 2, 1, 2, 2, 1, 2, 2, 2),)
        """
        if self.base_ring().characteristic():
            raise ValueError("characteristic must be 0")
        # cut down until we get a 1-d space using Hecke operators T_p
        # with p coprime to the level.
        V = self
        p = Integer(2)
        N = self.level()
        while V.dimension() >= 2:
            while not(N % p):
                p = p.next_prime()
            A = V.T(p) - (p + 1)
            V = A.kernel()
        return V

    def is_cuspidal(self):
        r"""
        Return whether ``self`` is cuspidal, i.e. has no Eisenstein part.

        EXAMPLES::

            sage: B = BrandtModule(3, 4)
            sage: B.is_cuspidal()
            False
            sage: B.eisenstein_subspace()
            Brandt module of dimension 1 of level 3*4 of weight 2 over Rational Field
        """
        return not self.eisenstein_subspace().dimension()

    @cached_method
    def monodromy_weights(self):
        r"""
        Return the weights for the monodromy pairing on this Brandt
        module.

        The weights are associated to each ideal class in our
        fixed choice of basis. The weight of an ideal class `[I]` is
        half the number of units of the right order `I`.

        NOTE: The base ring must be `\QQ` or `\ZZ`.

        EXAMPLES::

            sage: BrandtModule(11).monodromy_weights()
            (2, 3)
            sage: BrandtModule(37).monodromy_weights()
            (1, 1, 1)
            sage: BrandtModule(43).monodromy_weights()
            (2, 1, 1, 1)
            sage: BrandtModule(7,10).monodromy_weights()
            (1, 1, 1, 2, 1, 1, 2, 1, 1, 1)
            sage: BrandtModule(5,13).monodromy_weights()
            (1, 3, 1, 1, 1, 3)
            sage: BrandtModule(2).monodromy_weights()
            (12,)
            sage: BrandtModule(2,7).monodromy_weights()
            (3, 3)
        """
        # Before normalization,
        #
        #     theta(R) = 1 + e*q + ....
        #
        # where e is the number of units in the order R.
        #
        # Since the theta series may be normalized as
        #
        #     c * theta(R) = a[0] + a[1]*q + ...
        #
        # we recover e = a[1]/a[0] regardless of normalization.
        orders = self._ideal_products(diagonal_only=True)
        thetas = (R.theta_series_vector(2) for R in orders)
        return tuple(a[1] / a[0] / 2 for a in thetas)


#############################################################################
# Benchmarking
#############################################################################
def benchmark_magma(levels, silent=False):
    """
    INPUT:

    - ``levels`` -- list of pairs `(p,M)` where `p` is a prime not
      dividing `M`
    - ``silent`` -- bool, default ``False``; if ``True`` suppress
      printing during computation

    OUTPUT:

    list of 4-tuples ('magma', p, M, tm), where tm is the
    CPU time in seconds to compute T2 using Magma

    EXAMPLES::

        sage: a = sage.modular.quatalg.brandt.benchmark_magma([(11,1), (37,1), (43,1), (97,1)])  # optional - magma
        ('magma', 11, 1, ...)
        ('magma', 37, 1, ...)
        ('magma', 43, 1, ...)
        ('magma', 97, 1, ...)
        sage: a = sage.modular.quatalg.brandt.benchmark_magma([(11,2), (37,2), (43,2), (97,2)])  # optional - magma
        ('magma', 11, 2, ...)
        ('magma', 37, 2, ...)
        ('magma', 43, 2, ...)
        ('magma', 97, 2, ...)
    """
    ans = []
    from sage.interfaces.all import magma
    for p, M in levels:
        t = magma.cputime()
        magma.eval('HeckeOperator(BrandtModule(%s, %s),2)' % (p, M))
        tm = magma.cputime(t)
        v = ('magma', p, M, tm)
        if not silent:
            print(v)
        ans.append(v)
    return ans


def benchmark_sage(levels, silent=False):
    """
    INPUT:

    - ``levels`` -- list of pairs `(p,M)` where `p` is a prime
      not dividing `M`
    - ``silent`` -- bool, default ``False``; if ``True`` suppress
      printing during computation

    OUTPUT:

    list of 4-tuples ('sage', p, M, tm), where tm is the
    CPU time in seconds to compute T2 using Sage

    EXAMPLES::

        sage: a = sage.modular.quatalg.brandt.benchmark_sage([(11,1), (37,1), (43,1), (97,1)])
        ('sage', 11, 1, ...)
        ('sage', 37, 1, ...)
        ('sage', 43, 1, ...)
        ('sage', 97, 1, ...)
        sage: a = sage.modular.quatalg.brandt.benchmark_sage([(11,2), (37,2), (43,2), (97,2)])
        ('sage', 11, 2, ...)
        ('sage', 37, 2, ...)
        ('sage', 43, 2, ...)
        ('sage', 97, 2, ...)
    """
    from sage.misc.misc import cputime
    ans = []
    for p, M in levels:
        t = cputime()
        BrandtModule(p, M, use_cache=False).hecke_matrix(2)
        tm = cputime(t)
        v = ('sage', p, M, tm)
        if not silent:
            print(v)
        ans.append(v)
    return ans
