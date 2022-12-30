r"""
Key polynomials

Key polynomials (a.k.a. Demazure characters) are defined by applying an
operator `\pi_\sigma` to a monomial corresponding to an integer partition `\mu
\vdash n`. For a full definition, see
`SymmetricFunctions.com<https://www.symmetricfunctions.com/key.htm>`_. Key
polynomials are indexed by weak compositions with no trailing zeros, and
`\sigma` is the permutation of shortest length which sorts the indexing
composition into a partition.

.. SEEALSO::

    For Demazure characters in other types, see
    - :meth:`sage.combinat.root_system.weyl_characters.WeylCharacterRing.demazure_character`
    - :meth:`sage.categories.classical_crystals.ClassicalCrystals.ParentMethods.demazure_character`

AUTHORS:

- Trevor K. Karn (2022-08-17): initial version

"""

# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 (at) umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.composition import Composition
from sage.combinat.integer_vector import IntegerVector, IntegerVectors
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing, InfinitePolynomialRing_sparse
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial_sparse
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_commutative
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.rational_field import QQ
from sage.structure.element import parent

from collections.abc import Collection


class KeyPolynomial(CombinatorialFreeModule.Element):
    r"""
    Key polynomials are polynomials that form a basis for a polynomial ring
    and are indexed by weak compositions.

    Elements should be created by first creating the basis
    :class:`KeyPolynomialBasis` and passing a list representing the indexing
    composition.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
        sage: k = KeyPolynomialBasis(QQ)
        sage: f = k([4,3,2,1]) + k([1,2,3,4]); f
        k[1, 2, 3, 4] + k[4, 3, 2, 1]
        sage: f in k
        True
    """
    def _mul_(self, other):
        r"""
        Multiply the elements ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([4,3,2]) * k([1,1,1])
            k[5, 4, 3]

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k([4,3,2,0]) * k([1,1,1,0])
            k[5, 4, 3, 0]
        """
        return self.parent().from_polynomial(self.expand() * other.expand())

    def expand(self):
        r"""
        Return ``self`` written in the monomial basis.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: f = k([4,3,2,1])
            sage: f.expand()
            z_3*z_2^2*z_1^3*z_0^4

            sage: f = k([1,2,3])
            sage: f.expand()
            z_2^3*z_1^2*z_0 + z_2^3*z_1*z_0^2 + z_2^2*z_1^3*z_0
             + 2*z_2^2*z_1^2*z_0^2 + z_2^2*z_1*z_0^3 + z_2*z_1^3*z_0^2
             + z_2*z_1^2*z_0^3
        """
        R = self.parent()._polynomial_ring
        out = R.zero()
        z = self.parent().poly_gens()

        for m, c in self.monomial_coefficients().items():
            # get the partition itself
            mu = sorted(m, reverse=True)

            # create the monomial to apply
            monom = R.prod(z[i] ** mi for i, mi in enumerate(mu) if mi)

            # find the permutation sorting mu into m
            w = _sorting_word(m)

            out += c * _pi(self.parent(), w, monom)

        return out

    def pi(self, w):
        r"""
        Apply the operator `\pi_w` to ``self``.

        ``w`` may be either a ``Permutation`` or a list of indices of simple
        transpositions. In the case that it is a list of indices of simple
        transpositions, the convention follows that of the output of
        :meth:`sage.combinat.permutations.Permutation.reduced_word`
        namely that the transpositions are applied from left to right
        rather than from right to left.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([3,2,1]).pi(1)
            k[3, 1, 2]
            sage: k([3,2,1]).pi([1,0])
            k[1, 3, 2]
            sage: k([3,2,1]).pi(Permutation([3,2,1]))
            k[1, 2, 3]
            sage: f = k([3,2,1]) + k([3,2,1,1])
            sage: f.pi(1)
            k[3, 1, 2] + k[3, 1, 2, 1]

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k([3,2,1,0]).pi(1)
            k[3, 1, 2, 0]
            sage: k([3,2,1,0]).pi([1,0])
            k[1, 3, 2, 0]
            sage: k([3,2,1,0]).pi(Permutation([3,2,1,4]))
            k[1, 2, 3, 0]
            sage: f = k([3,2,1,0]) + k([3,2,1,1])
            sage: f.pi(1)
            k[3, 1, 2, 0] + k[3, 1, 2, 1]
        """
        P = self.parent()
        N = max(map(len, self.support()))
        if not isinstance(w, Permutation):
            if not isinstance(w, Collection):
                w = [w]
            from sage.combinat.permutation import Permutations
            w = Permutations(N).from_reduced_word([i+1 for i in w])

        total = P.zero()

        if P._k:
            for m, c in self.monomial_coefficients().items():
                if (n := len(m)) < N:
                    m = list(m) + [0]*(N-n)
                total += P._from_dict({P._indices(w.action(m)): c})
        else:
            for m, c in self.monomial_coefficients().items():
                if (n := len(m)) < N:
                    m = list(m) + [0]*(N-n)
                total += P._from_dict({P._indices(w.action(m)).trim(): c})

        return total

    def divided_difference(self, w):
        r"""
        Apply the divided difference operator `\partial_w` to ``self``.

        The convention is to apply from left to right so if
        ``w = [w0, w1, ..., wn]`` then we apply
        `\partial_{w_1 w_2 \cdots w_n} \circ \partial_{w_0}`

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([3,2,1]).divided_difference(1)
            k[3, 1, 1]
            sage: k([3,2,1]).divided_difference([1,2])
            k[3, 1]

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k([3,2,1,0]).divided_difference(1)
            k[3, 1, 1, 0]

        """
        if not isinstance(w, Collection):
            w = [w]
        P = self.parent()
        f = self.expand()
        for wi in w:
            f = _divided_difference(P, wi, f)
        return P.from_polynomial(f)

class KeyPolynomialBasis(CombinatorialFreeModule):
    r"""
    The key polynomial basis for a polynomial ring.

    EXAMPLES:

    Key polynomials are a basis, indexed by (weak) compositions,
    for polynomial rings::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
        sage: k = KeyPolynomialBasis(QQ)
        sage: k([3,0,1,2])
        k[3, 0, 1, 2]
        sage: R = k.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field

    We can expand them in the standard monomial basis::

        sage: k([3,0,1,2]).expand()
        z_3^2*z_2*z_0^3 + z_3^2*z_1*z_0^3 + z_3*z_2^2*z_0^3
         + 2*z_3*z_2*z_1*z_0^3 + z_3*z_1^2*z_0^3 + z_2^2*z_1*z_0^3
         + z_2*z_1^2*z_0^3

        sage: k([0,0,2]).expand()
        z_2^2 + z_2*z_1 + z_2*z_0 + z_1^2 + z_1*z_0 + z_0^2

    If we have a polynomial, we can express it in the key basis::

        sage: z = R.gen()
        sage: k.from_polynomial(z[2]^2*z[1]*z[0])
        k[1, 1, 2] - k[1, 2, 1]

        sage: f = z[3]^2*z[2]*z[0]^3 + z[3]^2*z[1]*z[0]^3 + z[3]*z[2]^2*z[0]^3 + \
        ....: 2*z[3]*z[2]*z[1]*z[0]^3 + z[3]*z[1]^2*z[0]^3 + z[2]^2*z[1]*z[0]^3 + \
        ....: z[2]*z[1]^2*z[0]^3
        sage: k.from_polynomial(f)
        k[3, 0, 1, 2]

    Since the ring of key polynomials may be regarded as a different choice of
    basis for a polynomial ring, it forms an algebra, so we have
    multiplication::

        sage: k([10,5,2])*k([1,1,1])
        k[11, 6, 3]

    We can also multiply by polynomials in the monomial basis::

        sage: k([10,9,1])*z[0]
        k[11, 9, 1]
        sage: z[0] * k([10,9,1])
        k[11, 9, 1]
        sage: k([10,9,1])*(z[0] + z[3])
        k[10, 9, 1, 1] + k[11, 9, 1]

    When the sorting permutation is the longest element, the key polynomial
    agrees with the Schur polynomial::

        sage: s = SymmetricFunctions(QQ).schur()
        sage: k([1,2,3]).expand()
        z_2^3*z_1^2*z_0 + z_2^3*z_1*z_0^2 + z_2^2*z_1^3*z_0
         + 2*z_2^2*z_1^2*z_0^2 + z_2^2*z_1*z_0^3 + z_2*z_1^3*z_0^2
         + z_2*z_1^2*z_0^3
        sage: s[3,2,1].expand(3)
        x0^3*x1^2*x2 + x0^2*x1^3*x2 + x0^3*x1*x2^2 + 2*x0^2*x1^2*x2^2
         + x0*x1^3*x2^2 + x0^2*x1*x2^3 + x0*x1^2*x2^3

    The polynomial expansions can be computed using crystals and expressed in
    terms of the key basis::

        sage: T = crystals.Tableaux(['A',3],shape=[2,1])
        sage: f = T.demazure_character([3,2,1])
        sage: k.from_polynomial(f)
        k[1, 0, 0, 2]
    
    The default behavior is to work in a polynomial ring with infinitely many
    variables. One can work in a specicfied number of variables::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
        sage: k = KeyPolynomialBasis(QQ, 4)
        sage: k([3,0,1,2]).expand()
        z_0^3*z_1^2*z_2 + z_0^3*z_1*z_2^2 + z_0^3*z_1^2*z_3 
         + 2*z_0^3*z_1*z_2*z_3 + z_0^3*z_2^2*z_3 + z_0^3*z_1*z_3^2 + z_0^3*z_2*z_3^2

        sage: k([0,0,2,0]).expand()
        z_0^2 + z_0*z_1 + z_1^2 + z_0*z_2  + z_1*z_2 + z_2^2

        sage: k([0,0,2,0]).expand().parent()
        Multivariate Polynomial Ring in z_0, z_1, z_2, z_3 over Rational Field

    If working in a specified number of variables, the length of the indexing
    composition must be the same as the number of variables::

        sage: k([0,0,2])
        Traceback (most recent call last):
         ...
        TypeError: do not know how to make x (= [0, 0, 2]) an element of self
         (=Key polynomial basis over Rational Field)

    One can also work in a specified polynomial ring::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
        sage: k = KeyPolynomialBasis(QQ['x0', 'x1', 'x2', 'x3'])
        sage: k([0,2,0,0])
        k[0, 2, 0, 0]
        sage: k([4,0,0,0]).expand()
        x0^4

    If one wishes to use a polynomial ring as coefficients for the key
    polynomials, pass the keyword argument ``poly_coeffs=True``::

        sage: k = KeyPolynomialBasis(QQ['q'], poly_coeffs=True)
        sage: R = k.base_ring(); R
        Univariate Polynomial Ring in q over Rational Field
        sage: R.inject_variables()
        Defining q
        sage: (q^2 + q + 1)*k([0,2,2,0,3,2])
        (q^2+q+1)*k[0, 2, 2, 0, 3, 2]
    """

    Element = KeyPolynomial

    @staticmethod
    def __classcall_private__(cls, R=None, k=None, poly_ring=None, poly_coeffs=False):
        r"""
        Normalize input.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: KeyPolynomialBasis(InfinitePolynomialRing(QQ, ['x', 'y']))
            Traceback (most recent call last):
            ...
            ValueError: Polynomial ring has too many generators

            sage: KeyPolynomialBasis(QQ['t0','t1','t2','t3'])
            Key polynomial basis over Rational Field
            
            sage: KeyPolynomialBasis(QQ['t'])
            Key polynomial basis over Rational Field
            
            sage: KeyPolynomialBasis(InfinitePolynomialRing(QQ['t'], 'z'))
            Key polynomial basis over Univariate Polynomial Ring in t over Rational Field

            sage: KeyPolynomialBasis(QQ)
            Key polynomial basis over Rational Field

            sage: KeyPolynomialBasis(QQ, 3)
            Key polynomial basis over Rational Field
        """

        poly_type = (PolynomialRing_commutative, 
                     MPolynomialRing_base,
                     InfinitePolynomialRing_sparse)

        if isinstance(R, poly_type):
            # if a polynomial ring is provided, we need to determine
            # if it is meant to be self.polynomial_ring() or self.base_ring()
            if isinstance(R, poly_type[0:2]):
                k = R.ngens()
            if isinstance(R, InfinitePolynomialRing_sparse) and R.ngens() > 1:
                raise ValueError(f"Polynomial ring has too many generators")
            if isinstance(R.base_ring(), poly_type[0:2]):
                # if R is of the form K[t_1, ..., t_n][z_*]
                # or K[t_1, ..., t_n][z_1, ..., z_k]
                return cls.__classcall__(cls, k=k, poly_ring=R)
            if poly_coeffs:
                # if R is a polynomial ring, but its base ring is not
                # and poly_coeffs is true, then we should interpret
                # R as the base ring
                return cls.__classcall__(cls, R=R)
            return cls.__classcall__(cls, k=k, poly_ring=R)
        else:
            # if R is not a polynomial ring, we know it is self.base_ring()
            return cls.__classcall__(cls, R=R, k=k)

    def __init__(self, R=None, k=None, poly_ring=None, poly_coeffs=False):
        """
        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k == loads(dumps(k))
            True
            sage: TestSuite(k).run()

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k == loads(dumps(k))
            True
            sage: TestSuite(k).run()
        """
        self._repr_option_bracket = False
        self._k = k

        if self._k:
            def build_index(m):
                return self._indices(m)
        else:
            def build_index(m):
                return self._indices(reversed(m)).trim()

        self._build_index = build_index

        if R:
            if poly_ring:
                raise ValueError("Specify only one of base_ring or poly_ring (not both)")
            if k:
                self._polynomial_ring = PolynomialRing(R, 'z_', k)
            else:
                self._polynomial_ring = InfinitePolynomialRing(R, 'z')
        if poly_ring:
            if R:
                raise ValueError("Specify only one of base_ring or poly_ring (not both)")
            R = poly_ring.base_ring()
            self._polynomial_ring = poly_ring

        self._name = f"Key polynomial basis"

        CombinatorialFreeModule.__init__(self, R, IntegerVectors(k=k),
                                         category=GradedAlgebrasWithBasis(R),
                                         prefix='k')

    def _coerce_map_from_(self, R):
        r"""
        EXAMPLES::
            
            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ) 
            sage: m1 = k([3,2,4,0]); m1
            k[3, 2, 4]
            sage: m2 = k(Composition([3, 2, 4])); m2
            k[3, 2, 4]
            sage: m1 == m2
            True

            sage: R = k.polynomial_ring()
            sage: z = R.gen()
            sage: z[0] * k([4, 3, 3, 2])
            k[5, 3, 3, 2]
        """
        P = self._polynomial_ring
        from sage.structure.coerce_maps import CallableConvertMap
        if R is P:
            return CallableConvertMap(R, self, self.from_polynomial)
        phi = P.coerce_map_from(R)
        if phi is not None: 
            return self.coerce_map_from(P) * phi
        return None

    def _monomial(self, x):
        r"""
        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k([3, 2, 3, 4])
            k[3, 2, 3, 4]
            sage: k = KeyPolynomialBasis(QQ, 5)
            sage: k([3, 2, 3, 4, 0])
            k[3, 2, 3, 4, 0]
        """
        if self._k:
            return self._from_dict({x: self.base_ring().one()}, remove_zeros=False)
        return self._from_dict({x.trim(): self.base_ring().one()}, remove_zeros=False)

    def one_basis(self):
        r"""
        Return the basis element indexing the identity.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k.one_basis()
            []

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k.one_basis()
            [0, 0, 0, 0]            
        """
        if self._k:
            return self._indices([0] * self._k)
        return self._indices([])

    def polynomial_ring(self):
        r"""
        Return the polynomial ring associated to ``self``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k.polynomial_ring()
            Infinite polynomial ring in z over Rational Field

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k.polynomial_ring()
            Multivariate Polynomial Ring in z_0, z_1, z_2, z_3 over Rational Field
        """
        return self._polynomial_ring

    def poly_gens(self):
        r"""
        Return the polynomial generators for the polynomial ring
        associated to ``self``.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: k.poly_gens()
            z_*

            sage: k = KeyPolynomialBasis(QQ, 4)
            sage: k.poly_gens()
            (z_0, z_1, z_2, z_3)
        """
        if self._k:
            return self._polynomial_ring.gens()
        return self._polynomial_ring.gen()

    def from_polynomial(self, f):
        r"""
        Expand a polynomial in terms of the key basis.

        EXAMPLES::

            sage: from sage.combinat.key_polynomial import KeyPolynomialBasis
            sage: k = KeyPolynomialBasis(QQ)
            sage: z = k.poly_gens(); z
            z_*
            sage: p = z[0]^4*z[1]^2*z[2]*z[3] + z[0]^4*z[1]*z[2]^2*z[3]
            sage: k.from_polynomial(p)
            k[4, 1, 2, 1]

            sage: all(k(c) == k.from_polynomial(k(c).expand()) for c in IntegerVectors(n=5, k=4))
            True

            sage: T = crystals.Tableaux(['A', 4], shape=[4,2,1,1])
            sage: k.from_polynomial(T.demazure_character([2]))
            k[4, 1, 2, 1]

        """
        if f not in self._polynomial_ring:
            try:  # to accept elements of SymbolicRing
                from sage.calculus.var import var
                f = f.substitute(list(d == var(f'z_{i}')
                                 for i, d in enumerate(f.variables())))
                f = self._polynomial_ring(f)
            except AttributeError:
                raise ValueError(f"f must be an element of {self._polynomial_ring}")

        out = self.zero()
        counter = 0

        while f:
            M = f.monomials()[0]
            c = f.monomial_coefficient(M)

            new_term = self._from_dict({self._build_index(*M.exponents()): c})

            f -= new_term.expand()
            out += new_term

        return out


def _divided_difference(P, i, f):
    r"""
    Apply the `i`th divided difference operator to the polynomial `f`.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _divided_difference
        sage: k = KeyPolynomialBasis(QQ)
        sage: z = k.poly_gens()
        sage: f = z[1]*z[2]^3 + z[1]*z[2]*z[3]
        sage: _divided_difference(k, 2, f)
        z_3^2*z_1 + z_3*z_2*z_1 + z_2^2*z_1

        sage: k = KeyPolynomialBasis(QQ, 4)
        sage: z = k.poly_gens()
        sage: f = z[1]*z[2]^3 + z[1]*z[2]*z[3]
        sage: _divided_difference(k, 2, f)
        z_1*z_2^2 + z_1*z_2*z_3 + z_1*z_3^2

        sage: k = KeyPolynomialBasis(QQ)
        sage: R = k.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field
        sage: z = R.gen()
        sage: _divided_difference(k, 1, z[1]*z[2]^3)
        -z_2^2*z_1 - z_2*z_1^2
        sage: _divided_difference(k, 2, z[1]*z[2]*z[3])
        0
        sage: _divided_difference(k, 3, z[1]*z[2]*z[3])
        z_2*z_1
        sage: _divided_difference(k, 3, z[1]*z[2]*z[4])
        -z_2*z_1

        sage: k = KeyPolynomialBasis(QQ, 5)
        sage: z = k.polynomial_ring().gens()
        sage: _divided_difference(k, 1, z[1]*z[2]^3)
        -z_1^2*z_2 - z_1*z_2^2
        sage: _divided_difference(k, 2, z[1]*z[2]*z[3])
        0
        sage: _divided_difference(k, 3, z[1]*z[2]*z[3])
        z_1*z_2
        sage: _divided_difference(k, 3, z[1]*z[2]*z[4])
        -z_1*z_2
    """
    R = P.polynomial_ring()
    z = P.poly_gens()

    si_f = f.subs({z[i+1]:z[i], z[i]:z[i+1]})

    return (si_f - f)//(z[i+1] - z[i])

def _pi(P, w, f):
    r"""
    Apply the operator `\pi_w` to the polynomial `f`.
    ``w`` may be either a single index or a list of
    indices of simple transpositions.

    .. WARNING::

        The simple transpositions should be applied from left to right. 

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _pi
        sage: k = KeyPolynomialBasis(QQ)
        sage: R = k.polynomial_ring()
        sage: z, = R.gens()
        sage: w = [2, 3]
        sage: _pi(k, w, z[1]^4*z[2]^2*z[3]*z[4])
        z_4^2*z_3*z_2*z_1^4 + z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: w = [3, 2]
        sage: _pi(k, w, z[1]^4*z[2]^2*z[3]*z[4])
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: w = [2, 1]
        sage: _pi(k, w, z[1]^2*z[2])
        z_3*z_2^2 + z_3*z_2*z_1 + z_3*z_1^2 + z_2^2*z_1 + z_2*z_1^2
    """
    if not hasattr(w, "__iter__"):  # this allows us to pass i instead of a word
        return _pi_i(P, w, f)
    for i in w:
        f = _pi_i(P, i, f)
    return f

def _pi_i(P, i, f):
    r"""
    Apply `\pi_i` for a single simple transposition `s_i = (i, i+1)`.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import KeyPolynomialBasis, _pi_i
        sage: k = KeyPolynomialBasis(QQ)
        sage: z = k.poly_gens()
        sage: _pi_i(k, 3, z[1]^4*z[2]^2*z[4])
        0

        sage: _pi_i(k, 2, z[1]^4*z[2]^2*z[3]*z[4])
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4
    """
    R = P.polynomial_ring()
    z = P.poly_gens()
    return _divided_difference(P, i, z[i] * f)

def _sorting_word(alpha):
    r"""
    Get a reduced word for the permutation which sorts ``alpha``
    into a partition. 

    The result is a list ``l = [i0, i1, i2, ...]`` where each ``ij``
    is a nonnegative integer such that it applies the simple
    transposition `(i_j, i_j+1)`. The transpositions are applied
    starting with ``i0``, then ``i1`` is applied, followed by ``i2``,
    and so on. See :meth:`sage.combinat.permutation.Permutation.reduced_words`
    for the convention used.

    EXAMPLES::

        sage: IV = IntegerVectors()
        sage: from sage.combinat.key_polynomial import _sorting_word
        sage: list(_sorting_word(IV([2,3,2])))
        [0]
        sage: list(_sorting_word(IV([5,6,7])))
        [0, 1, 0]
        sage: list(_sorting_word(IV([0,3,2])))
        [1, 0]
        sage: list(_sorting_word(IV([0,3,0,2])))
        [1, 2, 0]
        sage: list(_sorting_word(IV([3,2,1])))
        []
        sage: list(_sorting_word(IV([2,3,3])))
        [1, 0]
    """
    w = []
    n = len(alpha)

    # bubble sort to get the shortest sorting word
    with alpha.clone() as ac:
        for i in range(n-1):
            for j in range(n-i-1):
                if ac[j] < ac[j + 1]:
                    w.append(j)
                    ac[j], ac[j + 1] = ac[j + 1], ac[j]
    return reversed(w)
