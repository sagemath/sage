# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Key polynomials

Key polynomials (also known as type A Demazure characters) are defined by
applying the divided difference operator `\pi_\sigma`, where `\sigma` is
a permutation, to a monomial corresponding to an integer partition
`\mu \vdash n`.

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
from sage.misc.cachefunc import cached_method
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation
from sage.structure.element import parent
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing, InfinitePolynomialRing_sparse
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_commutative
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base

from collections.abc import Collection


class KeyPolynomial(CombinatorialFreeModule.Element):
    r"""
    A key polynomial.

    Key polynomials are polynomials that form a basis for a polynomial ring
    and are indexed by weak compositions.

    Elements should be created by first creating the basis
    :class:`KeyPolynomialBasis` and passing a list representing the indexing
    composition.

    EXAMPLES::

        sage: k = KeyPolynomials(QQ)
        sage: f = k([4,3,2,1]) + k([1,2,3,4]); f
        k[1, 2, 3, 4] + k[4, 3, 2, 1]
        sage: f in k
        True
    """
    def _mul_(self, other):
        r"""
        Multiply the elements ``self`` and ``other``.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k([4,3,2]) * k([1,1,1])
            k[5, 4, 3]

            sage: k = KeyPolynomials(QQ, 4)
            sage: k([4,3,2,0]) * k([1,1,1,0])
            k[5, 4, 3, 0]
        """
        return self.parent().from_polynomial(self.expand() * other.expand())

    def expand(self):
        r"""
        Return ``self`` written in the monomial basis (i.e., as an element
        in the corresponding polynomial ring).

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: f = k([4,3,2,1])
            sage: f.expand()
            z_3*z_2^2*z_1^3*z_0^4

            sage: f = k([1,2,3])
            sage: f.expand()
            z_2^3*z_1^2*z_0 + z_2^3*z_1*z_0^2 + z_2^2*z_1^3*z_0
             + 2*z_2^2*z_1^2*z_0^2 + z_2^2*z_1*z_0^3 + z_2*z_1^3*z_0^2
             + z_2*z_1^2*z_0^3
        """
        P = self.parent()
        R = P._polynomial_ring
        out = R.zero()
        z = P.poly_gens()

        for m, c in self.monomial_coefficients().items():
            # find the permutation sorting mu into m
            w, mu = sorting_word(m)

            # create the monomial to apply
            monom = R.prod(z[i] ** mi for i, mi in enumerate(mu) if mi)

            out += c * isobaric_divided_difference(monom, w)

        return out

    to_polynomial = expand

    def pi(self, w):
        r"""
        Apply the operator `\pi_w` to ``self``.

        ``w`` may be either a ``Permutation`` or a list of indices of simple
        transpositions (1-based).

        The convention is to apply from left to right so if
        ``w = [w1, w2, ..., wm]`` then we apply
        `\pi_{w_2 \cdots w_m} \circ \pi_{w_1}`

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k([3,2,1]).pi(2)
            k[3, 1, 2]
            sage: k([3,2,1]).pi([2,1])
            k[1, 3, 2]
            sage: k([3,2,1]).pi(Permutation([3,2,1]))
            k[1, 2, 3]
            sage: f = k([3,2,1]) + k([3,2,1,1])
            sage: f.pi(2)
            k[3, 1, 2] + k[3, 1, 2, 1]
            sage: k.one().pi(1)
            k[]

            sage: k([3,2,1,0]).pi(2).pi(2)
            k[3, 1, 2]
            sage: (-k([3,2,1,0]) + 4*k([3,1,2,0])).pi(2)
            3*k[3, 1, 2]

            sage: k = KeyPolynomials(QQ, 4)
            sage: k([3,2,1,0]).pi(2)
            k[3, 1, 2, 0]
            sage: k([3,2,1,0]).pi([2,1])
            k[1, 3, 2, 0]
            sage: k([3,2,1,0]).pi(Permutation([3,2,1,4]))
            k[1, 2, 3, 0]
            sage: f = k([3,2,1,0]) + k([3,2,1,1])
            sage: f.pi(2)
            k[3, 1, 2, 0] + k[3, 1, 2, 1]
            sage: k.one().pi(1)
            k[0, 0, 0, 0]

        TESTS:

        We check that this is consistent with the definition via the
        isobaric divided difference oerators::

            sage: from sage.combinat.key_polynomial import isobaric_divided_difference as idd
            sage: k = KeyPolynomials(QQ, 4)
            sage: S4 = Permutations(4)
            sage: f = k([4,2,2,0])
            sage: all(idd(f.expand(), w.reduced_word()) == f.pi(w).expand() for w in S4)
            True

            sage: f = k([4,2,0,1]) - 3 * k([2,0,1,2])
            sage: all(idd(f.expand(), w.reduced_word()) == f.pi(w).expand() for w in S4)
            True
        """
        P = self.parent()
        if isinstance(w, Permutation):
            w = w.reduced_word()
        if not isinstance(w, Collection):
            w = [w]

        if not w or not self:
            return self

        N = max(w) + 1

        if P._k is not None and N > P._k:
            raise ValueError(f"pi_{N-1} does not exist for this polynomial ring")

        ret = P.element_class(P, {})
        for m, c in self._monomial_coefficients.items():
            m = list(m)
            n = len(m)
            for i in w:
                if i > n:
                    continue
                if i == n:
                    m += [0]
                    n += 1
                if m[i-1] <= m[i]:
                    continue
                m[i-1], m[i] = m[i], m[i-1]
            m = P._indices(m)
            if P._k is None:
                m = m.trim()
            if m in ret._monomial_coefficients:
                ret._monomial_coefficients[m] += c
            else:
                ret._monomial_coefficients[m] = c
            if not ret._monomial_coefficients[m]:
                del ret._monomial_coefficients
        return ret

    isobaric_divided_difference = pi

    def divided_difference(self, w):
        r"""
        Apply the divided difference operator `\partial_w` to ``self``.

        The convention is to apply from left to right so if
        ``w = [w1, w2, ..., wm]`` then we apply
        `\partial_{w_2 \cdots w_m} \circ \partial_{w_1}`

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k([3,2,1]).divided_difference(2)
            k[3, 1, 1]
            sage: k([3,2,1]).divided_difference([2,3])
            k[3, 1]

            sage: k = KeyPolynomials(QQ, 4)
            sage: k([3,2,1,0]).divided_difference(2)
            k[3, 1, 1, 0]
        """
        if not isinstance(w, Collection):
            w = [w]
        f = self.expand()
        for wi in w:
            f = divided_difference(f, wi)
        return self.parent().from_polynomial(f)


class KeyPolynomialBasis(CombinatorialFreeModule):
    r"""
    The key polynomial basis for a polynomial ring.

    For a full definition, see
    `SymmetricFunctions.com <https://www.symmetricfunctions.com/key.htm>`_.
    Key polynomials are indexed by weak compositions with no trailing zeros,
    and `\sigma` is the permutation of shortest length which sorts the
    indexing composition into a partition.

    EXAMPLES:

    Key polynomials are a basis, indexed by (weak) compositions,
    for polynomial rings::

        sage: k = KeyPolynomials(QQ)
        sage: k([3,0,1,2])
        k[3, 0, 1, 2]
        sage: k([3,0,1,2])/2
        1/2*k[3, 0, 1, 2]
        sage: R = k.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field

        sage: K = KeyPolynomials(GF(5)); K
        Key polynomial basis over Finite Field of size 5
        sage: 2*K([3,0,1,2])
        2*k[3, 0, 1, 2]
        sage: 5*(K([3,0,1,2]) + K([3,1,1]))
        0

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

        sage: k = KeyPolynomials(QQ, 4)
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

        sage: k = KeyPolynomials(QQ['x0', 'x1', 'x2', 'x3'])
        sage: k([0,2,0,0])
        k[0, 2, 0, 0]
        sage: k([4,0,0,0]).expand()
        x0^4

    If one wishes to use a polynomial ring as coefficients for the key
    polynomials, pass the keyword argument ``poly_coeffs=True``::

        sage: k = KeyPolynomials(QQ['q'], poly_coeffs=True)
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

            sage: KeyPolynomials(InfinitePolynomialRing(QQ, ['x', 'y']))
            Traceback (most recent call last):
            ...
            ValueError: polynomial ring has too many generators

            sage: KeyPolynomials(QQ['t0','t1','t2','t3'])
            Key polynomial basis over Rational Field

            sage: KeyPolynomials(QQ['t'])
            Key polynomial basis over Rational Field

            sage: KeyPolynomials(InfinitePolynomialRing(QQ['t'], 'z'))
            Key polynomial basis over Univariate Polynomial Ring in t over Rational Field

            sage: KeyPolynomials(QQ)
            Key polynomial basis over Rational Field

            sage: KeyPolynomials(QQ, 3)
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
                raise ValueError("polynomial ring has too many generators")
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

    def __init__(self, R=None, k=None, poly_ring=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: R = GF(3)['t'].fraction_field()
            sage: k = KeyPolynomials(QQ)
            sage: TestSuite(k).run()
            sage: k = KeyPolynomials(R)
            sage: TestSuite(k).run()

            sage: k = KeyPolynomials(QQ, 4)
            sage: TestSuite(k).run()
            sage: k = KeyPolynomials(R, 4)
            sage: TestSuite(k).run()
        """
        self._k = k

        if self._k is not None:
            def build_index(m):
                return self._indices(m)
        else:
            def build_index(m):
                return self._indices(reversed(m)).trim()

        self._build_index = build_index

        if R is not None:
            if poly_ring:
                raise ValueError("specify only one of base_ring or poly_ring (not both)")
            if k:
                self._polynomial_ring = PolynomialRing(R, 'z_', k)
            else:
                self._polynomial_ring = InfinitePolynomialRing(R, 'z')
        if poly_ring is not None:
            if R is not None:
                raise ValueError("specify only one of base_ring or poly_ring (not both)")
            R = poly_ring.base_ring()
            self._polynomial_ring = poly_ring

        self._name = "Key polynomial basis"

        CombinatorialFreeModule.__init__(self, R, IntegerVectors(k=k),
                                         category=GradedAlgebrasWithBasis(R),
                                         prefix='k', bracket=False)

    def _coerce_map_from_(self, R):
        r"""
        Return the coercion map from ``R`` if it exists.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: m1 = k([3, 2, 4, 0]); m1
            k[3, 2, 4]
            sage: m2 = k(Composition([3, 2, 4])); m2
            k[3, 2, 4]
            sage: m1 == m2
            True

            sage: R = k.polynomial_ring()
            sage: z = R.gen()
            sage: z[0] * k([4, 3, 3, 2])
            k[5, 3, 3, 2]

            sage: X = SchubertPolynomialRing(QQ)
            sage: k(X([4, 3, 2, 1]))
            k[3, 2, 1]
        """
        P = self._polynomial_ring
        if R is P:
            return self.from_polynomial

        from sage.combinat.schubert_polynomial import SchubertPolynomialRing_xbasis
        if isinstance(R, SchubertPolynomialRing_xbasis):
            return self.from_schubert_polynomial

        phi = P.coerce_map_from(R)
        if phi is not None:
            return self.coerce_map_from(P) * phi
        return None

    def _monomial(self, x):
        r"""
        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k([3, 2, 3, 4, 0])
            k[3, 2, 3, 4]
            sage: k = KeyPolynomials(QQ, 5)
            sage: k([3, 2, 3, 4, 0])
            k[3, 2, 3, 4, 0]
        """
        if self._k:
            return self._from_dict({x: self.base_ring().one()}, remove_zeros=False)
        return self._from_dict({x.trim(): self.base_ring().one()}, remove_zeros=False)

    def __getitem__(self, c):
        """
        This method implements the abuses of notations ``k[2,1]``,
        ``k[[2,1]]``, etc.

        INPUT:

        - ``c`` -- anything that can represent an index of a basis element

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k[3]
            k[3]
            sage: k[3, 0, 2]
            k[3, 0, 2]
            sage: k[3, 1, 2, 0, 0]
            k[3, 1, 2]

            sage: k = KeyPolynomials(QQ, 4)
            sage: k[3, 0, 1, 0]
            k[3, 0, 1, 0]
            sage: k[3]
            Traceback (most recent call last):
            ...
            ValueError: [3] doesn't satisfy correct constraints
        """
        C = self._indices
        if not isinstance(c, C.element_class):
            if c in ZZ:
                c = C([c])
            else:
                c = C(c)
        return self._monomial(c)

    @cached_method
    def one_basis(self):
        r"""
        Return the basis element indexing the identity.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k.one_basis()
            []

            sage: k = KeyPolynomials(QQ, 4)
            sage: k.one_basis()
            [0, 0, 0, 0]
        """
        if self._k:
            return self._indices([0] * self._k)
        return self._indices([])

    def degree_on_basis(self, alpha):
        """
        Return the degree of the basis element indexed by ``alpha``.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k.degree_on_basis([2,1,0,2])
            5

            sage: k = KeyPolynomials(QQ, 5)
            sage: k.degree_on_basis([2,1,0,2,0])
            5
        """
        return ZZ(sum(alpha))

    def polynomial_ring(self):
        r"""
        Return the polynomial ring associated to ``self``.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k.polynomial_ring()
            Infinite polynomial ring in z over Rational Field

            sage: k = KeyPolynomials(QQ, 4)
            sage: k.polynomial_ring()
            Multivariate Polynomial Ring in z_0, z_1, z_2, z_3 over Rational Field
        """
        return self._polynomial_ring

    def poly_gens(self):
        r"""
        Return the polynomial generators for the polynomial ring
        associated to ``self``.

        EXAMPLES::

            sage: k = KeyPolynomials(QQ)
            sage: k.poly_gens()
            z_*

            sage: k = KeyPolynomials(QQ, 4)
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

            sage: k = KeyPolynomials(QQ)
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

        while f:
            M = f.monomials()[0]
            c = f.monomial_coefficient(M)

            new_term = self._from_dict({self._build_index(*M.exponents()): c})

            f -= new_term.expand()
            out += new_term

        return out

    def from_schubert_polynomial(self, x):
        r"""
        Expand a Schubert polynomial in the key basis.

        EXAMPLES::

            sage: k = KeyPolynomials(ZZ)
            sage: X = SchubertPolynomialRing(ZZ)
            sage: f = X([2,1,5,4,3])
            sage: k.from_schubert_polynomial(f)
            k[1, 0, 2, 1] + k[2, 0, 2] + k[3, 0, 0, 1]
            sage: k.from_schubert_polynomial(2)
            2*k[]
            sage: k(f)
            k[1, 0, 2, 1] + k[2, 0, 2] + k[3, 0, 0, 1]

            sage: k = KeyPolynomials(GF(7), 4)
            sage: k.from_schubert_polynomial(f)
            k[1, 0, 2, 1] + k[2, 0, 2, 0] + k[3, 0, 0, 1]

        TESTS::

            sage: k = KeyPolynomials(ZZ)
            sage: k.from_schubert_polynomial(k([3,2]))
            Traceback (most recent call last):
            ...
            ValueError: not a Schubert polynomial

            sage: k = KeyPolynomials(ZZ)
            sage: X = SchubertPolynomialRing(ZZ)
            sage: it = iter(Compositions())
            sage: for _ in range(50):
            ....:     C = next(it)
            ....:     assert k.from_schubert_polynomial(X(k[C])) == k[C], C

            sage: k = KeyPolynomials(ZZ, 4)
            sage: X = SchubertPolynomialRing(ZZ)
            sage: it = iter(k.basis().keys())
            sage: for _ in range(50):
            ....:     C = next(it)
            ....:     assert k.from_schubert_polynomial(X(k[C])) == k[C], C
        """
        if x in self.base_ring():
            return self(x)

        from sage.combinat.schubert_polynomial import SchubertPolynomial_class
        if not isinstance(x, SchubertPolynomial_class):
            raise ValueError('not a Schubert polynomial')

        from sage.combinat.diagram import RotheDiagram
        out = self.zero()
        if self._k is not None:
            def build_elt(wt):
                wt = list(wt)
                wt += [0] * (self._k - len(wt))
                return self[wt]
        else:
            def build_elt(wt):
                return self[wt]

        for m, c in x.monomial_coefficients().items():
            D = RotheDiagram(m)
            a = self.zero()
            for d in D.peelable_tableaux():
                a += build_elt(d.left_key_tableau().weight())
            out += c * a

        return out


def divided_difference(f, i):
    r"""
    Apply the ``i``-th divided difference operator to the polynomial ``f``.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import divided_difference
        sage: k = KeyPolynomials(QQ)
        sage: z = k.poly_gens()
        sage: f = z[1]*z[2]^3 + z[1]*z[2]*z[3]
        sage: divided_difference(f, 3)
        z_3^2*z_1 + z_3*z_2*z_1 + z_2^2*z_1

        sage: k = KeyPolynomials(QQ, 4)
        sage: z = k.poly_gens()
        sage: f = z[1]*z[2]^3 + z[1]*z[2]*z[3]
        sage: divided_difference(f, 3)
        z_1*z_2^2 + z_1*z_2*z_3 + z_1*z_3^2

        sage: k = KeyPolynomials(QQ)
        sage: R = k.polynomial_ring(); R
        Infinite polynomial ring in z over Rational Field
        sage: z = R.gen()
        sage: divided_difference(z[1]*z[2]^3, 2)
        -z_2^2*z_1 - z_2*z_1^2
        sage: divided_difference(z[1]*z[2]*z[3], 3)
        0
        sage: divided_difference(z[1]*z[2]*z[3], 4)
        z_2*z_1
        sage: divided_difference(z[1]*z[2]*z[4], 4)
        -z_2*z_1

        sage: k = KeyPolynomials(QQ, 5)
        sage: z = k.polynomial_ring().gens()
        sage: divided_difference(z[1]*z[2]^3, 2)
        -z_1^2*z_2 - z_1*z_2^2
        sage: divided_difference(z[1]*z[2]*z[3], 3)
        0
        sage: divided_difference(z[1]*z[2]*z[3], 4)
        z_1*z_2
        sage: divided_difference(z[1]*z[2]*z[4], 4)
        -z_1*z_2
    """
    P = parent(f)
    if isinstance(P, InfinitePolynomialRing_sparse):
        z = P.gen()
    else:
        z = P.gens()

    si_f = f.subs({z[i]: z[i-1], z[i-1]: z[i]})
    return (si_f - f) // (z[i] - z[i-1])


def isobaric_divided_difference(f, w):
    r"""
    Apply the isobaric divided difference operator `\pi_w` to the
    polynomial `f`.

    ``w`` may be either a single index or a list of
    indices of simple transpositions.

    .. WARNING::

        The simple transpositions should be applied from left to right.

    EXAMPLES::

        sage: from sage.combinat.key_polynomial import isobaric_divided_difference as idd
        sage: R.<z> = InfinitePolynomialRing(GF(3))
        sage: idd(z[1]^4*z[2]^2*z[4], 4)
        0

        sage: idd(z[1]^4*z[2]^2*z[3]*z[4], 3)
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: idd(z[1]^4*z[2]^2*z[3]*z[4], [3, 4])
        z_4^2*z_3*z_2*z_1^4 + z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: idd(z[1]^4*z[2]^2*z[3]*z[4], [4, 3])
        z_4*z_3^2*z_2*z_1^4 + z_4*z_3*z_2^2*z_1^4

        sage: idd(z[1]^2*z[2], [3, 2])
        z_3*z_2^2 + z_3*z_2*z_1 + z_3*z_1^2 + z_2^2*z_1 + z_2*z_1^2
    """
    P = parent(f)
    if isinstance(P, InfinitePolynomialRing_sparse):
        z = P.gen()
    else:
        z = P.gens()

    if not hasattr(w, "__iter__"):  # this allows us to pass i instead of a word
        w = [w]
    for i in w:
        fp = z[i-1] * f
        si_fp = fp.subs({z[i]: z[i-1], z[i-1]: z[i]})
        f = (si_fp - fp) // (z[i] - z[i-1])
    return f


def sorting_word(alpha):
    r"""
    Get a reduced word for the permutation which sorts ``alpha``
    into a partition.

    The result is a list ``l = [i0, i1, i2, ...]`` where each ``ij``
    is a positive integer such that it applies the simple
    transposition `(i_j, i_j+1)`. The transpositions are applied
    starting with ``i0``, then ``i1`` is applied, followed by ``i2``,
    and so on. See :meth:`sage.combinat.permutation.Permutation.reduced_words`
    for the convention used.

    EXAMPLES::

        sage: IV = IntegerVectors()
        sage: from sage.combinat.key_polynomial import sorting_word
        sage: list(sorting_word(IV([2,3,2]))[0])
        [1]
        sage: sorting_word(IV([2,3,2]))[1]
        [3, 2, 2]
        sage: list(sorting_word(IV([5,6,7]))[0])
        [1, 2, 1]
        sage: list(sorting_word(IV([0,3,2]))[0])
        [2, 1]
        sage: list(sorting_word(IV([0,3,0,2]))[0])
        [2, 3, 1]
        sage: list(sorting_word(IV([3,2,1]))[0])
        []
        sage: list(sorting_word(IV([2,3,3]))[0])
        [2, 1]
    """
    w = []
    L = list(alpha)
    n = len(L)

    # bubble sort to get the shortest sorting word
    for i in range(n-1):
        for j in range(n-i-1):
            if L[j] < L[j + 1]:
                w.append(j+1)
                L[j], L[j + 1] = L[j + 1], L[j]
    return reversed(w), L
