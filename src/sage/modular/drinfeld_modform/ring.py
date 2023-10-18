r"""
This module defines a class named :class:`DrinfeldModularForms`.

Currently, the implementation only supports the full group modular group
`\mathrm{GL}_r(A)` where `A = \mathbb{F}_q[T]`.

The implementation is based on the following identification:

.. MATH::

    M^r(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1, \ldots, g_{r-1}, g_{r}].

where `g_i` the `i`-th coefficient form of weight `q^{i} - 1`.

EXAMPLES::

    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 2)  # rank 2
    sage: M.gens()  # generators
    [g1, g2]
    sage: M.inject_variables()  # assign variables
    Defining g1, g2
    sage: g1.weight()
    2
    sage: g2.weight()
    8

::

    sage: M = DrinfeldModularForms(K, 3)  # rank 3
    sage: M.gens()
    [g1, g2, g3]
    sage: [g.weight() for g in M.gens()]  # list of the weights
    [2, 8, 26]
    sage: M.inject_variables()
    Defining g1, g2, g3
    sage: g1.weight() == 3 - 1
    True
    sage: g2.weight() == 3^2 - 1
    True
    sage: g3.weight() == 3^3 - 1
    True

It is possible to compute the coefficient forms::

    sage: M = DrinfeldModularForms(K, 2)
    sage: M.coefficient_form(1)
    g1
    sage: M.coefficient_form(2)
    g2
    sage: M.coefficient_form(2, T^2)
    g1^4 + (T^9 + T)*g2
    sage: M.coefficient_forms(T^3)
    [(T^6 + T^4 + T^2)*g1,
     (T^9 + T^3 + T)*g1^4 + (T^18 + T^10 + T^2)*g2,
     g1^13 + (T^27 + T^9 + T)*g1^9*g2 + (T^27 + T^3 + T)*g1*g2^3,
     g1^36*g2 + g1^28*g2^3 + g1^4*g2^9 + (T^81 + T^9 + T)*g2^10,
     g1^81*g2^10 + g1^9*g2^28 + g1*g2^30,
     g2^91]

One can compute basis for any subspace of given weight::

    sage: M = DrinfeldModularForms(K, 4)
    sage: M.basis_of_weight(q^3 - 1)
    [g3, g1*g2^3, g1^5*g2^2, g1^9*g2, g1^13]

We note that the elements of this ring may not be *modular forms* as
as they may have mixed weight components::

    sage: M = DrinfeldModularForms(K, 4)
    sage: M.inject_variables()
    Defining g1, g2, g3, g4
    sage: F = g1 + g2 + g3 + g4
    sage: F.is_drinfeld_modular_form()
    False

This is why we call these elements *graded Drinfeld modular forms*.

One can also consider the ring Drinfeld modular forms of arbitrary type::

    sage: A = GF(7)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 4, has_type=True)
    sage: M.inject_variables()
    Defining g1, g2, g3, h
    sage: h.weight()
    400
    sage: h.type_m()
    1
    sage: (g1*h^4).type_m()
    4

The last generator is known as Gekeler's `h` function.

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <davidayotte94@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.graded_algebras import GradedAlgebras

from sage.structure.parent import Parent

from sage.rings.fraction_field import FractionField_generic
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.integer_ring import ZZ

from sage.structure.unique_representation import UniqueRepresentation

from .element import DrinfeldModularFormsElement, DrinfeldModularFormsElement_rank_two

class DrinfeldModularForms(Parent, UniqueRepresentation):
    r"""
    Base class for the graded Drinfeld modular forms ring.

    INPUT:

    - ``base_ring`` -- The fraction field of a univariate polynomial
      ring over `\mathbb{F}_q`.
    - ``rank`` (integer, default: 2) -- the rank of the ring
    - ``group`` (NoneType) -- the group of self. The current
      implementation only supports the full group
      `\mathrm{GL}_r(A)`.
    - ``has_type`` (bool, default: ``False``) -- if set to True, returns
      the graded ring of arbitrary type.
    - ``names`` (string, default: ``'g'``) -- a single character or a
      comma seperated string of character representing the names of the
      generators.

    TESTS::

        sage: K = Frac(GF(3)['T'])
        sage: TestSuite(DrinfeldModularForms(K)).run()
        sage: TestSuite(DrinfeldModularForms(K, 3)).run()
        sage: TestSuite(DrinfeldModularForms(K, 4)).run()
        sage: K = Frac(GF(7)['T'])
        sage: TestSuite(DrinfeldModularForms(K)).run()
        sage: TestSuite(DrinfeldModularForms(K, 3)).run()
        sage: TestSuite(DrinfeldModularForms(K, 4)).run()

    ::

        sage: K = Frac(GF(2)['T'])
        sage: DrinfeldModularForms(K, rank=x)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert x to an integer
        sage: DrinfeldModularForms(GF(2)['T'])
        Traceback (most recent call last):
        ...
        TypeError: base ring must be a fraction field of a polynomial ring
        sage: DrinfeldModularForms(Frac(ZZ['T']))
        Traceback (most recent call last):
        ...
        ValueError: base ring characteristic must be finite
        sage: R = GF(2)['u']
        sage: K = Frac(R['T'])
        sage: DrinfeldModularForms(K)
        Traceback (most recent call last):
        ...
        ValueError: the ring of constants must be a field
        sage: K = Frac(Frac(R)['T'])
        sage: DrinfeldModularForms(K)
        Traceback (most recent call last):
        ...
        ValueError: the ring of constants must be finite
        sage: DrinfeldModularForms(Frac(GF(2)['T']), group=2)
        Traceback (most recent call last):
        ...
        NotImplementedError: Drinfeld modular forms are currently only implemented for the full group
        sage: DrinfeldModularForms(Frac(GF(2)['T']), names=1)
        Traceback (most recent call last):
        ...
        TypeError: names must be a string
        sage: DrinfeldModularForms(Frac(GF(2)['T']), rank=3, names='f1, f2, f3, f4')
        Traceback (most recent call last):
        ...
        ValueError: the the number of generators must be equal to the rank (=3)
    """

    Element = DrinfeldModularFormsElement

    @staticmethod
    def __classcall_private__(cls, base_ring, rank=2, group=None,
                              has_type=False, names='g'):
        rank = ZZ(rank)  # check the type of rank
        if not isinstance(base_ring, FractionField_generic):
            raise TypeError("base ring must be a fraction field of a "
                            "polynomial ring")
        if not isinstance(base_ring.base(), PolynomialRing_general):  # not sure if this test is relevant
            raise TypeError("the base of the base ring must be a "
                            "polynomial ring")
        if not base_ring.characteristic():
            raise ValueError("base ring characteristic must be finite")
        if not base_ring.base().base().is_field():
            raise ValueError("the ring of constants must be a field")
        if not base_ring.base().base().is_finite():
            raise ValueError("the ring of constants must be finite")
        if group is not None:  # placeholder
            raise NotImplementedError("Drinfeld modular forms are currently "
                                      "only implemented for the full group")
        if not isinstance(names, str):
            raise TypeError("names must be a string")
        nb_names = len(names.split())
        if nb_names == 1:
            n0 = names
            names += "1, "
            for i in range(2, rank, 1):
                names += n0 + str(i) + ", "
            if has_type:
                names += 'h'
            else:
                names += n0 + str(rank)
        else:
            if nb_names != rank:
                raise ValueError("the the number of generators "
                                 f"must be equal to the rank (={rank})")
        if rank == 2:
            return DrinfeldModularForms_rank_two(base_ring, group,
                                                 has_type, names)
        return cls.__classcall__(cls, base_ring, rank, group, has_type, names)

    def __init__(self, base_ring, rank=2, group=None, has_type=False,
                 names='g'):
        self._has_type = has_type
        self._rank = rank
        self._base_ring = base_ring
        q = base_ring.base_ring().cardinality()
        if has_type:
            degs = [q**i - 1 for i in range(1, rank, 1)]
            degs.append((q**rank - 1)/(q - 1))
        else:
            degs = [q**i - 1 for i in range(1, rank + 1, 1)]
        self._poly_ring = PolynomialRing(base_ring, rank, names=names,
                                         order=TermOrder('wdeglex', degs))
        self._assign_names(names)

        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def _an_element_(self):
        r"""
        Return an element of self.

        TESTS::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K)
            sage: M.an_element()
            g1
        """
        return self.element_class(self, self._poly_ring.an_element())

    def _repr_(self):
        r"""
        Return the string representation of self.

        TESTS::

            sage: A = GF(2)['T']
            sage: K = Frac(A)
            sage: M = DrinfeldModularForms(K, 3)
            sage: M._repr_()
            'Ring of Drinfeld modular forms of rank 3 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 2 (using GF2X)'
        """
        return ("Ring of Drinfeld modular forms of rank %s over %s"
                % (self._rank, self._base_ring))

    def _generator_coefficient_form(self, i):
        r"""
        Return the i-th coefficient form at `T`.

        For internal use only, the user should use
        :meth:`coefficient_form` instead.

        TESTS::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 3)
            sage: M._generator_coefficient_form(1)
            g1
            sage: M._generator_coefficient_form(2)
            g2
            sage: M._generator_coefficient_form(3)
            g3

        ::

            sage: M = DrinfeldModularForms(K, 3, has_type=True)
            sage: M._generator_coefficient_form(1)
            g1
            sage: M._generator_coefficient_form(2)
            g2
            sage: M._generator_coefficient_form(3)
            h^2
            sage: M._generator_coefficient_form(0)
            Traceback (most recent call last):
            ...
            ValueError: index (=0) must be >= 1 and <= rank (=3)
            sage: M._generator_coefficient_form(4)
            Traceback (most recent call last):
            ...
            ValueError: index (=0) must be >= 1 and <= rank (=3)
        """
        i = ZZ(i)
        r = self.rank()
        if i < 1 or i > r:
            raise ValueError(f"index (={i}) must be >= 1 and <= rank (={r})")
        if self._has_type and i == r:
            q = self._base_ring.base_ring().cardinality()
            return self.gen(i-1)**(q - 1)
        return self.gen(i - 1)

    def coefficient_form(self, i, a=None):
        r"""
        Return the `i`-th coefficient form of the universal Drinfeld
        module over `\Omega^r(\mathbb{C}_{\infty})`:

        ..MATH::

            \phi_{w, a} = a + g_{1, a}\tau + \cdots + g_{r d_a, a}\tau^{r d_a}

        where `d_a := \mathrm{deg}(a)`.

        INPUT:

        - ``i`` -- an integer between 1 and `r d_a`;

        - ``a`` -- (default: ``None``) an element in the ring of regular
          functions. If `a` is ``None``, then the method returns the
          `i`-th coefficient form of `\phi_{w, T}`.

        EXAMPLES::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 3)
            sage: M.coefficient_form(1)
            g1
            sage: M.coefficient_form(2)
            g2
            sage: M.coefficient_form(3)
            g3
            sage: M.coefficient_form(5, T^2)
            g2^27*g3 + g2*g3^9

        ::

            sage: M = DrinfeldModularForms(K, 2, has_type=True)
            sage: M.coefficient_form(1)
            g1
            sage: M.coefficient_form(2)
            h^2
            sage: M.coefficient_form(2, T^3 + T^2 + T)
            (T^9 + T^3 + T + 1)*g1^4 + (T^18 + T^10 + T^9 + T^2 + T + 1)*h^2
        """
        if i not in ZZ:
            raise TypeError("i must be an integer")
        i = ZZ(i)
        if i < 1:
            raise ValueError("i must be >= 1")
        if a is None:
            return self._generator_coefficient_form(i)
        if a not in self._base_ring:
            raise TypeError("a should be an element of the base ring")
        a = self._base_ring(a)
        if not a.denominator().is_one():
            raise ValueError("a should be in the ring of regular functions")
        a = a.numerator()
        poly_ring = PolynomialRing(self._base_ring, self.rank(), 'g')
        poly_ring_gens = poly_ring.gens()
        Frob = poly_ring.frobenius_endomorphism()
        gen = [self._base_ring.gen()]
        for g in poly_ring_gens:
            gen.append(g)
        ore_pol_ring = OrePolynomialRing(poly_ring, Frob, 't')
        gen = ore_pol_ring(gen)
        f = sum(c*(gen**idx) for idx, c in enumerate(a.coefficients(sparse=False)))
        form = f[i]
        coeff_form = form.subs({g: self._generator_coefficient_form(j+1) for j, g in enumerate(poly_ring_gens)})
        return coeff_form

    def coefficient_forms(self, a=None):
        r"""
        Return the list of all coefficient forms at `a`.

        See also :meth:`coefficient_form`.

        EXAMPLE::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.coefficient_forms()
            [g1, g2]
            sage: M.coefficient_forms(T^2)
            [(T^3 + T)*g1, g1^4 + (T^9 + T)*g2, g1^9*g2 + g1*g2^3, g2^10]
            sage: M.coefficient_forms(T^3)
            [(T^6 + T^4 + T^2)*g1,
             (T^9 + T^3 + T)*g1^4 + (T^18 + T^10 + T^2)*g2,
             g1^13 + (T^27 + T^9 + T)*g1^9*g2 + (T^27 + T^3 + T)*g1*g2^3,
             g1^36*g2 + g1^28*g2^3 + g1^4*g2^9 + (T^81 + T^9 + T)*g2^10,
             g1^81*g2^10 + g1^9*g2^28 + g1*g2^30,
             g2^91]
        """
        if a is None:
            return [self._generator_coefficient_form(i) for i in range(1, self.rank() + 1)]
        a = self._base_ring(a)
        if not a.denominator().is_one():
            raise ValueError("the input should be in the ring of regular"
                             " functions")
        d = a.numerator().degree()
        return [self.coefficient_form(i, a) for i in range(1, self.rank()*d + 1)]

    def gen(self, n):
        r"""
        Return the `n`-th generator of the ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0
            g1
            sage: M.1
            g2
        """
        return self(self._poly_ring.gen(n))

    def gens(self):
        r"""
        Return a list of generators for this ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 5)
            sage: M.gens()
            [g1, g2, g3, g4, g5]
        """
        return [self(g) for g in self._poly_ring.gens()]

    def ngens(self):
        r"""
        Return the number of generators of the ring.

        Note that the number of generators is equal to the rank.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 5)
            sage: M.ngens()
            5
        """
        return self._rank

    def _element_constructor_(self, polynomial):
        r"""
        Return the element corresponding to the given polynomial.

        TESTS::

            sage: A = GF(3)['T']
            sage: K = Frac(A)
            sage: M = DrinfeldModularForms(K, 3)
            sage: g1, g2, g3 = polygens(K, 3, 'g1, g2, g3')
            sage: M(g1*g2 + g3^4)
            g3^4 + g1*g2
        """
        return self.element_class(self, polynomial)

    def rank(self):
        r"""
        Return the rank of the ring of Drinfeld modular forms.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A);
            sage: DrinfeldModularForms(K, 2).rank()
            2
            sage: DrinfeldModularForms(K, 3).rank()
            3
            sage: DrinfeldModularForms(K, 4).rank()
            4
        """
        return self._rank

    def one(self):
        r"""
        Return the multiplicative identity of the ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.one()
            1
            sage: M.one() * M.0
            g1
            sage: M.one().is_one()
            True
        """
        return self(self._poly_ring.one())

    def zero(self):
        r"""
        Return the additive identity of the ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.zero()
            0
            sage: M.zero() + M.1
            g2
            sage: M.zero() * M.1
            0
            sage: M.zero().is_zero()
            True
        """
        return self(self._poly_ring.zero())

    def basis_of_weight(self, k):
        r"""
        Return a list of Drinfeld modular forms which forms a basis for the
        subspace of weight `k`.

        Note that if `k\not\equiv 0` modulo `q-1`, then the subspace is 0.

        An alias of this method is ``basis``.

        INPUT:

        - ``k`` -- an integer.

        EXAMPLES::

            sage: q = 3; A = GF(q)['T']; K = Frac(A);
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.basis_of_weight(q - 1)
            [g1]
            sage: M.basis_of_weight(q^2 - 1)
            [g2, g1^4]
            sage: M.basis_of_weight(q^3 - 1)
            [g1*g2^3, g1^5*g2^2, g1^9*g2, g1^13]
            sage: M.basis_of_weight(19*(q-1))
            [g1^3*g2^4, g1^7*g2^3, g1^11*g2^2, g1^15*g2, g1^19]
        """
        return [self(mon) for mon in self._poly_ring.monomials_of_degree(k)]

    basis = basis_of_weight  # alias

    def polynomial_ring(self):
        r"""
        Return the multivariate polynomial ring over the base ring where
        each variable corresponds to a generator of a ring.

        EXAMPLES::

            sage: q = 3; A = GF(q)['T']; K = Frac(A);
            sage: M = DrinfeldModularForms(K, 2)
            sage: P = M.polynomial_ring()
            sage: P
            Multivariate Polynomial Ring in g1, g2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

        The degree of the variables corresponds to the weight of the
        associated generator::

            sage: P.inject_variables()
            Defining g1, g2
            sage: g1.degree()
            2
            sage: g2.degree()
            8
        """
        return self._poly_ring


class DrinfeldModularForms_rank_two(DrinfeldModularForms):

    Element = DrinfeldModularFormsElement_rank_two

    def __init__(self, base_ring, group, has_type, names):
        super().__init__(base_ring, 2, group, has_type, names)

    def from_expansion(self, expansion, weight):
        return

    def from_petrov_expansion(self, k, n, name='t'):
        return

    def sturm_bound(self, weight):
        return

    def eisenstein_series(self, weight):
        return
