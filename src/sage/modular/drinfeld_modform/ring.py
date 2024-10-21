r"""
Graded rings of Drinfeld modular forms

This module defines a class named :class:`DrinfeldModularForms`.
Currently, the implementation only supports the full modular group
`\mathrm{GL}_r(A)` where `A = \mathbb{F}_q[T]`.

The implementation is based on the following identification:

.. MATH::

    M^r(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1, \ldots, g_{r-1}, g_{r}].

where `g_i` is the `i`-th coefficient form of weight `q^{i} - 1`.

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <da.ayotte@outlook.com>
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

from .element import DrinfeldModularFormsElement


class DrinfeldModularForms(Parent, UniqueRepresentation):
    r"""
    Base class for the graded ring of Drinfeld modular forms.

    If `K = \mathrm{Frac}(A)` where `A = \mathbb{F}_q[T]`, then the
    ring of Drinfeld modular forms over `K` of rank `r` and type zero
    for `\mathrm{GL}_r(A)` is

    .. MATH::

        M^{r, 0}(\mathrm{GL}_r(A)) = K[g_1, \ldots, g_{r-1}, g_{r}].

    where `g_i` the `i`-th coefficient form of weight `q^{i} - 1` at
    `T`.

    Similarly, the ring of Drinfeld modular forms over `K` of rank `r`
    and arbitrary type is

    .. MATH::

        M^{r}(\mathrm{GL}_r(A)) = K[g_1, \ldots, g_{r-1}, h_{r}].

    where `h_r` is a form of weight `(q^r - 1)/(q - 1)` and type `1`.

    We will see the elements of this ring as formal objects given by
    algebraic combination of the generator of the ring. See the class
    :class:`~sage.modular.drinfeld_modform.element.DrinfeldModularFormsElement`
    for more details about their implementation.

    INPUT:

    - ``base_ring`` -- the fraction field of a univariate polynomial
      ring over `\mathbb{F}_q`

    - ``rank`` integer (default: ``None``); the rank of the ring. If
      the rank is ``None``, then the names of the generators must be
      specified.

    - ``group`` -- (not implemented, default: ``None``) the group of the
      ring. The current implementation only supports the full modular
      group `\mathrm{GL}_r(A)`.

    - ``has_type`` -- boolean (default: ``False``); if set to ``True``,
      returns the graded ring of arbitrary type

    - ``names`` -- string, tuple or list (default: ``None``); a single
      character, a tuple or list of character, or comma separated string
      of character representing the names of the generators. If this
      parameter is set to ``None`` and the rank is specified, then the
      default names for the generators will be:

      * ``g1, g2, ..., gr`` for the type zero forms

      * ``g1, g2, ..., hr`` for the arbitrary type forms.

      If this parameter is a single character, for example ``f``, and a
      rank is specified, then the names will be of the form
      ``f1, f2, ..., fr``. Finally, if this parameter is a list, a tupe
      or a string of comma separated characters, then each character
      will corresponds to a generator. Note that in this case, it not
      necessary to specify the rank.

    EXAMPLES::

        sage: q = 3
        sage: A = GF(q)['T']
        sage: K.<T> = Frac(A)
        sage: M = DrinfeldModularForms(K, 3)
        sage: M
        Ring of Drinfeld modular forms of rank 3 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

    Use the :meth:`gens` method to obtain the generators of the ring::

        sage: M.gens()
        [g1, g2, g3]
        sage: M.inject_variables()  # assign the variable g1, g2, g3
        Defining g1, g2, g3
        sage: T*g1*g2 + g3
        g3 + T*g1*g2

    When creating the ring, one can name the generators in various
    ways::

        sage: M.<F, G, H> = DrinfeldModularForms(K)
        sage: M.gens()
        [F, G, H]
        sage: M = DrinfeldModularForms(K, 5, names='f')  # must specify the rank
        sage: M.gens()
        [f1, f2, f3, f4, f5]
        sage: M = DrinfeldModularForms(K, names='u, v, w, x')
        sage: M.gens()
        [u, v, w, x]
        sage: M = DrinfeldModularForms(K, names=['F', 'G', 'H'])
        sage: M.gens()
        [F, G, H]

    Set the keyword parameter ``has_type`` to ``True`` in order to create
    the ring of Drinfeld modular forms of arbitrary type::

        sage: M = DrinfeldModularForms(K, 4, has_type=True)
        sage: M.gens()
        [g1, g2, g3, h4]
        sage: h4 = M.3
        sage: h4.type()
        1

    To obtain a generating set of the subspace of forms of a fixed
    weight, use the methode :meth:`basis_of_weight`::

        sage: M = DrinfeldModularForms(K, 2)
        sage: M.basis_of_weight(q^3 - 1)
        [g1*g2^3, g1^5*g2^2, g1^9*g2, g1^13]

    In order to compute the coefficient forms, use the methods
    :meth:`coefficient_form` and :meth:`coefficient_forms`::

        sage: M = DrinfeldModularForms(K, 3)
        sage: M.coefficient_form(1)
        g1
        sage: M.coefficient_form(2)
        g2
        sage: M.coefficient_form(3)
        g3
        sage: M.coefficient_forms(T)
        [g1, g2, g3]
        sage: M.coefficient_forms(T^2)
        [(T^3 + T)*g1,
         g1^4 + (T^9 + T)*g2,
         g1^9*g2 + g1*g2^3 + (T^27 + T)*g3,
         g1^27*g3 + g1*g3^3 + g2^10,
         g2^27*g3 + g2*g3^9,
         g3^28]

    REFERENCE:

    For a quick introduction to Drinfeld modular forms, see the
    :ref:`tutorial <sage.modular.drinfeld_modform.tutorial>`. For more
    extensive references, see [Gek1988]_ and [BRP2018]_.
    """

    Element = DrinfeldModularFormsElement

    @staticmethod
    def __classcall_private__(cls, base_ring, rank=None, group=None,
                              has_type=False, names=None):
        r"""
        Check input validity and return a ``DrinfeldModularForms``
        object.

        TESTS::

            sage: K = Frac(GF(2)['T'])
            sage: DrinfeldModularForms(K, rank=x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to an integer

            sage: DrinfeldModularForms(GF(2)['T'])
            Traceback (most recent call last):
            ...
            TypeError: base ring must be a fraction field of a polynomial ring

            sage: k.<x, y> = GF(5)[]
            sage: K = k.quotient(x+y).fraction_field()
            sage: DrinfeldModularForms(K, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Drinfeld modular forms are currently only implemented for A = Fq[T]

            sage: DrinfeldModularForms(Frac(ZZ['T']))
            Traceback (most recent call last):
            ...
            ValueError: base ring characteristic must be finite

            sage: R = GF(2)['u']
            sage: K = Frac(R['T'])
            sage: DrinfeldModularForms(K, 2)
            Traceback (most recent call last):
            ...
            ValueError: the ring of constants must be a field

            sage: K = Frac(Frac(R)['T'])
            sage: DrinfeldModularForms(K, 2)
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
            TypeError: names must be None, a comma separated string or a list of string

            sage: DrinfeldModularForms(Frac(GF(2)['T']), rank=3, names='f1, f2, f3, f4')
            Traceback (most recent call last):
            ...
            ValueError: the number of generators (=4) must be equal to the rank (=3)

            sage: DrinfeldModularForms(Frac(GF(2)['T']), rank=3, names=['f1', 'f2', 'f3', 'f4'])
            Traceback (most recent call last):
            ...
            ValueError: the number of generators (=4) must be equal to the rank (=3)

            sage: DrinfeldModularForms(Frac(GF(2)['T']), rank=3, names=('f1', 'f2', 'f3', 'f4'))
            Traceback (most recent call last):
            ...
            ValueError: the number of generators (=4) must be equal to the rank (=3)

            sage: DrinfeldModularForms(Frac(GF(2)['T']))
            Traceback (most recent call last):
            ...
            TypeError: rank or names must be specified
        """
        if not isinstance(base_ring, FractionField_generic):
            raise TypeError("base ring must be a fraction field of a "
                            "polynomial ring")
        if not isinstance(base_ring.base(), PolynomialRing_general):
            raise NotImplementedError("Drinfeld modular forms are currently "
                                      "only implemented for A = Fq[T]")
        if not base_ring.characteristic():
            raise ValueError("base ring characteristic must be finite")
        if not base_ring.base().base().is_field():
            raise ValueError("the ring of constants must be a field")
        if not base_ring.base().base().is_finite():
            raise ValueError("the ring of constants must be finite")
        if group is not None:  # placeholder
            raise NotImplementedError("Drinfeld modular forms are currently "
                                      "only implemented for the full group")
        if names is None: # default names
            if rank is None:
                raise TypeError("rank or names must be specified")
            rank = ZZ(rank)  # check the type of rank
            names = [f'g{i}' for i in range(1, rank, 1)]
            last = f"h{rank}" if has_type else f"g{rank}"
            names.append(last)
        elif isinstance(names, (str, list, tuple)):
            if isinstance(names, str):
                names = names.replace(' ', '').split(',')
            nb_names = len(names)
            if rank is None:
                rank = nb_names
            else:
                rank = ZZ(rank)
                if nb_names == 1 and rank > 1:
                    g = names[0]
                    names = [f'{g}{i}' for i in range(1, rank + 1)]
                elif nb_names != rank:
                    raise ValueError(f"the number of generators (={nb_names}) "
                                     f"must be equal to the rank (={rank})")
        else:
            raise TypeError("names must be None, a comma separated string "
                            "or a list of string")
        return cls.__classcall__(cls, base_ring, rank, group, has_type,
                                 tuple(names))

    def __init__(self, base_ring, rank, group, has_type, names):
        r"""
        Initialize the ``DrinfeldModularForms`` class.

        TESTS::

            sage: K = Frac(GF(3)['T'])
            sage: TestSuite(DrinfeldModularForms(K, 2)).run()
            sage: TestSuite(DrinfeldModularForms(K, 3)).run()
            sage: TestSuite(DrinfeldModularForms(K, 4)).run()

            sage: K = Frac(GF(7)['T'])
            sage: TestSuite(DrinfeldModularForms(K, 2)).run()
            sage: TestSuite(DrinfeldModularForms(K, 3)).run()
            sage: TestSuite(DrinfeldModularForms(K, 4)).run()

            sage: K = Frac(GF(2^3)['T'])
            sage: TestSuite(DrinfeldModularForms(K, 2)).run()
            sage: TestSuite(DrinfeldModularForms(K, 3)).run()
            sage: TestSuite(DrinfeldModularForms(K, 4)).run()
        """
        self._has_type = has_type
        self._rank = rank
        self._base_ring = base_ring
        q = base_ring.base_ring().cardinality()
        if has_type:
            degs = [q**i - 1 for i in range(1, rank, 1)]
            degs.append((q**rank - 1) / (q - 1))
        else:
            degs = [q**i - 1 for i in range(1, rank + 1, 1)]
        self._poly_ring = PolynomialRing(base_ring, rank, names=names,
                                         order=TermOrder('wdeglex', degs))
        self._assign_names(names)
        cat = GradedAlgebras(base_ring).Commutative()
        super().__init__(base=base_ring, category=cat)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        TESTS::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.an_element()
            g1
        """
        return self.element_class(self, self._poly_ring.an_element())

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

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
        Return the `i`-th coefficient form at `T`.

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
            h3^2
        """
        if self._has_type and i == self.rank():
            q = self._base_ring.base_ring().cardinality()
            return self.gen(i-1)**(q - 1)
        return self.gen(i - 1)

    def _coefficient_forms(self, a):
        r"""
        Return the list of all coefficients of the universal Drinfeld
        module at `a`.

        This method is used in the methods :meth:`coefficient_form` and
        :meth:`coefficient_forms`. The main difference is that we don't
        check the input here.

        INPUT:

        - ``a`` -- an element in the ring of regular functions

        OUTPUT: list of Drinfeld modular forms. The `i`-th element of
        that list corresponds to the `(i+1)`-th coefficient form at `a`.

        TESTS::

            sage: K.<T> = Frac(GF(2)['T'])
            sage: M = DrinfeldModularForms(K, 2)
            sage: M._coefficient_forms(T)
            [g1, g2]
            sage: M._coefficient_forms(T^2)
            [(T^2 + T)*g1, g1^3 + (T^4 + T)*g2, g1^4*g2 + g1*g2^2, g2^5]
        """
        a = a.numerator()
        d = a.degree()
        poly_ring = PolynomialRing(self._base_ring, self.rank(), 'g')
        poly_ring_gens = poly_ring.gens()
        Frob = poly_ring.frobenius_endomorphism()
        gen = [self._base_ring.gen()]
        for g in poly_ring_gens:
            gen.append(g)
        ore_pol_ring = OrePolynomialRing(poly_ring, Frob, 't')
        gen = ore_pol_ring(gen)
        f = sum(c*(gen**idx) for idx, c in enumerate(a.coefficients(sparse=False)))
        coeff_forms = []
        for i in range(1, a.degree()*self.rank()+1):
            form = f[i]
            coeff_forms.append(form.subs({g: self._generator_coefficient_form(j+1)
                                          for j, g in enumerate(poly_ring_gens)}))
        return coeff_forms

    def coefficient_form(self, i, a=None):
        r"""
        Return the `i`-th coefficient form of the universal Drinfeld
        module over `\Omega^r(\mathbb{C}_{\infty})`:

        .. MATH::

            \phi_{w, a} = a + g_{1, a}\tau + \cdots + g_{r d_a, a}\tau^{r d_a}

        where `d_a := \mathrm{deg}(a)`.

        INPUT:

        - ``i`` -- integer between 1 and `r d_a`

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
            sage: M.coefficient_form(3, T^2)
            g1^9*g2 + g1*g2^3 + (T^27 + T)*g3

        ::

            sage: M = DrinfeldModularForms(K, 2, has_type=True)
            sage: M.coefficient_form(1)
            g1
            sage: M.coefficient_form(2)
            h2^2
            sage: M.coefficient_form(2, T^3 + T^2 + T)
            (T^9 + T^3 + T + 1)*g1^4 + (T^18 + T^10 + T^9 + T^2 + T + 1)*h2^2

        TESTS::

            sage: M.coefficient_form(0)
            Traceback (most recent call last):
            ...
            ValueError: index (=0) must be >= 1 and <= rank (=2)
            sage: M.coefficient_form(3)
            Traceback (most recent call last):
            ...
            ValueError: index (=3) must be >= 1 and <= rank (=2)
            sage: M.coefficient_form(9, T^2)
            Traceback (most recent call last):
            ...
            ValueError: index (=9) must be >= 1 and <= deg(a)*rank (=4)
            sage: M.coefficient_form(1, 1/(T+1))
            Traceback (most recent call last):
            ...
            ValueError: a must be an integral element
            sage: M.coefficient_form(1, 'x')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert a to an element in Fq[T]
        """
        i = ZZ(i)
        if a is None:
            if i < 1 or i > self.rank():
                raise ValueError(f"index (={i}) must be >= 1 and <= rank "
                                 f"(={self.rank()})")
            return self._generator_coefficient_form(i)
        try:
            A = self._base_ring.base()
            a = A(a)
        except TypeError:
            raise TypeError("unable to convert a to an element in Fq[T]")
        except ValueError:
            raise ValueError("a must be an integral element")
        if i < 1 or i > a.degree()*self.rank():
            raise ValueError(f"index (={i}) must be >= 1 and <= deg(a)*rank "
                             f"(={a.degree()*self.rank()})")
        coeff_forms = self._coefficient_forms(a)
        return coeff_forms[i - 1]

    def coefficient_forms(self, a=None):
        r"""
        Return the list of all coefficients of the universal Drinfeld
        module at `a`.

        See also :meth:`coefficient_form` for definitions.

        INPUT:

        - ``a`` -- (default: ``None``) an element in the ring of regular
          functions. If `a` is ``None``, then the method returns the
          coefficients forms at `a = T`.

        OUTPUT: list of Drinfeld modular forms. The `i`-th element of
        that list corresponds to the `(i+1)`-th coefficient form at `a`.

        EXAMPLES::

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

        TESTS::

            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.coefficient_forms(1/T)
            Traceback (most recent call last):
            ...
            ValueError: a must be an integral element
            sage: M.coefficient_forms('x')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert a to an element in Fq[T]
        """
        K = self._base_ring
        T = K.gen()
        if a is None:
            return [self._generator_coefficient_form(i)
                    for i in range(1, self.rank() + 1)]
        try:
            A = self._base_ring.base()
            a = A(a)
        except TypeError:
            raise TypeError("unable to convert a to an element in Fq[T]")
        except ValueError:
            raise ValueError("a must be an integral element")
        return self._coefficient_forms(a)

    def gen(self, n):
        r"""
        Return the `n`-th generator of this ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.gen(0)
            g1
            sage: M.1  # equivalent to M.gen(1)
            g2

        .. NOTE::

            Recall that the ring of Drinfeld modular forms is generated
            by the `r` coefficient forms of the universal Drinfeld
            module at `T`, `g_1, g_2, \ldots, g_r`, see
            :meth:`coefficient_forms`. We highlight however that we make
            a shift in the indexing so that the `i`-th generator
            corresponds to the `i+1`-th coefficient form for
            `0\leq i \leq r-1`.
        """
        return self(self._poly_ring.gen(n))

    def gens(self):
        r"""
        Return a list of generators of this ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 5)
            sage: M.gens()
            [g1, g2, g3, g4, g5]
        """
        return [self(g) for g in self._poly_ring.gens()]

    def ngens(self):
        r"""
        Return the number of generators of this ring.

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
        Return the rank of this ring of Drinfeld modular forms.

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
        Return the multiplicative unit.

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
        Return the additive identity.

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
        Return a list of Drinfeld modular forms which forms a basis for
        the subspace of weight `k`.

        Note that if `k\not\equiv 0` modulo `q-1`, then the subspace is
        0.

        An alias of this method is ``basis``.

        INPUT:

        - ``k`` -- integer

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
        each variable corresponds to a generator of this Drinfeld
        modular forms ring.

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
