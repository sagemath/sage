r"""
Univariate skew polynomial rings

This module provides the
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`.
In the class hierarchy in Sage, the locution *Skew Polynomial* is used
for a Ore polynomial without twisting derivation.

This module also provides:

- the class :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order`,
  which is a specialized class for skew polynomial rings over fields
  equipped with an automorphism of finite order. It inherits from
  :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`
  but contains more methods and provides better algorithms.

- the class :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_field`,
  which is a specialized class for skew polynomial rings over finite fields.

.. SEEALSO::

    :class:`~sage.rings.polynomial.ore_polynomial_ring.OrePolynomialRing`

AUTHOR:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests
  and refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and
  doctest errors
"""

# ***************************************************************************
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.category_object import normalize_names

from sage.rings.ring import Field
from sage.matrix.matrix_space import MatrixSpace

from sage.rings.morphism import RingHomomorphism
from sage.categories.homset import Hom
from sage.categories.map import Section

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

WORKING_CENTER_MAX_TRIES = 1000


# Helper functions

def _base_ring_to_fraction_field(S):
    r"""
    Return the unique skew polynomial ring over the fraction field of
    ``S.base_ring()`` which has ``S`` a sub-ring (internal method).

    INPUT:

    - ``S`` -- a skew polynomial ring

    OUTPUT:

    - ``Q`` -- the skew polynomial ring over the fraction field
      of ``S.base_ring``

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _base_ring_to_fraction_field
        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t + 1])
        sage: S.<x> = R['x', sigma]
        sage: _base_ring_to_fraction_field(S)
        Ore Polynomial Ring in x over Fraction Field of Univariate Polynomial Ring in t over Integer Ring twisted by t |-->  t + 1
    """
    R = S.base_ring()
    if isinstance(R, Field):
        return S
    else:
        Q = R.fraction_field()
        gens = R.gens()
        sigmaS = S.twisting_morphism()
        # try:
        sigmaQ = Q.hom([Q(sigmaS(g)) for g in gens])
        return Q[S.variable_name(), sigmaQ]
        # except Exception, e:
        #     raise ValueError("unable to lift the twisting morphism to a twisting morphism over %s (error was: %s)" % (Q, e))


def _minimal_vanishing_polynomial(R, eval_pts):
    r"""
    Return the minimal vanishing polynomial (internal function).

    See the documentation for
    :meth:`SkewPolynomialRing.minimal_vanishing_polynomial` for a description.

    INPUT:

    - ``R`` -- a skew polynomial ring over a field

    - ``eval_pts`` -- a list of evaluation points

    OUTPUT:

    The minimal vanishing polynomial.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.polynomial.skew_polynomial_ring import _minimal_vanishing_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [1, t, t^2]
        sage: b = _minimal_vanishing_polynomial(S, eval_pts); b
        doctest:...: FutureWarning: This class/method/function is marked as experimental.
         It, its functionality or its interface might change without a formal deprecation.
        See https://github.com/sagemath/sage/issues/13215 for details.
        x^3 + 4
    """
    l = len(eval_pts)
    if l == 0:
        return R.one()
    elif l == 1:
        e = eval_pts[0]
        if e.is_zero():
            return R.one()
        else:
            return R.gen() - R.twisting_morphism()(e) / e
    else:
        t = l // 2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        B_moved = M_A.multi_point_evaluation(B)
        M_at_B_moved = _minimal_vanishing_polynomial(R, B_moved)
        return M_at_B_moved * M_A


def _lagrange_polynomial(R, eval_pts, values):
    r"""
    Return the Lagrange polynomial of the given points if it exists.

    Otherwise return an unspecified polynomial (internal method).

    See the documentation for
    :meth:`SkewPolynomialRing.lagrange_polynomial` for a description
    of Lagrange polynomial.

    INPUT:

    - ``R`` -- a skew polynomial ring over a field

    - ``eval_pts`` -- list of evaluation points

    - ``values`` -- list of values that the Lagrange polynomial takes
        at the respective ``eval_pts``

    OUTPUT:

    - the Lagrange polynomial

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.polynomial.skew_polynomial_ring import _lagrange_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x', Frob]
        sage: eval_pts = [t , t^2]
        sage: values = [3*t^2 + 4*t + 4, 4*t]
        sage: d = _lagrange_polynomial(S, eval_pts, values); d
        x + t
        sage: d.multi_point_evaluation(eval_pts) == values
        True

    The following restrictions are impossible to satisfy because the evaluation
    points are linearly dependent over the fixed field of the twisting morphism, and the
    corresponding values do not match::

        sage: eval_pts = [t, 2*t]                                                       # needs sage.rings.finite_rings
        sage: values = [1, 3]
        sage: _lagrange_polynomial(S, eval_pts, values)                                 # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: the given evaluation points are linearly dependent over the fixed field of the twisting morphism,
        so a Lagrange polynomial could not be determined (and might not exist)
    """
    l = len(eval_pts)
    if l == 1:
        if eval_pts[0].is_zero():
            # This is due to linear dependence among the eval_pts.
            raise ValueError("the given evaluation points are linearly dependent"
                             " over the fixed field of the twisting morphism,"
                             " so a Lagrange polynomial could not be determined"
                             " (and might not exist)")
        return (values[0] / eval_pts[0]) * R.one()
    else:
        t = l // 2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        M_B = _minimal_vanishing_polynomial(R, B)
        A_ = M_B.multi_point_evaluation(A)
        B_ = M_A.multi_point_evaluation(B)
        I_1 = _lagrange_polynomial(R, A_, values[:t])
        I_2 = _lagrange_polynomial(R, B_, values[t:])
        return I_1 * M_B + I_2 * M_A


# Generic implementation of skew polynomial rings
#################################################

class SkewPolynomialRing(OrePolynomialRing):
    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``twisting_morphism`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``category`` -- a category

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t + 1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.category()
            Category of algebras over Univariate Polynomial Ring in t over Integer Ring
            sage: S([1]) + S([-1])
            0
            sage: TestSuite(S).run()
        """
        if derivation is not None:
            raise NotImplementedError
        if self.Element is None:
            import sage.rings.polynomial.skew_polynomial_element
            self.Element = sage.rings.polynomial.skew_polynomial_element.SkewPolynomial_generic_dense
        OrePolynomialRing.__init__(self, base_ring, morphism, None, name, sparse, category)

    def minimal_vanishing_polynomial(self, eval_pts):
        r"""
        Return the minimal-degree, monic skew polynomial which vanishes at all
        the given evaluation points.

        The degree of the vanishing polynomial is at most the length of
        ``eval_pts``. Equality holds if and only if the elements of
        ``eval_pts`` are linearly independent over the fixed field of
        ``self.twisting_morphism()``.

        - ``eval_pts`` -- list of evaluation points which are linearly
          independent over the fixed field of the twisting morphism of
          the associated skew polynomial ring

        OUTPUT:

        The minimal vanishing polynomial.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: eval_pts = [1, t, t^2]
            sage: b = S.minimal_vanishing_polynomial(eval_pts); b
            x^3 + 4

        The minimal vanishing polynomial evaluates to 0 at each of
        the evaluation points::

            sage: eval = b.multi_point_evaluation(eval_pts); eval                       # needs sage.rings.finite_rings
            [0, 0, 0]

        If the evaluation points are linearly dependent over the fixed
        field of the twisting morphism, then the returned polynomial has
        lower degree than the number of evaluation points::

            sage: S.minimal_vanishing_polynomial([t])                                   # needs sage.rings.finite_rings
            x + 3*t^2 + 3*t
            sage: S.minimal_vanishing_polynomial([t, 3*t])                              # needs sage.rings.finite_rings
            x + 3*t^2 + 3*t
        """
        return _minimal_vanishing_polynomial(_base_ring_to_fraction_field(self), eval_pts)

    def lagrange_polynomial(self, points):
        r"""
        Return the minimal-degree polynomial which interpolates the given
        points.

        More precisely, given `n` pairs `(x_1, y_1), \ldots, (x_n, y_n) \in R^2`,
        where `R` is ``self.base_ring()``, compute a skew polynomial `p(x)`
        such that `p(x_i) = y_i` for each `i`, under the condition that
        the `x_i` are linearly independent over the fixed field of
        ``self.twisting_morphism()``.

        If the `x_i` are linearly independent over the fixed field of
        ``self.twisting_morphism()`` then such a polynomial is guaranteed
        to exist. Otherwise, it might exist depending on the `y_i`, but
        the algorithm used in this implementation does not support that,
        and so an error is always raised.

        INPUT:

        - ``points`` -- a list of pairs `(x_1, y_1), \ldots, (x_n, y_n)` of
          elements of the base ring of ``self``; the `x_i` should be linearly
          independent over the fixed field of ``self.twisting_morphism()``

        OUTPUT:

        The Lagrange polynomial.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: points = [(t, 3*t^2 + 4*t + 4), (t^2, 4*t)]
            sage: d = S.lagrange_polynomial(points); d
            x + t

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t + 1])
            sage: T.<x> = R['x', sigma]
            sage: points = [(1, t^2 + 3*t + 4), (t, 2*t^2 + 3*t + 1), (t^2, t^2 + 3*t + 4)]
            sage: p = T.lagrange_polynomial(points); p
            ((-t^4 - 2*t - 3)/-2)*x^2 + (-t^4 - t^3 - t^2 - 3*t - 2)*x
             + (-t^4 - 2*t^3 - 4*t^2 - 10*t - 9)/-2
            sage: p.multi_point_evaluation([1, t, t^2]) == [t^2 + 3*t + 4, 2*t^2 + 3*t + 1, t^2 + 3*t + 4]
            True

        If the `x_i` are linearly dependent over the fixed field of
        ``self.twisting_morphism()``, then an error is raised::

            sage: T.lagrange_polynomial([(t, 1), (2*t, 3)])
            Traceback (most recent call last):
            ...
            ValueError: the given evaluation points are linearly dependent over the fixed field of the twisting morphism,
            so a Lagrange polynomial could not be determined (and might not exist)
        """
        l = len(points)
        if not all(len(pair) == 2 for pair in points):
            raise TypeError("supplied points must be pairs of elements of base ring")
        eval_pts = [x for (x, _) in points]
        values = [y for (_, y) in points]

        if l > len(set(eval_pts)):
            raise TypeError("the evaluation points must be distinct")
        zero_i = [i for i in range(l) if eval_pts[i].is_zero()]
        if zero_i and not values[zero_i[0]].is_zero():
            raise TypeError("a skew polynomial always evaluates to 0 at 0, but a non-zero value was requested")

        return _lagrange_polynomial(_base_ring_to_fraction_field(self), eval_pts, values)


# Special classes for twisting morphisms with finite order
##########################################################

class SectionSkewPolynomialCenterInjection(Section):
    r"""
    Section of the canonical injection of the center of a skew
    polynomial ring into this ring.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.center()
        sage: iota = S.convert_map_from(Z)
        sage: sigma = iota.section()
        sage: TestSuite(sigma).run(skip=['_test_category'])
    """
    def _call_(self, x):
        r"""
        Return `x` viewed as an element of the center.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(x^3)
            z
            sage: sigma(x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in the center
        """
        order = self.inverse()._order
        section = self.inverse()._embed.section()
        lx = x.list()
        l = []
        mod = 0
        for c in lx:
            if mod == 0:
                l.append(section(c))
            else:
                if not c.is_zero():
                    raise ValueError("%s is not in the center" % x)
            mod += 1
            if mod == order:
                mod = 0
        return self.codomain()(l)

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: sigma = iota.section()

            sage: # needs sage.rings.finite_rings
            sage: s = loads(dumps(sigma))
            sage: s == sigma
            True
            sage: s != sigma
            False
            sage: s is sigma
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        if op == op_NE:
            return (self.domain() is not right.domain()) or (self.codomain() is not right.codomain())
        return NotImplemented


class SkewPolynomialCenterInjection(RingHomomorphism):
    r"""
    Canonical injection of the center of a skew polynomial ring
    into this ring.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.center()
        sage: iota = S.convert_map_from(Z)
        sage: TestSuite(iota).run(skip=['_test_category'])
    """
    def __init__(self, domain, codomain, embed, order):
        r"""
        Initialize this morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: S.convert_map_from(Z)   # indirect doctest
            Embedding of the center of Ore Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
        """
        RingHomomorphism.__init__(self, Hom(domain, codomain))
        self._embed = embed
        self._order = order
        self._codomain = codomain
        self._section = SectionSkewPolynomialCenterInjection(self)

    def _repr_(self):
        r"""
        Return a string representation of this morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: iota
            Embedding of the center of Ore Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
            sage: iota._repr_()
            'Embedding of the center of Ore Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring'
        """
        return "Embedding of the center of %s into this ring" % self._codomain

    def _call_(self, x):
        r"""
        Return the image of `x` by this morphism.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z.<z> = S.center()
            sage: iota = S.convert_map_from(Z)

            sage: iota(z)                                                               # needs sage.rings.finite_rings
            x^3
        """
        k = self._codomain.base_ring()
        l = []
        lz = [k(0)] * (self._order - 1)
        for c in x.list():
            l += [self._embed(c)] + lz
        return self._codomain(l)

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)

            sage: # needs sage.rings.finite_rings
            sage: i = loads(dumps(iota))
            sage: i == iota
            True
            sage: i != iota
            False
            sage: i is iota
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        if op == op_NE:
            return (self.domain() is not right.domain()) or (self.codomain() is not right.codomain())
        return NotImplemented

    def section(self):
        r"""
        Return a section of this morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(x^3)
            z
        """
        return self._section


class SkewPolynomialRing_finite_order(SkewPolynomialRing):
    """
    A specialized class for skew polynomial rings whose twising morphism
    has finite order.

    .. SEEALSO::

        - :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`
        - :mod:`sage.rings.polynomial.skew_polynomial_finite_order`
    """
    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        r"""
        Initialize this skew polynomial ring.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]; S
            Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: S.category()
            Category of algebras over Finite Field in t of size 5^3

            sage: TestSuite(S).run()                                                    # needs sage.rings.finite_rings

        We check that a call to the method
        :meth:`sage.rings.polynomial.skew_polynomial_finite_order.SkewPolynomial_finite_order.is_central`
        does not affect the behaviour of default central variable names::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(7^4)
            sage: phi = k.frobenius_endomorphism()
            sage: S.<x> = k['x', phi]
            sage: (x^4).is_central()
            True
            sage: Z.<u> = S.center()
            sage: S.center() is Z
            True
        """
        if self.Element is None:
            import sage.rings.polynomial.skew_polynomial_finite_order
            self.Element = sage.rings.polynomial.skew_polynomial_finite_order.SkewPolynomial_finite_order_dense
        if self._fraction_field_class is None:
            from sage.rings.polynomial.ore_function_field import OreFunctionField_with_large_center
            self._fraction_field_class = OreFunctionField_with_large_center
        SkewPolynomialRing.__init__(self, base_ring, morphism, derivation, name, sparse, category)
        self._order = morphism.order()
        (self._constants, self._embed_constants) = morphism.fixed_field()

        # Configure and create center
        self._center = {}
        self._center_variable_name = 'z'
        for i in range(WORKING_CENTER_MAX_TRIES):
            try:
                self._working_center = self.center()
                self._center_variable_name = None
                break
            except ValueError:
                self._center_variable_name = "z%s_" % i
        if self._center_variable_name is not None:
            raise NotImplementedError("unable to create the center")

    def center(self, name=None, names=None, default=False):
        r"""
        Return the center of this skew polynomial ring.

        .. NOTE::

            If `F` denotes the subring of `R` fixed by `\sigma` and `\sigma`
            has order `r`, the center of `K[x,\sigma]` is `F[x^r]`, that
            is a univariate polynomial ring over `F`.

        INPUT:

        - ``name`` -- a string or ``None`` (default: ``None``);
          the name for the central variable (namely `x^r`)

        - ``default`` -- a boolean (default: ``False``); if ``True``,
          set the default variable name for the center to ``name``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]; S
            Ore Polynomial Ring in x over Finite Field in t of size 5^3
             twisted by t |--> t^5
            sage: Z = S.center(); Z
            Univariate Polynomial Ring in z over Finite Field of size 5
            sage: Z.gen()
            z

        We can pass in another variable name::

            sage: S.center(name='y')                                                    # needs sage.rings.finite_rings
            Univariate Polynomial Ring in y over Finite Field of size 5

        or use the bracket notation::

            sage: Zy.<y> = S.center(); Zy                                               # needs sage.rings.finite_rings
            Univariate Polynomial Ring in y over Finite Field of size 5
            sage: y.parent() is Zy                                                      # needs sage.rings.finite_rings
            True

        A coercion map from the center to the skew polynomial ring is set::

            sage: # needs sage.rings.finite_rings
            sage: S.has_coerce_map_from(Zy)
            True
            sage: P = y + x; P
            x^3 + x
            sage: P.parent()
            Ore Polynomial Ring in x over Finite Field in t of size 5^3
             twisted by t |--> t^5
            sage: P.parent() is S
            True

        together with a conversion map in the reverse direction::

            sage: Zy(x^6 + 2*x^3 + 3)                                                   # needs sage.rings.finite_rings
            y^2 + 2*y + 3

            sage: Zy(x^2)                                                               # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in the center

        Two different skew polynomial rings can share the same center::

            sage: S1.<x1> = k['x1', Frob]                                               # needs sage.rings.finite_rings
            sage: S2.<x2> = k['x2', Frob]                                               # needs sage.rings.finite_rings
            sage: S1.center() is S2.center()                                            # needs sage.rings.finite_rings
            True

        .. RUBRIC:: About the default name of the central variable

        A priori, the default is ``z``.

        However, a variable name is given the first time this method is
        called, the given name become the default for the next calls::

            sage: # needs sage.rings.finite_rings
            sage: K.<t> = GF(11^3)
            sage: phi = K.frobenius_endomorphism()
            sage: A.<X> = K['X', phi]
            sage: C.<u> = A.center()  # first call
            sage: C
            Univariate Polynomial Ring in u over Finite Field of size 11
            sage: A.center()  # second call: the variable name is still u
            Univariate Polynomial Ring in u over Finite Field of size 11
            sage: A.center() is C
            True

        We can update the default variable name by passing in the argument
        ``default=True``::

            sage: # needs sage.rings.finite_rings
            sage: D.<v> = A.center(default=True)
            sage: D
            Univariate Polynomial Ring in v over Finite Field of size 11
            sage: A.center()
            Univariate Polynomial Ring in v over Finite Field of size 11
            sage: A.center() is D
            True

        TESTS::

            sage: C.<a,b> = S.center()                                                  # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            IndexError: the number of names must equal the number of generators
        """
        if name is not None and names is not None:
            raise ValueError("you must specify the name of the variable")
        if names is None:
            if name is None:
                name = self._center_variable_name
            if name is None:
                name = 'z'
            names = (name,)
        names = normalize_names(1, names)
        name = names[0]
        if name in self._center:
            center = self._center[name]
        else:
            center = PolynomialRing(self._constants, names)
            embed = SkewPolynomialCenterInjection(center, self, self._embed_constants, self._order)
            try:
                assert not self.has_coerce_map_from(center)
                self.register_coercion(embed)
                center.register_conversion(embed.section())
            except AssertionError:
                raise ValueError("creation of coercion map fails; consider using another variable name")
            self._center[name] = center
        if default or (self._center_variable_name is None):
            self._center_variable_name = name
        return center


# Special class for skew polynomial over finite fields
######################################################

class SkewPolynomialRing_finite_field(SkewPolynomialRing_finite_order):
    r"""
    A specialized class for skew polynomial rings over finite fields.

    .. SEEALSO::

        - :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`
        - :mod:`sage.rings.polynomial.skew_polynomial_finite_field`

    .. TODO::

        Add methods related to center of skew polynomial ring, irreducibility, karatsuba
        multiplication and factorization.
    """
    def __init__(self, base_ring, morphism, derivation, names, sparse, category=None):
        r"""
        This method is a constructor for a general, dense univariate skew polynomial ring
        over a finite field.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to be used in ring

        ..NOTE::

            Multivariate and Sparse rings are not implemented.

        EXAMPLES::

            sage: k.<t> = GF(5^3)                                                       # needs sage.rings.finite_rings
            sage: Frob = k.frobenius_endomorphism()                                     # needs sage.rings.finite_rings
            sage: T.<x> = k['x', Frob]; T                                               # needs sage.rings.finite_rings
            Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        if self.Element is None:
            import sage.rings.polynomial.skew_polynomial_finite_field
            self.Element = sage.rings.polynomial.skew_polynomial_finite_field.SkewPolynomial_finite_field_dense
        SkewPolynomialRing_finite_order.__init__(self, base_ring, morphism, derivation, names, sparse, category)
        self._matrix_retraction = None

    def _new_retraction_map(self, seed=None):
        r"""
        Create a retraction map from the ring of coefficient
        of this skew polynomial ring to its fixed subfield under
        the twisting morphism

        This is an internal function used in factorization.

        INPUT:

        - ``seed`` -- an element of the base ring or ``None``
          (default: ``None``); it ``None``, a random element
          is picked

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(11^4)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: S._new_retraction_map()
            sage: S._matrix_retraction   # random
            [ 9  4 10  4]

        We can specify a seed::

            sage: # needs sage.rings.finite_rings
            sage: S._new_retraction_map(seed=a)
            sage: S._matrix_retraction
            [ 0  6  3 10]
            sage: S._new_retraction_map(seed=a)
            sage: S._matrix_retraction
            [ 0  6  3 10]
        """
        k = self.base_ring()
        section = self._embed_constants.section()
        if seed is None:
            seed = k.random_element()
        self._seed_retraction = seed
        trace = [ ]
        elt = seed
        for _ in range(k.degree()):
            x = elt
            tr = elt
            for _ in range(1, self._order):
                x = self._morphism(x)
                tr += x
            elt *= k.gen()
            trace.append(section(tr))
        self._matrix_retraction = MatrixSpace(self._constants, 1, k.degree())(trace)

    def _retraction(self, x, newmap=False, seed=None):
        r"""
        Return the image of `x` under the retraction map
        (see also :meth:`_new_retraction_map`)

        This is an internal function used in factorization.

        INPUT:

        - ``newmap`` -- a boolean (default: ``False``); whether we
          first create and use a new retraction map

        - ``seed`` -- an element of the base ring or ``None`` (default:
          ``None``); if given, first create a new random retraction map
          with given seed

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(11^4)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: S._retraction(a)   # random
            6

        Note that a retraction map has been automatically created::

            sage: S._matrix_retraction   # random                                       # needs sage.rings.finite_rings
            [ 0  6  3 10]

        If we call again the method :meth:`_retraction`,
        the same retraction map is used::

            sage: S._retraction(a)   # random                                           # needs sage.rings.finite_rings
            6

        We can specify a seed::

            sage: S._retraction(a^2, seed=a)  # random                                  # needs sage.rings.finite_rings
            10
        """
        # Better to return the retraction map but more difficult
        if newmap or seed is not None or self._matrix_retraction is None:
            self._new_retraction_map()
        return (self._matrix_retraction*self.base_ring()(x)._vector_())[0]
