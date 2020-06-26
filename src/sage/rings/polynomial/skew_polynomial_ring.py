r"""
Skew Univariate Polynomial Rings

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`
which constructs a general dense skew univariate polynomials over commutative base rings with
automorphisms over the base rings. This is the set of formal polynomials where the coefficients
are written on the left of the variable of the skew polynomial ring. The modified multiplication
operation over elements of the base ring is extended to all elements of the skew polynomial ring
by associativity and distributivity.

This module also provides :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order`
which is a specialized class for skew polynomial rings over fields equipped with an automorphism of
finite order. It inherits from
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing` but contains more
methods and provides better algorithms.

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

import sage

from sage.structure.richcmp import op_EQ
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.structure.category_object import normalize_names

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Algebra, Field
from sage.rings.integer import Integer

from sage.categories.commutative_rings import CommutativeRings
from sage.categories.algebras import Algebras
from sage.categories.fields import Fields

from sage.categories.morphism import Morphism, IdentityMorphism
from sage.rings.morphism import RingHomomorphism
from sage.categories.homset import Hom
from sage.categories.map import Section

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection

WORKING_CENTER_MAX_TRIES = 1000


# Helper functions

def _base_ring_to_fraction_field(S):
    """
    Return the unique skew polynomial ring over the fraction field of
    ``S.base_ring()`` which has ``S`` a sub-ring (internal method).

    INPUT:

    - ``S`` -- a skew polynomial ring.

    OUTPUT:

    - ``Q`` -- the skew polynomial ring over the fraction field of
      ``S.base_ring``.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _base_ring_to_fraction_field
        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x', sigma]
        sage: _base_ring_to_fraction_field(S)
        Skew Polynomial Ring in x over Fraction Field of Univariate Polynomial Ring in t over Integer Ring twisted by t |-->  t + 1
    """
    R = S.base_ring()
    if isinstance(R, Field):
        return S
    else:
        Q = R.fraction_field()
        gens = R.gens()
        sigmaS = S.twist_map()
        # try:
        sigmaQ = Q.hom([Q(sigmaS(g)) for g in gens])
        return Q[S.variable_name(), sigmaQ]
        # except Exception, e:
        #     raise ValueError("unable to lift the twist map to a twist map over %s (error was: %s)" % (Q, e))


def _minimal_vanishing_polynomial(R, eval_pts):
    """
    Return the minimal vanishing polynomial (internal function).

    See the documentation for
    :meth:`SkewPolynomialRing.minimal_vanishing_polynomial` for a description.

    INPUT:

    - ``R`` -- A skew polynomial ring over a field.

    - ``eval_pts`` -- a list of evaluation points

    OUTPUT:

    The minimal vanishing polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _minimal_vanishing_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [1, t, t^2]
        sage: b = _minimal_vanishing_polynomial(S, eval_pts); b
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/13215 for details.
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
            return R.gen() - R.twist_map()(e) / e
    else:
        t = l // 2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        B_moved = M_A.multi_point_evaluation(B)
        M_at_B_moved = _minimal_vanishing_polynomial(R, B_moved)
        return M_at_B_moved * M_A


def _lagrange_polynomial(R, eval_pts, values):
    """
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

    - the Lagrange polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_ring import _lagrange_polynomial
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [ t , t^2 ]
        sage: values = [ 3*t^2 + 4*t + 4 , 4*t ]
        sage: d = _lagrange_polynomial(S, eval_pts, values); d
        x + t
        sage: d.multi_point_evaluation(eval_pts) == values
        True

    The following restrictions are impossible to satisfy because the evaluation
    points are linearly dependent over the fixed field of the twist map, and the
    corresponding values do not match::

        sage: eval_pts = [ t, 2*t ]
        sage: values = [ 1, 3 ]
        sage: _lagrange_polynomial(S, eval_pts, values)
        Traceback (most recent call last):
        ...
        ValueError: the given evaluation points are linearly dependent over the fixed field of the twist map,
        so a Lagrange polynomial could not be determined (and might not exist).
    """
    l = len(eval_pts)
    if l == 1:
        if eval_pts[0].is_zero():
            # This is due to linear dependence among the eval_pts.
            raise ValueError("the given evaluation points are linearly dependent over the fixed field of the twist map, so a Lagrange polynomial could not be determined (and might not exist).")
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

class SkewPolynomialRing(Algebra, UniqueRepresentation):
    r"""
    Construct and return the globally unique skew polynomial ring with the
    given properties and variable names.

    Given a ring `R` and a ring automorphism `\sigma` of `R`, the ring of
    skew polynomials `R[X, \sigma]` is the usual abelian group polynomial
    `R[X]` equipped with the modification multiplication deduced from the
    rule `X a = \sigma(a) X`.
    We refer to [Ore1933]_ for more material on skew polynomials.

    .. SEEALSO::

        - :class:`sage.rings.polynomial.skew_polynomial_element.SkewPolynomial`

    INPUT:

    - ``base_ring`` -- a commutative ring

    - ``twist_map`` -- an automorphism of the base ring

    - ``names`` -- a string or a list of strings

    - ``sparse`` -- a boolean (default: ``False``). Currently not supported.

    .. NOTE::

        The current implementation of skew polynomial rings does not
        support derivations. Sparse skew polynomials and multivariate skew
        polynomials are also not implemented.

    OUTPUT:

    A univariate skew polynomial ring over ``base_ring`` twisted by
    ``twist_map`` when ``names`` is a string with no
    commas (``,``) or a list of length 1. Otherwise we raise a
    ``NotImplementedError`` as multivariate skew polynomial rings are
    not yet implemented.

    UNIQUENESS and IMMUTABILITY:

    In Sage, there is exactly one skew polynomial ring for each
    triple (base ring, twisting map, name of the variable).

    EXAMPLES of VARIABLE NAME CONTEXT::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = SkewPolynomialRing(R, sigma); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    The names of the variables defined above cannot be arbitrarily
    modified because each skew polynomial ring is unique in Sage and other
    objects in Sage could have pointers to that skew polynomial ring.

    However, the variable can be changed within the scope of a ``with``
    block using the localvars context::

        sage: with localvars(S, ['y']):
        ....:     print(S)
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    SQUARE BRACKETS NOTATION:

    You can alternatively create a skew polynomial ring over `R`
    twisted by ``twist_map`` by writing
    ``R['varname', twist_map]``.

    EXAMPLES:

    We first define the base ring::

        sage: R.<t> = ZZ[]; R
        Univariate Polynomial Ring in t over Integer Ring

    and the twisting map::

        sage: twist_map = R.hom([t+1]); twist_map
        Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
          Defn: t |--> t + 1

    Now, we are ready to define the skew polynomial ring::

        sage: S = SkewPolynomialRing(R, twist_map, names='x'); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    Use the diamond brackets notation to make the variable ready
    for use after you define the ring::

        sage: S.<x> = SkewPolynomialRing(R, twist_map)
        sage: (x + t)^2
        x^2 + (2*t + 1)*x + t^2

    Here is an example with the square bracket notations::

        sage: S.<x> = R['x', twist_map]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    Rings with different variables names are different::

        sage: R['x', twist_map] == R['y', twist_map]
        False

    Of course, skew polynomial rings with different twist maps are not
    equal either::

        sage: R['x',sigma] == R['x',sigma^2]
        False

    TESTS:

    You must specify a variable name::

        sage: SkewPolynomialRing(R, twist_map)
        Traceback (most recent call last):
        ...
        TypeError: you must specify the name of the variable

    With this syntax, it is not possible to omit the name of the
    variable neither in LHS nor in RHS. If we omit it in LHS, the
    variable is not created::

        sage: Sy = R['y', twist_map]; Sy
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1
        sage: y.parent()
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined

    If we omit it in RHS, sage tries to create a polynomial ring and fails::

        sage: Sz.<z> = R[twist_map]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n  Defn: t |--> t + 1' is not alphanumeric

    Multivariate skew polynomial rings are not supported::

        sage: S = SkewPolynomialRing(R, twist_map,names=['x','y'])
        Traceback (most recent call last):
        ...
        NotImplementedError: multivariate skew polynomials rings not supported

    Sparse skew polynomial rings are not implemented::

        sage: S = SkewPolynomialRing(R, twist_map, names='x', sparse=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: sparse skew polynomial rings are not implemented

    Saving and loading of polynomial rings works::

        sage: loads(dumps(R['x',sigma])) == R['x',sigma]
        True

    .. TODO::

        - Sparse Skew Polynomial Ring
        - Multivariate Skew Polynomial Ring
        - Add derivations.
    """
    Element = sage.rings.polynomial.skew_polynomial_element.SkewPolynomial_generic_dense

    @staticmethod
    def __classcall_private__(cls, base_ring, twist_map=None, names=None, sparse=False):
        r"""
        Construct the skew polynomial ring associated to the given parameters

        TESTS::

            sage: k.<t> = ZZ[]
            sage: theta = k.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(k, theta)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: type(S)
            <class 'sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_with_category'>

        We check the uniqueness property of parents::

            sage: sigma = k.hom([t+1])
            sage: T.<x> = SkewPolynomialRing(k, sigma)
            sage: S is T
            True

        When the twisting morphism has finite order, a special class
        is used::

            sage: k.<a> = GF(7^5)
            sage: Frob = k.frobenius_endomorphism(2)
            sage: S.<x> = SkewPolynomialRing(k, Frob)
            sage: type(S)
            <class 'sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_order_with_category'>
        """
        if base_ring not in CommutativeRings():
            raise TypeError('base_ring must be a commutative ring')
        if twist_map is None:
            twist_map = IdentityMorphism(base_ring)
        else:
            if (not isinstance(twist_map, Morphism)
                or twist_map.domain() is not base_ring
                or twist_map.codomain() is not base_ring):
                raise TypeError("the twist map must be a ring automorphism of base_ring (=%s)" % base_ring)
        if sparse:
            raise NotImplementedError("sparse skew polynomial rings are not implemented")
        if names is None:
            raise TypeError("you must specify the name of the variable")
        try:
            names = normalize_names(1, names)[0]
        except IndexError:
            raise NotImplementedError("multivariate skew polynomials rings not supported")

        # We check if the twisting morphism has finite order
        if base_ring in Fields():
            try:
                order = twist_map.order()
                if order is not Infinity:
                    return SkewPolynomialRing_finite_order(base_ring, twist_map, names, sparse)
            except (AttributeError, NotImplementedError):
                pass

        # We fallback to generic implementation
        return cls.__classcall__(cls, base_ring, twist_map, names, sparse)

    def __init__(self, base_ring, twist_map, name, sparse, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``twist_map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``category`` -- a category

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.category()
            Category of algebras over Univariate Polynomial Ring in t over Integer Ring
            sage: S([1]) + S([-1])
            0
            sage: TestSuite(S).run()
        """
        self.__is_sparse = sparse
        self._map = twist_map
        self._maps = {0: IdentityMorphism(base_ring), 1: self._map}
        category = Algebras(base_ring).or_subcategory(category)
        Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)

    def _element_constructor_(self, a=None, check=True, construct=False, **kwds):
        r"""
        Convert a base ring element ``a`` into a constant of this univariate
        skew polynomial ring, possibly non-canonically.

        INPUT:

        - ``a`` -- (default: ``None``) an element of the base ring
          of ``self`` or a ring that has a coerce map from ``self``

        - ``check`` -- boolean (default: ``True``)

        - ``construct`` -- boolean (default: ``False``)

        OUTPUT:

        An zero-degree skew polynomial in ``self``, equal to ``a``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S(1 + x + x^2 + x^3)
            x^3 + x^2 + x + 1
            sage: S(1 + t)
            t + 1
            sage: S(1 + t).degree()
            0
            sage: S(0).list()
            []

        TESTS::

            sage: S(x, check=True)
            x
        """
        C = self.Element
        if isinstance(a, list):
            return C(self, a, check=check, construct=construct)
        if isinstance(a, sage.structure.element.Element):
            P = a.parent()

            def build(check):
                if a.is_zero():
                    return P.zero()
                else:
                    return C(self, [a], check=check, construct=construct)
            if P is self:
                return a
            elif P is self.base_ring():
                build(False)
            elif P == self.base_ring() or self.base_ring().has_coerce_map_from(P):
                build(True)
        try:
            return a._polynomial_(self)
        except AttributeError:
            pass
        if isinstance(a, str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(a))
            except NameError:
                raise TypeError("unable to coerce string")
        return C(self, a, check, construct=construct, **kwds)

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: S.<x> = SkewPolynomialRing(R, R.hom([t + 1]))
            sage: S.coerce_map_from(R)
            Skew Polynomial base injection morphism:
              From: Univariate Polynomial Ring in t over Integer Ring
              To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: x.parent()
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: t.parent()
            Univariate Polynomial Ring in t over Integer Ring
            sage: y = x + t  # indirect doctest
            sage: y
            x + t
            sage: y.parent() is S
            True
        """
        return SkewPolynomialBaseringInjection(self.base_ring(), self)

    def _coerce_map_from_(self, P):
        r"""
        Check whether ``self`` has a coerce map from ``P``.

        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

        INPUT:

        - ``P`` -- a ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.has_coerce_map_from(S)
            True
            sage: S.has_coerce_map_from(R)
            True
            sage: S.has_coerce_map_from(ZZ)
            True
            sage: S.has_coerce_map_from(GF(5^3))
            False

            sage: S.coerce_map_from(ZZ)
            Composite map:
                From: Integer Ring
                To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
                Defn:   Polynomial base injection morphism:
                        From: Integer Ring
                        To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                        Skew Polynomial base injection morphism:
                        From: Univariate Polynomial Ring in t over Integer Ring
                        To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S.coerce_map_from(S)
            Identity endomorphism of Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        base_ring = self.base_ring()
        try:
            connecting = base_ring.coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(base_ring) * connecting
        except TypeError:
            pass
        if isinstance(P, SkewPolynomialRing):
            if self.__is_sparse and not P.is_sparse():
                return False
            if P.variable_name() == self.variable_name():
                return base_ring.has_coerce_map_from(P.base_ring())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        s = "Skew Polynomial Ring in %s over %s twisted by %s" % (self.variable_name(),
                                                                  self.base_ring(),
                                                                  self._map._repr_short())
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: latex(S)
            \Bold{Z}[t][x,\begin{array}{l}
            \text{\texttt{Ring{ }endomorphism...}}
            \end{array}]
        """
        from sage.misc.latex import latex
        return "%s[%s,%s]" % (latex(self.base_ring()), self.latex_variable_names()[0],
                              latex(self._map))

    def change_var(self, var):
        r"""
        Return the skew polynomial ring in variable ``var`` with the same base
        ring and twist map as ``self``.

        INPUT:

        - ``var`` -- a string representing the name of the new variable.

        OUTPUT:

        ``self`` with variable name changed to ``var``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = SkewPolynomialRing(k,Frob); R
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry = R.change_var('y'); Ry
            Skew Polynomial Ring in y over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry is R.change_var('y')
            True
        """
        return SkewPolynomialRing(self.base_ring(), self._map, names=var,
                                  sparse=self.__is_sparse)

    def characteristic(self):
        r"""
        Return the characteristic of the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: R['x',sigma].characteristic()
            0

            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: k['y',Frob].characteristic()
            5
        """
        return self.base_ring().characteristic()

    @cached_method
    def twist_map(self, n=1):
        r"""
        Return the twist map, the automorphism of the base ring of
        ``self``, iterated ``n`` times.

        INPUT:

        -  ``n`` - an integer (default: 1)

        OUTPUT:

        ``n``-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.twist_map()
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1
            sage: S.twist_map() == sigma
            True
            sage: S.twist_map(10)
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 10

        If ``n`` in negative, Sage tries to compute the inverse of the
        twist map::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

        Sometimes it fails, even if the twist map is actually invertible::

            sage: S.twist_map(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1
        """
        try:
            return self._map ** n
        except TypeError as e:
            if n < 0:
                raise NotImplementedError("inversion of the twist map %s" % self._map)
            else:
                raise ValueError("Unexpected error in iterating the twist map: %s", e)

    @cached_method
    def gen(self, n=0):
        r"""
        Return the indeterminate generator of this skew polynomial ring.

        INPUT:

        - ``n`` -- index of generator to return (default: 0). Exists for
          compatibility with other polynomial rings.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            sage: y = S.gen(); y
            x
            sage: y == x
            True
            sage: y is x
            True
            sage: S.gen(0)
            x

        This is also known as the parameter::

            sage: S.parameter() is S.gen()
            True
        """
        if n != 0:
            raise IndexError("generator %s not defined" % n)
        return self.Element(self, [0, 1])

    parameter = gen

    def gens_dict(self):
        r"""
        Return a {name: variable} dictionary of the generators of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self):
        r"""
        Return ``False`` since skew polynomial rings are not finite
        (unless the base ring is `0`.)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: k.is_finite()
            True
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_finite()
            False
        """
        R = self.base_ring()
        return R.is_finite() and R.order() == 1

    def is_exact(self):
        r"""
        Return ``True`` if elements of this skew polynomial ring are exact.
        This happens if and only if elements of the base ring are exact.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_exact()
            True
            sage: S.base_ring().is_exact()
            True

            sage: R.<u> = k[[]]
            sage: sigma = R.hom([u+u^2])
            sage: T.<y> = R['y',sigma]
            sage: T.is_exact()
            False
            sage: T.base_ring().is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_sparse(self):
        r"""
        Return ``True`` if the elements of this polynomial ring are sparsely
        represented.

        .. WARNING::

            Since sparse skew polynomials are not yet implemented, this
            function always returns ``False``.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.is_sparse()
            False
        """
        return self.__is_sparse

    def ngens(self):
        r"""
        Return the number of generators of this skew polynomial ring,
        which is 1.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, monic=False, *args, **kwds):
        r"""
        Return a random skew polynomial in ``self``.

        INPUT:

        - ``degree`` -- (default: 2) integer with degree
          or a tuple of integers with minimum and maximum degrees

        - ``monic`` -- (default: ``False``) if ``True``, return a monic
          skew polynomial

        - ``*args, **kwds`` -- passed on to the ``random_element`` method
          for the base ring

        OUTPUT:

        Skew polynomial such that the coefficients of `x^i`, for `i` up
        to ``degree``, are random elements from the base ring, randomized
        subject to the arguments ``*args`` and ``**kwds``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: S.random_element()  # random
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: S.random_element(monic=True)  # random
            x^2 + (2*t^2 + t + 1)*x + 3*t^2 + 3*t + 2

        Use ``degree`` to obtain polynomials of higher degree

            sage: p = S.random_element(degree=5)   # random
            (t^2 + 3*t)*x^4 + (4*t + 4)*x^3 + (4*t^2 + 4*t)*x^2 + (2*t^2 + 1)*x + 3

        When ``monic`` is ``False``, the returned skew polynomial may have
        a degree less than ``degree`` (it happens when the random leading
        coefficient is zero). However, if ``monic`` is ``True``, this can't
        happen::

            sage: p = S.random_element(degree=4, monic=True)
            sage: p.leading_coefficient() == S.base_ring().one()
            True
            sage: p.degree() == 4
            True

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree, both inclusive::

            sage: S.random_element(degree=(2,7))  # random
            (3*t^2 + 1)*x^4 + (4*t + 2)*x^3 + (4*t + 1)*x^2
             + (t^2 + 3*t + 3)*x + 3*t^2 + 2*t + 2

        If the first tuple element is greater than the second, a a
        ``ValueError`` is raised::

            sage: S.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: first degree argument must be less or equal to the second
        """
        R = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("first degree argument must be less or equal to the second")
            degree = randint(*degree)
        if monic:
            return self([R.random_element(*args, **kwds)
                         for _ in range(degree)] + [R.one()])
        else:
            return self([R.random_element(*args, **kwds)
                         for _ in range(degree + 1)])

    def is_commutative(self):
        r"""
        Return ``True`` if this skew polynomial ring is commutative, i.e. if the
        twist map is the identity.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_commutative()
            False

            sage: T.<y> = k['y',Frob^3]
            sage: T.is_commutative()
            True
        """
        return self.twist_map().is_identity()

    def minimal_vanishing_polynomial(self, eval_pts):
        """
        Return the minimal-degree, monic skew polynomial which vanishes at all
        the given evaluation points.

        The degree of the vanishing polynomial is at most the length of
        ``eval_pts``. Equality holds if and only if the elements of ``eval_pts``
        are linearly independent over the fixed field of ``self.twist_map()``.

        - ``eval_pts`` -- list of evaluation points which are linearly
          independent over the fixed field of the twist map of the associated
          skew polynomial ring

        OUTPUT:

        The minimal vanishing polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: eval_pts = [1, t, t^2]
            sage: b = S.minimal_vanishing_polynomial(eval_pts); b
            x^3 + 4

        The minimal vanishing polynomial evaluates to 0 at each of the evaluation points::

            sage: eval = b.multi_point_evaluation(eval_pts); eval
            [0, 0, 0]

        If the evaluation points are linearly dependent over the fixed field of
        the twist map, then the returned polynomial has lower degree than the
        number of evaluation points::

            sage: S.minimal_vanishing_polynomial([t])
            x + 3*t^2 + 3*t
            sage: S.minimal_vanishing_polynomial([t, 3*t])
            x + 3*t^2 + 3*t
        """
        return _minimal_vanishing_polynomial(_base_ring_to_fraction_field(self), eval_pts)

    def lagrange_polynomial(self, points):
        r"""
        Return the minimal-degree polynomial which interpolates the given
        points.

        More precisely, given `n` pairs `(x_1, y_1), ..., (x_n, y_n) \in R^2`,
        where `R` is ``self.base_ring()``, compute a skew polynomial `p(x)` such
        that `p(x_i) = y_i` for each `i`, under the condition that the `x_i` are
        linearly independent over the fixed field of ``self.twist_map()``.

        If the `x_i` are linearly independent over the fixed field of
        ``self.twist_map()`` then such a polynomial is guaranteed to exist.
        Otherwise, it might exist depending on the `y_i`, but the algorithm used
        in this implementation does not support that, and so an error is always
        raised.

        INPUT:

        - ``points`` -- a list of pairs ``(x_1, y_1),..., (x_n, y_n)`` of
          elements of the base ring of ``self``. The `x_i` should be linearly
          independent over the fixed field of ``self.twist_map()``.

        OUTPUT:

        The Lagrange polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: points = [(t, 3*t^2 + 4*t + 4), (t^2, 4*t)]
            sage: d = S.lagrange_polynomial(points); d
            x + t

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: T.<x> = R['x', sigma]
            sage: points = [ (1, t^2 + 3*t + 4), (t, 2*t^2 + 3*t + 1), (t^2, t^2 + 3*t + 4) ]
            sage: p = T.lagrange_polynomial(points); p
            ((-t^4 - 2*t - 3)/-2)*x^2 + (-t^4 - t^3 - t^2 - 3*t - 2)*x + (-t^4 - 2*t^3 - 4*t^2 - 10*t - 9)/-2
            sage: p.multi_point_evaluation([1, t, t^2]) == [ t^2 + 3*t + 4, 2*t^2 + 3*t + 1, t^2 + 3*t + 4 ]
            True

        If the `x_i` are linearly dependent over the fixed field of
        ``self.twist_map()``, then an error is raised::

            sage: T.lagrange_polynomial([ (t, 1), (2*t, 3) ])
            Traceback (most recent call last):
            ...
            ValueError: the given evaluation points are linearly dependent over the fixed field of the twist map,
            so a Lagrange polynomial could not be determined (and might not exist).
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
            raise TypeError("a skew polynomial always evaluates to 0 at 0, but a non-zero value was requested.")

        return _lagrange_polynomial(_base_ring_to_fraction_field(self), eval_pts, values)


# Special classes for twisting morphisms with finite order
##########################################################

class SectionSkewPolynomialCenterInjection(Section):
    r"""
    Section of the canonical injection of the center of a skew
    polynomial ring into this ring

    TESTS::

        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.center()
        sage: iota = S.convert_map_from(Z)
        sage: sigma = iota.section()
        sage: TestSuite(sigma).run(skip=['_test_category'])
    """
    def _call_(self, x):
        r"""
        Return `x` viewed as an element of the center

        EXAMPLES::

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
        Compare this morphism with ``right``

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: sigma = iota.section()

            sage: s = loads(dumps(sigma))
            sage: s == sigma
            True
            sage: s is sigma
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented


class SkewPolynomialCenterInjection(RingHomomorphism):
    r"""
    Canonical injection of the center of a skew polynomial ring
    into this ring

    TESTS::

        sage: k.<a> = GF(5^3)
        sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
        sage: Z = S.center()
        sage: iota = S.convert_map_from(Z)
        sage: TestSuite(iota).run(skip=['_test_category'])
    """
    def __init__(self, domain, codomain, embed, order):
        r"""
        Initialize this morphism

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: S.convert_map_from(Z)   # indirect doctest
            Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
        """
        RingHomomorphism.__init__(self, Hom(domain, codomain))
        self._embed = embed
        self._order = order
        self._codomain = codomain
        self._section = SectionSkewPolynomialCenterInjection(self)

    def _repr_(self):
        r"""
        Return a string representation of this morphism

        EXAMPLES::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)
            sage: iota
            Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring
            sage: iota._repr_()
            'Embedding of the center of Skew Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this ring'
        """
        return "Embedding of the center of %s into this ring" % self._codomain

    def _call_(self, x):
        r"""
        Return the image of `x` by this morphism

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z.<z> = S.center()
            sage: iota = S.convert_map_from(Z)

            sage: iota(z)
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
        Compare this morphism with ``right``

        TESTS::

            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: Z = S.center()
            sage: iota = S.convert_map_from(Z)

            sage: i = loads(dumps(iota))
            sage: i == iota
            True
            sage: i is iota
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented

    def section(self):
        r"""
        Return a section of this morphism

        EXAMPLES::

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
    A specialized class for skew polynomial rings over finite fields.

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing`
        :mod:`sage.rings.polynomial.skew_polynomial_finite_order`
    """
    import sage.rings.polynomial.skew_polynomial_finite_order
    Element = sage.rings.polynomial.skew_polynomial_finite_order.SkewPolynomial_finite_order_dense

    def __init__(self, base_ring, twist_map, name, sparse, category=None):
        r"""
        Initialize this skew polynomial

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]; S
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: S.category()
            Category of algebras over Finite Field in t of size 5^3

            sage: TestSuite(S).run()

        We check that a call to the method
        :meth:`sage.rings.polynomial.skew_polynomial_finite_order.SkewPolynomial_finite_order.is_central`
        does not affect the behaviour of default central variable names::

            sage: k.<a> = GF(7^4)
            sage: phi = k.frobenius_endomorphism()
            sage: S.<x> = k['x', phi]
            sage: (x^4).is_central()
            True
            sage: Z.<u> = S.center()
            sage: S.center() is Z
            True
        """
        SkewPolynomialRing.__init__(self, base_ring, twist_map, name, sparse, category)
        self._order = twist_map.order()
        (self._constants, self._embed_constants) = twist_map.fixed_field()

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

            If F denotes the subring of R fixed by `\sigma` and `\sigma`
            has order `r`, the center of `K[x,\sigma]` is `F[x^r]`, that
            is a univariate polynomial ring over `F`.

        INPUT:

        - ``name`` -- a string or ``None`` (default: ``None``);
          the name for the central variable (namely `x^r`)

        - ``default`` -- a boolean (default: ``False``); if ``True``,
          set the default variable name for the center to ``name``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]; S
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

            sage: Z = S.center(); Z
            Univariate Polynomial Ring in z over Finite Field of size 5
            sage: Z.gen()
            z

        We can pass in another variable name::

            sage: S.center(name='y')
            Univariate Polynomial Ring in y over Finite Field of size 5

        or use the bracket notation::

            sage: Zy.<y> = S.center(); Zy
            Univariate Polynomial Ring in y over Finite Field of size 5
            sage: y.parent() is Zy
            True

        A coercion map from the center to the skew polynomial ring is set::

            sage: S.has_coerce_map_from(Zy)
            True

            sage: P = y + x; P
            x^3 + x
            sage: P.parent()
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: P.parent() is S
            True

        together with a conversion map in the reverse direction::

            sage: Zy(x^6 + 2*x^3 + 3)
            y^2 + 2*y + 3

            sage: Zy(x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in the center

        Two different skew polynomial rings can share the same center::

            sage: S1.<x1> = k['x1', Frob]
            sage: S2.<x2> = k['x2', Frob]
            sage: S1.center() is S2.center()
            True

        ABOUT THE DEFAULT NAME OF THE CENTRAL VARIABLE:

        A priori, the default is ``z``.

        However, a variable name is given the first time this method is
        called, the given name become the default for the next calls::

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

            sage: D.<v> = A.center(default=True)
            sage: D
            Univariate Polynomial Ring in v over Finite Field of size 11
            sage: A.center()
            Univariate Polynomial Ring in v over Finite Field of size 11
            sage: A.center() is D
            True

        TESTS::

            sage: C.<a,b> = S.center()
            Traceback (most recent call last):
            ...
            IndexError: the number of names must equal the number of generators
        """
        if name is not None and names is not None:
            raise ValueError
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
