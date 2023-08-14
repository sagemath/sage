r"""
Wehler K3 Surfaces

AUTHORS:

- Ben Hutz (11-2012)
- Joao Alberto de Faria (10-2013)

.. TODO::

    Hasse-Weil Zeta Function

    Picard Number

    Number Fields

REFERENCES: [FH2015]_, [CS1996]_, [Weh1998]_, [Hutz2007]

"""
# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.calculus.functions import jacobian
from sage.categories.fields import Fields
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.number_fields import NumberFields
from sage.misc.functional import sqrt
from sage.misc.cachefunc import cached_method
from sage.misc.mrange import xmrange
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.order import is_NumberFieldOrder
from sage.rings.padics.factory import Qp
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.schemes.product_projective.subscheme import AlgebraicScheme_subscheme_product_projective
from sage.schemes.product_projective.space import ProductProjectiveSpaces
from copy import copy

_NumberFields = NumberFields()
_Fields = Fields()
_CommutativeRings = CommutativeRings()


def WehlerK3Surface(polys):
    r"""
    Defines a K3 Surface over `\mathbb{P}^2 \times \mathbb{P}^2` defined as
    the intersection of a bilinear and biquadratic form. [Weh1998]_

    INPUT: Bilinear and biquadratic polynomials as a tuple or list

    OUTPUT: :class:`WehlerK3Surface_ring`

    EXAMPLES::

        sage: PP.<x0,x1, x2, y0, y1, y2> = ProductProjectiveSpaces([2, 2],QQ)
        sage: L = x0*y0 + x1*y1 - x2*y2
        sage: Q = x0*x1*y1^2 + x2^2*y0*y2
        sage: WehlerK3Surface([L, Q])
        Closed subscheme of Product of projective spaces P^2 x P^2 over Rational
        Field defined by:
         x0*y0 + x1*y1 - x2*y2,
         x0*x1*y1^2 + x2^2*y0*y2
    """
    if not isinstance(polys, (list, tuple)):
        raise TypeError("polys must be a list or tuple of polynomials")

    R = polys[0].parent().base_ring()
    if R in _Fields:
        if R in _Fields.Finite():
            return WehlerK3Surface_finite_field(polys)
        return WehlerK3Surface_field(polys)
    if R in _CommutativeRings:
        return WehlerK3Surface_ring(polys)
    raise TypeError("R (= %s) must be a commutative ring" % R)


def random_WehlerK3Surface(PP):
    r"""
    Produces a random K3 surface in `\mathbb{P}^2 \times \mathbb{P}^2` defined as the
    intersection of a bilinear and biquadratic form. [Weh1998]_

    INPUT: Projective space cartesian product

    OUTPUT: :class:`WehlerK3Surface_ring`

    EXAMPLES::

        sage: PP.<x0, x1, x2, y0, y1, y2> = ProductProjectiveSpaces([2, 2], GF(3))
        sage: w = random_WehlerK3Surface(PP)
        sage: type(w)
        <class 'sage.dynamics.arithmetic_dynamics.wehlerK3.WehlerK3Surface_finite_field_with_category'>
    """

    CR = PP.coordinate_ring()
    BR = PP.base_ring()
    Q = 0
    for a in xmrange([3,3]):
        for b in xmrange([3,3]):
            Q += BR.random_element() * CR.gen(a[0]) * CR.gen(a[1]) * CR.gen(3+b[0]) * CR.gen(3+b[1])
    #We can always change coordinates to make L diagonal
    L = CR.gen(0) * CR.gen(3) + CR.gen(1) * CR.gen(4) + CR.gen(2) * CR.gen(5)
    return WehlerK3Surface([L, Q])


class WehlerK3Surface_ring(AlgebraicScheme_subscheme_product_projective):
    r"""
    A K3 surface in `\mathbb{P}^2 \times \mathbb{P}^2` defined as the
    intersection of a bilinear and biquadratic form. [Weh1998]_

    EXAMPLES::

        sage: R.<x,y,z,u,v,w> = PolynomialRing(QQ, 6)
        sage: L = x*u - y*v
        sage: Q = x*y*v^2 + z^2*u*w
        sage: WehlerK3Surface([L, Q])
        Closed subscheme of Product of projective spaces P^2 x P^2 over Rational
        Field defined by:
          x*u - y*v,
          x*y*v^2 + z^2*u*w
    """
    def __init__(self, polys):
        if not isinstance(polys, (list, tuple)):
            raise TypeError("polys must be a list or tuple of polynomials")
        R = polys[0].parent()
        vars = R.variable_names()
        A = ProductProjectiveSpaces([2, 2],R.base_ring(),vars)
        CR = A.coordinate_ring()
        #Check for following:
        #    Is the user calling in 2 polynomials from a list or tuple?
        #    Is there one biquadratic and one bilinear polynomial?
        if len(polys) != 2:
            raise AttributeError("there must be 2 polynomials")

        if (all(((e[0] + e[1] + e[2]) == 1 and (e[3] + e[4] + e[5]) == 1) for e in polys[0].exponents())):
            self.L = CR(polys[0])
        elif (all(((e[0] + e[1] + e[2]) == 1 and (e[3] + e[4] + e[5]) == 1) for e in polys[1].exponents())):
            self.L = CR(polys[1])
        else:
            raise AttributeError("there must be one bilinear polynomial")

        if (all(((e[0] + e[1] + e[2]) == 2 and (e[3] + e[4] + e[5]) == 2) for e in polys[0].exponents())):
            self.Q = CR(polys[0])
        elif (all(((e[0] + e[1] + e[2]) == 2 and (e[3] + e[4] + e[5]) == 2) for e in polys[1].exponents())):
            self.Q = CR(polys[1])
        else:
            raise AttributeError("there must be one biquadratic polynomial")
        AlgebraicScheme_subscheme.__init__(self, A, polys)

    def change_ring(self, R):
        r"""
        Changes the base ring on which the Wehler K3 Surface is defined.

        INPUT: ``R`` - ring

        OUTPUT: K3 Surface defined over input ring

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], GF(3))
            sage: L = x0*y0 + x1*y1 - x2*y2
            sage: Q = x0*x1*y1^2 + x2^2*y0*y2
            sage: W = WehlerK3Surface([L, Q])
            sage: W.base_ring()
            Finite Field of size 3
            sage: T = W.change_ring(GF(7))
            sage: T.base_ring()
            Finite Field of size 7
        """

        LR = self.L.change_ring(R)
        LQ = self.Q.change_ring(R)
        return (WehlerK3Surface( [LR,LQ]))

    def _check_satisfies_equations(self, P):
        r"""
        Function checks to see if point ``P`` lies on the K3 Surface.

        INPUT: ``P`` - point in `\mathbb{P}^2 \times \mathbb{P}^2`

        OUTPUT: AttributeError True if the point is not on the surface

        EXAMPLES::

            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 \
            + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - \
            4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0 * y0 + x1 * y1 + x2 * y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._check_satisfies_equations([0, 0, 1, 1, 0, 0])
            True

        ::

            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 \
            + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - \
            4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._check_satisfies_equations([0, 1, 1, 1, 0, 0])
            Traceback (most recent call last):
            ...
            AttributeError: point not on surface
        """
        Point = list(P)
        if self.L(*Point) == 0 and self.Q(*Point) == 0:
            return True
        else:
            raise AttributeError("point not on surface")

    def _Lcoeff(self, component, i):
        r"""
        Return the polynomials  `L^x_i` or `L^y_i`.

        These polynomials are defined as:

        `L^x_i` = the coefficients of `y_i` in `L(x, y)` (Component = 0)

        `L^y_i` = the coefficients of `x_i` in `L(x, y)` (Component = 1)

        Definition and Notation from: [CS1996]_

        INPUT:

        - ``component`` - Integer: 0 or 1

        - ``i`` - Integer: 0, 1 or 2

        OUTPUT: polynomial in terms of either y (Component = 0) or x (Component = 1)

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._Lcoeff(0, 0)
            y0

        ::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z =x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._Lcoeff(1, 0)
            x0
        """
        # Error Checks for Passed in Values
        if component not in [0, 1]:
            raise ValueError("component can only be 1 or 0")
        if i not in [0, 1, 2]:
            raise ValueError("index must be 0, 1, or 2")
        R = self.ambient_space().coordinate_ring()
        return self.L.coefficient(R.gen(component*3 + i))

    def _Qcoeff(self, component, i, j):
        r"""
        Return the polynomials `Q^x_{ij}` or `Q^y_{ij}`.

        These polynomials are defined as:

        `Q^x_{ij}` = the coefficients of `y_{i}y_{j}` in `Q(x, y)` (Component = 0).

        `Q^y_{ij}` = the coefficients of `x_{i}x_{j}` in `Q(x, y)` (Component = 1).

        Definition and Notation from: [CS1996]_.

        INPUT:

        - ``component`` - Integer: 0 or 1

        - ``i`` - Integer: 0, 1 or 2

        - ``j`` - Integer: 0, 1 or 2

        OUTPUT: polynomial in terms of either y (Component = 0) or x (Component = 1)

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._Qcoeff(0, 0, 0)
            y0*y1 + y2^2

        ::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X._Qcoeff(1, 1, 0)
            x0^2
        """
        # Check Errors in Passed in Values
        if component not in [0, 1]:
            raise ValueError("component can only be 1 or 0")

        if i not in [0, 1, 2] or j not in [0, 1, 2]:
            raise ValueError("the two indexes must be either 0, 1, or 2")

        R = self.ambient_space().coordinate_ring()
        return self.Q.coefficient(R.gen(component * 3 + i) * R.gen(component * 3 + j))

    @cached_method
    def Gpoly(self, component, k):
        r"""
        Return the G polynomials  `G^*_k`.

        They are defined as:
        `G^*_k = \left(L^*_j\right)^2Q^*_{ii}-L^*_iL^*_jQ^*_{ij}+\left(L^*_i\right)^2Q^*_{jj}`\
        where {i, j, k} is some permutation of (0, 1, 2) and * is either
        x (Component = 1) or y (Component = 0).

        INPUT:

        - ``component`` - Integer: 0 or 1

        - ``k`` - Integer: 0, 1 or 2

        OUTPUT: polynomial in terms of either y (Component = 0) or x (Component = 1)

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.Gpoly(1, 0)
            x0^2*x1^2 + x1^4 - x0*x1^2*x2 + x1^3*x2 + x1^2*x2^2 + x2^4
        """
        # Check Errors in passed in values
        if component not in [0, 1]:
            raise ValueError("component can only be 1 or 0")

        if k not in [0, 1, 2]:
            raise ValueError("index must be either 0, 1, or 2")

        Indices = [0, 1, 2]
        Indices.remove( k)
        i = Indices[0]
        j = Indices[1]

        return (self._Lcoeff(component, j)**2) * (self._Qcoeff(component, i, i)) - (self._Lcoeff(component, i))* \
            (self._Lcoeff(component, j)) * (self._Qcoeff(component, i, j)) + (self._Lcoeff( component, i)**2)* \
            (self._Qcoeff( component, j, j))

    @cached_method
    def Hpoly(self, component, i, j):
        r"""
        Return the H polynomials defined as `H^*_{ij}`.

        This polynomial is defined by:

        `H^*_{ij} = 2L^*_iL^*_jQ^*_{kk}-L^*_iL^*_kQ^*_{jk} - L^*_jL^*_kQ^*_{ik}+\left(L^*_k\right)^2Q^*_{ij}`
        where {i, j, k} is some permutation of (0, 1, 2) and * is either y (Component = 0) or x (Component = 1).

        INPUT:

        - ``component`` - Integer: 0 or 1

        - ``i`` - Integer: 0, 1 or 2

        - ``j`` - Integer: 0, 1 or 2

        OUTPUT: polynomial in terms of either y (Component = 0) or x (Component = 1)

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 \
            + x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.Hpoly(0, 1, 0)
             2*y0*y1^3 + 2*y0*y1*y2^2 - y1*y2^3
        """
        #Check Errors in Passed in Values
        if component not in [0, 1]:
            raise ValueError("component can only be 1 or 0")

        if (i not in [0, 1, 2]) or (j not in [0, 1, 2]):
            raise ValueError("the two indexes must be either 0, 1, or 2")

        Indices = [0, 1, 2]
        Indices.remove(i)
        Indices.remove(j)

        k = Indices[0]

        return 2*(self._Lcoeff(component, i)) * (self._Lcoeff(component, j)) * (self._Qcoeff(component, k, k)) -\
             (self._Lcoeff(component, i)) * (self._Lcoeff( component, k)) * (self._Qcoeff(component, j, k)) -\
              (self._Lcoeff(component, j)) * (self._Lcoeff(component, k)) * (self._Qcoeff( component, i, k)) +\
               (self._Lcoeff(component, k)**2) * (self._Qcoeff(component, i, j))

    def Lxa(self, a):
        r"""
        Function will return the L polynomial defining the fiber, given by `L^{x}_{a}`.

        This polynomial is defined as:

        `L^{x}_{a} = \{(a, y) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon L(a, y) = 0\}`.

        Notation and definition from: [CS1996]_

        INPUT: ``a`` - Point in `\mathbb{P}^2`

        OUTPUT: A polynomial representing the fiber

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 \
            + 3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - \
            x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 \
            + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 \
            + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP(1, 1, 0, 1, 0, 0)
            sage: X.Lxa(T[0])
            y0 + y1
        """
        if a not in self.ambient_space()[0]:
            raise TypeError("point must be in projective space of dimension 2")
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[1]
        PSYC = PSY.coordinate_ring()
        #Define projection homomorphism
        p = ASC.hom([a[0],a[1],a[2]] + list(PSY.gens()), PSYC)
        return p(self.L)

    def Qxa(self, a):
        r"""
        Function will return the Q polynomial defining a fiber given by `Q^{x}_{a}`.

        This polynomial is defined as:

        `Q^{x}_{a} = \{(a,y) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon Q(a,y) = 0\}`.

        Notation and definition from: [CS1996]_

        INPUT: ``a`` - Point in `\mathbb{P}^2`

        OUTPUT: A polynomial representing the fiber

        EXAMPLES::

           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
           sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 \
           - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 \
           + 5*x0*x2*y0*y2 \
           - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
           sage: Y = x0*y0 + x1*y1 + x2*y2
           sage: X = WehlerK3Surface([Z, Y])
           sage: T = PP(1, 1, 0, 1, 0, 0)
           sage: X.Qxa(T[0])
           5*y0^2 + 7*y0*y1 + y1^2 + 11*y1*y2 + y2^2
        """
        if a not in self.ambient_space()[0]:
            raise TypeError("point must be in Projective Space of dimension 2")
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[1]
        PSYC = PSY.coordinate_ring()
        #Define projection homomorphism
        p = ASC.hom([a[0], a[1], a[2]] + list(PSY.gens()), PSYC)
        return p(self.Q)

    def Sxa(self, a):
        r"""
        Function will return fiber by `S^{x}_{a}`.

        This function is defined as:

        `S^{x}_{a} = L^{x}_{a} \cap Q^{x}_{a}`.

        Notation and definition from: [CS1996]_

        INPUT: ``a`` - Point in `\mathbb{P}^2`

        OUTPUT: A subscheme representing the fiber

        EXAMPLES::

           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
           sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 \
           + 3*x0*x1*y0*y1 \
           - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 \
           + 5*x0*x2*y0*y2 \
           - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
           sage: Y = x0*y0 + x1*y1 + x2*y2
           sage: Y = x0*y0 + x1*y1 + x2*y2
           sage: X = WehlerK3Surface([Z, Y])
           sage: T = PP(1, 1, 0, 1, 0, 0)
           sage: X.Sxa(T[0])
             Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              y0 + y1,
              5*y0^2 + 7*y0*y1 + y1^2 + 11*y1*y2 + y2^2
        """
        if a not in self.ambient_space()[0]:
            raise TypeError("point must be in projective space of dimension 2")
        PSY = self.ambient_space()[1]
        return PSY.subscheme([self.Lxa(a),self.Qxa(a)])

    def Lyb(self, b):
        r"""
        Function will return a fiber by `L^{y}_{b}`.

        This polynomial is defined as:

        `L^{y}_{b} = \{(x,b) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon L(x,b) = 0\}`.

        Notation and definition from: [CS1996]_

        INPUT: ``b`` - Point in projective space

        OUTPUT: A polynomial representing the fiber

        EXAMPLES::

           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
           sage: Z =x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 \
           + 3*x0*x1*y0*y1 \
           - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 \
           + 5*x0*x2*y0*y2 \
           - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
           sage: Y = x0*y0 + x1*y1 + x2*y2
           sage: Y = x0*y0 + x1*y1 + x2*y2
           sage: X = WehlerK3Surface([Z, Y])
           sage: T = PP(1, 1, 0, 1, 0, 0)
           sage: X.Lyb(T[1])
             x0
        """
        if b not in self.ambient_space()[1]:
            raise TypeError("point must be in projective space of dimension 2")
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[0]
        PSYC = PSY.coordinate_ring()
        p = ASC.hom(list(PSY.gens()) + [b[0], b[1], b[2]], PSYC)
        return (p(self.L))

    def Qyb(self, b):
        r"""

        Function will return a fiber by `Q^{y}_{b}`.

        This polynomial is defined as:

        `Q^{y}_{b} = \{(x,b) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon Q(x,b) = 0\}`.

        Notation and definition from: [CS1996]_

        INPUT: ``b`` - Point in projective space

        OUTPUT: A polynomial representing the fiber

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 \
            + 3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 \
            + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP(1, 1, 0, 1, 0, 0)
            sage: X.Qyb(T[1])
            x0^2 + 3*x0*x1 + x1^2
        """
        if b not in self.ambient_space()[1]:
            raise TypeError("point must be in projective space of dimension 2")
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[0]
        PSYC = PSY.coordinate_ring()
        p = ASC.hom(list(PSY.gens()) + [b[0], b[1], b[2]], PSYC)
        return (p(self.Q))

    def Syb(self, b):
        r"""
        Function will return fiber by `S^{y}_{b}`.

        This function is defined by:

        `S^{y}_{b} = L^{y}_{b} \cap Q^{y}_{b}`.

        Notation and definition from: [CS1996]_

        INPUT: ``b`` - Point in `\mathbb{P}^2`

        OUTPUT: A subscheme representing the fiber

        EXAMPLES::

           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
           sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
           3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
           - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 \
           + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
           sage: Y = x0 * y0 + x1 * y1 + x2 * y2
           sage: X = WehlerK3Surface([Z, Y])
           sage: T = PP(1, 1, 0, 1, 0, 0)
           sage: X.Syb(T[1])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x0,
              x0^2 + 3*x0*x1 + x1^2

        """
        if b not in self.ambient_space()[1]:
            raise TypeError("point must be in projective space of dimension 2")
        AS = self.ambient_space()
        PSY = AS[0]
        return PSY.subscheme([self.Lyb(b), self.Qyb(b)])

    def Ramification_poly(self, i):
        r"""
        Function will return the Ramification polynomial  `g^*`.

        This polynomial is defined by:

        `g^* = \frac{\left(H^*_{ij}\right)^2 - 4G^*_iG^*_j}{\left(L^*_k\right)^2}`.

        The roots of this polynomial will either be degenerate fibers or fixed points
        of the involutions `\sigma_x` or `\sigma_y` for more information, see [CS1996]_.

        INPUT: ``i`` - Integer, either 0 (polynomial in y) or 1 (polynomial in x)

        OUTPUT: Polynomial in the coordinate ring of the ambient space

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1\
            - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2\
            - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.Ramification_poly(0)
            8*y0^5*y1 - 24*y0^4*y1^2 + 48*y0^2*y1^4 - 16*y0*y1^5 + y1^6 + 84*y0^3*y1^2*y2
            + 46*y0^2*y1^3*y2 - 20*y0*y1^4*y2 + 16*y1^5*y2 + 53*y0^4*y2^2 + 56*y0^3*y1*y2^2
            - 32*y0^2*y1^2*y2^2 - 80*y0*y1^3*y2^2 - 92*y1^4*y2^2 - 12*y0^2*y1*y2^3
            - 168*y0*y1^2*y2^3 - 122*y1^3*y2^3 + 14*y0^2*y2^4 + 8*y0*y1*y2^4 - 112*y1^2*y2^4 + y2^6
        """
        return ((self._Lcoeff(i, 0))**2)*(self._Qcoeff(i, 1, 2))**2 + \
            ((self._Lcoeff(i, 1))**2)*(self._Qcoeff(i, 0, 2)**2)+ \
            ((self._Lcoeff(i, 2))**2)*(self._Qcoeff(i, 0, 1)**2)- \
            2*(self._Lcoeff(i, 0))*(self._Lcoeff(i, 1))*(self._Qcoeff(i, 0, 2))*(self._Qcoeff(i, 1, 2))\
            -2*(self._Lcoeff(i, 0))*(self._Lcoeff(i, 2))*(self._Qcoeff(i, 0, 1))*(self._Qcoeff(i, 1, 2))\
            -2*(self._Lcoeff(i, 1))*(self._Lcoeff(i, 2))*(self._Qcoeff(i, 0, 1))*(self._Qcoeff(i, 0, 2)) + \
             4*(self._Lcoeff(i, 0))*(self._Lcoeff(i, 1))*(self._Qcoeff(i, 0, 1))*(self._Qcoeff(i, 2, 2)) + \
             4*(self._Lcoeff(i, 0))*(self._Lcoeff(i, 2))*(self._Qcoeff(i, 0, 2))*(self._Qcoeff(i, 1, 1)) + \
             4*(self._Lcoeff(i, 1))*(self._Lcoeff(i, 2))*(self._Qcoeff(i, 1, 2))*(self._Qcoeff(i, 0, 0)) - \
             4*((self._Lcoeff(i, 0))**2)*(self._Qcoeff(i, 1, 1))*(self._Qcoeff(i, 2, 2)) - \
             4*((self._Lcoeff(i, 1))**2)*(self._Qcoeff(i, 0, 0))*(self._Qcoeff(i, 2, 2)) - \
             4*((self._Lcoeff(i, 2))**2)*(self._Qcoeff(i, 1, 1))*(self._Qcoeff(i, 0, 0))

    @cached_method
    def is_degenerate(self):
        r"""
        Function will return True if there is a fiber (over the algebraic closure of the
        base ring) of dimension greater than 0 and False otherwise.

        OUTPUT: boolean

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 + x2^2*y2^2 + \
            x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.is_degenerate()
            True

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - \
            4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.is_degenerate()
            False

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], GF(3))
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - \
            4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.is_degenerate()
            True
        """
        PP = self.ambient_space()
        K = FractionField(PP[0].base_ring())
        R = PP.coordinate_ring()
        PS = PP[0] #check for x fibers
        vars = list(PS.gens())
        R0 = PolynomialRing(K, 3, vars) #for dimension calculation to work,
            #must be done with Polynomial ring over a field
        #Degenerate is equivalent to a common zero, see Prop 1.4 in [CS1996]_
        I = R.ideal(self.Gpoly(1, 0), self.Gpoly(1, 1), self.Gpoly(1, 2), self.Hpoly(1, 0, 1),
                    self.Hpoly(1, 0, 2), self.Hpoly(1, 1, 2))
        phi = R.hom(vars + [0, 0, 0], R0)
        I = phi(I)
        if I.dimension() != 0:
            return True

        PS = PP[1] #check for y fibers
        vars = list(PS.gens())
        R0 = PolynomialRing(K,3,vars) #for dimension calculation to work,
        #must be done with Polynomial ring over a field
        #Degenerate is equivalent to a common zero, see Prop 1.4 in [CS1996]_
        I = R.ideal(self.Gpoly(0, 0), self.Gpoly(0, 1), self.Gpoly(0, 2), self.Hpoly(0, 0, 1),
                    self.Hpoly(0, 0, 2), self.Hpoly(0, 1, 2))
        phi = R.hom([0, 0, 0] + vars, R0)
        I = phi(I)
        return I.dimension() != 0

    def degenerate_fibers(self):
        r"""
        Return the (rational) degenerate fibers of the surface defined over
        the base ring, or the fraction field of the base ring if it is not a field.

        ALGORITHM:

        The criteria for degeneracy by the common vanishing of the polynomials
        ``self.Gpoly(1, 0)``, ``self.Gpoly(1, 1)``, ``self.Gpoly(1, 2)``,
        ``self.Hpoly(1, 0, 1)``,``self.Hpoly(1, 0, 2)``,
        ``self.Hpoly(1, 1, 2)`` (for the first component), is from Proposition 1.4
        in the following article: [CS1996]_.

        This function finds the common solution through elimination via Groebner bases
        by using the .variety() function on the three affine charts in each component.

        OUTPUT: The output is a list of lists where the elements of lists are
        points in the appropriate projective space.
        The first list is the points whose pullback by the projection to the
        first component (projective space) is dimension greater than 0.
        The second list is points in the second component

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 + x2^2*y2^2\
            + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.degenerate_fibers()
            [[], [(1 : 0 : 0)]]

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1\
            - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2\
            - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.degenerate_fibers()
            [[], []]

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: R = PP.coordinate_ring()
            sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
            sage: q = (y1*y0 + y2^2)*x0^2 + ((y0^2 - y2*y1)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0 \
            + (y2*y0 + y1^2)*x1^2 + (y0^2 + (-y1^2 + y2^2))*x2*x1
            sage: X = WehlerK3Surface([l,q])
            sage: X.degenerate_fibers()
            [[(-1 : 1 : 1), (0 : 0 : 1)], [(-1 : -1 : 1), (0 : 0 : 1)]]
        """
        PP = self.ambient_space()
        R = PP.coordinate_ring()
        PSX = PP[0]
        vars = list(PSX.gens())
        K = FractionField(PSX.base_ring())
        R0 = PolynomialRing(K, 3, vars)
        I = R.ideal(self.Gpoly(1, 0), self.Gpoly(1, 1), self.Gpoly(1, 2), self.Hpoly(1, 0,1 ),
                    self.Hpoly(1, 0, 2), self.Hpoly(1, 1, 2))
        phi = R.hom(vars + [0, 0, 0], R0)
        I = phi(I)
        xFibers = []
        #check affine charts
        for n in range(3):
            affvars = list(R0.gens())
            del affvars[n]
            R1 = PolynomialRing(K, 2, affvars, order='lex')
            mapvars = list(R1.gens())
            mapvars.insert(n,1)
            phi1 = R0.hom(mapvars, R1)
            J = phi1(I)
            if (J.dimension() == 0):
                Var = J.variety()
                if Var != [{}]:
                    for d in Var: #iterate through dictionaries
                        P = [] #new point
                        for z in mapvars: #assign coordinate values
                            if (z == 1):
                                P.append(1)
                            else:
                                P.append(d[z])
                        MP = PSX(P) #make projective point
                        if MP not in xFibers:
                            xFibers.append(MP)
        PSY = PP[1]
        vars = list(PSY.gens())
        K = FractionField(PSY.base_ring())
        R0 = PolynomialRing(K, 3, vars)
        I = R.ideal(self.Gpoly(0, 0), self.Gpoly(0, 1), self.Gpoly(0, 2), self.Hpoly(0, 0, 1),
                    self.Hpoly(0, 0, 2), self.Hpoly(0, 1, 2))
        phi = PP.coordinate_ring().hom([0, 0, 0] + vars, R0)
        I = phi(I)
        yFibers = []
        #check affine charts
        for n in range(3):
            affvars = list(R0.gens())
            del affvars[n]
            R1 = PolynomialRing(K, 2, affvars, order='lex')
            mapvars = list(R1.gens())
            mapvars.insert(n, 1)
            phi1 = R0.hom(mapvars,R1)
            J = phi1(I)
            if (J.dimension() == 0):
                Var = J.variety()
                if Var != [{}]:
                    for d in Var: #iterate through dictionaries
                        P = [] #new point
                        for z in mapvars: #assign coordinate values
                            if (z == 1):
                                P.append(1)
                            else:
                                P.append(d[z])
                        MP = PSY(P) #make projective point
                        if MP not in yFibers:
                            yFibers.append(MP)
        return [xFibers,yFibers]

    @cached_method
    def degenerate_primes(self,check=True):
        r"""
        Determine which primes `p` self has degenerate fibers over `GF(p)`.

        If check is False, then may return primes that do not have degenerate fibers.
        Raises an error if the surface is degenerate.
        Works only for ``ZZ`` or ``QQ``.

        INPUT: ``check`` -- (default: True) boolean, whether the primes are verified

        ALGORITHM:

        `p` is a prime of bad reduction if and only if the defining
        polynomials of self plus the G and H polynomials have a common
        zero. Or stated another way, `p` is a prime of bad reduction if
        and only if the radical of the ideal defined by the defining
        polynomials of self plus the G and H polynomials is not
        `(x_0,x_1,\ldots,x_N)`.  This happens if and only if some
        power of each `x_i` is not in the ideal defined by the
        defining polynomials of self (with G and H). This last condition
        is what is checked. The lcm of the coefficients of the monomials `x_i` in
        a groebner basis is computed. This may return extra primes.

        OUTPUT: List of primes.

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ, 6)
            sage: L =  y0*x0 + (y1*x1 + y2*x2)
            sage: Q = (2*y0^2 + y2*y0 + (2*y1^2 + y2^2))*x0^2 + ((y0^2 + y1*y0 + \
            (y1^2 + 2*y2*y1 + y2^2))*x1 + (2*y1^2 + y2*y1 + y2^2)*x2)*x0 + ((2*y0^2\
            + (y1 + 2*y2)*y0 + (2*y1^2 + y2*y1))*x1^2 + ((2*y1 + 2*y2)*y0 + (y1^2 + \
            y2*y1 + 2*y2^2))*x2*x1 + (2*y0^2 + y1*y0 + (2*y1^2 + y2^2))*x2^2)
            sage: X = WehlerK3Surface([L, Q])
            sage: X.degenerate_primes()
            [2, 3, 5, 11, 23, 47, 48747691, 111301831]
        """
        PP = self.ambient_space()
        if PP.base_ring() in _NumberFields or is_NumberFieldOrder(PP.base_ring()):
            if PP.base_ring() != ZZ and PP.base_ring() != QQ:
                raise NotImplementedError("must be ZZ or QQ")
        else:
            raise TypeError("must be over a number field")
        if self.is_degenerate():
            raise TypeError("surface is degenerate at all primes")
        RR = PP.coordinate_ring()

        #x-fibers
        PSX = PP[0]
        vars = list(PSX.gens())
        K = PSX.base_ring()
        R = PolynomialRing(K, 3, vars)
        I = RR.ideal(self.Gpoly(1, 0), self.Gpoly(1, 1), self.Gpoly(1, 2), self.Hpoly(1, 0, 1),
                     self.Hpoly(1, 0, 2), self.Hpoly(1, 1, 2))
        phi = PP.coordinate_ring().hom(vars + [0, 0, 0], R)
        I = phi(I)
        bad_primes = []

        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S = PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            I = S.ideal(I.gens())
        GB = I.groebner_basis()
        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT = GB[i].lt().degrees()
            power = 0
            for j in range(R.ngens()):
                if LT[j] != 0:
                    power += 1
            if power == 1:
                bad_primes = bad_primes+GB[i].lt().coefficients()[0].support()

        #y-fibers
        PSY = PP[1]
        vars = list(PSY.gens())
        K = PSY.base_ring()
        R = PolynomialRing(K, 3, vars)
        I = RR.ideal(self.Gpoly(0, 0), self.Gpoly(0, 1), self.Gpoly(0, 2), self.Hpoly(0, 0, 1),
                     self.Hpoly(0, 0, 2), self.Hpoly(0, 1, 2))
        phi = PP.coordinate_ring().hom([0, 0, 0] + vars, R)
        I = phi(I)
        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S = PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            I = S.ideal(I.gens())
        GB = I.groebner_basis()
        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT = GB[i].lt().degrees()
            power = 0
            for j in range(R.ngens()):
                if LT[j] != 0:
                    power += 1
            if power == 1:
                bad_primes = bad_primes+GB[i].lt().coefficients()[0].support()
        bad_primes = sorted(set(bad_primes))
        #check to return only the truly bad primes
        if check:
            for p in bad_primes:
                X = self.change_ring(GF(p))
                if not X.is_degenerate():
                    bad_primes.remove(p)
        return bad_primes

    def is_smooth(self):
        r"""
        Function will return the status of the smoothness of the surface.

        ALGORITHM:

        Checks to confirm that all of the 2x2 minors of the Jacobian generated from
        the biquadratic and bilinear forms have no common vanishing points.

        OUTPUT:  Boolean

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +\
             x2^2*y2^2 + x2^2*y1^2 + x1^2*y2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.is_smooth()
            False

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.is_smooth()
            True
        """
        vars = list(self.ambient_space().gens())
        M = jacobian([self.L, self.Q], vars)
        R = self.ambient_space().coordinate_ring()
        I = R.ideal(M.minors(2) + [self.L,self.Q])
        T = PolynomialRing(self.ambient_space().base_ring().fraction_field(), 4, 'h')
        #check the 9 affine charts for a singular point
        for l in xmrange([3, 3]):
            vars = list(T.gens())
            vars.insert(l[0], 1)
            vars.insert(3 + l[1], 1)
            phi = R.hom(vars, T)
            J = phi(I)
            if J.dimension() != -1:
                return False
        return True

    def sigmaX(self, P, **kwds):
        r"""
        Function returns the involution on the  Wehler K3 surface induced by the double covers.

        In particular, it fixes the projection to the first coordinate and swaps the
        two points in the fiber, i.e. `(x, y) \to (x, y')`.
        Note that in the degenerate case, while we can split fiber into pairs of points,
        it is not always possibleto distinguish them, using this algorithm.

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute `\sigma_x`, `\sigma_y`, `\phi`,
        and `\psi`" in [CS1996FH2015.
        For the degenerate case refer to [FH2015]_.

        INPUT:

        - ``P`` - a point in `\mathbb{P}^2 \times \mathbb{P}^2`

        kwds:

        - ``check`` - (default: ``True``) boolean checks to see if point is on the surface

        - ``normalize`` -- (default: ``True``) boolean normalizes the point

        OUTPUT: A point on the K3 surface

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 +\
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -\
            4*x1*x2*y1^2 + 5*x0*x2*y0*y2 -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 +\
            4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP(0, 0, 1, 1, 0, 0)
            sage: X.sigmaX(T)
            (0 : 0 : 1 , 0 : 1 : 0)

        degenerate examples::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
            sage: q = (y1*y0)*x0^2 + ((y0^2)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0\
            + (y2*y0 + y1^2)*x1^2 + (y0^2 + (-y1^2 + y2^2))*x2*x1
            sage: X = WehlerK3Surface([l, q])
            sage: X.sigmaX(X([1, 0, 0, 0, 1, -2]))
            (1 : 0 : 0 , 0 : 1/2 : 1)
            sage: X.sigmaX(X([1, 0, 0, 0, 0, 1]))
            (1 : 0 : 0 , 0 : 0 : 1)
            sage: X.sigmaX(X([-1, 1, 1, -1, -1, 1]))
            (-1 : 1 : 1 , 2 : 2 : 1)
            sage: X.sigmaX(X([0, 0, 1, 1, 1, 0]))
            (0 : 0 : 1 , 1 : 1 : 0)
            sage: X.sigmaX(X([0, 0, 1, 1, 1, 1]))
            (0 : 0 : 1 , -1 : -1 : 1)

        Case where we cannot distinguish the two points::

            sage: PP.<y0,y1,y2,x0,x1,x2>=ProductProjectiveSpaces([2, 2], GF(3))
            sage: l = x0*y0 + x1*y1 + x2*y2
            sage: q=-3*x0^2*y0^2 + 4*x0*x1*y0^2 - 3*x0*x2*y0^2 - 5*x0^2*y0*y1 - \
            190*x0*x1*y0*y1- 5*x1^2*y0*y1 + 5*x0*x2*y0*y1 + 14*x1*x2*y0*y1 + \
            5*x2^2*y0*y1 - x0^2*y1^2 - 6*x0*x1*y1^2- 2*x1^2*y1^2 + 2*x0*x2*y1^2 - \
            4*x2^2*y1^2 + 4*x0^2*y0*y2 - x1^2*y0*y2 + 3*x0*x2*y0*y2+ 6*x1*x2*y0*y2 - \
            6*x0^2*y1*y2 - 4*x0*x1*y1*y2 - x1^2*y1*y2 + 51*x0*x2*y1*y2 - 7*x1*x2*y1*y2 - \
            9*x2^2*y1*y2 - x0^2*y2^2 - 4*x0*x1*y2^2 + 4*x1^2*y2^2 - x0*x2*y2^2 + 13*x1*x2*y2^2 - x2^2*y2^2
            sage: X = WehlerK3Surface([l, q])
            sage: P = X([1, 0, 0, 0, 1, 1])
            sage: X.sigmaX(X.sigmaX(P))
            Traceback (most recent call last):
            ...
            ValueError: cannot distinguish points in the degenerate fiber
        """
        check = kwds.get("check", True)
        normalize = kwds.get("normalize", True)

        if check:
            if self != P.codomain():
                try:
                    P = self(list(P))
                except (TypeError, NotImplementedError, AttributeError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self))
        pt = list(P[0]) + [0, 0, 0]
        if P[1][0] != 0:
            [a,b,c] = [P[1][0]*self.Gpoly(1, 0)(*pt),
                       -1*P[1][0]*self.Hpoly(1, 0, 1)(*pt) - P[1][1]*self.Gpoly(1, 0)(*pt),
                       -P[1][0]*self.Hpoly(1, 0, 2)(*pt) - P[1][2]*self.Gpoly(1, 0)(*pt)]
        elif P[1][1] != 0:
            [a,b,c] = [-1*P[1][1]*self.Hpoly(1, 0, 1)(*pt)-P[1][0]*self.Gpoly(1, 1)(*pt),
                        P[1][1]*self.Gpoly(1, 1)(*pt),
                       -P[1][1]*self.Hpoly(1, 1, 2)(*pt)-P[1][2]*self.Gpoly(1, 1)(*pt)]
        else:
            [a,b,c] = [-1*P[1][2]*self.Hpoly(1, 0, 2)(*pt) - P[1][0]*self.Gpoly(1, 2)(*pt),
                       -P[1][2]*self.Hpoly(1, 1, 2)(*pt) - P[1][1]*self.Gpoly(1, 2)(*pt),
                       P[1][2]*self.Gpoly(1, 2)(*pt)]
        Point = [P[0][0], P[0][1], P[0][2], a, b, c]

        if [a,b,c] != [0,0,0]:
            if normalize:
                Point = self.point(Point,False)
                Point.normalize_coordinates()
                return Point
            return self.point(Point,False)
        #Start of the degenerate case
        R = self.ambient_space().coordinate_ring()
        BR = self.ambient_space().base_ring()
        S = PolynomialRing(BR, 6, 's0, s1, w1, z0, z1, z2')
        s0,s1,w1,z0,z1,z2 = S.gens()
        #Define the blow-up map with (s0,s1) the new `\mathbb{P}^1` coordinates
        #so that the points on the fiber come in pairs on the lines defined by `(s0,s1)`
        #this allows us to extend the involution to degenerate fibers
        if P[0][0]!= 0:
            t1 = BR(P[0][1]/P[0][0])
            t = w1 - t1
            phi = R.hom([s0, s0*w1, s1*t + s0*P[0][2]/P[0][0], z0, z1, z2], S)
        elif P[0][1]!= 0:
            t1 = BR(P[0][0]/P[0][1])
            t = w1 - t1
            phi = R.hom([s0*w1, s0, s1*t + s0*P[0][2]/P[0][1], z0, z1, z2], S)
        else:
            t1 = BR(P[0][1]/P[0][2])
            t = w1 - t1
            phi = R.hom([s1*(t) + s0*P[0][0]/P[0][2], s0*w1, s0, z0, z1, z2], S)

        # Blow-up the fiber
        T = [phi(self.L),phi(self.Q),
             phi(self.Gpoly(1, 0)),
             phi(self.Gpoly(1, 1)),
             phi(self.Gpoly(1, 2)),
            -phi(self.Hpoly(1, 0, 1)),
            -phi(self.Hpoly(1, 0, 2)),
            -phi(self.Hpoly(1, 1, 2))]
        maxexp = []

        #Find highest exponent that we can divide out by to get a non zero answer
        for i in range(2,len(T)):
            e = 0
            while (T[i]/t**e).subs({w1:t1}) == 0:
                e += 1
            maxexp.append(e)

        e = min(maxexp)

        #Fix L and Q
        for i in range(0,2):
            while T[i].subs({w1:t1}) == 0:
                T[i] = T[i]/t
            T[i] = T[i].subs({w1:t1})

        #Fix G and H polys
        for i in range(2,len(T)):
            T[i] = T[i]/t**e
            T[i] = T[i].subs({w1:t1})

        #Defines the ideal whose solution gives `(s0, s1)` and the two points
        #on the fiber
        RR = PolynomialRing(BR, 5,'s0, s1, z0, z1, z2',order='lex')
        s0, s1, z0, z1, z2 = RR.gens()
        I = RR.ideal([RR(T[0]),
                      RR(T[1]),
                      RR(T[2]) - P[1][0]*z0, RR(T[3]) - P[1][1]*z1, RR(T[4])-P[1][2]*z2,
                      RR(T[5]) - (P[1][0]*z1 + P[1][1]*z0),
                      RR(T[6]) - (P[1][0]*z2 + P[1][2]*z0),
                      RR(T[7]) - (P[1][1]*z2 + P[1][2]*z1)])

        #Find the points
        SS = PolynomialRing(BR, 4,'s, z0, z1, z2', order='lex')
        s, z0, z1, z2 = SS.gens()
        phi = RR.hom([s, 1, z0, z1, z2], SS)
        J = phi(I)
        if J.dimension() > 0:
            raise ValueError("cannot distinguish points in the degenerate fiber")
        V = J.variety()
        #Our blow-up point has more than one line passing through it, thus we cannot find
        #the corresponding point on the surface
        if len(V) > 2:
            raise ValueError("cannot distinguish points in the degenerate fiber")
        #We always expect to have the trivial solution (0, 0, 0)
        if len(V) == 2:
            for D in V:
                if D[s] != 0:
                    [a,b,c] = [D[z0], D[z1], D[z2]]
        else:
            newT = [phi(tee) for tee in T]
            for i in range(2):
                while newT[i] != 0 and s.divides(newT[i]):
                    newT[i] = SS(newT[i]/s)
            maxexp = []

            for i in range(2, len(T)):
                e = 0
                if newT[i] != 0:
                    while (newT[i]/s**e).subs({s:0})==0:
                        e += 1
                    maxexp.append(e)
            e = min(maxexp)

            #Cancel the powers of s
            for i in range(2,len(T)):
                newT[i] = newT[i]/s**e
            #Create the new ideal
            II = SS.ideal([SS(newT[0]),
                           SS(newT[1]),
                           SS(newT[2]) - P[1][0]*z0,
                           SS(newT[3]) - P[1][1]*z1,
                           SS(newT[4]) - P[1][2]*z2,
                           SS(newT[5]) - (P[1][0]*z1 + P[1][1]*z0),
                           SS(newT[6]) - (P[1][0]*z2 + P[1][2]*z0),
                           SS(newT[7]) - (P[1][1]*z2 + P[1][2]*z1)])

            #Find the points
            SSS = PolynomialRing(BR, 3, 'z0, z1, z2', order='lex')
            z0,z1,z2 = SSS.gens()
            phi = SS.hom([0, z0, z1, z2], SSS)
            J2 = phi(II)
            if J2.dimension() > 0:
                raise ValueError("cannot distinguish points in the degenerate fiber")
            V = J2.variety()
            if len(V) > 1:
                raise ValueError("cannot distinguish points in the degenerate fiber")

            if len(V) == 1:
                [a, b, c] = [V[0][z0], V[0][z1], V[0][z2]]

            if len(V) == 0 or [a,b,c] == [0, 0, 0]:
                SS = PolynomialRing(BR, 3, 'z0, z1, z2', order='lex')
                z0,z1,z2 = SS.gens()
                phi = RR.hom([1, 0, z0, z1, z2], SS)
                J = phi(I)
                if J.dimension() > 0:
                    raise ValueError( "cannot distinguish points in the degenerate fiber")
                V = phi(I).variety()
                if len(V) > 1:
                    raise ValueError( "cannot distinguish points in the degenerate fiber")
                [a,b,c] = [V[0][z0], V[0][z1], V[0][z2]]

        Point = [P[0][0], P[0][1], P[0][2], a, b, c]
        if normalize:
            Point = self.point(Point, False)
            Point.normalize_coordinates()
            return Point
        return self.point(Point, False)

    def sigmaY(self,P, **kwds):
        r"""
        Function returns the involution on the Wehler K3 surfaces induced by the double covers.

        In particular,it fixes the projection to the second coordinate and swaps
        the two points in the fiber, i.e. `(x,y) \to (x',y)`.
        Note that in the degenerate case, while we can split the fiber into two points,
        it is not always possibleto distinguish them, using this algorithm.

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute `\sigma_x`, `\sigma_y`, `\phi`,
        and `\psi`" in [CS1996]_.
        For the degenerate case refer to [FH2015]_.

        INPUT:

        - ``P`` - a point in `\mathbb{P}^2 \times \mathbb{P}^2`

        kwds:

        - ``check`` - (default: ``True``) boolean checks to see if point is on the surface

        - ``normalize`` -- (default: ``True``) boolean normalizes the point

        OUTPUT: A point on the K3 surface

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP(0, 0, 1, 1, 0, 0)
            sage: X.sigmaY(T)
            (0 : 0 : 1 , 1 : 0 : 0)

        degenerate examples::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
            sage: q = (y1*y0)*x0^2 + ((y0^2)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0 +\
             (y2*y0 + y1^2)*x1^2 + (y0^2 + (-y1^2 + y2^2))*x2*x1
            sage: X = WehlerK3Surface([l, q])
            sage: X.sigmaY(X([1, -1, 0 ,-1, -1, 1]))
            (1/10 : -1/10 : 1 , -1 : -1 : 1)
            sage: X.sigmaY(X([0, 0, 1, -1, -1, 1]))
            (-4 : 4 : 1 , -1 : -1 : 1)
            sage: X.sigmaY(X([1, 2, 0, 0, 0, 1]))
            (-3 : -3 : 1 , 0 : 0 : 1)
            sage: X.sigmaY(X([1, 1, 1, 0, 0, 1]))
            (1 : 0 : 0 , 0 : 0 : 1)

        Case where we cannot distinguish the two points::

            sage: PP.<x0,x1,x2,y0,y1,y2>=ProductProjectiveSpaces([2, 2], GF(3))
            sage: l = x0*y0 + x1*y1 + x2*y2
            sage: q=-3*x0^2*y0^2 + 4*x0*x1*y0^2 - 3*x0*x2*y0^2 - 5*x0^2*y0*y1 - 190*x0*x1*y0*y1 \
            - 5*x1^2*y0*y1 + 5*x0*x2*y0*y1 + 14*x1*x2*y0*y1 + 5*x2^2*y0*y1 - x0^2*y1^2 - 6*x0*x1*y1^2 \
            - 2*x1^2*y1^2 + 2*x0*x2*y1^2 - 4*x2^2*y1^2 + 4*x0^2*y0*y2 - x1^2*y0*y2 + 3*x0*x2*y0*y2 \
            + 6*x1*x2*y0*y2 - 6*x0^2*y1*y2 - 4*x0*x1*y1*y2 - x1^2*y1*y2 + 51*x0*x2*y1*y2 - 7*x1*x2*y1*y2 \
            - 9*x2^2*y1*y2 - x0^2*y2^2 - 4*x0*x1*y2^2 + 4*x1^2*y2^2 - x0*x2*y2^2 + 13*x1*x2*y2^2 - x2^2*y2^2
            sage: X = WehlerK3Surface([l ,q])
            sage: P = X([0, 1, 1, 1, 0, 0])
            sage: X.sigmaY(X.sigmaY(P))
            Traceback (most recent call last):
            ...
            ValueError: cannot distinguish points in the degenerate fiber
        """
        check = kwds.get("check", True)
        normalize = kwds.get("normalize", True)

        if check:
            if self != P.codomain():
                try:
                    P = self(list(P))
                except (TypeError, NotImplementedError, AttributeError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self))
        pt = [0, 0, 0] + list(P[1])
        if P[0][0] != 0:
            [a, b, c] = [P[0][0]*self.Gpoly(0, 0)(*pt),
                      -1*P[0][0]*self.Hpoly(0, 0, 1)(*pt) - P[0][1]*self.Gpoly(0, 0)(*pt),
                      -P[0][0]*self.Hpoly(0, 0, 2)(*pt) - P[0][2]*self.Gpoly(0, 0)(*pt)]
        elif P[0][1] != 0:
            [a, b, c] = [-1*P[0][1]*self.Hpoly(0, 0, 1)(*pt) - P[0][0]*self.Gpoly(0, 1)(*pt),
                       P[0][1]*self.Gpoly(0, 1)(*pt),
                       -P[0][1]*self.Hpoly(0, 1, 2)(*pt) - P[0][2]*self.Gpoly(0, 1)(*pt)]
        else:
            [a, b, c] = [-1*P[0][2]*self.Hpoly(0, 0, 2)(*pt) - P[0][0]*self.Gpoly(0, 2)(*pt),
                       - P[0][2]*self.Hpoly(0, 1, 2)(*pt) - P[0][1]*self.Gpoly(0, 2)(*pt),
                       P[0][2]*self.Gpoly(0, 2)(*pt)]
        Point = [a, b, c, P[1][0], P[1][1], P[1][2]]
        if [a, b, c] != [0, 0, 0]:
            if normalize:
                Point = self.point(Point, False)
                Point.normalize_coordinates()
                return Point
            return self.point(Point, False)

        #Start of the degenerate case
        R = self.ambient_space().coordinate_ring()
        BR = self.ambient_space().base_ring()
        S = PolynomialRing(BR, 6, 'z0, z1, z2, s0, s1, w1')
        z0, z1, z2, s0, s1, w1 = S.gens()
        #Define the blow-up map with (s0,s1) the new `\mathbb{P}^1` coordinates
        #so that the points on the fiber come in pairs on the lines defined by `(s0,s1)`
        #this allows us to extend the involution to degenerate fibers
        if P[1][0] != 0:
            t1 = BR(P[1][1]/P[1][0])
            t = w1 - t1
            phi = R.hom([z0, z1, z2, s0, s0*w1, s1*t + s0*P[1][2]/P[1][0]], S)
        elif P[1][1]!= 0:
            t1 = BR(P[1][0]/P[1][1])
            t = w1 - t1
            phi = R.hom([z0, z1, z2, s0*w1, s0, s1*t + s0*P[1][2]/P[1][1]], S)
        else:
            t1 = BR(P[1][1]/P[1][2])
            t = w1 - t1
            phi = R.hom([z0, z1, z2, s1*(t) + s0*P[1][0]/P[1][2], s0*w1, s0], S)

        #Blow-up the fiber
        T = [phi(self.L),
         phi(self.Q),
         phi(self.Gpoly(0, 0)),
         phi(self.Gpoly(0, 1)),
         phi(self.Gpoly(0, 2)),
        -phi(self.Hpoly(0, 0, 1)),
        -phi(self.Hpoly(0, 0, 2)),
        -phi(self.Hpoly(0, 1, 2))]
        maxexp = []

        #Find highest exponent that we can divide out by to get a non zero answer
        for i in range(2, len(T)):
            e = 0
            while (T[i]/t**e).subs({w1:t1}) == 0:
                e += 1
            maxexp.append(e)

        e = min(maxexp)

        for i in range(0, 2):
            while T[i].subs({w1:t1}) == 0:
                T[i] = T[i]/t
            T[i] = T[i].subs({w1:t1})
        for i in range(2, len(T)):
            T[i] = T[i]/t**e
            T[i] = T[i].subs({w1:t1})

        # Defines the ideal whose solution gives `(s0,s1)` and the two points
        # on the fiber
        RR = PolynomialRing(BR, 5, 's0, s1, z0, z1, z2', order='lex')
        s0, s1, z0, z1, z2 = RR.gens()
        I = RR.ideal([RR(T[0]),
                      RR(T[1]),
                      RR(T[2]) - P[0][0]*z0,
                      RR(T[3]) - P[0][1]*z1,
                      RR(T[4]) - P[0][2]*z2,
                      RR(T[5]) - (P[0][0]*z1 + P[0][1]*z0),
                      RR(T[6]) - (P[0][0]*z2 + P[0][2]*z0),
                      RR(T[7]) - (P[0][1]*z2 + P[0][2]*z1)])
        #Find the points
        SS = PolynomialRing(BR, 4, 's, z0, z1, z2', order='lex')
        s, z0, z1, z2 = SS.gens()
        phi = RR.hom([s, 1, z0, z1, z2], SS)
        J = phi(I)
        if J.dimension() > 0:
            raise ValueError("cannot distinguish points in the degenerate fiber")
        V = J.variety()

        #Our blow-up point has more than one line passing through it, thus we cannot find
        #the corresponding point on the surface
        if len(V) > 2:
            raise ValueError("cannot distinguish points in the degenerate fiber")
        # We always expect to have the trivial solution (0, 0, 0)
        if len(V) == 2:
            for D in V:
                if D[s] != 0:
                    a, b, c = [D[z0], D[z1], D[z2]]
        else:
            newT = [phi(tee) for tee in T]
            for i in range(2):
                while newT[i] != 0 and s.divides(newT[i]):
                    newT[i] = SS(newT[i]/s)
            maxexp = []
            for i in range(2, len(T)):
                e = 0
                if newT[i] != 0:
                    while (newT[i]/s**e).subs({s:0}) == 0:
                        e += 1
                    maxexp.append(e)
            e = min(maxexp)
            #Cancel out the powers of s
            for i in range(2,len(T)):
                newT[i] = newT[i]/s**e
            #Create the new ideal
            II = SS.ideal([SS(newT[0]),
                           SS(newT[1]),
                           SS(newT[2]) - P[0][0]*z0,
                           SS(newT[3]) - P[0][1]*z1,
                           SS(newT[4]) - P[0][2]*z2,
                           SS(newT[5]) - (P[0][0]*z1 + P[0][1]*z0),
                           SS(newT[6]) - (P[0][0]*z2 + P[0][2]*z0),
                           SS(newT[7]) - (P[0][1]*z2 + P[0][2]*z1)])
            # Find the points
            SSS = PolynomialRing(BR, 3, 'z0, z1, z2', order='lex')
            z0, z1, z2 = SSS.gens()
            phi = SS.hom([0, z0, z1, z2], SSS)
            J2 = phi(II)
            if J2.dimension() > 0:
                raise ValueError("cannot distinguish points in the degenerate fiber")
            V = J2.variety()

            if len(V) > 1:
                raise ValueError("cannot distinguish points in the degenerate fiber")
            if len(V) == 1:
                [a, b, c] = [V[0][z0], V[0][z1], V[0][z2]]
            if len(V) == 0 or [a,b,c] == [0, 0, 0]:
                SS = PolynomialRing(BR, 3, 'z0, z1, z2', order='lex')
                z0,z1,z2 = SS.gens()
                phi = RR.hom([1, 0, z0, z1, z2], SS)
                J = phi(I)
                if J.dimension() > 0:
                    raise ValueError("cannot distinguish points in the degenerate fiber")
                V = phi(I).variety()
                if len(V) > 1:
                    raise ValueError("cannot distinguish points in the degenerate fiber")
                [a,b,c] = [V[0][z0], V[0][z1], V[0][z2]]

        Point = [a, b, c, P[1][0], P[1][1], P[1][2]]
        if normalize:
            Point = self.point(Point, False)
            Point.normalize_coordinates()
            return Point
        return self.point(Point, False)

    def phi(self, a, **kwds):
        r"""
        Evaluates the function `\phi = \sigma_y \circ \sigma_x`.

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute `\sigma_x`, `\sigma_y`,
        `\phi`, and `\psi`" in [CS1996]_.

        For the degenerate case refer to [FH2015]_.

        INPUT:

        - ``a`` - Point in `\mathbb{P}^2 \times \mathbb{P}^2`

        kwds:

        - ``check`` - (default: ``True``) boolean checks to see if point is on the surface

        - ``normalize`` -- (default: ``True``) boolean normalizes the point

        OUTPUT: A point on this surface

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP([0, 0, 1, 1 ,0, 0])
            sage: X.phi(T)
            (-1 : 0 : 1 , 0 : 1 : 0)
       """
        A = self.sigmaX(a, **kwds)
        kwds.update({"check":False})
        return self.sigmaY(A, **kwds)

    def psi(self,a, **kwds):
        r"""
        Evaluates the function `\psi = \sigma_x \circ \sigma_y`.

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute `\sigma_x`, `\sigma_y`,
        `\phi`, and `\psi`" in [CS1996]_.

        For the degenerate case refer to [FH2015]_.

        INPUT:

        - ``a`` - Point in `\mathbb{P}^2 \times \mathbb{P}^2`

        kwds:

        - ``check`` - (default: ``True``) boolean checks to see if point is on the surface

        - ``normalize`` -- (default: ``True``) boolean normalizes the point

        OUTPUT: A point on this surface

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP([0, 0, 1, 1, 0, 0])
            sage: X.psi(T)
            (0 : 0 : 1 , 0 : 1 : 0)
        """
        A = self.sigmaY(a, **kwds)
        kwds.update({"check":False})
        return self.sigmaX(A, **kwds)

    def lambda_plus(self, P, v, N, m, n, prec=100):
        r"""
        Evaluates the  local canonical height plus function of Call-Silverman at
        the place ``v`` for ``P`` with ``N`` terms of the series.

        Use ``v = 0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        Sum over local heights using convergent series, for more details,
        see section 4 of [CS1996]_.

        INPUT:

        - ``P`` -- a surface point

        - ``N`` -- positive integer. number of terms of the series to use

        - ``v`` -- non-negative integer. a place, use v = 0 for the Archimedean place

        - ``m,n`` -- positive integers, We compute the local height for the divisor `E_{mn}^{+}`.
                     These must be indices of non-zero coordinates of the point ``P``.

        - ``prec`` -- (default: 100) float point or p-adic precision

        OUTPUT: A real number

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1\
            - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2\
            - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: P = X([0, 0, 1, 1, 0, 0])
            sage: X.lambda_plus(P, 0, 10, 2, 0)
            0.89230705169161608922595928129
        """
        if not (v == 0 or v.is_prime()):
            raise ValueError("invalid valuation (= %s) entered"%v)
        R = RealField(prec)
        if v == 0:
            K = R
        else:
            K = Qp(v, prec)
        PK = P.change_ring(K)
        W = self.change_ring(K)
        Rx = W.ambient_space().coordinate_ring().hom(
             list(W.ambient_space()[0].coordinate_ring().gens()) + [0, 0, 0],
             W.ambient_space()[0].coordinate_ring())
        Ry = W.ambient_space().coordinate_ring().hom(
             [0, 0, 0] + list(W.ambient_space()[1].coordinate_ring().gens()),
             W.ambient_space()[1].coordinate_ring())
        beta = R(2 + sqrt(3))
        L = [x.abs() for x in list(PK[0])]
        i = L.index(max(L))
        L = [y.abs() for y in list(PK[1])]
        j = L.index(max(L))

        #Compute the local height wrt the divisor E_{mn}^{+}
        local_height = beta*R((PK[0][i]/PK[0][m]).abs()).log() - R((PK[1][j]/PK[1][n]).abs()).log()

        for e in range(N):
            #Take next iterate
            Q = W.phi(PK, check=False)
            L = [x.abs() for x in list(Q[0])]
            k = L.index(max(L))
            L = [y.abs() for y in list(Q[1])]
            l = L.index(max(L))
            newP = copy(PK)
            #normalize PK
            newP.scale_by([1/PK[0][i],1])

            #Find B and A, helper values for the local height
            if PK[1][j].abs() <= PK[1][l].abs():
                B = Rx(W.Gpoly(1, l))(tuple(newP[0]))*PK[1][j]/PK[1][l]
            else:
                B = -Rx(W.Gpoly(1, j))(tuple(newP[0]))*PK[1][l]/PK[1][j]
                B = B - Rx(W.Hpoly(1, j, l))(tuple(newP[0]))

            #Normalize Q
            newQ = copy(Q)
            newQ.scale_by([1,1/Q[1][l]])

            if PK[0][i].abs() <= PK[0][k].abs():
                A = Ry(W.Gpoly(0, k))(tuple(newQ[1]))*PK[0][i]/PK[0][k]
            else:
                A = -Ry(W.Gpoly(0, i))(tuple(newQ[1]))*PK[0][k]/PK[0][i]
                A = A - Ry(W.Hpoly(0, i, k))(tuple(newQ[1]))
            #Compute the new local height
            local_height += beta**(-2*R(e)-1)*R(A.abs()).log() + beta**(-2*R(e))*R(B.abs()).log()

            i = k
            j = l
            newQ.scale_by([1/Q[0][k], 1])
            PK = newQ
        return local_height

    def lambda_minus(self, P, v, N, m, n, prec=100):
        r"""
        Evaluates the local canonical height minus function of Call-Silverman
        at the place ``v`` for ``P`` with ``N`` terms of the series.

        Use ``v = 0`` for the Archimedean place. Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        Sum over local heights using convergent series, for more details,
        see section 4 of [CS1996]_.

        INPUT:

        - ``P`` -- a projective point

        - ``N`` -- positive integer. number of terms of the series to use

        - ``v`` -- non-negative integer. a place, use v = 0 for the Archimedean place

        - ``m,n`` -- positive integers, We compute the local height for the divisor `E_{mn}^{+}`.
                    These must be indices of non-zero coordinates of the point ``P``.

        - ``prec`` -- (default: 100) float point or p-adic precision

        OUTPUT: A real number

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 \
            - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2\
            - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: P = X([0, 0, 1, 1, 0, 0])
            sage: X.lambda_minus(P, 2, 20, 2, 0, 200)
            -0.18573351672047135037172805779671791488351056677474271893705
        """
        R = RealField(prec)
        if v == 0:
            K = R
        else:
            K = Qp(v, prec)
        PK = P.change_ring(K)
        W = self.change_ring(K)
        Rx = W.ambient_space().coordinate_ring().hom(list(W.ambient_space()[0].coordinate_ring().gens())
            + [0, 0, 0], W.ambient_space()[0].coordinate_ring())
        Ry = W.ambient_space().coordinate_ring().hom([0, 0, 0] +
                list(W.ambient_space()[1].coordinate_ring().gens()),
                W.ambient_space()[1].coordinate_ring())
        beta = R(2 + sqrt(3))
        L = [x.abs() for x in list(PK[0])]
        j = L.index(max(L))
        L = [y.abs() for y in list(PK[1])]
        i = L.index(max(L))

        ##Compute the local height wrt the divisor E_{mn}^{-}
        local_height = beta*R((PK[1][i]/PK[1][n]).abs()).log() - R((PK[0][j]/PK[0][m]).abs()).log()
        for e in range(N):
            #Take the next iterate
            Q = W.psi(PK, check=False)
            L = [x.abs() for x in list(Q[0])]
            l = L.index(max(L))
            L = [y.abs() for y in list(Q[1])]
            k = L.index(max(L))
            #Normalize the point
            newP = copy(PK)
            newP.scale_by([1, 1/PK[1][i]])
            #Find A and B, helper functions for computing local height
            if PK[0][j].abs() <= PK[0][l].abs():
                B = Ry(W.Gpoly(0, l))(tuple(newP[1]))*PK[0][j]/PK[0][l]
            else:
                B = -Ry(W.Gpoly(0, j))(tuple(newP[1]))*PK[0][l]/PK[0][j]
                B = B - Ry(W.Hpoly(0, j, l))(tuple(newP[1]))

            #Normalize Q
            newQ = copy(Q)
            newQ.scale_by([1/Q[0][l], 1])

            if PK[1][i].abs()<= PK[1][k].abs():
                A = Rx(W.Gpoly(1, k))(tuple(newQ[0]))*PK[1][i]/PK[1][k]
            else:
                A = -Rx(W.Gpoly(1, i))(tuple(newQ[0]))*PK[1][k]/PK[1][i]
                A = A-Rx(W.Hpoly(1, i, k))(tuple(newQ[0]))

            #Compute the local height
            local_height += beta**(-2*R(e)-1)*R(A.abs()).log() + beta**(-2*R(e))*R(B.abs()).log()
            i = k
            j = l
            newQ.scale_by([1, 1/Q[1][k]])
            PK = newQ
        return local_height

    def canonical_height_plus(self, P, N, badprimes=None, prec=100):
        r"""
        Evaluates the canonical height plus function of Call-Silverman
        for ``P`` with ``N`` terms of the series of the local heights.

        Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        Sum over the lambda plus heights (local heights) in a convergent series,
        for more detail see section 7 of [CS1996]_.

        INPUT:

        - ``P`` -- a surface point

        - ``N`` -- positive integer. Number of terms of the series to use

        - ``badprimes`` -- (optional) list of integer primes (where the surface is degenerate)

        - ``prec`` -- (default: 100) float point or p-adic precision

        OUTPUT: A real number

        EXAMPLES::

            sage: set_verbose(None)
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ, 6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + \
            (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 \
            + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L, Q])
            sage: P = X([1, 0, -1, 1, -1, 0]) #order 16
            sage: X.canonical_height_plus(P, 5)  # long time
            0.00000000000000000000000000000

        Call-Silverman Example::

            sage: set_verbose(None)
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: P = X([0, 1, 0, 0, 0, 1])
            sage: X.canonical_height_plus(P, 4) # long time
            0.14752753298983071394400412161
        """
        if badprimes is None:
            badprimes = self.degenerate_primes()
        m = 2
        while P[0][m] == 0:
            m = m - 1
        n = 2
        while P[1][n] == 0:
            n = n-1
        h = self.lambda_plus(P, 0, N, m, n, prec)
        for p in badprimes:
            h += self.lambda_plus(P, p, N, m, n, prec)
        return h

    def canonical_height_minus(self, P, N, badprimes=None, prec=100):
        r"""
        Evaluates the canonical height minus function of Call-Silverman
        for ``P`` with ``N`` terms of the series of the local heights.

        Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        Sum over the lambda minus heights (local heights) in a convergent series,
        for more detail see section 7 of [CS1996]_.

        INPUT:

        - ``P`` -- a surface point

        - ``N`` -- positive integer (number of terms of the series to use)

        - ``badprimes`` -- (optional) list of integer primes (where the surface is degenerate)

        - ``prec`` -- (default: 100) float point or p-adic precision

        OUTPUT: A real number

        EXAMPLES::

            sage: set_verbose(None)
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ, 6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1\
             + (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1\
              + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L, Q])
            sage: P = X([1, 0, -1, 1, -1, 0]) #order 16
            sage: X.canonical_height_minus(P, 5)  # long time
            0.00000000000000000000000000000

        Call-Silverman example::

            sage: set_verbose(None)
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 +\
             3*x0*x1*y0*y1 - 2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - \
             4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + \
             x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: P = X([0, 1, 0, 0, 0, 1])
            sage: X.canonical_height_minus(P, 4) # long time
            0.55073705369676788175590206734
        """
        if badprimes is None:
            badprimes = self.degenerate_primes()
        m = 2
        while P[0][m] == 0:
            m = m - 1
        n = 2
        while P[1][n] == 0:
            n = n-1
        h = self.lambda_minus(P, 0, N, m, n, prec)
        for p in badprimes:
            h += self.lambda_minus(P, p, N, m, n, prec)
        return h

    def canonical_height(self, P, N, badprimes=None, prec=100):
        r"""
        Evaluates the canonical height for ``P`` with ``N`` terms of the series of the local
        heights.

        ALGORITHM:

        The sum of the canonical height minus and canonical height plus,
        for more info see section 4 of [CS1996]_.

        INPUT:

        - ``P`` -- a surface point

        - ``N`` -- positive integer (number of terms of the series to use)

        - ``badprimes`` -- (optional) list of integer primes (where the surface is degenerate)

        - ``prec`` -- (default: 100) float point or p-adic precision

        OUTPUT: A real number

        EXAMPLES::

            sage: set_verbose(None)
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ, 6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + \
            (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 \
            + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L, Q])
            sage: P = X([1, 0, -1, 1,- 1, 0]) #order 16
            sage: X.canonical_height(P, 5)  # long time
            0.00000000000000000000000000000

        Call-Silverman example::

            sage: set_verbose(None)
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 \
            -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: P = X(0, 1, 0, 0, 0, 1)
            sage: X.canonical_height(P, 4)
            0.69826458668659859569990618895
        """
        if badprimes is None:
            badprimes = self.degenerate_primes()
        return (self.canonical_height_plus(P, N, badprimes, prec) +
                self.canonical_height_minus(P, N, badprimes, prec))

    def fiber(self, p, component):
        r"""
        Return the fibers [y (component = 1) or x (Component = 0)] of
        a point on a K3 Surface.

        This will work for nondegenerate fibers only.

        For algorithm, see [Hutz2007]_.

        INPUT:

        -``p`` - a point in `\mathbb{P}^2`

        OUTPUT: The corresponding fiber (as a list)

        EXAMPLES::

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ, 6)
            sage: Y = x0*y0 + x1*y1 - x2*y2
            sage: Z = y0^2*x0*x1 + y0^2*x2^2 - y0*y1*x1*x2 + y1^2*x2*x1 + y2^2*x2^2 +\
            y2^2*x1^2 + y1^2*x2^2
            sage: X = WehlerK3Surface([Z, Y])
            sage: Proj = ProjectiveSpace(QQ, 2)
            sage: P = Proj([1, 0, 0])
            sage: X.fiber(P, 1)
            Traceback (most recent call last):
            ...
            TypeError: fiber is degenerate

        ::

            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 - \
            2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - \
            4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: Proj = P[0]
            sage: T = Proj([0, 0, 1])
            sage: X.fiber(T, 1)
            [(0 : 0 : 1 , 0 : 1 : 0), (0 : 0 : 1 , 2 : 0 : 0)]

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], GF(7))
            sage: L = x0*y0 + x1*y1 - 1*x2*y2
            sage: Q=(2*x0^2 + x2*x0 + (2*x1^2 + x2^2))*y0^2 + ((x0^2 + x1*x0 +(x1^2 + 2*x2*x1 + x2^2))*y1 + \
            (2*x1^2 + x2*x1 + x2^2)*y2)*y0 + ((2*x0^2+ (x1 + 2*x2)*x0 + (2*x1^2 + x2*x1))*y1^2 + ((2*x1 + 2*x2)*x0 + \
            (x1^2 +x2*x1 + 2*x2^2))*y2*y1 + (2*x0^2 + x1*x0 + (2*x1^2 + x2^2))*y2^2)
            sage: W = WehlerK3Surface([L, Q])
            sage: W.fiber([4, 0, 1], 0)
            [(0 : 1 : 0 , 4 : 0 : 1), (4 : 0 : 2 , 4 : 0 : 1)]
        """
        R = self.base_ring()
        Zero = R(0)
        One = R(1)
        P = []
        for i in list(p):
            j = R(i)
            P.append(j)
        if component == 1:
            P0 = P + [Zero, Zero, Zero]
        else:
            P0 = [Zero, Zero, Zero] + P
        Points = []

        if (self.Gpoly(component,0)(P0)!= 0):
             #We are using the quadratic formula, we need this check to ensure that the points
             #will be rational
            T0 = (self.Hpoly(component, 0, 1)(P0)**2 -4*self.Gpoly(component, 0)(P0)*self.Gpoly(component, 1)(P0))
            T1 = (self.Hpoly(component, 0, 2)(P0)**2 -4*self.Gpoly(component, 0)(P0)*self.Gpoly(component, 2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = T0.sqrt()
                T1 = T1.sqrt()
                B1 = (-self.Hpoly(component, 0, 1)(P0)+T0)/(2*self.Gpoly(component, 0)(P0))
                B2 = (-self.Hpoly(component, 0, 1)(P0)-T0)/(2*self.Gpoly(component, 0)(P0))
                C1 = (-self.Hpoly(component, 0, 2)(P0)+T1)/(2*self.Gpoly(component, 0)(P0))
                C2 = (-self.Hpoly(component, 0, 2)(P0)-T1)/(2*self.Gpoly(component, 0)(P0))
                if component == 1:
                    Points.append(P+[One, B1, C1])
                    Points.append(P+[One, B2, C1])
                    Points.append(P+[One, B1, C2])
                    Points.append(P+[One, B2, C2])
                else:
                    Points.append([One, B1, C1]+P)
                    Points.append([One, B2, C1]+P)
                    Points.append([One, B1, C2]+P)
                    Points.append([One, B2, C2]+P)
            else:
                return []
        elif (self.Gpoly(component, 1)(P0) != 0):
            T0 = (self.Hpoly(component, 0, 1)(P0)**2 - 4*self.Gpoly(component, 0)(P0)*self.Gpoly(component, 1)(P0))
            T1 = (self.Hpoly(component, 1, 2)(P0)**2 - 4*self.Gpoly(component, 1)(P0)*self.Gpoly(component, 2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = T0.sqrt()
                T1 = T1.sqrt()
                A1 = (-self.Hpoly(component, 0, 1)(P0)+T0)/(2*self.Gpoly(component, 1)(P0))
                A2 = (-self.Hpoly(component, 0, 1)(P0)-T0)/(2*self.Gpoly(component, 1)(P0))
                C1 = (-self.Hpoly(component, 1, 2)(P0)+T1)/(2*self.Gpoly(component, 1)(P0))
                C2 = (-self.Hpoly(component, 1, 2)(P0)-T1)/(2*self.Gpoly(component, 1)(P0))
                if component == 1:
                    Points.append(P + [A1, One, C1])
                    Points.append(P + [A1, One, C2])
                    Points.append(P + [A2, One, C1])
                    Points.append(P + [A2, One, C2])
                else:
                    Points.append([A1, One, C1] + P)
                    Points.append([A1, One, C2] + P)
                    Points.append([A2, One, C1] + P)
                    Points.append([A2, One, C2] + P)
            else:
                return []
        elif self.Gpoly(component, 2)(P0) != 0:
            T0 = (self.Hpoly(component, 0, 2)(P0)**2 - 4*self.Gpoly(component, 0)(P0)*self.Gpoly(component, 2)(P0))
            T1 = (self.Hpoly(component, 1, 2)(P0)**2 - 4*self.Gpoly(component, 1)(P0)*self.Gpoly(component, 2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = T0.sqrt()
                T1 = T1.sqrt()
                A1 = (-self.Hpoly(component, 0, 2)(P0)+T0)/(2*self.Gpoly(component, 2)(P0))
                A2 = (-self.Hpoly(component, 0, 2)(P0)-T0)/(2*self.Gpoly(component, 2)(P0))
                B1 = (-self.Hpoly(component, 1, 2)(P0)+T1)/(2*self.Gpoly(component, 2)(P0))
                B2 = (-self.Hpoly(component, 1, 2)(P0)-T1)/(2*self.Gpoly(component, 2)(P0))
                if component == 1:
                    Points.append(P + [A1, B1, One])
                    Points.append(P + [A1, B2, One])
                    Points.append(P + [A2, B1, One])
                    Points.append(P + [A2, B2, One])
                else:
                    Points.append([A1, B1, One] + P)
                    Points.append([A1, B2, One] + P)
                    Points.append([A2, B1, One] + P)
                    Points.append([A2, B2, One] + P)
            else:
                return []
        elif self.Hpoly(component, 0, 1)(P0) != 0:
            if component == 1:
                Points.append(P+[Zero, One, Zero])
                Points.append(P+[-self.Hpoly(component, 0, 1)(P0),Zero,
                                 -self.Hpoly(component, 1, 2)(P0)])
                Points.append(P+[One,Zero,Zero])
                Points.append(P+[Zero,-self.Hpoly(component, 0, 1)(P0),
                                 -self.Hpoly(component, 0, 2)(P0)])
            else:
                Points.append([Zero,One,Zero]+P)
                Points.append([-self.Hpoly(component, 0, 1)(P0),Zero,
                               -self.Hpoly(component, 1, 2)(P0)] + P)
                Points.append([One,Zero,Zero]+P)
                Points.append([Zero,-self.Hpoly(component, 0, 1)(P0),
                               -self.Hpoly(component, 0, 2)(P0)] + P)
        elif self.Hpoly(component, 0, 2)(P0) != 0:
            if component == 1:
                Points.append(P+[Zero, Zero, One])
                Points.append(P+[-self.Hpoly(component, 0, 2)(P0),
                                 -self.Hpoly(component, 1, 2)(P0), Zero])
            else:
                Points.append([Zero, Zero, One]+P)
                Points.append([-self.Hpoly(component, 0, 2)(P0),
                               -self.Hpoly(component, 1, 2)(P0), Zero] + P)
        elif self.Hpoly(component, 1, 2)(P0) != 0:
            if component == 1:
                Points.append(P + [Zero, Zero, One])
                Points.append(P + [Zero, One, Zero])
            else:
                Points.append([Zero, Zero, One] + P)
                Points.append([Zero, One, Zero] + P)
        else:
            raise TypeError("fiber is degenerate")

        fiber = []
        for x in Points:
            if (self.L(x) == 0) and (self.Q(x) == 0):
                Y = self.point(x, False)
                if Y not in fiber:
                    fiber.append(Y)
        return fiber

    def nth_iterate_phi(self, P, n, **kwds):
        r"""
        Computes the nth iterate for the phi function.

        INPUT:

        - ``P`` -- - a point in `\mathbb{P}^2 \times \mathbb{P}^2`

        - ``n`` -- an integer

        kwds:

        - ``check`` - (default: ``True``) boolean checks to see if point is on the surface

        - ``normalize`` -- (default: ``False``) boolean normalizes the point

        OUTPUT: The nth iterate of the point given the phi function (if ``n`` is positive), or the
        psi function (if ``n`` is negative)

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 - x0*x1*y2^2
            sage: W = WehlerK3Surface([L ,Q])
            sage: T = W([-1, -1, 1, 1, 0, 1])
            sage: W.nth_iterate_phi(T, 7)
            (-1 : 0 : 1 , 1 : -2 : 1)

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 - x0*x1*y2^2
            sage: W = WehlerK3Surface([L, Q])
            sage: T = W([-1, -1, 1, 1, 0, 1])
            sage: W.nth_iterate_phi(T, -7)
            (1 : 0 : 1 , -1 : 2 : 1)

        ::

            sage: R.<x0,x1,x2,y0,y1,y2>=PolynomialRing(QQ, 6)
            sage: L = (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + (-y0^2 - y2*y1)*x2)*x0 \
            + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L, Q])
            sage: P = X([1, 0, -1, 1, -1, 0])
            sage: X.nth_iterate_phi(P, 8) == X.nth_iterate_psi(P, 8)
            True
        """
        try:
            n = ZZ(n)
        except TypeError:
            raise TypeError("iterate number must be an integer")
        #Since phi and psi are inverses and automorphisms
        if n < 0:
            return self.nth_iterate_psi(P, abs(n), **kwds)
        if n == 0:
            return self
        else:
            Q = self.phi(P, **kwds)
            for i in range(2, n+1):
                Q = self.phi(Q, **kwds)
            return Q

    def nth_iterate_psi(self, P, n, **kwds):
        r"""
        Computes the nth iterate for the psi function.

        INPUT:

        - ``P`` -- - a point in `\mathbb{P}^2 \times \mathbb{P}^2`

        - ``n`` -- an integer

        kwds:

        - ``check`` -- (default: ``True``) boolean, checks to see if point is on the surface

        - ``normalize`` -- (default: ``False``) boolean, normalizes the point

        OUTPUT: The nth iterate of the point given the psi function (if ``n`` is positive),
        or the phi function (if ``n`` is negative)

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 - x0*x1*y2^2
            sage: W = WehlerK3Surface([L, Q])
            sage: T = W([-1, -1, 1, 1, 0, 1])
            sage: W.nth_iterate_psi(T, -7)
            (-1 : 0 : 1 , 1 : -2 : 1)

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 - x0*x1*y2^2
            sage: W = WehlerK3Surface([L, Q])
            sage: T = W([-1, -1, 1, 1, 0, 1])
            sage: W.nth_iterate_psi(T, 7)
            (1 : 0 : 1 , -1 : 2 : 1)
        """
        try:
            n = ZZ(n)
        except TypeError:
            raise TypeError("iterate number must be an integer")
        #Since phi and psi and inverses
        if n < 0:
            return self.nth_iterate_phi(P, abs(n), **kwds)
        if n == 0:
            return self
        else:
            Q = self.psi(P, **kwds)
            for i in range(2, n+1):
                Q = self.psi(Q, **kwds)
            return Q

    def orbit_phi(self, P, N, **kwds):
        r"""
        Return the orbit of the `\phi` function defined by
        `\phi = \sigma_y \circ \sigma_x`.

        This function is defined in [CS1996]_.

        INPUT:

        - ``P`` - Point on the K3 surface

        - ``N`` - a non-negative integer or list or tuple of two non-negative integers

        kwds:

        - ``check`` -- (default: ``True``) boolean, checks to see if point is on the surface

        - ``normalize`` -- (default: ``False``) boolean, normalizes the point

        OUTPUT: List of points in the orbit

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - \
            4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + \
            x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP(0, 0, 1, 1, 0, 0)
            sage: X.orbit_phi(T,2, normalize = True)
            [(0 : 0 : 1 , 1 : 0 : 0), (-1 : 0 : 1 , 0 : 1 : 0), (-12816/6659 : 55413/6659 : 1 , 1 : 1/9 : 1)]
            sage: X.orbit_phi(T,[2,3], normalize = True)
            [(-12816/6659 : 55413/6659 : 1 , 1 : 1/9 : 1),
            (7481279673854775690938629732119966552954626693713001783595660989241/18550615454277582153932951051931712107449915856862264913424670784695
            : 3992260691327218828582255586014718568398539828275296031491644987908/18550615454277582153932951051931712107449915856862264913424670784695 :
            1 , -117756062505511/54767410965117 : -23134047983794359/37466994368025041 : 1)]
        """
        if not isinstance(N, (list, tuple)):
            N = [0, N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError("orbit bounds must be integers")
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be non-negative")
        if N[0] > N[1]:
            return []
        Q = self(copy(P))
        for i in range(1, N[0] + 1):
            Q = self.phi(Q, **kwds)
        Orb = [Q]
        for i in range(N[0] + 1, N[1] + 1):
            Q = self.phi(Q, **kwds)
            Orb.append(Q)
        return Orb

    def orbit_psi(self, P, N, **kwds):
        r"""
        Return the orbit of the `\psi` function defined by
        `\psi = \sigma_x \circ \sigma_y`.

        This function is defined in [CS1996]_.

        INPUT:

        - ``P`` - a point on the K3 surface

        - ``N`` - a non-negative integer or list or tuple of two non-negative integers

        kwds:

        - ``check`` - (default: ``True``) boolean, checks to see if point is on the surface

        - ``normalize`` -- (default: ``False``) boolean, normalizes the point

        OUTPUT: a list of points in the orbit

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
             3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 - \
             4*x1*x2*y1^2 + 5*x0*x2*y0*y2 -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + \
              x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = X(0, 0, 1, 1, 0, 0)
            sage: X.orbit_psi(T, 2, normalize = True)
            [(0 : 0 : 1 , 1 : 0 : 0), (0 : 0 : 1 , 0 : 1 : 0), (-1 : 0 : 1 , 1 : 1/9 : 1)]
            sage: X.orbit_psi(T,[2,3], normalize = True)
            [(-1 : 0 : 1 , 1 : 1/9 : 1),
            (-12816/6659 : 55413/6659 : 1 , -117756062505511/54767410965117 : -23134047983794359/37466994368025041 : 1)]
        """
        if not isinstance(N, (list, tuple)):
            N = [0, N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError("orbit bounds must be integers")
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be non-negative")
        if N[0] > N[1]:
            return []
        Q = self(copy(P))
        for i in range(1, N[0] + 1):
            Q = self.psi(Q, **kwds)
        Orb = [Q]
        for i in range(N[0] + 1, N[1] + 1):
            Q = self.psi(Q, **kwds)
            Orb.append(Q)
        return Orb

    def is_isomorphic(self, right):
        r"""
        Checks to see if two K3 surfaces have the same defining ideal.

        INPUT:

        - ``right`` - the K3 surface to compare to the original

        OUTPUT: Boolean

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            -4*x1*x2*y1^2 + 5*x0*x2*y0*y2 - 4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: W = WehlerK3Surface([Z + Y^2, Y])
            sage: X.is_isomorphic(W)
            True

        ::

            sage: R.<x,y,z,u,v,w> = PolynomialRing(QQ, 6)
            sage: L = x*u-y*v
            sage: Q = x*y*v^2 + z^2*u*w
            sage: W1 = WehlerK3Surface([L, Q])
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 -x0*x1*y2^2
            sage: W2 = WehlerK3Surface([L, Q])
            sage: W1.is_isomorphic(W2)
            False
        """
        return self.defining_ideal() == right.defining_ideal()

    def is_symmetric_orbit(self,orbit):
        r"""
        Checks to see if the orbit is symmetric (i.e. if one of the points on the
        orbit is fixed by '\sigma_x' or '\sigma_y').

        INPUT:

        - ``orbit``- a periodic cycle of either psi or phi

        OUTPUT: Boolean

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], GF(7))
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + 3*x0*x1*y0*y1 \
            -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 -4*x1*x2*y1^2 + 5*x0*x2*y0*y2 \
            -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: T = PP([0, 0, 1, 1, 0, 0])
            sage: orbit = X.orbit_psi(T, 4)
            sage: X.is_symmetric_orbit(orbit)
            True

        ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: L = x0*y0 + x1*y1 + x2*y2
            sage: Q = x1^2*y0^2 + 2*x2^2*y0*y1 + x0^2*y1^2 - x0*x1*y2^2
            sage: W = WehlerK3Surface([L, Q])
            sage: T = W([-1, -1, 1, 1, 0, 1])
            sage: Orb = W.orbit_phi(T, 7)
            sage: W.is_symmetric_orbit(Orb)
            False
        """
        N = len(orbit)
        if self.nth_iterate_phi(orbit[0], N) != orbit[0] and self.nth_iterate_psi(orbit[0], N) != orbit[0]:
            raise ValueError("must be an orbit of phi or psi functions")
        sym = False
        i = 0
        while i < len(orbit) - 1 and not sym:
            P = orbit[i]
            Q = orbit[i + 1]
            if P[0] == Q[0] or P[1] == Q[1]:
                sym = True
            i += 1
        return sym


class WehlerK3Surface_field( WehlerK3Surface_ring):
    pass


class WehlerK3Surface_finite_field( WehlerK3Surface_field):
    def cardinality( self):
        r"""
        Counts the total number of points on the K3 surface.

        ALGORITHM:

        Enumerate points over `\mathbb{P}^2`, and then count the points on the fiber of
        each of those points.

        OUTPUT: Integer - total number of points on the surface

        EXAMPLES::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], GF(7))
            sage: Z = x0^2*y0^2 + 3*x0*x1*y0^2 + x1^2*y0^2 + 4*x0^2*y0*y1 + \
            3*x0*x1*y0*y1 -2*x2^2*y0*y1 - x0^2*y1^2 + 2*x1^2*y1^2 - x0*x2*y1^2 \
            - 4*x1*x2*y1^2 + 5*x0*x2*y0*y2 -4*x1*x2*y0*y2 + 7*x0^2*y1*y2 + 4*x1^2*y1*y2 \
            + x0*x1*y2^2 + 3*x2^2*y2^2
            sage: Y = x0*y0 + x1*y1 + x2*y2
            sage: X = WehlerK3Surface([Z, Y])
            sage: X.cardinality()
            55
        """
        def getPx1():
            return ([x, y, 1] for x in self.base_ring() for y in self.base_ring())

        def getPx2():
            return ([x, 1, 0] for x in self.base_ring())
        Count = 0
        Xpoint = [1, 0, 0]
        Ypoint = [1, 0, 0]
        # Create all possible Px1 Values
        for i in getPx1():
            for j in getPx1():
                A = i + j
                if self.L(A) == 0 and self.Q(A) == 0:
                    Count += 1
            for k in getPx2():
                A = i + k
                if self.L(A) == 0 and self.Q(A) == 0:
                    Count += 1
            B = i + Ypoint
            if self.L(B) == 0 and self.Q(B) == 0:
                Count += 1
        #Create all possible Px2 Values
        for i in getPx2():
            for j in getPx1():
                A = i + j
                if self.L(A) == 0 and self.Q(A) == 0:
                    Count += 1
            for k in getPx2():
                A = i + k
                if (self.L(A) == 0 and self.Q(A) == 0):
                    Count += 1
            B = i + Ypoint
            if (self.L(B) == 0 and self.Q(B) == 0):
                Count += 1
        #Create all Xpoint values
        for j in getPx1():
            A = Xpoint+j
            if (self.L(A) == 0 and self.Q(A) == 0):
                Count += 1
        for k in getPx2():
            B = Xpoint + k
            if (self.L(B) == 0 and self.Q(B) == 0):
                Count += 1
            C = Xpoint + Ypoint
        if (self.L(C) == 0 and self.Q(C) == 0):
            Count += 1
        return Count
