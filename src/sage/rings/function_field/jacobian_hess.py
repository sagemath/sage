r"""
Jacobians in Hess model

This module implements Jacobian arithmetic based on divisor representation by
ideals. This approach to Jacobian arithmetic implementation is attributed to
Hess [Hes2002]_.

Jacobian
--------

To create a Jacobian in Hess model, specify ``'hess'`` model and provide a base divisor
of degree `g`, which is the genus of the function field::

    sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
    sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
    sage: C.geometric_genus()
    1
    sage: B = C([0,1,0]).place()
    sage: B.degree()
    1
    sage: J = C.jacobian(model='hess', base_div=B)
    sage: J
    Jacobian of Projective Plane Curve over Finite Field of size 29
     defined by x^3 - y^2*z + 5*z^3 (Hess model)

Group of rational points
------------------------

The group of rational points of a Jacobian is created from the Jacobian. A
point of the Jacobian group is determined by a divisor (of degree zero) of the
form `D - B` where `D` is an effective divisor of degree `g` and `B` is the base
divisor. Hence a point of the Jacobian group is represented by `D`.

::

    sage: G = J.group()
    sage: P1 = C([1,8,1]).place()
    sage: P2 = C([2,10,1]).place()
    sage: p1 = G(P1)
    sage: p2 = G(P2)
    sage: p1
    [Place (y + 21, z + 28)]
    sage: p2
    [Place (y + 24, z + 14)]
    sage: p1 + p2
    [Place (y + 8, z + 28)]

AUTHORS:

- Kwankyu Lee (2022-01-24): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import op_EQ, richcmp

from sage.categories.map import Map
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.homset import Hom

from sage.arith.misc import integer_ceil
from sage.arith.functions import lcm

from sage.rings.integer import Integer
from sage.matrix.constructor import matrix

from sage.combinat.integer_vector_weighted import WeightedIntegerVectors

from .place import FunctionFieldPlace
from .divisor import FunctionFieldDivisor

from .jacobian_base import (Jacobian_base,
                            JacobianGroup_base,
                            JacobianGroup_finite_field_base,
                            JacobianPoint_base,
                            JacobianPoint_finite_field_base)


class JacobianPoint(JacobianPoint_base):
    """
    Points of Jacobians represented by a pair of ideals.

    If a point of Jacobian is determined by `D`, then the divisor `D` is
    represented by a pair of ideals in the finite maximal order and the
    infinite maximal order of the function field.

    For efficiency reasons, the actual ideals stored are the inverted ideals
    of the ideals representing the divisor `D`.

    INPUT:

    - ``parent`` -- Jacobian group

    - ``dS`` -- an ideal of the finite maximal order of a function field

    - ``ds`` -- an ideal of infinite maximal order of a function field

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: G = C.jacobian(model='hess', base_div=b).group()
        sage: pl = C([1,8,1]).place()
        sage: p = G.point(pl - b)
        sage: dS, ds = p._data
        sage: -(dS.divisor() + ds.divisor()) == pl
        True
    """
    def __init__(self, parent, dS, ds):
        """
        Initialize.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl = C([1,8,1]).place()
            sage: p = G.point(pl - b)
            sage: TestSuite(p).run(skip=['_test_category','_test_pickling'])
        """
        super().__init__(parent)
        self._data = (dS, ds)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: G.zero()
            [Place (1/y, 1/y*z)]
        """
        dS, ds = self._data
        divisor = (~dS).divisor() + (~ds).divisor()
        return f'[{divisor}]'

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: f = x/(y + 1)
            sage: d = f.divisor()
            sage: {d: 1}
            {Place (1/x, 1/x^4*y^2 + 1/x^2*y + 1)
              + Place (1/x, 1/x^2*y + 1)
              + 3*Place (x, (1/(x^3 + x^2 + x))*y^2)
              - 6*Place (x + 1, y + 1): 1}
        """
        return hash(self._data)

    def _richcmp_(self, other, op):
        """
        Compare ``self`` with ``other`` with respect to operator ``op``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([2,10,1]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 == p1
            True
            sage: p1 != p2
            True
            sage: p1 > p1
            False
            sage: p1 > p2
            False
            sage: p1 < p2
            True
        """
        if op is op_EQ:
            J = self.parent()
            idS, ids = self._data
            jdS, jds = other._data
            return J._normalize(idS / jdS, ids / jds) is not None
        else:
            return richcmp(self._data, other._data, op)

    def _add_(self, other):
        """
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([2,10,1]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 + p2
            [Place (y + 8, z + 3)]
            sage: p1 + p2 == p2 + p1
            True
        """
        G = self.parent()
        idS, ids = self._data
        jdS, jds = other._data
        bdS, bds = G._base_point
        dS, ds = G._normalize(idS * jdS * bdS, ids * jds * bds)
        return G.element_class(self.parent(), dS, ds)

    def _neg_(self):
        """
        Return the negative of this point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: -p
            [Place (y + 27, z + 1)]
            sage: -(-p) == p
            True
        """
        G = self.parent()
        idS, ids = self._data
        bdS2, bds2 = G._base_point_double
        dS, ds = G._normalize(~(idS * bdS2), ~(ids * bds2))
        return G.element_class(self.parent(), dS, ds)

    def multiple(self, n):
        """
        Return the ``n``-th multiple of this point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: p.multiple(100)
            [Place (1/y, 1/y*z + 8)]
        """
        if n == 0:
            return self.parent().zero()

        G = self.parent()
        idS, ids = self._data
        bdS, bds = G._base_point
        bdS2, bds2 = G._base_point_double
        idSbdS2 = idS * bdS2
        idsbds2 = ids * bds2

        if n < 0:
            bits = Integer(-n).digits(2)
        else:
            bits = Integer(n).digits(2)
        bits.pop()

        dS = idS
        ds = ids
        for i in range(len(bits)):
            b = bits.pop()
            if b > 0:
                dS, ds = G._normalize(dS * dS * idSbdS2 , ds * ds * idsbds2)
            else:
                dS, ds = G._normalize(dS * dS * bdS, ds * ds * bds)
        if n < 0:
            dS, ds = G._normalize(~(dS * bdS2), ~(ds * bds2))

        return G.element_class(self.parent(), dS, ds)

    def addflip(self, other):
        """
        Return the addflip of this and ``other`` point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([2,19,1]).place()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1.addflip(p2)
            [Place (y + 8, z + 27)]
            sage: _ == -(p1 + p2)
            True
        """
        return -(self + other)

    def defining_divisor(self):
        """
        Return the effective divisor that defines this point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: p.defining_divisor() == pl
            True
        """
        dS, ds = self._data
        return (~dS).divisor() + (~ds).divisor()

    def order(self, bound=None):
        """
        Return the order of this point.

        ALGORITHM: Shanks' Baby Step Giant Step

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: p = C([-1,2,1]).place()
            sage: pt = G.point(p - b)
            sage: pt.order()
            30
        """
        if bound is None:  # naive
            J = self.parent()
            zero = J.zero()

            m = self
            r = 1
            while m != zero:
                m = m + self
                r += 1
            return r

        # if bound is given, deploy Shanks' Baby Step Giant Step

        J = self.parent()
        B = J.bound_on_order()
        q = integer_ceil(B.sqrt())
        zero = J.zero()

        # baby steps
        b = [zero]
        g = self
        for i in range(q - 1):
            if g == zero:
                return i + 1
            b.append(g)
            g = g + self

        # giant steps
        g0 = (-q)*(self)
        g = g0
        for i in range(q - 1):
            for r in range(q):
                if g == b[r]:
                    return q * (i + 1) + r
            g = g + g0

        # order is neither smaller or nor larger than this
        return q**2

    def divisor(self):
        """
        Return the divisor representing this point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: G.point(p.divisor()) == p
            True
        """
        J = self.parent()
        dS, ds = self._data
        return (~dS).divisor() + (~ds).divisor() - J._base_div


class JacobianPoint_finite_field(JacobianPoint, JacobianPoint_finite_field_base):
    """
    Points of Jacobians over finite fields
    """
    pass


class JacobianGroupEmbedding(Map):
    """
    Embeddings between Jacobian groups.

    INPUT:

    - ``base_group`` -- Jacobian group over a base field

    - ``extension_group`` -- Jacobian group over an extension field

    EXAMPLES::

        sage: k = GF(17)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: J = C.jacobian(model='hess', base_div=b)
        sage: G1 = J.group()
        sage: K = k.extension(3)
        sage: G3 = J.group(K)
        sage: G3.coerce_map_from(G1)
        Jacobian group embedding map:
          From: Group of rational points of Jacobian
           over Finite Field of size 17 (Hess model)
          To:   Group of rational points of Jacobian
           over Finite Field in z3 of size 17^3 (Hess model)
    """
    def __init__(self, base_group, extension_group):
        """
        Initialize.

        TESTS::

            sage: k = GF(17)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G1 = J.group()
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: map = G3.coerce_map_from(G1)
            sage: TestSuite(map).run(skip=['_test_category', '_test_pickling'])
        """
        F = base_group._function_field
        F_base = F.base_field()
        K = F.constant_base_field()

        F_ext = extension_group._function_field
        F_ext_base = F_ext.base_field()
        K_ext = F_ext.constant_base_field()

        # construct embedding of F into F_ext
        embedK = K_ext.coerce_map_from(K)
        embedF_base = F_base.hom(F_ext_base.gen(), embedK)
        if F.degree() > 1:
            embedF = F.hom(F_ext.gen(), embedF_base)
        else:
            embedF = embedF_base

        self._embedF = embedF
        self._O_ext = F_ext.maximal_order()
        self._Oinf_ext = F_ext.maximal_order_infinite()

        Map.__init__(self, Hom(base_group, extension_group, CommutativeAdditiveGroups()))

    def _repr_type(self):
        """
        Return string representation of ``self``.

        TESTS::

            sage: k = GF(17)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G1 = J.group()
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: G3.coerce_map_from(G1)  # indirect doctest
            Jacobian group embedding map:
              From: Group of rational points of Jacobian
               over Finite Field of size 17 (Hess model)
              To:   Group of rational points of Jacobian
               over Finite Field in z3 of size 17^3 (Hess model)
        """
        return 'Jacobian group embedding'

    def _call_(self, x):
        """
        Conorm map from F to F_ext.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^5 - x^3 - 2*x - 1).projective_closure()
            sage: J = C.jacobian(model='hess')
            sage: G1 = J.group()
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: m = G3.coerce_map_from(G1)
            sage: m(G1.zero()) == G3.zero()
            True
        """
        embedF = self._embedF
        O_ext = self._O_ext
        Oinf_ext = self._Oinf_ext

        idS,ids = x._data
        dS = O_ext.ideal([embedF(g) for g in idS.gens()])
        ds = Oinf_ext.ideal([embedF(g) for g in ids.gens()])
        return self.codomain().element_class(self.codomain(), dS, ds)


class JacobianGroup(UniqueRepresentation, JacobianGroup_base):
    """
    Groups of rational points of a Jacobian.

    INPUT:

    - ``parent`` -- a Jacobian

    - ``function_field`` -- a function field

    - ``base_div`` -- an effective divisor of the function field

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: J = C.jacobian(model='hess', base_div=b)
        sage: J.group()
        Group of rational points of Jacobian
         over Finite Field of size 17 (Hess model)
    """
    Element = JacobianPoint
    _embedding_map_class = JacobianGroupEmbedding

    def __init__(self, parent, function_field, base_div):
        """
        Initialize.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: TestSuite(G).run(skip=['_test_elements', '_test_pickling'])
        """
        super().__init__(parent, function_field, base_div)

        bdS, bds = self._get_dS_ds(-base_div)
        try:
            bdS._gens_two()  # speed up multiplication with these ideals
            bds._ideal._gens_two()  # by storing vector forms of two generators
        except AttributeError:
            pass
        self._base_point = (bdS, bds)
        self._base_point_double = (bdS * bdS, bds * bds)

        self._base_place = None

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: J.group()
            Group of rational points of Jacobian
             over Finite Field of size 17 (Hess model)
        """
        r = super()._repr_()
        return r + ' (Hess model)'

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        If ``x`` is an effective divisor, then it is assumed to be of
        degree `g`, the genus of the function field.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='hess', base_div=b).group()
            sage: G(0)
            [Place (1/y, 1/y*z)]
        """
        if x == 0:
            return self.zero()

        if isinstance(x, FunctionFieldPlace):
            if (self._base_place is not None
                and x in self._function_field.place_set()
                and x.degree() == 1):
                x = x - self._base_place
            else:
                x = x.divisor()

        if (isinstance(x, FunctionFieldDivisor)
            and x in self._function_field.divisor_group()):
            if x.degree() == 0:
                return self.point(x)
            if x.is_effective():
                if x.degree() != self._genus:
                    raise ValueError(f"effective divisor is not of degree {self._genus}")
                return self.element_class(self, *self._get_dS_ds(x))

        raise ValueError(f"cannot construct a point from {x}")

    def _get_dS_ds(self, divisor):
        """
        Return (dS,ds) representation of the divisor.

        TESTS::

            sage: k = GF(17)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: pl = C([2,8,1]).place()
            sage: dS, ds = G._get_dS_ds(2*pl)
            sage: (~dS).divisor() + (~ds).divisor() == 2*pl
            True
        """
        F = self._function_field
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(1)
        J = Oinf.ideal(1)
        for p in divisor._data:
            m = divisor._data[p]
            if p.is_infinite_place():
                J *= p.prime_ideal() ** (-m)
            else:
                I *= p.prime_ideal() ** (-m)

        return I, J

    def _normalize(self, I, J):
        """
        Return a pair of normalized ideals from `I` and `J`.

        INPUT:

        - ``I`` -- an ideal of the finite maximal order

        - ``J`` -- an ideal of the infinite maximal order

        The output represents an effective divisor linearly equivalent to the
        divisor represented by the given ideals `I` and `J`.

        ALGORITHM:

        Computes a function `f` in the Riemann-Roch space of the divisor `D`
        represented by the (inverted) ideals `I` and `J`. The output is the
        pair of the (inverted) ideals representing the effective divisor `(f) + D`,
        which is linearly equivalent to `D`.

        TESTS::

            sage: k = GF(17)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: pl = C([2,8,1]).place()
            sage: p = G.point(pl - b)
            sage: dS, ds = (p + p)._data  # indirect doctest
            sage: G.point((~dS).divisor() + (~ds).divisor() - b) == p + p
            True
        """
        F = self._function_field
        n = F.degree()

        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        # Step 1: construct matrix M of rational functions in x such that
        # M * B == C where B = [b1,b1,...,bn], C =[v1,v2,...,vn]
        V,fr,to = F.free_module(map=True)
        B = matrix([to(b) for b in J.gens_over_base()])
        C = matrix([to(v) for v in I.gens_over_base()])
        M = C * B.inverse()

        # Step 2: get the denominator d of M and set mat = d * M
        den = lcm([e.denominator() for e in M.list()])
        R = den.parent() # polynomial ring
        one = R.one()
        mat = matrix(R, n, [e.numerator() for e in (den*M).list()])
        gens = list(I.gens_over_base())

        # Step 3: transform mat to a weak Popov form, together with gens

        # initialise pivot_row and conflicts list
        found = None
        pivot_row = [[] for i in range(n)]
        conflicts = []
        for i in range(n):
            bestp = -1
            best = -1
            for c in range(n):
                d = mat[i,c].degree()
                if d >= best:
                    bestp = c
                    best = d

            if best <= den.degree():
                found = i
                break

            if best >= 0:
                pivot_row[bestp].append((i,best))
                if len(pivot_row[bestp]) > 1:
                    conflicts.append(bestp)

        if found is None:
            # while there is a conflict, do a simple transformation
            while conflicts:
                c = conflicts.pop()
                row = pivot_row[c]
                i,ideg = row.pop()
                j,jdeg = row.pop()

                if jdeg > ideg:
                    i,j = j,i
                    ideg,jdeg = jdeg,ideg

                coeff = - mat[i,c].lc() / mat[j,c].lc()
                s = coeff * one.shift(ideg - jdeg)

                mat.add_multiple_of_row(i, j, s)
                gens[i] += s * gens[j]

                row.append((j,jdeg))

                bestp = -1
                best = -1
                for c in range(n):
                    d = mat[i,c].degree()
                    if d >= best:
                        bestp = c
                        best = d

                if best <= den.degree():
                    found = i
                    break

                if best >= 0:
                    pivot_row[bestp].append((i,best))
                    if len(pivot_row[bestp]) > 1:
                        conflicts.append(bestp)
            else:
                return None

        f = gens[found]
        return (O.ideal(~f) * I, Oinf.ideal(~f) * J)

    def point(self, divisor):
        """
        Return the point represented by the divisor of degree zero.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: p = C([-1,2,1]).place()
            sage: G.point(p - b)
            [Place (y + 2, z + 1)]
        """
        if divisor.degree() != 0:
            raise ValueError('divisor not of degree zero')
        c = divisor + self._base_div
        f = c.basis_function_space()[0]
        d = f.divisor() + c
        dS, ds = self._get_dS_ds(d)
        return self.element_class(self, dS, ds)

    @cached_method
    def zero(self):
        """
        Return the zero element of this group.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: G.zero()
            [Place (1/y, 1/y*z)]
        """
        bdS,bds = self._base_point
        return self.element_class(self, ~bdS, ~bds)


class JacobianGroup_finite_field(JacobianGroup, JacobianGroup_finite_field_base):
    """
    Jacobian groups of function fields over finite fields.

    INPUT:

    - ``parent`` -- a Jacobian

    - ``function_field`` -- a function field

    - ``base_div`` -- an effective divisor of the function field

    EXAMPLES::

        sage: k = GF(17)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: J = C.jacobian(model='hess', base_div=b)
        sage: G1 = J.group()
        sage: K = k.extension(3)
        sage: G3 = J.group(K)
        sage: G3.coerce_map_from(G1)
        Jacobian group embedding map:
          From: Group of rational points of Jacobian
           over Finite Field of size 17 (Hess model)
          To:   Group of rational points of Jacobian
           over Finite Field in z3 of size 17^3 (Hess model)
    """
    Element = JacobianPoint_finite_field

    def __init__(self, parent, function_field, base_div):
        """
        Initialize.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: TestSuite(G).run(skip=['_test_elements','_test_pickling'])
        """
        super().__init__(parent, function_field, base_div)

        F = self._function_field
        K = F.constant_base_field()

        R = F.base_field()  # base rational function field
        x = R.gen()
        y = F.gen()

        r = self._parent._function_field.constant_base_field().degree()
        frob_K = K.frobenius_endomorphism(r)
        frob_R = R.hom(x, base_morphism=frob_K)
        frob_F = F.hom(y, base_morphism=frob_R)

        self._frobenius = frob_F

    def __iter__(self):
        """
        Return generator of points of this group.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='hess')
            sage: G = J.group()
            sage: len([pt for pt in G])
            11
        """
        g = self._parent._function_field.genus()
        F = self._function_field
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        deg = 1
        support = []
        degrees = []
        multiples = []
        lst = []

        places_infinite = F.places_infinite()
        generators = [iter(places_infinite)]
        num_of_infinite_places = len(places_infinite)
        while True:
            while True:
                try:
                    new_pl = next(generators[-1])
                    break
                except StopIteration:
                    if deg > g:
                        return
                    generators.append(F._places_finite(deg))
                    deg += 1
            multiples.append((g + 1)*[None])
            P = ~new_pl.prime_ideal()
            dn = new_pl.degree()
            I0 = O.ideal(1)
            J0 = Oinf.ideal(1)
            dr = 0
            for r in range(1, g // new_pl.degree() + 1):
                if new_pl.is_infinite_place():
                    J0 = J0 * P
                else:
                    I0 = I0 * P
                multiples[-1][r] = (I0, J0)
                dr = dr + dn
                for weights in WeightedIntegerVectors(g - dr, degrees):
                    I = I0
                    J = J0
                    for i in range(len(support)):
                        w = weights[i]
                        if w > 0:
                            dS, ds = multiples[i][w]
                            if i < num_of_infinite_places:
                                J *= ds  # dS is the unit ideal
                            else:
                                I *= dS  # ds is the unit ideal
                    pt = self.element_class(self, I, J)
                    if pt not in lst:
                        lst.append(pt)
                        yield pt
            support.append(new_pl)
            degrees.append(new_pl.degree())

    def _frobenius_on(self, pt):
        """
        Return the image of ``pt`` acted by the Frobenius automorphism.

        EXAMPLES::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='hess')
            sage: G1 = J.group()
            sage: G1.order()
            11
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: pts1 = G1.get_points(11)
            sage: pts3 = G3.get_points(12)
            sage: pt = next(pt for pt in pts3 if pt not in pts1)
            sage: pt.frobenius().frobenius().frobenius() == pt  # indirect doctest
            True
            sage: pt.frobenius() == pt
            False
        """
        frob_F = self._frobenius

        F = self._function_field
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        idS,ids = pt._data
        dS = O.ideal([frob_F(g) for g in idS.gens()])
        ds = Oinf.ideal([frob_F(g) for g in ids.gens()])
        return self.element_class(self, dS, ds)


class Jacobian(Jacobian_base, UniqueRepresentation):
    """
    Jacobians of function fields.

    EXAMPLES::

        sage: k = GF(17)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: C.jacobian(model='hess', base_div=b)
        Jacobian of Projective Plane Curve over Finite Field of size 17
         defined by x^3 - y^2*z + 5*z^3 (Hess model)
    """
    def __init__(self, function_field, base_div, **kwds):
        """
        Initialize.

        TESTS::

            sage: k = GF(17)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: TestSuite(J).run(skip=['_test_elements','_test_pickling'])
        """
        super().__init__(function_field, base_div, **kwds)

        if function_field.constant_base_field().is_finite():
            self._group_class = JacobianGroup_finite_field
        else:
            self._group_class = JacobianGroup

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: C.jacobian(model='hess', base_div=b)
            Jacobian of Projective Plane Curve over Finite Field of size 17
             defined by x^3 - y^2*z + 5*z^3 (Hess model)
        """
        r = super()._repr_()
        return r + ' (Hess model)'
