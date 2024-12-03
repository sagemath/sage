r"""
Jacobians of function fields

This module provides base classes for Jacobians of function fields.

Jacobian
--------

The Jacobian of a function field is created by default in the Hess model, with
a base divisor of degree `g` the genus of the function field. The base divisor
is automatically chosen if not given. ::

    sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
    sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
    sage: F = C.function_field()
    sage: J = F.jacobian()
    sage: J
    Jacobian of Function field in z defined by z^3 + 23*y^2*z + 6 (Hess model)
    sage: J.base_divisor().degree()
    1

Explicitly specify a model if you want Jacobians in different models. ::

    sage: J_km = F.jacobian(model='km_large')
    sage: J_km
    Jacobian of Function field in z defined by z^3 + 23*y^2*z + 6 (Khuri-Makdisi large model)

Group of rational points
------------------------

The group of rational points of a Jacobian is created from the Jacobian. A
point of the Jacobian group is determined by a divisor of degree zero. To
represent the point, a divisor of the form `D-B` is selected where `D` is an
effective divisor of the same degree with the base divisor `B`. Hence the point
is simply represented by the divisor `D`. ::

    sage: G = J.group()
    sage: G.order()
    30
    sage: pl1 = C([1,8,1]).place()
    sage: pl2 = C([2,10,1]).place()
    sage: p1 = G.point(pl1 - pl2)
    sage: p1
    [Place (y + 1, z + 6)]
    sage: p2 = G.point(pl2 - pl1)
    sage: p2
    [Place (y + 28, z + 6)]
    sage: p1 + p2 == G.zero()
    True
    sage: p1.order()
    5

We can get the corresponding point in the Jacobian in a different model. ::

    sage: # long time
    sage: p1km = J_km(p1)
    sage: p1km.order()
    5
    sage: p1km
    Point of Jacobian determined by
    [ 1  0  0  0  0  0 11  0  0]
    [ 0  1  0  0  0  0 18  0  0]
    [ 0  0  1  0  0  0 11  0  0]
    [ 0  0  0  1  0  0 18  1  0]
    [ 0  0  0  0  1  0 25  0 19]
    [ 0  0  0  0  0  1  8  8  0]

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

import math

from sage.arith.misc import integer_floor, integer_ceil

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement

from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.schemes import Jacobians
from sage.categories.pushout import ConstructionFunctor, pushout

from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer


class JacobianPoint_base(ModuleElement):
    """
    Abstract base class of points of Jacobian groups.
    """
    pass


class JacobianPoint_finite_field_base(JacobianPoint_base):
    """
    Points of Jacobians over finite fields.
    """
    def order(self):
        """
        Return the order of this point.

        EXAMPLES::

            sage: # long time
            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: F = C.function_field()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = F.get_place(1)
            sage: pl = C([-1,2,1]).place()
            sage: p = G.point(pl - b)
            sage: p.order()
            15

        ALGORITHM: Shanks' Baby Step Giant Step
        """
        G = self.parent()
        B = G._bound_on_order()
        q = integer_ceil(B.sqrt())
        zero = G.zero()

        # baby steps
        b = [zero]
        g = self
        for i in range(q - 1):
            if g == zero:
                return i + 1
            b.append(g)
            g = g + self

        # giant steps
        g0 = self.multiple(-q)
        g = g0
        for i in range(q - 1):
            for r in range(q):
                if g == b[r]:
                    return q * (i + 1) + r
            g = g + g0

        # order is neither smaller or nor larger than this
        return q**2

    def frobenius(self):
        """
        Return the image of the point acted by the Frobenius automorphism.

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
            sage: pt.frobenius() == pt
            False
            sage: pt.frobenius().frobenius().frobenius() == pt
            True
        """
        G = self.parent()
        return G._frobenius_on(self)


class JacobianGroupFunctor(ConstructionFunctor):
    """
    A construction functor for Jacobian groups.

    EXAMPLES::

        sage: k = GF(7)
        sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: J = C.jacobian(model='hess')
        sage: G = J.group()
        sage: F, obj = G.construction()
        sage: F
        JacobianGroupFunctor
    """
    rank = 20

    def __init__(self, base_field, field):
        """
        Initialize.

        TESTS::

            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: K = k.extension(2)
            sage: G = J.group(K)
            sage: F, obj = G.construction()
            sage: TestSuite(F).run()
        """
        super().__init__(Jacobians(base_field), CommutativeAdditiveGroups())

        self._field = field

    def _apply_functor(self, jacobian):
        """
        Apply this functor to ``jacobian``.

        INPUT:

        - ``jacobian`` -- a Jacobian

        EXAMPLES::

            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: K = k.extension(2)
            sage: G = J.group(K)
            sage: F, obj = G.construction()
            sage: F(obj) is G  # indirect doctest
            True
        """
        return jacobian.group(self._field)

    def merge(self, other):
        """
        Return the functor merging ``self`` and ``other``.

        INPUT:

        - ``other`` -- a functor

        EXAMPLES::

            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: K2 = k.extension(2)
            sage: G2 = J.group(K2)
            sage: K3 = k.extension(3)
            sage: G3 = J.group(K3)
            sage: sage.categories.pushout.pushout(G2, G3)  # indirect doctest
            Group of rational points of Jacobian over Finite Field in z6 of size 7^6 (Hess model)
        """
        if not isinstance(other, JacobianGroupFunctor):
            return None
        if not self.domain() == other.domain():
            return None
        K = pushout(self._field, other._field)
        return JacobianGroupFunctor(self.domain().base(), K)


class JacobianGroup_base(Parent):
    """
    Groups of rational points of Jacobians.

    INPUT:

    - ``parent`` -- a Jacobian

    - ``function_field`` -- a function field

    - ``base_div`` -- an effective divisor of the function field

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: J = C.jacobian(model='hess')
        sage: J.group()
        Group of rational points of Jacobian over Finite Field of size 7 (Hess model)
    """
    _embedding_map_class = None

    def __init__(self, parent, function_field, base_div):
        """
        Initialize.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: G = J.group()
            sage: TestSuite(G).run(skip=['_test_elements', '_test_pickling'])
        """
        super().__init__(base=IntegerRing(), category=CommutativeAdditiveGroups())

        self._parent = parent
        self._function_field = function_field
        self._genus = parent._function_field.genus()  # equals function_field.genus()
        self._base_div = base_div

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: J.group()
            Group of rational points of Jacobian over Finite Field of size 7 (Hess model)
        """
        F = self._function_field
        k = F.constant_base_field()
        return f'Group of rational points of Jacobian over {k}'

    def _coerce_map_from_(self, S):
        """
        Return the coerce map from ``S`` if ``S`` is embedded to ``self``.

        EXAMPLES::

            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G1 = J.group()
            sage: K = k.extension(3)
            sage: G3 = J.group(K)
            sage: G3.has_coerce_map_from(G1)
            True
        """
        if isinstance(S, JacobianGroup_base) and S.parent() is self.parent():
            K = self._function_field.constant_base_field()
            k = S._function_field.constant_base_field()
            if K.has_coerce_map_from(k):
                return self._embedding_map_class(S, self)
        return None

    def construction(self):
        """
        Return the data for a functorial construction of this Jacobian group.

        EXAMPLES::

            sage: k = GF(7)
            sage: P2.<x,y,z> = ProjectiveSpace(k, 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: K2 = k.extension(2)
            sage: G2 = J.group(K2)
            sage: K3= k.extension(3)
            sage: G3 = J.group(K3)
            sage: p1, p2 = G2.get_points(2)
            sage: q1, q2 = G3.get_points(2)
            sage: (p1 + q1).parent() is (p2 + q2).parent()
            True
        """
        k = self._parent._function_field.constant_base_field()
        K = self._function_field.constant_base_field()
        return (JacobianGroupFunctor(k, K), self._parent)

    def parent(self):
        """
        Return the Jacobian to which this Jacobian group belongs.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: G = J.group()
            sage: G.parent()
            Jacobian of Projective Plane Curve over Finite Field of size 7
             defined by x^3 - y^2*z - 2*z^3 (Hess model)
        """
        return self._parent

    def function_field(self):
        """
        Return the function field to which this Jacobian group attached.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='hess')
            sage: G = J.group()
            sage: G.function_field()
            Function field in z defined by z^3 + 4*y^2*z + 3
        """
        return self._function_field

    def base_divisor(self):
        """
        Return the base divisor that is used to represent points of this group.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: G.base_divisor()
            Place (1/y, 1/y*z)
            sage: _ == 1*b
            True

        The base divisor is the denominator (negative part) of the divisor of
        degree zero that represents a point. ::

            sage: p = C([-1,2,1]).place()
            sage: G.point(p - b).divisor()
            - Place (1/y, 1/y*z)
             + Place (y + 2, z + 1)
        """
        return self._base_div


class JacobianGroup_finite_field_base(JacobianGroup_base):
    """
    Jacobian groups of function fields over finite fields.

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
        sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: J = C.jacobian(model='hess', base_div=b)
        sage: J.group()
        Group of rational points of Jacobian over Finite Field of size 7 (Hess model)
    """
    def _bound_on_order(self):
        """
        Return an upper bound on the order of the abelian group.

        This bound depends on the genus and the order of the constant field
        of the function field. This simple bound is from [Hes2004]_.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: G._bound_on_order()
            23
        """
        F = self._function_field
        q = F.constant_base_field().order()
        g = self._genus

        c = 2*g/(q.sqrt() - 1)
        return integer_floor(math.exp(c)*q**g)

    def order(self, algorithm='numeric'):
        """
        Return the order of the Jacobian group.

        INPUT:

        - ``algorithm`` -- ``'numeric'`` (default) or ``'algebraic'``

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: J = C.jacobian(model='hess', base_div=b)
            sage: G = J.group()
            sage: G.order()
            7
        """
        F = self._parent._function_field
        g = F.genus()
        b = self._function_field.constant_base_field().degree() // F.constant_base_field().degree()

        f = F.L_polynomial()

        if algorithm == 'numeric':
            # numeric method - fast but might be inaccurate by numerical noise
            from sage.rings.qqbar import AlgebraicField
            h = Integer(math.prod([(1-a**(-b))**m for a, m in f.change_ring(AlgebraicField()).roots()]))
            return h

        # algebraic method - slow

        es = []
        s = -1
        for i in range(1, 2*g + 1):
            es.append(s*f[i])
            s = -s
        es

        ps = [es[0]]
        for i in range(1, 2*g):
            p = 0
            s = 1
            for j in range(i):
                p = p + s*es[j]*ps[-j-1]
                s = -s
            ps.append(p + s*(i + 1)*es[i])

        while len(ps) < b*2*g:
            p = 0
            s = 1
            for j in range(2*g):
                p = p + s*es[j]*ps[-j-1]
                s = -s
            ps.append(p)

        qs = [ps[b*(i + 1) - 1] for i in range(2*g)]

        fs = [qs[0]]
        for i in range(1, 2*g):
            k = qs[i]
            s = -1
            for j in range(i):
                k = k + s*fs[j]*qs[i - j - 1]
                s = -s
            fs.append(-s*k // (i + 1))

        bs = [1]
        s = -1
        for i in range(2*g):
            bs.append(s*fs[i])
            s = -s

        return sum(bs)

    def get_points(self, n):
        """
        Return `n` points of the Jacobian group.

        If `n` is greater than the order of the group, then returns
        all points of the group.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='hess')
            sage: G = J.group()
            sage: pts = G.get_points(G.order())
            sage: len(pts)
            11
        """
        lst = []
        S = iter(self)
        try:
            for i in range(n):
                lst.append(next(S))
        except StopIteration:
            pass

        return lst


class Jacobian_base(Parent):
    """
    Jacobians of function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
        sage: F.jacobian()
        Jacobian of Function field in y defined by y^2 + y + (x^2 + 1)/x (Hess model)
    """
    def __init__(self, function_field, base_div, **kwds):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: TestSuite(J).run()
        """
        self._function_field = function_field
        self._base_div = base_div
        self._system = {}
        self._base_place = None
        self._curve = kwds.get('curve')
        super().__init__(category=Jacobians(function_field.constant_base_field()),
                         base=function_field.constant_base_field(),
                         facade=True)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: F.jacobian()
            Jacobian of Function field in y defined by y^2 + y + (x^2 + 1)/x (Hess model)
        """
        return f'Jacobian of {self.base_curve()}'

    def _an_element_(self):
        """
        Return an element of ``self``.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: J.an_element()
            [Place (1/x, 1/x*y)]
        """
        return next(iter(self.group()))

    def __call__(self, x):
        """
        Return the point of ``self`` constructed from ``x``.

        It is assumed that ``self`` and ``x`` are points of the Jacobians
        attached to the same function field.

        TESTS::

            sage: # long time
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J_hess = F.jacobian(model='hess')
            sage: G = J_hess.group()
            sage: p = G.get_points(3)[2]
            sage: Jkm = F.jacobian(model='km_large')
            sage: q = Jkm(p)
            sage: p.order() == q.order()
            True
            sage: J_hess(q) == p
            True

        If ``x`` is an effective divisor, it is checked that the degree
        is equal to the degree of the base divisor. See :issue:38623.

            sage: K.<x> = FunctionField(GF(7))
            sage: _.<t> = K[]
            sage: F.<y> = K.extension(t^2 - x^6 - 3)
            sage: O = F.maximal_order()
            sage: D1 = (O.ideal(x + 1, y + 2)*O.ideal(x + 2, y + 2)).divisor()
            sage: I = O.ideal(x + 3, y+5)*O.ideal(x + 4, y + 5)*O.ideal(x + 5, y + 5)
            sage: D2 = I.divisor()
            sage: J = F.jacobian(model='hess')
            sage: J(D1)
            [Place (x + 1, y + 2)
             + Place (x + 2, y + 2)]
            sage: J(D2)
            Traceback (most recent call last):
            ...
            ValueError: effective divisor is not of degree 2
            sage: J.base_divisor().degree()
            2
        """
        F = self._function_field
        if isinstance(x, JacobianPoint_base):
            Gx = x.parent()
            Jx = Gx.parent()
            if Jx._function_field is F:
                k = Gx._function_field.constant_base_field()
                G = self.group(k)
                K = G._function_field
                return G.point(K.divisor_group()(x.divisor()))
        if x in F.place_set():
            return self(x - x.degree()*self._base_place)
        if x == 0:
            return self.group().zero()
        if x in F.divisor_group():
            return self.group()(x)
        raise ValueError(f"cannot create a point of the Jacobian from {x}")

    def curve(self):
        """
        Return the projective curve to which this Jacobian is attached.

        If the Jacobian was constructed from a function field, then returns nothing.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: J.curve()
        """
        return self._curve

    def base_curve(self):
        """
        Return the base curve or the function field that abstractly defines a curve.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: J.base_curve()
            Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        return self._function_field if self._curve is None else self._curve

    def facade_for(self):
        """
        Return the system of groups that this Jacobian is a facade for.

        The Jacobian can be seen as a facade for all groups of rational points
        over field extensions of the base constant field of the function field.
        This method returns only the internally constructed system of such
        groups.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: J.facade_for()
            [Group of rational points of Jacobian over Finite Field of size 2 (Hess model)]
        """
        if not self._system:
            return [self.group()]
        return list(self.group(k) for k in self._system)

    def base_divisor(self):
        """
        Return the base divisor used to construct the Jacobian.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: b = F.get_place(1)
            sage: J = F.jacobian(base_div=b)
            sage: J.base_divisor() == b
            True
        """
        return self._base_div

    def group(self, k_ext=None):
        """
        Return the group of rational points of the Jacobian.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: b = F.get_place(1)
            sage: J = F.jacobian(base_div=b)
            sage: J.group()
            Group of rational points of Jacobian over Finite Field of size 2 (Hess model)
        """
        F = self._function_field
        k = F.constant_base_field()

        if k_ext in self._system:
            return self._system[k_ext][0]

        if k_ext is None or k_ext is k:
            ext = F.extension_constant_field(k)
            grp = self._group_class(self, F, self._base_div)
            if self._base_place is not None:
                grp._base_place = self._base_place
            self._system[k] = (grp, ext)
        else:
            ext = F.extension_constant_field(k_ext)
            base_div = ext.conorm_divisor(self._base_div)
            grp = self._group_class(self, ext.top(), base_div)
            if self._base_place is not None:
                grp._base_place = ext.conorm_place(self._base_place)
            self._system[k_ext] = (grp, ext)

        return grp

    def set_base_place(self, place):
        """
        Set ``place`` as the base place.

        INPUT:

        - ``place`` -- a rational place of the function field

        The base place `B` is used to map a rational place `P` of the function
        field to the point of the Jacobian defined by the divisor `P - B`.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: J = F.jacobian()
            sage: B = F.get_place(1)
            sage: J.set_base_place(B)
            sage: Q = F.places()[-1]
            sage: J(Q)
            [Place (x + 1, x*y + 1)]
            sage: J(Q).parent()
            Group of rational points of Jacobian over Finite Field of size 2 (Hess model)
            sage: J(B)
            [Place (x, x*y)]
            sage: J(B).is_zero()
            True
        """
        self._base_place = place

        for k in self._system:
            grp, ext = self._system[k]
            grp._base_place = ext.conorm_place(place)
