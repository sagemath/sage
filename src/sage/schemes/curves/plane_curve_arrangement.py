# -*- coding: utf-8 -*-
r"""
Affine Plane Curve Arrangements

We create an :class:`OrderedAffinePlaneCurveArrangements`
object following the properties of :class:`HyperplaneArrangements`

    sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
    sage: C = H(3*x + 2*y - x^2 + y^3 - 7);  C
    Arrangement (Affine Plane Curve over Rational Field defined by y^3 - x^2 + 3*x + 2*y - 7)

The individual curves will be in  :class:`AffinePlaneCurve`::

    sage: C[0].parent()
    <class 'sage.schemes.curves.affine_curve.IntegralAffinePlaneCurve_with_category'>

The default base field is `\QQ`, the rational numbers.
Number fields are also possible (also with fixed embeddings in ``QQbar``)::

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ, 'x')
    sage: NF.<a> = NumberField(x^4 - 5 * x^2 + 5, embedding=1.90)
    sage: H.<y,z> = OrderedAffinePlaneCurveArrangements(NF)
    sage: A = H(y^2 - a * z, y^2 + a * z); A
    Arrangement (Affine Plane Curve over Number Field in a with defining polynomial
                 x^4 - 5*x^2 + 5 with a = 1.902113032590308? defined by y^2 + (-a)*z,
                 Affine Plane Curve over Number Field in a with defining polynomial
                 x^4 - 5*x^2 + 5 with a = 1.902113032590308? defined by y^2 + a*z)
    sage: A.base_ring()
    Number Field in a with defining polynomial x^4 - 5*x^2 + 5
    with a = 1.902113032590308?


AUTHORS:

- Enrique Artal (2023-10): initial version
"""

# *****************************************************************************
#       Copyright (C) 2023 Enrique Artal <artal@unizar.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.sets_cat import Sets
from sage.combinat.combination import Combinations
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.qqbar import QQbar
from sage.rings.ring import _Fields
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.curves.affine_curve import AffinePlaneCurve
from sage.schemes.curves.constructor import Curve
from sage.schemes.curves.zariski_vankampen import braid_monodromy
from sage.schemes.curves.zariski_vankampen import fundamental_group_arrangement
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


class OrderedAffinePlaneCurveArrangementsElement(Element):
    """
    An ordered affine plane curve arrangement.
    #
    # .. WARNING::
    #
    #     You should never create
    #     :class:`HyperplaneArrangementElement` instances directly,
    #     always use the parent.
    """

    def __init__(self, parent, curves, check=True):
        """
        Construct an ordered affine plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`HyperplaneArrangements`

        - ``curves`` -- a tuple of curves

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement (Affine Plane Curve over Rational Field defined by x,
                         Affine Plane Curve over Rational Field defined by y)
            sage: TestSuite(elt).run()
        """
        super().__init__(parent)
        self._curves = curves
        if check:
            if not isinstance(curves, tuple):
                raise ValueError("the curves must be given as a tuple")
            if not all(isinstance(h, AffinePlaneCurve) for h in curves):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is self.parent().ambient_space() for h in curves):
                raise ValueError("not all curves are in the same ambient space")
        self._braid_monodromy = None
        self._vertical_braid_monodromy = None
        self._strands = dict()
        self._vertical_strands = dict()
        self._fundamental_group = None
        self._meridians = dict()
        self._infinity = None

    def __getitem__(self, i):
        """
        Return the `i`-th curve.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        The `i`-th curve.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H(y^2 - x, y^3 + 2 * x^2, x^4 + y^4 + 1); h
            Arrangement (Affine Plane Curve over Rational Field defined by y^2 - x,
                         Affine Plane Curve over Rational Field defined by y^3 + 2*x^2,
                         Affine Plane Curve over Rational Field defined by x^4 + y^4 + 1)
        """
        return self._curves[i]

    def __hash__(self):
        r"""
        TESTS::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H((x * y, x + y +1))
            sage: len_dict = {h: len(h)}
        """
        return hash(self.curves())

    def n_curves(self):
        r"""
        Return the number of curves in the arrangement.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H((x * y, x + y +1))
            sage: h.n_curves()
            2
            sage: len(h)    # equivalent
            2
        """
        return len(self._curves)

    __len__ = n_curves

    def curves(self):
        r"""
        Return the curves in the arrangement as a tuple.

        OUTPUT:

        A tuple.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H((x * y, x + y + 1))
            sage: h.curves()
            (Affine Plane Curve over Rational Field defined by x*y,
             Affine Plane Curve over Rational Field defined by x + y + 1)

        Note that the hyperplanes can be indexed as if they were a list::

            sage: h[1]
            Affine Plane Curve over Rational Field defined by x + y + 1
        """
        return self._curves

    def _repr_(self):
        r"""
        String representation for a curve arrangement.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: h
            Arrangement of 5 curves
            sage: H(())
            Empty curve arrangement
        """
        if len(self) == 0:
            return 'Empty curve arrangement'
        elif len(self) < 5:
            curves = ', '.join(h._repr_() for h in self._curves)
            return 'Arrangement ({0})'.format(curves)
        return 'Arrangement of {0} curves'.format(len(self))

    def _richcmp_(self, other, op):
        """
        Compare two curve arrangements.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: H(x) == H(y)
            False
            sage: H(x) == H(2 * x)
            True

        TESTS::

            sage: H(x) == 0
            False
        """
        return richcmp(self._curves, other._curves, op)

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a curve arrangement or something that can
          be converted into a curve arrangement

        OUTPUT:

        A new curve arrangement.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: C = Curve(x^8 - y^8 -x^4 * y^4)
            sage: h1 = h.union(C); h1
            Arrangement of 6 curves
            sage: h1 == h1.union(C)
            Repeated curve
            True
        """
        P = self.parent()
        other_h = P(other)
        curves0 = self._curves + other_h._curves
        curves = ()
        for h in curves0:
            if h not in curves:
                curves += (h, )
            else:
                print("Repeated curve")
        result = P(*curves)
        return result

    add_curve = union

    __or__ = union

    def deletion(self, curves):
        r"""
        Return the curve arrangement obtained by removing ``h``.

        INPUT:

        - ``h`` -- a curve or curve arrangement

        OUTPUT:

        A new curve arrangement with the given curve(s)
        ``h`` removed.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: C = h[-1]
            sage: h.deletion(C)
            Arrangement (Affine Plane Curve over Rational Field defined by x*y,
                         Affine Plane Curve over Rational Field defined by x + y + 1,
                         Affine Plane Curve over Rational Field defined by -y^5 + x^3,
                         Affine Plane Curve over Rational Field defined by x^5 + y^5 + x^2*y^2)
            """
        parent = self.parent()
        curves = parent(curves)
        planes = list(self)
        for curve in curves:
            try:
                planes.remove(curve)
            except ValueError:
                raise ValueError('curve is not in the arrangement')
        return parent(planes)

    def change_ring(self, base_ring):
        """
        Return curve arrangement over the new base ring.

        INPUT:

        - ``base_ring`` -- the new base ring; must be a field for
          curve arrangements

        OUTPUT:

        The curve arrangement obtained by changing the base
        field, as a new curve arrangement.

        EXAMPLES::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: K.<a> = CyclotomicField(3)  # optional - sage.rings.number_field
            sage: A.change_ring(K)            # optional - sage.rings.number_field
            Arrangement (Affine Plane Curve over Cyclotomic Field of order 3 and degree 2 defined by -x^3 + y^2,
                         Affine Plane Curve over Cyclotomic Field of order 3 and degree 2 defined by x,
                         Affine Plane Curve over Cyclotomic Field of order 3 and degree 2 defined by y,
                         Affine Plane Curve over Cyclotomic Field of order 3 and degree 2 defined by x^2 + x*y + y^2)
        """
        parent = self.parent().change_ring(base_ring)
        curves = tuple(c.change_ring(base_ring) for c in self)
        return parent(curves)

    def defining_polynomials(self):
        r"""
        Return the defining polynomials of the elements of``self``.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: A.defining_polynomials()
            (-x^3 + y^2, x, y, x^2 + x*y + y^2)
        """
        return tuple(h.defining_polynomial() for h in self)

    def defining_polynomial(self, simplified=True):
        r"""
        Return the defining polynomial of the union of the curves in ``self``.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y ** 2 + x ** 2, x, y)
            sage: prod(A.defining_polynomials()) == A.defining_polynomial()
            True
        """
        return prod(self.defining_polynomials())

    def have_common_factors(self):
        r"""
        Check if th curves have common factors.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(x * y, x^2 + x* y^3)
            sage: A.have_common_factors()
            True
        """
        L = [c.defining_polynomial() for c in self]
        C = Combinations(L, 2)
        for f1, f2 in C:
            if f1.gcd(f2).degree() > 0:
                return True
        return False

    def reduce(self):
        r"""
        Replace the curves by their reduction.

        EXAMPLES::

            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2, (x + y)^3 * (x^2 + x * y + y^2))
            sage: A.reduce()
            Arrangement (Affine Plane Curve over Rational Field defined by y,
            Affine Plane Curve over Rational Field defined by x^3 + 2*x^2*y + 2*x*y^2 + y^3)
        """
        P = self.parent()
        L = [c.defining_polynomial().radical() for c in self]
        return P(*L)

    def fundamental_group(self, simplified=True, vertical=False):
        r"""
        It computes the fundamental group of the complement of the union
        of affine plane curves in `\mathbb{C}^2`.

        INPUT:

        - ``vertical`` -- boolean (default: False). If it is ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified braid braid_monodromy is used.

        OUTPUT:

        A group finitely presented with the assignation of each curve to
        a set of meridians, including the line at infinity.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.fundamental_group()
            (Finitely presented group < x0, x1, x2 | x2^-1*x1^-1*x2*x1,
                                                     x1*x0*x1^-1*x0^-1,
                                                     (x0*x2)^2*(x0^-1*x2^-1)^2 >,
             {0: [x2, x0*x2*x0^-1], 1: [x1], 2: [x0],
              3: [x0*x2^-1*x0^-1*x2^-1*x1^-1*x0^-1]})
            sage: A.fundamental_group(vertical=True)
            (Finitely presented group < x0, x1, x2 | x1*x0^-1*x1^-1*x0,
                                                     x2*x1*x2^-1*x1^-1, x2*x0*x2^-1*x0^-1 >,
             {0: [x1], 1: [x0], 2: [x2], 3: [x2^-1*x1^-2*x0^-1]})

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        C = self.reduce()
        L = C.defining_polynomials()
        G, dic = fundamental_group_arrangement(L, simplified=simplified, puiseux=True, vertical=vertical)
        return (G, dic)

    def braid_monodromy(self, vertical=False):
        r"""
        It computes the braid monodromy of the complement of the union
        of affine plane curves in `\mathbb{C}^2`. If there are vertical
        asymptotes a change of variable is done.

        INPUT:

        - ``vertical`` -- boolean (default: False). If it is ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified braid braid_monodromy is computed.

        OUTPUT:

        A braid monodromy with dictionnaries identifying strans with components
        and braids with vertical lines..

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.braid_monodromy()
            ([s1*s0*(s1*s2*s1)^2*s2*(s1^-1*s2^-1)^2*s1^-1*s0^-1*s1^-1,
              s1*s0*(s1*s2)^2*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
              s1*s0*s1*s2*(s1*s2^-1)^2*s0*s1*s2*s1*s0*s2^-1*s1^-3*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
              s1*s0*s1*s2*s1*s2^-1*s1^4*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
              s1*s0*s1*s2*s1*s2^-1*s1^-1*s2*s0^-1*s1^-1],
              {0: 2, 1: 1, 2: 0, 3: 0}, {})
            sage: A.braid_monodromy(vertical=True)
            ([s1*s0*s1*s0^-1*s1^-1*s0, s0^-1*s1*s0*s1^-1*s0, s0^-1*s1^2*s0],
             {0: 1, 1: 0, 2: 0}, {1: 2})

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if self._braid_monodromy:
            return self._braid_monodromy
        L = self.defining_polynomials()
        return braid_monodromy(prod(L), arrangement=L, vertical=vertical)[:-1]


class OrderedAffinePlaneCurveArrangements(Parent, UniqueRepresentation):
    """
    Curve arrangements.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
        sage: H(x, y^2, x-1, y-1)
        Arrangement (Affine Plane Curve over Rational Field defined by x,
                     Affine Plane Curve over Rational Field defined by y^2,
                     Affine Plane Curve over Rational Field defined by x - 1,
                     Affine Plane Curve over Rational Field defined by y - 1)
    """
    Element = OrderedAffinePlaneCurveArrangementsElement

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: K = OrderedAffinePlaneCurveArrangements(QQ, names=('x', 'y'))
            sage: H is K
            True
            sage: type(K)
            <class 'sage.schemes.curves.plane_curve_arrangement.OrderedAffinePlaneCurveArrangements_with_category'>

        TESTS::

            sage: H.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: TestSuite(H).run()
        """
        if base_ring not in _Fields:
            raise ValueError('base ring must be a field')
        super().__init__(category=Sets())
        self._base_ring = base_ring
        self._names = names

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        The base ring of the curve arrangement.

        EXAMPLES::

            sage: L.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: L.base_ring()
            Rational Field
        """
        return self._base_ring

    def change_ring(self, base_ring):
        """
        Return curve arrangements over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring.

        OUTPUT:

        A new :class:`OrderedAffinePlaneCurveArrangements` instance over the new
        base ring.

        EXAMPLES::

            sage: L.<x,y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: L.gen(0)
            x
            sage: L.change_ring(RR).base_ring()
            Real Field with 53 bits of precision

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return OrderedAffinePlaneCurveArrangements(base_ring, names=self._names)

    @cached_method
    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Affine Space of dimension 2 over Rational Field
        """
        return AffineSpace(self.base_ring(), 2, self._names)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ);  L
            Curve arrangements in Affine Space of dimension 2 over Rational Field
        """
        return 'Curve arrangements in {0}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a curve

        EXAMPLES::

            sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: A = L._element_constructor_(x, y); A
            Arrangement (Affine Plane Curve over Rational Field defined by x,
            Affine Plane Curve over Rational Field defined by y)
            sage: L._element_constructor_([x, y]) == A
            True
            sage: L._element_constructor_(Curve(x), Curve(y)) == A
            True
            sage: L._element_constructor_(y, x) == A
            False
       """
        if len(args) == 1 and not (isinstance(args[0], (tuple, list))):
            arg = (args[0], )
        elif len(args) == 1:
            arg = tuple(args[0])
        else:
            arg = tuple(args)
        # process keyword arguments
        if len(kwds) > 0:
            raise ValueError('unknown keyword argument')
        # process positional arguments
        AA = self.ambient_space()
        R = AA.coordinate_ring()
        curves = ()
        for h in arg:
            try:
                ambient = h.ambient_space()
                if ambient == AA:
                    curves += (h, )
                else:
                    return None
            except AttributeError:
                try:
                    h = R(h)
                    curves += (Curve(h), )
                except TypeError:
                    return None
        return self.element_class(self, curves)

    def _an_element_(self):
        """
        Dirty trick to avoid test run failure.

        TESTS::

            sage: H.<t, s> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: H._an_element_()
            Arrangement (Affine Plane Curve over Rational Field defined by t)
        """
        x = self.gen(0)
        return self(x)

    @cached_method
    def ngens(self):
        """
        Return the number of variables, i.e. 2, kept for completness.

        OUTPUT:

        An integer (2).

        EXAMPLES::

            sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: L.ngens()
            2
        """
        return len(self._names)

    @cached_method
    def gens(self):
        """
        Return the coordinates.

        OUTPUT:

        A tuple of linear expressions, one for each linear variable.

        EXAMPLES::

            sage: L = OrderedAffinePlaneCurveArrangements(QQ, ('x', 'y'))
            sage: L.gens()
            (x, y)
        """
        return self.ambient_space().gens()

    def gen(self, i):
        """
        Return the `i`-th coordinate.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        A variable.

        EXAMPLES::

            sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ)
            sage: L.gen(1)
            y
        """
        return self.gens()[i]

    # def _coerce_map_from_(self, P):
    #     """
    #     Return whether there is a coercion.
    #
    #     TESTS::
    #
    #         sage: L.<x, y> = OrderedAffinePlaneCurveArrangements(QQ);  L
    #         Curve arrangements in Affine Space of dimension 2 over Rational Field
    #         sage: M.<x, y> = HyperplaneArrangements(RR);  M
    #         Hyperplane arrangements in 2-dimensional linear space over Real Field with 53 bits of precision with coordinates x, y
    #
    #         sage: L.coerce_map_from(ZZ)
    #         Coercion map:
    #           From: Integer Ring
    #           To:   Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
    #         sage: M.coerce_map_from(L)
    #         Coercion map:
    #           From: Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
    #           To:   Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y
    #         sage: L.coerce_map_from(M)
    #     """
    #     if self.ambient_space().has_coerce_map_from(P):
    #         return True
    #     if isinstance(P, OrderedAffinePlaneCurveArrangements):
    #         return self.base_ring().has_coerce_map_from(P.base_ring())
    #     return super()._coerce_map_from_(P)
