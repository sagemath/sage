# -*- coding: utf-8 -*-
r"""
Affine and Projective Plane Curve Arrangements

We create classes :class:`AffinePlaneCurveArrangements`
and :class:`ProjectivePlaneCurveArrangements`
following the properties of :class:`HyperplaneArrangements`::

    sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
    sage: C = H(3*x + 2*y - x^2 + y^3 - 7);  C
    Arrangement (y^3 - x^2 + 3*x + 2*y - 7) in Affine Space of dimension 2 over Rational Field

The individual curves will be in  :class:`AffinePlaneCurve` or in :class:`ProjectivePlaneCurve`::

    sage: C[0].parent()
    <class 'sage.schemes.curves.affine_curve.IntegralAffinePlaneCurve_with_category'>

The default base field is `\mathbb{Q}`, the rational numbers.
Number fields are also possible (also with fixed embeddings in
`\overline{\mathbb{Q}}`)::

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ, 'x')
    sage: NF.<a> = NumberField(x^4 - 5 * x^2 + 5, embedding=1.90)
    sage: H.<y,z> = AffinePlaneCurveArrangements(NF)
    sage: A = H(y^2 - a * z, y^2 + a * z); A
    Arrangement (y^2 + (-a)*z, y^2 + a*z) in Affine Space of dimension 2
    over Number Field in a with defining polynomial
    x^4 - 5*x^2 + 5 with a = 1.902113032590308?
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
from sage.groups.free_group import FreeGroup
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar
from sage.rings.ring import _Fields
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.curves.affine_curve import AffinePlaneCurve
from sage.schemes.curves.constructor import Curve
from sage.schemes.curves.projective_curve import ProjectiveSpace
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve
from sage.schemes.curves.zariski_vankampen import braid_monodromy
from sage.schemes.curves.zariski_vankampen import fundamental_group_arrangement
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


class AffinePlaneCurveArrangementsElement(Element):
    """
    An ordered affine plane curve arrangement.
    """

    def __init__(self, parent, curves, check=True):
        """
        Construct an ordered affine plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`AffinePlaneCurveArrangements`

        - ``curves`` -- a tuple of curves

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement (x, y) in Affine Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
        """
        super().__init__(parent)
        self._curves = curves
        if check:
            if not isinstance(curves, tuple):
                raise ValueError("the curves must be given as a tuple")
            if not all(isinstance(h, AffinePlaneCurve) for h in curves):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is self.parent().ambient_space()
                       for h in curves):
                raise ValueError("not all curves are in the same ambient space")
        self._braid_monodromy = None
        self._braid_monodromy_with_vertical = None
        self._strands = dict()
        self._strands_with_vertical = dict()
        self._vertical_lines_in_braid_mon = None
        self._fundamental_group_vertical_simplified = None
        self._meridians_vertical_simplified = dict()
        self._fundamental_group_vertical = None
        self._meridians_vertical = dict()
        self._fundamental_group_simplified = None
        self._meridians_simplified = dict()
        self._fundamental_group = None
        self._meridians = dict()

    def __getitem__(self, i):
        """
        Return the `i`-th curve.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        The `i`-th curve.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H(y^2 - x, y^3 + 2 * x^2, x^4 + y^4 + 1); h
            Arrangement (y^2 - x, y^3 + 2*x^2, x^4 + y^4 + 1)
            in Affine Space of dimension 2 over Rational Field
        """
        return self._curves[i]

    def __hash__(self):
        r"""
        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
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

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
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

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H((x * y, x + y + 1))
            sage: h.curves()
            (Affine Plane Curve over Rational Field defined by x*y,
             Affine Plane Curve over Rational Field defined by x + y + 1)

        Note that the curves can be indexed as if they were a list::

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

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: h
            Arrangement of 5 curves in Affine Space of dimension 2 over Rational Field
            sage: H(())
            Empty curve arrangement in Affine Space of dimension 2 over Rational Field
        """
        if len(self) == 0:
            return 'Empty curve arrangement in {0}'.format(self.parent().ambient_space())
        elif len(self) < 5:
            curves = ', '.join(h.defining_polynomial()._repr_()
                               for h in self._curves)
            return 'Arrangement ({0}) in {1}'.format(curves,
                                                     self.parent().ambient_space())
        return 'Arrangement of {0} curves in {1}'.format(len(self),
                                                         self.parent().ambient_space())

    def _richcmp_(self, other, op):
        """
        Compare two curve arrangements.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
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

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: C = Curve(x^8 - y^8 -x^4 * y^4)
            sage: h1 = h.union(C); h1
            Arrangement of 6 curves in Affine Space of dimension 2 over Rational Field
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

    add_curves = union

    __or__ = union

    def deletion(self, curves):
        r"""
        Return the curve arrangement obtained by removing ``curves``.

        INPUT:

        - ``curves`` -- a curve or curve arrangement

        OUTPUT:

        A new curve arrangement with the given curve(s)
        ``h`` removed.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: C = h[-1]
            sage: h.deletion(C)
            Arrangement (x*y, x + y + 1, -y^5 + x^3, x^5 + y^5 + x^2*y^2)
            in Affine Space of dimension 2 over Rational Field
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

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: K.<a> = CyclotomicField(3)  # optional - sage.rings.number_field
            sage: A.change_ring(K)            # optional - sage.rings.number_field
            Arrangement (-x^3 + y^2, x, y, x^2 + x*y + y^2) in Affine Space of
            dimension 2 over Cyclotomic Field of order 3 and degree 2
        """
        parent = self.parent().change_ring(base_ring)
        curves = tuple(c.change_ring(base_ring) for c in self)
        return parent(curves)

    def coordinate_ring(self):
        """
        Return the coordinate ring of ``self``.

        OUTPUT:

        The coordinate ring of the curve arrangement.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: C = L(x, y)
            sage: C.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.curves()[0].defining_polynomial().parent()

    def defining_polynomials(self):
        r"""
        Return the defining polynomials of the elements of ``self``.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: A.defining_polynomials()
            (-x^3 + y^2, x, y, x^2 + x*y + y^2)
        """
        return tuple(h.defining_polynomial() for h in self)

    def defining_polynomial(self, simplified=True):
        r"""
        Return the defining polynomial of the union of the curves in ``self``.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y ** 2 + x ** 2, x, y)
            sage: prod(A.defining_polynomials()) == A.defining_polynomial()
            True
        """
        return prod(self.defining_polynomials())

    def have_common_factors(self):
        r"""
        Check if the curves have common factors.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
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

    def reduce(self, clean=False):
        r"""
        Replace the curves by their reduction.

        INPUT:

        - ``clean`` -- boolean (default: False); if ``False``
          and there are common factors it returns ``None`` and
          a warning message. If ``True``, the common factors are kept
          only in the first occurance.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2, (x + y)^3 * (x^2 + x * y + y^2))
            sage: A.reduce()
            Arrangement (y, x^3 + 2*x^2*y + 2*x*y^2 + y^3) in Affine Space
            of dimension 2 over Rational Field
            sage: C = H(x*y, x*(y + 1))
            sage: C.reduce()
            Some curves have common components
            sage: C.reduce(clean=True)
            Arrangement (x*y, y + 1) in Affine Space of dimension 2
            over Rational Field
            sage: C = H(x*y, x)
            sage: C.reduce(clean=True)
            Arrangement (x*y) in Affine Space of dimension 2 over Rational Field
        """
        P = self.parent()
        R = self.coordinate_ring()
        L = [self[0].defining_polynomial().radical()]
        for c in self[1:]:
            g = c.defining_polynomial().radical()
            for f in L:
                d = g.gcd(f)
                if d.degree() > 0 and not clean:
                    print("Some curves have common components")
                    return None
                g = R(g / d)
            if g.degree() > 0:
                L.append(g)
        return P(*L)

    def fundamental_group(self, simplified=True, vertical=True,
                          projective=False):
        r"""
        The fundamental group of the complement of the union
        of affine plane curves in `\mathbb{C}^2`.

        INPUT:

        - ``vertical`` -- boolean (default: True); if it is ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified braid braid_monodromy is used.

        - ``simplified`` -- boolean (default: True); if it is ``True``, the
          group is simplified.

        - ``projective`` -- boolean (default: False); to be used in the
          method for projective curves.

        OUTPUT:

        A finitely presented group.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.fundamental_group()
            Finitely presented group
            < x0, x1, x2 | x2*x0*x2^-1*x0^-1, x1*x0*x1^-1*x0^-1, (x2*x1)^2*(x2^-1*x1^-1)^2 >
            sage: A.meridians()
            {0: [x1, x2*x1*x2^-1], 1: [x0], 2: [x2],
             3: [x1^-1*x2^-1*x1^-1*x0^-1]}
            sage: A.fundamental_group(simplified=False)
            Finitely presented group
            < x0, x1, x2, x3 | x1*x0*x1^-1*x0^-1, x2*x0*x2^-1*x0^-1,
                               x3*x0*x1*x0^-1*x2^-1*x0*x2*x0*x1^-1*x0^-1*x3^-1*x0^-1,
                               x3*x0*x1*x0^-1*x2^-1*x0^-1*(x2*x0)^2*x1^-1*x0^-1*x3^-1*x1^-1,
                               x3*x0*x1*x0^-1*x3^-1*x2^-1 >
            sage: A.meridians(simplified=False)
            {0: [x1, x2], 1: [x0], 2: [x3], 3: [x3^-1*x2^-1*x1^-1*x0^-1]}
            sage: A.fundamental_group(vertical=False)
            Finitely presented group
            < x0, x1, x2 | x2^-1*x1^-1*x2*x1, x1*x0*x1^-1*x0^-1, (x0*x2)^2*(x0^-1*x2^-1)^2 >
            sage: A.meridians(vertical=False)
            {0: [x2, x0*x2*x0^-1], 1: [x1], 2: [x0], 3: [x0*x2^-1*x0^-1*x2^-1*x1^-1*x0^-1]}
            sage: A.fundamental_group(simplified=False, vertical=False)
            Finitely presented group
            < x0, x1, x2, x3 | x3*x2^-1*x1*x2*x3^-1*x2^-1*x1^-1*x2,
                              x1*x0*x1^-1*x0^-1,
                              x3^-1*x2^-1*x0^-1*x2*x3*x2^-1*x0*x2*x3*x2^-1,
                              x3^-1*(x2^-1*x0*x2*x3)^2*x2^-1*x0^-1*x2*x3^-1*x2^-1*x0^-1*x2,
                              x3*x2^-1*x1*x2*x3^-1*x2^-1*x1^-1*x2 >
            sage: A.meridians(simplified=False, vertical=False)
            {0: [x2, x3], 1: [x1], 2: [x0], 3: [x3^-1*x2^-1*x1^-1*x0^-1]}
            sage: A = H(x * y^2 + x + y, y + x -1, x, y)
            sage: A.fundamental_group()
            Finitely presented group
            < x0, x1, x2, x3 | x3*x0^-1*x3^-1*x0, x3*x1*x3^-1*x1^-1,
                               x3*x2*x3^-1*x2^-1, x2*x0*x2^-1*x0^-1,
                               x1*x0^-1*x1^-1*x0, x1*x2*x1^-1*x2^-1 >

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if self._fundamental_group and not vertical and not simplified:
            return self._fundamental_group
        if self._fundamental_group_simplified and simplified and not vertical:
            return self._fundamental_group_simplified
        if self._fundamental_group_vertical and not simplified and vertical:
            return self._fundamental_group_vertical
        if self._fundamental_group_vertical_simplified and simplified \
           and vertical:
            return self._fundamental_group_vertical_simplified
        K = self.base_ring()
        R = self.coordinate_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        C = self.reduce()
        L = C.defining_polynomials()
        if not vertical and self._braid_monodromy is not None:
            d1 = prod(L).degree()
            bd = (self._braid_monodromy, self._strands, dict(), d1)
        elif vertical and self._braid_monodromy_with_vertical is not None:
            d1 = prod(L).degree(R.gen(1))
            bd = (self._braid_monodromy_with_vertical, self._strands_with_vertical,
                  self._vertical_lines_in_braid_mon, d1)
        else:
            bd = None
        G, dic = fundamental_group_arrangement(L, simplified=simplified,
                                               puiseux=True,
                                               projective=projective,
                                               vertical=vertical,
                                               braid_data=bd)
        if not vertical and not simplified:
            self._fundamental_group = G
            self._meridians = dic
        elif not vertical:
            self._fundamental_group_simplified = G
            self._meridians_simplified = dic
        elif not simplified:
            self._fundamental_group_vertical = G
            self._meridians_vertical = dic
        else:
            self._fundamental_group_vertical_simplified = G
            self._meridians_vertical_simplified = dic
        return G

    def meridians(self, simplified=True, vertical=True):
        r"""
        Meridians of each irreducible component.

        OUTPUT:

        A dictionary which associates the index of each curve with its meridians,
        including the line at infinity if it can be omputed

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(x-1, y, x, y - 1)
            sage: A.fundamental_group()
            Finitely presented group
            < x0, x1, x2, x3 | x2*x0*x2^-1*x0^-1, x2*x1*x2^-1*x1^-1,
                               x3*x0*x3^-1*x0^-1, x3*x1*x3^-1*x1^-1 >
            sage: A.meridians()
            {0: [x2], 1: [x0], 2: [x3], 3: [x1], 4: [x3^-1*x2^-1*x1^-1*x0^-1]}

        This function needs
        :func:`AffinePlaneCurveArrangements.fundamental_group` with the same options,
        where some examples are shown.

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not vertical and not simplified:
            computed = bool(self._meridians)
        elif not vertical:
            computed = bool(self._meridians_simplified)
        elif not simplified:
            computed = bool(self._meridians_vertical)
        else:
            computed = bool(self._meridians_vertical_simplified)
        if not computed:
            self.fundamental_group(simplified=simplified, vertical=vertical)
        if not vertical and not simplified:
            return self._meridians
        if simplified and not vertical:
            return self._meridians_simplified
        if not simplified and vertical:
            return self._meridians_vertical
        if simplified and vertical:
            return self._meridians_vertical_simplified

    def braid_monodromy(self, vertical=True):
        r"""
        It computes the braid monodromy of the compleme+nt of the union
        of affine plane curves in `\mathbb{C}^2`. If there are vertical
        asymptotes a change of variable is done.

        INPUT:

        - ``vertical`` -- boolean (default: True). If it is ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified braid_monodromy is computed.

        OUTPUT:

        A braid monodromy with dictionnaries identifying strans with components
        and braids with vertical lines.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.braid_monodromy(vertical=False)
            [s1*s0*(s1*s2*s1)^2*s2*(s1^-1*s2^-1)^2*s1^-1*s0^-1*s1^-1,
             s1*s0*(s1*s2)^2*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
             s1*s0*s1*s2*(s1*s2^-1)^2*s0*s1*s2*s1*s0*s2^-1*s1^-3*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
             s1*s0*s1*s2*s1*s2^-1*s1^4*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
             s1*s0*s1*s2*s1*s2^-1*s1^-1*s2*s0^-1*s1^-1]
            sage: A.braid_monodromy(vertical=True)
            [s1*s0*s1*s0^-1*s1^-1*s0, s0^-1*s1*s0*s1^-1*s0, s0^-1*s1^2*s0]

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not vertical and self._braid_monodromy is not None:
            return self._braid_monodromy
        if vertical and self._braid_monodromy_with_vertical is not None:
            return self._braid_monodromy_with_vertical
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        L = self.defining_polynomials()
        bm, dic, dv, d1 = braid_monodromy(prod(L), arrangement=L,
                                          vertical=vertical)
        if vertical:
            self._braid_monodromy_with_vertical = bm
            self._strands_with_vertical = dic
            self._vertical_lines_in_braid_mon = dv
        else:
            self._braid_monodromy = bm
            self._strands = dic
        return bm

    def strands(self):
        r"""
        Strands for each member of the arrangement.

        OUTPUT:

        A dictionary which associates to the index of each strand
        its associated component if the braid monodromy has been
        calculated with ``vertical=False``.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: bm = A.braid_monodromy()
            sage: A.strands()
            {0: 2, 1: 1, 2: 0, 3: 0}

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not self._strands:
            self._braid_monodromy = self.braid_monodromy(vertical=False)
        return self._strands

    def vertical_strands(self):
        r"""
        Vertical strands for each member of the arrangement.

        OUTPUT:

        A dictionary which associates to the index of each strand
        its associated component if the braid monodromy has been
        calculated with ``vertical=True``.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: bm = A.braid_monodromy(vertical=True)
            sage: A.vertical_strands()
            {0: 1, 1: 0, 2: 0}

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not self._strands_with_vertical:
            self.braid_monodromy(vertical=True)
        return self._strands_with_vertical

    def vertical_lines_in_braid_mon(self):
        r"""
        Vertical lines in the arrangement.

        OUTPUT:

        A dictionary which associates the index of a braid
        to the index of the vertical line associated to the braid.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.braid_monodromy(vertical=True)
            [s1*s0*s1*s0^-1*s1^-1*s0, s0^-1*s1*s0*s1^-1*s0, s0^-1*s1^2*s0]
            sage: A.vertical_lines_in_braid_mon()
            {1: 2}

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not self._vertical_lines_in_braid_mon:
            self.braid_monodromy(vertical=True)
        return self._vertical_lines_in_braid_mon


class ProjectivePlaneCurveArrangementsElement(AffinePlaneCurveArrangementsElement):
    """
    An ordered projective plane curve arrangement.
    """

    def __init__(self, parent, curves, check=True):
        """
        Construct an ordered projective plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`ProjectivePlaneCurveArrangements`

        - ``curves`` -- a tuple of curves

        EXAMPLES::

            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: elt = H(x, y, z); elt
            Arrangement (x, y, z) in Projective Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
        """
        Element.__init__(self, parent)
        self._curves = curves
        if check:
            if not isinstance(curves, tuple):
                raise ValueError("the curves must be given as a tuple")
            if not all(isinstance(h, ProjectivePlaneCurve) for h in curves):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is self.parent().ambient_space()
                       for h in curves):
                raise ValueError("not all curves are in the same ambient space")
        self._fundamental_group_simplified = None
        self._meridians_simplified = dict()
        self._fundamental_group = None
        self._meridians = dict()

    def fundamental_group(self, simplified=True):
        r"""
        The fundamental group of the complement of the union
        of projective plane curves in `\mathbb{P}^2`.

        INPUT:

        - ``simplified`` -- boolean (default: True); set if the group
          is simplified..

        OUTPUT:

        A finitely presented group.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x*z, y + x - z, x)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1 | x1^-1*x0*x1*x0^-1 >
            sage: A.meridians()
            {0: [x1], 1: [x1^-1*x0^-1*x1^-1], 2: [x0]}
            sage: A.fundamental_group(simplified=False)
            Finitely presented group
            < x0, x1, x2, x3 | x3*x2^-1*x1*x2*x3^-1*x2^-1*x1^-1*x2,
                               x1*x0*x1^-1*x0^-1, x3^-1*x2^-1*x0^-1*x2*x3*x2^-1*x0*x2*x3*x2^-1,
                                x3^-1*(x2^-1*x0*x2*x3)^2*x2^-1*x0^-1*x2*x3^-1*x2^-1*x0^-1*x2,
                                x3*x2^-1*x1*x2*x3^-1*x2^-1*x1^-1*x2, x0*x1*x2*x3 >
            sage: A.meridians(simplified=False)
            {0: [x2, x3], 1: [x1], 2: [x0]}
            sage: A = H(y^2 + x*z, z, x)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1 | (x1*x0)^2*(x1^-1*x0^-1)^2 >
            sage: A = H(y^2 + x*z, z*x, y)
            sage: A.fundamental_group()
            Finitely presented group
            < x0, x1, x2 | x2*x0*x1*x0^-1*x2^-1*x1^-1,
                           x1*(x2*x0)^2*x2^-1*x1^-1*x0^-1*x2^-1*x0^-1 >

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if self._fundamental_group and not simplified:
            return self._fundamental_group
        if self._fundamental_group_simplified and simplified:
            return self._fundamental_group_simplified
        H = self.parent()
        K = self.base_ring()
        R = self.coordinate_ring()
        x, y, z = R.gens()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        C = self.reduce()
        n = len(C)
        infinity = Curve(z)
        infinity_in_C = infinity in C
        if infinity_in_C and n == 1:
            G = FreeGroup(0) / []
            self._fundamental_group = G
            self._meridians = {0: 0}
            self._fundamental_group_simplified = G
            self._meridians_simplified = {0: [G.one()]}
            return G
        if infinity_in_C:
            j = C.curves().index(infinity)
            C = H(C.curves()[:j] + C.curves()[j + 1:])
        infinity_divides = False
        for j, c in enumerate(C):
            g = c.defining_polynomial()
            infinity_divides = z.divides(g)
            if infinity_divides:
                h = R(g / z)
                C = H(C.curves()[:j] + (h, ) + C.curves()[j + 1:])
                break
        affine = AffinePlaneCurveArrangements(K, names=('u', 'v'))
        u, v = affine.gens()
        affines = [f.defining_polynomial().subs({x: u, y: v, z: 1}) for f in C]
        changes = any([g.degree(v) < g.degree() > 1 for g in affines])
        while changes:
            affines = [f.subs({u: u + v}) for f in affines]
            changes = any([g.degree(v) < g.degree() > 1 for g in affines])
        C_affine = affine(affines)
        proj = not (infinity_divides or infinity_in_C)
        G = C_affine.fundamental_group(simplified=simplified, vertical=True,
                                       projective=proj)
        dic = C_affine.meridians(simplified=simplified, vertical=True)
        if infinity_in_C:
            dic1 = dict()
            for k in range(j):
                dic1[k] = dic[k]
            dic1[j] = dic[n - 1]
            for k in range(j + 1, n):
                dic1[k] = dic[k - 1]
        elif infinity_divides:
            dic1 = {k: dic[k] for k in range(n)}
            dic1[j] += dic[n]
        else:
            dic1 = dic
        if not simplified:
            self._fundamental_group = G
            self._meridians = dic1
        else:
            self._fundamental_group_simplified = G
            self._meridians_simplified = dic1
        return G

    def meridians(self, simplified=True):
        r"""
        Meridians of each irreducible component.

        OUTPUT:

        A dictionary which associates the index of each curve with
        its meridians, including the line at infinity if it can be
        computed

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x*z, y + x - z, x)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1 | x1^-1*x0*x1*x0^-1 >
            sage: A.meridians()
            {0: [x1], 1: [x1^-1*x0^-1*x1^-1], 2: [x0]}
            sage: A = H(y^2 + x*z, z, x)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1 | (x1*x0)^2*(x1^-1*x0^-1)^2 >
            sage: A.meridians()
            {0: [x0, x1*x0*x1^-1], 1: [x0^-1*x1^-1*x0^-1], 2: [x1]}
            sage: A = H(y^2 + x*z, z*x, y)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1, x2 | x2*x0*x1*x0^-1*x2^-1*x1^-1,
                                                    x1*(x2*x0)^2*x2^-1*x1^-1*x0^-1*x2^-1*x0^-1 >
            sage: A.meridians()
            {0: [x0, x2*x0*x2^-1], 1: [x2, x0^-1*x2^-1*x1^-1*x0^-1], 2: [x1]}

        This function needs
        :func:`ProjectivePlaneCurveArrangements.fundamental_group`
        with the same options, where some examples are shown.

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        if not simplified:
            computed = bool(self._meridians)
        else:
            computed = bool(self._meridians_simplified)
        if not computed:
            _ = self._fundamental_group(simplified=simplified)
        if not simplified:
            return self._meridians
        return self._meridians_simplified


class AffinePlaneCurveArrangements(Parent, UniqueRepresentation):
    """
    Affine curve arrangements.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
        sage: H(x, y^2, x-1, y-1)
        Arrangement (x, y^2, x - 1, y - 1) in Affine Space
        of dimension 2 over Rational Field
    """
    Element = AffinePlaneCurveArrangementsElement

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: K = AffinePlaneCurveArrangements(QQ, names=('x', 'y'))
            sage: H is K
            True
            sage: type(K)
            <class 'sage.schemes.curves.plane_curve_arrangement.AffinePlaneCurveArrangements_with_category'>

        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
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

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.base_ring()
            Rational Field
        """
        return self._base_ring

    def coordinate_ring(self):
        """
        Return the coordinate ring.

        OUTPUT:

        The coordinate ring of the curve arrangement.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return PolynomialRing(self._base_ring, self._names)

    def change_ring(self, base_ring):
        """
        Return curve arrangements over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring.

        OUTPUT:

        A new :class:`AffinePlaneCurveArrangements` instance over the new
        base ring.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.gen(0)
            x
            sage: L.change_ring(RR).base_ring()
            Real Field with 53 bits of precision

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return AffinePlaneCurveArrangements(base_ring, names=self._names)

    @cached_method
    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
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

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ);  L
            Curve arrangements in Affine Space of dimension 2 over Rational Field
        """
        return 'Curve arrangements in {0}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a curve

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = L._element_constructor_(x, y); A
            Arrangement (x, y) in Affine Space of dimension 2 over Rational Field
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
        ambient_space = self.ambient_space()
        R = ambient_space.coordinate_ring()
        curves = ()
        for h in arg:
            try:
                ambient = h.ambient_space()
                if ambient == ambient_space:
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
        Return an element of ``self``.

        TESTS::

            sage: H.<t, s> = AffinePlaneCurveArrangements(QQ)
            sage: H._an_element_()
            Arrangement (t) in Affine Space of dimension 2 over Rational Field
            sage: H.<t, s, r> = ProjectivePlaneCurveArrangements(QQ)
            sage: H._an_element_()
            Arrangement (t) in Projective Space of dimension 2 over Rational Field
        """
        x = self.gen(0)
        return self(x)

    @cached_method
    def ngens(self):
        """
        Return the number of variables, i.e. 2 or 3, kept for completness.

        OUTPUT:

        An integer, 2 or 3, depending if the arrangement is projective or affine.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.ngens()
            2
            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.ngens()
            3
        """
        return len(self._names)

    @cached_method
    def gens(self):
        """
        Return the coordinates.

        OUTPUT:

        A tuple of linear expressions, one for each linear variable.

        EXAMPLES::

            sage: L = AffinePlaneCurveArrangements(QQ, ('x', 'y'))
            sage: L.gens()
            (x, y)
            sage: L = ProjectivePlaneCurveArrangements(QQ, ('x', 'y', 'z'))
            sage: L.gens()
            (x, y, z)
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

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.gen(1)
            y
            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.gen(2)
            z
        """
        return self.gens()[i]


class ProjectivePlaneCurveArrangements(AffinePlaneCurveArrangements):
    """
    Projective curve arrangements.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
        sage: H(x, y^2, x-z, y-z)
        Arrangement (x, y^2, x - z, y - z) in Projective Space
        of dimension 2 over Rational Field
    """
    Element = ProjectivePlaneCurveArrangementsElement

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: K = ProjectivePlaneCurveArrangements(QQ, names=('x', 'y', 'z'))
            sage: H is K
            True
            sage: type(K)
            <class 'sage.schemes.curves.plane_curve_arrangement.ProjectivePlaneCurveArrangements_with_category'>

        TESTS::

            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: TestSuite(H).run()
        """
        if base_ring not in _Fields:
            raise ValueError('base ring must be a field')
        super().__init__(base_ring, names=names)
        self._base_ring = base_ring
        self._names = names

    @cached_method
    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Projective Space of dimension 2 over Rational Field
        """
        return ProjectiveSpace(self.base_ring(), 2, self._names)
