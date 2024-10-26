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

The default base field is `\QQ`, the rational numbers.
Number fields are also possible (also with fixed embeddings in
`\QQbar`)::

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

from itertools import combinations
from sage.categories.sets_cat import Sets
from sage.groups.free_group import FreeGroup
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar
from sage.rings.ring import _Fields
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.curves.affine_curve import AffinePlaneCurve
from sage.schemes.curves.constructor import Curve
from sage.schemes.curves.projective_curve import ProjectiveSpace, ProjectivePlaneCurve
from sage.schemes.curves.zariski_vankampen import braid_monodromy, fundamental_group_arrangement
from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


class PlaneCurveArrangementElement(Element):
    """
    An ordered plane curve arrangement.
    """
    def __init__(self, parent, curves, check=True):
        """
        Construct a plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`PlaneCurveArrangements`

        - ``curves`` -- tuple of curves

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement (x, y) in Affine Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement (x, y) in Projective Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
        """
        super().__init__(parent)
        self._curves = tuple(curves)
        if check:
            affine = all(isinstance(h, AffinePlaneCurve) for h in curves)
            projective = all(isinstance(h, ProjectivePlaneCurve) for h in curves)
            if not (affine or projective):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is parent.ambient_space()
                       for h in curves):
                raise ValueError("not all curves are in the same ambient space")

    def __getitem__(self, i):
        """
        Return the `i`-th curve.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: H(y^2 - x, y^3 + 2 * x^2, x^4 + y^4 + 1)
            Arrangement (y^2 - x, y^3 + 2*x^2, x^4 + y^4 + 1)
            in Affine Space of dimension 2 over Rational Field
        """
        return self._curves[i]

    def __hash__(self):
        r"""
        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: H((x * y, x + y +1)).__hash__()   # random
            -4938643871296220686
        """
        return hash(self.curves())

    def ncurves(self):
        r"""
        Return the number of curves in the arrangement.

        OUTPUT: integer

        EXAMPLES::

            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: h = H((x * y, x + y + z))
            sage: h.ncurves()
            2
            sage: len(h)    # equivalent
            2
        """
        return len(self._curves)

    __len__ = ncurves

    def curves(self):
        r"""
        Return the curves in the arrangement as a tuple.

        OUTPUT: a tuple

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

        OUTPUT: string

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: h
            Arrangement of 5 curves in Affine Space of dimension 2 over Rational Field
            sage: H(())
            Arrangement () in Affine Space of dimension 2 over Rational Field
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + z, x^3 * z^2 - y^5, x^2 * y^2 * z + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - z^3)^2])
            sage: h
            Arrangement of 5 curves in Projective Space of dimension 2 over Rational Field
        """
        if not self:
            return 'Empty curve arrangement in {}'.format(self.parent().ambient_space())
        elif len(self) < 5:
            curves = ', '.join(h.defining_polynomial()._repr_()
                               for h in self._curves)
            return 'Arrangement ({}) in {}'.format(curves,
                                                     self.parent().ambient_space())
        return 'Arrangement of {} curves in {}'.format(len(self),
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

        OUTPUT: a new curve arrangement

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: h = H([x * y, x + y + 1, x^3 - y^5, x^2 * y^2 + x^5 + y^5, (x^2 + y^2)^3 + (x^3 + y^3 - 1)^2])
            sage: C = Curve(x^8 - y^8 -x^4 * y^4)
            sage: h1 = h.union(C); h1
            Arrangement of 6 curves in Affine Space of dimension 2 over Rational Field
            sage: h1 == h1.union(C)
            True
        """
        P = self.parent()
        other_h = P(other)
        curves0 = self._curves + other_h._curves
        curves = []
        for h in curves0:
            if h not in curves:
                curves.append(h)
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
            sage: h.deletion(x)
            Traceback (most recent call last):
            ...
            ValueError: curve is not in the arrangement
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

            sage: # needs sage.rings.number_field
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: K.<a> = CyclotomicField(3)
            sage: A.change_ring(K)
            Arrangement (-x^3 + y^2, x, y, x^2 + x*y + y^2) in Affine Space of
            dimension 2 over Cyclotomic Field of order 3 and degree 2
        """
        parent = self.parent().change_ring(base_ring)
        curves = tuple(c.change_ring(base_ring) for c in self)
        return parent(curves)

    @cached_method
    def coordinate_ring(self):
        """
        Return the coordinate ring of ``self``.

        OUTPUT: the coordinate ring of the curve arrangement

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: C = L(x, y)
            sage: C.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: P.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: C = P(x, y)
            sage: C.coordinate_ring()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self._curves[0].defining_polynomial().parent()

    def defining_polynomials(self):
        r"""
        Return the defining polynomials of the elements of ``self``.

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 - x^3, x, y, y^2 + x * y + x^2)
            sage: A.defining_polynomials()
            (-x^3 + y^2, x, y, x^2 + x*y + y^2)
        """
        return tuple([h.defining_polynomial() for h in self._curves])

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
            sage: H(x * y, x + y^3).have_common_factors()
            False
        """
        L = [c.defining_polynomial() for c in self._curves]
        C = combinations(L, 2)
        return any(f1.gcd(f2).degree() > 0 for f1, f2 in C)

    def reduce(self, clean=False, verbose=False):
        r"""
        Replace the curves by their reduction.

        INPUT:

        - ``clean`` -- boolean (default: ``False``); if ``False``
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
            sage: C.reduce(verbose=True)
            Some curves have common components
            sage: C.reduce(clean=True)
            Arrangement (x*y, y + 1) in Affine Space of dimension 2
            over Rational Field
            sage: C = H(x*y, x)
            sage: C.reduce(clean=True)
            Arrangement (x*y) in Affine Space of dimension 2 over Rational Field
        """
        P = self.parent()
        L = [self._curves[0].defining_polynomial().radical()]
        for c in self._curves[1:]:
            g = c.defining_polynomial().radical()
            for f in L:
                d = g.gcd(f)
                if d.degree() > 0 and not clean:
                    if verbose:
                        print("Some curves have common components")
                    return None
                g //= d
            if g.degree() > 0:
                L.append(g)
        return P(*L)


class AffinePlaneCurveArrangementElement(PlaneCurveArrangementElement):
    """
    An ordered affine plane curve arrangement.
    """
    def __init__(self, parent, curves, check=True):
        """
        Construct an ordered affine plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`AffinePlaneCurveArrangements`

        - ``curves`` -- tuple of curves

        EXAMPLES::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: elt = H(x, y); elt
            Arrangement (x, y) in Affine Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
        """
        Element.__init__(self, parent)
        self._curves = tuple(curves)
        if check:
            if not all(isinstance(h, AffinePlaneCurve) for h in curves):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is self.parent().ambient_space()
                       for h in curves):
                raise ValueError("not all curves are in the same ambient space")
        self._braid_monodromy_non_vertical = None
        self._braid_monodromy_vertical = None
        self._strands_nonvertical = None
        self._strands_vertical = None
        self._fundamental_group_nonsimpl_nonvertical = None
        self._fundamental_group_nonsimpl_vertical = None
        self._fundamental_group_simpl_nonvertical = None
        self._fundamental_group_simpl_vertical = None
        self._meridians_nonsimpl_nonvertical = None
        self._meridians_nonsimpl_vertical = None
        self._meridians_simpl_nonvertical = None
        self._meridians_simpl_vertical = None
        self._vertical_lines_in_braid_mon = None

    def fundamental_group(self, simplified=True, vertical=True,
                          projective=False):
        r"""
        Return the fundamental group of the complement of the union
        of affine plane curves in `\CC^2`.

        INPUT:

        - ``vertical`` -- boolean (default: ``True``); if ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified braid :func:`braid_monodromy` is used

        - ``simplified`` -- boolean (default: ``True``); if it is ``True``, the
          group is simplified

        - ``projective`` -- boolean (default: ``False``); to be used in the
          method for projective curves

        OUTPUT: a finitely presented group

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

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
            sage: G = A.fundamental_group(simplified=False)
            sage: G.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3 | x3^-1*x2^-1*x3*x0*x1*x0^-1,
                               x3^-1*x1^-1*x3*x0*x1*x0^-1*x2^-1*x0^-1*(x2*x0)^2*x1^-1*x0^-1,
                               x3^-1*x0^-1*x3*x0*x1*x0^-1*x2^-1*x0*x2*x0*x1^-1*x0^-1,
                               x2^-1*x0^-1*x2*x0, x1^-1*x0^-1*x1*x0 >
            sage: A.meridians(simplified=False)
            {0: [x1, x2], 1: [x0], 2: [x3], 3: [x3^-1*x2^-1*x1^-1*x0^-1]}
            sage: A.fundamental_group(vertical=False)
            Finitely presented group
            < x0, x1, x2 | x2^-1*x1^-1*x2*x1, x1*x0*x1^-1*x0^-1, (x0*x2)^2*(x0^-1*x2^-1)^2 >
            sage: A.meridians(vertical=False)
            {0: [x2, x0*x2*x0^-1], 1: [x1], 2: [x0], 3: [x0*x2^-1*x0^-1*x2^-1*x1^-1*x0^-1]}
            sage: G = A.fundamental_group(simplified=False, vertical=False)
            sage: G.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3 | x3^-1*x2^-1*x1^-1*x2*x3*x2^-1*x1*x2,
                               x3^-1*x2^-1*x1^-1*x2*x3*x2^-1*x1*x2,
                               (x3^-1*x2^-1*x0^-1*x2)^2*(x3*x2^-1*x0*x2)^2,
                               x3^-1*x2^-1*x0^-1*x2*x3^-1*x2^-1*x0*x2*x3*x2,
                               x1^-1*x0^-1*x1*x0 >
            sage: A.meridians(simplified=False, vertical=False)
            {0: [x2, x3], 1: [x1], 2: [x0], 3: [x3^-1*x2^-1*x1^-1*x0^-1]}
            sage: A = H(x * y^2 + x + y, y + x -1, x, y)
            sage: G = A.fundamental_group()
            sage: G.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3 | x3^-1*x2^-1*x3*x2, x3^-1*x1^-1*x3*x1,
                               x3^-1*x0^-1*x3*x0, x2^-1*x1^-1*x2*x1,
                               x2^-1*x0^-1*x2*x0, x1^-1*x0^-1*x1*x0 >
        """
        if simplified and vertical:
            computed = self._fundamental_group_simpl_vertical
        elif simplified and not vertical:
            computed = self._fundamental_group_simpl_nonvertical
        elif not simplified and vertical:
            computed = self._fundamental_group_nonsimpl_vertical
        else:
            computed = self._fundamental_group_nonsimpl_nonvertical
        if computed:
            return computed
        K = self.base_ring()
        R = self.coordinate_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        C = self.reduce()
        L = C.defining_polynomials()
        if vertical:
            bm = self._braid_monodromy_vertical
        else:
            bm = self._braid_monodromy_non_vertical
        if bm is not None:  # bm could be []
            if not vertical:
                st = self._strands_nonvertical
                d1 = prod(L).degree()
                bd = (bm, st, {}, d1)
            else:
                st = self._strands_vertical
                d1 = prod(L).degree(R.gen(1))
                bd = (bm, st, self._vertical_lines_in_braid_mon, d1)
        else:
            bd = None
        G, dic = fundamental_group_arrangement(L, simplified=simplified,
                                               puiseux=True,
                                               projective=projective,
                                               vertical=vertical,
                                               braid_data=bd)
        if simplified and vertical:
            self._fundamental_group_simpl_vertical = G
            self._meridians_simpl_vertical = dic
        elif simplified and not vertical:
            self._fundamental_group_simpl_nonvertical = G
            self._meridians_simpl_nonvertical = dic
        elif not simplified and vertical:
            self._fundamental_group_nonsimpl_vertical = G
            self._meridians_nonsimpl_vertical = dic
        else:
            self._fundamental_group_nonsimpl_nonvertical = G
            self._meridians_nonsimpl_nonvertical = dic
        return G

    def meridians(self, simplified=True, vertical=True):
        r"""
        Return the meridians of each irreducible component.

        OUTPUT:

        A dictionary which associates the index of each curve with its meridians,
        including the line at infinity if it can be omputed

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed
           and :meth:`AffinePlaneCurveArrangements.fundamental_group` with the same options,
           where some examples are shown.

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
        """
        if simplified and vertical:
            computed = self._meridians_simpl_vertical
        elif simplified and not vertical:
            computed = self._meridians_simpl_nonvertical
        elif not simplified and vertical:
            computed = self._meridians_nonsimpl_vertical
        else:
            computed = self._meridians_nonsimpl_nonvertical
        if computed:
            return dict(computed)
        self.fundamental_group(simplified=simplified, vertical=vertical)
        if simplified and vertical:
            return dict(self._meridians_simpl_vertical)
        elif simplified and not vertical:
            return dict(self._meridians_group_simpl_nonvertical)
        elif not simplified and vertical:
            return dict(self._meridians_nonsimpl_vertical)
        else:
            return dict(self._meridians_nonsimpl_nonvertical)

    def braid_monodromy(self, vertical=True):
        r"""
        Return the braid monodromy of the complement of the union
        of affine plane curves in `\CC^2`. If there are vertical
        asymptotes a change of variable is done.

        INPUT:

        - ``vertical`` -- boolean (default: ``True``); if it is ``True``, there
          are no vertical asymptotes, and there are vertical lines, then a
          simplified :func:`braid_monodromy` is computed.

        OUTPUT:

        A braid monodromy with dictionnaries identifying strands with components
        and braids with vertical lines.

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

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
        """
        if vertical:
            computed = self._braid_monodromy_vertical
        else:
            computed = self._braid_monodromy_non_vertical
        if computed is not None:
            return computed
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        L = self.defining_polynomials()
        bm, dic, dv, d1 = braid_monodromy(prod(L), arrangement=L,
                                          vertical=vertical)
        if vertical:
            self._braid_monodromy_vertical = bm
            self._strands_vertical = dic
            self._vertical_lines_in_braid_mon = dv
        else:
            self._braid_monodromy_non_vertical = bm
            self._strands_nonvertical = dic
        return bm

    def strands(self):
        r"""
        Return the strands for each member of the arrangement.

        OUTPUT:

        A dictionary which associates to the index of each strand
        its associated component if the braid monodromy has been
        calculated with ``vertical=False``.

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: bm = A.braid_monodromy()
            sage: A.strands()
            {0: 2, 1: 1, 2: 0, 3: 0}
        """
        if not self._strands_nonvertical:
            self._braid_monodromy = self.braid_monodromy(vertical=False)
        return self._strands_nonvertical

    def vertical_strands(self):
        r"""
        Return the strands if the braid monodromy has been computed with
        the vertical option.

        OUTPUT:

        A dictionary which associates to the index of each strand
        its associated component if the braid monodromy has been
        calculated with ``vertical=True``.

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.vertical_strands()
            {0: 1, 1: 0, 2: 0}
            sage: A.braid_monodromy(vertical=True)
            [s1*s0*s1*s0^-1*s1^-1*s0, s0^-1*s1*s0*s1^-1*s0, s0^-1*s1^2*s0]
        """
        if not self._strands_vertical:
            self.braid_monodromy(vertical=True)
        return self._strands_vertical

    def vertical_lines_in_braid_monodromy(self):
        r"""
        Return the vertical lines in the arrangement.

        OUTPUT:

        A dictionary which associates the index of a braid
        to the index of the vertical line associated to the braid.

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x, y + x - 1, x)
            sage: A.vertical_lines_in_braid_monodromy()
            {1: 2}
            sage: A.braid_monodromy(vertical=True)
            [s1*s0*s1*s0^-1*s1^-1*s0, s0^-1*s1*s0*s1^-1*s0, s0^-1*s1^2*s0]
        """
        if not self._vertical_lines_in_braid_mon:
            self.braid_monodromy(vertical=True)
        return self._vertical_lines_in_braid_mon


class ProjectivePlaneCurveArrangementElement(PlaneCurveArrangementElement):
    """
    An ordered projective plane curve arrangement.
    """
    def __init__(self, parent, curves, check=True):
        """
        Construct an ordered projective plane curve arrangement.

        INPUT:

        - ``parent`` -- the parent :class:`ProjectivePlaneCurveArrangements`

        - ``curves`` -- tuple of curves

        EXAMPLES::

            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: elt = H(x, y, z); elt
            Arrangement (x, y, z) in Projective Space of dimension 2 over Rational Field
            sage: TestSuite(elt).run()
        """
        Element.__init__(self, parent)
        self._curves = tuple(curves)
        if check:
            if not all(isinstance(h, ProjectivePlaneCurve) for h in curves):
                raise ValueError("not all elements are curves")
            if not all(h.ambient_space() is self.parent().ambient_space()
                       for h in curves):
                raise ValueError("not all curves are in the same ambient space")
        self._fundamental_group_nonsimpl = None
        self._fundamental_group_simpl = None
        self._meridians_nonsimpl = None
        self._meridians_simpl = None

    def fundamental_group(self, simplified=True):
        r"""
        Return the fundamental group of the complement of the union
        of an arragnement of projective plane curves
        in the projective plane.

        INPUT:

        - ``simplified`` -- boolean (default: ``True``); set if the group
          is simplified

        OUTPUT: a finitely presented group

        .. NOTE::

           This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: H(z).fundamental_group()
            Finitely presented group <  |  >
            sage: H(x*y).fundamental_group()
            Finitely presented group < x |  >
            sage: A = H(y^2 + x*z, y + x - z, x)
            sage: A.fundamental_group().sorted_presentation()
            Finitely presented group < x0, x1 | x1^-1*x0^-1*x1*x0 >
            sage: A.meridians()
            {0: [x1], 1: [x0], 2: [x1^-1*x0^-1*x1^-1]}
            sage: G = A.fundamental_group(simplified=False)
            sage: G.sorted_presentation()
            Finitely presented group
            < x0, x1, x2, x3 | x3^-1*x2^-1*x1^-1*x0^-1, x3^-1*x2^-1*x3*x0*x1*x0^-1,
                               x3^-1*x1^-1*x3*x0*x1*x0^-1*x2^-1*x0^-1*(x2*x0)^2*x1^-1*x0^-1,
                               x3^-1*x0^-1*x3*x0*x1*x0^-1*x2^-1*x0*x2*x0*x1^-1*x0^-1,
                               x2^-1*x0^-1*x2*x0, x1^-1*x0^-1*x1*x0 >
            sage: A.meridians(simplified=False)
            {0: [x1, x2], 1: [x0], 2: [x3]}
            sage: A = H(y^2 + x*z, z, x)
            sage: A.fundamental_group()
            Finitely presented group < x0, x1 | (x1*x0)^2*(x1^-1*x0^-1)^2 >
            sage: A = H(y^2 + x*z, z*x, y)
            sage: A.fundamental_group()
            Finitely presented group
            < x0, x1, x2 | x2*x0*x1*x0^-1*x2^-1*x1^-1,
                           x1*(x2*x0)^2*x2^-1*x1^-1*x0^-1*x2^-1*x0^-1 >
        """
        if simplified:
            computed = self._fundamental_group_simpl
        else:
            computed = self._fundamental_group_nonsimpl
        if computed:
            return computed
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
            if simplified:
                self._fundamental_group_simpl = G
                self._meridians_simpl = {0: [G.one()]}
            else:
                self._fundamental_group_nonsimpl = G
                self._meridians_nonsimpl = {0: [G.one()]}
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
        changes = any(g.degree(v) < g.degree() > 1 for g in affines)
        while changes:
            affines = [f.subs({u: u + v}) for f in affines]
            changes = any(g.degree(v) < g.degree() > 1 for g in affines)
        C_affine = affine(affines)
        proj = not (infinity_divides or infinity_in_C)
        G = C_affine.fundamental_group(simplified=simplified, vertical=True,
                                       projective=proj)
        dic = C_affine.meridians(simplified=simplified, vertical=True)
        if infinity_in_C:
            dic1 = {}
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
        if simplified:
            self._fundamental_group_simpl = G
            self._meridians_simpl = dic1
        else:
            self._fundamental_group_nonsimpl = G
            self._meridians_nonsimpl = dic1
        return G

    def meridians(self, simplified=True):
        r"""
        Return the meridians of each irreducible component.

        OUTPUT:

        A dictionary which associates the index of each curve with
        its meridians, including the line at infinity if it can be
        computed

        .. NOTE::

           This function requires the ``sirocco`` package to be installed and
           :func:`ProjectivePlaneCurveArrangements.fundamental_group`
           with the same options, where some examples are shown.

        EXAMPLES::

            sage: # needs sirocco
            sage: H.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: A = H(y^2 + x*z, y + x - z, x)
            sage: A.fundamental_group().sorted_presentation()
            Finitely presented group < x0, x1 | x1^-1*x0^-1*x1*x0 >
            sage: A.meridians()
            {0: [x1], 1: [x0], 2: [x1^-1*x0^-1*x1^-1]}
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
        """
        if simplified:
            computed = self._meridians_simpl
        else:
            computed = self._meridians_nonsimpl
        if computed:
            return dict(computed)
        self._fundamental_group(simplified=simplified)
        if simplified:
            return dict(self._meridians_simpl)
        else:
            return dict(self._meridians_nonsimpl)


class PlaneCurveArrangements(UniqueRepresentation, Parent):
    """
    Plane curve arrangements.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
        sage: H(x, y^2, x-1, y-1)
        Arrangement (x, y^2, x - 1, y - 1) in Affine Space
        of dimension 2 over Rational Field
    """
    Element = PlaneCurveArrangementElement

    @staticmethod
    def __classcall__(cls, base, names=()):
        """
        Normalize the inputs to ensure a unique representation.

        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: K = AffinePlaneCurveArrangements(QQ, names=('x', 'y'))
            sage: H is K
            True
        """
        names = normalize_names(len(names), names)
        return super().__classcall__(cls, base, names)

    def __init__(self, base_ring, names=()):
        """
        Initialize ``self``.

        TESTS::

            sage: H.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: TestSuite(H).run()
        """
        if base_ring not in _Fields:
            raise ValueError('base ring must be a field')
        super().__init__(base_ring, names=names, category=Sets())
        self._embedded = None
        if len(names) == 2:
            self._embedded = 'affine'
        elif len(names) == 3:
            self._embedded = 'projective'

    def coordinate_ring(self):
        """
        Return the coordinate ring.

        OUTPUT: the coordinate ring of the curve arrangement

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return PolynomialRing(self.base_ring(), self.variable_names())

    def change_ring(self, base_ring):
        """
        Return curve arrangements over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring; the new base ring

        OUTPUT:

        A new :class:`PlaneCurveArrangements` instance over the new
        base ring

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.change_ring(RR).base_ring()
            Real Field with 53 bits of precision

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return self.__reduce__()[1][0](base_ring, names=self.variable_names())

    @abstract_method
    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y> = PlaneCurveArrangements(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method ambient_space at  0x...>
            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Affine Space of dimension 2 over Rational Field
            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Projective Space of dimension 2 over Rational Field
        """

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ);  L
            Curve arrangements in Affine Space of dimension 2 over Rational Field
        """
        return 'Curve arrangements in {}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a curve

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: A = L(x, y); A
            Arrangement (x, y) in Affine Space of dimension 2 over Rational Field
            sage: L([x, y]) == A
            True
            sage: L(Curve(x), Curve(y)) == A
            True
            sage: L(y, x) == A
            False
       """
        if len(args) == 1:
            if not isinstance(args[0], (tuple, list)):
                arg = (args[0], )
            else:
                arg = tuple(args[0])
        else:
            arg = tuple(args)
        ambient_space = self.ambient_space()
        R = ambient_space.coordinate_ring()
        curves = ()
        for h in arg:
            try:
                ambient = h.ambient_space()
                if ambient == ambient_space:
                    curves += (h, )
                else:
                    raise TypeError('the curves do not have the same ambient space')
            except AttributeError:
                try:
                    h = R(h)
                    curves += (Curve(h), )
                except TypeError:
                    raise TypeError('elements are not curves')
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
        return self(self.gen(0))

    @cached_method
    def ngens(self):
        """
        Return the number of variables, i.e. 2 or 3, kept for completeness.

        OUTPUT: integer, 2 or 3, depending if the arrangement is projective or affine

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.ngens()
            2
            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.ngens()
            3
        """
        return len(self.variable_names())

    def gens(self):
        """
        Return the coordinates.

        OUTPUT: a tuple of linear expressions, one for each linear variable

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

        OUTPUT: a variable

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.gen(1)
            y
            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.gen(2)
            z
        """
        return self.gens()[i]


class AffinePlaneCurveArrangements(PlaneCurveArrangements):
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
    Element = AffinePlaneCurveArrangementElement

    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y> = AffinePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Affine Space of dimension 2 over Rational Field
        """
        return AffineSpace(self.base_ring(), 2, self._names)


class ProjectivePlaneCurveArrangements(PlaneCurveArrangements):
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
    Element = ProjectivePlaneCurveArrangementElement

    def ambient_space(self):
        """
        Return the ambient space.

        EXAMPLES::

            sage: L.<x, y, z> = ProjectivePlaneCurveArrangements(QQ)
            sage: L.ambient_space()
            Projective Space of dimension 2 over Rational Field
        """
        return ProjectiveSpace(self.base_ring(), 2, self._names)
