# -*- coding: utf-8 -*-
r"""
Affine Plane Curve Arrangements

We create an :class:`OrderedAffinePlaneCurveArrangements`
object
# to define the variables `x`, `y`, `z`::
#
#     sage: H.<x,y,z> = HyperplaneArrangements(QQ)
#     sage: h = 3*x + 2*y - 5*z - 7;  h
#     Hyperplane 3*x + 2*y - 5*z - 7
#     sage: h.normal()
#     (3, 2, -5)
#     sage: h.constant_term()
#     -7

The individual curves will be in  :class:`AffinePlaneCurve`::

    # sage: -2*h
    # Hyperplane -6*x - 4*y + 10*z + 14
    # sage: x, y, z
    # (Hyperplane x + 0*y + 0*z + 0,
    #  Hyperplane 0*x + y + 0*z + 0,
    #  Hyperplane 0*x + 0*y + z + 0)

The default base field is `\QQ`, the rational numbers.
Number fields are also possible (also with fixed embeddings in ``QQbar``)::

    # sage: # needs sage.rings.number_field
    # sage: x = polygen(QQ, 'x')
    # sage: NF.<a> = NumberField(x**4 - 5*x**2 + 5, embedding=1.90)
    # sage: H.<y,z> = HyperplaneArrangements(NF)
    # sage: A = H([[(-a**3 + 3*a, -a**2 + 4), 1], [(a**3 - 4*a, -1), 1],
    # ....:        [(0, 2*a**2 - 6), 1], [(-a**3 + 4*a, -1), 1],
    # ....:        [(a**3 - 3*a, -a**2 + 4), 1]])
    # sage: A
    # Arrangement of 5 hyperplanes of dimension 2 and rank 2
    # sage: A.base_ring()
    # Number Field in a with defining polynomial x^4 - 5*x^2 + 5
    #  with a = 1.902113032590308?


AUTHORS:

- Enrique Artal (2023-10): initial version
"""

# *****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Possible extensions:


# from sage.categories.fields import Fields
# from sage.categories.homset import Hom, End, hom
# from sage.categories.number_fields import NumberFields
from sage.categories.sets_cat import Sets
from sage.combinat.combination import Combinations
# from sage.combinat.permutation import Permutation
# from sage.groups.free_group import FreeGroup
# from sage.interfaces.singular import singular
# from sage.matrix.constructor import matrix, vector
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
# from sage.rings.infinity import infinity
# from sage.misc.lazy_attribute import lazy_attribute
# from sage.rings.number_field.number_field import NumberField
# from sage.rings.polynomial.multi_polynomial_element import degree_lowest_rational_function
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.rings.qqbar import number_field_elements_from_algebraics
from sage.rings.qqbar import QQbar
# from sage.rings.rational_field import QQ, is_RationalField
# from sage.rings.rational_field import is_RationalField
from sage.rings.ring import _Fields
from sage.schemes.affine.affine_space import AffineSpace
# from sage.schemes.affine.affine_space import is_AffineSpace
# from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
# from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine_field
from sage.schemes.curves.affine_curve import AffinePlaneCurve
from sage.schemes.curves.constructor import Curve
from sage.schemes.curves.zariski_vankampen import fundamental_group_arrangement
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation

# from .closed_point import IntegralAffineCurveClosedPoint
# from .curve import Curve_generic
# from .point import (AffineCurvePoint_field,
#                     AffinePlaneCurvePoint_field,
#                     AffinePlaneCurvePoint_finite_field,
#                     IntegralAffineCurvePoint,
#                     IntegralAffineCurvePoint_finite_field,
#                     IntegralAffinePlaneCurvePoint,
#                     IntegralAffinePlaneCurvePoint_finite_field)


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

        - ``hyperplanes`` -- a tuple of hyperplanes

        - ``check`` -- boolean (optional; default ``True``); whether
          to check input

        - ``backend`` -- string (optional; default: ``None``); the backend to
          use for the related polyhedral objects

        # EXAMPLES::
        #
        #     sage: H.<x,y> = HyperplaneArrangements(QQ)
        #     sage: elt = H(x, y); elt
        #     Arrangement <y | x>
        #     sage: TestSuite(elt).run()
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

    def _first_ngens(self, n):
        """
        Workaround to support the construction with names.

        INPUT/OUTPUT:

        See :meth:`OrderedAffinePlaneCurveArrangements._first_ngens`.

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.braid(3)   # indirect doctest     # needs sage.graphs
            sage: (x, y) == a._first_ngens(2)                                           # needs sage.graphs
            True
        """
        return self.parent()._first_ngens(n)

    def __getitem__(self, i):
        """
        Return the `i`-th curve.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        The `i`-th curve.

        EXAMPLES::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: h = x|y;  h
            # Arrangement <y | x>
            # sage: h[0]
            # Hyperplane 0*x + y + 0
            # sage: h[1]
            # Hyperplane x + 0*y + 0
        """
        return self._curves[i]

    def __hash__(self):
        r"""
        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = x|y; h
            Arrangement <y | x>
            sage: len_dict = {h: len(h)}
        """
        return hash(self.curves())

    def n_curves(self):
        r"""
        Return the number of curves in the arrangement.

        OUTPUT:

        An integer.

        EXAMPLES::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            # sage: A.n_hyperplanes()
            # 3
            # sage: len(A)    # equivalent
            # 3
        """
        return len(self._curves)

    __len__ = n_curves

    def curves(self):
        r"""
        Return the curves in the arrangement as a tuple.

        OUTPUT:

        An tuple.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            sage: A.hyperplanes()
            (Hyperplane x + 0*y + 1, Hyperplane 3*x - y + 2, Hyperplane 5*x + 3*y + 4)

        Note that the hyperplanes can be indexed as if they were a list::

            sage: A[0]
            Hyperplane x + 0*y + 1
        """
        return self._curves

    def _repr_(self):
        r"""
        String representation for a curve arrangement.

        OUTPUT:

        A string.

        EXAMPLES::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: H(x, y, x-1, y-1)
            # Arrangement <y - 1 | y | x - 1 | x>
            # sage: x | y | x - 1 | y - 1 | x + y | x - y
            # Arrangement of 6 hyperplanes of dimension 2 and rank 2
            # sage: H()
            # Empty hyperplane arrangement of dimension 2
        """
        if len(self) == 0:
            return 'Empty curve arrangement'
        elif len(self) < 5:
            curves = ', '.join(h._repr_() for h in self._curves)
            return 'Arrangement ({0})'.format(curves)
        return 'Arrangement of {0} curves'.format(len(self))

    def _richcmp_(self, other, op):
        """
        Compare two hyperplane arrangements.

        EXAMPLES::

            # sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            # sage: H(x) == H(y)
            # False

        TESTS::

            # sage: H(x) == 0
            # False
        """
        return richcmp(self._curves, other._curves, op)

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a curvee arrangement or something that can
          be converted into a curve arrangement

        OUTPUT:

        A new curve arrangement.

        EXAMPLES::

        #     sage: H.<x,y> = HyperplaneArrangements(QQ)
        #     sage: H1.<x,y> = OrderedHyperplaneArrangements(QQ)
        #     sage: A = H([1,2,3], [0,1,1], [0,1,-1], [1,-1,0], [1,1,0])
        #     sage: B = H([1,1,1], [1,-1,1], [1,0,-1])
        #     sage: C = A.union(B); C
        #     Arrangement of 8 hyperplanes of dimension 2 and rank 2
        #     sage: C == A | B   # syntactic sugar
        #     True
        #     sage: A1 = H1(A)
        #     sage: B1 = H1(B)
        #     sage: C1 = A1.union(B1); C1
        #     Arrangement of 8 hyperplanes of dimension 2 and rank 2
        #     sage: [C1.hyperplanes().index(h) for h in C.hyperplanes()]
        #     [0, 5, 6, 1, 2, 3, 7, 4]
        #
        # A single hyperplane is coerced into a hyperplane arrangement
        # if necessary::
        #
        #     sage: A.union(x+y-1)
        #     Arrangement of 6 hyperplanes of dimension 2 and rank 2
        #     sage: A.add_hyperplane(x+y-1)    # alias
        #     Arrangement of 6 hyperplanes of dimension 2 and rank 2
        #
        #     sage: P.<x,y> = HyperplaneArrangements(RR)
        #     sage: C = P(2*x + 4*y + 5)
        #     sage: C.union(A)
        #     Arrangement of 6 hyperplanes of dimension 2 and rank 2
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

    def plot(self, **kwds):
        """
        Plot an arrangement of curves.

        OUTPUT:

        A graphics object.

        EXAMPLES::

            # sage: L.<x, y> = HyperplaneArrangements(QQ)
            # sage: L(x, y, x+y-2).plot()                                                 # needs sage.plot
            # Graphics object consisting of 3 graphics primitives
        """
        from sage.geometry.hyperplane_arrangement.plot import plot
        return plot(self, **kwds)

    def deletion(self, curves):
        r"""
        Return the curve arrangement obtained by removing ``h``.

        INPUT:

        - ``h`` -- a curve or curve arrangement

        OUTPUT:

        A new curve arrangement with the given curve(s)
        ``h`` removed.

        EXAMPLES::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: A = H([0,1,0], [1,0,1], [-1,0,1], [0,1,-1], [0,1,1]);  A
            # Arrangement of 5 hyperplanes of dimension 2 and rank 2
            # sage: A.deletion(x)
            # Arrangement <y - 1 | y + 1 | x - y | x + y>
            # sage: h = H([0,1,0], [0,1,1])
            # sage: A.deletion(h)
            # Arrangement <y - 1 | y + 1 | x - y>
        """
        parent = self.parent()
        curves = parent(curves)
        planes = list(self)
        for curve in curves:
            try:
                planes.remove(curve)
            except ValueError:
                raise ValueError('curve is not in the arrangement')
        return parent(*planes)

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

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: A = H([(1,1), 0], [(2,3), -1])
            # sage: A.change_ring(FiniteField(2))
            # Arrangement <y + 1 | x + y>
        """
        parent = self.parent().change_ring(base_ring)
        return parent(self)

    def defining_polynomials(self):
        r"""
        Return the defining polynomials of ``A``.

        # Let `A = (H_i)_i` be a hyperplane arrangement in a vector space `V`
        # corresponding to the null spaces of `\alpha_{H_i} \in V^*`. Then
        # the *defining polynomial* of `A` is given by
        #
        # .. MATH::
        #
        #     Q(A) = \prod_i \alpha_{H_i} \in S(V^*).

        EXAMPLES::

            # sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            # sage: A = H([2*x + y - z, -x - 2*y + z])
            # sage: p = A.defining_polynomial(); p
            # -2*x^2 - 5*x*y - 2*y^2 + 3*x*z + 3*y*z - z^2
            # sage: p.factor()
            # (-1) * (x + 2*y - z) * (2*x + y - z)
        """
        return tuple(h.defining_polynomial() for h in self)

    def defining_polynomial(self, simplified=True):
        r"""
        Return the defining polynomial of ``A``.

        # Let `A = (H_i)_i` be a hyperplane arrangement in a vector space `V`
        # corresponding to the null spaces of `\alpha_{H_i} \in V^*`. Then
        # the *defining polynomial* of `A` is given by
        #
        # .. MATH::
        #
        #     Q(A) = \prod_i \alpha_{H_i} \in S(V^*).

        EXAMPLES::

            # sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            # sage: A = H([2*x + y - z, -x - 2*y + z])
            # sage: p = A.defining_polynomial(); p
            # -2*x^2 - 5*x*y - 2*y^2 + 3*x*z + 3*y*z - z^2
            # sage: p.factor()
            # (-1) * (x + 2*y - z) * (2*x + y - z)
        """
        return prod(self.defining_polynomials())

    def have_common_factors(self):
        L = [c.defining_polynomial() for c in self]
        C = Combinations(L, 2)
        for f1, f2 in C:
            if f1.gcd(f2).degree() > 0:
                return True
        return False

    def reduce(self):
        P = self.parent()
        L = [c.defining_polynomial().radical() for c in self]
        return P(*L)

    def fundamental_group(self, vertical=True):
        r"""
        It computes the fundamental group of the complement of the union
        of affine projective curves in `\mathbb{C}^2`.

        OUTPUT:

        A group finitely presented with the assignation of each curve to
        a set of meridians, including the line at infinity.

        EXAMPLES::

            # sage: A.<x, y> = OrderedHyperplaneArrangements(QQ)
            # sage: L = [y + x, y + x - 1]
            # sage: H = A(L)
            # sage: G, dic = H.fundamental_group(); G                                   # optional - sirocco
            # Finitely presented group < x0, x1 |  >
            # sage: L = [x, y, x + 1, y + 1, x - y]
            # sage: H = A(L); list(H)
            # [Hyperplane x + 0*y + 0,
            #  Hyperplane 0*x + y + 0,
            #  Hyperplane x + 0*y + 1,
            #  Hyperplane 0*x + y + 1,
            #  Hyperplane x - y + 0]
            # sage: G, dic = H.fundamental_group()                                      # optional - sirocco
            # sage: G.simplified()                                                        # optional - sirocco
            # Finitely presented group < x0, x1, x2, x3, x4 | x3*x2*x3^-1*x2^-1,
            #                            x2^-1*x0^-1*x2*x4*x0*x4^-1,
            #                            x0*x1*x3*x0^-1*x3^-1*x1^-1,
            #                            x0*x2*x4*x2^-1*x0^-1*x4^-1,
            #                            x0*x1^-1*x0^-1*x3^-1*x1*x3,
            #                            x4^-1*x3^-1*x1*x3*x4*x3^-1*x1^-1*x3 >
            # sage: dic                                                                   # optional - sirocco
            #     {0: [x2], 1: [x4], 2: [x1], 3: [x3], 4: [x0], 5: [x4^-1*x3^-1*x2^-1*x1^-1*x0^-1]}
            # sage: H=A(x,y,x+y)
            # sage: H._fundamental_group_()                                               # optional - sirocco
            # (Finitely presented group < x0, x1, x2 | x1*x2*x0*x2^-1*x1^-1*x0^-1, x1*x0^-1*x2^-1*x1^-1*x2*x0 >,
            #  {0: [x0], 1: [x2], 2: [x1], 3: [x2^-1*x1^-1*x0^-1]})
            # sage: H._fundamental_group_(proj=True)                                      # optional - sirocco
            # (Finitely presented group < x0, x1 |  >, {1: (1,), 2: (2,), 3: (-2, -1)})
            # sage: A3.<x, y, z> = OrderedHyperplaneArrangements(QQ)
            # sage: H = A3(hyperplane_arrangements.braid(4).essentialization())               # optional - sage.graphs
            # sage: G, dic = H._fundamental_group_(proj=True)                             # optional - sage.graphs, sirocco
            # sage: h = G.simplification_isomorphism()                                    # optional - sage.graphs, sirocco
            # sage: G.simplified()                                                        # optional - sage.graphs, sirocco
            # Finitely presented group < x0, x1, x3, x4, x5 | x0*x3*x0^-1*x3^-1,
            #                                                 x1*x4*x1^-1*x4^-1,
            #                                                 x1*x5*x1^-1*x0^-1*x5^-1*x0,
            #                                                 x5*x3*x4*x3^-1*x5^-1*x4^-1,
            #                                                 x5^-1*x1^-1*x0*x1*x5*x0^-1,
            #                                                 x4*x5^-1*x4^-1*x3^-1*x5*x3 >
            # sage: {j: h(dic[j][0]) for j in dic.keys()}                                 # optional - sage.graphs, sirocco
            # {0: x5, 1: x0, 2: x1, 3: x3, 4: x4, 5: x0^-1*x5^-1*x4^-1*x1^-1*x3^-1}

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        K = self.base_ring()
        if not K.is_subring(QQbar):
            raise TypeError('the base field is not in QQbar')
        C = self.reduce()
        L = C.defining_polynomials()
        G, dic = fundamental_group_arrangement(L, puiseux=True, vertical=vertical)
        return (G, dic)


class OrderedAffinePlaneCurveArrangements(Parent, UniqueRepresentation):
    """
    Curve arrangements.

    INPUT:

    - ``base_ring`` -- ring; the base ring

    - ``names`` -- tuple of strings; the variable names

    EXAMPLES::

        sage: H.<x,y> = HyperplaneArrangements(QQ)
        sage: x
        Hyperplane x + 0*y + 0
        sage: x + y
        Hyperplane x + y + 0
        sage: H(x, y, x-1, y-1)
        Arrangement <y - 1 | y | x - 1 | x>
    """
    Element = OrderedAffinePlaneCurveArrangementsElement

    def __init__(self, base_ring, names=tuple()):
        """
        Initialize ``self``.

        TESTS::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: K = HyperplaneArrangements(QQ, names=('x', 'y'))
            # sage: H is K
            # True
            # sage: type(K)
            # <class 'sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements_with_category'>
            # sage: K.change_ring(RR).gen(0)
            # Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000

        TESTS::

            # sage: H.<x,y> = HyperplaneArrangements(QQ)
            # sage: TestSuite(H).run()
            # sage: K = HyperplaneArrangements(QQ)
            # sage: TestSuite(K).run()
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

            # sage: L.<x,y> = HyperplaneArrangements(QQ)
            # sage: L.base_ring()
            # Rational Field
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
            #
            # sage: L.<x,y> = HyperplaneArrangements(QQ)
            # sage: L.gen(0)
            # Hyperplane x + 0*y + 0
            # sage: L.change_ring(RR).gen(0)
            # Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000

        TESTS::

            # sage: L.change_ring(QQ) is L
            # True
        """
        return OrderedAffinePlaneCurveArrangements(base_ring, names=self._names)

    @cached_method
    def ambient_space(self):
        """
        Return the ambient space.

        # The ambient space is the parent of hyperplanes. That is, new
        # hyperplanes are always constructed internally from the ambient
        # space instance.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L.ambient_space()([(1,0), 0])
            Hyperplane x + 0*y + 0
            sage: L.ambient_space()([(1,0), 0]) == x
            True
        """
        return AffineSpace(self.base_ring(), 2, self._names)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            # sage: L.<x, y> = HyperplaneArrangements(QQ);  L
            # Hyperplane arrangements in 2-dimensional linear space over Rational Field with coordinates x, y
        """
        return 'Curve arrangements in {0}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of ``self``.

        INPUT:

        - ``*args`` -- positional arguments, each defining a
          hyperplane; alternatively, a single polytope or a single
          hyperplane arrangement

        - ``warn_duplicates`` -- boolean (optional, default: ``False``);
          whether to issue a warning if duplicate hyperplanes were
          passed -- note that duplicate hyperplanes are always removed,
          whether or not there is a warning shown


        EXAMPLES::

        #     sage: L.<x, y> = HyperplaneArrangements(QQ)
        #     sage: L._element_constructor_(x, y)
        #     Arrangement <y | x>
        #     sage: L._element_constructor_([x, y])
        #     Arrangement <y | x>
        #     sage: L._element_constructor_([0, 1, 0], [0, 0, 1])
        #     Arrangement <y | x>
        #     sage: L._element_constructor_([[0, 1, 0], [0, 0, 1]])
        #     Arrangement <y | x>
        #
        #     sage: L._element_constructor_(polytopes.hypercube(2))
        #     Arrangement <-x + 1 | -y + 1 | y + 1 | x + 1>
        #
        #     sage: L(x, x, warn_duplicates=True)
        #     doctest:...: UserWarning: Input contained 2 hyperplanes, but only 1 are distinct.
        #     Arrangement <x>
        #     sage: L(-x, x + y - 1, signed=False)
        #     Arrangement <-x - y + 1 | x>
        #
        # TESTS::
        #
        #     sage: L()
        #     Empty hyperplane arrangement of dimension 2
        #     sage: L(0)        # zero is equivalent to no argument, Issue #8648
        #     Empty hyperplane arrangement of dimension 2
        #     sage: L(0*x)      # degenerate hyperplane is NOT allowed
        #     Traceback (most recent call last):
        #     ...
        #     ValueError: linear expression must be non-constant to define a hyperplane
        #     sage: L(0*x, y)   # ditto
        #     Traceback (most recent call last):
        #     ...
        #     ValueError: linear expression must be non-constant to define a hyperplane
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

    @cached_method
    def ngens(self):
        """
        Return the number of linear variables.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space
             over Rational Field with coordinates x, y, z
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

            # sage: L = HyperplaneArrangements(QQ, ('x', 'y', 'z'))
            # sage: L.gens()
            # (Hyperplane x + 0*y + 0*z + 0,
            #  Hyperplane 0*x + y + 0*z + 0,
            #  Hyperplane 0*x + 0*y + z + 0)
        """
        return self.ambient_space().gens()

    def gen(self, i):
        """
        Return the `i`-th coordinate.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        A linear expression.

        EXAMPLES::

            # sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            # Hyperplane arrangements in
            #  3-dimensional linear space over Rational Field with coordinates x, y, z
            # sage: L.gen(0)
            # Hyperplane x + 0*y + 0*z + 0
        """
        return self.gens()[i]

    def _coerce_map_from_(self, P):
        """
        Return whether there is a coercion.

        TESTS::

            sage: L.<x> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.<y> = HyperplaneArrangements(RR);  M
            Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y

            sage: L.coerce_map_from(ZZ)
            Coercion map:
              From: Integer Ring
              To:   Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.coerce_map_from(L)
            Coercion map:
              From: Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
              To:   Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y
            sage: L.coerce_map_from(M)
        """
        if self.ambient_space().has_coerce_map_from(P):
            return True
        if isinstance(P, OrderedAffinePlaneCurveArrangements):
            return self.base_ring().has_coerce_map_from(P.base_ring())
        return super()._coerce_map_from_(P)
