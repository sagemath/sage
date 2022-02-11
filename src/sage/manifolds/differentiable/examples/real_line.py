r"""
The Real Line and Open Intervals

The class :class:`OpenInterval` implement open intervals as 1-dimensional
differentiable manifolds over `\RR`. The derived class :class:`RealLine` is
devoted to `\RR` itself, as the open interval `(-\infty, +\infty)`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee2013]_

"""
#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.rings.infinity import infinity, minus_infinity
from sage.symbolic.ring import SR
from sage.rings.real_mpfr import RR
from sage.typeset.unicode_characters import unicode_mathbbR
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.structure import RealDifferentialStructure
from sage.categories.manifolds import Manifolds

class OpenInterval(DifferentiableManifold):
    r"""
    Open interval as a 1-dimensional differentiable manifold over `\RR`.

    INPUT:

    - ``lower`` -- lower bound of the interval (possibly ``-Infinity``)
    - ``upper`` -- upper bound of the interval (possibly ``+Infinity``)
    - ``ambient_interval`` -- (default: ``None``) another open interval,
      to which the constructed interval is a subset of
    - ``name`` -- (default: ``None``) string; name (symbol) given to
      the interval; if ``None``, the name is constructed from ``lower``
      and ``upper``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      interval; if ``None``, the LaTeX symbol is constructed from ``lower``
      and ``upper`` if ``name`` is ``None``, otherwise, it is set to ``name``
    - ``coordinate`` -- (default: ``None``) string defining the symbol of the
      canonical coordinate set on the interval; if none is provided and
      ``names`` is ``None``, the symbol 't' is used
    - ``names`` -- (default: ``None``) used only when ``coordinate`` is
      ``None``: it must be a single-element tuple containing the canonical
      coordinate symbol (this is guaranteed if the shortcut ``<names>`` is
      used, see examples below)
    - ``start_index`` -- (default: 0) unique value of the index for vectors
      and forms on the interval manifold

    EXAMPLES:

    The interval `(0,\pi)`::

        sage: I = manifolds.OpenInterval(0, pi); I
        Real interval (0, pi)
        sage: latex(I)
        \left(0, \pi\right)

    ``I`` is a 1-dimensional smooth manifold over `\RR`::

        sage: I.category()
        Category of smooth connected manifolds over Real Field with 53 bits of
         precision
        sage: I.base_field()
        Real Field with 53 bits of precision
        sage: dim(I)
        1

    It is infinitely differentiable (smooth manifold)::

        sage: I.diff_degree()
        +Infinity

    The instance is unique (as long as the constructor arguments
    are the same)::

        sage: I is manifolds.OpenInterval(0, pi)
        True
        sage: I is manifolds.OpenInterval(0, pi, name='I')
        False

    The display of the interval can be customized::

        sage: I  # default display
        Real interval (0, pi)
        sage: latex(I)  # default LaTeX display
        \left(0, \pi\right)
        sage: I1 = manifolds.OpenInterval(0, pi, name='I'); I1
        Real interval I
        sage: latex(I1)
        I
        sage: I2 = manifolds.OpenInterval(0, pi, name='I', latex_name=r'\mathcal{I}'); I2
        Real interval I
        sage: latex(I2)
        \mathcal{I}

    ``I`` is endowed with a canonical chart::

        sage: I.canonical_chart()
        Chart ((0, pi), (t,))
        sage: I.canonical_chart() is I.default_chart()
        True
        sage: I.atlas()
        [Chart ((0, pi), (t,))]

    The canonical coordinate is returned by the method
    :meth:`canonical_coordinate`::

        sage: I.canonical_coordinate()
        t
        sage: t = I.canonical_coordinate()
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

    However, it can be obtained in the same step as the interval construction
    by means of the shortcut ``I.<names>``::

        sage: I.<t> = manifolds.OpenInterval(0, pi)
        sage: t
        t
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

    The trick is performed by the Sage preparser::

        sage: preparse("I.<t> = manifolds.OpenInterval(0, pi)")
        "I = manifolds.OpenInterval(Integer(0), pi, names=('t',)); (t,) = I._first_ngens(1)"

    In particular the shortcut can be used to set a canonical
    coordinate symbol different from ``'t'``::

        sage: J.<x> = manifolds.OpenInterval(0, pi)
        sage: J.canonical_chart()
        Chart ((0, pi), (x,))
        sage: J.canonical_coordinate()
        x

    The LaTeX symbol of the canonical coordinate can be adjusted via
    the same syntax as a chart declaration (see
    :class:`~sage.manifolds.chart.RealChart`)::

        sage: J.<x> = manifolds.OpenInterval(0, pi, coordinate=r'x:\xi')
        sage: latex(x)
        {\xi}
        sage: latex(J.canonical_chart())
        \left(\left(0, \pi\right),({\xi})\right)

    An element of the open interval ``I``::

        sage: x = I.an_element(); x
        Point on the Real interval (0, pi)
        sage: x.coord() # coordinates in the default chart = canonical chart
        (1/2*pi,)

    As for any manifold subset, a specific element of ``I`` can be created
    by providing a tuple containing its coordinate(s) in a given chart::

        sage: x = I((2,)) # (2,) = tuple of coordinates in the canonical chart
        sage: x
        Point on the Real interval (0, pi)

    But for convenience, it can also be created directly from the coordinate::

        sage: x = I(2); x
        Point on the Real interval (0, pi)
        sage: x.coord()
        (2,)
        sage: I(2) == I((2,))
        True

    By default, the coordinates passed for the element ``x`` are those
    relative to the canonical chart::

        sage: I(2) ==  I((2,), chart=I.canonical_chart())
        True

    The lower and upper bounds of the interval ``I``::

        sage: I.lower_bound()
        0
        sage: I.upper_bound()
        pi

    One of the endpoint can be infinite::

        sage: J = manifolds.OpenInterval(1, +oo); J
        Real interval (1, +Infinity)
        sage: J.an_element().coord()
        (2,)

    The construction of a subinterval can be performed via the argument
    ``ambient_interval`` of ``OpenInterval``::

        sage: J = manifolds.OpenInterval(0, 1, ambient_interval=I); J
        Real interval (0, 1)

    However, it is recommended to use the method :meth:`open_interval`
    instead::

        sage: J = I.open_interval(0, 1); J
        Real interval (0, 1)
        sage: J.is_subset(I)
        True
        sage: J.manifold() is I
        True

    A subinterval of a subinterval::

        sage: K = J.open_interval(1/2, 1); K
        Real interval (1/2, 1)
        sage: K.is_subset(J)
        True
        sage: K.is_subset(I)
        True
        sage: K.manifold() is I
        True

    We have::

        sage: list(I.subset_family())
        [Real interval (0, 1), Real interval (0, pi), Real interval (1/2, 1)]
        sage: list(J.subset_family())
        [Real interval (0, 1), Real interval (1/2, 1)]
        sage: list(K.subset_family())
        [Real interval (1/2, 1)]

    As any open subset of a manifold, open subintervals are created in a
    category of subobjects of smooth manifolds::

        sage: J.category()
        Join of Category of subobjects of sets and Category of smooth manifolds
         over Real Field with 53 bits of precision and Category of connected
         manifolds over Real Field with 53 bits of precision
        sage: K.category()
        Join of Category of subobjects of sets and Category of smooth manifolds
         over Real Field with 53 bits of precision and Category of connected
         manifolds over Real Field with 53 bits of precision

    On the contrary, ``I``, which has not been created as a subinterval,
    is in the category of smooth manifolds (see
    :class:`~sage.categories.manifolds.Manifolds`)::

        sage: I.category()
        Category of smooth connected manifolds over Real Field with 53 bits of
         precision

    and we have::

        sage: J.category() is I.category().Subobjects()
        True

    All intervals are parents::

        sage: x = J(1/2); x
        Point on the Real interval (0, pi)
        sage: x.parent() is J
        True
        sage: y = K(3/4); y
        Point on the Real interval (0, pi)
        sage: y.parent() is K
        True

    We have::

        sage: x in I, x in J, x in K
        (True, True, False)
        sage: y in I, y in J, y in K
        (True, True, True)

    The canonical chart of subintervals is inherited from the canonical chart
    of the parent interval::

        sage: XI = I.canonical_chart(); XI
        Chart ((0, pi), (t,))
        sage: XI.coord_range()
        t: (0, pi)
        sage: XJ = J.canonical_chart(); XJ
        Chart ((0, 1), (t,))
        sage: XJ.coord_range()
        t: (0, 1)
        sage: XK = K.canonical_chart(); XK
        Chart ((1/2, 1), (t,))
        sage: XK.coord_range()
        t: (1/2, 1)

    """
    @staticmethod
    def __classcall_private__(cls, lower, upper, ambient_interval=None,
                              name=None, latex_name=None, coordinate=None,
                              names=None, start_index=0):
        r"""
        Determine the correct interval to return based upon the input.

        TESTS:

        Check whether :trac:`30830` is fixed::

            sage: I = manifolds.OpenInterval(0,2)
            sage: J = manifolds.OpenInterval(0,1, ambient_interval=I, coordinate='t')
            sage: I.open_interval(0,1)
            Real interval (0, 1)

        """
        if ambient_interval:
            # cope the UniqueRepresentation framework for subintervals and
            # reset irrelevant information only:
            coordinate = None
            names = None
            start_index = 0
        return super(cls, OpenInterval).__classcall__(cls, lower, upper,
                          ambient_interval=ambient_interval, name=name,
                          latex_name=latex_name, coordinate=coordinate,
                          names=names, start_index=start_index)

    def __init__(self, lower, upper, ambient_interval=None,
                 name=None, latex_name=None,
                 coordinate=None, names=None, start_index=0):
        r"""
        Construct an open interval.

        TESTS::

            sage: I = manifolds.OpenInterval(-1,1); I
            Real interval (-1, 1)
            sage: TestSuite(I).run(skip='_test_elements')  # pickling of elements fails
            sage: J = manifolds.OpenInterval(-oo, 2); J
            Real interval (-Infinity, 2)
            sage: TestSuite(J).run(skip='_test_elements')  # pickling of elements fails

        """
        if latex_name is None:
            if name is None:
                latex_name = r"\left({}, {}\right)".format(latex(lower), latex(upper))
            else:
                latex_name = name
        if name is None:
            name = "({}, {})".format(lower, upper)
        if ambient_interval is None:
            ambient_manifold = None
        else:
            if not isinstance(ambient_interval, OpenInterval):
                raise TypeError("the argument ambient_interval must be an open interval")
            ambient_manifold = ambient_interval.manifold()
        field = 'real'
        structure = RealDifferentialStructure()
        category = Manifolds(RR).Smooth().Connected()
        DifferentiableManifold.__init__(self, 1, name, field, structure,
                                        base_manifold=ambient_manifold,
                                        latex_name=latex_name,
                                        start_index=start_index,
                                        category=category)
        if ambient_interval is None:
            if coordinate is None:
                if names is None:
                    coordinate = 't'
                else:
                    coordinate = names[0]
        else:
            if lower < ambient_interval.lower_bound():
                raise ValueError("the lower bound is smaller than that of "
                                 + "the containing interval")
            if upper > ambient_interval.upper_bound():
                raise ValueError("the upper bound is larger than that of "
                                 + "the containing interval")
            self.declare_subset(ambient_interval)
            ambient_interval._top_subsets.add(self)
        if lower != minus_infinity:
            if upper != infinity:
                restrictions = lambda t: [t > lower, t < upper]
            else:
                restrictions = lambda t: t > lower
        else:
            if upper != infinity:
                restrictions = lambda t: t < upper
            else:
                restrictions = None
        if ambient_interval is None:
            self._canon_chart = self.chart(coordinates=coordinate,
                                           coord_restrictions=restrictions)
        else:
            self._canon_chart = ambient_interval.canonical_chart().restrict(self,
                                                                            restrictions=restrictions)
        self._lower = lower
        self._upper = upper

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: I = manifolds.OpenInterval(-1, pi)
            sage: I
            Real interval (-1, pi)
            sage: I = manifolds.OpenInterval(-1, +oo)
            sage: I
            Real interval (-1, +Infinity)
            sage: I = manifolds.OpenInterval(-oo,0)
            sage: I
            Real interval (-Infinity, 0)

        """
        return "Real interval " + self._name

    def _first_ngens(self, n):
        r"""
        Return the coordinate of the canonical chart.

        This is useful only for the use of Sage preparser.

        INPUT:

        - ``n`` -- the number of coordinates: must be 1

        TESTS::

            sage: I = manifolds.OpenInterval(-1, 1)
            sage: I._first_ngens(1)
            (t,)
            sage: I = manifolds.OpenInterval(-1, 1, coordinate='x')
            sage: I._first_ngens(1)
            (x,)
            sage: I = manifolds.OpenInterval(-1, 1, names=('x',))
            sage: I._first_ngens(1)
            (x,)

        """
        return self._canon_chart[:]

    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct an element of ``self``.

        This is a redefinition of
        :meth:`sage.manifolds.differentiable.DifferentiableManifold._element_constructor_`
        to allow for construction from a number (considered as the canonical
        coordinate).

        INPUT:

        - ``coords`` -- (default: ``None``) either (i) the point coordinates
          (as a single-element tuple or list) in the chart ``chart``, (ii) the
          value of the point coordinate in the chart ``chart``, or (iii)
          another point in the interval
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the interval's default chart
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          point; if none is provided, the LaTeX symbol is set to ``name``
        - ``check_coords`` -- (default: ``True``) determines whether ``coords``
          are valid coordinates for the chart ``chart``; for symbolic
          coordinates, it is recommended to set ``check_coords`` to ``False``

        OUTPUT:

        - :class:`~sage.manifolds.point.TopologicalManifoldPoint`
          representing a point in the current interval

        EXAMPLES::

            sage: I = manifolds.OpenInterval(-1, 4)
            sage: I((2,))  # standard used of TopologicalManifoldSubset._element_constructor_
            Point on the Real interval (-1, 4)
            sage: I(2)  # specific use with a single coordinate
            Point on the Real interval (-1, 4)
            sage: I(2).coord()
            (2,)
            sage: I(2) == I((2,))
            True
            sage: I(pi)
            Point on the Real interval (-1, 4)
            sage: I(pi).coord()
            (pi,)
            sage: I(8)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates (8,) are not valid on the Chart
             ((-1, 4), (t,))

        """
        if coords in SR:
            coords = (coords,)
        return super(OpenInterval, self)._element_constructor_(coords=coords,
                                 chart=chart, name=name, latex_name=latex_name,
                                 check_coords=check_coords)

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of curves in ``other`` with parameter in ``self``.

        INPUT:

        - ``other`` -- a differentiable manifold `M`
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the set of curves `I \to M`,  where `I` is ``self``

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
            for more documentation.

        TESTS::

            sage: I = manifolds.OpenInterval(-1,1)
            sage: M = Manifold(3, 'M')
            sage: H = I._Hom_(M); H
            Set of Morphisms from Real interval (-1, 1) to 3-dimensional
             differentiable manifold M in Category of smooth manifolds over Real
             Field with 53 bits of precision
            sage: H is Hom(I, M)
            True

        """
        from sage.manifolds.differentiable.manifold_homset import DifferentiableCurveSet
        return DifferentiableCurveSet(self, other)

    def canonical_chart(self):
        r"""
        Return the canonical chart defined on ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES:

        Canonical chart on the interval `(0, \pi)`::

            sage: I = manifolds.OpenInterval(0, pi)
            sage: I.canonical_chart()
            Chart ((0, pi), (t,))
            sage: I.canonical_chart().coord_range()
            t: (0, pi)

        The symbol used for the coordinate of the canonical chart is that
        defined during the construction of the interval::

            sage: I.<x> = manifolds.OpenInterval(0, pi)
            sage: I.canonical_chart()
            Chart ((0, pi), (x,))

        """
        return self._canon_chart

    def canonical_coordinate(self):
        r"""
        Return the canonical coordinate defined on the interval.

        OUTPUT:

        - the symbolic variable representing the canonical coordinate

        EXAMPLES:

        Canonical coordinate on the interval `(0, \pi)`::

            sage: I = manifolds.OpenInterval(0, pi)
            sage: I.canonical_coordinate()
            t
            sage: type(I.canonical_coordinate())
            <class 'sage.symbolic.expression.Expression'>
            sage: I.canonical_coordinate().is_real()
            True

        The canonical coordinate is the first (unique) coordinate of the
        canonical chart::

            sage: I.canonical_coordinate() is I.canonical_chart()[0]
            True

        Its default symbol is `t`; but it can be customized during the
        creation of the interval::

            sage: I = manifolds.OpenInterval(0, pi, coordinate='x')
            sage: I.canonical_coordinate()
            x
            sage: I.<x> = manifolds.OpenInterval(0, pi)
            sage: I.canonical_coordinate()
            x

        """
        return self._canon_chart._xx[0]

    def lower_bound(self):
        r"""
        Return the lower bound (infimum) of the interval.

        EXAMPLES::

            sage: I = manifolds.OpenInterval(1/4, 3)
            sage: I.lower_bound()
            1/4
            sage: J = manifolds.OpenInterval(-oo, 2)
            sage: J.lower_bound()
            -Infinity

        An alias of :meth:`lower_bound` is :meth:`inf`::

            sage: I.inf()
            1/4
            sage: J.inf()
            -Infinity

        """
        return self._lower

    inf = lower_bound

    def upper_bound(self):
        r"""
        Return the upper bound (supremum) of the interval.

        EXAMPLES::

            sage: I = manifolds.OpenInterval(1/4, 3)
            sage: I.upper_bound()
            3
            sage: J = manifolds.OpenInterval(1, +oo)
            sage: J.upper_bound()
            +Infinity

        An alias of :meth:`upper_bound` is :meth:`sup`::

            sage: I.sup()
            3
            sage: J.sup()
            +Infinity

        """
        return self._upper

    sup = upper_bound

    def open_interval(self, lower, upper, name=None, latex_name=None):
        r"""
        Define an open subinterval of ``self``.

        INPUT:

        - ``lower`` -- lower bound of the subinterval (possibly ``-Infinity``)
        - ``upper`` -- upper bound of the subinterval (possibly ``+Infinity``)
        - ``name`` -- (default: ``None``) string; name (symbol) given to the
          subinterval; if ``None``, the name is constructed from ``lower`` and
          ``upper``
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the subinterval; if ``None``, the LaTeX symbol is constructed from
          ``lower`` and ``upper`` if ``name`` is ``None``, otherwise, it is set
          to ``name``

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.examples.real_line.OpenInterval`
          representing the open interval (``lower``, ``upper``)

        EXAMPLES:

        The interval `(0, \pi)` as a subinterval of `(-4, 4)`::

            sage: I = manifolds.OpenInterval(-4, 4)
            sage: J = I.open_interval(0, pi); J
            Real interval (0, pi)
            sage: J.is_subset(I)
            True
            sage: list(I.subset_family())
            [Real interval (-4, 4), Real interval (0, pi)]

        ``J`` is considered as an open submanifold of ``I``::

            sage: J.manifold() is I
            True

        The subinterval `(-4, 4)` is ``I`` itself::

            sage: I.open_interval(-4, 4) is I
            True

        """
        if lower == self._lower and upper == self._upper:
            return self
        # To cope with the unique representation framework, we have to
        # distinguish several cases, instead of performing a mere
        # return OpenInterval(lower, upper, ambient_interval=self, name=name,
        #                     latex_name=latex_name)
        if name is None:
            if latex_name is None:
                return OpenInterval(lower, upper, ambient_interval=self)
            return OpenInterval(lower, upper, ambient_interval=self,
                                latex_name=latex_name)
        if latex_name is None:
            return OpenInterval(lower, upper, ambient_interval=self, name=name)
        return OpenInterval(lower, upper, ambient_interval=self, name=name,
                            latex_name=latex_name)


#******************************************************************************

class RealLine(OpenInterval):
    r"""
    Field of real numbers, as a differentiable manifold of dimension 1 (real
    line) with a canonical coordinate chart.

    INPUT:

    - ``name`` -- (default: ``'R'``) string; name (symbol) given to
      the real line
    - ``latex_name`` -- (default: ``r'\Bold{R}'``) string; LaTeX symbol to
      denote the real line
    - ``coordinate`` -- (default: ``None``) string defining the symbol of the
      canonical coordinate set on the real line; if none is provided and
      ``names`` is ``None``, the symbol 't' is used
    - ``names`` -- (default: ``None``) used only when ``coordinate`` is
      ``None``: it must be a single-element tuple containing the canonical
      coordinate symbol (this is guaranteed if the shortcut ``<names>`` is
      used, see examples below)
    - ``start_index`` -- (default: 0) unique value of the index for vectors
      and forms on the real line manifold

    EXAMPLES:

    Constructing the real line without any argument::

        sage: R = manifolds.RealLine() ; R
        Real number line ℝ
        sage: latex(R)
        \Bold{R}

    ``R`` is a 1-dimensional real smooth manifold::

        sage: R.category()
        Category of smooth connected manifolds over Real Field with 53 bits of
         precision
        sage: isinstance(R, sage.manifolds.differentiable.manifold.DifferentiableManifold)
        True
        sage: dim(R)
        1

    It is endowed with a canonical chart::

        sage: R.canonical_chart()
        Chart (ℝ, (t,))
        sage: R.canonical_chart() is R.default_chart()
        True
        sage: R.atlas()
        [Chart (ℝ, (t,))]

    The instance is unique (as long as the constructor arguments are the
    same)::

        sage: R is manifolds.RealLine()
        True
        sage: R is manifolds.RealLine(latex_name='R')
        False

    The canonical coordinate is returned by the method
    :meth:`~sage.manifolds.differentiable.examples.real_line.OpenInterval.canonical_coordinate`::

        sage: R.canonical_coordinate()
        t
        sage: t = R.canonical_coordinate()
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

    However, it can be obtained in the same step as the real line construction
    by means of the shortcut ``R.<names>``::

        sage: R.<t> = manifolds.RealLine()
        sage: t
        t
        sage: type(t)
        <class 'sage.symbolic.expression.Expression'>

    The trick is performed by Sage preparser::

        sage: preparse("R.<t> = manifolds.RealLine()")
        "R = manifolds.RealLine(names=('t',)); (t,) = R._first_ngens(1)"

    In particular the shortcut is to be used to set a canonical
    coordinate symbol different from 't'::

        sage: R.<x> = manifolds.RealLine()
        sage: R.canonical_chart()
        Chart (ℝ, (x,))
        sage: R.atlas()
        [Chart (ℝ, (x,))]
        sage: R.canonical_coordinate()
        x

    The LaTeX symbol of the canonical coordinate can be adjusted via the same
    syntax as a chart declaration (see
    :class:`~sage.manifolds.chart.RealChart`)::

        sage: R.<x> = manifolds.RealLine(coordinate=r'x:\xi')
        sage: latex(x)
        {\xi}
        sage: latex(R.canonical_chart())
        \left(\Bold{R},({\xi})\right)

    The LaTeX symbol of the real line itself can also be customized::

        sage: R.<x> = manifolds.RealLine(latex_name=r'\mathbb{R}')
        sage: latex(R)
        \mathbb{R}

    Elements of the real line can be constructed directly from a number::

        sage: p = R(2) ; p
        Point on the Real number line ℝ
        sage: p.coord()
        (2,)
        sage: p = R(1.742) ; p
        Point on the Real number line ℝ
        sage: p.coord()
        (1.74200000000000,)

    Symbolic variables can also be used::

        sage: p = R(pi, name='pi') ; p
        Point pi on the Real number line ℝ
        sage: p.coord()
        (pi,)
        sage: a = var('a')
        sage: p = R(a) ; p
        Point on the Real number line ℝ
        sage: p.coord()
        (a,)

    The real line is considered as the open interval `(-\infty, +\infty)`::

        sage: isinstance(R, sage.manifolds.differentiable.examples.real_line.OpenInterval)
        True
        sage: R.lower_bound()
        -Infinity
        sage: R.upper_bound()
        +Infinity

    A real interval can be created from ``R`` means of the method
    :meth:`~sage.manifolds.differentiable.examples.real_line.OpenInterval.open_interval`::

        sage: I = R.open_interval(0, 1); I
        Real interval (0, 1)
        sage: I.manifold()
        Real number line ℝ
        sage: list(R.subset_family())
        [Real interval (0, 1), Real number line ℝ]

    """
    @staticmethod
    def __classcall__(cls, name=unicode_mathbbR, latex_name=r'\Bold{R}',
                      coordinate=None, names=None, start_index=0):
        r"""
        Determine the correct interval to return based upon the input.

        TESTS::

            sage: R = manifolds.RealLine(); R
            Real number line ℝ
            sage: R1 = manifolds.RealLine('ℝ'); R1
            Real number line ℝ
            sage: R is R1
            True

        """
        return super(cls, RealLine).__classcall__(cls, name=name,
                                           latex_name=latex_name,
                                           coordinate=coordinate,
                                           names=names, start_index=start_index)

    def __init__(self, name=unicode_mathbbR, latex_name=r'\Bold{R}',
                 coordinate=None, names=None, start_index=0):
        r"""
        Construct the real line manifold.

        TESTS::

            sage: R = manifolds.RealLine() ; R
            Real number line ℝ
            sage: R.category()
            Category of smooth connected manifolds over Real Field with 53 bits
             of precision
            sage: TestSuite(R).run(skip='_test_elements')  # pickling of elements fails

        """
        OpenInterval.__init__(self, minus_infinity, infinity, name=name,
                              latex_name=latex_name, coordinate=coordinate,
                              names=names, start_index=start_index)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: R = manifolds.RealLine()
            sage: R._repr_()
            'Real number line ℝ'
            sage: R = manifolds.RealLine(name='r')
            sage: R._repr_()
            'Real number line r'

        """
        return "Real number line " + self._name

