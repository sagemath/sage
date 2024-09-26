r"""
Coordinate Charts

The class :class:`Chart` implements coordinate charts on a topological
manifold over a topological field `K`. The subclass :class:`RealChart`
is devoted to the case `K=\RR`, for which the concept of coordinate
range is meaningful.
Moreover, :class:`RealChart` is endowed with some plotting
capabilities (cf. method :meth:`~sage.manifolds.chart.RealChart.plot`).

Transition maps between charts are implemented via the class
:class:`CoordChange`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2015): review tweaks
- Eric Gourgoulhon (2019): periodic coordinates,
  add :meth:`~Chart.calculus_method`

REFERENCES:

- Chap. 2 of [Lee2011]_
- Chap. 1 of [Lee2013]_
"""

# ****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh@umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.symbolic.ring import SR
from sage.rings.infinity import Infinity
from sage.misc.latex import latex
from sage.misc.decorators import options
from sage.manifolds.chart_func import ChartFunctionRing
from sage.manifolds.calculus_method import CalculusMethod
from sage.symbolic.expression import Expression
from sage.ext.fast_callable import fast_callable


class Chart(UniqueRepresentation, SageObject):
    r"""
    Chart on a topological manifold.

    Given a topological manifold `M` of dimension `n` over a topological
    field `K`, a *chart* on `M` is a pair `(U, \varphi)`, where `U` is an
    open subset of `M` and `\varphi : U \rightarrow V \subset K^n` is a
    homeomorphism from `U` to an open subset `V` of `K^n`.

    The components `(x^1, \ldots, x^n)` of `\varphi`, defined by
    `\varphi(p) = (x^1(p), \ldots, x^n(p)) \in K^n` for any point
    `p \in U`, are called the *coordinates* of the chart `(U, \varphi)`.

    INPUT:

    - ``domain`` -- open subset `U` on which the chart is defined (must be
      an instance of :class:`~sage.manifolds.manifold.TopologicalManifold`)
    - ``coordinates`` -- (default: '' (empty string)) single string defining
      the coordinate symbols, with ``' '`` (whitespace) as a separator; each
      item has at most three fields, separated by a colon (``:``):

      1. the coordinate symbol (a letter or a few letters)
      2. (optional) the period of the coordinate if the coordinate is
         periodic; the period field must be written as ``period=T``, where
         ``T`` is the period (see examples below)
      3. (optional) the LaTeX spelling of the coordinate; if not provided the
         coordinate symbol given in the first field will be used

      The order of fields 2 and 3 does not matter and each of them can be
      omitted. If it contains any LaTeX expression, the string ``coordinates``
      must be declared with the prefix 'r' (for "raw") to allow for a proper
      treatment of LaTeX's backslash character (see examples below).
      If no period and no LaTeX spelling are to be set for any coordinate, the
      argument ``coordinates`` can be omitted when the shortcut operator
      ``<,>`` is used to declare the chart (see examples below).
    - ``calc_method`` -- (default: ``None``) string defining the calculus
      method for computations involving coordinates of the chart; must be
      one of

      - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
      - ``'sympy'``: SymPy
      - ``None``: the default of
        :class:`~sage.manifolds.calculus_method.CalculusMethod` will be
        used
    - ``names`` -- (default: ``None``) unused argument, except if
      ``coordinates`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``coord_restrictions`` -- additional restrictions on the coordinates.
      A restriction can be any symbolic equality or inequality involving
      the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
      of the list (or set or frozenset) ``coord_restrictions`` are combined
      with the ``and`` operator; if some restrictions are to be combined with
      the ``or`` operator instead, they have to be passed as a tuple in some
      single item of the list (or set or frozenset) ``coord_restrictions``.
      For example::

        coord_restrictions=[x > y, (x != 0, y != 0), z^2 < x]

      means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
      If the list ``coord_restrictions`` contains only one item, this
      item can be passed as such, i.e. writing ``x > y`` instead
      of the single element list ``[x > y]``.  If the chart variables have
      not been declared as variables yet, ``coord_restrictions`` must
      be ``lambda``-quoted.

    EXAMPLES:

    A chart on a complex 2-dimensional topological manifold::

        sage: M = Manifold(2, 'M', field='complex', structure='topological')
        sage: X = M.chart('x y'); X
        Chart (M, (x, y))
        sage: latex(X)
        \left(M,(x, y)\right)
        sage: type(X)
        <class 'sage.manifolds.chart.Chart'>

    To manipulate the coordinates `(x,y)` as global variables,
    one has to set::

        sage: x,y = X[:]

    However, a shortcut is to use the declarator ``<x,y>`` in the left-hand
    side of the chart declaration (there is then no need to pass the string
    ``'x y'`` to ``chart()``)::

        sage: M = Manifold(2, 'M', field='complex', structure='topological')
        sage: X.<x,y> = M.chart(); X
        Chart (M, (x, y))

    The coordinates are then immediately accessible::

        sage: y
        y
        sage: x is X[0] and y is X[1]
        True

    Note that ``x`` and ``y`` declared in ``<x,y>`` are mere Python variable
    names and do not have to coincide with the coordinate symbols;
    for instance, one may write::

        sage: M = Manifold(2, 'M', field='complex', structure='topological')
        sage: X.<x1,y1> = M.chart('x y'); X
        Chart (M, (x, y))

    Then ``y`` is not known as a global Python variable and the
    coordinate `y` is accessible only through the global variable ``y1``::

        sage: y1
        y
        sage: latex(y1)
        y
        sage: y1 is X[1]
        True

    However, having the name of the Python variable coincide with the
    coordinate symbol is quite convenient; so it is recommended to declare::

        sage: M = Manifold(2, 'M', field='complex', structure='topological')
        sage: X.<x,y> = M.chart()

    In the above example, the chart X covers entirely the manifold ``M``::

        sage: X.domain()
        Complex 2-dimensional topological manifold M

    Of course, one may declare a chart only on an open subset of ``M``::

        sage: U = M.open_subset('U')
        sage: Y.<z1, z2> = U.chart(r'z1:\zeta_1 z2:\zeta_2'); Y
        Chart (U, (z1, z2))
        sage: Y.domain()
        Open subset U of the Complex 2-dimensional topological manifold M

    In the above declaration, we have also specified some LaTeX writing
    of the coordinates different from the text one::

        sage: latex(z1)
        {\zeta_1}

    Note the prefix ``r`` in front of the string ``r'z1:\zeta_1 z2:\zeta_2'``;
    it makes sure that the backslash character is treated as an ordinary
    character, to be passed to the LaTeX interpreter.

    Periodic coordinates are declared through the keyword ``period=`` in the
    coordinate field::

        sage: N = Manifold(2, 'N', field='complex', structure='topological')
        sage: XN.<Z1,Z2> = N.chart('Z1:period=1+2*I Z2')
        sage: XN.periods()
        (2*I + 1, None)

    Coordinates are Sage symbolic variables (see
    :mod:`sage.symbolic.expression`)::

        sage: type(z1)
        <class 'sage.symbolic.expression.Expression'>

    In addition to the Python variable name provided in the operator ``<.,.>``,
    the coordinates are accessible by their indices::

        sage: Y[0], Y[1]
        (z1, z2)

    The index range is that declared during the creation of the manifold. By
    default, it starts at 0, but this can be changed via the parameter
    ``start_index``::

        sage: M1 = Manifold(2, 'M_1', field='complex', structure='topological',
        ....:               start_index=1)
        sage: Z.<u,v> = M1.chart()
        sage: Z[1], Z[2]
        (u, v)

    The full set of coordinates is obtained by means of the slice
    operator ``[:]``::

        sage: Y[:]
        (z1, z2)

    Some partial sets of coordinates::

        sage: Y[:1]
        (z1,)
        sage: Y[1:]
        (z2,)

    Each constructed chart is automatically added to the manifold's user
    atlas::

        sage: M.atlas()
        [Chart (M, (x, y)), Chart (U, (z1, z2))]

    and to the atlas of the chart's domain::

        sage: U.atlas()
        [Chart (U, (z1, z2))]

    Manifold subsets have a *default chart*, which, unless changed via the
    method
    :meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart`,
    is the first defined chart on the subset (or on a open subset of it)::

        sage: M.default_chart()
        Chart (M, (x, y))
        sage: U.default_chart()
        Chart (U, (z1, z2))

    The default charts are not privileged charts on the manifold, but rather
    charts whose name can be skipped in the argument list of functions having
    an optional ``chart=`` argument.

    The chart map `\varphi` acting on a point is obtained by passing
    it as an input to the map::

        sage: p = M.point((1+i, 2), chart=X); p
        Point on the Complex 2-dimensional topological manifold M
        sage: X(p)
        (I + 1, 2)
        sage: X(p) == p.coord(X)
        True

    Setting additional coordinate restrictions::

        sage: M = Manifold(2, 'M', field='complex', structure='topological')
        sage: X.<x,y> = M.chart(coord_restrictions=lambda x,y: abs(x) > 1)
        sage: X.valid_coordinates(2+i, 1)
        True
        sage: X.valid_coordinates(i, 1)
        False

    .. SEEALSO::

        :class:`sage.manifolds.chart.RealChart` for charts on topological
        manifolds over `\RR`.
    """

    @staticmethod
    def __classcall__(cls, domain, coordinates='',
                      calc_method=None, names=None,
                      coord_restrictions=None, **coordinate_options):
        r"""
        Normalize init args and implement unique representation behavior.

        TESTS::

            sage: from sage.manifolds.chart import Chart
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: var("u v")
            (u, v)
            sage: Chart(M, (u, v)) is Chart(M, "u v")
            True
        """
        if isinstance(coordinates, str):
            if coordinates == '':
                for x in names:
                    coordinates += x + ' '
                coordinates = coordinates[:-1]
            coordinates, parsed_options = cls._parse_coordinates(domain,
                                                                 coordinates)
            if not coordinate_options:
                coordinate_options = parsed_options

        coord_string = ' '.join(str(x) for x in coordinates)

        try:
            return domain._charts_by_coord[coord_string]
        except KeyError:
            # Make coord_restrictions hashable
            coord_restrictions = cls._normalize_coord_restrictions(coordinates,
                                                                   coord_restrictions)
            self = super().__classcall__(cls, domain, coordinates, calc_method,
                                         coord_restrictions=coord_restrictions,
                                         **coordinate_options)
            domain._charts_by_coord[coord_string] = self
            return self

    def __init__(self, domain, coordinates, calc_method=None, periods=None,
                 coord_restrictions=None):
        r"""
        Construct a chart.

        TESTS::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X
            Chart (M, (x, y))
            sage: type(X)
            <class 'sage.manifolds.chart.Chart'>
            sage: assumptions() # no assumptions on x,y set
            []
            sage: TestSuite(X).run()

        Check that :issue:`32112` has been fixed::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: XU = U.chart('x y')
            sage: XV = V.chart('x y')
            sage: M.top_charts()
            [Chart (U, (x, y)), Chart (V, (x, y))]
        """
        from sage.manifolds.manifold import TopologicalManifold
        if not isinstance(domain, TopologicalManifold):
            raise TypeError("the first argument must be an open subset of " +
                            "a topological manifold")
        self._manifold = domain.manifold()
        self._domain = domain
        self._sindex = self._manifold.start_index()
        # Handling of calculus methods available on this chart:
        self._calc_method = CalculusMethod(current=calc_method,
                             base_field_type=self.manifold().base_field_type())
        self.simplify = self._calc_method.simplify

        # Treatment of the coordinates:
        self._periods = periods

        if len(coordinates) != self._manifold.dim():
            raise ValueError("the list of coordinates must contain " +
                             "{} elements".format(self._manifold.dim()))
        self._xx = coordinates
        #
        # Additional restrictions on the coordinates.
        self._restrictions = sorted(coord_restrictions, key=str)
        #
        # The chart is added to the domain's atlas, as well as to all the
        # atlases of the domain's supersets; moreover the first defined chart
        # is considered as the default chart
        for sd in domain.open_supersets():
            # the chart is added in the top charts iff its coordinates have
            # not been used on a domain including the chart's domain:
            for chart in sd._atlas:
                if (domain.is_subset(chart._domain)
                    and self._xx == chart._xx):
                    break
            else:
                sd._top_charts.append(self)
            sd._atlas.append(self)
            if sd._def_chart is None:
                sd._def_chart = self
        # The chart is added to the list of the domain's covering charts:
        domain._covering_charts.append(self)
        # Initialization of the set of charts that are restrictions of the
        # current chart to subsets of the chart domain:
        self._subcharts = set([self])
        # Initialization of the set of charts which the current chart is a
        # restriction of:
        self._supercharts = set([self])

        self._dom_restrict = {}  # dict. of the restrictions of self to
                                 # subsets of self._domain, with the
                                 # subsets as keys
        # The null and one functions of the coordinates:
        # Expression in self of the zero and one scalar fields of open sets
        # containing the domain of self:
        for dom in domain.open_supersets():
            dom._zero_scalar_field._express[self] = self.function_ring().zero()
            dom._one_scalar_field._express[self] = self.function_ring().one()

    @classmethod
    def _parse_coordinates(cls, domain, coordinates):
        r"""
        Initialization of the coordinates as symbolic variables.

        INPUT:

        - ``coord_list`` -- list (or space-separated concatenation) of
          coordinate fields.  Each field is a string of at most 3 items,
          separated by ":". These items are: the coordinate symbol, the
          (optional) indicator of the periodic character of the
          coordinate, and the (optional) coordinate LaTeX symbol

        OUTPUT:

        - a tuple of variables (as elements of ``SR``)
        - a dictionary with possible keys:

          - ``'periods'`` -- a tuple of periods

        TESTS::

            sage: from sage.manifolds.chart import Chart
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: Chart._parse_coordinates(M, ['z1', 'z2'])
            ((z1, z2), {'periods': (None, None)})
            sage: Chart._parse_coordinates(M, 'z1 z2')
            ((z1, z2), {'periods': (None, None)})
            sage: Chart._parse_coordinates(M, [r'z1:\zeta_1', r'z2:\zeta_2'])
            ((z1, z2), {'periods': (None, None)})
        """
        if isinstance(coordinates, str):
            coord_list = coordinates.split()
        else:
            coord_list = coordinates
        xx_list = [] # will contain the coordinates as Sage symbolic variables
        period_list = []
        for coord_index, coord_field in enumerate(coord_list):
            coord_properties = coord_field.split(':')
            coord_symb = coord_properties[0].strip() # the coordinate symbol
            coord_latex = None # possibly redefined below
            period = None      # possibly redefined below
            # scan of the properties other than the symbol:
            for prop in coord_properties[1:]:
                prop1 = prop.strip()
                if prop1[0:6] == 'period':
                    if domain.base_field_type() in ['real', 'complex']:
                        period = SR(prop1[7:])
                    else:
                        period = domain.base_field()(prop1[7:])
                else:
                    # prop1 is the coordinate's LaTeX symbol
                    coord_latex = prop1
            # Construction of the coordinate as a Sage symbolic variable:
            coord_var = SR.var(coord_symb, latex_name=coord_latex)
            xx_list.append(coord_var)
            period_list.append(period)
        return tuple(xx_list), dict(periods=tuple(period_list))

    @staticmethod
    def _normalize_coord_restrictions(coordinates, coord_restrictions):
        r"""
        Rewrite ``coord_restrictions`` as a ``frozenset``, representing a logical "and", of other clauses.

        Also replace ``list`` by ``frozenset`` , making the result hashable.

        EXAMPLES::

            sage: from sage.manifolds.chart import Chart
            sage: coordinates = var("x y z")
            sage: Chart._normalize_coord_restrictions(coordinates, None)
            frozenset()
            sage: Chart._normalize_coord_restrictions(coordinates, x > y)
            frozenset({x > y})
            sage: Chart._normalize_coord_restrictions(coordinates, (x != 0, y != 0))
            frozenset({(x != 0, y != 0)})
            sage: Chart._normalize_coord_restrictions(coordinates, [x > y, (x != 0, y != 0), z^2 < x])
            frozenset({(x != 0, y != 0), x > y, z^2 < x})
        """
        def normalize(r):
            if isinstance(r, tuple):  # or
                return tuple(normalize(x) for x in r)
            elif isinstance(r, (list, set, frozenset)):  # and
                return frozenset(normalize(x) for x in r)
            else:
                return r

        if coord_restrictions is None:
            return frozenset()

        if callable(coord_restrictions) and not isinstance(coord_restrictions,
                                                           Expression):
            # lambda-quoted
            coord_restrictions = coord_restrictions(*coordinates)

        if not isinstance(coord_restrictions, (list, set, frozenset)):
            # case of a single condition or conditions to be combined by "or"
            coord_restrictions = [coord_restrictions]

        return normalize(coord_restrictions)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X
            Chart (M, (x, y))
        """
        return 'Chart ({}, {})'.format(self.domain()._name, self._xx)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X._latex_()
            '\\left(M,(x, y)\\right)'
            sage: Y.<z1, z2> = M.chart(r'z1:\zeta_1 z2:\zeta2')
            sage: Y._latex_()
            '\\left(M,({\\zeta_1}, {\\zeta2})\\right)'
            sage: latex(Y)
            \left(M,({\zeta_1}, {\zeta2})\right)
        """
        description = r'\left(' + latex(self.domain()).strip() + ',('
        n = len(self._xx)
        for i in range(n-1):
            description += latex(self._xx[i]).strip() + ', '
        description += latex(self._xx[n-1]).strip() + r')\right)'
        return description

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates.

        This is useful only for the use of Sage preparser::

            sage: preparse("c_cart.<x,y,z> = M.chart()")
            "c_cart = M.chart(names=('x', 'y', 'z',)); (x, y, z,) = c_cart._first_ngens(3)"
        """
        return self[:]

    def __getitem__(self, i):
        r"""
        Access to the coordinates.

        INPUT:

        - ``i`` -- index of the coordinate; if the slice ``[:]``, then all
            the coordinates are returned

        OUTPUT:

        - the coordinate of index ``i`` or all the coordinates (as a tuple)
          if ``i`` is the slice ``[:]``

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The index range is controlled by the parameter ``start_index``::

            sage: M = Manifold(2, 'M', field='complex', structure='topological',
            ....:              start_index=1)
            sage: X.<x,y> = M.chart()
            sage: X[1]
            x
            sage: X[2]
            y
            sage: X[:]
            (x, y)

        We check that slices are properly shifted as well::

            sage: X[2:]
            (y,)
            sage: X[:2]
            (x,)
        """
        if isinstance(i, slice):
            start, stop = i.start, i.stop
            if start is not None:
                start -= self._sindex
            if stop is not None:
                stop -= self._sindex
            return self._xx[start:stop:i.step]
        return self._xx[i-self._sindex]

    def __call__(self, point):
        r"""
        Return the coordinates of a given point.

        INPUT:

        - ``point`` -- point in the domain of the chart

        OUTPUT: tuple of the coordinates of the point

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1+i, 2-i), chart=X)
            sage: X(p)
            (I + 1, -I + 2)
            sage: X(M.an_element())
            (0, 0)
        """
        return point.coord(self)

    def domain(self):
        r"""
        Return the open subset on which the chart is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.domain()
            2-dimensional topological manifold M
            sage: U = M.open_subset('U')
            sage: Y.<u,v> = U.chart()
            sage: Y.domain()
            Open subset U of the 2-dimensional topological manifold M
        """
        return self._domain

    def manifold(self):
        r"""
        Return the manifold on which the chart is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart()
            sage: X.manifold()
            2-dimensional topological manifold M
            sage: X.domain()
            Open subset U of the 2-dimensional topological manifold M
        """
        return self._manifold

    def periods(self):
        r"""
        Return the coordinate periods.

        OUTPUT:

        - a tuple containing the period of each coordinate, with the
          value ``None`` if the coordinate is not periodic

        EXAMPLES:

        A chart without any periodic coordinate::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.periods()
            (None, None)

        Charts with a periodic coordinate::

            sage: Y.<u,v> = M.chart("u v:(0,2*pi):periodic")
            sage: Y.periods()
            (None, 2*pi)
            sage: Z.<a,b> = M.chart(r"a:period=sqrt(2):\alpha b:\beta")
            sage: Z.periods()
            (sqrt(2), None)

        Complex manifold with a periodic coordinate::

            sage: M = Manifold(2, 'M', field='complex', structure='topological',
            ....:              start_index=1)
            sage: X.<x,y> = M.chart("x y:period=1+I")
            sage: X.periods()
            (None, I + 1)

        TESTS::

            sage: M = Manifold(2, 'M', field=QQ, structure='topological')
            sage: X.<xq,yq> = M.chart(r"xq:period=3/2 yq:\zeta:period=2")
            sage: X.periods()
            (3/2, 2)
        """
        return self._periods

    def add_restrictions(self, restrictions):
        r"""
        Add some restrictions on the coordinates.

        This is deprecated; provide the restrictions at the time of creating
        the chart.

        INPUT:

        - ``restrictions`` -- list of restrictions on the
          coordinates, in addition to the ranges declared by the intervals
          specified in the chart constructor

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list ``restrictions`` are combined with the ``and`` operator;
        if some restrictions are to be combined with the ``or`` operator
        instead, they have to be passed as a tuple in some single item
        of the list ``restrictions``. For example::

          restrictions = [x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.add_restrictions(abs(x) > 1)
            doctest:warning...
            DeprecationWarning: Chart.add_restrictions is deprecated; provide the
            restrictions at the time of creating the chart
            See https://github.com/sagemath/sage/issues/32102 for details.
            sage: X.valid_coordinates(2+i, 1)
            True
            sage: X.valid_coordinates(i, 1)
            False
        """
        from sage.misc.superseded import deprecation
        deprecation(32102, "Chart.add_restrictions is deprecated; provide the restrictions at the time of creating the chart")
        self._restrictions.extend(self._normalize_coord_restrictions(self._xx, restrictions))

    def restrict(self, subset, restrictions=None):
        r"""
        Return the restriction of ``self`` to some open subset of its domain.

        If the current chart is `(U,\varphi)`, a *restriction* (or *subchart*)
        is a chart `(V,\psi)` such that `V\subset U` and `\psi = \varphi |_V`.

        If such subchart has not been defined yet, it is constructed here.

        The coordinates of the subchart bare the same names as the coordinates
        of the current chart.

        INPUT:

        - ``subset`` -- open subset `V` of the chart domain `U` (must be an
          instance of :class:`~sage.manifolds.manifold.TopologicalManifold`)
        - ``restrictions`` -- (default: ``None``) list of coordinate
          restrictions defining the subset `V`

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list ``restrictions`` are combined with the ``and`` operator;
        if some restrictions are to be combined with the ``or`` operator
        instead, they have to be passed as a tuple in some single item
        of the list ``restrictions``. For example::

          restrictions = [x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.

        OUTPUT:

        - chart `(V, \psi)` as a :class:`Chart`

        EXAMPLES:

        Coordinates on the unit open ball of  `\CC^2` as a subchart
        of the global coordinates of `\CC^2`::

            sage: M = Manifold(2, 'C^2', field='complex', structure='topological')
            sage: X.<z1, z2> = M.chart()
            sage: B = M.open_subset('B')
            sage: X_B = X.restrict(B, abs(z1)^2 + abs(z2)^2 < 1); X_B
            Chart (B, (z1, z2))
        """
        if subset == self.domain():
            return self
        if subset not in self._dom_restrict:
            if not subset.is_subset(self.domain()):
                raise ValueError("the specified subset is not a subset " +
                                 "of the domain of definition of the chart")
            coordinates = ""
            for coord in self._xx:
                coordinates += repr(coord) + ' '
            res_coord_restrictions = set(self._restrictions)
            res_coord_restrictions.update(self._normalize_coord_restrictions(self._xx,
                                                                             restrictions))
            res = type(self)(subset, coordinates,
                             calc_method=self._calc_method._current,
                             periods=self._periods,
                             # The coordinate restrictions are added
                             # to the result chart
                             coord_restrictions=res_coord_restrictions)
            # Update of supercharts and subcharts:
            res._supercharts.update(self._supercharts)
            for schart in self._supercharts:
                schart._subcharts.add(res)
                schart._dom_restrict[subset] = res
            # Update of domain restrictions:
            self._dom_restrict[subset] = res
        return self._dom_restrict[subset]

    def valid_coordinates(self, *coordinates, **kwds):
        r"""
        Check whether a tuple of coordinates can be the coordinates of a
        point in the chart domain.

        INPUT:

        - ``*coordinates`` -- coordinate values
        - ``**kwds`` -- options:

          - ``parameters=None``, dictionary to set numerical values to
            some parameters (see example below)

        OUTPUT:

        - ``True`` if the coordinate values are admissible in the chart
          image, ``False`` otherwise

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart(coord_restrictions=lambda x,y: [abs(x)<1, y!=0])
            sage: X.valid_coordinates(0, i)
            True
            sage: X.valid_coordinates(i, 1)
            False
            sage: X.valid_coordinates(i/2, 1)
            True
            sage: X.valid_coordinates(i/2, 0)
            False
            sage: X.valid_coordinates(2, 0)
            False

        Example of use with the keyword ``parameters`` to set a specific value
        to a parameter appearing in the coordinate restrictions::

            sage: var('a')  # the parameter is a symbolic variable
            a
            sage: Y.<u,v> = M.chart(coord_restrictions=lambda u,v: abs(v)<a)
            sage: Y.valid_coordinates(1, i, parameters={a: 2})  # setting a=2
            True
            sage: Y.valid_coordinates(1, 2*i, parameters={a: 2})
            False
        """
        if len(coordinates) != self.domain()._dim:
            return False
        if 'parameters' in kwds:
            parameters = kwds['parameters']
        else:
            parameters = None
        # Check of restrictions:
        if self._restrictions:
            substitutions = dict(zip(self._xx, coordinates))
            if parameters:
                substitutions.update(parameters)
            return self._check_restrictions(self._restrictions, substitutions)
        return True

    def _check_restrictions(self, restrict, substitutions):
        r"""
        Recursive helper function to check the validity of coordinates
        given some restrictions

        INPUT:

        - ``restrict`` -- a tuple of conditions (combined with 'or'), a list of
          conditions (combined with 'and') or a single coordinate condition
        - ``substitutions`` -- dictionary (keys: coordinates of ``self``) giving the
          value of each coordinate

        OUTPUT: boolean stating whether the conditions are fulfilled by the
        coordinate values

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X._check_restrictions(x>0, {x: pi, y: 0})
            True
            sage: X._check_restrictions(x>0, {x: -sqrt(2), y: 0})
            False
            sage: X._check_restrictions((x>0, [x<y, y<0]), {x: 1, y: 2})
            True
            sage: X._check_restrictions((x>0, [x<y, y<0]), {x: -1, y: 2})
            False
            sage: X._check_restrictions((x>0, [x<y, y<0]), {x: -1, y: -1/2})
            True
            sage: X._check_restrictions([(x<y, y<0), x>0], {x: 1, y: 2})
            True
            sage: X._check_restrictions([(x<y, y<0), x>0], {x: -1, y: 2})
            False
            sage: X._check_restrictions([(x<y, y<0), x>0], {x: 1, y: -2})
            True
            sage: X._check_restrictions([(x<y, y<0), x>0], {x: 2, y: 1})
            False
        """
        if isinstance(restrict, tuple): # case of 'or' conditions
            return any(self._check_restrictions(cond, substitutions)
                       for cond in restrict)
        elif isinstance(restrict, (list, set, frozenset)): # case of 'and' conditions
            return all(self._check_restrictions(cond, substitutions)
                       for cond in restrict)
        # Case of a single condition:
        return bool(restrict.subs(substitutions))

    def codomain(self):
        r"""
        Return the codomain of ``self`` as a set.

        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.codomain()
            Vector space of dimension 2 over Complex Field with 53 bits of precision
        """
        from sage.modules.free_module import VectorSpace
        ambient = VectorSpace(self.manifold().base_field(), self.manifold().dimension())
        if self._restrictions:
            return self._restrict_set(ambient, self._restrictions)
        else:
            return ambient

    def _restrict_set(self, universe, coord_restrictions):
        """
        Return a set corresponding to coordinate restrictions.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: universe = RR^2
            sage: X._restrict_set(universe, x>0)
            { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x > 0 }
            sage: X._restrict_set(universe, x>0)
            { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x > 0 }
            sage: X._restrict_set(universe, (x>0, [x<y, y<0]))
            Set-theoretic union of
             { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x > 0 } and
             { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x < y, y < 0 }
            sage: X._restrict_set(universe, [(x<y, y<0), x>0])
            Set-theoretic intersection of
             Set-theoretic union of
              { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x < y } and
              { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : y < 0 } and
             { (x, y) ∈ Vector space of dimension 2 over Real Field with 53 bits of precision : x > 0 }
        """
        if isinstance(coord_restrictions, tuple): # case of 'or' conditions
            A = self._restrict_set(universe, coord_restrictions[0])
            if len(coord_restrictions) == 1:
                return A
            else:
                return A.union(self._restrict_set(universe, coord_restrictions[1:]))
        elif isinstance(coord_restrictions, (list, set, frozenset)): # case of 'and' conditions
            A = self._restrict_set(universe, coord_restrictions[0])
            if len(coord_restrictions) == 1:
                return A
            else:
                return A.intersection(self._restrict_set(universe, coord_restrictions[1:]))
        # Case of a single condition:
        from sage.sets.condition_set import ConditionSet
        return ConditionSet(universe, coord_restrictions, vars=self._xx)

    def transition_map(self, other, transformations, intersection_name=None,
                       restrictions1=None, restrictions2=None):
        r"""
        Construct the transition map between the current chart,
        `(U, \varphi)` say, and another one, `(V, \psi)` say.

        If `n` is the manifold's dimension, the *transition map*
        is the map

        .. MATH::

            \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset K^n
            \rightarrow \psi(U\cap V) \subset K^n,

        where `K` is the manifold's base field. In other words, the
        transition map expresses the coordinates `(y^1, \ldots, y^n)` of
        `(V, \psi)` in terms of the coordinates `(x^1, \ldots, x^n)` of
        `(U, \varphi)` on the open subset where the two charts intersect,
        i.e. on `U \cap V`.

        INPUT:

        - ``other`` -- the chart `(V, \psi)`
        - ``transformations`` -- tuple (or list) `(Y_1, \ldots, Y_n)`, where
          `Y_i` is the symbolic expression of the coordinate `y^i` in terms
          of the coordinates `(x^1, \ldots, x^n)`
        - ``intersection_name`` -- (default: ``None``) name to be given to the
          subset `U \cap V` if the latter differs from `U` or `V`
        - ``restrictions1`` -- (default: ``None``) list of conditions on the
          coordinates of the current chart that define `U \cap V` if the
          latter differs from `U`
        - ``restrictions2`` -- (default: ``None``) list of conditions on the
          coordinates of the chart `(V,\psi)` that define `U \cap V` if the
          latter differs from `V`

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list ``restrictions`` are combined with the ``and`` operator;
        if some restrictions are to be combined with the ``or`` operator
        instead, they have to be passed as a tuple in some single item
        of the list ``restrictions``. For example::

          restrictions = [x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.

        OUTPUT:

        - the transition map `\psi \circ \varphi^{-1}` defined on
          `U \cap V` as a :class:`CoordChange`

        EXAMPLES:

        Transition map between two stereographic charts on the circle `S^1`::

            sage: M = Manifold(1, 'S^1', structure='topological')
            sage: U = M.open_subset('U') # Complement of the North pole
            sage: cU.<x> = U.chart() # Stereographic chart from the North pole
            sage: V = M.open_subset('V') # Complement of the South pole
            sage: cV.<y> = V.chart() # Stereographic chart from the South pole
            sage: M.declare_union(U,V)   # S^1 is the union of U and V
            sage: trans = cU.transition_map(cV, 1/x, intersection_name='W',
            ....:                           restrictions1= x!=0, restrictions2 = y!=0)
            sage: trans
            Change of coordinates from Chart (W, (x,)) to Chart (W, (y,))
            sage: trans.display()
            y = 1/x

        The subset `W`, intersection of `U` and `V`, has been created by
        ``transition_map()``::

            sage: F = M.subset_family(); F
            Set {S^1, U, V, W} of open subsets of the 1-dimensional topological manifold S^1
            sage: W = F['W']
            sage: W is U.intersection(V)
            True
            sage: M.atlas()
            [Chart (U, (x,)), Chart (V, (y,)), Chart (W, (x,)), Chart (W, (y,))]

        Transition map between the spherical chart and the Cartesian
        one on `\RR^2`::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart()
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x >= 0}
            sage: c_spher.<r,phi> = U.chart(r'r:(0,+oo) phi:(0,2*pi):\phi')
            sage: trans = c_spher.transition_map(c_cart, (r*cos(phi), r*sin(phi)),
            ....:                                restrictions2=(y!=0, x<0))
            sage: trans
            Change of coordinates from Chart (U, (r, phi)) to Chart (U, (x, y))
            sage: trans.display()
            x = r*cos(phi)
            y = r*sin(phi)

        In this case, no new subset has been created since `U \cap M = U`::

            sage: M.subset_family()
            Set {R^2, U} of open subsets of the 2-dimensional topological manifold R^2

        but a new chart has been created: `(U, (x, y))`::

            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (r, phi)), Chart (U, (x, y))]
        """
        dom1 = self.domain()
        dom2 = other.domain()
        dom = dom1.intersection(dom2, name=intersection_name)
        if dom is dom1:
            chart1 = self
        else:
            chart1 = self.restrict(dom, restrictions1)
        if dom is dom2:
            chart2 = other
        else:
            chart2 = other.restrict(dom, restrictions2)
        if not isinstance(transformations, (tuple, list)):
                transformations = [transformations]
        return CoordChange(chart1, chart2, *transformations)

    def preimage(self, codomain_subset, name=None, latex_name=None):
        """
        Return the preimage (pullback) of ``codomain_subset`` under ``self``.

        It is the subset of the domain of ``self`` formed by the points
        whose coordinate vectors lie in ``codomain_subset``.

        INPUT:

        - ``codomain_subset`` -- an instance of
          :class:`~sage.geometry.convex_set.ConvexSet_base` or another
          object with a ``__contains__`` method that accepts coordinate
          vectors
        - ``name`` -- string; name (symbol) given to the subset
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to
          denote the subset; if none are provided, it is set to ``name``

        OUTPUT:

        - either a :class:`~sage.manifolds.manifold.TopologicalManifold` or
          a :class:`~sage.manifolds.subsets.pullback.ManifoldSubsetPullback`

        EXAMPLES::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

        Pulling back a polytope under a chart::

            sage: # needs sage.geometry.polyhedron
            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [2, 1]]); P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: McP = c_cart.preimage(P); McP
            Subset x_y_inv_P of the 2-dimensional topological manifold R^2
            sage: M((1, 2)) in McP
            True
            sage: M((2, 0)) in McP
            False

        Pulling back the interior of a polytope under a chart::

            sage: # needs sage.geometry.polyhedron
            sage: int_P = P.interior(); int_P
            Relative interior of
             a 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: McInt_P = c_cart.preimage(int_P, name='McInt_P'); McInt_P
            Open subset McInt_P of the 2-dimensional topological manifold R^2
            sage: M((0, 0)) in McInt_P
            False
            sage: M((1, 1)) in McInt_P
            True

        Pulling back a point lattice::

            sage: W = span([[1, 0], [3, 5]], ZZ); W
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 5]
            sage: McW = c_cart.pullback(W, name='McW'); McW
            Subset McW of the 2-dimensional topological manifold R^2
            sage: M((4, 5)) in McW
            True
            sage: M((4, 4)) in McW
            False

        Pulling back a real vector subspaces::

            sage: V = span([[1, 2]], RR); V
            Vector space of degree 2 and dimension 1 over Real Field with 53 bits of precision
            Basis matrix:
            [1.00000000000000 2.00000000000000]
            sage: McV = c_cart.pullback(V, name='McV'); McV
            Subset McV of the 2-dimensional topological manifold R^2
            sage: M((2, 4)) in McV
            True
            sage: M((1, 0)) in McV
            False

        Pulling back a finite set of points::

            sage: F = Family([vector(QQ, [1, 2], immutable=True),
            ....:             vector(QQ, [2, 3], immutable=True)])
            sage: McF = c_cart.pullback(F, name='McF'); McF
            Subset McF of the 2-dimensional topological manifold R^2
            sage: M((2, 3)) in McF
            True
            sage: M((0, 0)) in McF
            False

        Pulling back the integers::

            sage: R = manifolds.RealLine(); R
            Real number line ℝ
            sage: McZ = R.canonical_chart().pullback(ZZ, name='ℤ'); McZ
            Subset ℤ of the Real number line ℝ
            sage: R((3/2,)) in McZ
            False
            sage: R((-2,)) in McZ
            True
        """
        from .subsets.pullback import ManifoldSubsetPullback
        return ManifoldSubsetPullback(self, codomain_subset,
                                      name=name, latex_name=latex_name)

    pullback = preimage

    def function_ring(self):
        """
        Return the ring of coordinate functions on ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.function_ring()
            Ring of chart functions on Chart (M, (x, y))
        """

        return ChartFunctionRing(self)

    def function(self, expression, calc_method=None, expansion_symbol=None,
                 order=None):
        r"""
        Define a coordinate function to the base field.

        If the current chart belongs to the atlas of a `n`-dimensional manifold
        over a topological field `K`, a *coordinate function* is a map

        .. MATH::

            \begin{array}{cccc}
                f:&  V\subset K^n & \longrightarrow & K \\
                  &  (x^1,\ldots, x^n) & \longmapsto & f(x^1,\ldots, x^n),
            \end{array}

        where `V` is the chart codomain and `(x^1, \ldots, x^n)` are the
        chart coordinates.

        INPUT:

        - ``expression`` -- a symbolic expression involving the chart
          coordinates, to represent `f(x^1,\ldots, x^n)`

        - ``calc_method`` -- string (default: ``None``); the calculus method
          with respect to which the internal expression of the function must be
          initialized from ``expression``; one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the chart current calculus method is assumed

        - ``expansion_symbol`` -- (default: ``None``) symbolic variable (the
          "small parameter") with respect to which the coordinate expression is
          expanded in power series (around the zero value of this variable)

        - ``order`` -- integer (default: ``None``); the order of the expansion
          if ``expansion_symbol`` is not ``None``; the *order* is defined as
          the degree of the polynomial representing the truncated power series
          in ``expansion_symbol``.

          .. WARNING::

             The value of ``order`` is `n-1`, where `n` is the order of the
             big `O` in the power series expansion

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.chart_func.ChartFunction`
          representing the coordinate function `f`

        EXAMPLES:

        A symbolic coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(sin(x*y))
            sage: f
            sin(x*y)
            sage: type(f)
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
            sage: f.display()
            (x, y) ↦ sin(x*y)
            sage: f(2,3)
            sin(6)

        Using SymPy for the internal representation of the function (dictionary
        ``_express``)::

            sage: g = X.function(x^2 + x*cos(y), calc_method='sympy')
            sage: g._express
            {'sympy': x**2 + x*cos(y)}

        On the contrary, for ``f``, only the ``SR`` part has been initialized::

            sage: f._express
            {'SR': sin(x*y)}

        See :class:`~sage.manifolds.chart_func.ChartFunction` for more examples.
        """
        parent = self.function_ring()
        return parent.element_class(parent, expression, calc_method=calc_method,
                                    expansion_symbol=expansion_symbol,
                                    order=order)

    def zero_function(self):
        r"""
        Return the zero function of the coordinates.

        If the current chart belongs to the atlas of a `n`-dimensional manifold
        over a topological field `K`, the zero coordinate function is the map

        .. MATH::

            \begin{array}{cccc}
                f:&  V\subset K^n & \longrightarrow & K \\
                  &  (x^1,\ldots, x^n) & \longmapsto & 0,
            \end{array}

        where `V` is the chart codomain.

        See class :class:`~sage.manifolds.chart_func.ChartFunction`
        for a complete documentation.

        OUTPUT:

        - a :class:`~sage.manifolds.chart_func.ChartFunction`
          representing the zero coordinate function `f`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.zero_function()
            0
            sage: X.zero_function().display()
            (x, y) ↦ 0
            sage: type(X.zero_function())
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>

        The result is cached::

            sage: X.zero_function() is X.zero_function()
            True

        Zero function on a `p`-adic manifold::

            sage: # needs sage.rings.padics
            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5)); M
            2-dimensional topological manifold M over the 5-adic Field with
             capped relative precision 20
            sage: X.<x,y> = M.chart()
            sage: X.zero_function()
            0
            sage: X.zero_function().display()
            (x, y) ↦ 0
        """
        return self.function_ring().zero()

    def one_function(self):
        r"""
        Return the constant function of the coordinates equal to one.

        If the current chart belongs to the atlas of a `n`-dimensional manifold
        over a topological field `K`, the "one" coordinate function is the map

        .. MATH::

            \begin{array}{cccc}
                f:&  V\subset K^n & \longrightarrow & K \\
                  &  (x^1,\ldots, x^n) & \longmapsto & 1,
            \end{array}

        where `V` is the chart codomain.

        See class :class:`~sage.manifolds.chart_func.ChartFunction`
        for a complete documentation.

        OUTPUT:

        - a :class:`~sage.manifolds.chart_func.ChartFunction`
          representing the one coordinate function `f`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.one_function()
            1
            sage: X.one_function().display()
            (x, y) ↦ 1
            sage: type(X.one_function())
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>

        The result is cached::

            sage: X.one_function() is X.one_function()
            True

        One function on a `p`-adic manifold::

            sage: # needs sage.rings.padics
            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5)); M
            2-dimensional topological manifold M over the 5-adic Field with
             capped relative precision 20
            sage: X.<x,y> = M.chart()
            sage: X.one_function()
            1 + O(5^20)
            sage: X.one_function().display()
            (x, y) ↦ 1 + O(5^20)
        """
        return self.function_ring().one()

    def calculus_method(self):
        r"""
        Return the interface governing the calculus engine for expressions
        involving coordinates of this chart.

        The calculus engine can be one of the following:

        - Sage's symbolic engine (Pynac + Maxima), implemented via the
          Symbolic Ring ``SR``
        - SymPy

        .. SEEALSO::

            :class:`~sage.manifolds.calculus_method.CalculusMethod` for a
            complete documentation.

        OUTPUT: an instance of :class:`~sage.manifolds.calculus_method.CalculusMethod`

        EXAMPLES:

        The default calculus method relies on Sage's Symbolic Ring::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.calculus_method()
            Available calculus methods (* = current):
             - SR (default) (*)
             - sympy

        Accordingly the method
        :meth:`~sage.manifolds.chart_func.ChartFunction.expr` of a function
        ``f`` defined on the chart ``X`` returns a Sage symbolic expression::

            sage: f = X.function(x^2 + cos(y)*sin(x))
            sage: f.expr()
            x^2 + cos(y)*sin(x)
            sage: type(f.expr())
            <class 'sage.symbolic.expression.Expression'>
            sage: parent(f.expr())
            Symbolic Ring
            sage: f.display()
            (x, y) ↦ x^2 + cos(y)*sin(x)

        Changing to SymPy::

            sage: X.calculus_method().set('sympy')
            sage: f.expr()
            x**2 + sin(x)*cos(y)
            sage: type(f.expr())
            <class 'sympy.core.add.Add'>
            sage: parent(f.expr())
            <class 'sympy.core.add.Add'>
            sage: f.display()
            (x, y) ↦ x**2 + sin(x)*cos(y)

        Back to the Symbolic Ring::

            sage: X.calculus_method().set('SR')
            sage: f.display()
            (x, y) ↦ x^2 + cos(y)*sin(x)
        """
        return self._calc_method

    def multifunction(self, *expressions):
        r"""
        Define a coordinate function to some Cartesian power of the base field.

        If `n` and `m` are two positive integers and `(U, \varphi)` is a
        chart on a topological manifold `M` of dimension `n` over a
        topological field `K`, a *multi-coordinate function* associated
        to `(U,\varphi)` is a map

        .. MATH::

            \begin{array}{llcl}
            f:& V \subset K^n & \longrightarrow & K^m \\
              & (x^1, \ldots, x^n) & \longmapsto & (f_1(x^1, \ldots, x^n),
                \ldots, f_m(x^1, \ldots, x^n)),
            \end{array}

        where `V` is the codomain of `\varphi`. In other words, `f` is a
        `K^m`-valued function of the coordinates associated to the chart
        `(U, \varphi)`.

        See :class:`~sage.manifolds.chart_func.MultiCoordFunction` for a
        complete documentation.

        INPUT:

        - ``expressions`` -- list (or tuple) of `m` elements to construct the
          coordinate functions `f_i` (`1\leq i \leq m`); for
          symbolic coordinate functions, this must be symbolic expressions
          involving the chart coordinates, while for numerical coordinate
          functions, this must be data file names

        OUTPUT:

        - a :class:`~sage.manifolds.chart_func.MultiCoordFunction`
          representing `f`

        EXAMPLES:

        Function of two coordinates with values in `\RR^3`::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.multifunction(x+y, sin(x*y), x^2 + 3*y); f
            Coordinate functions (x + y, sin(x*y), x^2 + 3*y) on the Chart (M, (x, y))
            sage: f(2,3)
            (5, sin(6), 13)

        TESTS::

            sage: type(f)
            <class 'sage.manifolds.chart_func.MultiCoordFunction'>
        """
        from sage.manifolds.chart_func import MultiCoordFunction
        return MultiCoordFunction(self, expressions)


# *****************************************************************************

class RealChart(Chart):
    r"""
    Chart on a topological manifold over `\RR`.

    Given a topological manifold `M` of dimension `n` over `\RR`, a *chart*
    on `M` is a pair `(U,\varphi)`, where `U` is an open subset of `M` and
    `\varphi : U \to V \subset \RR^n` is a homeomorphism from `U` to
    an open subset `V` of `\RR^n`.

    The components `(x^1, \ldots, x^n)` of `\varphi`, defined by
    `\varphi(p) = (x^1(p), \ldots, x^n(p))\in \RR^n` for any point
    `p \in U`, are called the *coordinates* of the chart `(U, \varphi)`.

    INPUT:

    - ``domain`` -- open subset `U` on which the chart is defined
    - ``coordinates`` -- (default: '' (empty string)) single string defining
      the coordinate symbols, with ``' '`` (whitespace) as a separator; each
      item has at most four fields, separated by a colon (``:``):

      1. the coordinate symbol (a letter or a few letters)
      2. (optional) the interval `I` defining the coordinate range: if not
         provided, the coordinate is assumed to span all `\RR`; otherwise
         `I` must be provided in the form ``(a,b)`` (or equivalently
         ``]a,b[``); the bounds ``a`` and ``b`` can be ``+/-Infinity``,
         ``Inf``, ``infinity``, ``inf`` or ``oo``; for *singular*
         coordinates, non-open intervals such as ``[a,b]`` and ``(a,b]``
         (or equivalently ``]a,b]``) are allowed; note that the interval
         declaration must not contain any whitespace
      3. (optional) indicator of the periodic character of the coordinate,
         either as ``period=T``, where ``T`` is the period, or as the keyword
         ``periodic`` (the value of the period is then deduced from the
         interval `I` declared in field 2; see examples below)
      4. (optional) the LaTeX spelling of the coordinate; if not provided the
         coordinate symbol given in the first field will be used

      The order of fields 2 to 4 does not matter and each of them can be
      omitted. If it contains any LaTeX expression, the string ``coordinates``
      must be declared with the prefix 'r' (for "raw") to allow for a proper
      treatment of LaTeX's backslash character (see examples below).
      If interval range, no period and no LaTeX spelling are to be set for any
      coordinate, the argument ``coordinates`` can be omitted when the shortcut
      operator ``<,>`` is used to declare the chart (see examples below).
    - ``calc_method`` -- (default: ``None``) string defining the calculus
      method for computations involving coordinates of the chart; must be
      one of

      - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
      - ``'sympy'``: SymPy
      - ``None``: the default of
        :class:`~sage.manifolds.calculus_method.CalculusMethod` will be
        used
    - ``names`` -- (default: ``None``) unused argument, except if
      ``coordinates`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``coord_restrictions`` -- additional restrictions on the coordinates.
      A restriction can be any symbolic equality or inequality involving
      the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
      of the list (or set or frozenset) ``coord_restrictions`` are combined
      with the ``and`` operator; if some restrictions are to be combined with
      the ``or`` operator instead, they have to be passed as a tuple in some
      single item of the list (or set or frozenset) ``coord_restrictions``.
      For example::

        coord_restrictions=[x > y, (x != 0, y != 0), z^2 < x]

      means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
      If the list ``coord_restrictions`` contains only one item, this
      item can be passed as such, i.e. writing ``x > y`` instead
      of the single element list ``[x > y]``.  If the chart variables have
      not been declared as variables yet, ``coord_restrictions`` must
      be ``lambda``-quoted.

    EXAMPLES:

    Cartesian coordinates on `\RR^3`::

        sage: M = Manifold(3, 'R^3', r'\RR^3', structure='topological',
        ....:              start_index=1)
        sage: c_cart = M.chart('x y z'); c_cart
        Chart (R^3, (x, y, z))
        sage: type(c_cart)
        <class 'sage.manifolds.chart.RealChart'>

    To have the coordinates accessible as global variables, one has to set::

        sage: (x,y,z) = c_cart[:]

    However, a shortcut is to use the declarator ``<x,y,z>`` in the left-hand
    side of the chart declaration (there is then no need to pass the string
    ``'x y z'`` to  ``chart()``)::

        sage: M = Manifold(3, 'R^3', r'\RR^3', structure='topological',
        ....:              start_index=1)
        sage: c_cart.<x,y,z> = M.chart(); c_cart
        Chart (R^3, (x, y, z))

    The coordinates are then immediately accessible::

        sage: y
        y
        sage: y is c_cart[2]
        True

    Note that ``x, y, z`` declared in ``<x,y,z>`` are mere Python variable
    names and do not have to coincide with the coordinate symbols; for
    instance, one may write::

        sage: M = Manifold(3, 'R^3', r'\RR^3', structure='topological',
        ....:              start_index=1)
        sage: c_cart.<x1,y1,z1> = M.chart('x y z'); c_cart
        Chart (R^3, (x, y, z))

    Then ``y`` is not known as a global variable and the coordinate `y`
    is accessible only through the global variable ``y1``::

        sage: y1
        y
        sage: y1 is c_cart[2]
        True

    However, having the name of the Python variable coincide with the
    coordinate symbol is quite convenient; so it is recommended to declare::

        sage: forget()   # for doctests only
        sage: M = Manifold(3, 'R^3', r'\RR^3', structure='topological', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()

    Spherical coordinates on the subset `U` of `\RR^3` that is the
    complement of the half-plane `\{y=0, x \geq 0\}`::

        sage: U = M.open_subset('U')
        sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: c_spher
        Chart (U, (r, th, ph))

    Note the prefix 'r' for the string defining the coordinates in the
    arguments of ``chart``.

    Coordinates are Sage symbolic variables (see
    :mod:`sage.symbolic.expression`)::

        sage: type(th)
        <class 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        {\theta}
        sage: assumptions(th)
        [th is real, th > 0, th < pi]

    Coordinate are also accessible by their indices::

        sage: x1 = c_spher[1]; x2 = c_spher[2]; x3 = c_spher[3]
        sage: [x1, x2, x3]
        [r, th, ph]
        sage: (x1, x2, x3) == (r, th, ph)
        True

    The full set of coordinates is obtained by means of the slice ``[:]``::

        sage: c_cart[:]
        (x, y, z)
        sage: c_spher[:]
        (r, th, ph)

    Let us check that the declared coordinate ranges have been taken into
    account::

        sage: c_cart.coord_range()
        x: (-oo, +oo); y: (-oo, +oo); z: (-oo, +oo)
        sage: c_spher.coord_range()
        r: (0, +oo); th: (0, pi); ph: (0, 2*pi)
        sage: bool(th>0 and th<pi)
        True
        sage: assumptions()  # list all current symbolic assumptions
        [x is real, y is real, z is real, r is real, r > 0, th is real,
         th > 0, th < pi, ph is real, ph > 0, ph < 2*pi]

    The coordinate ranges are used for simplifications::

        sage: simplify(abs(r)) # r has been declared to lie in the interval (0,+oo)
        r
        sage: simplify(abs(x)) # no positive range has been declared for x
        abs(x)

    A coordinate can be declared periodic by adding the keyword ``periodic``
    to its range::

        sage: V = M.open_subset('V')
        sage: c_spher1.<r,th,ph1> = \
        ....: V.chart(r'r:(0,+oo) th:(0,pi):\theta ph1:(0,2*pi):periodic:\phi_1')
        sage: c_spher1.periods()
        (None, None, 2*pi)
        sage: c_spher1.coord_range()
        r: (0, +oo); th: (0, pi); ph1: [0, 2*pi] (periodic)

    It is equivalent to give the period as ``period=2*pi``, skipping the
    coordinate range::

        sage: c_spher2.<r,th,ph2> = \
        ....: V.chart(r'r:(0,+oo) th:(0,pi):\theta ph2:period=2*pi:\phi_2')
        sage: c_spher2.periods()
        (None, None, 2*pi)
        sage: c_spher2.coord_range()
        r: (0, +oo); th: (0, pi); ph2: [0, 2*pi] (periodic)

    Each constructed chart is automatically added to the manifold's
    user atlas::

        sage: M.atlas()
        [Chart (R^3, (x, y, z)), Chart (U, (r, th, ph)),
         Chart (V, (r, th, ph1)), Chart (V, (r, th, ph2))]

    and to the atlas of its domain::

        sage: U.atlas()
        [Chart (U, (r, th, ph))]

    Manifold subsets have a *default chart*, which, unless changed
    via the method
    :meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart`,
    is the first defined chart on the subset (or on a open subset of it)::

        sage: M.default_chart()
        Chart (R^3, (x, y, z))
        sage: U.default_chart()
        Chart (U, (r, th, ph))

    The default charts are not privileged charts on the manifold, but rather
    charts whose name can be skipped in the argument list of functions having
    an optional ``chart=`` argument.

    The chart map `\varphi` acting on a point is obtained by means of the
    call operator, i.e. the operator ``()``::

        sage: p = M.point((1,0,-2)); p
        Point on the 3-dimensional topological manifold R^3
        sage: c_cart(p)
        (1, 0, -2)
        sage: c_cart(p) == p.coord(c_cart)
        True
        sage: q = M.point((2,pi/2,pi/3), chart=c_spher) # point defined by its spherical coordinates
        sage: c_spher(q)
        (2, 1/2*pi, 1/3*pi)
        sage: c_spher(q) == q.coord(c_spher)
        True
        sage: a = U.point((1,pi/2,pi)) # the default coordinates on U are the spherical ones
        sage: c_spher(a)
        (1, 1/2*pi, pi)
        sage: c_spher(a) == a.coord(c_spher)
        True

    Cartesian coordinates on `U` as an example of chart construction with
    coordinate restrictions: since `U` is the complement of the half-plane
    `\{y = 0, x \geq 0\}`, we must have `y \neq 0` or `x < 0` on U.
    Accordingly, we set::

        sage: c_cartU.<x,y,z> = U.chart(coord_restrictions=lambda x,y,z: (y!=0, x<0))
        sage: U.atlas()
        [Chart (U, (r, th, ph)), Chart (U, (x, y, z))]
        sage: M.atlas()
        [Chart (R^3, (x, y, z)), Chart (U, (r, th, ph)),
         Chart (V, (r, th, ph1)), Chart (V, (r, th, ph2)),
         Chart (U, (x, y, z))]
        sage: c_cartU.valid_coordinates(-1,0,2)
        True
        sage: c_cartU.valid_coordinates(1,0,2)
        False
        sage: c_cart.valid_coordinates(1,0,2)
        True

    Note that, as an example, the following would have meant `y \neq 0`
    *and* `x < 0`::

        c_cartU.<x,y,z> = U.chart(coord_restrictions=lambda x,y,z: [y!=0, x<0])

    Chart grids can be drawn in 2D or 3D graphics thanks to the method
    :meth:`plot`.
    """
    def __init__(self, domain, coordinates, calc_method=None, bounds=None,
                 periods=None, coord_restrictions=None):
        r"""
        Construct a chart on a real topological manifold.

        TESTS::

            sage: forget()  # for doctests only
            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X
            Chart (M, (x, y))
            sage: type(X)
            <class 'sage.manifolds.chart.RealChart'>
            sage: assumptions()  # assumptions set in X._init_coordinates
            [x is real, y is real]
            sage: TestSuite(X).run()
        """
        super().__init__(domain, coordinates, calc_method=calc_method,
                         periods=periods, coord_restrictions=coord_restrictions)
        self._bounds = bounds
        self._tighten_bounds()
        self._fast_valid_coordinates = None

    @classmethod
    def _parse_coordinates(cls, domain, coordinates):
        r"""
        Initialization of the coordinates as symbolic variables.

        This method must be redefined by derived classes in order to take
        into account specificities (e.g. enforcing real coordinates).

        INPUT:

        - ``coord_list`` -- list (or space-separated concatenation) of
          coordinate fields.  Each field is a string of at most 3 items,
          separated by ":". These items are: the coordinate symbol, the
          (optional) coordinate range or indicator of the periodic
          character of the coordinate, and the (optional) coordinate
          LaTeX symbol

        OUTPUT:

        - a tuple of variables (as elements of ``SR``)
        - a dictionary with possible keys:

          - ``'periods'`` -- a tuple of periods
          - ``'bounds'`` -- a tuple of coordinate ranges

        TESTS::

            sage: from sage.manifolds.chart import RealChart
            sage: M = Manifold(2, 'M', structure='topological')
            sage: RealChart._parse_coordinates(M, ['x', 'y'])
            ((x, y),
            {'bounds': (((-Infinity, False), (+Infinity, False)),
            ((-Infinity, False), (+Infinity, False))),
            'periods': (None, None)})
            sage: RealChart._parse_coordinates(M, [r'x1:\xi:(0,1)', r'y1:\eta'])
            ((x1, y1),
            {'bounds': (((0, False), (1, False)),
            ((-Infinity, False), (+Infinity, False))),
            'periods': (None, None)})
        """
        from sage.symbolic.assumptions import assume
        if isinstance(coordinates, str):
            coord_list = coordinates.split()
        else:
            coord_list = coordinates
        xx_list = [] # will contain the coordinates as Sage symbolic variables
        bounds_list = [] # will contain the coordinate bounds
        period_list = []
        for coord_index, coord_field in enumerate(coord_list):
            coord_properties = coord_field.split(':')
            coord_symb = coord_properties[0].strip() # the coordinate symbol
            # default values, possibly redefined below:
            coord_latex = None
            xmin = -Infinity
            xmin_included = False
            xmax = +Infinity
            xmax_included = False
            period = None
            # scan of the properties other than the symbol:
            is_periodic = False
            for prop in coord_properties[1:]:
                prop1 = prop.strip()
                delim_min = prop1[0]
                if delim_min in ['[', ']', '(']:
                    # prop1 is the coordinate's range
                    xmin_str, xmax_str = prop1[1:len(prop1)-1].split(',')
                    if xmin_str not in ['-inf', '-Inf', '-infinity',
                                        '-Infinity', '-oo']:
                        xmin = SR(xmin_str)
                        xmin_included = ( delim_min == '[' )
                    if xmax_str not in ['inf', '+inf', 'Inf', '+Inf',
                                        'infinity', '+infinity', 'Infinity',
                                        '+Infinity', 'oo', '+oo']:
                        xmax = SR(xmax_str)
                        xmax_included = ( prop1[-1] == ']' )
                elif prop1[0:6] == 'period':
                    # prop1 indicates a periodic coordinate
                    is_periodic = True
                    if prop1[6:8] != 'ic':
                        # case prop1 = 'period=value'
                        xmin = 0
                        xmax = SR(prop1[7:])
                else:
                    # prop1 is the coordinate's LaTeX symbol
                    coord_latex = prop1
            # Construction of the coordinate as a Sage symbolic variable:
            coord_var = SR.var(coord_symb, domain='real',
                               latex_name=coord_latex)
            assume(coord_var, 'real')
            if is_periodic:
                period = xmax - xmin
                xmin_included = 'periodic'
                xmax_included = 'periodic'
            else:
                if not (xmin == -Infinity):
                    if xmin_included:
                        assume(coord_var >= xmin)
                    else:
                        assume(coord_var > xmin)
                if not (xmax == Infinity):
                    if xmax_included:
                        assume(coord_var <= xmax)
                    else:
                        assume(coord_var < xmax)
            xx_list.append(coord_var)
            bounds_list.append(((xmin, xmin_included), (xmax, xmax_included)))
            period_list.append(period)
        return tuple(xx_list), dict(bounds=tuple(bounds_list),
                                    periods=tuple(period_list))

    def coord_bounds(self, i=None):
        r"""
        Return the lower and upper bounds of the range of a coordinate.

        For a nicely formatted output, use :meth:`coord_range` instead.

        INPUT:

        - ``i`` -- (default: ``None``)  index of the coordinate; if ``None``,
          the bounds of all the coordinates are returned

        OUTPUT:

        - the coordinate bounds as the tuple
          ``((xmin, min_included), (xmax, max_included))`` where

          - ``xmin`` is the coordinate lower bound
          - ``min_included`` is a boolean, indicating whether the coordinate
            can take the value ``xmin``, i.e. ``xmin`` is a strict lower
            bound iff ``min_included`` is ``False``
          - ``xmin`` is the coordinate upper bound
          - ``max_included`` is a boolean, indicating whether the coordinate
            can take the value ``xmax``, i.e. ``xmax`` is a strict upper
            bound iff ``max_included`` is ``False``

        EXAMPLES:

        Some coordinate bounds on a 2-dimensional manifold::

            sage: forget()  # for doctests only
            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart('x y:[0,1)')
            sage: c_xy.coord_bounds(0)  # x in (-oo,+oo) (the default)
            ((-Infinity, False), (+Infinity, False))
            sage: c_xy.coord_bounds(1)  # y in [0,1)
            ((0, True), (1, False))
            sage: c_xy.coord_bounds()
            (((-Infinity, False), (+Infinity, False)), ((0, True), (1, False)))
            sage: c_xy.coord_bounds() == (c_xy.coord_bounds(0), c_xy.coord_bounds(1))
            True

        The coordinate bounds can also be recovered via the method
        :meth:`coord_range`::

            sage: c_xy.coord_range()
            x: (-oo, +oo); y: [0, 1)
            sage: c_xy.coord_range(y)
            y: [0, 1)

        or via Sage's function
        :func:`sage.symbolic.assumptions.assumptions`::

            sage: assumptions(x)
            [x is real]
            sage: assumptions(y)
            [y is real, y >= 0, y < 1]
        """
        if i is None:
            return self._bounds
        else:
            return self._bounds[i-self._sindex]

    def codomain(self):
        r"""
        Return the codomain of ``self`` as a set.

        EXAMPLES::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x >= 0}
            sage: c_spher.<r,phi> = U.chart(r'r:(0,+oo) phi:(0,2*pi):\phi')
            sage: c_spher.codomain()
            The Cartesian product of ((0, +oo), (0, 2*pi))

            sage: M = Manifold(3, 'R^3', r'\RR^3', structure='topological', start_index=1)
            sage: c_cart.<x,y,z> = M.chart()
            sage: c_cart.codomain()
            Vector space of dimension 3 over Real Field with 53 bits of precision

        In the current implementation, the codomain of periodic coordinates are represented
        by a fundamental domain::

            sage: V = M.open_subset('V')
            sage: c_spher1.<r,th,ph1> = \
            ....: V.chart(r'r:(0,+oo) th:(0,pi):\theta ph1:(0,2*pi):periodic:\phi_1')
            sage: c_spher1.codomain()
            The Cartesian product of ((0, +oo), (0, pi), [0, 2*pi))
        """
        from sage.sets.real_set import RealSet
        from sage.modules.free_module import VectorSpace
        from sage.categories.cartesian_product import cartesian_product
        intervals = tuple(RealSet.interval(xmin, xmax,
                                           lower_closed=(min_included == 'periodic' or min_included),
                                           upper_closed=(max_included != 'periodic' and max_included))
                          for ((xmin, min_included), (xmax, max_included)) in self._bounds)
        if all(interval.is_universe()
               for interval in intervals):
            ambient = VectorSpace(self.manifold().base_field(), self.manifold().dimension())
        else:
            ambient = cartesian_product(intervals)
        if self._restrictions:
            return self._restrict_set(ambient, self._restrictions)
        else:
            return ambient

    def coord_range(self, xx=None):
        r"""
        Display the range of a coordinate (or all coordinates), as an
        interval.

        INPUT:

        - ``xx`` -- (default: ``None``) symbolic expression corresponding
          to a coordinate of the current chart; if ``None``, the ranges of
          all coordinates are displayed

        EXAMPLES:

        Ranges of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.coord_range()
            x: (-oo, +oo); y: (-oo, +oo)
            sage: X.coord_range(x)
            x: (-oo, +oo)
            sage: U = M.open_subset('U', coord_def={X: [x>1, y<pi]})
            sage: XU = X.restrict(U)  # restriction of chart X to U
            sage: XU.coord_range()
            x: (1, +oo); y: (-oo, pi)
            sage: XU.coord_range(x)
            x: (1, +oo)
            sage: XU.coord_range(y)
            y: (-oo, pi)

        The output is LaTeX-formatted for the notebook::

            sage: latex(XU.coord_range(y))
            y :\ \left( -\infty, \pi \right)
        """
        from sage.tensor.modules.format_utilities import FormattedExpansion

        def _display_coord_range(self, xx, rtxt, rlatex):
            ind = self._xx.index(xx)
            bounds = self._bounds[ind]
            rtxt += "{}: ".format(xx)
            rlatex += latex(xx) + r":\ "
            if bounds[0][1]:
                rtxt += "["
                rlatex += r"\left["
            else:
                rtxt += "("
                rlatex += r"\left("
            xmin = bounds[0][0]
            if xmin == -Infinity:
                rtxt += "-oo, "
                rlatex += r"-\infty,"
            else:
                rtxt += "{}, ".format(xmin)
                rlatex += latex(xmin) + ","
            xmax = bounds[1][0]
            if xmax == Infinity:
                rtxt += "+oo"
                rlatex += r"+\infty"
            else:
                rtxt += "{}".format(xmax)
                rlatex += latex(xmax)
            if bounds[1][1]:
                rtxt += "]"
                rlatex += r"\right]"
                if bounds[1][1] == 'periodic':
                    rtxt += " (periodic)"
                    rlatex += r"\text{(periodic)}"
            else:
                rtxt += ")"
                rlatex += r"\right)"
            return rtxt, rlatex

        resu_txt = ""
        resu_latex = ""
        if xx is None:
            for x in self._xx:
                if resu_txt != "":
                    resu_txt += "; "
                    resu_latex += r";\quad "
                resu_txt, resu_latex = _display_coord_range(self, x, resu_txt,
                                                            resu_latex)
        else:
            resu_txt, resu_latex = _display_coord_range(self, xx, resu_txt,
                                                        resu_latex)
        return FormattedExpansion(resu_txt, resu_latex)

    def add_restrictions(self, restrictions):
        r"""
        Add some restrictions on the coordinates.

        This is deprecated; provide the restrictions at the time of creating
        the chart.

        INPUT:

        - ``restrictions`` -- list of restrictions on the
          coordinates, in addition to the ranges declared by the intervals
          specified in the chart constructor

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list ``restrictions`` are combined with the ``and`` operator;
        if some restrictions are to be combined with the ``or`` operator
        instead, they have to be passed as a tuple in some single item
        of the list ``restrictions``. For example::

          restrictions = [x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.

        EXAMPLES:

        Cartesian coordinates on the open unit disc in `\RR^2`::

            sage: M = Manifold(2, 'M', structure='topological') # the open unit disc
            sage: X.<x,y> = M.chart()
            sage: X.add_restrictions(x^2+y^2<1)
            doctest:warning...
            DeprecationWarning: Chart.add_restrictions is deprecated; provide the
            restrictions at the time of creating the chart
            See https://github.com/sagemath/sage/issues/32102 for details.
            sage: X.valid_coordinates(0,2)
            False
            sage: X.valid_coordinates(0,1/3)
            True

        The restrictions are transmitted to subcharts::

            sage: A = M.open_subset('A') # annulus 1/2 < r < 1
            sage: X_A = X.restrict(A, x^2+y^2 > 1/4)
            sage: X_A._restrictions
            [x^2 + y^2 < 1, x^2 + y^2 > (1/4)]
            sage: X_A.valid_coordinates(0,1/3)
            False
            sage: X_A.valid_coordinates(2/3,1/3)
            True

        If appropriate, the restrictions are transformed into bounds on
        the coordinate ranges::

            sage: U = M.open_subset('U')
            sage: X_U = X.restrict(U)
            sage: X_U.coord_range()
            x: (-oo, +oo); y: (-oo, +oo)
            sage: X_U.add_restrictions([x<0, y>1/2])
            sage: X_U.coord_range()
            x: (-oo, 0); y: (1/2, +oo)
        """
        super().add_restrictions(restrictions)
        self._tighten_bounds()

    def _tighten_bounds(self):
        """
        Update coordinate bounds from the coordinate restrictions.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological') # the open unit disc
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U')
            sage: X_U = X.restrict(U, restrictions=[x<0, y>1/2])
            sage: X_U.coord_range()
            x: (-oo, 0); y: (1/2, +oo)
        """
        import operator
        bounds = list(self._bounds) # convert to a list for modifications
        new_restrictions = []
        for restrict in self._restrictions:
            restrict_used = False # determines whether restrict is used
                                  # to set some coordinate bound
            if not isinstance(restrict, (tuple, list, set, frozenset)): # case of combined
                                                                        # conditions excluded
                operands = restrict.operands()
                left = operands[0]
                right = operands[1]
                right_var = right.variables()
                if left in self._xx:
                    # the l.h.s. of the restriction is a single
                    # coordinate
                    right_coord = [coord for coord in self._xx
                                   if coord in right_var]
                    if not right_coord:
                        # there is no other coordinate in the r.h.s.
                        ind = self._xx.index(left)
                        left_bounds = list(bounds[ind])
                        oper = restrict.operator()
                        oinf = left_bounds[0][0] # old coord inf
                        osup = left_bounds[1][0] # old coord sup
                        if oper == operator.lt:
                            if osup == Infinity or right <= osup:
                                left_bounds[1] = (right, False)
                                restrict_used = True
                        elif oper == operator.le:
                            if osup == Infinity or right < osup:
                                left_bounds[1] = (right, True)
                                restrict_used = True
                        elif oper == operator.gt:
                            if oinf == -Infinity or right >= oinf:
                                left_bounds[0] = (right, False)
                                restrict_used = True
                        elif oper == operator.ge:
                            if oinf == -Infinity or right > oinf:
                                left_bounds[0] = (right, True)
                                restrict_used = True
                        bounds[ind] = tuple(left_bounds)
            if not restrict_used:
                # if restrict has not been used to set a coordinate bound
                # it is maintained in the list of restrictions:
                new_restrictions.append(restrict)
        self._bounds = tuple(bounds)
        self._restrictions = new_restrictions
        self._fast_valid_coordinates = None

    def restrict(self, subset, restrictions=None):
        r"""
        Return the restriction of the chart to some open subset of its domain.

        If the current chart is `(U, \varphi)`, a *restriction* (or *subchart*)
        is a chart `(V, \psi)` such that `V \subset U` and `\psi = \varphi|_V`.

        If such subchart has not been defined yet, it is constructed here.

        The coordinates of the subchart bare the same names as the coordinates
        of the current chart.

        INPUT:

        - ``subset`` -- open subset `V` of the chart domain `U` (must be an
          instance of :class:`~sage.manifolds.manifold.TopologicalManifold`)
        - ``restrictions`` -- (default: ``None``) list of coordinate
          restrictions defining the subset `V`

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list ``restrictions`` are combined with the ``and`` operator;
        if some restrictions are to be combined with the ``or`` operator
        instead, they have to be passed as a tuple in some single item
        of the list ``restrictions``. For example::

          restrictions = [x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.

        OUTPUT:

        - the chart `(V, \psi)` as a :class:`RealChart`

        EXAMPLES:

        Cartesian coordinates on the unit open disc in `\RR^2` as a subchart
        of the global Cartesian coordinates::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: p = M.point((1/2, 0))
            sage: p in D
            True
            sage: q = M.point((1, 2))
            sage: q in D
            False

        Cartesian coordinates on the annulus `1 < \sqrt{x^2 + y^2} < 2`::

            sage: A = M.open_subset('A')
            sage: c_cart_A = c_cart.restrict(A, [x^2+y^2>1, x^2+y^2<4])
            sage: p in A, q in A
            (False, False)
            sage: a = M.point((3/2,0))
            sage: a in A
            True

        TESTS:

        Check that :issue:`32929` is fixed::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart(r"x:(0,+oo) y:(0,2):periodic")
            sage: U = M.open_subset('U', coord_def={X: x<1})
            sage: XU = X.restrict(U)
            sage: XU.coord_range()
            x: (0, 1); y: [0, 2] (periodic)
            sage: XU.periods()
            (None, 2)
        """
        if subset == self.domain():
            return self
        if subset not in self._dom_restrict:
            if not subset.is_subset(self.domain()):
                raise ValueError("the specified subset is not a subset " +
                                 "of the domain of definition of the chart")
            coordinates = ""
            for coord in self._xx:
                coordinates += repr(coord) + ' '
            res_coord_restrictions = set(self._restrictions)
            res_coord_restrictions.update(self._normalize_coord_restrictions(self._xx,
                                                                             restrictions))
            res = type(self)(subset, coordinates,
                             calc_method=self._calc_method._current,
                             bounds=self._bounds, periods=self._periods,
                             # The coordinate restrictions are added
                             # to the result chart and possibly
                             # transformed into coordinate bounds:
                             coord_restrictions=res_coord_restrictions)
            # Update of supercharts and subcharts:
            res._supercharts.update(self._supercharts)
            for schart in self._supercharts:
                schart._subcharts.add(res)
                schart._dom_restrict[subset] = res
            # Update of domain restrictions:
            self._dom_restrict[subset] = res
        return self._dom_restrict[subset]

    def valid_coordinates(self, *coordinates, **kwds):
        r"""
        Check whether a tuple of coordinates can be the coordinates of a
        point in the chart domain.

        INPUT:

        - ``*coordinates`` -- coordinate values
        - ``**kwds`` -- options:

          - ``tolerance=0``, to set the absolute tolerance in the test of
            coordinate ranges
          - ``parameters=None``, to set some numerical values to parameters

        OUTPUT:

        - ``True`` if the coordinate values are admissible in the chart range
          and ``False`` otherwise

        EXAMPLES:

        Cartesian coordinates on a square interior::

            sage: forget()  # for doctest only
            sage: M = Manifold(2, 'M', structure='topological')  # the square interior
            sage: X.<x,y> = M.chart('x:(-2,2) y:(-2,2)')
            sage: X.valid_coordinates(0,1)
            True
            sage: X.valid_coordinates(-3/2,5/4)
            True
            sage: X.valid_coordinates(0,3)
            False

        The unit open disk inside the square::

            sage: D = M.open_subset('D', coord_def={X: x^2+y^2<1})
            sage: XD = X.restrict(D)
            sage: XD.valid_coordinates(0,1)
            False
            sage: XD.valid_coordinates(-3/2,5/4)
            False
            sage: XD.valid_coordinates(-1/2,1/2)
            True
            sage: XD.valid_coordinates(0,0)
            True

        Another open subset of the square, defined by `x^2+y^2<1` or
        (`x>0` and `|y|<1`)::

            sage: B = M.open_subset('B',
            ....:                   coord_def={X: (x^2+y^2<1,
            ....:                                  [x>0, abs(y)<1])})
            sage: XB = X.restrict(B)
            sage: XB.valid_coordinates(-1/2, 0)
            True
            sage: XB.valid_coordinates(-1/2, 3/2)
            False
            sage: XB.valid_coordinates(3/2, 1/2)
            True
        """
        n = len(coordinates)
        if n != self._manifold._dim:
            return False
        if 'tolerance' in kwds:
            tolerance = kwds['tolerance']
        else:
            tolerance = 0
        if 'parameters' in kwds:
            parameters = kwds['parameters']
        else:
            parameters = None
        # Check of the coordinate ranges:
        for x, bounds in zip(coordinates, self._bounds):
            xmin = bounds[0][0] - tolerance
            min_included = bounds[0][1]
            if min_included == 'periodic':
                continue  # no range to check for a periodic coordinate
            xmax = bounds[1][0] + tolerance
            max_included = bounds[1][1]
            if parameters:
                xmin = xmin.subs(parameters)
                xmax = xmax.subs(parameters)
            if min_included:
                if x < xmin:
                    return False
            elif x <= xmin:
                return False
            if max_included:
                if x > xmax:
                    return False
            elif x >= xmax:
                return False
        # Check of additional restrictions:
        if self._restrictions:
            substitutions = dict(zip(self._xx, coordinates))
            if parameters:
                substitutions.update(parameters)
            return self._check_restrictions(self._restrictions, substitutions)
        return True

    def valid_coordinates_numerical(self, *coordinates):
        r"""
        Check whether a tuple of float coordinates can be the coordinates
        of a point in the chart domain.

        This version is optimized for float numbers, and cannot accept
        parameters nor tolerance. The chart restriction must also be
        specified in CNF (i.e. a list of tuples).

        INPUT:

        - ``*coordinates`` -- coordinate values

        OUTPUT:

        - ``True`` if the coordinate values are admissible in the chart
          range and ``False`` otherwise

        EXAMPLES:

        Cartesian coordinates on a square interior::

            sage: forget()  # for doctest only
            sage: M = Manifold(2, 'M', structure='topological')  # the square interior
            sage: X.<x,y> = M.chart('x:(-2,2) y:(-2,2)')
            sage: X.valid_coordinates_numerical(0,1)
            True
            sage: X.valid_coordinates_numerical(-3/2,5/4)
            True
            sage: X.valid_coordinates_numerical(0,3)
            False

        The unit open disk inside the square::

            sage: D = M.open_subset('D', coord_def={X: x^2+y^2<1})
            sage: XD = X.restrict(D)
            sage: XD.valid_coordinates_numerical(0,1)
            False
            sage: XD.valid_coordinates_numerical(-3/2,5/4)
            False
            sage: XD.valid_coordinates_numerical(-1/2,1/2)
            True
            sage: XD.valid_coordinates_numerical(0,0)
            True

        Another open subset of the square, defined by `x^2 + y^2 < 1` or
        (`x > 0` and `|y| < 1`)::

            sage: B = M.open_subset('B',coord_def={X: [(x^2+y^2<1, x>0),
            ....:                   (x^2+y^2<1,  abs(y)<1)]})
            sage: XB = X.restrict(B)
            sage: XB.valid_coordinates_numerical(-1/2, 0)
            True
            sage: XB.valid_coordinates_numerical(-1/2, 3/2)
            False
            sage: XB.valid_coordinates_numerical(3/2, 1/2)
            True
        """
        # case fast callable already computed
        if self._fast_valid_coordinates is not None:
            return self._fast_valid_coordinates(*coordinates)

        # case fast callable has to be computed
        from operator import lt, gt

        if not isinstance(self._restrictions, (list, set, frozenset)):
            if isinstance(self._restrictions, tuple):
                self._restrictions = [self._restrictions]
            elif isinstance(self._restrictions, Expression):
                self._restrictions = [(self._restrictions,)]
            else:
                raise ValueError("restrictions must be in CNF (list of tuples)")

        list_of_clause = []
        for clause in self._restrictions:
            if not isinstance(clause, tuple):
                if isinstance(clause, Expression):
                    clause = (clause,)
                else:
                    raise ValueError("restrictions must be in CNF (list of tuples)")
            list_of_fast_callable = []
            for literal in clause:
                if not isinstance(literal, Expression):
                    raise ValueError("Restrictions must be in CNF (list of tuples)")
                # End of checks

                fl = fast_callable(literal.lhs(), vars=self[:], domain=float)
                fr = fast_callable(literal.rhs(), vars=self[:], domain=float)
                op = literal.operator()
                list_of_fast_callable.append((fl, fr, op))
            list_of_clause.append(list_of_fast_callable)

        # adding bounds as restrictions
        for x, bounds in zip(self[:], self._bounds):
            if bounds[0][1] == 'periodic':
                continue  # no range to check for a periodic coordinate
            xmin = bounds[0][0]
            xmax = bounds[1][0]

            if x <= xmin:
                return False
            if x >= xmax:
                return False

            if xmin is not -Infinity:
                fl = fast_callable(x, vars=self[:], domain=float)
                fr = fast_callable(SR(xmin), vars=self[:], domain=float)
                list_of_clause.append(((fl, fr, gt),))
            if xmax is not Infinity:
                fl = fast_callable(x, vars=self[:], domain=float)
                fr = fast_callable(SR(xmax), vars=self[:], domain=float)
                list_of_clause.append(((fl, fr, lt),))

        # final call
        def evaluate_fast_callable(*coordinates):
            for clause in list_of_clause:
                temp = False
                for fl, fr, op in clause:
                    temp = temp or op(fl(*coordinates), fr(*coordinates))
                if not temp:
                    return False
            return True

        self._fast_valid_coordinates = evaluate_fast_callable
        return self._fast_valid_coordinates(*coordinates)

    @options(max_range=8, color='red', style='-', thickness=1, plot_points=75,
             label_axes=True)
    def plot(self, chart=None, ambient_coords=None, mapping=None,
             fixed_coords=None, ranges=None, number_values=None,
             steps=None, parameters=None, **kwds):
        r"""
        Plot ``self`` as a grid in a Cartesian graph based on
        the coordinates of some ambient chart.

        The grid is formed by curves along which a chart coordinate
        varies, the other coordinates being kept fixed. It is drawn in
        terms of two (2D graphics) or three (3D graphics) coordinates
        of another chart, called hereafter the *ambient chart*.

        The ambient chart is related to the current chart either by
        a transition map if both charts are defined on the same manifold,
        or by the coordinate expression of some continuous map (typically an
        immersion). In the latter case, the two charts may be defined on two
        different manifolds.

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, the ambient chart is set to the current chart
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2
          or 3 coordinates of the ambient chart in terms of which the plot
          is performed; if ``None``, all the coordinates of the ambient
          chart are considered
        - ``mapping`` -- (default: ``None``)
          :class:`~sage.manifolds.continuous_map.ContinuousMap`; continuous
          manifold map providing the link between the current chart and the
          ambient chart (cf. above); if ``None``, both charts are supposed
          to be defined on the same manifold and related by some transition
          map (see :meth:`~sage.manifolds.chart.Chart.transition_map`)
        - ``fixed_coords`` -- (default: ``None``) dictionary with keys the
          chart coordinates that are not drawn and with values the fixed
          value of these coordinates; if ``None``, all the coordinates of the
          current chart are drawn
        - ``ranges`` -- (default: ``None``) dictionary with keys the
          coordinates to be drawn and values tuples ``(x_min, x_max)``
          specifying the coordinate range for the plot; if ``None``, the
          entire coordinate range declared during the chart construction
          is considered (with ``-Infinity`` replaced by ``-max_range``
          and ``+Infinity`` by ``max_range``)
        - ``number_values`` -- (default: ``None``) either an integer or a
          dictionary with keys the coordinates to be drawn and values the
          number of constant values of the coordinate to be considered; if
          ``number_values`` is a single integer, it represents the number of
          constant values for all coordinates; if ``number_values`` is ``None``,
          it is set to 9 for a 2D plot and to 5 for a 3D plot
        - ``steps`` -- (default: ``None``) dictionary with keys the coordinates
          to be drawn and values the step between each constant value of
          the coordinate; if ``None``, the step is computed from the coordinate
          range (specified in ``ranges``) and ``number_values``. On the contrary
          if the step is provided for some coordinate, the corresponding
          number of constant values is deduced from it and the coordinate range.
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the relation between
          the two coordinate systems
        - ``max_range`` -- (default: 8) numerical value substituted to
          +Infinity if the latter is the upper bound of the range of a
          coordinate for which the plot is performed over the entire coordinate
          range (i.e. for which no specific plot range has been set in
          ``ranges``); similarly ``-max_range`` is the numerical valued
          substituted for ``-Infinity``
        - ``color`` -- (default: ``'red'``) either a single color or a
          dictionary of colors, with keys the coordinates to be drawn,
          representing the colors of the lines along which the coordinate
          varies, the other being kept constant; if ``color`` is a single
          color, it is used for all coordinate lines
        - ``style`` -- (default: ``'-'``) either a single line style or
          a dictionary of line styles, with keys the coordinates to be
          drawn, representing the style of the lines along which the
          coordinate varies, the other being kept constant; if ``style``
          is a single style, it is used for all coordinate lines;
          NB: ``style`` is effective only for 2D plots
        - ``thickness`` -- (default: 1) either a single line thickness or a
          dictionary of line thicknesses, with keys the coordinates to be drawn,
          representing the thickness of the lines along which the coordinate
          varies, the other being kept constant; if ``thickness`` is a single
          value, it is used for all coordinate lines
        - ``plot_points`` -- (default: 75) either a single number of points or
          a dictionary of integers, with keys the coordinates to be drawn,
          representing the number of points to plot the lines along which the
          coordinate varies, the other being kept constant; if ``plot_points``
          is a single integer, it is used for all coordinate lines
        - ``label_axes`` -- boolean (default: ``True``); determining whether the
          labels of the ambient coordinate axes shall be added to the graph;
          can be set to ``False`` if the graph is 3D and must be superposed
          with another graph

        OUTPUT:

        - a graphic object, either a :class:`~sage.plot.graphics.Graphics`
          for a 2D plot (i.e. based on 2 coordinates of the ambient chart)
          or a :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot
          (i.e. based on 3 coordinates of the ambient chart)

        EXAMPLES:

        A 2-dimensional chart plotted in terms of itself results in a
        rectangular grid::

            sage: R2 = Manifold(2, 'R^2', structure='topological')  # the Euclidean plane
            sage: c_cart.<x,y> = R2.chart()  # Cartesian coordinates
            sage: g = c_cart.plot(); g  # equivalent to c_cart.plot(c_cart)             # needs sage.plot
            Graphics object consisting of 18 graphics primitives

        .. PLOT::

            R2 = Manifold(2, 'R^2', structure='topological')
            c_cart = R2.chart('x y')
            g = c_cart.plot()
            sphinx_plot(g)

        Grid of polar coordinates in terms of Cartesian coordinates in the
        Euclidean plane::

            sage: U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)})  # the complement of the segment y=0 and x>0
            sage: c_pol.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')  # polar coordinates on U
            sage: pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            sage: g = c_pol.plot(c_cart); g                                             # needs sage.plot
            Graphics object consisting of 18 graphics primitives

        .. PLOT::

            R2 = Manifold(2, 'R^2', structure='topological')
            c_cart = R2.chart('x y'); x, y = c_cart[:]
            U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)})
            c_pol = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi'); r, ph = c_pol[:]
            pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            g = c_pol.plot(c_cart)
            sphinx_plot(g)

        Call with non-default values::

            sage: g = c_pol.plot(c_cart, ranges={ph:(pi/4,pi)},                         # needs sage.plot
            ....:                number_values={r:7, ph:17},
            ....:                color={r:'red', ph:'green'},
            ....:                style={r:'-', ph:'--'})

        .. PLOT::

            R2 = Manifold(2, 'R^2', structure='topological')
            c_cart = R2.chart('x y'); x, y = c_cart[:]
            U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)})
            c_pol = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi'); r, ph = c_pol[:]
            pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            g = c_pol.plot(c_cart, ranges={ph:(pi/4,pi)}, number_values={r:7, ph:17},
                           color={r:'red', ph:'green'}, style={r:'-', ph:'--'})
            sphinx_plot(g)

        A single coordinate line can be drawn::

            sage: g = c_pol.plot(c_cart,    # draw a circle of radius r=2               # needs sage.plot
            ....:                fixed_coords={r: 2})

        .. PLOT::

            R2 = Manifold(2, 'R^2', structure='topological')
            c_cart = R2.chart('x y'); x, y = c_cart[:]
            U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)})
            c_pol = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi'); r, ph = c_pol[:]
            pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            g = c_pol.plot(c_cart, fixed_coords={r: 2})
            sphinx_plot(g)

        ::

            sage: g = c_pol.plot(c_cart,    # draw a segment at phi=pi/4                # needs sage.plot
            ....:                fixed_coords={ph: pi/4})

        .. PLOT::

            R2 = Manifold(2, 'R^2', structure='topological')
            c_cart = R2.chart('x y'); x, y = c_cart[:]
            U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)})
            c_pol = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi'); r, ph = c_pol[:]
            pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            g = c_pol.plot(c_cart, fixed_coords={ph: pi/4})
            sphinx_plot(g)

        An example with the ambient chart lying in an another manifold (the
        plot is then performed via some manifold map passed as the
        argument ``mapping``): 3D plot of the stereographic charts on the
        2-sphere::

            sage: S2 = Manifold(2, 'S^2', structure='topological')  # the 2-sphere
            sage: U = S2.open_subset('U'); V = S2.open_subset('V')  # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart()  # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart()  # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: R3 = Manifold(3, 'R^3', structure='topological')  # the Euclidean space R^3
            sage: c_cart.<X,Y,Z> = R3.chart()  # Cartesian coordinates on R^3
            sage: Phi = S2.continuous_map(R3, {(c_xy, c_cart): [2*x/(1+x^2+y^2),
            ....:                          2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
            ....:                          (c_uv, c_cart): [2*u/(1+u^2+v^2),
            ....:                          2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]},
            ....:                         name='Phi', latex_name=r'\Phi')  # Embedding of S^2 in R^3
            sage: g = c_xy.plot(c_cart, mapping=Phi); g                                 # needs sage.plot
            Graphics3d Object

        .. PLOT::

            S2 = Manifold(2, 'S^2', structure='topological')
            U = S2.open_subset('U') ; V = S2.open_subset('V')
            S2.declare_union(U,V)
            c_xy = U.chart('x y'); x, y = c_xy[:]
            c_uv = V.chart('u v'); u, v = c_uv[:]
            xy_to_uv = c_xy.transition_map(c_uv, (x/(x**2+y**2), y/(x**2+y**2)),
                             intersection_name='W', restrictions1= x**2+y**2!=0,
                             restrictions2= u**2+v**2!=0)
            uv_to_xy = xy_to_uv.inverse()
            R3 = Manifold(3, 'R^3', structure='topological')
            c_cart = R3.chart('X Y Z')
            Phi = S2.continuous_map(R3, {(c_xy, c_cart): [2*x/(1+x**2+y**2),
                              2*y/(1+x**2+y**2), (x**2+y**2-1)/(1+x**2+y**2)],
                              (c_uv, c_cart): [2*u/(1+u**2+v**2),
                              2*v/(1+u**2+v**2), (1-u**2-v**2)/(1+u**2+v**2)]},
                              name='Phi', latex_name=r'\Phi')
            sphinx_plot(c_xy.plot(c_cart, mapping=Phi))

        NB: to get a better coverage of the whole sphere, one should increase
        the coordinate sampling via the argument ``number_values`` or the
        argument ``steps`` (only the default value, ``number_values = 5``, is
        used here, which is pretty low).

        The same plot without the ``(X,Y,Z)`` axes labels::

            sage: g = c_xy.plot(c_cart, mapping=Phi, label_axes=False)                  # needs sage.plot

        The North and South stereographic charts on the same plot::

            sage: g2 = c_uv.plot(c_cart, mapping=Phi, color='green')                    # needs sage.plot
            sage: g + g2                                                                # needs sage.plot
            Graphics3d Object

        .. PLOT::

            S2 = Manifold(2, 'S^2', structure='topological')
            U = S2.open_subset('U') ; V = S2.open_subset('V')
            S2.declare_union(U,V)
            c_xy = U.chart('x y'); x, y = c_xy[:]
            c_uv = V.chart('u v'); u, v = c_uv[:]
            xy_to_uv = c_xy.transition_map(c_uv, (x/(x**2+y**2), y/(x**2+y**2)),
                             intersection_name='W', restrictions1= x**2+y**2!=0,
                             restrictions2= u**2+v**2!=0)
            uv_to_xy = xy_to_uv.inverse()
            R3 = Manifold(3, 'R^3', structure='topological')
            c_cart = R3.chart('X Y Z')
            Phi = S2.continuous_map(R3, {(c_xy, c_cart): [2*x/(1+x**2+y**2),
                              2*y/(1+x**2+y**2), (x**2+y**2-1)/(1+x**2+y**2)],
                              (c_uv, c_cart): [2*u/(1+u**2+v**2),
                              2*v/(1+u**2+v**2), (1-u**2-v**2)/(1+u**2+v**2)]},
                              name='Phi', latex_name=r'\Phi')
            g = c_xy.plot(c_cart, mapping=Phi, label_axes=False)
            g2 = c_uv.plot(c_cart, mapping=Phi, color='green')
            sphinx_plot(g+g2)

        South stereographic chart drawn in terms of the North one (we split
        the plot in four parts to avoid the singularity at `(u,v)=(0,0)`)::

            sage: # long time, needs sage.plot
            sage: W = U.intersection(V) # the subset common to both charts
            sage: c_uvW = c_uv.restrict(W) # chart (W,(u,v))
            sage: gSN1 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[-6.,-0.02]})
            sage: gSN2 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[0.02,6.]})
            sage: gSN3 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[-6.,-0.02]})
            sage: gSN4 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[0.02,6.]})
            sage: show(gSN1+gSN2+gSN3+gSN4, xmin=-1.5, xmax=1.5, ymin=-1.5, ymax=1.5)

        .. PLOT::

            S2 = Manifold(2, 'S^2', structure='topological')
            U = S2.open_subset('U'); V = S2.open_subset('V'); S2.declare_union(U,V)
            c_xy = U.chart('x y'); x, y = c_xy[:]
            c_uv = V.chart('u v'); u, v = c_uv[:]
            xy_to_uv = c_xy.transition_map(c_uv, (x/(x**2+y**2), y/(x**2+y**2)),
                              intersection_name='W', restrictions1= x**2+y**2!=0,
                              restrictions2= u**2+v**2!=0)
            uv_to_xy = xy_to_uv.inverse()
            c_uvW = c_uv.restrict(U.intersection(V))
            gSN1 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[-6.,-0.02]})
            gSN2 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[0.02,6.]})
            gSN3 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[-6.,-0.02]})
            gSN4 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[0.02,6.]})
            g = gSN1+gSN2+gSN3+gSN4; g.set_axes_range(-1.5, 1.5, -1.5, 1.5)
            sphinx_plot(g)

        The coordinate line `u = 1` (red) and the coordinate line `v = 1`
        (green) on the same plot::

            sage: # long time, needs sage.plot
            sage: gu1 = c_uvW.plot(c_xy, fixed_coords={u: 1}, max_range=20,
            ....:                  plot_points=300)
            sage: gv1 = c_uvW.plot(c_xy, fixed_coords={v: 1}, max_range=20,
            ....:                  plot_points=300, color='green')
            sage: gu1 + gv1
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            S2 = Manifold(2, 'S^2', structure='topological')
            U = S2.open_subset('U'); V = S2.open_subset('V'); S2.declare_union(U,V)
            c_xy = U.chart('x y'); x, y = c_xy[:]
            c_uv = V.chart('u v'); u, v = c_uv[:]
            xy_to_uv = c_xy.transition_map(c_uv, (x/(x**2+y**2), y/(x**2+y**2)),
                              intersection_name='W', restrictions1= x**2+y**2!=0,
                              restrictions2= u**2+v**2!=0)
            uv_to_xy = xy_to_uv.inverse()
            c_uvW = c_uv.restrict(U.intersection(V))
            gu1 = c_uvW.plot(c_xy, fixed_coords={u: 1}, max_range=20, plot_points=300)
            gv1 = c_uvW.plot(c_xy, fixed_coords={v: 1}, max_range=20, plot_points=300,
                             color='green')
            sphinx_plot(gu1+gv1)

        Note that we have set ``max_range=20`` to have a wider range for
        the coordinates `u` and `v`, i.e. to have `[-20, 20]` instead of
        the default `[-8, 8]`.

        A 3-dimensional chart plotted in terms of itself results in a 3D
        rectangular grid::

            sage: # long time, needs sage.plot
            sage: g = c_cart.plot()  # equivalent to c_cart.plot(c_cart)
            sage: g
            Graphics3d Object

        .. PLOT::

            R3 = Manifold(3, 'R^3', structure='topological')
            c_cart = R3.chart('X Y Z')
            sphinx_plot(c_cart.plot())

        A 4-dimensional chart plotted in terms of itself (the plot is
        performed for at most 3 coordinates, which must be specified via
        the argument ``ambient_coords``)::

            sage: # needs sage.plot
            sage: M = Manifold(4, 'M', structure='topological')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = X.plot(ambient_coords=(t,x,y))  # the coordinate z is not depicted  # long time
            sage: g                                                                       # long time
            Graphics3d Object

        .. PLOT::

            M = Manifold(4, 'M', structure='topological')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            g = X.plot(ambient_coords=(t,x,y))
            sphinx_plot(g)

        ::

            sage: # needs sage.plot
            sage: g = X.plot(ambient_coords=(t,y))  # the coordinates x and z are not depicted
            sage: g
            Graphics object consisting of 18 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M', structure='topological')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            g = X.plot(ambient_coords=(t,y))
            sphinx_plot(g)

        Note that the default values of some arguments of the method ``plot``
        are stored in the dictionary ``plot.options``::

            sage: X.plot.options  # random (dictionary output)
            {'color': 'red', 'label_axes': True, 'max_range': 8,
             'plot_points': 75, 'style': '-', 'thickness': 1}

        so that they can be adjusted by the user::

            sage: X.plot.options['color'] = 'blue'

        From now on, all chart plots will use blue as the default color.
        To restore the original default options, it suffices to type::

            sage: X.plot.reset()
        """
        from sage.misc.functional import numerical_approx
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from sage.manifolds.continuous_map import ContinuousMap
        from sage.manifolds.utilities import set_axes_labels

        # Extract the kwds options
        max_range = kwds['max_range']
        color = kwds['color']
        style = kwds['style']
        thickness = kwds['thickness']
        plot_points = kwds['plot_points']
        label_axes = kwds['label_axes']

        def _plot_xx_list(xx_list, rem_coords, ranges, steps, number_values):
            r"""
            Helper function to plot the coordinate grid.
            """
            coord = rem_coords[0]
            xmin = ranges[coord][0]
            sx = steps[coord]
            resu = []
            for xx in xx_list:
                xc = xmin
                for i in range(number_values[coord]):
                    nxx = list(xx)
                    nxx[self._xx.index(coord)] = xc
                    resu.append(nxx)
                    xc += sx
            if len(rem_coords) == 1:
                return resu
            else:
                rem_coords.remove(coord)
                return _plot_xx_list(resu, rem_coords, ranges, steps,
                                     number_values)

        if chart is None:
            chart = self
        elif not isinstance(chart, Chart):
            raise TypeError("the argument 'chart' must be a coordinate chart")
        #
        # 1/ Determination of the relation between self and chart
        #    ------------------------------------------------------------
        nc = self._manifold.dimension()
        if chart is self:
            transf = self.multifunction(*(self._xx))
            if nc > 3:
                if ambient_coords is None:
                    raise TypeError("the argument 'ambient_coords' must be provided")
                if len(ambient_coords) > 3:
                    raise ValueError("too many ambient coordinates")
                fixed_coords = {}
                for coord in self._xx:
                    if coord not in ambient_coords:
                        fixed_coords[coord] = 0
        else:
            transf = None # to be the MultiCoordFunction object relating self
                          # to the ambient chart
            if mapping is None:
                if not self.domain().is_subset(chart.domain()):
                    raise ValueError("the domain of {} is not ".format(self) +
                                     "included in that of {}".format(chart))
                coord_changes = chart.domain()._coord_changes
                for chart_pair in coord_changes:
                    if chart_pair == (self, chart):
                        transf = coord_changes[chart_pair]._transf
                        break
                else:
                    # Search for a subchart
                    for chart_pair in coord_changes:
                        for schart in chart._subcharts:
                            if chart_pair == (self, schart):
                                transf = coord_changes[chart_pair]._transf
            else:
                if not isinstance(mapping, ContinuousMap):
                    raise TypeError("the argument 'mapping' must be a "
                                    "continuous manifold map")
                if not self.domain().is_subset(mapping.domain()):
                    raise ValueError("the domain of {} is not ".format(self) +
                                     "included in that of {}".format(mapping))
                if not chart.domain().is_subset(mapping._codomain):
                    raise ValueError("the domain of {} is not ".format(chart) +
                                     "included in the codomain of {}".format(
                                                                      mapping))
                try:
                    transf = mapping.coord_functions(chart1=self, chart2=chart)
                except ValueError:
                    pass
            if transf is None:
                raise ValueError("no relation has been found between " +
                                 "{} and {}".format(self, chart))
        #
        # 2/ Treatment of input parameters
        #    -----------------------------
        if fixed_coords is None:
            coords = self._xx
        else:
            fixed_coord_list = fixed_coords.keys()
            coords = []
            for coord in self._xx:
                if coord not in fixed_coord_list:
                    coords.append(coord)
            coords = tuple(coords)
        if ambient_coords is None:
            ambient_coords = chart._xx
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca != 3:
            raise ValueError("bad number of ambient coordinates: {}".format(nca))
        if ranges is None:
            ranges = {}
        ranges0 = {}
        for coord in coords:
            if coord in ranges:
                ranges0[coord] = (numerical_approx(ranges[coord][0]),
                                  numerical_approx(ranges[coord][1]))
            else:
                bounds = self._bounds[self._xx.index(coord)]
                if bounds[0][0] == -Infinity:
                    xmin = numerical_approx(-max_range)
                elif bounds[0][1]:
                    xmin = numerical_approx(bounds[0][0])
                else:
                    xmin = numerical_approx(bounds[0][0] + 1.e-3)
                if bounds[1][0] == Infinity:
                    xmax = numerical_approx(max_range)
                elif bounds[1][1]:
                    xmax = numerical_approx(bounds[1][0])
                else:
                    xmax = numerical_approx(bounds[1][0] - 1.e-3)
                ranges0[coord] = (xmin, xmax)
        ranges = ranges0
        if number_values is None:
            if nca == 2: # 2D plot
                number_values = 9
            else:   # 3D plot
                number_values = 5
        if not isinstance(number_values, dict):
            number_values0 = {}
            for coord in coords:
                number_values0[coord] = number_values
            number_values = number_values0
        if steps is None:
            steps = {}
        for coord in coords:
            if coord not in steps:
                steps[coord] = ((ranges[coord][1] - ranges[coord][0])
                                / (number_values[coord]-1))
            else:
                from sage.functions.other import floor
                number_values[coord] = 1 + floor((ranges[coord][1] - ranges[coord][0])
                                             / steps[coord])
        if not isinstance(color, dict):
            color0 = {}
            for coord in coords:
                color0[coord] = color
            color = color0
        if not isinstance(style, dict):
            style0 = {}
            for coord in coords:
                style0[coord] = style
            style = style0
        if not isinstance(thickness, dict):
            thickness0 = {}
            for coord in coords:
                thickness0[coord] = thickness
            thickness = thickness0
        if not isinstance(plot_points, dict):
            plot_points0 = {}
            for coord in coords:
                plot_points0[coord] = plot_points
            plot_points = plot_points0
        #
        # 3/ Plots
        #    -----
        xx0 = [0] * nc
        if fixed_coords is not None:
            if len(fixed_coords) != nc - len(coords):
                raise ValueError("bad number of fixed coordinates")
            for fc, val in fixed_coords.items():
                xx0[self._xx.index(fc)] = val
        ind_a = [chart._xx.index(ac) for ac in ambient_coords]
        resu = Graphics()
        for coord in coords:
            color_c, style_c = color[coord], style[coord]
            thickness_c = thickness[coord]
            rem_coords = list(coords)
            rem_coords.remove(coord)
            xx_list = [xx0]
            if len(rem_coords) >= 1:
                xx_list = _plot_xx_list(xx_list, rem_coords, ranges, steps,
                                        number_values)
            xmin, xmax = ranges[coord]
            nbp = plot_points[coord]
            dx = (xmax - xmin) / (nbp-1)
            ind_coord = self._xx.index(coord)
            for xx in xx_list:
                curve = []
                first_invalid = False # initialization
                xc = xmin
                xp = list(xx)
                if parameters is None:
                    for i in range(nbp):
                        xp[ind_coord] = xc
                        if self.valid_coordinates(*xp, tolerance=1e-13):
                            yp = transf(*xp, simplify=False)
                            curve.append( [numerical_approx(yp[j])
                                           for j in ind_a] )
                            first_invalid = True # next invalid point will be
                                                 # the first one
                        else:
                            if first_invalid:
                                # the curve is stopped at previous point and
                                # added to the graph:
                                resu += line(curve, color=color_c,
                                             linestyle=style_c,
                                             thickness=thickness_c)
                                curve = [] # a new curve will start at the
                                           # next valid point
                            first_invalid = False # next invalid point will not
                                                  # be the first one
                        xc += dx
                else:
                    for i in range(nbp):
                        xp[ind_coord] = xc
                        if self.valid_coordinates(*xp, tolerance=1e-13,
                                                  parameters=parameters):
                            yp = transf(*xp, simplify=False)
                            curve.append([numerical_approx(yp[j].substitute(parameters))
                                          for j in ind_a])
                            first_invalid = True # next invalid point will be
                                                 # the first one
                        else:
                            if first_invalid:
                                # the curve is stopped at previous point and
                                # added to the graph:
                                resu += line(curve, color=color_c,
                                             linestyle=style_c,
                                             thickness=thickness_c)
                                curve = [] # a new curve will start at the
                                           # next valid point
                            first_invalid = False # next invalid point will not
                                                  # be the first one
                        xc += dx
                if curve:
                    resu += line(curve, color=color_c,
                                 linestyle=style_c,
                                 thickness=thickness_c)
        if nca == 2:  # 2D graphic
            resu.set_aspect_ratio(1)
            if label_axes:
                # We update the dictionary _extra_kwds (options to be passed
                # to show()), instead of using the method
                # Graphics.axes_labels() since the latter is not robust w.r.t.
                # graph addition
                resu._extra_kwds['axes_labels'] = [r'$'+latex(ac)+r'$'
                                                   for ac in ambient_coords]
        else: # 3D graphic
            resu.aspect_ratio(1)
            if label_axes:
                labels = [str(ac) for ac in ambient_coords]
                resu = set_axes_labels(resu, *labels)
        return resu

# *****************************************************************************

class CoordChange(SageObject):
    r"""
    Transition map between two charts of a topological manifold.

    Giving two coordinate charts `(U, \varphi)` and `(V, \psi)` on a
    topological manifold `M` of dimension `n` over a topological field `K`,
    the *transition map from* `(U, \varphi)` *to* `(V, \psi)` is the map

    .. MATH::

        \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset K^n
        \rightarrow \psi(U\cap V) \subset K^n.

    In other words, the transition map `\psi \circ \varphi^{-1}` expresses
    the coordinates `(y^1, \ldots, y^n)` of `(V, \psi)` in terms of the
    coordinates `(x^1, \ldots, x^n)` of `(U, \varphi)` on the open subset
    where the two charts intersect, i.e. on `U \cap V`.

    INPUT:

    - ``chart1`` -- chart `(U, \varphi)`
    - ``chart2`` -- chart `(V, \psi)`
    - ``transformations`` -- tuple (or list) `(Y_1, \ldots, Y_2)`, where
      `Y_i` is the symbolic expression of the coordinate `y^i` in terms
      of the coordinates `(x^1, \ldots, x^n)`

    EXAMPLES:

    Transition map on a 2-dimensional topological manifold::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: Y.<u,v> = M.chart()
        sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
        sage: X_to_Y
        Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        sage: type(X_to_Y)
        <class 'sage.manifolds.chart.CoordChange'>
        sage: X_to_Y.display()
        u = x + y
        v = x - y
    """
    def __init__(self, chart1, chart2, *transformations):
        r"""
        Construct a transition map.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: type(X_to_Y)
            <class 'sage.manifolds.chart.CoordChange'>
            sage: TestSuite(X_to_Y).run()
        """
        self._n1 = len(chart1._xx)
        self._n2 = len(chart2._xx)
        if len(transformations) != self._n2:
            raise ValueError("{} coordinate transformations ".format(self._n2)
                             + "must be provided")
        self._chart1 = chart1
        self._chart2 = chart2
        # The coordinate transformations are implemented via the class
        # MultiCoordFunction:
        self._transf = chart1.multifunction(*transformations)
        self._inverse = None
        # If the two charts are on the same open subset, the coordinate change
        # is added to the subset (and supersets) dictionary:
        if chart1.domain() == chart2.domain():
            domain = chart1.domain()
            for sdom in domain.open_supersets():
                sdom._coord_changes[(chart1, chart2)] = self

    def _repr_(self):
        r"""
        String representation of the transition map.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y._repr_()
            'Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))'
            sage: repr(X_to_Y)  # indirect doctest
            'Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))'
            sage: X_to_Y  # indirect doctest
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        """
        return "Change of coordinates from {} to {}".format(self._chart1,
                                                            self._chart2)

    def _latex_(self):
        r"""
        LaTeX representation of the transition map.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y._latex_()
            \left(M,(x, y)\right) \rightarrow \left(M,(u, v)\right)
            sage: latex(X_to_Y)  # indirect doctest
            \left(M,(x, y)\right) \rightarrow \left(M,(u, v)\right)
        """
        return latex(self._chart1) + r' \rightarrow ' + latex(self._chart2)

    def __eq__(self, other):
        r"""
        Equality operator.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y == X_to_Y
            True
            sage: X_to_Y1 = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y == X_to_Y1
            True
            sage: X_to_Y2 = X.transition_map(Y, [2*y, -x])
            sage: X_to_Y == X_to_Y2
            False
            sage: Z.<w,z> = M.chart()
            sage: X_to_Z = X.transition_map(Z, [x+y, x-y])
            sage: X_to_Y == X_to_Z
            False
        """
        if other is self:
            return True
        if not isinstance(other, CoordChange):
            return False
        return ((self._chart1 == other._chart1)
                and (self._chart2 == other._chart2)
                and (self._transf == other._transf))

    def __ne__(self, other):
        r"""
        Non-equality operator.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y2 = X.transition_map(Y, [2*y, -x])
            sage: X_to_Y != X_to_Y2
            True
        """
        return not (self == other)

    def __call__(self, *coords):
        r"""
        Compute the new coordinates from old ones.

        INPUT:

        - ``coords`` -- values of coordinates of ``chart1``

        OUTPUT: tuple of values of coordinates of ``chart2``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y(1,2)
            (3, -1)
        """
        return self._transf(*coords)

    def inverse(self):
        r"""
        Return the inverse coordinate transformation.

        If the inverse is not already known, it is computed here. If the
        computation fails, the inverse can be set by hand via the method
        :meth:`set_inverse`.

        OUTPUT:

        - an instance of :class:`CoordChange` representing the inverse of
          the current coordinate transformation

        EXAMPLES:

        Inverse of a coordinate transformation corresponding to a rotation
        in the Cartesian plane::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: phi = var('phi', domain='real')
            sage: xy_to_uv = c_xy.transition_map(c_uv,
            ....:                                [cos(phi)*x + sin(phi)*y,
            ....:                                 -sin(phi)*x + cos(phi)*y])
            sage: M.coord_changes()
            {(Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: uv_to_xy = xy_to_uv.inverse(); uv_to_xy
            Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
            sage: uv_to_xy.display()
            x = u*cos(phi) - v*sin(phi)
            y = v*cos(phi) + u*sin(phi)
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}

        The result is cached::

            sage: xy_to_uv.inverse() is uv_to_xy
            True

        We have as well::

            sage: uv_to_xy.inverse() is xy_to_uv
            True
        """
        from sage.symbolic.relation import solve
        if self._inverse is not None:
            return self._inverse
        # The computation is necessary:
        x1 = self._chart1._xx  # list of coordinates in chart1
        x2 = self._chart2._xx  # list of coordinates in chart2
        n1 = self._n1
        n2 = self._n2
        if n1 != n2:
            raise ValueError("the change of coordinates is not invertible " +
                             "(different number of coordinates in the two " +
                             "charts)")
        # New symbolic variables (different from x2 to allow for a
        #  correct solution even when chart2 = chart1):
        base_field = self._chart1.domain().base_field_type()
        if base_field == 'real':
            coord_domain = ['real' for i in range(n2)]
        elif base_field == 'complex':
            coord_domain = ['complex' for i in range(n2)]
        else:
            coord_domain = [None for i in range(n2)]
        for i in range(n2):
            if x2[i].is_positive():
                coord_domain[i] = 'positive'
        xp2 = [ SR.temp_var(domain=coord_domain[i]) for i in range(n2) ]
        xx2 = self._transf.expr()
        equations = [xp2[i] == xx2[i] for i in range(n2)]
        try:
            solutions = solve(equations, *x1, solution_dict=True)
        except RuntimeError:
            raise RuntimeError("the system could not be solved; use " +
                               "set_inverse() to set the inverse manually")
        substitutions = dict(zip(xp2, x2))
        if len(solutions) == 1:
            x2_to_x1 = [solutions[0][x1[i]].subs(substitutions)
                                                            for i in range(n1)]
            x2_to_x1_simpl = [] # to store simplified transformations
            for transf in x2_to_x1:
                try:
                    transf = self._chart2.simplify(transf)
                except AttributeError:
                    pass
                x2_to_x1_simpl.append(transf)
            x2_to_x1 = x2_to_x1_simpl
        else:
            list_x2_to_x1 = []
            for sol in solutions:
                if x2[0] in sol:
                    raise ValueError("the system could not be solved; use " +
                                     "set_inverse() to set the inverse " +
                                     "manually")
                try:
                    x2_to_x1 = [sol[x1[i]].subs(substitutions) for i in range(n1)]
                except KeyError: # sol is not a valid solution
                    continue
                x2_to_x1_simpl = [] # to store simplified transformations
                for transf in x2_to_x1:
                    try:
                        transf = self._chart2.simplify(transf)
                    except AttributeError:
                        pass
                    x2_to_x1_simpl.append(transf)
                x2_to_x1 = x2_to_x1_simpl
                if self._chart1.valid_coordinates(*x2_to_x1):
                    list_x2_to_x1.append(x2_to_x1)
            if len(list_x2_to_x1) == 0:
                raise ValueError("no solution found; use set_inverse() to " +
                                 "set the inverse manually")
            if len(list_x2_to_x1) > 1:
                print("Multiple solutions found: ")
                print(list_x2_to_x1)
                raise ValueError(
                   "non-unique solution to the inverse coordinate " +
                   "transformation; use set_inverse() to set the inverse " +
                   "manually")
            x2_to_x1 = list_x2_to_x1[0]
        self._inverse = type(self)(self._chart2, self._chart1, *x2_to_x1)
        self._inverse._inverse = self
        SR.cleanup_var(xp2)
        return self._inverse

    def set_inverse(self, *transformations, **kwds):
        r"""
        Set the inverse of the coordinate transformation.

        This is useful when the automatic computation via :meth:`inverse()`
        fails.

        INPUT:

        - ``transformations`` -- the inverse transformations expressed as a
          list of the expressions of the "old" coordinates in terms of the
          "new" ones
        - ``kwds`` -- optional arguments; valid keywords are

          - ``check`` -- boolean (default: ``True``); whether the
            provided transformations are checked to be indeed the inverse
            coordinate transformations
          - ``verbose`` -- boolean (default: ``False``); whether
            some details of the check are printed out; if ``False``, no
            output is printed if the check is passed (see example below)

        EXAMPLES:

        From spherical coordinates to Cartesian ones in the plane::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: U = M.open_subset('U') # complement of the half line {y=0, x>= 0}
            sage: c_cart.<x,y> = U.chart()
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = c_spher.transition_map(c_cart,
            ....:                                        [r*cos(ph), r*sin(ph)])
            sage: spher_to_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))
            Check of the inverse coordinate transformation:
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
              x == x  *passed*
              y == y  *passed*
            NB: a failed report can reflect a mere lack of simplification.

        As indicated, the failure for ``ph`` is due to a lack of simplification
        of the ``arctan2`` term, not to any error in the provided inverse
        formulas.

        We have now::

            sage: spher_to_cart.inverse()
            Change of coordinates from Chart (U, (x, y)) to Chart (U, (r, ph))
            sage: spher_to_cart.inverse().display()
            r = sqrt(x^2 + y^2)
            ph = arctan2(y, x)
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (U, (r, ph)),
              Chart (U, (x, y))): Change of coordinates from Chart (U, (r, ph))
               to Chart (U, (x, y)),
             (Chart (U, (x, y)),
              Chart (U, (r, ph))): Change of coordinates from Chart (U, (x, y))
               to Chart (U, (r, ph))}

        One can suppress the check of the provided formulas by means of the
        optional argument ``check=False``::

            sage: spher_to_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x),
            ....:                           check=False)

        However, it is not recommended to do so, the check being (obviously)
        useful to avoid some mistake. For instance, if the term
        ``sqrt(x^2+y^2)`` contains a typo (``x^3`` instead of ``x^2``),
        we get::

            sage: spher_to_cart.set_inverse(sqrt(x^3+y^2), atan2(y,x))
            Check of the inverse coordinate transformation:
              r == sqrt(r*cos(ph)^3 + sin(ph)^2)*r  **failed**
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
              x == sqrt(x^3 + y^2)*x/sqrt(x^2 + y^2)  **failed**
              y == sqrt(x^3 + y^2)*y/sqrt(x^2 + y^2)  **failed**
            NB: a failed report can reflect a mere lack of simplification.

        If the check is passed, no output is printed out::

            sage: M = Manifold(2, 'M')
            sage: X1.<x,y> = M.chart()
            sage: X2.<u,v> = M.chart()
            sage: X1_to_X2 = X1.transition_map(X2, [x+y, x-y])
            sage: X1_to_X2.set_inverse((u+v)/2, (u-v)/2)

        unless the option ``verbose`` is set to ``True``::

            sage: X1_to_X2.set_inverse((u+v)/2, (u-v)/2, verbose=True)
            Check of the inverse coordinate transformation:
              x == x  *passed*
              y == y  *passed*
              u == u  *passed*
              v == v  *passed*

        TESTS:

        Check that :issue:`31923` is fixed::

            sage: X1_to_X2.inverse().inverse() is X1_to_X2
            True

        Check of keyword arguments::

            sage: X1_to_X2.set_inverse((u+v)/2, (u-v)/2, bla=3)
            Traceback (most recent call last):
            ...
            TypeError: bla is not a valid keyword argument
        """
        check = kwds.pop('check', True)
        verbose = kwds.pop('verbose', False)
        for unknown_key in kwds:
            raise TypeError("{} is not a valid keyword "
                            "argument".format(unknown_key))
        self._inverse = type(self)(self._chart2, self._chart1,
                                   *transformations)
        self._inverse._inverse = self
        if check:
            infos = ["Check of the inverse coordinate transformation:"]
            x1 = self._chart1._xx
            x2 = self._chart2._xx
            x1_to_x1 = self._inverse(*(self(*x1)))
            x2_to_x2 = self(*(self._inverse(*x2)))
            any_failure = False  # a priori
            for x, xc in zip(x1, x1_to_x1):
                eq = x == self._chart1.simplify(xc)
                if bool(eq):
                    resu = '*passed*'
                else:
                    resu = '**failed**'
                    any_failure = True
                infos.append("  {}  {}".format(eq, resu))
            for x, xc in zip(x2, x2_to_x2):
                eq = x == self._chart2.simplify(xc)
                if bool(eq):
                    resu = '*passed*'
                else:
                    resu = '**failed**'
                    any_failure = True
                infos.append("  {}  {}".format(eq, resu))
            if any_failure:
                infos.append("NB: a failed report can reflect a mere lack of "
                             "simplification.")
            if verbose or any_failure:
                for li in infos:
                    print(li)

    def __mul__(self, other):
        r"""
        Composition with another change of coordinates.

        INPUT:

        - ``other`` -- another change of coordinate, the final chart of
          it is the initial chart of ``self``

        OUTPUT:

        - the change of coordinates `X_1 \to X_3`, where `X_1` is the initial
          chart of ``other`` and `X_3` is the final chart of ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: U.<u,v> = M.chart()
            sage: X_to_U = X.transition_map(U, (x+y, x-y))
            sage: W.<w,z> = M.chart()
            sage: U_to_W = U.transition_map(W, (u+cos(u)/2, v-sin(v)/2))
            sage: X_to_W = U_to_W * X_to_U; X_to_W
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (w, z))
            sage: X_to_W.display()
            w = 1/2*cos(x)*cos(y) - 1/2*sin(x)*sin(y) + x + y
            z = -1/2*cos(y)*sin(x) + 1/2*cos(x)*sin(y) + x - y
        """
        if not isinstance(other, CoordChange):
            raise TypeError("{} is not a change of coordinate".format(other))
        if other._chart2 != self._chart1:
            raise ValueError("composition not possible: " +
                             "{} is different from {}".format(other._chart2,
                                                              other._chart1))
        transf = self._transf(*(other._transf.expr()))
        return type(self)(other._chart1, self._chart2, *transf)

    def restrict(self, dom1, dom2=None):
        r"""
        Restriction to subsets.

        INPUT:

        - ``dom1`` -- open subset of the domain of ``chart1``
        - ``dom2`` -- (default: ``None``) open subset of the domain of
          ``chart2``; if ``None``, ``dom1`` is assumed

        OUTPUT:

        - the transition map between the charts restricted to the
          specified subsets

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: U = M.open_subset('U', coord_def={X: x>0, Y: u+v>0})
            sage: X_to_Y_U = X_to_Y.restrict(U); X_to_Y_U
            Change of coordinates from Chart (U, (x, y)) to Chart (U, (u, v))
            sage: X_to_Y_U.display()
            u = x + y
            v = x - y

        The result is cached::

            sage: X_to_Y.restrict(U) is X_to_Y_U
            True
        """
        if dom2 is None:
            dom2 = dom1
        ch1 = self._chart1.restrict(dom1)
        ch2 = self._chart2.restrict(dom2)
        if (ch1, ch2) in dom1.coord_changes():
            return dom1.coord_changes()[(ch1, ch2)]
        return type(self)(self._chart1.restrict(dom1),
                          self._chart2.restrict(dom2), *(self._transf.expr()))

    def display(self):
        r"""
        Display of the coordinate transformation.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        From spherical coordinates to Cartesian ones in the plane::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x>= 0}
            sage: c_cart.<x,y> = U.chart()
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = c_spher.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            sage: spher_to_cart.display()
            x = r*cos(ph)
            y = r*sin(ph)
            sage: latex(spher_to_cart.display())
            \left\{\begin{array}{lcl} x & = & r \cos\left({\phi}\right) \\
             y & = & r \sin\left({\phi}\right) \end{array}\right.

        A shortcut is ``disp()``::

            sage: spher_to_cart.disp()
            x = r*cos(ph)
            y = r*sin(ph)
        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        coords2 = self._chart2[:]
        n2 = len(coords2)
        expr = self._transf.expr('SR')
        rtxt = ""
        if n2 == 1:
            rlatex = r"\begin{array}{lcl}"
        else:
            rlatex = r"\left\{\begin{array}{lcl}"
        for i in range(n2):
            x2 = coords2[i]
            x2f = expr[i]
            rtxt += repr(x2) + " = " + repr(x2f) + "\n"
            rlatex += latex(x2) + r" & = & " + latex(x2f) + r"\\"
        rtxt = rtxt[:-1]  # remove the last new line
        rlatex = rlatex[:-2] + r"\end{array}"
        if n2 > 1:
            rlatex += r"\right."
        return FormattedExpansion(rtxt, rlatex)

    disp = display
