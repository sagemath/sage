# sage_setup: distribution = sagemath-symbolics
r"""
Coordinate Charts on Differentiable Manifolds

The class :class:`DiffChart` implements coordinate charts on a differentiable
manifold over a topological field `K` (in most applications, `K = \RR` or
`K = \CC`).

The subclass :class:`RealDiffChart` is devoted
to the case `K=\RR`, for which the concept of coordinate range is meaningful.
Moreover, :class:`RealDiffChart` is endowed with some plotting
capabilities (cf. method :meth:`~sage.manifolds.chart.RealChart.plot`).

Transition maps between charts are implemented via the class
:class:`DiffCoordChange`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- Chap. 1 of [Lee2013]_
"""
# ****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.manifolds.chart import Chart, RealChart, CoordChange
from sage.manifolds.differentiable.vectorframe import CoordFrame


class DiffChart(Chart):
    r"""
    Chart on a differentiable manifold.

    Given a differentiable manifold `M` of dimension `n` over a topological
    field `K`, a *chart* is a member `(U,\varphi)` of the manifold's
    differentiable atlas; `U` is then an open subset of `M` and
    `\varphi: U \rightarrow V \subset K^n` is a homeomorphism from
    `U` to an open subset `V` of `K^n`.

    The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
    `\varphi(p) = (x^1(p),\ldots,x^n(p))\in K^n` for any point `p\in U`, are
    called the *coordinates* of the chart `(U,\varphi)`.

    INPUT:

    - ``domain`` -- open subset `U` on which the chart is defined
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

      - ``'SR'`` -- Sage's default symbolic engine (Symbolic Ring)
      - ``'sympy'`` -- SymPy
      - ``None`` -- the default of :class:`~sage.manifolds.calculus_method.CalculusMethod`
        will be used
    - ``names`` -- (default: ``None``) unused argument, except if
      ``coordinates`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used).
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

    A chart on a complex 2-dimensional differentiable manifold::

        sage: M = Manifold(2, 'M', field='complex')
        sage: X = M.chart('x y'); X
        Chart (M, (x, y))
        sage: latex(X)
        \left(M,(x, y)\right)
        sage: type(X)
        <class 'sage.manifolds.differentiable.chart.DiffChart'>

    To manipulate the coordinates `(x,y)` as global variables, one has to set::

        sage: x,y = X[:]

    However, a shortcut is to use the declarator ``<x,y>`` in the left-hand
    side of the chart declaration (there is then no need to pass the string
    ``'x y'`` to ``chart()``)::

        sage: M = Manifold(2, 'M', field='complex')
        sage: X.<x,y> = M.chart(); X
        Chart (M, (x, y))

    The coordinates are then immediately accessible::

        sage: y
        y
        sage: x is X[0] and y is X[1]
        True

    The trick is performed by Sage preparser::

        sage: preparse("X.<x,y> = M.chart()")
        "X = M.chart(names=('x', 'y',)); (x, y,) = X._first_ngens(2)"

    Note that ``x`` and ``y`` declared in ``<x,y>`` are mere Python variable
    names and do not have to coincide with the coordinate symbols;
    for instance, one may write::

        sage: M = Manifold(2, 'M', field='complex')
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

        sage: M = Manifold(2, 'M', field='complex')
        sage: X.<x,y> = M.chart()

    In the above example, the chart X covers entirely the manifold M::

        sage: X.domain()
        2-dimensional complex manifold M

    Of course, one may declare a chart only on an open subset of M::

        sage: U = M.open_subset('U')
        sage: Y.<z1, z2> = U.chart(r'z1:\zeta_1 z2:\zeta_2'); Y
        Chart (U, (z1, z2))
        sage: Y.domain()
        Open subset U of the 2-dimensional complex manifold M

    In the above declaration, we have also specified some LaTeX writing
    of the coordinates different from the text one::

        sage: latex(z1)
        {\zeta_1}

    Note the prefix ``r`` in front of the string ``r'z1:\zeta_1 z2:\zeta_2'``;
    it makes sure that the backslash character is treated as an ordinary
    character, to be passed to the LaTeX interpreter.

    Periodic coordinates are declared through the keyword ``period=`` in the
    coordinate field::

        sage: N = Manifold(2, 'N', field='complex')
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

        sage: M1 = Manifold(2, 'M_1', field='complex', start_index=1)
        sage: Z.<u,v> = M1.chart()
        sage: Z[1], Z[2]
        (u, v)

    The full set of coordinates is obtained by means of the operator
    ``[:]``::

        sage: Y[:]
        (z1, z2)

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

    The action of the chart map `\varphi` on a point is obtained by means of
    the call operator, i.e. the operator ``()``::

        sage: p = M.point((1+i, 2), chart=X); p
        Point on the 2-dimensional complex manifold M
        sage: X(p)
        (I + 1, 2)
        sage: X(p) == p.coord(X)
        True

    A vector frame is naturally associated to each chart::

        sage: X.frame()
        Coordinate frame (M, (∂/∂x,∂/∂y))
        sage: Y.frame()
        Coordinate frame (U, (∂/∂z1,∂/∂z2))

    as well as a dual frame (basis of 1-forms)::

        sage: X.coframe()
        Coordinate coframe (M, (dx,dy))
        sage: Y.coframe()
        Coordinate coframe (U, (dz1,dz2))

    .. SEEALSO::

        :class:`~sage.manifolds.differentiable.chart.RealDiffChart` for charts
        on differentiable manifolds over `\RR`.
    """
    def __init__(self, domain, coordinates, calc_method=None, periods=None, coord_restrictions=None):
        r"""
        Construct a chart.

        TESTS::

            sage: M = Manifold(2, 'M', field='complex')
            sage: X.<x,y> = M.chart()
            sage: X
            Chart (M, (x, y))
            sage: type(X)
            <class 'sage.manifolds.differentiable.chart.DiffChart'>
            sage: assumptions() # no assumptions on x,y set by X._init_coordinates
            []
            sage: TestSuite(X).run()
        """
        super().__init__(domain, coordinates, calc_method=calc_method,
                         periods=periods, coord_restrictions=coord_restrictions)
        # Construction of the coordinate frame associated to the chart:
        self._frame = CoordFrame(self)
        self._coframe = self._frame._coframe

    def transition_map(self, other, transformations, intersection_name=None,
                       restrictions1=None, restrictions2=None):
        r"""
        Construct the transition map between the current chart,
        `(U,\varphi)` say, and another one, `(V,\psi)` say.

        If `n` is the manifold's dimension, the *transition map* is the
        map

        .. MATH::

            \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset K^n
            \rightarrow \psi(U\cap V) \subset K^n,

        where `K` is the manifold's base field. In other words, the
        transition map expresses the coordinates `(y^1,\ldots,y^n)` of
        `(V,\psi)` in terms of the coordinates `(x^1,\ldots,x^n)` of
        `(U,\varphi)` on the open subset where the two charts intersect, i.e.
        on `U\cap V`.

        By definition, the transition map `\psi\circ\varphi^{-1}` must be
        of class `C^k`, where `k` is the degree of differentiability of the
        manifold (cf.
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.diff_degree`).

        INPUT:

        - ``other`` -- the chart `(V,\psi)`
        - ``transformations`` -- tuple (or list) `(Y_1,\ldots,Y_2)`, where
          `Y_i` is the symbolic expression of the coordinate `y^i` in terms
          of the coordinates `(x^1,\ldots,x^n)`
        - ``intersection_name`` -- (default: ``None``) name to be given to the
          subset `U\cap V` if the latter differs from `U` or `V`
        - ``restrictions1`` -- (default: ``None``) list of conditions on the
          coordinates of the current chart that define `U\cap V` if the
          latter differs from `U`. ``restrictions1`` must be a list of
          of symbolic equalities or inequalities involving the
          coordinates, such as x>y or x^2+y^2 != 0. The items of the list
          ``restrictions1`` are combined with the ``and`` operator; if some
          restrictions are to be combined with the ``or`` operator instead,
          they have to be passed as a tuple in some single item of the list
          ``restrictions1``. For example, ``restrictions1`` = [x>y,
          (x!=0, y!=0), z^2<x] means (x>y) and ((x!=0) or (y!=0)) and (z^2<x).
          If the list ``restrictions1`` contains only one item, this item can
          be passed as such, i.e. writing x>y instead of the single-element
          list [x>y].
        - ``restrictions2`` -- (default: ``None``) list of conditions on the
          coordinates of the chart `(V,\psi)` that define `U\cap V` if the
          latter differs from `V` (see ``restrictions1`` for the syntax)

        OUTPUT:

        - The transition map `\psi\circ\varphi^{-1}` defined on `U\cap V`, as an
          instance of :class:`DiffCoordChange`.

        EXAMPLES:

        Transition map between two stereographic charts on the circle `S^1`::

            sage: M = Manifold(1, 'S^1')
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
            Set {S^1, U, V, W} of open subsets of the 1-dimensional differentiable manifold S^1
            sage: W = F['W']
            sage: W is U.intersection(V)
            True
            sage: M.atlas()
            [Chart (U, (x,)), Chart (V, (y,)), Chart (W, (x,)), Chart (W, (y,))]

        Transition map between the polar chart and the Cartesian one on
        `\RR^2`::

            sage: M = Manifold(2, 'R^2')
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

        In this case, no new subset has been created since `U\cap M = U`::

            sage: M.subset_family()
            Set {R^2, U} of open subsets of the 2-dimensional differentiable manifold R^2

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
        return DiffCoordChange(chart1, chart2, *transformations)

    def frame(self):
        r"""
        Return the vector frame (coordinate frame) associated with ``self``.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.CoordFrame`
          representing the coordinate frame

        EXAMPLES:

        Coordinate frame associated with some chart on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_xy.frame()
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: type(c_xy.frame())
            <class 'sage.manifolds.differentiable.vectorframe.CoordFrame_with_category'>

        Check that ``c_xy.frame()`` is indeed the coordinate frame associated
        with the coordinates `(x,y)`::

            sage: ex = c_xy.frame()[0] ; ex
            Vector field ∂/∂x on the 2-dimensional differentiable manifold M
            sage: ey = c_xy.frame()[1] ; ey
            Vector field ∂/∂y on the 2-dimensional differentiable manifold M
            sage: ex(M.scalar_field(x)).display()
            1: M → ℝ
               (x, y) ↦ 1
            sage: ex(M.scalar_field(y)).display()
            zero: M → ℝ
               (x, y) ↦ 0
            sage: ey(M.scalar_field(x)).display()
            zero: M → ℝ
               (x, y) ↦ 0
            sage: ey(M.scalar_field(y)).display()
            1: M → ℝ
               (x, y) ↦ 1
        """
        return self._frame

    def coframe(self):
        r"""
        Return the coframe (basis of coordinate differentials) associated
        with ``self``.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.CoordCoFrame`
          representing the coframe

        EXAMPLES:

        Coordinate coframe associated with some chart on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_xy.coframe()
            Coordinate coframe (M, (dx,dy))
            sage: type(c_xy.coframe())
            <class 'sage.manifolds.differentiable.vectorframe.CoordCoFrame_with_category'>

        Check that ``c_xy.coframe()`` is indeed the coordinate coframe
        associated with the coordinates `(x, y)`::

            sage: dx = c_xy.coframe()[0] ; dx
            1-form dx on the 2-dimensional differentiable manifold M
            sage: dy = c_xy.coframe()[1] ; dy
            1-form dy on the 2-dimensional differentiable manifold M
            sage: ex = c_xy.frame()[0] ; ex
            Vector field ∂/∂x on the 2-dimensional differentiable manifold M
            sage: ey = c_xy.frame()[1] ; ey
            Vector field ∂/∂y on the 2-dimensional differentiable manifold M
            sage: dx(ex).display()
            dx(∂/∂x): M → ℝ
               (x, y) ↦ 1
            sage: dx(ey).display()
            dx(∂/∂y): M → ℝ
               (x, y) ↦ 0
            sage: dy(ex).display()
            dy(∂/∂x): M → ℝ
               (x, y) ↦ 0
            sage: dy(ey).display()
            dy(∂/∂y): M → ℝ
               (x, y) ↦ 1
        """
        return self._coframe

    def restrict(self, subset, restrictions=None):
        r"""
        Return the restriction of ``self`` to some subset.

        If the current chart is `(U, \varphi)`, a *restriction* (or
        *subchart*) is a chart `(V, \psi)` such that `V \subset U`
        and `\psi = \varphi |_V`.

        If such subchart has not been defined yet, it is constructed here.

        The coordinates of the subchart bare the same names as the
        coordinates of the original chart.

        INPUT:

        - ``subset`` -- open subset `V` of the chart domain `U`
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

        - a :class:`DiffChart` `(V, \psi)`

        EXAMPLES:

        Coordinates on the unit open ball of  `\CC^2` as a subchart
        of the global coordinates of `\CC^2`::

            sage: M = Manifold(2, 'C^2', field='complex')
            sage: X.<z1, z2> = M.chart()
            sage: B = M.open_subset('B')
            sage: X_B = X.restrict(B, abs(z1)^2 + abs(z2)^2 < 1); X_B
            Chart (B, (z1, z2))
        """
        if subset == self.domain():
            return self
        if subset not in self._dom_restrict:
            resu = Chart.restrict(self, subset, restrictions=restrictions)
            # Update of superframes and subframes:
            resu._frame._superframes.update(self._frame._superframes)
            for sframe in self._frame._superframes:
                sframe._subframes.add(resu._frame)
                sframe._restrictions[subset] = resu._frame
            # The subchart frame is not a "top frame" in the supersets
            # (including self.domain()):
            for dom in self.domain().open_supersets():
                if resu._frame in dom._top_frames:
                    # it was added by the Chart constructor invoked in
                    # Chart.restrict above
                    dom._top_frames.remove(resu._frame)
        return self._dom_restrict[subset]

    def symbolic_velocities(self, left='D', right=None):
        r"""
        Return a list of symbolic variables ready to be used by the
        user as the derivatives of the coordinate functions with respect
        to a curve parameter (i.e. the velocities along the curve).
        It may actually serve to denote anything else than velocities,
        with a name including the coordinate functions.
        The choice of strings provided as 'left' and 'right' arguments
        is not entirely free since it must comply with Python
        prescriptions.

        INPUT:

        - ``left`` -- (default: ``D``) string to concatenate to the left
          of each coordinate functions of the chart
        - ``right`` -- (default: ``None``) string to concatenate to the
          right of each coordinate functions of the chart

        OUTPUT: list of symbolic expressions with the desired names

        EXAMPLES:

        Symbolic derivatives of the Cartesian coordinates of the
        3-dimensional Euclidean space::

            sage: R3 = Manifold(3, 'R3', start_index=1)
            sage: cart.<X,Y,Z> = R3.chart()
            sage: D = cart.symbolic_velocities(); D
            [DX, DY, DZ]
            sage: D = cart.symbolic_velocities(left='d', right='/dt'); D
            Traceback (most recent call last):
            ...
            ValueError: The name "dX/dt" is not a valid Python
             identifier.
            sage: D = cart.symbolic_velocities(left='d', right='_dt'); D
            [dX_dt, dY_dt, dZ_dt]
            sage: D = cart.symbolic_velocities(left='', right="'"); D
            Traceback (most recent call last):
            ...
            ValueError: The name "X'" is not a valid Python
             identifier.
            sage: D = cart.symbolic_velocities(left='', right='_dot'); D
            [X_dot, Y_dot, Z_dot]
            sage: R.<t> = manifolds.RealLine()
            sage: canon_chart = R.default_chart()
            sage: D = canon_chart.symbolic_velocities() ; D
            [Dt]
        """

        from sage.symbolic.ring import var

        # The case len(self[:]) = 1 is treated apart due to the
        # following fact.
        # In the case of several coordinates, the argument of 'var' (as
        # implemented below after the case len(self[:]) = 1) is a list
        # of strings of the form ['Dx1', 'Dx2', ...] and not a unique
        # string of the form 'Dx1 Dx2 ...'.
        # Although 'var' is supposed to accept both syntaxes, the first
        # one causes an error when it contains only one argument, due to
        # line 784 of sage/symbolic/ring.pyx :
        # "return self.symbol(name, latex_name=formatted_latex_name, domain=domain)"
        # In this line, the first argument 'name' of 'symbol' is a list
        # and not a string if the argument of 'var' is a list of one
        # string (of the type ['Dt']), which causes error in 'symbol'.
        # This might be corrected.
        if len(self[:]) == 1:
            string_vel = left + format(self[:][0]) # will raise an error
            # in case left is not a string
            if right is not None:
                string_vel += right # will raise an error in case right
                # is not a string

            # If the argument of 'var' contains only one word, for
            # instance::
            # var('Dt')
            # then 'var' does not return a tuple containing one symbolic
            # expression, but the symbolic expression itself.
            # This is taken into account below in order to return a list
            # containing one symbolic expression.
            return [var(string_vel)]

        list_strings_velocities = [left + format(coord_func)
                                   for coord_func in self[:]] # will
        # raise an error in case left is not a string

        if right is not None:
            list_strings_velocities = [str_vel + right for str_vel
                                       in list_strings_velocities] # will
            # raise an error in case right is not a string

        return list(var(list_strings_velocities))


#*****************************************************************************

class RealDiffChart(DiffChart, RealChart):
    r"""
    Chart on a differentiable manifold over `\RR`.

    Given a differentiable manifold `M` of dimension `n` over `\RR`,
    a *chart* is a member `(U,\varphi)` of the manifold's
    differentiable atlas; `U` is then an open subset of `M` and
    `\varphi: U \rightarrow V \subset \RR^n` is a homeomorphism from
    `U` to an open subset `V` of `\RR^n`.

    The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
    `\varphi(p) = (x^1(p),\ldots,x^n(p))\in \RR^n` for any point `p\in U`, are
    called the *coordinates* of the chart `(U,\varphi)`.

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

      - ``'SR'`` -- Sage's default symbolic engine (Symbolic Ring)
      - ``'sympy'`` -- SymPy
      - ``None`` -- the default of :class:`~sage.manifolds.calculus_method.CalculusMethod`
        will be used
    - ``names`` -- (default: ``None``) unused argument, except if
      ``coordinates`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used).
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

        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart = M.chart('x y z'); c_cart
        Chart (R^3, (x, y, z))
        sage: type(c_cart)
        <class 'sage.manifolds.differentiable.chart.RealDiffChart'>

    To have the coordinates accessible as global variables, one has to set::

        sage: (x,y,z) = c_cart[:]

    However, a shortcut is to use the declarator ``<x,y,z>`` in the left-hand
    side of the chart declaration (there is then no need to pass the string
    ``'x y z'`` to  ``chart()``)::

        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart(); c_cart
        Chart (R^3, (x, y, z))

    The coordinates are then immediately accessible::

        sage: y
        y
        sage: y is c_cart[2]
        True

    The trick is performed by Sage preparser::

        sage: preparse("c_cart.<x,y,z> = M.chart()")
        "c_cart = M.chart(names=('x', 'y', 'z',)); (x, y, z,) = c_cart._first_ngens(3)"

    Note that ``x, y, z`` declared in ``<x,y,z>`` are mere Python variable
    names and do not have to coincide with the coordinate symbols; for instance,
    one may write::

        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
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
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()

    Spherical coordinates on the subset `U` of `\RR^3` that is the
    complement of the half-plane `\{y=0, x\geq 0\}`::

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

    The full set of coordinates is obtained by means of the operator [:]::

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

    Manifold subsets have a *default chart*, which, unless changed via the
    method
    :meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart`,
    is the first defined chart on the subset (or on a open subset of it)::

        sage: M.default_chart()
        Chart (R^3, (x, y, z))
        sage: U.default_chart()
        Chart (U, (r, th, ph))

    The default charts are not privileged charts on the manifold, but rather
    charts whose name can be skipped in the argument list of functions having
    an optional ``chart=`` argument.

    The action of the chart map `\varphi` on a point is obtained by means of
    the call operator, i.e. the operator ``()``::

        sage: p = M.point((1,0,-2)); p
        Point on the 3-dimensional differentiable manifold R^3
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
    `\{y=0, x\geq 0\}`, we must have `y\not=0` or `x<0` on U. Accordingly,
    we set::

        sage: c_cartU.<x,y,z> = U.chart(coord_restrictions=lambda x,y,z: (y!=0, x<0))
        ....:    # the tuple (y!=0, x<0) means y!=0 or x<0
        ....:    #           [y!=0, x<0] would have meant y!=0 AND x<0
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

    A vector frame is naturally associated to each chart::

        sage: c_cart.frame()
        Coordinate frame (R^3, (∂/∂x,∂/∂y,∂/∂z))
        sage: c_spher.frame()
        Coordinate frame (U, (∂/∂r,∂/∂th,∂/∂ph))

    as well as a dual frame (basis of 1-forms)::

        sage: c_cart.coframe()
        Coordinate coframe (R^3, (dx,dy,dz))
        sage: c_spher.coframe()
        Coordinate coframe (U, (dr,dth,dph))

    Chart grids can be drawn in 2D or 3D graphics thanks to the method
    :meth:`~sage.manifolds.chart.RealChart.plot`.
    """
    def __init__(self, domain, coordinates, calc_method=None,
                 bounds=None, periods=None, coord_restrictions=None):
        r"""
        Construct a chart on a real differentiable manifold.

        TESTS::

            sage: forget()  # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: X
            Chart (M, (x, y))
            sage: type(X)
            <class 'sage.manifolds.differentiable.chart.RealDiffChart'>
            sage: assumptions()  # assumptions set in X._init_coordinates
            [x is real, y is real]
            sage: TestSuite(X).run()
        """
        RealChart.__init__(self, domain, coordinates, calc_method=calc_method,
                           bounds=bounds, periods=periods, coord_restrictions=coord_restrictions)
        # Construction of the coordinate frame associated to the chart:
        self._frame = CoordFrame(self)
        self._coframe = self._frame._coframe

    def restrict(self, subset, restrictions=None):
        r"""
        Return the restriction of the chart to some subset.

        If the current chart is `(U, \varphi)`, a *restriction* (or
        *subchart*) is a chart `(V, \psi)` such that `V \subset U`
        and `\psi = \varphi |_V`.

        If such subchart has not been defined yet, it is constructed here.

        The coordinates of the subchart bare the same names as the
        coordinates of the original chart.

        INPUT:

        - ``subset`` -- open subset `V` of the chart domain `U`
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

        - a :class:`RealDiffChart` `(V, \psi)`

        EXAMPLES:

        Cartesian coordinates on the unit open disc in `\RR^2` as a subchart
        of the global Cartesian coordinates::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: p = M.point((1/2, 0))
            sage: p in D
            True
            sage: q = M.point((1, 2))
            sage: q in D
            False

        Cartesian coordinates on the annulus `1 < \sqrt{x^2+y^2} < 2`::

            sage: A = M.open_subset('A')
            sage: c_cart_A = c_cart.restrict(A, [x^2+y^2>1, x^2+y^2<4])
            sage: p in A, q in A
            (False, False)
            sage: a = M.point((3/2,0))
            sage: a in A
            True
        """
        if subset == self.domain():
            return self
        if subset not in self._dom_restrict:
            resu = RealChart.restrict(self, subset, restrictions=restrictions)
            # Update of superframes and subframes:
            resu._frame._superframes.update(self._frame._superframes)
            for sframe in self._frame._superframes:
                sframe._subframes.add(resu._frame)
                sframe._restrictions[subset] = resu._frame
            # The subchart frame is not a "top frame" in the supersets
            # (including self.domain()):
            for dom in self.domain().open_supersets():
                if resu._frame in dom._top_frames:
                    # it was added by the Chart constructor invoked in
                    # Chart.restrict above
                    dom._top_frames.remove(resu._frame)
        return self._dom_restrict[subset]

#******************************************************************************

class DiffCoordChange(CoordChange):
    r"""
    Transition map between two charts of a differentiable manifold.

    Giving two coordinate charts `(U,\varphi)` and `(V,\psi)` on a
    differentiable manifold `M` of dimension `n` over a topological field `K`,
    the *transition map from* `(U,\varphi)` *to* `(V,\psi)` is the map

    .. MATH::

        \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset K^n
        \rightarrow \psi(U\cap V) \subset K^n,

    In other words, the transition map `\psi\circ\varphi^{-1}` expresses the
    coordinates `(y^1,\ldots,y^n)` of `(V,\psi)` in terms of the coordinates
    `(x^1,\ldots,x^n)` of `(U,\varphi)` on the open subset where the two
    charts intersect, i.e. on `U\cap V`.

    By definition, the transition map `\psi\circ\varphi^{-1}` must be
    of class `C^k`, where `k` is the degree of differentiability of the
    manifold (cf.
    :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.diff_degree`).

    INPUT:

    - ``chart1`` -- chart `(U,\varphi)`
    - ``chart2`` -- chart `(V,\psi)`
    - ``transformations`` -- tuple (or list) `(Y_1,\ldots,Y_2)`, where
      `Y_i` is the symbolic expression of the coordinate `y^i` in terms
      of the coordinates `(x^1,\ldots,x^n)`

    EXAMPLES:

    Transition map on a 2-dimensional differentiable manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: Y.<u,v> = M.chart()
        sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
        sage: X_to_Y
        Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        sage: type(X_to_Y)
        <class 'sage.manifolds.differentiable.chart.DiffCoordChange'>
        sage: X_to_Y.display()
        u = x + y
        v = x - y
    """
    def __init__(self, chart1, chart2, *transformations):
        r"""
        Construct a transition map.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: X_to_Y
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: type(X_to_Y)
            <class 'sage.manifolds.differentiable.chart.DiffCoordChange'>
            sage: TestSuite(X_to_Y).run(skip='_test_pickling')

        .. TODO::

            fix _test_pickling
        """
        CoordChange.__init__(self, chart1, chart2, *transformations)
        # Jacobian matrix:
        self._jacobian = self._transf.jacobian()
        # If the two charts are on the same open subset, the Jacobian matrix is
        # added to the dictionary of changes of frame:
        if chart1.domain() == chart2.domain():
            domain = chart1.domain()
            frame1 = chart1._frame
            frame2 = chart2._frame
            vf_module = domain.vector_field_module()
            ch_basis = vf_module.automorphism()
            ch_basis.add_comp(frame1)[:, chart1] = self._jacobian
            ch_basis.add_comp(frame2)[:, chart1] = self._jacobian
            vf_module._basis_changes[(frame2, frame1)] = ch_basis
            for sdom in domain.open_supersets():
                sdom._frame_changes[(frame2, frame1)] = ch_basis
            # The inverse is computed only if it does not exist already
            # (because if it exists it may have a simpler expression than that
            #  obtained from the matrix inverse)
            if (frame1, frame2) not in vf_module._basis_changes:
                ch_basis_inv = ch_basis.inverse()
                vf_module._basis_changes[(frame1, frame2)] = ch_basis_inv
                for sdom in domain.open_supersets():
                    sdom._frame_changes[(frame1, frame2)] = ch_basis_inv

    def jacobian(self):
        r"""
        Return the Jacobian matrix of ``self``.

        If ``self`` corresponds to the change of coordinates

        .. MATH::

            y^i = Y^i(x^1,\ldots,x^n)\qquad 1\leq i \leq n

        the Jacobian matrix `J` is given by

        .. MATH::

            J_{ij} = \frac{\partial Y^i}{\partial x^j}

        where `i` is the row index and `j` the column one.

        OUTPUT:

        - Jacobian matrix `J`, the elements `J_{ij}` of which being
          coordinate functions
          (cf. :class:`~sage.manifolds.chart_func.ChartFunction`)

        EXAMPLES:

        Jacobian matrix of a 2-dimensional transition map::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y^2, 3*x-y])
            sage: X_to_Y.jacobian()
            [  1 2*y]
            [  3  -1]

        Each element of the Jacobian matrix is a coordinate function::

            sage: parent(X_to_Y.jacobian()[0,0])
            Ring of chart functions on Chart (M, (x, y))
        """
        return self._jacobian  # has been computed in __init__

    @cached_method
    def jacobian_det(self):
        r"""
        Return the Jacobian determinant of ``self``.

        The Jacobian determinant is the determinant of the Jacobian
        matrix (see :meth:`jacobian`).

        OUTPUT:

        - determinant of the Jacobian matrix `J` as a coordinate
          function
          (cf. :class:`~sage.manifolds.chart_func.ChartFunction`)

        EXAMPLES:

        Jacobian determinant of a 2-dimensional transition map::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y^2, 3*x-y])
            sage: X_to_Y.jacobian_det()
            -6*y - 1
            sage: X_to_Y.jacobian_det() == det(X_to_Y.jacobian())
            True

        The Jacobian determinant is a coordinate function::

            sage: parent(X_to_Y.jacobian_det())
            Ring of chart functions on Chart (M, (x, y))
        """
        return self._transf.jacobian_det()
