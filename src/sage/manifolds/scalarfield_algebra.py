# sage_setup: distribution = sagemath-symbolics
r"""
Algebra of Scalar Fields

The class :class:`ScalarFieldAlgebra` implements the commutative algebra
`C^0(M)` of scalar fields on a topological manifold `M` over a topological
field `K`. By *scalar field*, it
is meant a continuous function `M \to K`. The set
`C^0(M)` is an algebra over `K`, whose ring product is the pointwise
multiplication of `K`-valued functions, which is clearly commutative.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee2011]_
- [KN1963]_
"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.topological_spaces import TopologicalSpaces
from sage.symbolic.ring import SymbolicRing, SR
from sage.manifolds.scalarfield import ScalarField

class ScalarFieldAlgebra(UniqueRepresentation, Parent):
    r"""
    Commutative algebra of scalar fields on a topological manifold.

    If `M` is a topological manifold over a topological field `K`, the
    commutative algebra of scalar fields on `M` is the set `C^0(M)` of all
    continuous maps `M \to K`. The set `C^0(M)` is an algebra over `K`,
    whose ring product is the pointwise multiplication of `K`-valued
    functions, which is clearly commutative.

    If `K = \RR` or `K = \CC`, the field `K` over which the
    algebra `C^0(M)` is constructed is represented by the :class:`Symbolic
    Ring <sage.symbolic.ring.SymbolicRing>` ``SR``, since there is no exact
    representation of `\RR` nor `\CC`.

    INPUT:

    - ``domain`` -- the topological manifold `M` on which the scalar fields
      are defined

    EXAMPLES:

    Algebras of scalar fields on the sphere `S^2` and on some open
    subsets of it::

        sage: M = Manifold(2, 'M', structure='topological') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U')  # complement of the North pole
        sage: c_xy.<x,y> = U.chart()  # stereographic coordinates from the North pole
        sage: V = M.open_subset('V')  # complement of the South pole
        sage: c_uv.<u,v> = V.chart()  # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)    # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                                intersection_name='W',
        ....:                                restrictions1= x^2+y^2!=0,
        ....:                                restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: CM = M.scalar_field_algebra(); CM
        Algebra of scalar fields on the 2-dimensional topological manifold M
        sage: W = U.intersection(V)  # S^2 minus the two poles
        sage: CW = W.scalar_field_algebra(); CW
        Algebra of scalar fields on the Open subset W of the
         2-dimensional topological manifold M

    `C^0(M)` and `C^0(W)` belong to the category of commutative
    algebras over `\RR` (represented here by
    :class:`~sage.symbolic.ring.SymbolicRing`)::

        sage: CM.category()
        Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces
        sage: CM.base_ring()
        Symbolic Ring
        sage: CW.category()
        Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces
        sage: CW.base_ring()
        Symbolic Ring

    The elements of `C^0(M)` are scalar fields on `M`::

        sage: CM.an_element()
        Scalar field on the 2-dimensional topological manifold M
        sage: CM.an_element().display()  # this sample element is a constant field
        M → ℝ
        on U: (x, y) ↦ 2
        on V: (u, v) ↦ 2

    Those of `C^0(W)` are scalar fields on `W`::

        sage: CW.an_element()
        Scalar field on the Open subset W of the 2-dimensional topological
         manifold M
        sage: CW.an_element().display()  # this sample element is a constant field
        W → ℝ
        (x, y) ↦ 2
        (u, v) ↦ 2

    The zero element::

        sage: CM.zero()
        Scalar field zero on the 2-dimensional topological manifold M
        sage: CM.zero().display()
        zero: M → ℝ
        on U: (x, y) ↦ 0
        on V: (u, v) ↦ 0

    ::

        sage: CW.zero()
        Scalar field zero on the Open subset W of the 2-dimensional
         topological manifold M
        sage: CW.zero().display()
        zero: W → ℝ
           (x, y) ↦ 0
           (u, v) ↦ 0

    The unit element::

        sage: CM.one()
        Scalar field 1 on the 2-dimensional topological manifold M
        sage: CM.one().display()
        1: M → ℝ
        on U: (x, y) ↦ 1
        on V: (u, v) ↦ 1

    ::

        sage: CW.one()
        Scalar field 1 on the Open subset W of the 2-dimensional topological
         manifold M
        sage: CW.one().display()
        1: W → ℝ
          (x, y) ↦ 1
          (u, v) ↦ 1

    A generic element can be constructed by using a dictionary of
    the coordinate expressions defining the scalar field::

        sage: f = CM({c_xy: atan(x^2+y^2), c_uv: pi/2 - atan(u^2+v^2)}); f
        Scalar field on the 2-dimensional topological manifold M
        sage: f.display()
        M → ℝ
        on U: (x, y) ↦ arctan(x^2 + y^2)
        on V: (u, v) ↦ 1/2*pi - arctan(u^2 + v^2)
        sage: f.parent()
        Algebra of scalar fields on the 2-dimensional topological manifold M

    Specific elements can also be constructed in this way::

        sage: CM(0) == CM.zero()
        True
        sage: CM(1) == CM.one()
        True

    Note that the zero scalar field is cached::

        sage: CM(0) is CM.zero()
        True

    Elements can also be constructed by means of the method
    :meth:`~sage.manifolds.manifold.TopologicalManifold.scalar_field` acting
    on the domain (this allows one to set the name of the scalar field at the
    construction)::

        sage: f1 = M.scalar_field({c_xy: atan(x^2+y^2), c_uv: pi/2 - atan(u^2+v^2)},
        ....:                     name='f')
        sage: f1.parent()
        Algebra of scalar fields on the 2-dimensional topological manifold M
        sage: f1 == f
        True
        sage: M.scalar_field(0, chart='all') == CM.zero()
        True

    The algebra `C^0(M)` coerces to `C^0(W)` since `W` is an open
    subset of `M`::

        sage: CW.has_coerce_map_from(CM)
        True

    The reverse is of course false::

        sage: CM.has_coerce_map_from(CW)
        False

    The coercion map is nothing but the restriction to `W` of scalar fields
    on `M`::

        sage: fW = CW(f) ; fW
        Scalar field on the Open subset W of the
         2-dimensional topological manifold M
        sage: fW.display()
        W → ℝ
          (x, y) ↦ arctan(x^2 + y^2)
          (u, v) ↦ 1/2*pi - arctan(u^2 + v^2)

    ::

        sage: CW(CM.one()) == CW.one()
        True

    The coercion map allows for the addition of elements of `C^0(W)`
    with elements of `C^0(M)`, the result being an element of
    `C^0(W)`::

        sage: s = fW + f
        sage: s.parent()
        Algebra of scalar fields on the Open subset W of the
         2-dimensional topological manifold M
        sage: s.display()
        W → ℝ
          (x, y) ↦ 2*arctan(x^2 + y^2)
          (u, v) ↦ pi - 2*arctan(u^2 + v^2)

    Another coercion is that from the Symbolic Ring.
    Since the Symbolic Ring is the base ring for the algebra ``CM``, the
    coercion of a symbolic expression ``s`` is performed by the operation
    ``s*CM.one()``, which invokes the (reflected) multiplication operator.
    If the symbolic expression does not involve any chart coordinate,
    the outcome is a constant scalar field::

        sage: h = CM(pi*sqrt(2)) ; h
        Scalar field on the 2-dimensional topological manifold M
        sage: h.display()
        M → ℝ
        on U: (x, y) ↦ sqrt(2)*pi
        on V: (u, v) ↦ sqrt(2)*pi
        sage: a = var('a')
        sage: h = CM(a); h.display()
        M → ℝ
        on U: (x, y) ↦ a
        on V: (u, v) ↦ a

    If the symbolic expression involves some coordinate of one of the
    manifold's charts, the outcome is initialized only on the chart domain::

        sage: h = CM(a+x); h.display()
        M → ℝ
        on U: (x, y) ↦ a + x
        on W: (u, v) ↦ (a*u^2 + a*v^2 + u)/(u^2 + v^2)
        sage: h = CM(a+u); h.display()
        M → ℝ
        on W: (x, y) ↦ (a*x^2 + a*y^2 + x)/(x^2 + y^2)
        on V: (u, v) ↦ a + u

    If the symbolic expression involves coordinates of different charts,
    the scalar field is created as a Python object, but is not initialized,
    in order to avoid any ambiguity::

        sage: h = CM(x+u); h.display()
        M → ℝ

    TESTS:

    Ring laws::

        sage: h = CM(pi*sqrt(2))
        sage: s = f + h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ sqrt(2)*pi + arctan(x^2 + y^2)
        on V: (u, v) ↦ 1/2*pi*(2*sqrt(2) + 1) - arctan(u^2 + v^2)

    ::

        sage: s = f - h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ -sqrt(2)*pi + arctan(x^2 + y^2)
        on V: (u, v) ↦ -1/2*pi*(2*sqrt(2) - 1) - arctan(u^2 + v^2)

    ::

        sage: s = f*h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ sqrt(2)*pi*arctan(x^2 + y^2)
        on V: (u, v) ↦ 1/2*sqrt(2)*(pi^2 - 2*pi*arctan(u^2 + v^2))

    ::

        sage: s = f/h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ 1/2*sqrt(2)*arctan(x^2 + y^2)/pi
        on V: (u, v) ↦ 1/4*sqrt(2)*(pi - 2*arctan(u^2 + v^2))/pi

    ::

        sage: f*(h+f) == f*h + f*f
        True

    Ring laws with coercion::

        sage: f - fW == CW.zero()
        True
        sage: f/fW == CW.one()
        True
        sage: s = f*fW ; s
        Scalar field on the Open subset W of the 2-dimensional topological
         manifold M
        sage: s.display()
        W → ℝ
        (x, y) ↦ arctan(x^2 + y^2)^2
        (u, v) ↦ 1/4*pi^2 - pi*arctan(u^2 + v^2) + arctan(u^2 + v^2)^2
        sage: s/f == fW
        True

    Multiplication by a real number::

        sage: s = 2*f ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ 2*arctan(x^2 + y^2)
        on V: (u, v) ↦ pi - 2*arctan(u^2 + v^2)

    ::

        sage: 0*f == CM.zero()
        True
        sage: 1*f == f
        True
        sage: 2*(f/2) == f
        True
        sage: (f+2*f)/3 == f
        True
        sage: 1/3*(f+2*f) == f
        True

    The Sage test suite for algebras is passed::

        sage: TestSuite(CM).run()

    It is passed also for `C^0(W)`::

        sage: TestSuite(CW).run()
    """

    Element = ScalarField

    def __init__(self, domain):
        r"""
        Construct an algebra of scalar fields.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra(); CM
            Algebra of scalar fields on the 2-dimensional topological
             manifold M
            sage: type(CM)
            <class 'sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra_with_category'>
            sage: type(CM).__base__
            <class 'sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra'>
            sage: TestSuite(CM).run()
        """
        base_field = domain.base_field()
        if domain.base_field_type() in ['real', 'complex']:
            base_field = SR
        Parent.__init__(self, base=base_field,
                        category=CommutativeAlgebras(base_field) & TopologicalSpaces().Homsets())
        self._domain = domain
        self._populate_coercion_lists_()

    #### Methods required for any Parent
    def _element_constructor_(self, coord_expression=None, chart=None,
                              name=None, latex_name=None):
        r"""
        Construct a scalar field.

        INPUT:

        - ``coord_expression`` -- (default: ``None``) element(s) to construct
          the scalar field; this can be either

          - a scalar field defined on a domain that encompass ``self._domain``;
            then ``_element_constructor_`` return the restriction of
            the scalar field to ``self._domain``
          - a dictionary of coordinate expressions in various charts on the
            domain, with the charts as keys
          - a single coordinate expression; if the argument ``chart`` is
            ``'all'``, this expression is set to all the charts defined
            on the open set; otherwise, the expression is set in the
            specific chart provided by the argument ``chart``

        - ``chart`` -- (default: ``None``) chart defining the coordinates used
          in ``coord_expression`` when the latter is a single coordinate
          expression; if none is provided (default), the default chart of the
          open set is assumed. If ``chart=='all'``, ``coord_expression`` is
          assumed to be independent of the chart (constant scalar field).

        - ``name`` -- (default: ``None``) string; name (symbol) given to the
          scalar field

        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the scalar field; if none is provided, the LaTeX symbol is set to
          ``name``

        If ``coord_expression`` is ``None`` or incomplete, coordinate
        expressions can be added after the creation of the object, by means
        of the methods :meth:`add_expr`, :meth:`add_expr_by_continuation` and
        :meth:`set_expr`

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: f = CM({X: x+y^2}); f
            Scalar field on the 2-dimensional topological manifold M
            sage: f.display()
            M → ℝ
            (x, y) ↦ y^2 + x
            sage: f = CM({X: x+y^2}, name='f'); f
            Scalar field f on the 2-dimensional topological manifold M
            sage: f.display()
            f: M → ℝ
               (x, y) ↦ y^2 + x
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: CU = U.scalar_field_algebra()
            sage: fU = CU(f); fU
            Scalar field f on the Open subset U of the 2-dimensional topological
             manifold M
            sage: fU.display()
            f: U → ℝ
               (x, y) ↦ y^2 + x
        """
        try:
            if coord_expression.is_trivial_zero():
                return self.zero()
            elif (coord_expression - 1).is_trivial_zero():
                return self.one()
        except AttributeError:
            if coord_expression == 0:
                return self.zero()
            if coord_expression == 1:
                return self.one()
        if isinstance(coord_expression, ScalarField):
            if self._domain.is_subset(coord_expression._domain):
                # restriction of the scalar field to self._domain:
                return coord_expression.restrict(self._domain)
        else:
            # Anything going wrong here should produce a readable error:
            try:
                # generic constructor:
                resu = self.element_class(self,
                                          coord_expression=coord_expression,
                                          name=name, latex_name=latex_name,
                                          chart=chart)
            except TypeError:
                raise TypeError("cannot convert " +
                                "{} to a scalar ".format(coord_expression) +
                                "field on {}".format(self._domain))
        return resu

    def _an_element_(self):
        r"""
        Construct some element of the algebra.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: f = CM._an_element_(); f
            Scalar field on the 2-dimensional topological manifold M
            sage: f.display()
            M → ℝ
            (x, y) ↦ 2
        """
        return self.element_class(self, coord_expression=2, chart='all')

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from ``other``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: CM._coerce_map_from_(SR)
            True
            sage: SCR = SR.subring(no_variables=True); SCR
            Symbolic Constants Subring
            sage: CM._coerce_map_from_(SCR)
            True
            sage: CM._coerce_map_from_(X.function_ring())
            True
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: CU = U.scalar_field_algebra()
            sage: CM._coerce_map_from_(CU)
            False
            sage: CU._coerce_map_from_(CM)
            True
        """
        from .chart_func import ChartFunctionRing
        if isinstance(other, SymbolicRing):
            return True  # coercion from the base ring (multiplication by the
                         # algebra unit, i.e. self.one())
                         # cf. ScalarField._lmul_() for the implementation of
                         # the coercion map
        elif isinstance(other, ScalarFieldAlgebra):
            return self._domain.is_subset(other._domain)
        elif isinstance(other, ChartFunctionRing):
            return self._domain.is_subset(other._chart.domain())
        else:
            return False

    #### End of methods required for any Parent

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: CM = M.scalar_field_algebra()
            sage: CM._repr_()
            'Algebra of scalar fields on the 2-dimensional topological manifold M'
            sage: CM
            Algebra of scalar fields on the 2-dimensional topological manifold M
        """
        return "Algebra of scalar fields on the {}".format(self._domain)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: CM = M.scalar_field_algebra()
            sage: CM._latex_()
            'C^0 \\left(M\\right)'
            sage: latex(CM)
            C^0 \left(M\right)
        """
        return r"C^0 \left(" + self._domain._latex_() + r"\right)"

    @cached_method
    def zero(self):
        r"""
        Return the zero element of the algebra.

        This is nothing but the constant scalar field `0` on the manifold,
        where `0` is the zero element of the base field.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: z = CM.zero(); z
            Scalar field zero on the 2-dimensional topological manifold M
            sage: z.display()
            zero: M → ℝ
               (x, y) ↦ 0

        The result is cached::

            sage: CM.zero() is z
            True
        """
        coord_express = {chart: chart.zero_function()
                         for chart in self._domain.atlas()}
        zero = self.element_class(self,
                                  coord_expression=coord_express,
                                  name='zero', latex_name='0')
        zero._is_zero = True
        zero.set_immutable()
        return zero

    @cached_method
    def one(self):
        r"""
        Return the unit element of the algebra.

        This is nothing but the constant scalar field `1` on the manifold,
        where `1` is the unit element of the base field.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: h = CM.one(); h
            Scalar field 1 on the 2-dimensional topological manifold M
            sage: h.display()
            1: M → ℝ
               (x, y) ↦ 1

        The result is cached::

            sage: CM.one() is h
            True
        """
        coord_express = {chart: chart.one_function()
                         for chart in self._domain.atlas()}
        one = self.element_class(self, coord_expression=coord_express,
                                 name='1', latex_name='1')
        one.set_immutable()
        return one
