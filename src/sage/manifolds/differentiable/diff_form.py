r"""
Differential Forms

Let `U` and `M` be two differentiable manifolds.
Given a positive integer `p` and a differentiable map `\Phi: U \rightarrow M`,
a *differential form of degree* `p`, or `p`-*form*,
*along* `U` *with values on* `M` is a field along `U` of alternating
multilinear forms of degree `p` in the tangent spaces to `M`.
The standard case of a differential form *on* a differentiable manifold
corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
`\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

Two classes implement differential forms, depending whether the manifold
`M` is parallelizable:

* :class:`DiffFormParal` when `M` is parallelizable
* :class:`DiffForm` when `M` is not assumed parallelizable.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version
- Joris Vankerschaver (2010): developed a previous class,
  ``DifferentialForm`` (cf. :issue:`24444`), which inspired the storage of the
  nonzero components as a dictionary whose keys are the indices.
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [KN1963]_
- [Lee2013]_
"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2010 Joris Vankerschaver <joris.vankerschaver@gmail.com>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import annotations
from typing import Optional, Union, TYPE_CHECKING
from sage.misc.cachefunc import cached_method
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal

if TYPE_CHECKING:
    from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
    from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
    from sage.manifolds.differentiable.symplectic_form import SymplecticForm


class DiffForm(TensorField):
    r"""
    Differential form with values on a generic (i.e. a priori not
    parallelizable) differentiable manifold.

    Given a differentiable manifold `U`, a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M` and a positive
    integer `p`, a *differential form of degree* `p` (or `p`-*form*)
    *along* `U` *with values on* `M\supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(0,p)}M

    (`T^{(0,p)}M` being the tensor bundle of type `(0,p)` over `M`) such that

    .. MATH::

        \forall x \in U,\quad a(x) \in \Lambda^p(T_{\Phi(x)}^* M) ,

    where `T_{\Phi(x)}^* M` is the dual of the tangent space to `M` at
    `\Phi(x)` and `\Lambda^p` stands for the exterior power of degree `p` (cf.
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule`).
    In other words, `a(x)` is an alternating multilinear form of degree `p` of
    the tangent vector space `T_{\Phi(x)} M`.

    The standard case of a differential form *on* a
    manifold `M` corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other
    common cases are `\Phi` being an immersion and `\Phi` being a curve in `M`
    (`U` is then an open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, the class :class:`DiffFormParal`
        must be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi`
    - ``degree`` -- the degree of the differential form (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the differential form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      differential form; if none is provided, the LaTeX symbol is set to
      ``name``

    EXAMPLES:

    Differential form of degree 2 on a non-parallelizable 2-dimensional
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
        ....:                                restrictions1= x>0, restrictions2= u+v>0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a = M.diff_form(2, name='a') ; a
        2-form a on the 2-dimensional differentiable manifold M
        sage: a.parent()
        Module Omega^2(M) of 2-forms on the 2-dimensional differentiable
         manifold M
        sage: a.degree()
        2

    Setting the components of ``a``::

        sage: a[eU,0,1] = x*y^2 + 2*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = (x*y^2 + 2*x) dx∧dy
        sage: a.display(eV)
        a = (-1/16*u^3 + 1/16*u*v^2 - 1/16*v^3
         + 1/16*(u^2 - 8)*v - 1/2*u) du∧dv

    A 1-form on ``M``::

        sage: a = M.one_form(name='a') ; a
        1-form a on the 2-dimensional differentiable manifold M
        sage: a.parent()
        Module Omega^1(M) of 1-forms on the 2-dimensional differentiable
         manifold M
        sage: a.degree()
        1

    Setting the components of the 1-form in a consistent way::

        sage: a[eU,:] = [-y, x]
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = -y dx + x dy
        sage: a.display(eV)
        a = 1/2*v du - 1/2*u dv

    It is also possible to set the components at the 1-form definition, via
    a dictionary whose keys are the vector frames::

        sage: a1 = M.one_form({eU: [-y, x], eV: [v/2, -u/2]}, name='a')
        sage: a1 == a
        True

    The exterior derivative of the 1-form is a 2-form::

        sage: da = a.exterior_derivative() ; da
        2-form da on the 2-dimensional differentiable manifold M
        sage: da.display(eU)
        da = 2 dx∧dy
        sage: da.display(eV)
        da = -du∧dv

    The exterior derivative can also be obtained by applying the function
    ``diff`` to a differentiable form::

        sage: diff(a) is a.exterior_derivative()
        True

    Another 1-form defined by its components in ``eU``::

        sage: b = M.one_form(1+x*y, x^2, frame=eU, name='b')

    Since ``eU`` is the default vector frame on ``M``, it can be omitted
    in the definition::

        sage: b = M.one_form(1+x*y, x^2, name='b')
        sage: b.add_comp_by_continuation(eV, W, c_uv)

    Adding two 1-forms results in another 1-form::

        sage: s = a + b ; s
        1-form a+b on the 2-dimensional differentiable manifold M
        sage: s.display(eU)
        a+b = ((x - 1)*y + 1) dx + (x^2 + x) dy
        sage: s.display(eV)
        a+b = (1/4*u^2 + 1/4*(u + 2)*v + 1/2) du
         + (-1/4*u*v - 1/4*v^2 - 1/2*u + 1/2) dv

    The exterior product of two 1-forms is a 2-form::

        sage: s = a.wedge(b) ; s
        2-form a∧b on the 2-dimensional differentiable manifold M
        sage: s.display(eU)
        a∧b = (-2*x^2*y - x) dx∧dy
        sage: s.display(eV)
        a∧b = (1/8*u^3 - 1/8*u*v^2 - 1/8*v^3 + 1/8*(u^2 + 2)*v + 1/4*u) du∧dv

    Multiplying a 1-form by a scalar field results in another 1-form::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
        sage: s = f*a ; s
        1-form f*a on the 2-dimensional differentiable manifold M
        sage: s.display(eU)
        f*a = (-x^2*y - 2*x*y^2 - y^3) dx + (x^3 + 2*x^2*y + x*y^2) dy
        sage: s.display(eV)
        f*a = 1/2*u^2*v du - 1/2*u^3 dv


    .. RUBRIC:: Examples with SymPy as the symbolic engine

    From now on, we ask that all symbolic calculus on manifold `M` are
    performed by SymPy::

        sage: M.set_calculus_method('sympy')

    We define a 2-form `a` as above::

        sage: a = M.diff_form(2, name='a')
        sage: a[eU,0,1] = x*y^2 + 2*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = (x*y**2 + 2*x) dx∧dy
        sage: a.display(eV)
        a = (-u**3/16 + u**2*v/16 + u*v**2/16 - u/2 - v**3/16 - v/2) du∧dv

    A 1-form on ``M``::

        sage: a = M.one_form(-y, x, name='a')
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = -y dx + x dy
        sage: a.display(eV)
        a = v/2 du - u/2 dv

    The exterior derivative of ``a``::

        sage: da = a.exterior_derivative()
        sage: da.display(eU)
        da = 2 dx∧dy
        sage: da.display(eV)
        da = -du∧dv

    Another 1-form::

        sage: b = M.one_form(1+x*y, x^2, name='b')
        sage: b.add_comp_by_continuation(eV, W, c_uv)

    Adding two 1-forms::

        sage: s = a + b
        sage: s.display(eU)
        a+b = (x*y - y + 1) dx + x*(x + 1) dy
        sage: s.display(eV)
        a+b = (u**2/4 + u*v/4 + v/2 + 1/2) du + (-u*v/4 - u/2 - v**2/4 + 1/2) dv

    The exterior product of two 1-forms::

        sage: s = a.wedge(b)
        sage: s.display(eU)
        a∧b = x*(-2*x*y - 1) dx∧dy
        sage: s.display(eV)
        a∧b = (u**3/8 + u**2*v/8 - u*v**2/8 + u/4 - v**3/8 + v/4) du∧dv

    Multiplying a 1-form by a scalar field::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
        sage: s = f*a
        sage: s.display(eU)
        f*a = y*(-x**2 - 2*x*y - y**2) dx + x*(x**2 + 2*x*y + y**2) dy
        sage: s.display(eV)
        f*a = u**2*v/2 du - u**3/2 dv
    """
    def __init__(self, vector_field_module, degree, name=None, latex_name=None):
        r"""
        Construct a differential form.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``DiffForm``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: A = M.diff_form_module(2)
            sage: XM = M.vector_field_module()
            sage: a = A.element_class(XM, 2, name='a'); a
            2-form a on the 2-dimensional differentiable manifold M
            sage: a[e_xy,0,1] = x+y
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.diff_form``::

            sage: a1 = M.diff_form(2, name='a'); a1
            2-form a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        .. TODO::

            Fix ``_test_pickling`` (in the superclass :class:`TensorField`).
        """
        TensorField.__init__(self, vector_field_module, (0, degree), name=name,
                             latex_name=latex_name, antisym=range(degree),
                             parent=vector_field_module.dual_exterior_power(degree))
        self._init_derived()  # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: a = M.diff_form(2, name='a')
            sage: a._repr_()
            '2-form a on the 3-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            '2-form a on the 3-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            2-form a on the 3-dimensional differentiable manifold M
            sage: b = M.diff_form(2)
            sage: b._repr_()
            '2-form on the 3-dimensional differentiable manifold M'
        """
        description = "{}-form ".format(self._tensor_rank)
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: a = M.diff_form(2, name='a')
            sage: a1 = a._new_instance(); a1
            2-form on the 3-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True
        """
        return type(self)(self._vmodule, self._tensor_rank)

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: a = M.diff_form(2, name='a')
            sage: a._del_derived()
        """
        TensorField._del_derived(self)
        self.exterior_derivative.clear_cache()

    @cached_method
    def exterior_derivative(self) -> DiffForm:
        r"""
        Compute the exterior derivative of ``self``.

        OUTPUT:

        - instance of :class:`DiffForm` representing the exterior derivative
          of the differential form

        EXAMPLES:

        Exterior derivative of a 1-form on the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()

        The 1-form::

            sage: a = M.one_form({e_xy: [-y^2, x^2]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a.display(e_xy)
            a = -y^2 dx + x^2 dy
            sage: a.display(e_uv)
            a = -(2*u^3*v - u^2*v^2 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du
             + (u^4 - u^2*v^2 + 2*u*v^3)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv

        Its exterior derivative::

            sage: da = a.exterior_derivative(); da
            2-form da on the 2-dimensional differentiable manifold M
            sage: da.display(e_xy)
            da = (2*x + 2*y) dx∧dy
            sage: da.display(e_uv)
            da = -2*(u + v)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) du∧dv

        The result is cached, i.e. is not recomputed unless ``a`` is changed::

            sage: a.exterior_derivative() is da
            True

        Instead of invoking the method :meth:`exterior_derivative`, one may
        use the global function ``diff``::

            sage: diff(a) is a.exterior_derivative()
            True

        Let us check Cartan's identity::

            sage: v = M.vector_field({e_xy: [-y, x]}, name='v')
            sage: v.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a.lie_der(v) == v.contract(diff(a)) + diff(a(v))  # long time
            True
        """
        from sage.tensor.modules.format_utilities import (
            format_unop_txt,
            format_unop_latex,
        )

        vmodule = self._vmodule  # shortcut
        rname = format_unop_txt("d", self._name)
        rlname = format_unop_latex(r"\mathrm{d}", self._latex_name)
        resu = vmodule.alternating_form(
            self._tensor_rank + 1, name=rname, latex_name=rlname
        )
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst.exterior_derivative()
        return resu

    derivative = exterior_derivative  # allows one to use functional notation,
                                      # e.g. diff(a) for a.exterior_derivative()

    def wedge(self, other: DiffForm) -> DiffForm:
        r"""
        Exterior product with another differential form.

        INPUT:

        - ``other`` -- another differential form (on the same manifold)

        OUTPUT:

        - instance of :class:`DiffForm` representing the exterior product
          ``self ∧ other``

        EXAMPLES:

        Exterior product of two 1-forms on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1) # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coord. (North and South)
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: a = M.one_form({e_xy: [y, x]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: b = M.one_form({e_xy: [x^2 + y^2, y]}, name='b')
            sage: b.add_comp_by_continuation(e_uv, W, c_uv)
            sage: c = a.wedge(b); c
            2-form a∧b on the 2-dimensional differentiable manifold S^2
            sage: c.display(e_xy)
            a∧b = (-x^3 - (x - 1)*y^2) dx∧dy
            sage: c.display(e_uv)
            a∧b = -(v^2 - u)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du∧dv

        If one of the two operands is unnamed, the result is unnamed too::

            sage: b1 = M.diff_form(1)  # no name set
            sage: b1[e_xy,:] = x^2 + y^2, y
            sage: b1.add_comp_by_continuation(e_uv, W, c_uv)
            sage: c1 = a.wedge(b1); c1
            2-form on the 2-dimensional differentiable manifold S^2
            sage: c1.display(e_xy)
            (-x^3 - (x - 1)*y^2) dx∧dy

        To give a name to the result, one shall use the method
        :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.set_name`::

            sage: c1.set_name('c');  c1
            2-form c on the 2-dimensional differentiable manifold S^2
            sage: c1.display(e_xy)
            c = (-x^3 - (x - 1)*y^2) dx∧dy

        Wedging with scalar fields yields the multiplication from right::

            sage: f = M.scalar_field(x, name='f')
            sage: f.add_expr_by_continuation(c_uv, W)
            sage: t = a.wedge(f)
            sage: t.display()
            f*a = x*y dx + x^2 dy
        """
        if other._tensor_rank == 0:
            return self * other
        from sage.typeset.unicode_characters import unicode_wedge
        from sage.tensor.modules.format_utilities import is_atomic
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior product")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior product")
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(other._ambient_domain)
        resu_degree = self._tensor_rank + other._tensor_rank
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        # Facilitate computations involving zero:
        if resu_degree > ambient_dom_resu._dim:
            return dom_resu.diff_form_module(resu_degree,
                                             dest_map=dest_map_resu).zero()
        if self._is_zero or other._is_zero:
            return dom_resu.diff_form_module(resu_degree,
                                             dest_map=dest_map_resu).zero()
        if self is other and (self._tensor_rank % 2) == 1:
            return dom_resu.diff_form_module(resu_degree,
                                             dest_map=dest_map_resu).zero()
        # Generic case:
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleAltForm version:
            return FreeModuleAltForm.wedge(self_r, other_r)
        # otherwise, the result is created here:
        resu_name = None
        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu_name = sname + unicode_wedge + oname
        resu_latex_name = None
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            resu_latex_name = slname + r'\wedge ' + olname
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        resu = vmodule.alternating_form(resu_degree, name=resu_name,
                                        latex_name=resu_latex_name)
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                resu._restrictions[dom] = self_r._restrictions[dom].wedge(
                                          other_r._restrictions[dom])
        return resu

    def degree(self) -> int:
        r"""
        Return the degree of ``self``.

        OUTPUT: integer `p` such that the differential form is a `p`-form

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: a = M.diff_form(2); a
            2-form on the 3-dimensional differentiable manifold M
            sage: a.degree()
            2
            sage: b = M.diff_form(1); b
            1-form on the 3-dimensional differentiable manifold M
            sage: b.degree()
            1
        """
        return self._tensor_rank

    def hodge_dual(
        self,
        nondegenerate_tensor: Union[
            PseudoRiemannianMetric, SymplecticForm, None
        ] = None,
        minus_eigenvalues_convention: bool = False,
    ) -> DiffForm:
        r"""
        Compute the Hodge dual of the differential form with respect to some non-degenerate
        bilinear form (Riemannian metric or symplectic form).

        If the differential form is a `p`-form `A`, its *Hodge dual* with
        respect to the non-degenerate form `g` is the
        `(n-p)`-form `*A` defined by

        .. MATH::

            *A_{i_1\ldots i_{n-p}} = \frac{1}{p!} A^{k_1\ldots k_p}
                \epsilon_{k_1\ldots k_p\, i_1\ldots i_{n-p}}

        where `n` is the manifold's dimension, `\epsilon` is the volume
        `n`-form associated with `g` (see
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.volume_form`)
        and the indices `k_1,\ldots, k_p` are raised with `g`.
        If `g` is a pseudo-Riemannian metric, sometimes an additional multiplicative
        factor of `(-1)^s` is introduced on the right-hand side,
        where `s` is the number of negative eigenvalues of `g`.
        This convention can be enforced by setting the option
        ``minus_eigenvalues_convention``.

        INPUT:

        - ``nondegenerate_tensor`` -- a non-degenerate bilinear form defined on the same manifold
          as the current differential form; must be an instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric` or
          :class:`~sage.manifolds.differentiable.symplectic_form.SymplecticForm`.
          If none is provided, the ambient domain of ``self`` is supposed to be endowed
          with a default metric and this metric is then used.
        - ``minus_eigenvalues_convention`` -- if `true`, a factor of `(-1)^s` is
            introduced with `s` being the number of negative eigenvalues of the
            ``nondegenerate_tensor``.

        OUTPUT:

        - the `(n-p)`-form `*A`

        EXAMPLES:

        Hodge dual of a 1-form on the 2-sphere equipped with the standard
        metric: we first construct `\mathbb{S}^2` and its metric `g`::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coord. (North and South)
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: g = M.metric('g')
            sage: g[eU,1,1], g[eU,2,2] = 4/(1+x^2+y^2)^2, 4/(1+x^2+y^2)^2
            sage: g[eV,1,1], g[eV,2,2] = 4/(1+u^2+v^2)^2, 4/(1+u^2+v^2)^2

        We endow `S^2` with the orientation defined by the stereographic
        frame from the North pole, i.e. ``eU``; ``eV`` is then left-handed and
        in order to define an orientation on the whole manifold, we introduce a
        vector frame on ``V`` by swapping ``eV``'s vectors::

            sage: f = V.vector_frame('f', (eV[2], eV[1]))
            sage: M.set_orientation([eU, f])

        Then we construct the 1-form and take its Hodge dual w.r.t. `g`::

            sage: a = M.one_form({eU: [-y, x]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, c_uv)
            sage: a.display(eU)
            a = -y dx + x dy
            sage: a.display(eV)
            a = -v/(u^4 + 2*u^2*v^2 + v^4) du + u/(u^4 + 2*u^2*v^2 + v^4) dv
            sage: sa = a.hodge_dual(g); sa
            1-form *a on the 2-dimensional differentiable manifold S^2
            sage: sa.display(eU)
            *a = -x dx - y dy
            sage: sa.display(eV)
            *a = u/(u^4 + 2*u^2*v^2 + v^4) du + v/(u^4 + 2*u^2*v^2 + v^4) dv

        Instead of calling the method :meth:`hodge_dual` on the differential
        form, one can invoke the method
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.hodge_star`
        of the metric::

            sage: a.hodge_dual(g) == g.hodge_star(a)
            True

        For a 1-form and a Riemannian metric in dimension 2, the Hodge dual
        applied twice is minus the identity::

            sage: ssa = sa.hodge_dual(g); ssa
            1-form **a on the 2-dimensional differentiable manifold S^2
            sage: ssa == -a
            True

        The Hodge dual of the metric volume 2-form is the constant scalar
        field 1 (considered as a 0-form)::

            sage: eps = g.volume_form(); eps
            2-form eps_g on the 2-dimensional differentiable manifold S^2
            sage: eps.display(eU)
            eps_g = 4/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) dx∧dy
            sage: eps.display(eV)
            eps_g = -4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) du∧dv
            sage: seps = eps.hodge_dual(g); seps
            Scalar field *eps_g on the 2-dimensional differentiable manifold S^2
            sage: seps.display()
            *eps_g: S^2 → ℝ
            on U: (x, y) ↦ 1
            on V: (u, v) ↦ 1

        Hodge dual of a 1-form in the Euclidean space `R^3`::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')  # the Euclidean metric
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: var('Ax Ay Az')
            (Ax, Ay, Az)
            sage: a = M.one_form(Ax, Ay, Az, name='A')
            sage: sa = a.hodge_dual(g) ; sa
            2-form *A on the 3-dimensional differentiable manifold M
            sage: sa.display()
            *A = Az dx∧dy - Ay dx∧dz + Ax dy∧dz
            sage: ssa = sa.hodge_dual(g) ; ssa
            1-form **A on the 3-dimensional differentiable manifold M
            sage: ssa.display()
            **A = Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Riemannian metric in dimension 3
            True

        See the documentation of
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.hodge_star`
        for more examples.

        TESTS:
        Fall back to use (ambient) metric::

            sage: M = Manifold(3, 'M', start_index=1, structure='Riemannian')
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric()
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: var('Ax Ay Az')
            (Ax, Ay, Az)
            sage: a = M.one_form(Ax, Ay, Az, name='A')
            sage: a.hodge_dual().display()
            *A = Az dx∧dy - Ay dx∧dz + Ax dy∧dz
        """
        from sage.functions.other import factorial
        from sage.tensor.modules.format_utilities import (
            format_unop_txt,
            format_unop_latex,
        )

        if nondegenerate_tensor is None:
            nondegenerate_tensor = self._vmodule._ambient_domain.metric()

        p = self.tensor_type()[1]
        # For performance reasons, we raise the indicies of the volume form
        # and not of the differential form; in the symplectic case this is wrong by
        # a factor of (-1)^p, which will be corrected below
        eps = nondegenerate_tensor.volume_form(p)
        if p == 0:
            common_domain = nondegenerate_tensor.domain().intersection(self.domain())
            result = self.restrict(common_domain) * eps.restrict(common_domain)
        else:
            result = self.contract(*range(p), eps, *range(p))
            if p > 1:
                result = result / factorial(p)
            if minus_eigenvalues_convention:
                from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
                if isinstance(nondegenerate_tensor, PseudoRiemannianMetric):
                    result = result * nondegenerate_tensor._indic_signat
            from sage.manifolds.differentiable.symplectic_form import SymplecticForm
            if isinstance(nondegenerate_tensor, SymplecticForm):
                # correction because we lifted the indicies of the volume (see above)
                result = result * (-1)**p

        result.set_name(
            name=format_unop_txt("*", self._name),
            latex_name=format_unop_latex(r"\star ", self._latex_name),
        )
        return result

    def interior_product(self, qvect):
        r"""
        Interior product with a multivector field.

        If ``self`` is a differential form `A` of degree `p` and `B` is a
        multivector field of degree `q\geq p` on the same manifold, the
        interior product of `A` by `B` is the multivector field `\iota_A B` of
        degree `q-p` defined by

        .. MATH::

            (\iota_A B)^{i_1\ldots i_{q-p}} = A_{k_1\ldots k_p}
                            B^{k_1\ldots k_p i_1\ldots i_{q-p}}

        .. NOTE::

            ``A.interior_product(B)`` yields the same result as
            ``A.contract(0,..., p-1, B, 0,..., p-1)`` (cf.
            :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.contract`),
            but ``interior_product`` is more efficient, the alternating
            character of `A` being not used to reduce the computation in
            :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.contract`

        INPUT:

        - ``qvect`` -- multivector field `B` (instance of
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`);
          the degree of `B` must be at least equal to the degree of ``self``

        OUTPUT:

        - scalar field (case `p=q`) or
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
          (case `p<q`) representing the interior product `\iota_A B`, where
          `A` is ``self``

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorField.interior_product`
            for the interior product of a multivector field with a differential
            form

        EXAMPLES:

        Interior product of a 1-form with a 2-vector field on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1) # the sphere S^2
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() # stereographic coord. North
            sage: c_uv.<u,v> = V.chart() # stereographic coord. South
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: a = M.one_form({e_xy: [y, x]}, name='a')
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: b = M.multivector_field(2, name='b')
            sage: b[e_xy,1,2] = x*y
            sage: b.add_comp_by_continuation(e_uv, W, c_uv)
            sage: s = a.interior_product(b); s
            Vector field i_a b on the 2-dimensional differentiable manifold S^2
            sage: s.display(e_xy)
            i_a b = -x^2*y ∂/∂x + x*y^2 ∂/∂y
            sage: s.display(e_uv)
            i_a b = (u^4*v - 3*u^2*v^3)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) ∂/∂u
             + (3*u^3*v^2 - u*v^4)/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6) ∂/∂v
            sage: s == a.contract(b)
            True

        Interior product of a 2-form with a 2-vector field::

            sage: a = M.diff_form(2, name='a')
            sage: a[e_xy,1,2] = 4/(x^2+y^2+1)^2   # the standard area 2-form
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: s = a.interior_product(b); s
            Scalar field i_a b on the 2-dimensional differentiable manifold S^2
            sage: s.display()
            i_a b: S^2 → ℝ
            on U: (x, y) ↦ 8*x*y/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
            on V: (u, v) ↦ 8*u*v/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)

        Some checks::

            sage: s == a.contract(0, 1, b, 0, 1)
            True
            sage: s.restrict(U) == 2 * a[[e_xy,1,2]] * b[[e_xy,1,2]]
            True
            sage: s.restrict(V) == 2 * a[[e_uv,1,2]] * b[[e_uv,1,2]]
            True
        """
        from sage.tensor.modules.format_utilities import is_atomic
        if self._domain.is_subset(qvect._domain):
            if not self._ambient_domain.is_subset(qvect._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        elif qvect._domain.is_subset(self._domain):
            if not qvect._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        dom_resu = self._domain.intersection(qvect._domain)
        ambient_dom_resu = self._ambient_domain.intersection(qvect._ambient_domain)
        self_r = self.restrict(dom_resu)
        qvect_r = qvect.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleAltForm version:
            return FreeModuleAltForm.interior_product(self_r, qvect_r)
        # Otherwise, the result is created here:
        # Name of the result
        resu_name = None
        if self._name is not None and qvect._name is not None:
            sname = self._name
            oname = qvect._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu_name = 'i_' + sname + ' ' + oname
        resu_latex_name = None
        if self._latex_name is not None and qvect._latex_name is not None:
            slname = self._latex_name
            olname = qvect._latex_name
            if not is_atomic(olname):
                olname = r'\left(' + olname + r'\right)'
            resu_latex_name = r'\iota_{' + slname + '} ' + olname
        # Domain and computation of the result
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        resu_degree = qvect._tensor_rank - self._tensor_rank
        resu = vmodule.alternating_contravariant_tensor(resu_degree,
                                    name=resu_name, latex_name=resu_latex_name)
        for dom in self_r._restrictions:
            if dom in qvect_r._restrictions:
                resu._restrictions[dom] = \
                    self_r._restrictions[dom].interior_product(
                                                    qvect_r._restrictions[dom])
        if resu_degree == 0:
            if not resu._express:  # only the restrictions to subdomains have
                                   # been initialized
                for chart in dom_resu.top_charts():
                    resu._express[chart] = \
                            resu.restrict(chart.domain()).coord_function(chart)
        return resu

# *****************************************************************************


class DiffFormParal(FreeModuleAltForm, TensorFieldParal, DiffForm):
    r"""
    Differential form with values on a parallelizable manifold.

    Given a differentiable manifold `U`,  a differentiable map
    `\Phi: U \rightarrow M` to a parallelizable manifold `M` and a positive
    integer `p`, a *differential form of degree* `p` (or `p`-*form*)
    *along* `U` *with values on* `M\supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(0,p)}M

    (`T^{(0,p)}M` being the tensor bundle of type `(0,p)` over `M`) such that

    .. MATH::

        \forall x \in U,\quad a(x) \in \Lambda^p(T_{\Phi(x)}^* M) ,

    where `T_{\Phi(x)}^* M` is the dual of the tangent space to `M` at
    `\Phi(x)` and `\Lambda^p` stands for the exterior power of degree `p` (cf.
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerDualFreeModule`).
    In other words, `a(x)` is an alternating multilinear form of degree `p` of
    the tangent vector space `T_{\Phi(x)} M`.

    The standard case of a differential form *on* a manifold `M` corresponds
    to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
    being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, the class :class:`DiffForm` must
        be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi`
    - ``degree`` -- the degree of the differential form (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the differential form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      differential form; if none is provided, the LaTeX symbol is set to
      ``name``

    EXAMPLES:

    A 2-form on a 4-dimensional manifold::

        sage: M = Manifold(4, 'M')
        sage: c_txyz.<t,x,y,z> = M.chart()
        sage: a = M.diff_form(2, name='a') ; a
        2-form a on the 4-dimensional differentiable manifold M
        sage: a.parent()
        Free module Omega^2(M) of 2-forms on the 4-dimensional differentiable
         manifold M

    A differential form is a tensor field of purely covariant type::

        sage: a.tensor_type()
        (0, 2)

    It is antisymmetric, its components being
    :class:`~sage.tensor.modules.comp.CompFullyAntiSym`::

        sage: a.symmetries()
        no symmetry;  antisymmetry: (0, 1)
        sage: a[0,1] = 2
        sage: a[1,0]
        -2
        sage: a.comp()
        Fully antisymmetric 2-indices components w.r.t. Coordinate frame (M, (∂/∂t,∂/∂x,∂/∂y,∂/∂z))
        sage: type(a.comp())
        <class 'sage.tensor.modules.comp.CompFullyAntiSym'>

    Setting a component with repeated indices to a nonzero value
    results in an error::

        sage: a[1,1] = 3
        Traceback (most recent call last):
        ...
        ValueError: by antisymmetry, the component cannot have a nonzero value
         for the indices (1, 1)
        sage: a[1,1] = 0  # OK, albeit useless
        sage: a[1,2] = 3  # OK

    The expansion of a differential form with respect to a given coframe is
    displayed via the method
    :meth:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.display`::

        sage: a.display() # expansion with respect to the default coframe (dt, dx, dy, dz)
        a = 2 dt∧dx + 3 dx∧dy
        sage: latex(a.display()) # output for the notebook
        a = 2 \mathrm{d} t\wedge \mathrm{d} x
         + 3 \mathrm{d} x\wedge \mathrm{d} y

    Differential forms can be added or subtracted::

        sage: b = M.diff_form(2)
        sage: b[0,1], b[0,2], b[0,3] = (1,2,3)
        sage: s = a + b ; s
        2-form on the 4-dimensional differentiable manifold M
        sage: a[:], b[:], s[:]
        (
        [ 0  2  0  0]  [ 0  1  2  3]  [ 0  3  2  3]
        [-2  0  3  0]  [-1  0  0  0]  [-3  0  3  0]
        [ 0 -3  0  0]  [-2  0  0  0]  [-2 -3  0  0]
        [ 0  0  0  0], [-3  0  0  0], [-3  0  0  0]
        )
        sage: s = a - b ; s
        2-form on the 4-dimensional differentiable manifold M
        sage: s[:]
        [ 0  1 -2 -3]
        [-1  0  3  0]
        [ 2 -3  0  0]
        [ 3  0  0  0]

    An example of 3-form is the volume element on `\RR^3` in Cartesian
    coordinates::

        sage: M = Manifold(3, 'R3', latex_name=r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()
        sage: eps = M.diff_form(3, name='epsilon', latex_name=r'\epsilon')
        sage: eps[1,2,3] = 1  # the only independent component
        sage: eps[:] # all the components are set from the previous line:
        [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]],
         [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
        sage: eps.display()
        epsilon = dx∧dy∧dz

    Spherical components of the volume element from the tensorial
    change-of-frame formula::

        sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: spher_to_cart = c_spher.transition_map(c_cart,
        ....:                [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
        sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2),
        ....:                              atan2(sqrt(x^2+y^2),z), atan2(y, x))
        Check of the inverse coordinate transformation:
          r == r  *passed*
          th == arctan2(r*sin(th), r*cos(th))  **failed**
          ph == arctan2(r*sin(ph)*sin(th), r*cos(ph)*sin(th))  **failed**
          x == x  *passed*
          y == y  *passed*
          z == z  *passed*
        NB: a failed report can reflect a mere lack of simplification.
        sage: eps.comp(c_spher.frame()) # computation of the components in the spherical frame
        Fully antisymmetric 3-indices components w.r.t. Coordinate frame
         (R3, (∂/∂r,∂/∂th,∂/∂ph))
        sage: eps.comp(c_spher.frame())[1,2,3, c_spher]
        r^2*sin(th)
        sage: eps.display(c_spher.frame())
        epsilon = sqrt(x^2 + y^2 + z^2)*sqrt(x^2 + y^2) dr∧dth∧dph
        sage: eps.display(c_spher.frame(), c_spher)
        epsilon = r^2*sin(th) dr∧dth∧dph

    As a shortcut of the above command, on can pass just the chart ``c_spher``
    to ``display``, the vector frame being then assumed to be the coordinate
    frame associated with the chart::

        sage: eps.display(c_spher)
        epsilon = r^2*sin(th) dr∧dth∧dph


    The exterior product of two differential forms is performed via the method
    :meth:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.wedge`::

        sage: a = M.one_form(x*y*z, -z*x, y*z, name='A')
        sage: b = M.one_form(cos(z), sin(x), cos(y), name='B')
        sage: ab = a.wedge(b) ; ab
        2-form A∧B on the 3-dimensional differentiable manifold R3
        sage: ab[:]
        [                         0  x*y*z*sin(x) + x*z*cos(z)  x*y*z*cos(y) - y*z*cos(z)]
        [-x*y*z*sin(x) - x*z*cos(z)                          0   -(x*cos(y) + y*sin(x))*z]
        [-x*y*z*cos(y) + y*z*cos(z)    (x*cos(y) + y*sin(x))*z                          0]
        sage: ab.display()
        A∧B = (x*y*z*sin(x) + x*z*cos(z)) dx∧dy + (x*y*z*cos(y) - y*z*cos(z)) dx∧dz
         - (x*cos(y) + y*sin(x))*z dy∧dz

    Let us check the formula relating the exterior product to the tensor
    product for 1-forms::

        sage: a.wedge(b) == a*b - b*a
        True

    The tensor product of a 1-form and a 2-form is not a 3-form but a tensor
    field of type `(0,3)` with less symmetries::

        sage: c = a*ab ; c
        Tensor field A⊗(A∧B) of type (0,3) on the 3-dimensional differentiable
         manifold R3
        sage: c.symmetries()  # the antisymmetry is only w.r.t. the last 2 arguments:
        no symmetry;  antisymmetry: (1, 2)
        sage: d = ab*a ; d
        Tensor field (A∧B)⊗A of type (0,3) on the 3-dimensional differentiable
         manifold R3
        sage: d.symmetries()  # the antisymmetry is only w.r.t. the first 2 arguments:
        no symmetry;  antisymmetry: (0, 1)

    The exterior derivative of a differential form is obtained by means
    of the method :meth:`exterior_derivative`::

        sage: da = a.exterior_derivative() ; da
        2-form dA on the 3-dimensional differentiable manifold R3
        sage: da.display()
        dA = -(x + 1)*z dx∧dy - x*y dx∧dz + (x + z) dy∧dz
        sage: db = b.exterior_derivative() ; db
        2-form dB on the 3-dimensional differentiable manifold R3
        sage: db.display()
        dB = cos(x) dx∧dy + sin(z) dx∧dz - sin(y) dy∧dz
        sage: dab = ab.exterior_derivative() ; dab
        3-form d(A∧B) on the 3-dimensional differentiable manifold R3

    or by applying the function ``diff`` to the differential form::

        sage: diff(a) is a.exterior_derivative()
        True

    As a 3-form over a 3-dimensional manifold, ``d(A∧B)`` is necessarily
    proportional to the volume 3-form::

        sage: dab == dab[[1,2,3]]/eps[[1,2,3]]*eps
        True

    We may also check that the classical anti-derivation formula is fulfilled::

        sage: dab == da.wedge(b) - a.wedge(db)
        True

    The Lie derivative of a 2-form is a 2-form::

        sage: v = M.vector_field(y*z, -x*z, x*y, name='v')
        sage: ab.lie_der(v)  # long time
        2-form on the 3-dimensional differentiable manifold R3

    Let us check Cartan formula, which expresses the Lie derivative in terms
    of exterior derivatives::

        sage: ab.lie_der(v) == (v.contract(ab.exterior_derivative())  # long time
        ....:                   + v.contract(ab).exterior_derivative())
        True

    A 1-form on a `\RR^3`::

        sage: om = M.one_form(name='omega', latex_name=r'\omega'); om
        1-form omega on the 3-dimensional differentiable manifold R3

    A 1-form is of course a differential form::

        sage: isinstance(om, sage.manifolds.differentiable.diff_form.DiffFormParal)
        True
        sage: om.parent()
        Free module Omega^1(R3) of 1-forms on the 3-dimensional differentiable
         manifold R3
        sage: om.tensor_type()
        (0, 1)

    Setting the components with respect to the manifold's default frame::

        sage: om[:] = (2*z, x, x-y)
        sage: om[:]
        [2*z, x, x - y]
        sage: om.display()
        omega = 2*z dx + x dy + (x - y) dz

    A 1-form acts on vector fields::

        sage: v = M.vector_field(x, 2*y, 3*z, name='V')
        sage: om(v)
        Scalar field omega(V) on the 3-dimensional differentiable manifold R3
        sage: om(v).display()
        omega(V): R3 → ℝ
           (x, y, z) ↦ 2*x*y + (5*x - 3*y)*z
           (r, th, ph) ↦ 2*r^2*cos(ph)*sin(ph)*sin(th)^2 + r^2*(5*cos(ph)
                            - 3*sin(ph))*cos(th)*sin(th)
        sage: latex(om(v))
        \omega\left(V\right)

    The tensor product of two 1-forms is a tensor field of type `(0,2)`::

        sage: a = M.one_form(1, 2, 3, name='A')
        sage: b = M.one_form(6, 5, 4, name='B')
        sage: c = a*b ; c
        Tensor field A⊗B of type (0,2) on the 3-dimensional differentiable
         manifold R3
        sage: c[:]
        [ 6  5  4]
        [12 10  8]
        [18 15 12]
        sage: c.symmetries()    # c has no symmetries:
        no symmetry;  no antisymmetry
    """
    def __init__(self, vector_field_module: VectorFieldModule, degree: int, name: Optional[str] = None,
                 latex_name: Optional[str] = None):
        r"""
        Construct a differential form.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``DiffFormParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: A = M.diff_form_module(2)
            sage: XM = M.vector_field_module()
            sage: a = A.element_class(XM, 2, name='a'); a
            2-form a on the 2-dimensional differentiable manifold M
            sage: a[0,1] = x*y
            sage: TestSuite(a).run()

        Construction via ``DifferentiableManifold.diff_form``::

            sage: a1 = M.diff_form(2, name='a'); a1
            2-form a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        Initializing the components at the construction::

            sage: a = M.diff_form(2, [[0, x*y], [-x*y, 0]], name='a')
            sage: a.display()
            a = x*y dx∧dy
        """
        FreeModuleAltForm.__init__(self, vector_field_module, degree,
                                   name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.diff_form(2, name='a')
            sage: a._repr_()
            '2-form a on the 3-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            '2-form a on the 3-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            2-form a on the 3-dimensional differentiable manifold M
            sage: b = M.diff_form(2)
            sage: b._repr_()
            '2-form on the 3-dimensional differentiable manifold M'
        """
        description = "{}-form ".format(self._tensor_rank)
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.diff_form(2, name='a')
            sage: a1 = a._new_instance(); a1
            2-form on the 3-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True
        """
        return type(self)(self._fmodule, self._tensor_rank)

    # This method is needed to redirect to the correct class (TensorFieldParal)
    def _init_derived(self):
        r"""
        Initialize the derived quantities of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.diff_form(2, name='a')
            sage: a._init_derived()
        """
        TensorFieldParal._init_derived(self)

    def _del_derived(self, del_restrictions: bool = True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- boolean (default: ``True``); determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.diff_form(2, name='a')
            sage: a._del_derived()
        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)
        self.exterior_derivative.clear_cache()

    def __call__(self, *args):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.__call__`
        to allow for domain treatment.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.diff_form(2, name='a')
            sage: a[0,1] = x*y
            sage: a.display()
            a = x*y dx∧dy
            sage: u = M.vector_field(1+x, 2-y, name='u')
            sage: v = M.vector_field(-y, x, name='v')
            sage: s = a.__call__(u,v); s
            Scalar field a(u,v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(u,v): M → ℝ
               (x, y) ↦ -x*y^3 + 2*x*y^2 + (x^3 + x^2)*y
            sage: s == a[[0,1]]*(u[[0]]*v[[1]] - u[[1]]*v[[0]])
            True
            sage: s == a(u,v)  # indirect doctest
            True
        """
        return TensorFieldParal.__call__(self, *args)

    @cached_method
    def exterior_derivative(self) -> DiffFormParal:
        r"""
        Compute the exterior derivative of ``self``.

        OUTPUT:

        - a :class:`DiffFormParal` representing the exterior
          derivative of the differential form

        EXAMPLES:

        Exterior derivative of a 1-form on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: c_txyz.<t,x,y,z> = M.chart()
            sage: a = M.one_form(t*x*y*z, z*y**2, x*z**2, x**2 + y**2, name='A')
            sage: da = a.exterior_derivative() ; da
            2-form dA on the 4-dimensional differentiable manifold M
            sage: da.display()
            dA = -t*y*z dt∧dx - t*x*z dt∧dy - t*x*y dt∧dz
             + (-2*y*z + z^2) dx∧dy + (-y^2 + 2*x) dx∧dz
             + (-2*x*z + 2*y) dy∧dz
            sage: latex(da)
            \mathrm{d}A

        The result is cached, i.e. is not recomputed unless ``a`` is changed::

            sage: a.exterior_derivative() is da
            True

        Instead of invoking the method :meth:`exterior_derivative`, one may
        use the global function ``diff``::

            sage: diff(a) is a.exterior_derivative()
            True

        The exterior derivative is nilpotent::

            sage: dda = da.exterior_derivative() ; dda
            3-form ddA on the 4-dimensional differentiable manifold M
            sage: dda.display()
            ddA = 0
            sage: dda == 0
            True

        Let us check Cartan's identity::

            sage: v = M.vector_field(-y, x, t, z, name='v')
            sage: a.lie_der(v) == v.contract(diff(a)) + diff(a(v)) # long time
            True
        """
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        from sage.tensor.modules.comp import CompFullyAntiSym
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        fmodule = self._fmodule # shortcut
        rname = format_unop_txt('d', self._name)
        rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
        resu = fmodule.alternating_form(self._tensor_rank + 1,
                                        name=rname,
                                        latex_name=rlname)
        # 1/ List of all coordinate frames in which the components of self
        # are known
        coord_frames = []
        for frame in self._components:
            if isinstance(frame, CoordFrame):
                coord_frames.append(frame)
        if not coord_frames:
            # A coordinate frame is searched, at the price of a change of
            # frame, privileging the frame of the domain's default chart
            dom = self._domain
            def_coordf = dom._def_chart._frame
            for frame in self._components:
                if (frame, def_coordf) in dom._frame_changes:
                    self.comp(def_coordf, from_basis=frame)
                    coord_frames = [def_coordf]
                    break
            if not coord_frames:
                for chart in dom._atlas:
                    if chart != dom._def_chart: # the case def_chart is
                                                # treated above
                        coordf = chart._frame
                        for frame in self._components:
                            if (frame, coordf) in dom._frame_changes:
                                self.comp(coordf, from_basis=frame)
                                coord_frames[coordf]
                                break
                        if coord_frames:
                            break
        # 2/ The computation:
        for frame in coord_frames:
            chart = frame._chart
            sc = self._components[frame]
            dc = CompFullyAntiSym(fmodule._ring, frame,
                                  self._tensor_rank+1,
                                  start_index=fmodule._sindex,
                                output_formatter=fmodule._output_formatter)
            for ind, val in sc._comp.items():
                for i in fmodule.irange():
                    ind_d = (i,) + ind
                    if len(ind_d) == len(set(ind_d)):
                        # all indices are different
                        dc[[ind_d]] += \
                           val.coord_function(chart).diff(i).scalar_field()
            resu._components[frame] = dc
        return resu

    derivative = exterior_derivative  # allows one to use functional notation,
                                      # e.g. diff(a) for a.exterior_derivative()

    def wedge(self, other):
        r"""
        Exterior product of ``self`` with another differential form.

        INPUT:

        - ``other`` -- another differential form

        OUTPUT:

        - instance of :class:`DiffFormParal` representing the exterior
          product ``self ∧ other``

        EXAMPLES:

        Exterior product of a 1-form and a 2-form on a 3-dimensional
        manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: a = M.one_form(2, 1+x, y*z, name='a')
            sage: b = M.diff_form(2, name='b')
            sage: b[1,2], b[1,3], b[2,3] = y^2, z+x, z^2
            sage: a.display()
            a = 2 dx + (x + 1) dy + y*z dz
            sage: b.display()
            b = y^2 dx∧dy + (x + z) dx∧dz + z^2 dy∧dz
            sage: s = a.wedge(b); s
            3-form a∧b on the 3-dimensional differentiable manifold M
            sage: s.display()
            a∧b = (-x^2 + (y^3 - x - 1)*z + 2*z^2 - x) dx∧dy∧dz

        Check::

            sage: s[1,2,3] == a[1]*b[2,3] + a[2]*b[3,1] + a[3]*b[1,2]
            True

        Wedging with scalar fields yields the multiplication from right::

            sage: f = M.scalar_field(x, name='f')
            sage: t = a.wedge(f)
            sage: t.display()
            f*a = 2*x dx + (x^2 + x) dy + x*y*z dz
        """
        if other._tensor_rank == 0:
            return self * other
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        dom_resu = self._domain.intersection(other._domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        return FreeModuleAltForm.wedge(self_r, other_r)

    def interior_product(self, qvect):
        r"""
        Interior product with a multivector field.

        If ``self`` is a differential form `A` of degree `p` and `B` is a
        multivector field of degree `q\geq p` on the same manifold, the
        interior product of `A` by `B` is the multivector field `\iota_A B` of
        degree `q-p` defined by

        .. MATH::

            (\iota_A B)^{i_1\ldots i_{q-p}} = A_{k_1\ldots k_p}
                            B^{k_1\ldots k_p i_1\ldots i_{q-p}}

        .. NOTE::

            ``A.interior_product(B)`` yields the same result as
            ``A.contract(0,..., p-1, B, 0,..., p-1)`` (cf.
            :meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.contract`),
            but ``interior_product`` is more efficient, the alternating
            character of `A` being not used to reduce the computation in
            :meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.contract`

        INPUT:

        - ``qvect`` -- multivector field `B` (instance of
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal`);
          the degree of `B` must be at least equal to the degree of
          ``self``

        OUTPUT:

        - scalar field (case `p=q`) or
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal`
          (case `p<q`) representing the interior product `\iota_A B`,
          where `A` is ``self``

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal.interior_product`
            for the interior product of a multivector field with a
            differential form

        EXAMPLES:

        Interior product of a 1-form with a 2-vector field on a
        3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: a = M.one_form(2, 1+x, y*z, name='a')
            sage: b = M.multivector_field(2, name='b')
            sage: b[1,2], b[1,3], b[2,3] = y^2, z+x, -z^2
            sage: s = a.interior_product(b); s
            Vector field i_a b on the 3-dimensional differentiable
             manifold M
            sage: s.display()
            i_a b = (-(x + 1)*y^2 - x*y*z - y*z^2) ∂/∂x
             + (y*z^3 + 2*y^2) ∂/∂y + (-(x + 1)*z^2 + 2*x + 2*z) ∂/∂z
            sage: s == a.contract(b)
            True

        Interior product of a 2-form with a 2-vector field::

            sage: a = M.diff_form(2, name='a')
            sage: a[1,2], a[1,3], a[2,3] = x*y, -3, z
            sage: s = a.interior_product(b); s
            Scalar field i_a b on the 3-dimensional differentiable manifold M
            sage: s.display()
            i_a b: M → ℝ
               (x, y, z) ↦ 2*x*y^3 - 2*z^3 - 6*x - 6*z
            sage: s == a.contract(0,1,b,0,1)
            True

        TESTS:

        Check that :issue:`33780` is fixed::

            sage: v = X.frame()[1]  # vector field d/dx
            sage: f = X.coframe()[2]  # 1-form dy
            sage: f.interior_product(v)
            Scalar field zero on the 3-dimensional differentiable manifold M
        """
        if self._domain.is_subset(qvect._domain):
            if not self._ambient_domain.is_subset(qvect._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        elif qvect._domain.is_subset(self._domain):
            if not qvect._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        dom_resu = self._domain.intersection(qvect._domain)
        self_r = self.restrict(dom_resu)
        qvect_r = qvect.restrict(dom_resu)
        return FreeModuleAltForm.interior_product(self_r, qvect_r)
