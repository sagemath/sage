r"""
Tensor Field Modules

The set of tensor fields along a differentiable manifold `U` with values on
a differentiable manifold `M` via a differentiable map `\Phi: U \rightarrow M`
(possibly `U = M` and `\Phi = \mathrm{Id}_M`) is a module over the algebra
`C^k(U)` of differentiable scalar fields on `U`. It is a free module if
and only if `M` is parallelizable. Accordingly, two classes are devoted
to tensor field modules:

- :class:`TensorFieldModule` for tensor fields with values on a generic (in
  practice, not parallelizable) differentiable manifold `M`,
- :class:`TensorFieldFreeModule` for tensor fields with values on a
  parallelizable manifold `M`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_
"""

# *****************************************************************************
#       Copyright (C) 2015-2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#                     2015      Michal Bejger <bejger@camk.edu.pl>
#                     2016      Travis Scrimshaw <tscrimsh@umn.edu>
#                     2020      Michael Jung
#                     2022      Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.tensor.modules.reflexive_module import ReflexiveModule_tensor
from sage.tensor.modules.tensor_free_module import TensorFreeModule
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.manifolds.differentiable.diff_form import (DiffForm,
                                                     DiffFormParal)
from sage.manifolds.differentiable.multivectorfield import (MultivectorField,
                                                            MultivectorFieldParal)
from sage.manifolds.differentiable.automorphismfield import (AutomorphismField,
                                                             AutomorphismFieldParal)


class TensorFieldModule(UniqueRepresentation, ReflexiveModule_tensor):
    r"""
    Module of tensor fields of a given type `(k,l)` along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `U \rightarrow M`.

    Given two nonnegative integers `k` and `l` and a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    the *tensor field module* `T^{(k,l)}(U,\Phi)` is the set of all tensor
    fields of the type

    .. MATH::

        t:\ U  \longrightarrow T^{(k,l)} M

    (where `T^{(k,l)} M` is the tensor bundle of type `(k,l)` over `M`) such
    that

    .. MATH::

        t(p) \in T^{(k,l)}(T_{\Phi(p)}M)

    for all `p \in U`, i.e. `t(p)` is a tensor of type `(k,l)` on the
    tangent vector space `T_{\Phi(p)} M`. The set `T^{(k,l)}(U,\Phi)`
    is a module over `C^k(U)`, the ring (algebra) of differentiable
    scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).

    The standard case of tensor fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`; we then denote
    `T^{(k,l)}(M,\mathrm{Id}_M)` by merely `T^{(k,l)}(M)`. Other common
    cases are `\Phi` being an immersion and `\Phi` being a curve in `M`
    (`U` is then an open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, the class :class:`TensorFieldFreeModule`
        should be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow M`
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant
      rank and `l` the covariant rank

    EXAMPLES:

    Module of type-`(2,0)` tensor fields on the 2-sphere::

        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: T20 = M.tensor_field_module((2,0)); T20
        Module T^(2,0)(M) of type-(2,0) tensors fields on the 2-dimensional
         differentiable manifold M

    `T^{(2,0)}(M)` is a module over the algebra `C^k(M)`::

        sage: T20.category()
        Category of tensor products of modules over Algebra of differentiable scalar fields
         on the 2-dimensional differentiable manifold M
        sage: T20.base_ring() is M.scalar_field_algebra()
        True

    `T^{(2,0)}(M)` is not a free module::

        sage: from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule_abstract
        sage: isinstance(T20, FiniteRankFreeModule_abstract)
        False

    because `M = S^2` is not parallelizable::

        sage: M.is_manifestly_parallelizable()
        False

    On the contrary, the module of type-`(2,0)` tensor fields on `U` is a
    free module, since `U` is parallelizable (being a coordinate domain)::

        sage: T20U = U.tensor_field_module((2,0))
        sage: isinstance(T20U, FiniteRankFreeModule_abstract)
        True
        sage: U.is_manifestly_parallelizable()
        True

    The zero element::

        sage: z = T20.zero() ; z
        Tensor field zero of type (2,0) on the 2-dimensional differentiable
         manifold M
        sage: z is T20(0)
        True
        sage: z[c_xy.frame(),:]
        [0 0]
        [0 0]
        sage: z[c_uv.frame(),:]
        [0 0]
        [0 0]

    The module `T^{(2,0)}(M)` coerces to any module of type-`(2,0)` tensor
    fields defined on some subdomain of `M`, for instance `T^{(2,0)}(U)`::

        sage: T20U.has_coerce_map_from(T20)
        True

    The reverse is not true::

        sage: T20.has_coerce_map_from(T20U)
        False

    The coercion::

        sage: T20U.coerce_map_from(T20)
        Coercion map:
          From: Module T^(2,0)(M) of type-(2,0) tensors fields on the 2-dimensional differentiable manifold M
          To:   Free module T^(2,0)(U) of type-(2,0) tensors fields on the Open subset U of the 2-dimensional differentiable manifold M

    The coercion map is actually the *restriction* of tensor fields defined
    on `M` to `U`::

        sage: t = M.tensor_field(2,0, name='t')
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: t[eU,:] = [[2,0], [0,-3]]
        sage: t.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: T20U(t)  # the conversion map in action
        Tensor field t of type (2,0) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: T20U(t) is t.restrict(U)
        True

    There is also a coercion map from fields of tangent-space automorphisms to
    tensor fields of type-`(1,1)`::

        sage: T11 = M.tensor_field_module((1,1)) ; T11
        Module T^(1,1)(M) of type-(1,1) tensors fields on the 2-dimensional
         differentiable manifold M
        sage: GL = M.automorphism_field_group() ; GL
        General linear group of the Module X(M) of vector fields on the
         2-dimensional differentiable manifold M
        sage: T11.has_coerce_map_from(GL)
        True

    Explicit call to the coercion map::

        sage: a = GL.one() ; a
        Field of tangent-space identity maps on the 2-dimensional
         differentiable manifold M
        sage: a.parent()
        General linear group of the Module X(M) of vector fields on the
         2-dimensional differentiable manifold M
        sage: ta = T11.coerce(a) ; ta
        Tensor field Id of type (1,1) on the 2-dimensional differentiable
         manifold M
        sage: ta.parent()
        Module T^(1,1)(M) of type-(1,1) tensors fields on the 2-dimensional
         differentiable manifold M
        sage: ta[eU,:]  # ta on U
        [1 0]
        [0 1]
        sage: ta[eV,:]  # ta on V
        [1 0]
        [0 1]

    TESTS::

        sage: T11.tensor_factors()
        [Module X(M) of vector fields on the 2-dimensional differentiable manifold M,
        Module Omega^1(M) of 1-forms on the 2-dimensional differentiable manifold M]
    """
    Element = TensorField

    def __init__(self, vector_field_module, tensor_type, category=None):
        r"""
        Construct a module of tensor fields taking values on a (a priori) not
        parallelizable differentiable manifold.

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: XM = M.vector_field_module()
            sage: from sage.manifolds.differentiable.tensorfield_module import TensorFieldModule
            sage: T21 = TensorFieldModule(XM, (2,1)); T21
            Module T^(2,1)(M) of type-(2,1) tensors fields on the 2-dimensional
             differentiable manifold M
            sage: T21 is M.tensor_field_module((2,1))
            True
            sage: TestSuite(T21).run(skip='_test_elements')

        In the above test suite, ``_test_elements`` is skipped because of the
        ``_test_pickling`` error of the elements (to be fixed in
        :class:`~sage.manifolds.differentiable.tensorfield.TensorField`)
        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "T^({},{})({}".format(kcon, lcov, domain._name)
        latex_name = r"\mathcal{{T}}^{{({},{})}}\left({}".format(kcon, lcov, domain._latex_name)
        if dest_map is not domain.identity_map():
            dm_name = dest_map._name
            dm_latex_name = dest_map._latex_name
            if dm_name is None:
                dm_name = "unnamed map"
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            name += "," + dm_name
            latex_name += "," + dm_latex_name
        self._name = name + ")"
        self._latex_name = latex_name + r"\right)"
        self._vmodule = vector_field_module
        self._tensor_type = tensor_type
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        category = Modules(self._ring).TensorProducts().or_subcategory(category)
        Parent.__init__(self, base=self._ring, category=category)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: T20 = M.tensor_field_module((2,0))
            sage: t = T20([[1+x, 2], [x*y, 3-y]], name='t'); t
            Tensor field t of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: t.display(c_xy.frame())
            t = (x + 1) ∂/∂x⊗∂/∂x + 2 ∂/∂x⊗∂/∂y + x*y ∂/∂y⊗∂/∂x
             + (-y + 3) ∂/∂y⊗∂/∂y
            sage: T20(0) is T20.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, DiffForm):
            # coercion of a p-form to a type-(0,p) tensor field:
            form = comp # for readability
            p = form.degree()
            if (self._tensor_type != (0,p) or
                self._vmodule != form.base_module()):
                raise TypeError("cannot convert the {}".format(form) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._vmodule, (0,p),
                                      name=form._name,
                                      latex_name=form._latex_name,
                                      antisym=asym)
            for dom, rst in form._restrictions.items():
                resu._restrictions[dom] = dom.tensor_field_module((0,p))(rst)
            return resu
        if isinstance(comp, MultivectorField):
            # coercion of a p-vector field to a type-(p,0) tensor:
            pvect = comp # for readability
            p = pvect.degree()
            if (self._tensor_type != (p,0) or
                self._vmodule != pvect.base_module()):
                raise TypeError("cannot convert the {}".format(pvect) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._vmodule, (p,0),
                                      name=pvect._name,
                                      latex_name=pvect._latex_name,
                                      antisym=asym)
            for dom, rst in pvect._restrictions.items():
                resu._restrictions[dom] = dom.tensor_field_module((p,0))(rst)
            return resu
        if isinstance(comp, AutomorphismField):
            # coercion of an automorphism to a type-(1,1) tensor:
            autom = comp # for readability
            if (self._tensor_type != (1,1) or
                self._vmodule != autom.base_module()):
                raise TypeError("cannot convert the {}".format(autom) +
                                " to an element of {}".format(self))
            resu = self.element_class(self._vmodule, (1,1),
                                      name=autom._name,
                                      latex_name=autom._latex_name)
            for dom, rest in autom._restrictions.items():
                resu._restrictions[dom] = dom.tensor_field_module((1,1))(rest)
            return resu
        if isinstance(comp, TensorField):
            # coercion by domain restriction
            if (self._tensor_type == comp._tensor_type
                and self._domain.is_subset(comp._domain)
                and self._ambient_domain.is_subset(comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {}".format(comp) +
                                " to an element of {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._vmodule, self._tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) tensor field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: T31 = M.tensor_field_module((3,1))
            sage: T31._an_element_()
            Tensor field of type (3,1) on the 2-dimensional differentiable
             manifold M
        """
        resu = self.element_class(self._vmodule, self._tensor_type)
        for oc in self._domain.open_covers(trivial=False):
            # the first non-trivial open cover is selected
            for dom in oc:
                vmodule_dom = dom.vector_field_module(dest_map=self._dest_map.restrict(dom))
                tmodule_dom = vmodule_dom.tensor_module(*(self._tensor_type))
                resu.set_restriction(tmodule_dom._an_element_())
            return resu
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: T02 = M.tensor_field_module((0,2))
            sage: T02U = U.tensor_field_module((0,2))
            sage: T02U._coerce_map_from_(T02)
            True
            sage: T02._coerce_map_from_(T02U)
            False
            sage: T02._coerce_map_from_(M.diff_form_module(2))
            True
            sage: T20 = M.tensor_field_module((2,0))
            sage: T20._coerce_map_from_(M.multivector_module(2))
            True
            sage: T11 = M.tensor_field_module((1,1))
            sage: T11._coerce_map_from_(M.automorphism_field_group())
            True
        """
        from sage.manifolds.differentiable.diff_form_module import \
                                                          DiffFormModule
        from sage.manifolds.differentiable.multivector_module import \
                                                       MultivectorModule
        from sage.manifolds.differentiable.automorphismfield_group \
                                           import AutomorphismFieldGroup
        if isinstance(other, (TensorFieldModule, TensorFieldFreeModule)):
            # coercion by domain restriction
            return (self._tensor_type == other._tensor_type
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(other._ambient_domain))
        if isinstance(other, DiffFormModule):
            # coercion of p-forms to type-(0,p) tensor fields
            return (self._vmodule is other.base_module()
                    and self._tensor_type == (0, other.degree()))
        if isinstance(other, MultivectorModule):
            # coercion of p-vector fields to type-(p,0) tensor fields
            return (self._vmodule is other.base_module()
                    and self._tensor_type == (other.degree(),0))
        if isinstance(other, AutomorphismFieldGroup):
            # coercion of automorphism fields to type-(1,1) tensor fields
            return (self._vmodule is other.base_module()
                    and self._tensor_type == (1,1))
        return False

    #### End of parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: T13 = M.tensor_field_module((1,3))
            sage: T13._repr_()
            'Module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional differentiable manifold M'
            sage: repr(T13)  # indirect doctest
            'Module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional differentiable manifold M'
            sage: T13  # indirect doctest
            Module T^(1,3)(M) of type-(1,3) tensors fields on the 2-dimensional
             differentiable manifold M
        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-({},{})".format(self._tensor_type[0],
                                                self._tensor_type[1])
        description += " tensors fields "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {}".format(self._domain) + \
                           " mapped into the {}".format(self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: T13 = M.tensor_field_module((1,3))
            sage: T13._latex_()
            '\\mathcal{T}^{(1,3)}\\left(M\\right)'
            sage: latex(T13)  # indirect doctest
            \mathcal{T}^{(1,3)}\left(M\right)
        """
        if self._latex_name is None:
            return r'\text{' + str(self) + r'}'
        else:
            return self._latex_name

    def base_module(self):
        r"""
        Return the vector field module on which ``self`` is constructed.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
          representing the module on which ``self`` is defined

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: T13 = M.tensor_field_module((1,3))
            sage: T13.base_module()
            Module X(M) of vector fields on the 2-dimensional differentiable
             manifold M
            sage: T13.base_module() is M.vector_field_module()
            True
            sage: T13.base_module().base_ring()
            Algebra of differentiable scalar fields on the 2-dimensional
             differentiable manifold M
        """
        return self._vmodule

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``.

        OUTPUT: pair `(k,l)` of nonnegative integers such that the tensor
        fields belonging to this module are of type `(k,l)`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: T13 = M.tensor_field_module((1,3))
            sage: T13.tensor_type()
            (1, 3)
            sage: T20 = M.tensor_field_module((2,0))
            sage: T20.tensor_type()
            (2, 0)
        """
        return self._tensor_type

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: T20 = M.tensor_field_module((2,0))
            sage: T20.zero()
            Tensor field zero of type (2,0) on the
             2-dimensional differentiable manifold M
        """
        resu = self._element_constructor_(name='zero', latex_name='0')
        for frame in self._domain._frames:
            if self._dest_map.restrict(frame._domain) == frame._dest_map:
                resu.add_comp(frame)
                # (since new components are initialized to zero)
        resu._is_zero = True  # This element is certainly zero
        resu.set_immutable()
        return resu

#***********************************************************************


class TensorFieldFreeModule(TensorFreeModule):
    r"""
    Free module of tensor fields of a given type `(k,l)` along a
    differentiable manifold `U` with values on a parallelizable manifold `M`,
    via a differentiable map `U \rightarrow M`.

    Given two nonnegative integers `k` and `l` and a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    the *tensor field module* `T^{(k,l)}(U, \Phi)` is the set of all tensor
    fields of the type

    .. MATH::

        t:\ U \longrightarrow T^{(k,l)} M

    (where `T^{(k,l)}M` is the tensor bundle of type `(k,l)` over `M`)
    such that

    .. MATH::

        t(p) \in T^{(k,l)}(T_{\Phi(p)}M)

    for all `p \in U`, i.e. `t(p)` is a tensor of type `(k,l)` on the
    tangent vector space `T_{\Phi(p)}M`. Since `M` is parallelizable,
    the set `T^{(k,l)}(U,\Phi)` is a free module over `C^k(U)`, the
    ring (algebra) of differentiable scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).

    The standard case of tensor fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`; we then denote
    `T^{(k,l)}(M,\mathrm{Id}_M)` by merely `T^{(k,l)}(M)`. Other common cases
    are `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then
    an open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, the class :class:`TensorFieldModule`
        should be used instead, for `T^{(k,l)}(U,\Phi)` is no longer a
        free module.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow M`
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank
      and `l` the covariant rank

    EXAMPLES:

    Module of type-`(2,0)` tensor fields on `\RR^3`::

        sage: M = Manifold(3, 'R^3')
        sage: c_xyz.<x,y,z> = M.chart()  # Cartesian coordinates
        sage: T20 = M.tensor_field_module((2,0)) ; T20
        Free module T^(2,0)(R^3) of type-(2,0) tensors fields on the
         3-dimensional differentiable manifold R^3

    `T^{(2,0)}(\RR^3)` is a module over the algebra `C^k(\RR^3)`::

        sage: T20.category()
        Category of tensor products of finite dimensional modules over
         Algebra of differentiable scalar fields on the 3-dimensional differentiable manifold R^3
        sage: T20.base_ring() is M.scalar_field_algebra()
        True

    `T^{(2,0)}(\RR^3)` is a free module::

        sage: from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule_abstract
        sage: isinstance(T20, FiniteRankFreeModule_abstract)
        True

    because `M = \RR^3` is parallelizable::

        sage: M.is_manifestly_parallelizable()
        True

    The zero element::

        sage: z = T20.zero() ; z
        Tensor field zero of type (2,0) on the 3-dimensional differentiable
         manifold R^3
        sage: z[:]
        [0 0 0]
        [0 0 0]
        [0 0 0]

    A random element::

        sage: t = T20.an_element() ; t
        Tensor field of type (2,0) on the 3-dimensional differentiable
         manifold R^3
        sage: t[:]
        [2 0 0]
        [0 0 0]
        [0 0 0]

    The module `T^{(2,0)}(\RR^3)` coerces to any module of type-`(2,0)`
    tensor fields defined on some subdomain of `\RR^3`::

        sage: U = M.open_subset('U', coord_def={c_xyz: x>0})
        sage: T20U = U.tensor_field_module((2,0))
        sage: T20U.has_coerce_map_from(T20)
        True
        sage: T20.has_coerce_map_from(T20U)  # the reverse is not true
        False
        sage: T20U.coerce_map_from(T20)
        Coercion map:
          From: Free module T^(2,0)(R^3) of type-(2,0) tensors fields on the 3-dimensional differentiable manifold R^3
          To:   Free module T^(2,0)(U) of type-(2,0) tensors fields on the Open subset U of the 3-dimensional differentiable manifold R^3

    The coercion map is actually the *restriction* of tensor fields defined
    on `\RR^3` to `U`.

    There is also a coercion map from fields of tangent-space automorphisms to
    tensor fields of type `(1,1)`::

        sage: T11 = M.tensor_field_module((1,1)) ; T11
        Free module T^(1,1)(R^3) of type-(1,1) tensors fields on the
         3-dimensional differentiable manifold R^3
        sage: GL = M.automorphism_field_group() ; GL
        General linear group of the Free module X(R^3) of vector fields on the
         3-dimensional differentiable manifold R^3
        sage: T11.has_coerce_map_from(GL)
        True

    An explicit call to this coercion map is::

        sage: id = GL.one() ; id
        Field of tangent-space identity maps on the 3-dimensional
         differentiable manifold R^3
        sage: tid = T11(id) ; tid
        Tensor field Id of type (1,1) on the 3-dimensional differentiable
         manifold R^3
        sage: tid[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    Element = TensorFieldParal

    def __init__(self, vector_field_module, tensor_type):
        r"""
        Construct a module of tensor fields taking values on a
        parallelizable differentiable manifold.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: from sage.manifolds.differentiable.tensorfield_module import TensorFieldFreeModule
            sage: T12 = TensorFieldFreeModule(XM, (1,2)); T12
            Free module T^(1,2)(M) of type-(1,2) tensors fields on the
             2-dimensional differentiable manifold M
            sage: T12 is M.tensor_field_module((1,2))
            True
            sage: TestSuite(T12).run()
        """
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "T^({},{})({}".format(kcon, lcov, domain._name)
        latex_name = r"\mathcal{{T}}^{{({}, {})}}\left({}".format(kcon,
                                               lcov, domain._latex_name)
        if dest_map is not domain.identity_map():
            dm_name = dest_map._name
            dm_latex_name = dest_map._latex_name
            if dm_name is None:
                dm_name = "unnamed map"
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            name += "," + dm_name
            latex_name += "," + dm_latex_name
        name += ")"
        latex_name += r"\right)"
        TensorFreeModule.__init__(self, vector_field_module, tensor_type,
                                  name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: T12 = M.tensor_field_module((1,2))
            sage: t = T12([[[x,-y], [2,y]], [[1+x,y^2], [x^2,3]],
            ....:          [[x*y, 1-x], [y^2, x]]], name='t'); t
            Tensor field t of type (1,2) on the 2-dimensional
             differentiable manifold M
            sage: t.display()
            t = x ∂/∂x⊗dx⊗dx - y ∂/∂x⊗dx⊗dy + 2 ∂/∂x⊗dy⊗dx + y ∂/∂x⊗dy⊗dy
             + (x + 1) ∂/∂y⊗dx⊗dx + y^2 ∂/∂y⊗dx⊗dy + x^2 ∂/∂y⊗dy⊗dx
             + 3 ∂/∂y⊗dy⊗dy
            sage: T12(0) is T12.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, DiffFormParal):
            # coercion of a p-form to a type-(0,p) tensor field:
            form = comp # for readability
            p = form.degree()
            if (self._tensor_type != (0,p) or
                self._fmodule != form.base_module()):
                raise TypeError("cannot convert the {}".format(form) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._fmodule, (0,p),
                                      name=form._name,
                                      latex_name=form._latex_name,
                                      antisym=asym)
            for frame, cp in form._components.items():
                resu._components[frame] = cp.copy()
            return resu
        if isinstance(comp, MultivectorFieldParal):
            # coercion of a p-vector field to a type-(p,0) tensor field:
            pvect = comp # for readability
            p = pvect.degree()
            if (self._tensor_type != (p,0) or
                self._fmodule != pvect.base_module()):
                raise TypeError("cannot convert the {}".format(pvect) +
                                " to an element of {}".format(self))
            if p == 1:
                asym = None
            else:
                asym = range(p)
            resu = self.element_class(self._fmodule, (p,0),
                                      name=pvect._name,
                                      latex_name=pvect._latex_name,
                                      antisym=asym)
            for frame, cp in pvect._components.items():
                resu._components[frame] = cp.copy()
            return resu
        if isinstance(comp, AutomorphismFieldParal):
            # coercion of an automorphism to a type-(1,1) tensor:
            autom = comp # for readability
            if (self._tensor_type != (1,1) or
                self._fmodule != autom.base_module()):
                raise TypeError("cannot convert the {}".format(autom) +
                                " to an element of {}".format(self))
            resu = self.element_class(self._fmodule, (1,1),
                                      name=autom._name,
                                      latex_name=autom._latex_name)
            for basis, comp in autom._components.items():
                resu._components[basis] = comp.copy()
            return resu
        if isinstance(comp, TensorField):
            # coercion by domain restriction
            if (self._tensor_type == comp._tensor_type
                and self._domain.is_subset(comp._domain)
                and self._ambient_domain.is_subset(
                                                 comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot convert the {}".format(comp) +
                                " to an element of {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # Standard construction
        resu = self.element_class(self._fmodule, self._tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass TensorFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: T02 = M.tensor_field_module((0,2))
            sage: T02U = U.tensor_field_module((0,2))
            sage: T02U._coerce_map_from_(T02)
            True
            sage: T02._coerce_map_from_(T02U)
            False
            sage: T02._coerce_map_from_(M.diff_form_module(2))
            True
            sage: T20 = M.tensor_field_module((2,0))
            sage: T20._coerce_map_from_(M.multivector_module(2))
            True
            sage: T11 = M.tensor_field_module((1,1))
            sage: T11._coerce_map_from_(M.automorphism_field_group())
            True
        """
        from sage.manifolds.differentiable.diff_form_module import \
                                                      DiffFormFreeModule
        from sage.manifolds.differentiable.multivector_module import \
                                                   MultivectorFreeModule
        from sage.manifolds.differentiable.automorphismfield_group \
                                      import AutomorphismFieldParalGroup
        if isinstance(other, (TensorFieldModule, TensorFieldFreeModule)):
            # coercion by domain restriction
            return (self._tensor_type == other._tensor_type
                    and self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(other._ambient_domain))
        if isinstance(other, DiffFormFreeModule):
            # coercion of p-forms to type-(0,p) tensor fields
            return (self._fmodule is other.base_module()
                    and self._tensor_type == (0, other.degree()))
        if isinstance(other, MultivectorFreeModule):
            # coercion of p-vector fields to type-(p,0) tensor fields
            return (self._fmodule is other.base_module()
                    and self._tensor_type == (other.degree(),0))
        if isinstance(other, AutomorphismFieldParalGroup):
            # coercion of automorphism fields to type-(1,1) tensor fields
            return (self._fmodule is other.base_module()
                    and self._tensor_type == (1,1))
        return False

    #### End of parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: T12 = M.tensor_field_module((1,2))
            sage: T12._repr_()
            'Free module T^(1,2)(M) of type-(1,2) tensors fields on the 2-dimensional differentiable manifold M'
            sage: repr(T12)  # indirect doctest
            'Free module T^(1,2)(M) of type-(1,2) tensors fields on the 2-dimensional differentiable manifold M'
            sage: T12  # indirect doctest
            Free module T^(1,2)(M) of type-(1,2) tensors fields on the
             2-dimensional differentiable manifold M
        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-({},{})".format(self._tensor_type[0],
                                                self._tensor_type[1])
        description += " tensors fields "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {}".format(self._domain) + \
                           " mapped into the {}".format(self._ambient_domain)
        return description
