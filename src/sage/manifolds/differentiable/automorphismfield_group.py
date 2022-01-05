r"""
Group of Tangent-Space Automorphism Fields

Given a differentiable manifold `U` and a differentiable map
`\Phi: U \rightarrow M` to a differentiable manifold `M` (possibly `U = M`
and `\Phi=\mathrm{Id}_M`), the *group of tangent-space automorphism fields*
associated with `U` and `\Phi` is the general linear group
`\mathrm{GL}(\mathfrak{X}(U,\Phi))` of the module `\mathfrak{X}(U,\Phi)` of
vector fields along `U` with values on `M\supset \Phi(U)` (see
:class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`).
Note that `\mathfrak{X}(U, \Phi)` is a module over
`C^k(U)`, the algebra of differentiable scalar fields on `U`.
Elements of `\mathrm{GL}(\mathfrak{X}(U, \Phi))` are fields along `U`
of automorphisms of tangent spaces to `M`.

Two classes implement `\mathrm{GL}(\mathfrak{X}(U, \Phi))` depending
whether `M` is parallelizable or not:
:class:`AutomorphismFieldParalGroup` and :class:`AutomorphismFieldGroup`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks
- Michael Jung (2019): improve treatment of the identity element

REFERENCES:

- Chap. 15 of [God1968]_

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.groups import Groups
from sage.misc.cachefunc import cached_method
from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
from sage.manifolds.differentiable.vectorfield_module import (VectorFieldModule,
                                                              VectorFieldFreeModule)
from sage.manifolds.differentiable.automorphismfield import (AutomorphismField,
                                                             AutomorphismFieldParal)

class AutomorphismFieldGroup(UniqueRepresentation, Parent):
    r"""
    General linear group of the module of vector fields along a differentiable
    manifold `U` with values on a differentiable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M` (possibly `U = M`
    and `\Phi = \mathrm{Id}_M`), the *group of tangent-space automorphism
    fields* associated with `U` and `\Phi` is the general linear group
    `\mathrm{GL}(\mathfrak{X}(U,\Phi))` of the module `\mathfrak{X}(U,\Phi)` of
    vector fields along `U` with values on `M \supset \Phi(U)` (see
    :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`).
    Note that `\mathfrak{X}(U,\Phi)` is a module over
    `C^k(U)`, the algebra of differentiable scalar fields on `U`.
    Elements of `\mathrm{GL}(\mathfrak{X}(U,\Phi))` are fields along `U` of
    automorphisms of tangent spaces to `M`.

    .. NOTE::

        If `M` is parallelizable, then :class:`AutomorphismFieldParalGroup`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` --
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`;
      module `\mathfrak{X}(U,\Phi)` of vector fields along `U` with values on `M`

    EXAMPLES:

    Group of tangent-space automorphism fields of the 2-sphere::

        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:             intersection_name='W',
        ....:             restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: G = M.automorphism_field_group() ; G
        General linear group of the Module X(M) of vector fields on the
         2-dimensional differentiable manifold M

    ``G`` is the general linear group of the vector field module
    `\mathfrak{X}(M)`::

        sage: XM = M.vector_field_module() ; XM
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: G is XM.general_linear_group()
        True

    ``G`` is a non-abelian group::

        sage: G.category()
        Category of groups
        sage: G in Groups()
        True
        sage: G in CommutativeAdditiveGroups()
        False

    The elements of ``G`` are tangent-space automorphisms::

        sage: a = G.an_element(); a
        Field of tangent-space automorphisms on the 2-dimensional
         differentiable manifold M
        sage: a.parent() is G
        True
        sage: a.restrict(U).display()
        2 ∂/∂x⊗dx + 2 ∂/∂y⊗dy
        sage: a.restrict(V).display()
        2 ∂/∂u⊗du + 2 ∂/∂v⊗dv

    The identity element of the group ``G``::

        sage: e = G.one() ; e
        Field of tangent-space identity maps on the 2-dimensional
         differentiable manifold M
        sage: eU = U.default_frame() ; eU
        Coordinate frame (U, (∂/∂x,∂/∂y))
        sage: eV = V.default_frame() ; eV
        Coordinate frame (V, (∂/∂u,∂/∂v))
        sage: e.display(eU)
        Id = ∂/∂x⊗dx + ∂/∂y⊗dy
        sage: e.display(eV)
        Id = ∂/∂u⊗du + ∂/∂v⊗dv

    """

    Element = AutomorphismField

    def __init__(self, vector_field_module):
        r"""
        See :class:`AutomorphismfieldGroup` for documentation and examples.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:  intersection_name='W', restrictions1= x>0,
            ....:  restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: from sage.manifolds.differentiable.automorphismfield_group \
            ....:                                 import AutomorphismFieldGroup
            sage: G = AutomorphismFieldGroup(M.vector_field_module()) ; G
            General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M
            sage: TestSuite(G).run(skip='_test_elements')

        ``_test_elements`` does not pass due to the failure
        of ``_test_pickling`` in
        :class:`sage.manifolds.differentiable.tensorfield.TensorField`.

        """
        if not isinstance(vector_field_module, VectorFieldModule):
            raise TypeError("{} is not a module of vector fields".format(
                            vector_field_module))
        Parent.__init__(self, category=Groups())
        self._vmodule = vector_field_module


    #### Parent methods ####

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a field of tangent-space automorphisms.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:  intersection_name='W', restrictions1= x>0,
            ....:  restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: G = M.automorphism_field_group()
            sage: a = G(1); a
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M
            sage: a = G([[1+x^2, 0], [0, 1+y^2]], frame=c_xy.frame(), name='a'); a
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M
            sage: a.display(c_xy.frame())
            a = (x^2 + 1) ∂/∂x⊗dx + (y^2 + 1) ∂/∂y⊗dy

        """
        if hasattr(comp, 'is_trivial_zero'):
            if (comp - 1).is_trivial_zero():
                return self.one()
        elif comp == 1:
            return self.one()
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._vmodule, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some specific field of tangent-space automorphisms.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:  intersection_name='W', restrictions1= x>0,
            ....:  restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: G = M.automorphism_field_group()
            sage: a = G.an_element() ; a
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: a.restrict(U).display()
            2 ∂/∂x⊗dx + 2 ∂/∂y⊗dy
            sage: a.restrict(V).display()
            2 ∂/∂u⊗du + 2 ∂/∂v⊗dv
            sage: a == G.an_element()  # indirect doctest
            True

        """
        resu = self.element_class(self._vmodule)
        for dom in resu.domain().subsets():
            if dom.is_manifestly_parallelizable():
                fmodule = dom.vector_field_module()
                idm = fmodule.identity_map()
                rst = fmodule.automorphism()
                for frame, comp in idm._components.items():
                    rst._components[frame] = 2 * comp
                resu._restrictions[dom] = rst
        return resu

    #### End of parent methods ####


    #### Monoid methods ####

    @cached_method
    def one(self):
        r"""
        Return identity element of ``self``.

        The group identity element is the field of tangent-space identity maps.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
          representing the identity element

        EXAMPLES:

        Identity element of the group of tangent-space automorphism fields of
        the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: G = M.automorphism_field_group()
            sage: G.one()
            Field of tangent-space identity maps on the 2-dimensional differentiable manifold M
            sage: G.one().restrict(U)[:]
            [1 0]
            [0 1]
            sage: G.one().restrict(V)[:]
            [1 0]
            [0 1]

        """
        # Specific initializations for the field of identity maps:
        resu = self._element_constructor_(name='Id', latex_name=r'\mathrm{Id}')
        resu._inverse = resu
        for dom in resu._domain._subsets:
            if dom.is_manifestly_parallelizable():
                fmodule = dom.vector_field_module()
                resu._restrictions[dom] = fmodule.identity_map(name='Id',
                                                     latex_name=r'\mathrm{Id}')
        resu._is_identity = True
        resu.set_immutable()
        return resu

    #### End of monoid methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: G = M.automorphism_field_group()
            sage: G._repr_()
            'General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M'
            sage: repr(G)  # indirect doctest
            'General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M'
            sage: G  # indirect doctest
            General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M

        """
        return "General linear group of the {}".format(self._vmodule)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: G = M.automorphism_field_group()
            sage: G._latex_()
            \mathrm{GL}\left( \mathfrak{X}\left(M\right) \right)
            sage: latex(G)  # indirect doctest
            \mathrm{GL}\left( \mathfrak{X}\left(M\right) \right)

        """
        from sage.misc.latex import latex
        return r"\mathrm{GL}\left("+ latex(self._vmodule)+ r"\right)"

    def base_module(self):
        r"""
        Return the vector-field module of which ``self`` is the general
        linear group.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`

        EXAMPLES:

        Base module of the group of tangent-space automorphism fields of
        the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: G = M.automorphism_field_group()
            sage: G.base_module()
            Module X(M) of vector fields on the 2-dimensional differentiable
             manifold M
            sage: G.base_module() is M.vector_field_module()
            True

        """
        return self._vmodule


#******************************************************************************

class AutomorphismFieldParalGroup(FreeModuleLinearGroup):
    r"""
    General linear group of the module of vector fields along a differentiable
    manifold `U` with values on a parallelizable manifold `M`.

    Given a differentiable manifold `U` and a differentiable map
    `\Phi: U \rightarrow M` to a parallelizable  manifold `M` (possibly `U = M`
    and `\Phi = \mathrm{Id}_M`), the *group of tangent-space automorphism
    fields* associated with `U` and `\Phi` is the general linear group
    `\mathrm{GL}(\mathfrak{X}(U, \Phi))` of the module `\mathfrak{X}(U, \Phi)`
    of vector fields along `U` with values on `M \supset \Phi(U)` (see
    :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`).
    Note that `\mathfrak{X}(U, \Phi)` is a free module over `C^k(U)`,
    the algebra of differentiable scalar fields on `U`.
    Elements of `\mathrm{GL}(\mathfrak{X}(U, \Phi))` are fields along `U` of
    automorphisms of tangent spaces to `M`.

    .. NOTE::

        If `M` is not parallelizable, the class
        :class:`AutomorphismFieldGroup` must be used instead.

    INPUT:

    - ``vector_field_module`` --
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`;
      free module `\mathfrak{X}(U,\Phi)` of vector fields along `U`
      with values on `M`

    EXAMPLES:

    Group of tangent-space automorphism fields of a 2-dimensional
    parallelizable manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: XM = M.vector_field_module() ; XM
        Free module X(M) of vector fields on the 2-dimensional differentiable
         manifold M
        sage: G = M.automorphism_field_group(); G
        General linear group of the Free module X(M) of vector fields on the
         2-dimensional differentiable manifold M
        sage: latex(G)
        \mathrm{GL}\left( \mathfrak{X}\left(M\right) \right)

    ``G`` is nothing but the general linear group of the module
    `\mathfrak{X}(M)`::

        sage: G is XM.general_linear_group()
        True

    ``G`` is a group::

        sage: G.category()
        Category of groups
        sage: G in Groups()
        True

    It is not an abelian group::

        sage: G in CommutativeAdditiveGroups()
        False

    The elements of ``G`` are tangent-space automorphisms::

        sage: G.Element
        <class 'sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal'>
        sage: a = G.an_element() ; a
        Field of tangent-space automorphisms on the 2-dimensional
         differentiable manifold M
        sage: a.parent() is G
        True

    As automorphisms of `\mathfrak{X}(M)`, the elements of ``G`` map a vector
    field to a vector field::

        sage: v = XM.an_element() ; v
        Vector field on the 2-dimensional differentiable manifold M
        sage: v.display()
        2 ∂/∂x + 2 ∂/∂y
        sage: a(v)
        Vector field on the 2-dimensional differentiable manifold M
        sage: a(v).display()
        2 ∂/∂x - 2 ∂/∂y

    Indeed the matrix of ``a`` with respect to the frame
    `(\partial_x, \partial_y)` is::

        sage: a[X.frame(),:]
        [ 1  0]
        [ 0 -1]

    The elements of ``G`` can also be considered as tensor fields of
    type `(1,1)`::

        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2
        sage: a.domain()
        2-dimensional differentiable manifold M
        sage: a.display()
        ∂/∂x⊗dx - ∂/∂y⊗dy

    The identity element of the group ``G`` is::

        sage: id = G.one() ; id
        Field of tangent-space identity maps on the 2-dimensional
         differentiable manifold M
        sage: id*a == a
        True
        sage: a*id == a
        True
        sage: a*a^(-1) == id
        True
        sage: a^(-1)*a == id
        True

    Construction of an element by providing its components with
    respect to the manifold's default frame (frame associated to
    the coordinates `(x,y)`)::

        sage: b = G([[1+x^2,0], [0,1+y^2]]) ; b
        Field of tangent-space automorphisms on the 2-dimensional
         differentiable manifold M
        sage: b.display()
        (x^2 + 1) ∂/∂x⊗dx + (y^2 + 1) ∂/∂y⊗dy
        sage: (~b).display()  # the inverse automorphism
        1/(x^2 + 1) ∂/∂x⊗dx + 1/(y^2 + 1) ∂/∂y⊗dy

    We check the group law on these elements::

        sage: (a*b)^(-1) == b^(-1) * a^(-1)
        True

    Invertible tensor fields of type `(1,1)` can be converted to
    elements of ``G``::

        sage: t = M.tensor_field(1, 1, name='t')
        sage: t[:] = [[1+exp(y), x*y], [0, 1+x^2]]
        sage: t1 = G(t) ; t1
        Field of tangent-space automorphisms t on the 2-dimensional
         differentiable manifold M
        sage: t1 in G
        True
        sage: t1.display()
        t = (e^y + 1) ∂/∂x⊗dx + x*y ∂/∂x⊗dy + (x^2 + 1) ∂/∂y⊗dy
        sage: t1^(-1)
        Field of tangent-space automorphisms t^(-1) on the 2-dimensional
         differentiable manifold M
        sage: (t1^(-1)).display()
        t^(-1) = 1/(e^y + 1) ∂/∂x⊗dx - x*y/(x^2 + (x^2 + 1)*e^y + 1) ∂/∂x⊗dy
         + 1/(x^2 + 1) ∂/∂y⊗dy

    Since any automorphism field can be considered as a tensor field of
    type-`(1,1)` on ``M``, there is a coercion map from ``G`` to the
    module `T^{(1,1)}(M)` of type-`(1,1)` tensor fields::

        sage: T11 = M.tensor_field_module((1,1)) ; T11
        Free module T^(1,1)(M) of type-(1,1) tensors fields on the
         2-dimensional differentiable manifold M
        sage: T11.has_coerce_map_from(G)
        True

    An explicit call of this coercion map is::

        sage: tt = T11(t1) ; tt
        Tensor field t of type (1,1) on the 2-dimensional differentiable
         manifold M
        sage: tt == t
        True

    An implicit call of the coercion map is performed to subtract an
    element of ``G`` from an element of `T^{(1,1)}(M)`::

        sage: s = t - t1 ; s
        Tensor field t-t of type (1,1) on
         the 2-dimensional differentiable manifold M
        sage: s.parent() is T11
        True
        sage: s.display()
        t-t = 0

    as well as for the reverse operation::

        sage: s = t1 - t ; s
        Tensor field t-t of type (1,1) on the 2-dimensional differentiable
         manifold M
        sage: s.display()
        t-t = 0

    TESTS::

        sage: TestSuite(G).run()

    """

    Element = AutomorphismFieldParal

    def __init__(self, vector_field_module):
        r"""
        See :class:`AutomorphismfieldParalGroup` for documentation and
        examples.

        TESTS::

            sage: M = Manifold(2, 'M') ; M
            2-dimensional differentiable manifold M
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: G = M.automorphism_field_group(); G
            General linear group of the Free module X(M) of vector fields on
             the 2-dimensional differentiable manifold M
            sage: TestSuite(G).run()

        """
        if not isinstance(vector_field_module, VectorFieldFreeModule):
            raise TypeError("{} is not a free module of vector fields".format(
                            vector_field_module))
        FreeModuleLinearGroup.__init__(self, vector_field_module)

