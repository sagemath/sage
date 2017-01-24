r"""
Vector Field Modules

The set of vector fields along a differentiable manifold `U` with values on
a differentiable manifold `M` via a differentiable map `\Phi: U \to M`
(possibly `U = M` and `\Phi=\mathrm{Id}_M`) is a module over the algebra
`C^k(U)` of differentiable scalar fields on `U`. If `\Phi` is the identity
map, this module is considered a Lie algebroid under the Lie bracket `[\ ,\ ]`
(cf. :wikipedia:`Lie_algebroid`). It is a free module if and only if `M` is
parallelizable. Accordingly, there are two classes for vector field modules:

- :class:`VectorFieldModule` for vector fields with values on a
  generic (in practice, not parallelizable) differentiable manifold `M`.
- :class:`VectorFieldFreeModule` for vector fields with values on a
  parallelizable manifold `M`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): structure of Lie algebroid (:trac:`20771`)

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.misc.cachefunc import cached_method
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.differentiable.vectorfield import (VectorField, VectorFieldParal)

class VectorFieldModule(UniqueRepresentation, Parent):
    r"""
    Module of vector fields along a differentiable manifold `U`
    with values on a differentiable manifold `M`, via a differentiable
    map `U \rightarrow M`.

    Given a differentiable map

    .. MATH::

        \Phi:\  U \longrightarrow M,

    the *vector field module* `\mathcal{X}(U,\Phi)` is the set of
    all vector fields of the type

    .. MATH::

        v:\ U  \longrightarrow TM

    (where `TM` is the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M,

    where `T_{\Phi(p)}M` is the tangent space to `M` at the point `\Phi(p)`.

    The set `\mathcal{X}(U,\Phi)` is a module over `C^k(U)`, the ring
    (algebra) of differentiable scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    Furthermore, it is a Lie algebroid under the Lie bracket (cf.
    :wikipedia:`Lie_algebroid`)

    .. MATH::

        [X, Y] = X \circ Y - Y \circ X

    over the scalarfields if `\Phi` is the identity map. That is to say
    the Lie bracket is antisymmetric, bilinear over the base field,
    satisfies the Jacobi identity, and `[X, fY] = X(f) Y + f[X, Y]`.

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`; we then denote
    `\mathcal{X}(M,\mathrm{Id}_M)` by merely `\mathcal{X}(M)`. Other common
    cases are `\Phi` being an immersion and `\Phi` being a curve in `M`
    (`U` is then an open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, the class :class:`VectorFieldFreeModule`
        should be used instead.

    INPUT:

    - ``domain`` -- differentiable manifold `U` along which the
      vector fields are defined
    - ``dest_map`` -- (default: ``None``) destination map
      `\Phi:\ U \rightarrow M`
      (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`);
      if ``None``, it is assumed that `U = M` and `\Phi` is the identity
      map of `M` (case of vector fields *on* `M`)

    EXAMPLES:

    Module of vector fields on the 2-sphere::

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
        sage: XM = M.vector_field_module() ; XM
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M

    `\mathcal{X}(M)` is a module over the algebra `C^k(M)`::

        sage: XM.category()
        Category of modules over Algebra of differentiable scalar fields on the
         2-dimensional differentiable manifold M
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    `\mathcal{X}(M)` is not a free module::

        sage: isinstance(XM, FiniteRankFreeModule)
        False

    because `M = S^2` is not parallelizable::

        sage: M.is_manifestly_parallelizable()
        False

    On the contrary, the module of vector fields on `U` is a free module,
    since `U` is parallelizable (being a coordinate domain)::

        sage: XU = U.vector_field_module()
        sage: isinstance(XU, FiniteRankFreeModule)
        True
        sage: U.is_manifestly_parallelizable()
        True

    The zero element of the module::

        sage: z = XM.zero() ; z
        Vector field zero on the 2-dimensional differentiable manifold M
        sage: z.display(c_xy.frame())
        zero = 0
        sage: z.display(c_uv.frame())
        zero = 0

    The module `\mathcal{X}(M)` coerces to any module of vector fields defined
    on a subdomain of `M`, for instance `\mathcal{X}(U)`::

        sage: XU.has_coerce_map_from(XM)
        True
        sage: XU.coerce_map_from(XM)
        Conversion map:
          From: Module X(M) of vector fields on the 2-dimensional differentiable manifold M
          To:   Free module X(U) of vector fields on the Open subset U of the 2-dimensional differentiable manifold M

    The conversion map is actually the restriction of vector fields defined
    on `M` to `U`.

    """
    Element = VectorField

    def __init__(self, domain, dest_map=None):
        r"""
        Construct the module of vector fields taking values on a (a priori)
        non-parallelizable differentiable manifold.

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
            sage: from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
            sage: XM = VectorFieldModule(M, dest_map=M.identity_map()); XM
            Module X(M) of vector fields on the 2-dimensional differentiable
             manifold M
            sage: XM is M.vector_field_module()
            True
            sage: TestSuite(XM).run(skip='_test_elements')

        In the above test suite, _test_elements is skipped because of the
        _test_pickling error of the elements (to be fixed in class
        TensorField)

        """
        self._domain = domain
        name = "X(" + domain._name
        latex_name = r"\mathcal{X}\left(" + domain._latex_name
        if dest_map is None:
            dest_map = domain.identity_map()
        self._dest_map = dest_map
        if dest_map is domain.identity_map():
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + self._dest_map._name + ")"
            latex_name += "," + self._dest_map._latex_name + r"\right)"
        self._ambient_domain = self._dest_map._codomain
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        # Dictionary of the tensor modules built on self
        #   (dict. keys = (k,l) --the tensor type)
        self._tensor_modules = {(1,0): self} # self is considered as the set of
                                             # tensors of type (1,0)
        # Dictionary of exterior powers of the dual of self
        #   (keys = p --the power degree) :
        self._dual_exterior_powers = {}

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of the module

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: XM = M.vector_field_module()
            sage: v = XM([-x,y], frame=c_xy.frame(), name='v'); v
            Vector field v on the 2-dimensional differentiable manifold M
            sage: v.display()
            v = -x d/dx + y d/dy
            sage: XM(0) is XM.zero()
            True

        """
        if comp == 0:
            return self.zero()
        if isinstance(comp, VectorField):
            if (self._domain.is_subset(comp._domain)
                   and self._ambient_domain.is_subset(comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {} ".format(comp) +
                                 "to a vector field in {}".format(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) element of the module.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: XM = M.vector_field_module()
            sage: XM._an_element_()
            Vector field on the 2-dimensional differentiable manifold M

        """
        resu = self.element_class(self)
        # Non-trivial open covers of the domain:
        open_covers = self._domain.open_covers()[1:]  # the open cover 0
                                                      # is trivial
        if open_covers != []:
            oc = open_covers[0]  # the first non-trivial open cover is selected
            for dom in oc:
                vmodule_dom = dom.vector_field_module(
                                         dest_map=self._dest_map.restrict(dom))
                resu.set_restriction(vmodule_dom._an_element_())
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: XM = M.vector_field_module()
            sage: XU = U.vector_field_module()
            sage: XM._coerce_map_from_(XU)
            False
            sage: XU._coerce_map_from_(XM)
            True

        """
        if isinstance(other, (VectorFieldModule, VectorFieldFreeModule)):
            return self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM._repr_()
            'Module X(M) of vector fields on the 2-dimensional differentiable manifold M'
            sage: repr(XM)  # indirect doctest
            'Module X(M) of vector fields on the 2-dimensional differentiable manifold M'
            sage: XM  # indirect doctest
            Module X(M) of vector fields on the 2-dimensional differentiable
             manifold M

        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of vector fields "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += ("along the {}".format(self._domain)
                            + " mapped into the {}".format(self._ambient_domain))
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM._latex_()
            '\\mathcal{X}\\left(M\\right)'
            sage: latex(XM)  # indirect doctest
            \mathcal{X}\left(M\right)

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def domain(self):
        r"""
        Return the domain of the vector fields in this module.

        If the module is `\mathcal{X}(U,\Phi)`, returns the domain `U` of
        `\Phi`.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the domain of the vector fields that belong to this
          module

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.domain()
            5-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Phi = U.diff_map(M, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.domain()
            2-dimensional differentiable manifold U

        """
        return self._domain

    def ambient_domain(self):
        r"""
        Return the manifold in which the vector fields of this module take
        their values.

        If the module is `\mathcal{X}(U,\Phi)`, returns the codomain `M` of
        `\Phi`.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the manifold in which the vector fields of this
          module take their values

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.ambient_domain()
            5-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Phi = U.diff_map(M, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.ambient_domain()
            5-dimensional differentiable manifold M

        """
        return self._ambient_domain

    def destination_map(self):
        r"""
        Return the differential map associated to this module.

        The differential map associated to this module is the map

        .. MATH::

            \Phi:\  U \longrightarrow M

        such that this module is the set `\mathcal{X}(U,\Phi)` of all
        vector fields of the type

        .. MATH::

            v:\ U  \longrightarrow TM

        (where `TM` is the tangent bundle of `M`) such that

        .. MATH::

            \forall p \in U,\ v(p) \in T_{\Phi(p)}M,

        where `T_{\Phi(p)}M` is the tangent space to `M` at the
        point `\Phi(p)`.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the differential map `\Phi`

        EXAMPLES::

            sage: M = Manifold(5, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.destination_map()
            Identity map Id_M of the 5-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Phi = U.diff_map(M, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.destination_map()
            Differentiable map Phi from the 2-dimensional differentiable
             manifold U to the 5-dimensional differentiable manifold M

        """
        return self._dest_map

    def tensor_module(self, k, l):
        r"""
        Return the module of type-`(k,l)` tensors on ``self``.

        INPUT:

        - ``k`` -- non-negative integer; the contravariant rank,
          the tensor type being `(k,l)`
        - ``l`` -- non-negative integer; the covariant rank,
          the tensor type being `(k,l)`

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`
          representing the module `T^{(k,l)}(U,\Phi)` of type-`(k,l)`
          tensors on the vector field module

        EXAMPLES:

        A tensor field module on a 2-dimensional differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.tensor_module(1,2)
            Module T^(1,2)(M) of type-(1,2) tensors fields on the 2-dimensional
             differentiable manifold M

        The special case of tensor fields of type (1,0)::

            sage: XM.tensor_module(1,0)
            Module X(M) of vector fields on the 2-dimensional differentiable
             manifold M

        The result is cached::

            sage: XM.tensor_module(1,2) is XM.tensor_module(1,2)
            True
            sage: XM.tensor_module(1,0) is XM
            True

        See
        :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`
        for more examples and documentation.

        """
        from sage.manifolds.differentiable.tensorfield_module import TensorFieldModule
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def dual_exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of the dual of the vector field
        module.

        If the vector field module is `\mathcal{X}(U,\Phi)`, the
        `p`-th exterior power of its dual is the set `\Lambda^p(U, \Phi)`
        of `p`-forms along `U` with values on `\Phi(U)`. It is a module
        over `C^k(U)`, the ring (algebra) of differentiable scalar
        fields on `U`.

        INPUT:

        - ``p`` -- non-negative integer

        OUTPUT:

        - for `p \geq 1`, instance of
          :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormModule`
          representing the module `\Lambda^p(U,\Phi)`; for `p=0`, the
          base ring, i.e. `C^k(U)`, is returned instead

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.dual_exterior_power(2)
            Module /\^2(M) of 2-forms on the 2-dimensional differentiable
             manifold M
            sage: XM.dual_exterior_power(1)
            Module /\^1(M) of 1-forms on the 2-dimensional differentiable
             manifold M
            sage: XM.dual_exterior_power(1) is XM.dual()
            True
            sage: XM.dual_exterior_power(0)
            Algebra of differentiable scalar fields on the 2-dimensional
             differentiable manifold M
            sage: XM.dual_exterior_power(0) is M.scalar_field_algebra()
            True

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormModule`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.diff_form_module import DiffFormModule
        if p == 0:
            return self._ring
        if p not in self._dual_exterior_powers:
            self._dual_exterior_powers[p] = DiffFormModule(self, p)
        return self._dual_exterior_powers[p]

    def dual(self):
        r"""
        Return the dual module.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.dual()
            Module /\^1(M) of 1-forms on the 2-dimensional differentiable
             manifold M

        """
        return self.dual_exterior_power(1)

    def general_linear_group(self):
        r"""
        Return the general linear group of ``self``.

        If the vector field module is `\mathcal{X}(U,\Phi)`, the *general
        linear group* is the group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
        automorphisms of `\mathcal{X}(U, \Phi)`. Note that an automorphism
        of `\mathcal{X}(U,\Phi)` can also be viewed as a *field* along `U`
        of automorphisms of the tangent spaces of `M \supset \Phi(U)`.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldGroup`
          representing `\mathrm{GL}(\mathcal{X}(U,\Phi))`

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.general_linear_group()
            General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldGroup`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.automorphismfield_group import AutomorphismFieldGroup
        return AutomorphismFieldGroup(self)

    def tensor(self, tensor_type, name=None, latex_name=None, sym=None,
               antisym=None, specific_type=None):
        r"""
        Construct a tensor on ``self``.

        The tensor is actually a tensor field on the domain of
        the vector field module.

        INPUT:

        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank
          and l the covariant rank
        - ``name`` -- (string; default: ``None``) name given to the tensor
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the tensor; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention position=0 for the first argument; for instance:

          * ``sym=(0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym=[(0,2),(1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments

        - ``antisym`` -- (default: ``None``) antisymmetry or list of
          antisymmetries among the arguments, with the same convention
          as for ``sym``
        - ``specific_type`` -- (default: ``None``) specific subclass of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField` for
          the output

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing the tensor defined on the vector field module with the
          provided characteristics

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.tensor((1,2), name='t')
            Tensor field t of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: XM.tensor((1,0), name='a')
            Vector field a on the 2-dimensional differentiable manifold M
            sage: XM.tensor((0,2), name='a', antisym=(0,1))
            2-form a on the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.automorphismfield import \
                                                              AutomorphismField
        from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
        if tensor_type==(1,0):
            return self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            return self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type == (1,1) and specific_type is not None:
            if issubclass(specific_type, AutomorphismField):
                return self.automorphism(name=name, latex_name=latex_name)
        elif (tensor_type[0] == 0 and tensor_type[1] > 1
              and antisym is not None and antisym != []):
            if isinstance(antisym, list):
                antisym0 = antisym[0]
            else:
                antisym0 = antisym
            if len(antisym0) == tensor_type[1]:
                return self.alternating_form(tensor_type[1], name=name,
                                             latex_name=latex_name)
            else:
                return self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name,
                                 sym=sym, antisym=antisym)
        elif tensor_type==(0,2) and specific_type is not None:
            if issubclass(specific_type, PseudoRiemannianMetric):
                return self.metric(name, latex_name=latex_name)
                # NB: the signature is not treated
        # Generic case
        return self.tensor_module(*tensor_type).element_class(self,
                        tensor_type, name=name, latex_name=latex_name, sym=sym,
                        antisym=antisym)

    def alternating_form(self, degree, name=None, latex_name=None):
        r"""
        Construct an alternating form on the vector field module.

        An alternating form on the vector field module is actually a
        differential form along the differentiable manifold `U` over which
        the vector field module is defined.

        INPUT:

        - ``degree`` -- the degree of the alternating form
          (i.e. its tensor rank)
        - ``name`` -- (string; optional) name given to the alternating form
        - ``latex_name`` -- (string; optional) LaTeX symbol to denote
          the alternating form; if none is provided, the
          LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.alternating_form(2, name='a')
            2-form a on the 2-dimensional differentiable manifold M
            sage: XM.alternating_form(1, name='a')
            1-form a on the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
            for more examples and documentation.

        """
        return self.dual_exterior_power(degree).element_class(self, degree,
                                              name=name, latex_name=latex_name)

    def linear_form(self, name=None, latex_name=None):
        r"""
        Construct a linear form on the vector field module.

        A linear form on the vector field module is actually a field
        of linear forms (i.e. a 1-form) along the differentiable manifold `U`
        over which the vector field module is defined.

        INPUT:

        - ``name`` -- (string; optional) name given to the linear form
        - ``latex_name`` -- (string; optional) LaTeX symbol to denote
          the linear form; if none is provided, the LaTeX symbol is
          set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.linear_form()
            1-form on the 2-dimensional differentiable manifold M
            sage: XM.linear_form(name='a')
            1-form a on the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
            for more examples and documentation.

        """
        return self.dual_exterior_power(1).element_class(self, 1, name=name,
                                                         latex_name=latex_name)

    def automorphism(self, name=None, latex_name=None):
        r"""
        Construct an automorphism of the vector field module.

        An automorphism of the vector field module is actually a field
        of tangent-space automorphisms along the differentiable manifold `U`
        over which the vector field module is defined.

        INPUT:

        - ``name`` -- (string; optional) name given to the automorphism
        - ``latex_name`` -- (string; optional) LaTeX symbol to denote
          the automorphism; if none is provided, the LaTeX symbol is
          set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.automorphism()
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: XM.automorphism(name='a')
            Field of tangent-space automorphisms a on the 2-dimensional
             differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
            for more examples and documentation.

        """
        return self.general_linear_group().element_class(self, name=name,
                                                         latex_name=latex_name)

    @cached_method
    def identity_map(self, name='Id', latex_name=None):
        r"""
        Construct the identity map on the vector field module.

        The identity map on the vector field module is actually a field
        of tangent-space identity maps along the differentiable manifold `U`
        over which the vector field module is defined.

        INPUT:

        - ``name`` -- (string; default: ``'Id'``) name given to the
          identity map
        - ``latex_name`` -- (string; optional) LaTeX symbol to denote
          the identity map;  if none is provided, the LaTeX symbol is
          set to ``'\mathrm{Id}'`` if ``name`` is ``'Id'`` and
          to ``name`` otherwise

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.identity_map()
            Field of tangent-space identity maps on the 2-dimensional
             differentiable manifold M

        """
        resu = self.general_linear_group().one()
        if latex_name is None:
            latex_name = name
        resu.set_name(name=name, latex_name=latex_name)
        return resu

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.zero()
            Vector field zero on the 2-dimensional differentiable manifold M
        """
        elt = self.element_class(self, name='zero', latex_name='0')
        for frame in self._domain._frames:
            if self._dest_map.restrict(frame._domain) == frame._dest_map:
                elt.add_comp(frame)
                # (since new components are initialized to zero)
        return elt

    def metric(self, name, signature=None, latex_name=None):
        r"""
        Construct a pseudo-Riemannian metric (nondegenerate symmetric bilinear
        form) on the current vector field module.

        A pseudo-Riemannian metric of the vector field module is actually a
        field of tangent-space non-degenerate symmetric bilinear forms along
        the manifold `U` on which the vector field module is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (integer; default: ``None``) signature `S` of the
          metric: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal
          writing of the metric components; if ``signature`` is not provided,
          `S` is set to the manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the metric; if ``None``, it is formed from ``name``

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined pseudo-Riemannian metric.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: XM.metric('g')
            Riemannian metric g on the 2-dimensional differentiable manifold M
            sage: XM.metric('g', signature=0)
            Lorentzian metric g on the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
            for more documentation.

        """
        from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
        return PseudoRiemannianMetric(self, name, signature=signature,
                                      latex_name=latex_name)


#******************************************************************************

class VectorFieldFreeModule(FiniteRankFreeModule):
    r"""
    Free module of vector fields along a differentiable manifold `U` with
    values on a parallelizable manifold `M`, via a differentiable map
    `U \rightarrow M`.

    Given a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M

    the *vector field module* `\mathcal{X}(U,\Phi)` is the set of all vector
    fields of the type

    .. MATH::

        v:\ U  \longrightarrow TM

    (where `TM` is the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)} M,

    where `T_{\Phi(p)} M` is the tangent space to `M` at the point `\Phi(p)`.

    Since `M` is parallelizable, the set `\mathcal{X}(U,\Phi)` is a
    free module over `C^k(U)`, the ring (algebra) of differentiable
    scalar fields on `U` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`).
    In fact, it carries the structure of a finite-dimensional Lie algebroid
    (cf. :wikipedia:`Lie_algebroid`).

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U=M` and `\Phi = \mathrm{Id}_M`; we then denote
    `\mathcal{X}(M,\mathrm{Id}_M)` by merely `\mathcal{X}(M)`. Other common
    cases are `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is
    then an open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, the class :class:`VectorFieldModule`
        should be used instead, for `\mathcal{X}(U,\Phi)` is no longer a
        free module.

    INPUT:

    - ``domain`` -- differentiable manifold `U` along which the vector fields
      are defined
    - ``dest_map`` -- (default: ``None``) destination map
      `\Phi:\ U \rightarrow M`
      (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`); if
      ``None``, it is assumed that `U=M` and `\Phi` is the identity map of
      `M` (case of vector fields *on* `M`)

    EXAMPLES:

    Module of vector fields on `\RR^2`::

        sage: M = Manifold(2, 'R^2')
        sage: cart.<x,y> = M.chart()  # Cartesian coordinates on R^2
        sage: XM = M.vector_field_module() ; XM
        Free module X(R^2) of vector fields on the 2-dimensional differentiable
         manifold R^2
        sage: XM.category()
        Category of finite dimensional modules
         over Algebra of differentiable scalar fields
         on the 2-dimensional differentiable manifold R^2
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    Since `\RR^2` is obviously parallelizable, ``XM`` is a free module::

        sage: isinstance(XM, FiniteRankFreeModule)
        True

    Some elements::

        sage: XM.an_element().display()
        2 d/dx + 2 d/dy
        sage: XM.zero().display()
        zero = 0
        sage: v = XM([-y,x]) ; v
        Vector field on the 2-dimensional differentiable manifold R^2
        sage: v.display()
        -y d/dx + x d/dy

    An example of module of vector fields with a destination map `\Phi`
    different from the identity map, namely a mapping
    `\Phi: I \rightarrow \RR^2`, where `I` is an open interval of `\RR`::

        sage: I = Manifold(1, 'I')
        sage: canon.<t> = I.chart('t:(0,2*pi)')
        sage: Phi = I.diff_map(M, coord_functions=[cos(t), sin(t)], name='Phi',
        ....:                      latex_name=r'\Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold
         I to the 2-dimensional differentiable manifold R^2
        sage: Phi.display()
        Phi: I --> R^2
           t |--> (x, y) = (cos(t), sin(t))
        sage: XIM = I.vector_field_module(dest_map=Phi) ; XIM
        Free module X(I,Phi) of vector fields along the 1-dimensional
         differentiable manifold I mapped into the 2-dimensional differentiable
         manifold R^2
        sage: XIM.category()
        Category of finite dimensional modules
         over Algebra of differentiable scalar fields
         on the 1-dimensional differentiable manifold I

    The rank of the free module `\mathcal{X}(I,\Phi)` is the dimension
    of the manifold `\RR^2`, namely two::

        sage: XIM.rank()
        2

    A basis of it is induced by the coordinate vector frame of `\RR^2`::

        sage: XIM.bases()
        [Vector frame (I, (d/dx,d/dy)) with values on the 2-dimensional
         differentiable manifold R^2]

    Some elements of this module::

        sage: XIM.an_element().display()
        2 d/dx + 2 d/dy
        sage: v = XIM([t, t^2]) ; v
        Vector field along the 1-dimensional differentiable manifold I with
         values on the 2-dimensional differentiable manifold R^2
        sage: v.display()
        t d/dx + t^2 d/dy

    The test suite is passed::

        sage: TestSuite(XIM).run()

    Let us now consider the module of vector fields on the circle `S^1`; we
    start by constructing the `S^1` manifold::

        sage: M = Manifold(1, 'S^1')
        sage: U = M.open_subset('U')  # the complement of one point
        sage: c_t.<t> =  U.chart('t:(0,2*pi)') # the standard angle coordinate
        sage: V = M.open_subset('V') # the complement of the point t=pi
        sage: M.declare_union(U,V)   # S^1 is the union of U and V
        sage: c_u.<u> = V.chart('u:(0,2*pi)') # the angle t-pi
        sage: t_to_u = c_t.transition_map(c_u, (t-pi,), intersection_name='W',
        ....:                     restrictions1 = t!=pi, restrictions2 = u!=pi)
        sage: u_to_t = t_to_u.inverse()
        sage: W = U.intersection(V)

    `S^1` cannot be covered by a single chart, so it cannot be covered by
    a coordinate frame. It is however parallelizable and we introduce a global
    vector frame as follows. We notice that on their common subdomain, `W`,
    the coordinate vectors `\partial/\partial t` and `\partial/\partial u`
    coincide, as we can check explicitly::

        sage: c_t.frame()[0].display(c_u.frame().restrict(W))
        d/dt = d/du

    Therefore, we can extend `\partial/\partial t` to all `V` and hence to all
    `S^1`, to form a vector field on `S^1` whose components w.r.t. both
    `\partial/\partial t` and `\partial/\partial u` are 1::

        sage: e = M.vector_frame('e')
        sage: U.set_change_of_frame(e.restrict(U), c_t.frame(),
        ....:                       U.tangent_identity_field())
        sage: V.set_change_of_frame(e.restrict(V), c_u.frame(),
        ....:                       V.tangent_identity_field())
        sage: e[0].display(c_t.frame())
        e_0 = d/dt
        sage: e[0].display(c_u.frame())
        e_0 = d/du

    Equipped with the frame `e`, the manifold `S^1` is manifestly
    parallelizable::

        sage: M.is_manifestly_parallelizable()
        True

    Consequently, the module of vector fields on `S^1` is a free module::

        sage: XM = M.vector_field_module() ; XM
        Free module X(S^1) of vector fields on the 1-dimensional differentiable
         manifold S^1
        sage: isinstance(XM, FiniteRankFreeModule)
        True
        sage: XM.category()
        Category of finite dimensional modules
         over Algebra of differentiable scalar fields
         on the 1-dimensional differentiable manifold S^1
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    The zero element::

        sage: z = XM.zero() ; z
        Vector field zero on the 1-dimensional differentiable manifold S^1
        sage: z.display()
        zero = 0
        sage: z.display(c_t.frame())
        zero = 0

    The module `\mathcal{X}(S^1)` coerces to any module of vector fields
    defined on a subdomain of `S^1`, for instance `\mathcal{X}(U)`::

        sage: XU = U.vector_field_module() ; XU
        Free module X(U) of vector fields on the Open subset U of the
         1-dimensional differentiable manifold S^1
        sage: XU.has_coerce_map_from(XM)
        True
        sage: XU.coerce_map_from(XM)
        Conversion map:
          From: Free module X(S^1) of vector fields on the 1-dimensional differentiable manifold S^1
          To:   Free module X(U) of vector fields on the Open subset U of the 1-dimensional differentiable manifold S^1

    The conversion map is actually the restriction of vector fields defined
    on `S^1` to `U`.

    The Sage test suite for modules is passed::

        sage: TestSuite(XM).run()

    """

    Element = VectorFieldParal

    def __init__(self, domain, dest_map=None):
        r"""
        Construct the free module of vector fields with values on a
        parallelizable manifold.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: from sage.manifolds.differentiable.vectorfield_module \
            ....:                                  import VectorFieldFreeModule
            sage: XM = VectorFieldFreeModule(M, dest_map=M.identity_map()); XM
            Free module X(M) of vector fields on the 2-dimensional
             differentiable manifold M
            sage: XM is M.vector_field_module()
            True
            sage: TestSuite(XM).run()

        """
        from sage.manifolds.differentiable.scalarfield import DiffScalarField
        self._domain = domain
        if dest_map is None:
            dest_map = domain.identity_map()
        self._dest_map = dest_map
        self._ambient_domain = self._dest_map._codomain
        name = "X(" + domain._name
        latex_name = r"\mathcal{X}\left(" + domain._latex_name
        if self._dest_map == domain.identity_map():
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + self._dest_map._name + ")"
            latex_name += "," + self._dest_map._latex_name + r"\right)"
        manif = self._ambient_domain.manifold()
        cat = Modules(domain.scalar_field_algebra()).FiniteDimensional()
        FiniteRankFreeModule.__init__(self, domain.scalar_field_algebra(),
                               manif._dim, name=name, latex_name=latex_name,
                               start_index=manif._sindex,
                               output_formatter=DiffScalarField.coord_function,
                               category=cat)
        #
        # Special treatment when self._dest_map != identity:
        # bases of self are created from vector frames of the ambient domain
        #
        self._induced_bases = {}
        if self._dest_map != self._domain.identity_map():
            for frame in self._ambient_domain._top_frames:
                basis = self.basis(from_frame=frame)
                self._induced_bases[frame] = basis

        # Initialization of the components of the zero element:
        zero = self.zero()
        for frame in self._domain._frames:
            if frame._dest_map == self._dest_map:
                zero.add_comp(frame) # since new components are
                                     # initialized to zero

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: v = XM([-y,x], name='v'); v
            Vector field v on the 2-dimensional differentiable manifold M
            sage: v.display()
            v = -y d/dx + x d/dy
            sage: XM(0) is XM.zero()
            True

        """
        if comp == 0:
            return self.zero()
        if isinstance(comp, VectorField):
            if (self._domain.is_subset(comp._domain)
                   and self._ambient_domain.is_subset(comp._ambient_domain)):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {}".format(comp) +
                                 "to a vector field in {}".format(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(basis)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass FiniteRankFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from parent ``other``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: U = M.open_subset('U')
            sage: XM = M.vector_field_module()
            sage: XU = U.vector_field_module()
            sage: XM._coerce_map_from_(XU)
            False
            sage: XU._coerce_map_from_(XM)
            True

        """
        if isinstance(other, (VectorFieldModule, VectorFieldFreeModule)):
            return (self._domain.is_subset(other._domain)
                    and self._ambient_domain.is_subset(other._ambient_domain))
        else:
            return False

    #### End of parent methods

    #### Methods to be redefined by derived classes of FiniteRankFreeModule ####

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM._repr_()
            'Free module X(M) of vector fields on the 2-dimensional differentiable manifold M'
            sage: repr(XM)  # indirect doctest
            'Free module X(M) of vector fields on the 2-dimensional differentiable manifold M'
            sage: XM  # indirect doctest
            Free module X(M) of vector fields on the 2-dimensional
             differentiable manifold M

        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of vector fields "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {}".format(self._domain) + \
                           " mapped into the {}".format(self._ambient_domain)
        return description

    def domain(self):
        r"""
        Return the domain of the vector fields in ``self``.

        If the module is `\mathcal{X}(U, \Phi)`, returns the domain `U`
        of `\Phi`.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the domain of the vector fields that belong to this
          module

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.domain()
            3-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Y.<u,v> = U.chart()
            sage: Phi = U.diff_map(M, {(Y,X): [u+v, u-v, u*v]}, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.domain()
            2-dimensional differentiable manifold U

        """
        return self._domain

    def ambient_domain(self):
        r"""
        Return the manifold in which the vector fields of ``self``
        take their values.

        If the module is `\mathcal{X}(U, \Phi)`, returns the codomain `M`
        of `\Phi`.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the manifold in which the vector fields of ``self``
          take their values

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.ambient_domain()
            3-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Y.<u,v> = U.chart()
            sage: Phi = U.diff_map(M, {(Y,X): [u+v, u-v, u*v]}, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.ambient_domain()
            3-dimensional differentiable manifold M

        """
        return self._ambient_domain

    def destination_map(self):
        r"""
        Return the differential map associated to ``self``.

        The differential map associated to this module is the map

        .. MATH::

            \Phi:\  U \longrightarrow M

        such that this module is the set `\mathcal{X}(U,\Phi)` of all vector
        fields of the type

        .. MATH::

            v:\ U  \longrightarrow TM

        (where `TM` is the tangent bundle of `M`) such that

        .. MATH::

            \forall p \in U,\ v(p) \in T_{\Phi(p)} M,

        where `T_{\Phi(p)} M` is the tangent space to `M` at the
        point `\Phi(p)`.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the differential map `\Phi`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.destination_map()
            Identity map Id_M of the 3-dimensional differentiable manifold M
            sage: U = Manifold(2, 'U')
            sage: Y.<u,v> = U.chart()
            sage: Phi = U.diff_map(M, {(Y,X): [u+v, u-v, u*v]}, name='Phi')
            sage: XU = U.vector_field_module(dest_map=Phi)
            sage: XU.destination_map()
            Differentiable map Phi from the 2-dimensional differentiable
             manifold U to the 3-dimensional differentiable manifold M

        """
        return self._dest_map

    def tensor_module(self, k, l):
        r"""
        Return the free module of all tensors of type `(k, l)` defined
        on ``self``.

        INPUT:

        - ``k`` -- non-negative integer; the contravariant rank,
          the tensor type being `(k, l)`
        - ``l`` -- non-negative integer; the covariant rank,
          the tensor type being `(k, l)`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`
          representing the free module of type-`(k,l)` tensors on the
          vector field module

        EXAMPLES:

        A tensor field module on a 2-dimensional differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.tensor_module(1,2)
            Free module T^(1,2)(M) of type-(1,2) tensors fields on the
             2-dimensional differentiable manifold M

        The special case of tensor fields of type (1,0)::

            sage: XM.tensor_module(1,0)
            Free module X(M) of vector fields on the 2-dimensional
             differentiable manifold M

        The result is cached::

            sage: XM.tensor_module(1,2) is XM.tensor_module(1,2)
            True
            sage: XM.tensor_module(1,0) is XM
            True

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.tensorfield_module import TensorFieldFreeModule
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldFreeModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def dual_exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of the dual of ``self``.

        If the vector field module is `\mathcal{X}(U,\Phi)`, the
        `p`-th exterior power of its dual is the set `\Lambda^p(U, \Phi)`
        of `p`-forms along `U` with values on `\Phi(U)`. It is a module
        over `C^k(U)`, the ring (algebra) of differentiable scalar
        fields on `U`.

        INPUT:

        - ``p`` -- non-negative integer

        OUTPUT:

        - for `p \geq 1`, a
          :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormFreeModule`
          representing the module `\Lambda^p(U,\Phi)`; for `p=0`, the
          base ring, i.e. `C^k(U)`, is returned instead

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.dual_exterior_power(2)
            Free module /\^2(M) of 2-forms on the 2-dimensional differentiable
             manifold M
            sage: XM.dual_exterior_power(1)
            Free module /\^1(M) of 1-forms on the 2-dimensional differentiable
             manifold M
            sage: XM.dual_exterior_power(1) is XM.dual()
            True
            sage: XM.dual_exterior_power(0)
            Algebra of differentiable scalar fields on the 2-dimensional
             differentiable manifold M
            sage: XM.dual_exterior_power(0) is M.scalar_field_algebra()
            True

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormFreeModule`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.diff_form_module import DiffFormFreeModule
        if p == 0:
            return self._ring
        if p not in self._dual_exterior_powers:
            self._dual_exterior_powers[p] = DiffFormFreeModule(self, p)
        return self._dual_exterior_powers[p]

    def general_linear_group(self):
        r"""
        Return the general linear group of ``self``.

        If the vector field module is `\mathcal{X}(U,\Phi)`, the *general
        linear group* is the group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
        automorphisms of `\mathcal{X}(U,\Phi)`. Note that an automorphism of
        `\mathcal{X}(U,\Phi)` can also be viewed as a *field* along `U` of
        automorphisms of the tangent spaces of `V=\Phi(U)`.

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldParalGroup`
          representing `\mathrm{GL}(\mathcal{X}(U,\Phi))`

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.general_linear_group()
            General linear group of the Free module X(M) of vector fields on
             the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldParalGroup`
            for more examples and documentation.

        """
        from sage.manifolds.differentiable.automorphismfield_group import AutomorphismFieldParalGroup
        return AutomorphismFieldParalGroup(self)

    def basis(self, symbol=None, latex_symbol=None, from_frame=None):
        r"""
        Define a basis of ``self``.

        A basis of the vector field module is actually a vector frame along
        the differentiable manifold `U` over which the vector field module
        is defined.

        If the basis specified by the given symbol already exists, it is
        simply returned.
        If no argument is provided the module's default basis is returned.

        INPUT:

        - ``symbol`` -- (string; default: ``None``) a letter (of a few
          letters) to denote a generic element of the basis; if ``None``
          and ``from_frame = None`` the module's default basis is returned
        - ``latex_symbol`` -- (string; default: ``None``) symbol to denote a
          generic element of the basis; if ``None``, the value of ``symbol``
          is used
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde{e}`
          on the codomain `M` of the destination map `\Phi` of ``self``;
          the returned basis `e` is then such that for all `p \in U`,
          we have `e(p) = \tilde{e}(\Phi(p))`

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          representing a basis on ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: e = XM.basis('e'); e
            Vector frame (M, (e_0,e_1))

        See :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
        for more examples and documentation.

        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        if symbol is None:
            if from_frame is None:
                return self.default_basis()
            else:
                symbol = from_frame._symbol
                latex_symbol = from_frame._latex_symbol
        for other in self._known_bases:
            if symbol == other._symbol:
                return other
        return VectorFrame(self, symbol=symbol, latex_symbol=latex_symbol,
                           from_frame=from_frame)

    def tensor(self, tensor_type, name=None, latex_name=None, sym=None,
               antisym=None, specific_type=None):
        r"""
        Construct a tensor on ``self``.

        The tensor is actually a tensor field along the differentiable
        manifold `U` over which ``self`` is defined.

        INPUT:

        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank
          and l the covariant rank
        - ``name`` -- (string; default: ``None``) name given to the tensor
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the tensor; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention position=0 for the first argument; for instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments

        - ``antisym`` -- (default: ``None``) antisymmetry or list of
          antisymmetries among the arguments, with the same convention
          as for ``sym``
        - ``specific_type`` -- (default: ``None``) specific subclass of
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          for the output

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          representing the tensor defined on ``self`` with the provided
          characteristics

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.tensor((1,2), name='t')
            Tensor field t of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: XM.tensor((1,0), name='a')
            Vector field a on the 2-dimensional differentiable manifold M
            sage: XM.tensor((0,2), name='a', antisym=(0,1))
            2-form a on the 2-dimensional differentiable manifold M

        See
        :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
        for more examples and documentation.

        """
        from sage.manifolds.differentiable.automorphismfield import (AutomorphismField,
                                                                     AutomorphismFieldParal)
        from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
        if tensor_type == (1,0):
            return self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            return self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type == (1,1) and specific_type is not None:
            if issubclass(specific_type,
                          (AutomorphismField, AutomorphismFieldParal)):
                return self.automorphism(name=name, latex_name=latex_name)
        elif (tensor_type[0] == 0 and tensor_type[1] > 1
              and antisym is not None and antisym != []):
            if isinstance(antisym, list):
                antisym0 = antisym[0]
            else:
                antisym0 = antisym
            if len(antisym0) == tensor_type[1]:
                return self.alternating_form(tensor_type[1], name=name,
                                             latex_name=latex_name)
            else:
                return self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name,
                                 sym=sym, antisym=antisym)
        elif tensor_type==(0,2) and specific_type is not None:
            if issubclass(specific_type, PseudoRiemannianMetric):
                return self.metric(name, latex_name=latex_name)
                # NB: the signature is not treated
        # Generic case
        return self.tensor_module(*tensor_type).element_class(self,
                        tensor_type, name=name, latex_name=latex_name, sym=sym,
                        antisym=antisym)

    def tensor_from_comp(self, tensor_type, comp, name=None, latex_name=None):
        r"""
        Construct a tensor on ``self`` from a set of components.

        The tensor is actually a tensor field along the differentiable
        manifold `U` over which the vector field module is defined.
        The tensor symmetries are deduced from those of the components.

        INPUT:

        - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant
          rank and `l` the covariant rank
        - ``comp`` -- :class:`~sage.tensor.modules.comp.Components`;
          the tensor components in a given basis
        - ``name`` -- string (default: ``None``); name given to the tensor
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the tensor; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          representing the tensor defined on the vector field module with the
          provided characteristics

        EXAMPLES:

        A 2-dimensional set of components transformed into a type-`(1,1)`
        tensor field::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: XM = M.vector_field_module()
            sage: from sage.tensor.modules.comp import Components
            sage: comp = Components(M.scalar_field_algebra(), X.frame(), 2,
            ....:                   output_formatter=XM._output_formatter)
            sage: comp[:] = [[1+x, -y], [x*y, 2-y^2]]
            sage: t = XM.tensor_from_comp((1,1), comp, name='t'); t
            Tensor field t of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: t.display()
            t = (x + 1) d/dx*dx - y d/dx*dy + x*y d/dy*dx + (-y^2 + 2) d/dy*dy

        The same set of components transformed into a type-`(0,2)`
        tensor field::

            sage: t = XM.tensor_from_comp((0,2), comp, name='t'); t
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: t.display()
            t = (x + 1) dx*dx - y dx*dy + x*y dy*dx + (-y^2 + 2) dy*dy

        """
        from sage.tensor.modules.comp import CompWithSym, CompFullySym, CompFullyAntiSym

        # 0/ Compatibility checks:
        if comp._ring is not self._ring:
             raise ValueError("the components are not defined on the same " +
                              "ring as the module")
        if comp._frame not in self._known_bases:
            raise ValueError("the components are not defined on a basis of " +
                             "the module")
        if comp._nid != tensor_type[0] + tensor_type[1]:
            raise ValueError("number of component indices not compatible " +
                             "with the tensor type")
        #
        # 1/ Construction of the tensor:
        if tensor_type == (1,0):
            resu = self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            resu = self.linear_form(name=name, latex_name=latex_name)
        elif (tensor_type[0] == 0 and tensor_type[1] > 1
              and isinstance(comp, CompFullyAntiSym)):
            resu = self.alternating_form(tensor_type[1], name=name,
                                         latex_name=latex_name)
        else:
            resu = self.tensor_module(*tensor_type).element_class(self,
                             tensor_type, name=name, latex_name=latex_name)
            # Tensor symmetries deduced from those of comp:
            if isinstance(comp, CompWithSym):
                resu._sym = comp._sym
                resu._antisym = comp._antisym

        # 2/ Tensor components set to comp:
        resu._components[comp._frame] = comp
        #
        return resu

    def sym_bilinear_form(self, name=None, latex_name=None):
        r"""
        Construct a symmetric bilinear form on ``self``.

        A symmetric bilinear form on the vector field module is
        actually a field of tangent-space symmetric bilinear forms
        along the differentiable manifold `U` over which the vector
        field module is defined.

        INPUT:

        - ``name`` -- string (default: ``None``); name given to the
          automorphism
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the automorphism; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
          of tensor type `(0,2)` and symmetric

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.sym_bilinear_form(name='a')
            Field of symmetric bilinear forms a on the 2-dimensional
             differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
            for more examples and documentation.

        """
        return self.tensor((0,2), name=name, latex_name=latex_name, sym=(0,1))

    #### End of methods to be redefined by derived classes of FiniteRankFreeModule ####

    def metric(self, name, signature=None, latex_name=None):
        r"""
        Construct a pseudo-Riemannian metric (nondegenerate symmetric bilinear
        form) on the current vector field module.

        A pseudo-Riemannian metric of the vector field module is actually a
        field of tangent-space non-degenerate symmetric bilinear forms along
        the manifold `U` on which the vector field module is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (integer; default: ``None``) signature `S` of the
          metric: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal
          writing of the metric components; if ``signature`` is not provided,
          `S` is set to the manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the metric; if ``None``, it is formed from ``name``

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetricParal`
          representing the defined pseudo-Riemannian metric.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()
            sage: XM.metric('g')
            Riemannian metric g on the 2-dimensional differentiable manifold M
            sage: XM.metric('g', signature=0)
            Lorentzian metric g on the 2-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetricParal`
            for more documentation.

        """
        from sage.manifolds.differentiable.metric import PseudoRiemannianMetricParal
        return PseudoRiemannianMetricParal(self, name, signature=signature,
                                           latex_name=latex_name)
