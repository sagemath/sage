# sage_setup: distribution = sagemath-symbolics
r"""
Differentiable Vector Bundles

Let `K` be a topological field. A `C^k`-differentiable *vector bundle* of rank
`n` over the field `K` and over a `C^k`-differentiable manifold `M` (base
space) is a `C^k`-differentiable manifold `E` (total space) together with a
`C^k` differentiable and surjective map `\pi: E \to M` such that for
every point `x \in M`:

- the set `E_x=\pi^{-1}(x)` has the vector space structure of `K^n`,
- there is a neighborhood `U \subset M` of `x` and a `C^k`-diffeomorphism
  `\varphi: \pi^{-1}(x) \to U \times K^n` such that
  `v \mapsto \varphi^{-1}(y,v)` is a linear isomorphism for any `y \in U`.

An important case of a differentiable vector bundle over a differentiable
manifold is the tensor bundle (see :class:`TensorBundle`)

AUTHORS:

- Michael Jung (2019) : initial version
"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.vector_bundles import VectorBundles
from sage.rings.cc import CC
from sage.rings.real_mpfr import RR
from sage.manifolds.vector_bundle import TopologicalVectorBundle
from sage.rings.infinity import infinity
from sage.misc.superseded import deprecated_function_alias
from sage.rings.rational_field import QQ

class DifferentiableVectorBundle(TopologicalVectorBundle):
    r"""
    An instance of this class represents a differentiable vector bundle
    `E \to M`

    INPUT:

    - ``rank`` -- positive integer; rank of the vector bundle
    - ``name`` -- string representation given to the total space
    - ``base_space`` -- the base space (differentiable manifold) `M` over which
      the vector bundle is defined
    - ``field`` -- field `K` which gives the fibers the structure of a
      vector space over `K`; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
        a vector bundle over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
        for a vector bundle over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of topological fields

    - ``latex_name`` -- (default: ``None``) LaTeX representation given to the
      total space
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``VectorBundles(base_space, c_field).Differentiable()`` is
      assumed (see the category
      :class:`~sage.categories.vector_bundles.VectorBundles`)

    EXAMPLES:

    A differentiable vector bundle of rank 2 over a 3-dimensional
    differentiable manifold::

        sage: M = Manifold(3, 'M')
        sage: E = M.vector_bundle(2, 'E', field='complex'); E
        Differentiable complex vector bundle E -> M of rank 2 over the base
         space 3-dimensional differentiable manifold M
        sage: E.category()
        Category of smooth vector bundles over Complex Field with 53 bits of
         precision with base space 3-dimensional differentiable manifold M

    At this stage, the differentiable vector bundle has the same
    differentiability degree as the base manifold::

        sage: M.diff_degree() == E.diff_degree()
        True
    """
    def __init__(self, rank, name, base_space, field='real', latex_name=None,
                 category=None, unique_tag=None):
        r"""
        Construct a differentiable vector bundle.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: from sage.manifolds.differentiable.vector_bundle import DifferentiableVectorBundle
            sage: DifferentiableVectorBundle(2, 'E', M)
            Differentiable real vector bundle E -> M of rank 2 over the base
             space 2-dimensional differentiable manifold M
        """
        diff_degree = base_space._diff_degree
        if category is None:
            if field == 'real':
                field_c = RR
            elif field == 'complex':
                field_c = CC
            else:
                field_c = field
            if diff_degree == infinity:
                category = VectorBundles(base_space, field_c).Smooth()
            else:
                category = VectorBundles(base_space, field_c).Differentiable()
        TopologicalVectorBundle.__init__(self, rank, name, base_space,
                                         field=field,
                                         latex_name=latex_name,
                                         category=category)
        self._diff_degree = diff_degree  # Override diff degree

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._repr_()
            'Differentiable real vector bundle E -> M of rank 1 over the base
             space 2-dimensional differentiable manifold M'
        """
        desc = "Differentiable "
        return desc + TopologicalVectorBundle._repr_object_name(self)

    def bundle_connection(self, name, latex_name=None):
        r"""
        Return a bundle connection on ``self``.

        OUTPUT:

        - a bundle connection on ``self`` as an instance of
          :class:`~sage.manifolds.differentiable.bundle_connection.BundleConnection`

        EXAMPLES::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla'); nab
            Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 2 over the base space 3-dimensional differentiable
             manifold M

        .. SEEALSO::

            Further examples can be found in
            :class:`~sage.manifolds.differentiable.bundle_connection.BundleConnection`.
        """
        from .bundle_connection import BundleConnection
        return BundleConnection(self, name, latex_name)

    def characteristic_cohomology_class_ring(self, base=QQ):
        r"""
        Return the characteristic cohomology class ring of ``self`` over
        a given base.

        INPUT:

        - ``base`` -- (default: ``QQ``) base over which the ring should be
          constructed; typically that would be `\ZZ`, `\QQ`, `\RR` or the
          symbolic ring

        EXAMPLES::

            sage: M = Manifold(4, 'M', start_index=1)
            sage: R = M.tangent_bundle().characteristic_cohomology_class_ring()
            sage: R
            Algebra of characteristic cohomology classes of the Tangent bundle
             TM over the 4-dimensional differentiable manifold M
            sage: p1 = R.gen(0); p1
            Characteristic cohomology class (p_1)(TM) of the Tangent bundle TM
             over the 4-dimensional differentiable manifold M
            sage: 1 + p1
            Characteristic cohomology class (1 + p_1)(TM) of the Tangent bundle
             TM over the 4-dimensional differentiable manifold M
        """
        from .characteristic_cohomology_class import CharacteristicCohomologyClassRing

        return CharacteristicCohomologyClassRing(base, self)

    def characteristic_cohomology_class(self, *args, **kwargs):
        r"""
        Return a characteristic cohomology class associated with the input
        data.

        INPUT:

        - ``val`` -- the input data associated with the characteristic class
          using the Chern-Weil homomorphism; this argument can be either a
          symbolic expression, a polynomial or one of the following predefined
          classes:

          - ``'Chern'`` -- total Chern class,
          - ``'ChernChar'`` -- Chern character,
          - ``'Todd'`` -- Todd class,
          - ``'Pontryagin'`` -- total Pontryagin class,
          - ``'Hirzebruch'`` -- Hirzebruch class,
          - ``'AHat'`` -- `\hat{A}` class,
          - ``'Euler'`` -- Euler class.

        - ``base_ring`` -- (default: ``QQ``) base ring over which the
          characteristic cohomology class ring shall be defined
        - ``name`` -- (default: ``None``) string representation given to the
          characteristic cohomology class; if ``None`` the default algebra
          representation or predefined name is used
        - ``latex_name`` -- (default: ``None``) LaTeX name given to the
          characteristic class; if ``None`` the value of ``name`` is used
        - ``class_type`` -- (default: ``None``) class type of the characteristic
          cohomology class; the following options are possible:

          - ``'multiplicative'`` -- returns a class of multiplicative type
          - ``'additive'`` -- returns a class of additive type
          - ``'Pfaffian'`` -- returns a class of Pfaffian type

          This argument must be stated if ``val`` is a polynomial or symbolic
          expression.

        EXAMPLES:

        Pontryagin class on the Minkowski space::

            sage: M = Manifold(4, 'M', structure='Lorentzian', start_index=1)
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.metric()
            sage: g[1,1] = -1
            sage: g[2,2] = 1
            sage: g[3,3] = 1
            sage: g[4,4] = 1
            sage: g.display()
            g = -dt⊗dt + dx⊗dx + dy⊗dy + dz⊗dz

        Let us introduce the corresponding Levi-Civita connection::

            sage: nab = g.connection(); nab
            Levi-Civita connection nabla_g associated with the Lorentzian
             metric g on the 4-dimensional Lorentzian manifold M
            sage: nab.set_immutable()  # make nab immutable

        Of course, `\nabla_g` is flat::

            sage: nab.display()

        Let us check the total Pontryagin class which must be the one
        element in the corresponding cohomology ring in this case::

            sage: TM = M.tangent_bundle(); TM
            Tangent bundle TM over the 4-dimensional Lorentzian manifold M
            sage: p = TM.characteristic_cohomology_class('Pontryagin'); p
            Characteristic cohomology class p(TM) of the Tangent bundle TM over
             the 4-dimensional Lorentzian manifold M
            sage: p_form = p.get_form(nab); p_form.display_expansion()
            p(TM, nabla_g) = 1

        .. SEEALSO::

            More examples can be found in
            :class:`~sage.manifolds.differentiable.characteristic_class.CharacteristicClass`.
        """
        base_ring = kwargs.get('base_ring', QQ)
        R = self.characteristic_cohomology_class_ring(base_ring)
        return R(*args, **kwargs)

    characteristic_class = deprecated_function_alias(29581, characteristic_cohomology_class)

    def diff_degree(self):
        r"""
        Return the vector bundle's degree of differentiability.

        The degree of differentiability is the integer `k` (possibly
        `k=\infty`) such that the vector bundle is of class `C^k` over
        its base field. The degree always corresponds to the degree of
        differentiability of it's base space.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.diff_degree()
            +Infinity
            sage: M = Manifold(2, 'M', structure='differentiable',
            ....:              diff_degree=3)
            sage: E = M.vector_bundle(2, 'E')
            sage: E.diff_degree()
            3
        """
        return self._diff_degree

    def total_space(self):
        r"""
        Return the total space of ``self``.

        .. NOTE::

            At this stage, the total space does not come with induced charts.

        OUTPUT:

        - the total space of ``self`` as an instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.total_space()
            6-dimensional differentiable manifold E
        """
        if self._total_space is None:
            from sage.manifolds.manifold import Manifold
            base_space = self._base_space
            dim = base_space._dim * self._rank
            sindex = base_space.start_index()
            self._total_space = Manifold(dim, self._name,
                                latex_name=self._latex_name,
                                field=self._field, structure='differentiable',
                                diff_degree=self._diff_degree,
                                start_index=sindex)

        # TODO: if update_atlas: introduce charts via self._atlas

        return self._total_space

# *****************************************************************************

class TensorBundle(DifferentiableVectorBundle):
    r"""
    Tensor bundle over a differentiable manifold along a differentiable map.

    An instance of this class represents the pullback tensor bundle
    `\Phi^* T^{(k,l)}N` along a differentiable map (called *destination map*)

    .. MATH::

        \Phi: M \longrightarrow N

    between two differentiable manifolds `M` and `N` over the topological field
    `K`.

    More precisely, `\Phi^* T^{(k,l)}N` consists of all pairs
    `(p,t) \in M \times T^{(k,l)}N` such that `t \in T_q^{(k,l)}N` for
    `q = \Phi(p)`, namely

    .. MATH::

        t:\ \underbrace{T_q^*N\times\cdots\times T_q^*N}_{k\ \; \text{times}}
            \times \underbrace{T_q N\times\cdots\times T_q N}_{l\ \; \text{times}}
            \longrightarrow K

    (`k` is called the *contravariant* and `l` the *covariant* rank of the
    tensor bundle).

    The trivializations are directly given by charts on the codomain (called
    *ambient domain*) of `\Phi`.
    In particular, let `(V, \varphi)` be a chart of `N` with components
    `(x^1, \dots, x^n)` such that `q=\Phi(p) \in V`. Then, the matrix entries
    of `t \in T_q^{(k,l)}N` are given by

    .. MATH::

        t^{a_1 \ldots a_k}_{\phantom{a_1 \ldots a_k} \, b_1 \ldots b_l} =
            t \left( \left.\frac{\partial}{\partial x^{a_1}}\right|_q, \dots,
            \left.\frac{\partial}{\partial x^{a_k}}\right|_q,
            \left.\mathrm{d}x^{b_1}\right|_q, \dots,
            \left.\mathrm{d}x^{b_l}\right|_q \right) \in K

    and a trivialization over `U=\Phi^{-1}(V) \subset M` is obtained via

    .. MATH::

        (p,t) \mapsto \left(p, t^{1 \ldots 1}_{\phantom{1 \ldots 1} \, 1 \ldots 1},
            \dots, t^{n \ldots n}_{\phantom{n \ldots n} \, n \ldots n} \right)
            \in U \times K^{n^{(k+l)}}.

    The standard case of a tensor bundle over a differentiable manifold
    corresponds to `M=N` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `N` (`M` is then an
    open interval of `\RR`).

    INPUT:

    - ``base_space`` -- the base space (differentiable manifold) `M` over which
      the tensor bundle is defined
    - ``k`` -- the contravariant rank of the corresponding tensor bundle
    - ``l`` -- the covariant rank of the corresponding tensor bundle
    - ``dest_map`` -- (default: ``None``) destination map
      `\Phi:\ M \rightarrow N`
      (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`); if
      ``None``, it is assumed that `M=M` and `\Phi` is the identity map of
      `M` (case of the standard tensor bundle over `M`)

    EXAMPLES:

    Pullback tangent bundle of `R^2` along a curve `\Phi`::

        sage: M = Manifold(2, 'M')
        sage: c_cart.<x,y> = M.chart()
        sage: R = Manifold(1, 'R')
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 2-dimensional differentiable manifold M
        sage: Phi.display()
        Phi: R → M
           t ↦ (x, y) = (cos(t), sin(t))
        sage: PhiTM = R.tangent_bundle(dest_map=Phi); PhiTM
        Tangent bundle Phi^*TM over the 1-dimensional differentiable manifold R
         along the Differentiable map Phi from the 1-dimensional differentiable
         manifold R to the 2-dimensional differentiable manifold M

    The section module is the corresponding tensor field module::

        sage: R_tensor_module = R.tensor_field_module((1,0), dest_map=Phi)
        sage: R_tensor_module is PhiTM.section_module()
        True
    """
    def __init__(self, base_space, k, l, dest_map=None):
        r"""
        Construct a tensor bundle.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: N = Manifold(2, 'N')
            sage: Phi = M.diff_map(N, name='Phi')
            sage: from sage.manifolds.differentiable.vector_bundle import TensorBundle
            sage: TensorBundle(M, 1, 2, dest_map=Phi)
            Tensor bundle Phi^*T^(1,2)N over the 2-dimensional differentiable
             manifold M along the Differentiable map Phi from the 2-dimensional
             differentiable manifold M to the 2-dimensional differentiable
             manifold N
        """
        if dest_map is None:
            self._dest_map = base_space.identity_map()
        else:
            self._dest_map = dest_map
        self._ambient_domain = self._dest_map._codomain
        self._tensor_type = (k, l)
        # Set total space name:
        if not self._dest_map.is_identity():
            if self._dest_map._name is None:
                name = "(unnamed map)^*"
            else:
                name = self._dest_map._name + "^*"
            if self._dest_map._latex_name is None:
                latex_name = r'\text{(unnamed map)}^* '
            else:
                latex_name = self._dest_map._latex_name + r'^* '
        else:
            name = ""
            latex_name = ""
        if self._tensor_type == (1, 0):
            name += "T{}".format(self._ambient_domain._name)
            latex_name += r'T{}'.format(self._ambient_domain._latex_name)
        elif self._tensor_type == (0, 1):
            name += "T*{}".format(self._ambient_domain._name)
            latex_name += r'T^*{}'.format(self._ambient_domain._latex_name)
        else:
            name += "T^({},{}){}".format(k, l, self._ambient_domain._name)
            latex_name += r'T^{(' + str(k) + r',' + str(l) + r')}' + \
                          self._ambient_domain._latex_name
        # Initialize differentiable vector bundle:
        rank = self._ambient_domain.dim() ** (k + l)
        DifferentiableVectorBundle.__init__(self, rank, name, base_space,
                                            field=base_space._field,
                                            latex_name=latex_name)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: TM = M.tangent_bundle()
            sage: TM._init_derived()
        """
        self._def_frame = None

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: TM = M.tangent_bundle()
            sage: TM # indirect doctest
            Tangent bundle TM over the 2-dimensional differentiable manifold M
            sage: repr(TM) # indirect doctest
            'Tangent bundle TM over the 2-dimensional differentiable manifold M'
            sage: TM._repr_()
            'Tangent bundle TM over the 2-dimensional differentiable manifold M'
            sage: cTM = M.cotangent_bundle()
            sage: cTM._repr_()
            'Cotangent bundle T*M over the 2-dimensional differentiable
             manifold M'
            sage: T12M = M.tensor_bundle(1, 2)
            sage: T12M._repr_()
            'Tensor bundle T^(1,2)M over the 2-dimensional differentiable
             manifold M'
        """
        if self._tensor_type == (1, 0):
            desc = "Tangent bundle "
        elif self._tensor_type == (0, 1):
            desc = "Cotangent bundle "
        else:
            desc = "Tensor bundle "
        desc += self._name + " over the {}".format(self._base_space)
        if not self._dest_map.is_identity():
            desc += " along the {}".format(self._dest_map)
        return desc

    def fiber(self, point):
        r"""
        Return the tensor bundle fiber over a point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` of the base manifold of ``self``

        OUTPUT:

        - an instance of :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
          representing the tensor bundle fiber over `p`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,2,1), name='p'); p
            Point p on the 3-dimensional differentiable manifold M
            sage: TM = M.tangent_bundle(); TM
            Tangent bundle TM over the 3-dimensional differentiable manifold M
            sage: TM.fiber(p)
            Tangent space at Point p on the 3-dimensional differentiable
             manifold M
            sage: TM.fiber(p) is M.tangent_space(p)
            True

        ::

            sage: T11M = M.tensor_bundle(1,1); T11M
            Tensor bundle T^(1,1)M over the 3-dimensional differentiable
             manifold M
            sage: T11M.fiber(p)
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 3-dimensional differentiable manifold M
            sage: T11M.fiber(p) is M.tangent_space(p).tensor_module(1,1)
            True
        """
        amb_point = self._dest_map(point)
        return self._ambient_domain.tangent_space(amb_point).tensor_module(*self._tensor_type)

    def atlas(self):
        r"""
        Return the list of charts that have been defined on the codomain of the
        destination map.

        .. NOTE::

            Since an atlas of charts gives rise to an atlas of trivializations,
            this method directly invokes
            :meth:`~sage.manifolds.manifold.TopologicalManifold.atlas`
            of the class
            :class:`~sage.manifolds.manifold.TopologicalManifold`.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: Y.<u,v> = M.chart()
            sage: TM = M.tangent_bundle()
            sage: TM.atlas()
            [Chart (M, (x, y)), Chart (M, (u, v))]
        """
        return self._base_space.atlas()

    def section_module(self, domain=None):
        r"""
        Return the section module on ``domain``, namely the corresponding
        tensor field module, of ``self`` on ``domain``.

        .. NOTE::

            This method directly invokes
            :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.tensor_field_module`
            of the class
            :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`.

        INPUT:

        - ``domain`` -- (default: ``None``) the domain of the corresponding
          section module; if ``None``, the base space is assumed

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`)
          representing the module `\mathcal{T}^{(k,l)}(U,\Phi)` of type-`(k,l)`
          tensor fields on the domain `U \subset M` taking values on
          `\Phi(U) \subset N`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U')
            sage: TM = M.tangent_bundle()
            sage: TUM = TM.section_module(domain=U); TUM
            Module X(U) of vector fields on the Open subset U of the
             2-dimensional differentiable manifold M
            sage: TUM is U.tensor_field_module((1,0))
            True
        """
        if domain is None:
            base_space = self.base_space()
            return base_space.tensor_field_module(self._tensor_type,
                                                  dest_map=self._dest_map)
        else:
            return domain.tensor_field_module(self._tensor_type,
                                      dest_map=self._dest_map.restrict(domain))

    def section(self, *args, **kwargs):
        r"""
        Return a section of ``self`` on ``domain``, namely a tensor field on
        the subset ``domain`` of the base space.

        .. NOTE::

            This method directly invokes
            :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.tensor_field`
            of the class
            :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`.

        INPUT:

        - ``comp`` -- (optional) either the components of the tensor field
          with respect to the vector frame specified by the argument
          ``frame`` or a dictionary of components, the keys of which are vector
          frames or pairs ``(f, c)`` where ``f`` is a vector frame and ``c``
          the chart in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``domain`` -- (default: ``None``) domain of the section; if ``None``,
          ``self.base_space()`` is assumed
        - ``name`` -- (default: ``None``) name given to the tensor field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          tensor field; if ``None``, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position=0`` for the first argument; for instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments

        - ``antisym`` -- (default: ``None``) antisymmetry or list of
          antisymmetries among the arguments, with the same convention as for
          ``sym``

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`)
          representing the defined tensor field on the domain `U \subset M`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                              intersection_name='W',
            ....:                              restrictions1= x>0,
            ....:                              restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: T11M = M.tensor_bundle(1, 1); T11M
            Tensor bundle T^(1,1)M over the 2-dimensional differentiable
             manifold M
            sage: t = T11M.section({eU: [[1, x], [0, 2]]}, name='t'); t
            Tensor field t of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: t.display()
            t = ∂/∂x⊗dx + x ∂/∂x⊗dy + 2 ∂/∂y⊗dy

        An example of use with the arguments ``comp`` and ``domain``::

            sage: TM = M.tangent_bundle()
            sage: w = TM.section([-y, x], domain=U); w
            Vector field on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: w.display()
            -y ∂/∂x + x ∂/∂y
        """
        nargs = [self._tensor_type[0], self._tensor_type[1]]
        nargs.extend(args)
        domain = kwargs.pop('domain', self._base_space)
        kwargs['dest_map'] = self._dest_map.restrict(domain)
        return domain.tensor_field(*nargs, **kwargs)

    def set_change_of_frame(self, frame1, frame2, change_of_frame,
                            compute_inverse=True):
        r"""
        Relate two vector frames by an automorphism.

        This updates the internal dictionary ``self._frame_changes`` of the
        base space `M`.

        .. SEEALSO::

            For further details on frames on ``self`` see
            :meth:`local_frame`.

        .. NOTE::

            Since frames on ``self`` are directly induced by vector frames on
            the base space, this method directly invokes
            :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_change_of_frame`
            of the class
            :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`.

        INPUT:

        - ``frame1`` -- frame 1, denoted `(e_i)` below
        - ``frame2`` -- frame 2, denoted `(f_i)` below
        - ``change_of_frame`` -- instance of class
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` -- boolean (default: ``True``); if set to ``True``, the inverse
          automorphism is computed and the change from basis `(f_i)` to `(e_i)`
          is set to it in the internal dictionary ``self._frame_changes``

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: f = M.vector_frame('f')
            sage: a = M.automorphism_field()
            sage: a[e,:] = [[1,2],[0,3]]
            sage: TM = M.tangent_bundle()
            sage: TM.set_change_of_frame(e, f, a)
            sage: f[0].display(e)
            f_0 = e_0
            sage: f[1].display(e)
            f_1 = 2 e_0 + 3 e_1
            sage: e[0].display(f)
            e_0 = f_0
            sage: e[1].display(f)
            e_1 = -2/3 f_0 + 1/3 f_1
            sage: TM.change_of_frame(e,f)[e,:]
            [1 2]
            [0 3]
        """
        if not frame1._domain.is_subset(self._ambient_domain):
            raise ValueError("the frames must be defined on a subset of "
                             "the {}".format(self._ambient_domain))
        frame1._domain.set_change_of_frame(frame1=frame1, frame2=frame2,
                                           change_of_frame=change_of_frame,
                                           compute_inverse=compute_inverse)

    def change_of_frame(self, frame1, frame2):
        r"""
        Return a change of vector frames defined on the base space of ``self``.

        .. SEEALSO::

            For further details on frames on ``self`` see
            :meth:`local_frame`.

        .. NOTE::

            Since frames on ``self`` are directly induced by vector frames on
            the base space, this method directly invokes
            :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.change_of_frame`
            of the class
            :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`.

        INPUT:

        - ``frame1`` -- local frame 1
        - ``frame2`` -- local frame 2

        OUTPUT:

        - a :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          representing, at each point, the vector space automorphism `P` that
          relates frame 1, `(e_i)` say, to frame 2, `(f_i)` say, according to
          `f_i = P(e_i)`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y))
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: TM = M.tangent_bundle()
            sage: TM.change_of_frame(c_xy.frame(), c_uv.frame())
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: TM.change_of_frame(c_xy.frame(), c_uv.frame())[:]
            [ 1/2  1/2]
            [ 1/2 -1/2]
            sage: TM.change_of_frame(c_uv.frame(), c_xy.frame())
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: TM.change_of_frame(c_uv.frame(), c_xy.frame())[:]
            [ 1  1]
            [ 1 -1]
            sage: TM.change_of_frame(c_uv.frame(), c_xy.frame()) == \
            ....:       M.change_of_frame(c_xy.frame(), c_uv.frame()).inverse()
            True
        """
        return self._base_space.change_of_frame(frame1=frame1, frame2=frame2)

    def changes_of_frame(self):
        r"""
        Return the changes of vector frames defined on the base space of
        ``self`` with respect to the destination map.

        .. SEEALSO::

            For further details on frames on ``self`` see
            :meth:`local_frame`.

        OUTPUT:

        - dictionary of automorphisms on the tangent bundle representing
          the changes of frames, the keys being the pair of frames

        EXAMPLES:

        Let us consider a first vector frame on a 2-dimensional
        differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: TM = M.tangent_bundle()
            sage: e = X.frame(); e
            Coordinate frame (M, (∂/∂x,∂/∂y))

        At this stage, the dictionary of changes of frame is empty::

            sage: TM.changes_of_frame()
            {}

        We introduce a second frame on the manifold, relating it to
        frame ``e`` by a field of tangent space automorphisms::

            sage: a = M.automorphism_field(name='a')
            sage: a[:] = [[-y, x], [1, 2]]
            sage: f = e.new_frame(a, 'f'); f
            Vector frame (M, (f_0,f_1))

        Then we have::

            sage: TM.changes_of_frame()  # random (dictionary output)
            {(Coordinate frame (M, (∂/∂x,∂/∂y)),
              Vector frame (M, (f_0,f_1))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M,
             (Vector frame (M, (f_0,f_1)),
              Coordinate frame (M, (∂/∂x,∂/∂y))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M}

        Some checks::

            sage: TM.changes_of_frame()[(e,f)] == a
            True
            sage: TM.changes_of_frame()[(f,e)] == a^(-1)
            True
        """
        base_cof = self._base_space.changes_of_frame()
        # Filter out all frames with respect to dest_map:
        cof = {}
        for frames in base_cof:
            if frames[0]._dest_map == self._dest_map:
                cof[(frames[0], frames[1])] = base_cof[frames]
        return cof

    def frames(self):
        r"""
        Return the list of all vector frames defined on the base space of
        ``self`` with respect to the destination map.

        .. SEEALSO::

            For further details on frames on ``self`` see
            :meth:`local_frame`.

        OUTPUT: list of local frames defined on ``self``

        EXAMPLES:

        Vector frames on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: TM = M.tangent_bundle()
            sage: TM.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y))]
            sage: e = TM.vector_frame('e')
            sage: TM.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y)),
             Vector frame (R^2, (e_0,e_1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1})
            sage: TU = U.tangent_bundle()
            sage: TU.frames()
            [Coordinate frame (U, (∂/∂x,∂/∂y))]
            sage: TM.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y)),
             Vector frame (R^2, (e_0,e_1)),
             Coordinate frame (U, (∂/∂x,∂/∂y))]

        List of vector frames of a tensor bundle of type `(1 ,1)` along a
        curve::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: e_cart = c_cart.frame() # standard basis
            sage: R = Manifold(1, 'R')
            sage: T.<t> = R.chart()  # canonical chart on R
            sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold R to the 2-dimensional differentiable manifold M
            sage: Phi.display()
            Phi: R → M
               t ↦ (x, y) = (cos(t), sin(t))
            sage: PhiT11 = R.tensor_bundle(1, 1, dest_map=Phi); PhiT11
            Tensor bundle Phi^*T^(1,1)M over the 1-dimensional differentiable
             manifold R along the Differentiable map Phi from the 1-dimensional
             differentiable manifold R to the 2-dimensional differentiable
             manifold M
            sage: f = PhiT11.local_frame(); f
            Vector frame (R, (∂/∂x,∂/∂y)) with values on the 2-dimensional
             differentiable manifold M
            sage: PhiT11.frames()
            [Vector frame (R, (∂/∂x,∂/∂y)) with values on the 2-dimensional
             differentiable manifold M]
        """
        if self._dest_map.is_identity():
            return self._base_space.frames()
        else:
            # Filter out all frames with respect to dest_map:
            frames = []
            for frame in self._base_space.frames():
                if frame._dest_map == self._dest_map:
                    frames.append(frame)
            return frames

    def coframes(self):
        r"""
        Return the list of coframes defined on the base manifold of ``self``
        with respect to the destination map.

        .. SEEALSO::

            For further details on frames on ``self`` see
            :meth:`local_frame`.

        OUTPUT: list of coframes defined on ``self``

        EXAMPLES:

        Coframes on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: TM = M.tangent_bundle()
            sage: TM.coframes()
            [Coordinate coframe (R^2, (dx,dy))]
            sage: e = TM.vector_frame('e')
            sage: M.coframes()
            [Coordinate coframe (R^2, (dx,dy)), Coframe (R^2, (e^0,e^1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1})
            sage: TU = U.tangent_bundle()
            sage: TU.coframes()
            [Coordinate coframe (U, (dx,dy))]
            sage: e.restrict(U)
            Vector frame (U, (e_0,e_1))
            sage: TU.coframes()
            [Coordinate coframe (U, (dx,dy)), Coframe (U, (e^0,e^1))]
            sage: TM.coframes()
            [Coordinate coframe (R^2, (dx,dy)),
             Coframe (R^2, (e^0,e^1)),
             Coordinate coframe (U, (dx,dy)),
             Coframe (U, (e^0,e^1))]
        """
        if self._dest_map.is_identity():
            return self._base_space.coframes()
        else:
            # Filter out all coframes with respect to dest_map:
            coframes = []
            for coframe in self._base_space.coframes():
                if coframe._dest_map == self._dest_map:
                    coframes.append(coframe)
            return coframes

    def trivialization(self, coordinates='', names=None, calc_method=None):
        r"""
        Return a trivialization of ``self`` in terms of a chart on the codomain
        of the destination map.

        .. NOTE::

            Since a chart gives direct rise to a trivialization, this method is
            nothing but an invocation of
            :meth:`~sage.manifolds.manifold.TopologicalManifold.chart` of the
            class
            :class:`~sage.manifolds.manifold.TopologicalManifold`.


        INPUT:

        - ``coordinates`` -- (default: ``''`` (empty string)) string
          defining the coordinate symbols, ranges and possible periodicities,
          see below
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)
        - ``calc_method`` -- (default: ``None``) string defining the calculus
          method to be used on this chart; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the current calculus method defined on the manifold is
            used (cf.
            :meth:`~sage.manifolds.manifold.TopologicalManifold.set_calculus_method`)

        The coordinates declared in the string ``coordinates`` are
        separated by ``' '`` (whitespace) and each coordinate has at most four
        fields, separated by a colon (``':'``):

        1. The coordinate symbol (a letter or a few letters).
        2. (optional, only for manifolds over `\RR`) The interval `I`
           defining the coordinate range: if not provided, the coordinate
           is assumed to span all `\RR`; otherwise `I` must be provided
           in the form ``(a,b)`` (or equivalently ``]a,b[``)
           The bounds ``a`` and ``b`` can be ``+/-Infinity``, ``Inf``,
           ``infinity``, ``inf`` or ``oo``. For *singular* coordinates,
           non-open intervals such as ``[a,b]`` and
           ``(a,b]`` (or equivalently ``]a,b]``) are allowed. Note that
           the interval declaration must not contain any space character.
        3. (optional) Indicator of the periodic character of the coordinate,
           either as ``period=T``, where ``T`` is the period, or, for manifolds
           over `\RR` only, as the keyword ``periodic`` (the value of the
           period is then deduced from the interval `I` declared in field 2;
           see the example below)
        4. (optional) The LaTeX spelling of the coordinate; if not provided
           the coordinate symbol given in the first field will be used.

        The order of fields 2 to 4 does not matter and each of them can
        be omitted. If it contains any LaTeX expression, the string
        ``coordinates`` must be declared with the prefix 'r' (for "raw") to
        allow for a proper treatment of the backslash character (see
        examples below). If no interval range, no period and no LaTeX spelling
        is to be set for any coordinate, the argument ``coordinates`` can be
        omitted when the shortcut operator ``<,>`` is used to declare the
        trivialization.

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.manifolds.chart.Chart` or one of its subclasses, like
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart` for
          differentiable manifolds over `\RR`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: TM = M.tangent_bundle()
            sage: X = TM.trivialization('x y'); X
            Chart (M, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)
        """
        return self._ambient_domain.chart(coordinates=coordinates, names=names,
                                          calc_method=calc_method)

    def transitions(self):
        r"""
        Return the transition maps between trivialization maps in terms of
        coordinate changes defined via charts on the codomain of the
        destination map.

        .. NOTE::

            Since a chart gives direct rise to a trivialization, this method is
            nothing but an invocation of
            :meth:`~sage.manifolds.manifold.TopologicalManifold.coord_changes`
            of the class
            :class:`~sage.manifolds.manifold.TopologicalManifold`.

        EXAMPLES:

        Various changes of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
            sage: TM = M.tangent_bundle()
            sage: TM.transitions()
            {(Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y))
               to Chart (M, (u, v))}
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: TM.transitions()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v))
              to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y))
               to Chart (M, (u, v))}
            sage: c_rs.<r,s> = M.chart()
            sage: uv_to_rs = c_uv.transition_map(c_rs, [-u+2*v, 3*u-v])
            sage: TM.transitions()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v))
               to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v))
              to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y))
               to Chart (M, (u, v))}
            sage: xy_to_rs = uv_to_rs * xy_to_uv
            sage: TM.transitions()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v))
               to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v))
              to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y))
               to Chart (M, (u, v)),
             (Chart (M, (x, y)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (x, y))
               to Chart (M, (r, s))}
        """
        return self._ambient_domain.coord_changes()

    def transition(self, chart1, chart2):
        r"""
        Return the change of trivializations in terms of a coordinate change
        between two differentiable charts defined on the codomain of the
        destination map.

        The differentiable chart must have been defined previously, for
        instance by the method
        :meth:`~sage.manifolds.chart.Chart.transition_map`.

        .. NOTE::

            Since a chart gives direct rise to a trivialization, this method is
            nothing but an invocation of
            :meth:`~sage.manifolds.manifold.TopologicalManifold.coord_change`
            of the class :class:`~sage.manifolds.manifold.TopologicalManifold`.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines coord. change
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: TM = M.tangent_bundle()
            sage: TM.transition(c_xy, c_uv) # returns the coord. change above
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        """
        return self._ambient_domain.coord_change(chart1, chart2)

    def is_manifestly_trivial(self):
        r"""
        Return ``True`` if ``self`` is known to be a trivial and ``False``
        otherwise.

        If ``False`` is returned, either the tensor bundle is not trivial
        or no vector frame has been defined on it yet.

        EXAMPLES:

        A just created manifold has a priori no manifestly trivial tangent
        bundle::

            sage: M = Manifold(2, 'M')
            sage: TM = M.tangent_bundle()
            sage: TM.is_manifestly_trivial()
            False

        Defining a vector frame on it makes it trivial::

            sage: e = TM.vector_frame('e')
            sage: TM.is_manifestly_trivial()
            True

        Defining a coordinate chart on the whole manifold also makes it
        trivial::

            sage: N = Manifold(4, 'N')
            sage: X.<t,x,y,z> = N.chart()
            sage: TN = N.tangent_bundle()
            sage: TN.is_manifestly_trivial()
            True

        The situation is not so clear anymore when a destination map to a
        non-parallelizable manifold is stated::

            sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereo coord from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereo coord from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2),
            ....:                                       y/(x^2+y^2)),
            ....:                                intersection_name='W',
            ....:                                restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: Phi = U.diff_map(M, {(c_xy, c_xy): [x, y]},
            ....:                  name='Phi') # inclusion map
            sage: PhiTU = U.tangent_bundle(dest_map=Phi); PhiTU
            Tangent bundle Phi^*TS^2 over the Open subset U of the
             2-dimensional differentiable manifold S^2 along the
             Differentiable map Phi from the Open subset U of the
             2-dimensional differentiable manifold S^2 to the 2-dimensional
             differentiable manifold S^2

        A priori, the pullback tangent bundle is not trivial::

            sage: PhiTU.is_manifestly_trivial()
            False

        But certainly, this bundle must be trivial since `U` is parallelizable.
        To ensure this, we need to define a local frame on `U` with values
        in `\Phi^*TS^2`::

            sage: PhiTU.local_frame('e', from_frame=c_xy.frame())
            Vector frame (U, (e_0,e_1)) with values on the 2-dimensional
             differentiable manifold S^2
            sage: PhiTU.is_manifestly_trivial()
            True
        """
        if self._dest_map.is_identity():
            # The standard case:
            return self._base_space.is_manifestly_parallelizable()
        else:
            # If the ambient domain is manifestly trivial, the pullback bundle
            # is certainly trivial:
            if self._ambient_domain.is_manifestly_parallelizable():
                return True
            # Otherwise check whether a global frame on the pullback bundle is
            # defined:
            for frame in self.frames():
                if frame._domain is self._base_space:
                    return True
            return False

    def local_frame(self, *args, **kwargs):
        r"""
        Define a vector frame on ``domain``, possibly with values in the
        tangent bundle of the ambient domain.

        If the basis specified by the given symbol already exists, it is
        simply returned.
        If no argument is provided the vector field module's default frame is
        returned.

        Notice, that a vector frame automatically induces a local frame on the
        tensor bundle ``self``. More precisely, if `e: U \to \Phi^*TN` is a
        vector frame on `U \subset M` with values in `\Phi^*TN` along the
        destination map

        .. MATH::

            \Phi: M \longrightarrow N

        then the map

        .. MATH::

            p \mapsto \Big(\underbrace{e^*(p), \dots, e^*(p)}_{k\ \; \text{times}},
            \underbrace{e(p), \dots, e(p)}_{l\ \; \text{times}}\Big) \in
            T^{(k,l)}_q N ,

        with `q=\Phi(p)`, defines a basis at each point `p \in U` and
        therefore gives rise to a local frame on `\Phi^* T^{(k,l)}N` on the
        domain `U`.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
            for complete documentation.

        INPUT:

        - ``symbol`` -- (default: ``None``) either a string, to be used as a
          common base for the symbols of the vector fields constituting the
          vector frame, or a list/tuple of strings, representing the individual
          symbols of the vector fields; can be ``None`` only if ``from_frame``
          is not ``None`` (see below)
        - ``vector_fields`` -- tuple or list of `n` linearly independent vector
          fields on ``domain`` (`n` being the dimension of  ``domain``)
          defining the vector frame; can be omitted if the vector frame is
          created from scratch or if ``from_frame`` is not ``None``
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the vector fields
          constituting the vector frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the vector fields;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde{e}`
          on the codomain `N` of the destination map `\Phi`; the returned
          frame `e` is then such that for all `p \in U`,
          we have `e(p) = \tilde{e}(\Phi(p))`
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the vector fields of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the vector fields;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe
        - ``domain`` -- (default: ``None``) domain on which the local frame is
          defined; if ``None`` is provided, the base space of ``self`` is
          assumed

        OUTPUT:

        - the vector frame corresponding to the above specifications; this is
          an instance of
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`.

        EXAMPLES:

        Defining a local frame for the tangent bundle of a 3-dimensional
        manifold::

            sage: M = Manifold(3, 'M')
            sage: TM = M.tangent_bundle()
            sage: e = TM.local_frame('e'); e
            Vector frame (M, (e_0,e_1,e_2))
            sage: e[0]
            Vector field e_0 on the 3-dimensional differentiable manifold M

        Specifying the domain of the vector frame::

            sage: U = M.open_subset('U')
            sage: f = TM.local_frame('f', domain=U); f
            Vector frame (U, (f_0,f_1,f_2))
            sage: f[0]
            Vector field f_0 on the Open subset U of the 3-dimensional
             differentiable manifold M

        .. SEEALSO::

            For more options, in particular for the choice of symbols and
            indices, see
            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`.
        """
        domain = kwargs.pop('domain', None)
        if domain is None:
            domain = self._base_space
        dest_map = self._dest_map.restrict(domain)
        if not args and not kwargs:
            # if no argument is provided, the default basis of the
            # base vector field module is returned:
            return domain.vector_field_module(dest_map=dest_map,
                                              force_free=True).basis()
        kwargs['dest_map'] = dest_map
        return domain.vector_frame(*args, **kwargs)

    vector_frame = local_frame

    def ambient_domain(self):
        r"""
        Return the codomain of the destination map.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing the codomain of the destination map

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: e_cart = c_cart.frame() # standard basis
            sage: R = Manifold(1, 'R')
            sage: T.<t> = R.chart()  # canonical chart on R
            sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold R to the 2-dimensional differentiable manifold M
            sage: Phi.display()
            Phi: R → M
               t ↦ (x, y) = (cos(t), sin(t))
            sage: PhiT11 = R.tensor_bundle(1, 1, dest_map=Phi)
            sage: PhiT11.ambient_domain()
            2-dimensional differentiable manifold M
        """
        return self._ambient_domain

    def destination_map(self):
        r"""
        Return the destination map.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.diff_map.DifferentialMap`
          representing the destination map

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: e_cart = c_cart.frame() # standard basis
            sage: R = Manifold(1, 'R')
            sage: T.<t> = R.chart()  # canonical chart on R
            sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold R to the 2-dimensional differentiable manifold M
            sage: Phi.display()
            Phi: R → M
               t ↦ (x, y) = (cos(t), sin(t))
            sage: PhiT11 = R.tensor_bundle(1, 1, dest_map=Phi)
            sage: PhiT11.destination_map()
            Differentiable map Phi from the 1-dimensional differentiable
             manifold R to the 2-dimensional differentiable manifold M
        """
        return self._dest_map

    def default_frame(self):
        r"""
        Return the default vector frame defined on ``self``.

        By *vector frame*, it is meant a field on the manifold that provides,
        at each point `p`, a vector basis of the pulled back tangent space at
        `p`.

        If the destination map is the identity map, the default frame is the
        the first one defined on the manifold, usually the coordinate frame,
        unless it is changed via :meth:`set_default_frame`.

        If the destination map is non-trivial, the default frame usually must
        be set via :meth:`set_default_frame`.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          representing the default vector frame

        EXAMPLES:

        The default vector frame is often the coordinate frame associated
        with the first chart defined on the manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: TM = M.tangent_bundle()
            sage: TM.default_frame()
            Coordinate frame (M, (∂/∂x,∂/∂y))
        """
        def_bframe = self._base_space.default_frame()
        if self._def_frame is None and def_bframe is not None:
            if def_bframe._dest_map == self._dest_map:
                self._def_frame = def_bframe
        return self._def_frame

    def set_default_frame(self, frame):
        r"""
        Changing the default vector frame on ``self``.

        .. NOTE::

            If the destination map is the identity, the default frame of the
            base manifold gets changed here as well.

        INPUT:

        - ``frame`` --
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          a vector frame defined on the base manifold

        EXAMPLES:

        Changing the default frame on the tangent bundle of a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: TM = M.tangent_bundle()
            sage: e = TM.vector_frame('e')
            sage: TM.default_frame()
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: TM.set_default_frame(e)
            sage: TM.default_frame()
            Vector frame (M, (e_0,e_1))
            sage: M.default_frame()
            Vector frame (M, (e_0,e_1))
        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        if not isinstance(frame, VectorFrame):
            raise TypeError("{} is not a vector frame".format(frame))
        if (not frame._domain.is_subset(self._base_space) or
                frame._dest_map != self._dest_map):
            raise ValueError("the frame must be defined on " +
                             "the {}".format(self))
        if self._dest_map.is_identity():
            self._base_space.set_default_frame(frame)
        else:
            frame._fmodule.set_default_basis(frame)
        self._def_frame = frame

    def set_orientation(self, orientation):
        r"""
        Set the preferred orientation of ``self``.

        INPUT:

        - ``orientation`` -- a vector frame or a list of vector frames, covering
          the base space of ``self``

        .. NOTE::

            If the destination map is the identity, the preferred orientation
            of the base manifold gets changed here as well.

        .. WARNING::

            It is the user's responsibility that the orientation set here
            is indeed an orientation. There is no check going on in the
            background. See :meth:`orientation` for the definition of an
            orientation.

        EXAMPLES:

        Set an orientation on a tensor bundle::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: T11 = M.tensor_bundle(1, 1)
            sage: e = T11.local_frame('e'); e
            Vector frame (M, (e_0,e_1))
            sage: T11.set_orientation(e)
            sage: T11.orientation()
            [Vector frame (M, (e_0,e_1))]

        Set an orientation in the non-trivial case::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: T12 = M.tensor_bundle(1, 2)
            sage: e = T12.local_frame('e', domain=U)
            sage: f = T12.local_frame('f', domain=V)
            sage: T12.set_orientation([e, f])
            sage: T12.orientation()
            [Vector frame (U, (e_0,e_1)), Vector frame (V, (f_0,f_1))]
        """
        if self._dest_map.is_identity():
            base_space = self._base_space
            base_space.set_orientation(orientation)
            self._orientation = base_space._orientation
        else:
            super().set_orientation(orientation)

    def orientation(self):
        r"""
        Get the preferred orientation of ``self`` if available.

        See :meth:`~sage.manifolds.vector_bundle.TopologicalVectorBundle.orientation`
        for details regarding orientations on vector bundles.

        The tensor bundle `\Phi^* T^{(k,l)}N` of a manifold is orientable if
        the manifold `\Phi(M)` is orientable. The converse does not
        necessarily hold true. The usual case corresponds to `\Phi`
        being the identity map, where the tensor bundle `T^{(k,l)}M` is
        orientable if and only if the manifold `M` is orientable.

        .. NOTE::

            Notice that the orientation of a general tensor bundle
            `\Phi^* T^{(k,l)}N` is canonically induced by the orientation of
            the tensor bundle `\Phi^* T^{(1,0)}N` as each local frame there
            induces the frames on `\Phi^* T^{(k,l)}N` in a canonical way.

        If no preferred orientation has been set before, and if the ambient
        space already admits a preferred orientation, the corresponding
        orientation is returned and henceforth fixed for the tensor bundle.

        EXAMPLES:

        In the trivial case, i.e. if the destination map is the identity
        and the tangent bundle is covered by one frame, the orientation is
        easily obtained::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: T11 = M.tensor_bundle(1, 1)
            sage: T11.orientation()
            [Coordinate frame (M, (∂/∂x,∂/∂y))]

        The same holds true if the ambient domain admits a trivial
        orientation::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: R = Manifold(1, 'R')
            sage: c_t.<t> = R.chart()
            sage: Phi = R.diff_map(M, name='Phi')
            sage: PhiT22 = R.tensor_bundle(2, 2, dest_map=Phi); PhiT22
            Tensor bundle Phi^*T^(2,2)M over the 1-dimensional differentiable
             manifold R along the Differentiable map Phi from the 1-dimensional
             differentiable manifold R to the 2-dimensional differentiable
             manifold M
            sage: PhiT22.local_frame()  # initialize frame
            Vector frame (R, (∂/∂x,∂/∂y)) with values on the 2-dimensional
             differentiable manifold M
            sage: PhiT22.orientation()
            [Vector frame (R, (∂/∂x,∂/∂y)) with values on the 2-dimensional
             differentiable manifold M]
            sage: PhiT22.local_frame() is PhiT22.orientation()[0]
            True

        In the non-trivial case, however, the orientation must be set
        manually by the user::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: T11 = M.tensor_bundle(1, 1); T11
            Tensor bundle T^(1,1)M over the 2-dimensional differentiable
             manifold M
            sage: T11.orientation()
            []
            sage: T11.set_orientation([c_xy.frame(), c_uv.frame()])
            sage: T11.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]

        If the destination map is the identity, the orientation is
        automatically set for the manifold, too::

            sage: M.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]

        Conversely, if one sets an orientation on the manifold,
        the orientation on its tensor bundles is set accordingly::

            sage: c_tz.<t,z> = U.chart()
            sage: M.set_orientation([c_tz, c_uv])
            sage: T11.orientation()
            [Coordinate frame (U, (∂/∂t,∂/∂z)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]
        """
        if self._dest_map.is_identity():
            self._orientation = self._base_space.orientation()
        if not self._orientation:
            if not self._dest_map.is_identity():
                # try to get orientation from ambient space:
                ambient_domain = self._ambient_domain
                amb_orient = ambient_domain.orientation()
                if amb_orient:
                    for frame in self.frames():
                        from_frame = frame._from_frame
                        if from_frame in amb_orient:
                            self._orientation.append(frame)
        return list(self._orientation)
