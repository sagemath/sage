# sage_setup: distribution = sagemath-symbolics
r"""
Section Modules

The set of sections over a vector bundle `E \to M` of class `C^k` on a domain
`U \in M` is a module over the algebra `C^k(U)` of scalar fields on `U`.

Depending on the domain, there are two classes of section modules:

- :class:`SectionModule` for local sections over a non-trivial part of a
  topological vector bundle
- :class:`SectionFreeModule` for local sections over a trivial part of a
  topological vector bundle

AUTHORS:

- Michael Jung (2019): initial version
"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.rings.infinity import infinity
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method
from sage.categories.modules import Modules
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.section import Section, TrivialSection

class SectionModule(UniqueRepresentation, Parent):
    r"""
    Module of sections over a vector bundle `E \to M` of class `C^k` on a domain
    `U \in M`.

    The *section module* `C^k(U;E)` is the set of all `C^k`-maps, called
    *sections*, of type

    .. MATH::

        s: U \longrightarrow E

    such that

    .. MATH::

        \forall p \in U,\ s(p) \in E_p,

    where `E_p` is the vector bundle fiber of `E` at the point `p`.

    `C^k(U;E)` is a module over `C^k(U)`, the algebra of `C^k` scalar fields on
    `U`.

    INPUT:

    - ``vbundle`` -- vector bundle `E` on which the sections takes its values
    - ``domain`` -- (default: ``None``) subdomain `U` of the base space on which
      the sections are defined

    EXAMPLES:

    Module of sections on the Möbius bundle::

        sage: M = Manifold(1, 'RP^1', structure='top', start_index=1)
        sage: U = M.open_subset('U')  # the complement of one point
        sage: c_u.<u> =  U.chart() # [1:u] in homogeneous coord.
        sage: V = M.open_subset('V') # the complement of the point u=0
        sage: M.declare_union(U,V)   # [v:1] in homogeneous coord.
        sage: c_v.<v> = V.chart()
        sage: u_to_v = c_u.transition_map(c_v, (1/u),
        ....:                             intersection_name='W',
        ....:                             restrictions1 = u!=0,
        ....:                             restrictions2 = v!=0)
        sage: v_to_u = u_to_v.inverse()
        sage: W = U.intersection(V)
        sage: E = M.vector_bundle(1, 'E')
        sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
        ....:                          domain=U)
        sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
        ....:                          domain=V)
        sage: transf = phi_U.transition_map(phi_V, [[u]])
        sage: C0 = E.section_module(); C0
        Module C^0(RP^1;E) of sections on the 1-dimensional topological manifold
         RP^1 with values in the real vector bundle E of rank 1

    `C^0(\RR P^1;E)` is a module over the algebra `C^0(\RR P^1)`::

        sage: C0.category()
        Category of modules over Algebra of scalar fields on the 1-dimensional
         topological manifold RP^1
        sage: C0.base_ring() is M.scalar_field_algebra()
        True

    However, `C^0(\RR P^1;E)` is not a free module::

        sage: isinstance(C0, FiniteRankFreeModule)
        False

    since the Möbius bundle is not trivial::

        sage: E.is_manifestly_trivial()
        False

    The section module over `U`, on the other hand, is a free module since
    `E|_U` admits a trivialization and therefore has a local frame::

        sage: C0_U = E.section_module(domain=U)
        sage: isinstance(C0_U, FiniteRankFreeModule)
        True

    The zero element of the module::

        sage: z = C0.zero() ; z
        Section zero on the 1-dimensional topological manifold RP^1 with values
         in the real vector bundle E of rank 1
        sage: z.display(phi_U.frame())
        zero = 0
        sage: z.display(phi_V.frame())
        zero = 0

    The module `C^0(M;E)` coerces to any module of sections defined
    on a subdomain of `M`, for instance `C^0(U;E)`::

        sage: C0_U.has_coerce_map_from(C0)
        True
        sage: C0_U.coerce_map_from(C0)
        Coercion map:
          From: Module C^0(RP^1;E) of sections on the 1-dimensional topological
           manifold RP^1 with values in the real vector bundle E of rank 1
          To:   Free module C^0(U;E) of sections on the Open subset U of the
           1-dimensional topological manifold RP^1 with values in the real vector
           bundle E of rank 1

    The conversion map is actually the restriction of sections defined
    on `M` to `U`.
    """
    Element = Section

    def __init__(self, vbundle, domain):
        r"""
        Construct the module of continuous sections over a vector bundle.

        TESTS::

            sage: M = Manifold(1, 'S^1', latex_name=r'S^1', start_index=1,
            ....:               structure='topological')
            sage: U = M.open_subset('U')
            sage: c_x.<x> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_u.<u> = V.chart()
            sage: M.declare_union(U, V)
            sage: x_to_u = c_x.transition_map(c_u, 1/x, intersection_name='W',
            ....:                   restrictions1= x!=0, restrictions2= u!=0)
            sage: W = U.intersection(V)
            sage: u_to_x = x_to_u.inverse()
            sage: E = M.vector_bundle(1, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V)
            sage: transf = phi_U.transition_map(phi_V, [[-1]])
            sage: C0 = E.section_module(); C0
            Module C^0(S^1;E) of sections on the 1-dimensional topological
             manifold S^1 with values in the real vector bundle E of rank 1
            sage: TestSuite(C0).run()
        """
        base_space = vbundle.base_space()
        if not domain.is_subset(base_space):
            raise ValueError("domain must be a subset of base space")
        if vbundle._diff_degree == infinity:
            repr_deg = "infinity" # to skip the "+" in repr(infinity)
            latex_deg = r"\infty" # to skip the "+" in latex(infinity)
        else:
            repr_deg = r"{}".format(vbundle._diff_degree)
            latex_deg = r"{}".format(vbundle._diff_degree)
        self._name = "C^{}({};{})".format(repr_deg, domain._name, vbundle._name)
        self._latex_name = r"C^{" + latex_deg + r"}" + \
                           r"({};{})".format(domain._latex_name,
                                             vbundle._latex_name)
        self._vbundle = vbundle
        self._domain = domain
        self._base_space = vbundle.base_space()
        self._ring = domain.scalar_field_algebra()
        self._def_frame = None
        Parent.__init__(self, base=self._ring,
                        category=Modules(self._ring))

    #### Begin of parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of the module.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xyz_U = c_xyz.restrict(U); c_xyz_V = c_xyz.restrict(V)
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: e = E.local_frame('e', domain=U)
            sage: f = E.local_frame('f', domain=V)
            sage: s = C0([-x,y], frame=e, name='s'); s
            Section s on the 3-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
            sage: s.display(e)
            s = -x e_0 + y e_1
            sage: C0(0) is C0.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, Section):
            if self._domain.is_subset(comp._domain):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {} ".format(comp) +
                                 "to a local section in {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unnamed) element of the module.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: C0._an_element_()
            Section on the 2-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
        """
        resu = self.element_class(self)
        for oc in self._domain.open_covers(trivial=False):
            # the first non-trivial open cover is selected
            for dom in oc:
                smodule_dom = self._vbundle.section_module(domain=dom)
                resu.set_restriction(smodule_dom._an_element_())
            return resu
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from ``other`` parent.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: C0_U = E.section_module(domain=U)
            sage: C0._coerce_map_from_(C0_U)
            False
            sage: C0_U._coerce_map_from_(C0)
            True
        """
        if isinstance(other, (SectionModule, SectionFreeModule)):
            return self._domain.is_subset(other._domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: C0._repr_()
            'Module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2'
            sage: repr(C0)  # indirect doctest
            'Module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2'
            sage: C0  # indirect doctest
            Module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
        """
        desc = "Module {} of sections on the {} with values in the {} vector " \
               "bundle {} of rank {}"
        desc = desc.format(self._name, self._domain,
                           self._vbundle.base_field_type(),
                           self._vbundle._name,
                           self._vbundle.rank())
        return desc

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: C = E.section_module()
            sage: C._latex_()
            'C^{\\infty}(M;E)'
            sage: latex(C) # indirect doctest
            C^{\infty}(M;E)
        """
        return self._latex_name

    def base_space(self):
        r"""
        Return the base space of the sections in this module.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = U.vector_bundle(2, 'E')
            sage: C0 = E.section_module(); C0
            Module C^0(U;E) of sections on the Open subset U of the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: C0.base_space()
            Open subset U of the 3-dimensional topological manifold M
        """
        return self._base_space

    def domain(self):
        r"""
        Return the domain of the section module.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0_U = E.section_module(domain=U); C0_U
            Module C^0(U;E) of sections on the Open subset U of the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: C0_U.domain()
            Open subset U of the 3-dimensional topological manifold M
        """
        return self._domain

    def vector_bundle(self):
        r"""
        Return the overlying vector bundle on which the section module is
        defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(); C0
            Module C^0(M;E) of sections on the 3-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
            sage: C0.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
             3-dimensional topological manifold M
            sage: E is C0.vector_bundle()
            True
        """
        return self._vbundle

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: z = C0.zero(); z
            Section zero on the 2-dimensional topological manifold M with values
             in the real vector bundle E of rank 2
            sage: z == 0
            True
        """
        res = self.element_class(self, name='zero', latex_name='0')
        for frame in self._vbundle._frames:
            if frame._domain.is_subset(self._domain):
                res.add_comp(frame)
                # (since new components are initialized to zero)
        res.set_immutable()
        return res

    def default_frame(self):
        r"""
        Return the default frame defined on ``self``.

        EXAMPLES:

        Get the default local frame of a non-trivial section module::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: e = E.local_frame('e', domain=U)
            sage: C0.default_frame()
            Local frame (E|_U, (e_0,e_1))

        The local frame is indeed the same, and not a copy::

            sage: e is C0.default_frame()
            True
        """
        return self._def_frame

    def set_default_frame(self, basis):
        r"""
        Set the default local frame on ``self``.

        EXAMPLES:

        Set a default frame of a non-trivial section module::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(); C0
            Module C^0(M;E) of sections on the 3-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
            sage: e = E.local_frame('e', domain=U)
            sage: C0.set_default_frame(e)
            sage: C0.default_frame()
            Local frame (E|_U, (e_0,e_1))

        The local frame is indeed the same, and not a copy::

            sage: e is C0.default_frame()
            True

        Notice, that the local frame is defined on a subset and is not part of
        the section module `C^k(M;E)`::

            sage: C0.default_frame().domain()
            Open subset U of the 3-dimensional topological manifold M
        """
        from .local_frame import LocalFrame
        if not isinstance(basis, LocalFrame):
            raise ValueError("the argument is not a local frame")
        elif not basis._domain.is_subset(self._domain):
            raise ValueError("local frame's domain must be a subset "
                             "of the {}".format(self._domain))
        self._def_frame = basis

#******************************************************************************

class SectionFreeModule(FiniteRankFreeModule):
    r"""
    Free module of sections over a vector bundle `E \to M` of class `C^k` on a
    domain `U \in M` which admits a trivialization or local frame.

    The *section module* `C^k(U;E)` is the set of all `C^k`-maps, called
    *sections*, of type

    .. MATH::

        s: U \longrightarrow E

    such that

    .. MATH::

        \forall p \in U,\ s(p) \in E_p,

    where `E_p` is the vector bundle fiber of `E` at the point `p`.

    Since the domain `U` admits a local frame, the corresponding vector bundle
    `E|_U \to U` is trivial and `C^k(U;E)` is a free module over `C^k(U)`.

    .. NOTE::

        If `E|_U` is not trivial, the class :class:`SectionModule` should be
        used instead, for `C^k(U;E)` is no longer a free module.

    INPUT:

    - ``vbundle`` -- vector bundle `E` on which the sections takes its values
    - ``domain`` -- (default: ``None``) subdomain `U` of the base space on which
      the sections are defined

    EXAMPLES:

    Module of sections on the 2-rank trivial vector bundle over the Euclidean
    plane `\RR^2`::

        sage: M = Manifold(2, 'R^2', structure='top')
        sage: c_cart.<x,y> = M.chart()
        sage: E = M.vector_bundle(2, 'E')
        sage: e = E.local_frame('e') # Trivializes the vector bundle
        sage: C0 = E.section_module(); C0
        Free module C^0(R^2;E) of sections on the 2-dimensional topological
         manifold R^2 with values in the real vector bundle E of rank 2
        sage: C0.category()
        Category of finite dimensional modules over Algebra of scalar fields on
         the 2-dimensional topological manifold R^2
        sage: C0.base_ring() is M.scalar_field_algebra()
        True

    The vector bundle admits a global frame and is therefore trivial::

        sage: E.is_manifestly_trivial()
        True

    Since the vector bundle is trivial, its section module of global sections
    is a free module::

        sage: isinstance(C0, FiniteRankFreeModule)
        True

    Some elements are::

        sage: C0.an_element().display()
        2 e_0 + 2 e_1
        sage: C0.zero().display()
        zero = 0
        sage: s = C0([-y,x]); s
        Section on the 2-dimensional topological manifold R^2 with values in the
         real vector bundle E of rank 2
        sage: s.display()
        -y e_0 + x e_1

    The rank of the free module equals the rank of the vector bundle::

        sage: C0.rank()
        2

    The basis is given by the definition above::

        sage: C0.bases()
        [Local frame (E|_R^2, (e_0,e_1))]

    The test suite is passed as well::

        sage: TestSuite(C0).run()
    """
    Element = TrivialSection

    def __init__(self, vbundle, domain):
        r"""
        Construct the free module of sections over a trivial part of a vector
        bundle.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: from sage.manifolds.section_module import SectionFreeModule
            sage: C0 = SectionFreeModule(E, M); C0
            Free module C^0(M;E) of sections on the 3-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
            sage: C0 is E.section_module(force_free=True)
            True
            sage: TestSuite(C0).run()
        """
        from .scalarfield import ScalarField
        self._domain = domain
        name = "C^0({};{})".format(domain._name, vbundle._name)
        latex_name = r'C^0({};{})'.format(domain._latex_name,
                                          vbundle._latex_name)
        base_space = vbundle.base_space()
        self._base_space = base_space
        self._vbundle = vbundle
        cat = Modules(domain.scalar_field_algebra()).FiniteDimensional()
        FiniteRankFreeModule.__init__(self, domain.scalar_field_algebra(),
                               vbundle.rank(), name=name,
                               latex_name=latex_name,
                               start_index=base_space._sindex,
                               output_formatter=ScalarField.coord_function,
                               category=cat)

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: C0 = E.section_module(); C0
            Free module C^0(M;E) of sections on the 3-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
            sage: s = C0([-x,y], basis=e, name='s'); s
            Section s on the 3-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
            sage: s.display(e)
            s = -x e_0 + y e_1
            sage: C0(0) is C0.zero()
            True
        """
        try:
            if comp.is_trivial_zero():
                return self.zero()
        except AttributeError:
            if comp == 0:
                return self.zero()
        if isinstance(comp, Section):
            if self._domain.is_subset(comp._domain):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {}".format(comp) +
                                 "to a local section in {}".format(self))
        if not isinstance(comp, (list, tuple)):
            raise TypeError("cannot convert the {} ".format(comp) +
                            "to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass FiniteRankFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from parent ``other``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: C0_U = E.section_module(domain=U, force_free=True)
            sage: C0._coerce_map_from_(C0_U)
            False
            sage: C0_U._coerce_map_from_(C0)
            True
        """
        if isinstance(other, (SectionModule, SectionFreeModule)):
            return self._domain.is_subset(other._domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: C0._repr_()
            'Free module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2'
            sage: repr(C0)  # indirect doctest
            'Free module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2'
            sage: C0  # indirect doctest
            Free module C^0(M;E) of sections on the 2-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
        """
        desc = "Free module {} of sections on the {} with values in the {} " \
               "vector bundle {} of rank {}"
        desc = desc.format(self._name, self._domain,
                           self._vbundle.base_field_type(),
                           self._vbundle._name,
                           self._vbundle.rank())
        return desc

    def domain(self):
        r"""
        Return the domain of the section module.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0_U = E.section_module(domain=U, force_free=True); C0_U
            Free module C^0(U;E) of sections on the Open subset U of the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: C0_U.domain()
            Open subset U of the 3-dimensional topological manifold M
        """
        return self._domain

    def base_space(self):
        r"""
        Return the base space of the sections in this module.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = U.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True); C0
            Free module C^0(U;E) of sections on the Open subset U of the
             3-dimensional topological manifold M with values in the real
             vector bundle E of rank 2
            sage: C0.base_space()
            Open subset U of the 3-dimensional topological manifold M
        """
        return self._base_space

    def vector_bundle(self):
        r"""
        Return the overlying vector bundle on which the section module is
        defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True); C0
            Free module C^0(M;E) of sections on the 3-dimensional topological
             manifold M with values in the real vector bundle E of rank 2
            sage: C0.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
             3-dimensional topological manifold M
            sage: E is C0.vector_bundle()
            True
        """
        return self._vbundle

    def basis(self, symbol=None, latex_symbol=None, from_frame=None,
              indices=None, latex_indices=None, symbol_dual=None,
              latex_symbol_dual=None):
        r"""
        Define a basis of ``self``.

        A basis of the section module is actually a local frame on the
        differentiable manifold `U` over which the section module is defined.

        If the basis specified by the given symbol already exists, it is
        simply returned.
        If no argument is provided the module's default basis is returned.

        INPUT:

        - ``symbol`` -- (default: ``None``) either a string, to be used as a
          common base for the symbols of the elements of the basis, or a
          tuple of strings, representing the individual symbols of the
          elements of the basis
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the elements of the basis,
          or a tuple of strings, representing the individual LaTeX symbols
          of the elements of the basis; if ``None``, ``symbol`` is used in
          place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices
          labelling the elements of the basis; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the elements of
          the basis; if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual basis; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual basis
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual basis

        OUTPUT:

        - a :class:`~sage.manifolds.local_frame.LocalFrame` representing a basis
          on ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: e = C0.basis('e'); e
            Local frame (E|_M, (e_0,e_1))

        See :class:`~sage.manifolds.local_frame.LocalFrame` for more examples
        and documentation.
        """
        from sage.manifolds.local_frame import LocalFrame
        if symbol is None:
            symbol = from_frame._symbol
            latex_symbol = from_frame._latex_symbol
            indices = from_frame._indices
            latex_indices = from_frame._latex_indices
            symbol_dual = from_frame._symbol_dual
            latex_symbol_dual = from_frame._latex_symbol_dual
        for other in self._known_bases:
            if symbol == other._symbol:
                return other
        return LocalFrame(self, symbol, latex_symbol=latex_symbol,
                          indices=indices,
                          latex_indices=latex_indices,
                          symbol_dual=symbol_dual,
                          latex_symbol_dual=latex_symbol_dual)

    set_default_frame = FiniteRankFreeModule.set_default_basis
    default_frame = FiniteRankFreeModule.default_basis
