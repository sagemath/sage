# sage_setup: distribution = sagemath-symbolics
r"""
Topological Vector Bundle

Let `K` be a topological field. A *vector bundle* of rank `n` over the field
`K` and over a topological manifold `B` (base space) is a topological manifold
`E` (total space) together with a continuous and surjective map `\pi: E \to B`
such that for every point `p \in B`, we have:

- the set `E_p=\pi^{-1}(p)` has the vector space structure of `K^n`,
- there is a neighborhood `U \subset B` of `p` and a homeomorphism
  (trivialization) `\varphi: \pi^{-1}(p) \to U \times K^n` such that `\varphi`
  is compatible with the fibers, namely `\pi \circ \varphi^{-1} = \mathrm{pr}_1`,
  and `v \mapsto \varphi^{-1}(q,v)` is a linear isomorphism between `K^n` and
  `E_q` for any `q \in U`.

AUTHORS:

- Michael Jung (2019) : initial version

REFERENCES:

- [Lee2013]_
- [Mil1974]_
"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.category_object import CategoryObject
from sage.categories.vector_bundles import VectorBundles
from sage.structure.unique_representation import UniqueRepresentation
import sage.rings.abc
from sage.rings.cc import CC
from sage.rings.real_mpfr import RR
from sage.rings.integer import Integer
from sage.manifolds.vector_bundle_fiber import VectorBundleFiber

class TopologicalVectorBundle(CategoryObject, UniqueRepresentation):
    r"""
    An instance of this class is a topological vector bundle `E \to B` over a
    topological field `K`.

    INPUT:

    - ``rank`` -- positive integer; rank of the vector bundle
    - ``name`` -- string representation given to the total space
    - ``base_space`` -- the base space (topological manifold) over which the
      vector bundle is defined
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
      ``None``, ``VectorBundles(base_space, c_field)`` is assumed (see the
      category :class:`~sage.categories.vector_bundles.VectorBundles`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior would return the previously constructed object corresponding to
      these arguments)

    EXAMPLES:

    A real line bundle over some 4-dimensional topological manifold::

        sage: M = Manifold(4, 'M', structure='top')
        sage: E = M.vector_bundle(1, 'E'); E
        Topological real vector bundle E -> M of rank 1 over the base space
         4-dimensional topological manifold M
        sage: E.base_space()
        4-dimensional topological manifold M
        sage: E.base_ring()
        Real Field with 53 bits of precision
        sage: E.rank()
        1

    For a more sophisticated example, let us define a non-trivial
    2-manifold to work with::

        sage: M = Manifold(2, 'M', structure='top')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:                    intersection_name='W', restrictions1= x>0,
        ....:                    restrictions2= u+v>0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: E = M.vector_bundle(2, 'E'); E
        Topological real vector bundle E -> M of rank 2 over the base space
         2-dimensional topological manifold M

    Now, there a two ways to go. Most effortlessly, we define
    trivializations similar to charts (see
    :class:`~sage.manifolds.trivialization.Trivialization`)::

        sage: phi_U = E.trivialization('phi_U', domain=U); phi_U
        Trivialization (phi_U, E|_U)
        sage: phi_V = E.trivialization('phi_V', domain=V); phi_V
        Trivialization (phi_V, E|_V)
        sage: transf = phi_U.transition_map(phi_V, [[0,x],[x,0]]) # transition map between trivializations
        sage: fU = phi_U.frame(); fU
        Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2)))
        sage: fV = phi_V.frame(); fV
        Trivialization frame (E|_V, ((phi_V^*e_1),(phi_V^*e_2)))
        sage: E.changes_of_frame() # random
        {(Local frame (E|_W, ((phi_U^*e_1),(phi_U^*e_2))),
         Local frame (E|_W, ((phi_V^*e_1),(phi_V^*e_2)))): Automorphism
         phi_U^(-1)*phi_V^(-1) of the Free module C^0(W;E) of sections on
         the Open subset W of the 2-dimensional topological manifold M with
         values in the real vector bundle E of rank 2,
         (Local frame (E|_W, ((phi_V^*e_1),(phi_V^*e_2))),
         Local frame (E|_W, ((phi_U^*e_1),(phi_U^*e_2)))): Automorphism
         phi_U^(-1)*phi_V of the Free module C^0(W;E) of sections on the
         Open subset W of the 2-dimensional topological manifold M with
         values in the real vector bundle E of rank 2}

    Then, the atlas of `E` consists of all known trivializations defined
    on E::

        sage: E.atlas() # a shallow copy of the atlas
        [Trivialization (phi_U, E|_U), Trivialization (phi_V, E|_V)]

    Or we just define frames, an automorphism on the free
    section module over the intersection domain `W` and declare the change
    of frame manually (for more details consult
    :class:`~sage.manifolds.local_frame.LocalFrame`)::

        sage: eU = E.local_frame('eU', domain=U); eU
        Local frame (E|_U, (eU_0,eU_1))
        sage: eUW = eU.restrict(W) # to trivialize E|_W
        sage: eV = E.local_frame('eV', domain=V); eV
        Local frame (E|_V, (eV_0,eV_1))
        sage: eVW = eV.restrict(W)
        sage: a = E.section_module(domain=W).automorphism(); a
        Automorphism of the Free module C^0(W;E) of sections on the Open
         subset W of the 2-dimensional topological manifold M with values in
         the real vector bundle E of rank 2
        sage: a[eUW,:] = [[0,x],[x,0]]
        sage: E.set_change_of_frame(eUW, eVW, a)
        sage: E.change_of_frame(eUW, eVW)
        Automorphism of the Free module C^0(W;E) of sections on the Open
         subset W of the 2-dimensional topological manifold M with values in
         the real vector bundle E of rank 2

    Now, the list of all known frames defined on `E` can be displayed via
    :meth:`frames`::

        sage: E.frames() # a shallow copy of all known frames on E
        [Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2))),
         Trivialization frame (E|_V, ((phi_V^*e_1),(phi_V^*e_2))),
         Local frame (E|_W, ((phi_U^*e_1),(phi_U^*e_2))),
         Local frame (E|_W, ((phi_V^*e_1),(phi_V^*e_2))),
         Local frame (E|_U, (eU_0,eU_1)),
         Local frame (E|_W, (eU_0,eU_1)),
         Local frame (E|_V, (eV_0,eV_1)),
         Local frame (E|_W, (eV_0,eV_1))]

    By definition `E` is a manifold, in this case of dimension 4 (notice
    that the induced charts are not implemented, yet)::

        sage: E.total_space()
        4-dimensional topological manifold E

    The method :meth:`section` returns a section while the method
    :meth:`section_module` returns the section module on the corresponding
    domain::

        sage: s = E.section(name='s'); s
        Section s on the 2-dimensional topological manifold M with values in
         the real vector bundle E of rank 2
        sage: s in E.section_module()
        True
    """
    def __init__(self, rank, name, base_space, field='real',
                 latex_name=None, category=None, unique_tag=None):
        r"""
        Construct a topological vector bundle.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: from sage.manifolds.vector_bundle import TopologicalVectorBundle
            sage: TopologicalVectorBundle(2, 'E', M)
            Topological real vector bundle E -> M of rank 2 over the base space
             2-dimensional topological manifold M
        """
        if base_space is None:
            raise ValueError("a base space must be provided")
        ###
        # Handle the field:
        if field == 'real':
            self._field = RR
            self._field_type = field
        elif field == 'complex':
            self._field = CC
            self._field_type = field
        else:
            self._field = field
            if isinstance(field, sage.rings.abc.RealField):
                self._field_type = 'real'
            elif isinstance(field, sage.rings.abc.ComplexField):
                self._field_type = 'complex'
            else:
                self._field_type = 'neither_real_nor_complex'
        bs_field = base_space.base_field()
        if not bs_field.is_subring(self._field):
            raise ValueError("for concrete implementation, manifold's base "
                             "field must be a subfield of the vector bundle's "
                             "base field")
        ###
        # Get the category:
        if category is None:
            category = VectorBundles(base_space, self._field)
        CategoryObject.__init__(self, base=self._field,
                                category=category)
        # Check rank:
        if not isinstance(rank, (int, Integer)):
            raise TypeError("the rank must be an integer")
        if rank < 1:
            raise ValueError("the rank must be strictly positive")
        ###
        # Define remaining attributes:
        self._rank = rank
        self._diff_degree = 0
        self._base_space = base_space
        self._total_space = None
        self._orientation = []  # set no orientation a priori
        ###
        # Set names:
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        ###
        # Initialize quantities like frames and trivializations:
        self._init_attributes()

    def _init_attributes(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E')
            sage: E._init_attributes()
        """
        self._section_modules = {} # dict of section modules with domains as
                                   # keys
        self._atlas = []  # list of trivializations defined on self
        self._transitions = {} # dictionary of transition maps (key: pair of
                               # of trivializations)
        self._frames = []  # list of local frames for self
        self._frame_changes = {}  # dictionary of changes of frames
        self._coframes = [] # list of local coframes for self
        self._trivial_parts = set() # subsets of base space on which self is
                                    # trivial
        self._def_frame = None

    def base_space(self):
        r"""
        Return the base space of the vector bundle.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.base_space()
            2-dimensional topological manifold M
        """
        return self._base_space

    def base_field_type(self):
        r"""
        Return the type of topological field on which the fibers are defined.

        OUTPUT:

        A string describing the field, with three possible values:

        - ``'real'`` for the real field `\RR`
        - ``'complex'`` for the complex field `\CC`
        - ``'neither_real_nor_complex'`` for a field different from `\RR` and `\CC`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E', field=CC)
            sage: E.base_field_type()
            'complex'
        """
        return self._field_type

    def base_field(self):
        r"""
        Return the field on which the fibers are defined.

        OUTPUT: a topological field

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E', field=CC)
            sage: E.base_field()
            Complex Field with 53 bits of precision
        """
        return self._field

    def rank(self):
        r"""
        Return the rank of the vector bundle.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(3, 'E')
            sage: E.rank()
            3
        """
        return self._rank

    def _repr_object_name(self):
        r"""
        String name of the object without structure.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._repr_object_name()
            'real vector bundle E -> M of rank 1 over the base space
             2-dimensional topological manifold M'
        """
        desc = self.base_field_type() + " "
        desc += "vector bundle "
        desc += self._name + " -> " + self.base_space()._name + " "
        desc += "of rank {} ".format(self._rank)
        desc += "over the base space {}".format(self.base_space())
        return desc

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._repr_()
            'Topological real vector bundle E -> M of rank 1 over the base space
             2-dimensional topological manifold M'
        """
        desc = "Topological "
        desc += self._repr_object_name()
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._latex_()
            'E\\to M'
        """
        latex = self._latex_name
        latex += r'\to '
        latex += self.base_space()._latex_name
        return latex

    def _add_local_frame(self, frame):
        r"""
        Helper method to add local frames to the vector bundle.

        INPUT:

        - ``frame`` -- the local frame that shall be added

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: E._add_local_frame(e)
            sage: E._frames
            [Local frame (E|_M, (e_0,e_1)), Local frame (E|_M, (e_0,e_1))]
        """
        self._trivial_parts.add(frame.domain())
        self._frames.append(frame)

    def trivialization(self, name, domain=None, latex_name=None):
        r"""
        Return a trivialization of ``self`` over the domain ``domain``.

        INPUT:

        - ``domain`` -- (default: ``None``) domain on which the trivialization
          is defined; if ``None`` the base space is assumed
        - ``name`` -- (default: ``None``) name given to the trivialization
        - ``latex_name`` -- (default: ``None``) LaTeX name given to the
          trivialization

        OUTPUT:

        - a :class:`~sage.manifolds.trivialization.Trivialization` representing
          a trivialization of `E`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U); phi
            Trivialization (phi, E|_U)
        """
        if domain is None:
            domain = self._base_space
        from .trivialization import Trivialization
        return Trivialization(self, name, domain=domain, latex_name=latex_name)

    def transitions(self):
        r"""
        Return the transition maps defined over subsets of the base space.

        OUTPUT: dictionary of transition maps, with pairs of trivializations as keys

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: X_UV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_U', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: E.transitions() # random
            {(Trivialization (phi_U, E|_U),
             Trivialization (phi_U, E|_V)): Transition map from Trivialization
             (phi_U, E|_U) to Trivialization (phi_U, E|_V),
             (Trivialization (phi_U, E|_V),
             Trivialization (phi_U, E|_U)): Transition map from Trivialization
             (phi_U, E|_V) to Trivialization (phi_U, E|_U)}
        """
        return self._transitions.copy()

    def transition(self, triv1, triv2):
        r"""
        Return the transition map between two trivializations defined over the
        manifold.

        The transition map must have been defined previously, for instance by
        the method
        :meth:`~sage.manifolds.trivialization.Trivialization.transition_map`.

        INPUT:

        - ``triv1`` -- trivialization 1
        - ``triv2`` -- trivialization 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.trivialization.TransitionMap`
          representing the transition map from trivialization 1 to
          trivialization 2

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: X_UV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: E.transition(phi_V, phi_U)
            Transition map from Trivialization (phi_V, E|_V) to Trivialization
             (phi_U, E|_U)
        """
        if (triv1, triv2) not in self._transitions:
            raise TypeError("the transition map from " +
                            "{} to {}".format(triv1, triv2) + " has not " +
                            "been defined on the {}".format(self))
        return self._transitions[(triv1, triv2)]

    def atlas(self):
        r"""
        Return the list of trivializations that have been defined for ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_M = E.trivialization('phi_M')
            sage: E.atlas()
            [Trivialization (phi_U, E|_U),
             Trivialization (phi_V, E|_V),
             Trivialization (phi_M, E|_M)]
        """
        return list(self._atlas) # Make a (shallow) copy

    def is_manifestly_trivial(self):
        r"""
        Return ``True`` if ``self`` is manifestly a trivial bundle, i.e. there
        exists a frame or a trivialization defined on the whole base space.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: phi_U = E.trivialization('phi_U', domain=U); phi_U
            Trivialization (phi_U, E|_U)
            sage: phi_V = E.trivialization('phi_V', domain=V); phi_V
            Trivialization (phi_V, E|_V)
            sage: E.is_manifestly_trivial()
            False
            sage: E.trivialization('phi_M', M)
            Trivialization (phi_M, E|_M)
            sage: E.is_manifestly_trivial()
            True
        """
        return self.base_space() in self._trivial_parts

    def section_module(self, domain=None, force_free=False):
        r"""
        Return the section module of continuous sections on ``self``.

        See :class:`~sage.manifolds.section_module.SectionModule` for a complete
        documentation.

        INPUT:

        - ``domain`` -- (default: ``None``) the domain on which the module is
          defined; if ``None`` the base space is assumed
        - ``force_free`` -- boolean (default: ``False``); if set to ``True``, force
          the construction of a *free* module (this implies that `E` is trivial)

        OUTPUT:

        - a
          :class:`~sage.manifolds.section_module.SectionModule`
          (or if `E` is trivial, a
          :class:`~sage.manifolds.section_module.SectionFreeModule`)
          representing the module of continuous sections on
          `U` taking values in `E`

        EXAMPLES:

        Module of sections on the Möbius bundle over the real-projective space
        `M=\RR P^1`::

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
            Module C^0(RP^1;E) of sections on the 1-dimensional topological
             manifold RP^1 with values in the real vector bundle E of rank 1

        `C^0(\RR P^1;E)` is a module over the algebra `C^0(\RR P^1)`::

            sage: C0.category()
            Category of modules over Algebra of scalar fields on the
             1-dimensional topological manifold RP^1
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

        The elements of `C^0(U)` are sections on `U`::

            sage: C0_U.an_element()
            Section on the Open subset U of the 1-dimensional topological
             manifold RP^1 with values in the real vector bundle E of rank 1
            sage: C0_U.an_element().display(phi_U.frame())
            2 (phi_U^*e_1)
        """
        if domain is None:
            domain = self._base_space
        from sage.manifolds.section_module import (SectionModule,
                                                   SectionFreeModule)
        if domain not in self._section_modules:
            if force_free or domain in self._trivial_parts:
                self._section_modules[domain] = SectionFreeModule(self, domain)
            else:
                self._section_modules[domain] = SectionModule(self, domain)

        return self._section_modules[domain]

    def fiber(self, point):
        r"""
        Return the vector bundle fiber over a point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` of the base space of ``self``

        OUTPUT:

        - instance of :class:`~sage.manifolds.vector_bundle_fiber.VectorBundleFiber`
          representing the fiber over `p`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,2,1), name='p'); p
            Point p on the 3-dimensional topological manifold M
            sage: E = M.vector_bundle(2, 'E'); E
            Topological real vector bundle E -> M of rank 2 over the base space
             3-dimensional topological manifold M
            sage: E.fiber(p)
            Fiber of E at Point p on the 3-dimensional topological manifold M
        """
        return VectorBundleFiber(self, point)

    def local_frame(self, *args, **kwargs):
        r"""
        Define a local frame on ``self``.

        A *local frame* is a section on a subset `U \subset M` in `E` that
        provides, at each point `p` of the base space, a vector basis of the
        fiber `E_p` at `p`.

        .. SEEALSO::

            :class:`~sage.manifolds.local_frame.LocalFrame` for complete
            documentation.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the sections constituting the local frame, or a list/tuple
          of strings, representing the individual symbols of the sections
        - ``sections`` -- tuple or list of `n` linearly independent sections on
          ``self`` (`n` being the rank of ``self``) defining the local
          frame; can be omitted if the local frame is created from scratch
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe
        - ``domain`` -- (default: ``None``) domain on which the local frame
          is defined; if ``None``, the whole base space is assumed

        OUTPUT:

        - a :class:`~sage.manifolds.local_frame.LocalFrame` representing the
          defined local frame

        EXAMPLES:


        Defining a local frame from two linearly independent sections on a
        real rank-2 vector bundle::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: X.<x,y,z> = U.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U)
            sage: s0 = E.section(name='s_0', domain=U)
            sage: s0[:] = 1+z^2, -2
            sage: s1 = E.section(name='s_1', domain=U)
            sage: s1[:] = 1, 1+x^2
            sage: e = E.local_frame('e', (s0, s1), domain=U); e
            Local frame (E|_U, (e_0,e_1))
            sage: (e[0], e[1]) == (s0, s1)
            True

        If the sections are not linearly independent, an error is raised::

            sage: e = E.local_frame('z', (s0, -s0), domain=U)
            Traceback (most recent call last):
            ...
            ValueError: the provided sections are not linearly independent

        It is also possible to create a local frame from scratch, without
        connecting it to previously defined local frames or sections
        (this can still be performed later via the method
        :meth:`set_change_of_frame`)::

            sage: f = E.local_frame('f', domain=U); f
            Local frame (E|_U, (f_0,f_1))

        For a global frame, the argument ``domain`` is omitted::

            sage: g = E.local_frame('g'); g
            Local frame (E|_M, (g_0,g_1))

        .. SEEALSO::

            For more options, in particular for the choice of symbols and
            indices, see :class:`~sage.manifolds.local_frame.LocalFrame`.
        """
        from sage.manifolds.local_frame import LocalFrame
        # Input processing
        n_args = len(args)
        if n_args < 1 or n_args > 2:
            raise TypeError("local_frame() takes one or two positional "
                            "arguments, not {}".format(n_args))
        symbol = args[0]
        sections = None
        if n_args == 2:
            sections = args[1]
        latex_symbol = kwargs.pop('latex_symbol', None)
        indices = kwargs.pop('indices', None)
        latex_indices = kwargs.pop('latex_indices', None)
        symbol_dual = kwargs.pop('symbol_dual', None)
        latex_symbol_dual = kwargs.pop('latex_symbol_dual', None)
        domain = kwargs.pop('domain', None)

        sec_module = self.section_module(domain=domain, force_free=True)
        resu = LocalFrame(sec_module, symbol=symbol, latex_symbol=latex_symbol,
                          indices=indices, latex_indices=latex_indices,
                          symbol_dual=symbol_dual,
                          latex_symbol_dual=latex_symbol_dual)
        if sections:
            linked = False
            try:
                resu._init_from_family(sections)
            except ArithmeticError as err:
                linked = str(err) in ["non-invertible matrix",
                                      "input matrix must be nonsingular"]
            if linked:
                raise ValueError("the provided sections are not linearly "
                                 "independent")
        return resu

    def section(self, *comp, **kwargs):
        r"""
        Return a continuous section of ``self``.

        INPUT:

        - ``domain`` -- (default: ``None``) domain on which the section shall be
          defined; if ``None``, the base space is assumed
        - ``name`` -- (default: ``None``) name of the local section
        - ``latex_name`` -- (default``None``) latex representation of the local
          section

        OUTPUT:

        - an instance of :class:`~sage.manifolds.section.Section` representing
          a continuous section of `M` with values on `E`

        EXAMPLES:

        A section on a non-trivial rank 2 vector bundle over a non-trivial
        2-manifold::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: E = M.vector_bundle(2, 'E') # define the vector bundle
            sage: phi_U = E.trivialization('phi_U', domain=U) # define trivializations
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: transf = phi_U.transition_map(phi_V, [[0,x],[x,0]]) # transition map between trivializations
            sage: fU = phi_U.frame(); fV = phi_V.frame() # define induced frames
            sage: s = E.section(name='s'); s
            Section s on the 2-dimensional topological manifold M with values in the
             real vector bundle E of rank 2
        """
        domain = kwargs.pop('domain', self._base_space)
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        smodule = self.section_module(domain=domain)  # the parent
        resu = smodule.element_class(smodule, name=name, latex_name=latex_name)
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    def total_space(self):
        r"""
        Return the total space of ``self``.

        .. NOTE::

            At this stage, the total space does not come with induced charts.

        OUTPUT:

        - the total space of ``self`` as an instance of
          :class:`~sage.manifolds.manifold.TopologicalManifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.total_space()
            6-dimensional topological manifold E
        """
        if self._total_space is None:
            from sage.manifolds.manifold import Manifold
            base_space = self._base_space
            dim = base_space._dim * self._rank
            sindex = base_space.start_index()
            self._total_space = Manifold(dim, self._name,
                                   latex_name=self._latex_name,
                                   field=self._field, structure='topological',
                                   start_index=sindex)

        # TODO: if update_atlas: introduce charts via self._atlas

        return self._total_space

    def set_change_of_frame(self, frame1, frame2, change_of_frame,
                            compute_inverse=True):
        r"""
        Relate two vector frames by an automorphism.

        This updates the internal dictionary ``self._frame_changes``.

        INPUT:

        - ``frame1`` -- frame 1, denoted `(e_i)` below
        - ``frame2`` -- frame 2, denoted `(f_i)` below
        - ``change_of_frame`` -- instance of class
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` -- boolean (default: ``True``); if set to True,
          the inverse automorphism is computed and the change from basis
          `(f_i)` to `(e_i)` is set to it in the internal dictionary
          ``self._frame_changes``

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: f = E.local_frame('f')
            sage: a = E.section_module().automorphism()
            sage: a[e,:] = [[1,2],[0,3]]
            sage: E.set_change_of_frame(e, f, a)
            sage: f[0].display(e)
            f_0 = e_0
            sage: f[1].display(e)
            f_1 = 2 e_0 + 3 e_1
            sage: e[0].display(f)
            e_0 = f_0
            sage: e[1].display(f)
            e_1 = -2/3 f_0 + 1/3 f_1
            sage: E.change_of_frame(e,f)[e,:]
            [1 2]
            [0 3]
        """
        from sage.tensor.modules.free_module_automorphism import \
            FreeModuleAutomorphism
        sec_module = frame1._fmodule
        if frame2._fmodule != sec_module:
            raise ValueError("the two frames are not defined on the same " +
                             "section module")
        if isinstance(change_of_frame, FreeModuleAutomorphism):
            auto = change_of_frame
        else: # Otherwise try to coerce the input
            auto_group = sec_module.general_linear_group()
            auto = auto_group(change_of_frame, basis=frame1)
        sec_module.set_change_of_basis(frame1, frame2, auto,
                                       compute_inverse=compute_inverse)
        self._frame_changes[(frame1, frame2)] = auto
        if compute_inverse:
            self._frame_changes[(frame2, frame1)] = ~auto

    def change_of_frame(self, frame1, frame2):
        r"""
        Return a change of local frames defined on ``self``.

        INPUT:

        - ``frame1`` -- local frame 1
        - ``frame2`` -- local frame 2

        OUTPUT:

        - a :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          representing, at each point, the vector space automorphism `P` that
          relates frame 1, `(e_i)` say, to frame 2, `(f_i)` say, according to
          `f_i = P(e_i)`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: a = E.section_module().automorphism() # Now, the section module is free
            sage: a[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: f = e.new_frame(a, 'f')
            sage: E.change_of_frame(e, f)
            Automorphism of the Free module C^0(M;E) of sections on the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: a == E.change_of_frame(e, f)
            True
            sage: a.inverse() == E.change_of_frame(f, e)
            True
        """
        if (frame1, frame2) not in self._frame_changes:
            raise ValueError("the change of frame from {} to {}".format(frame1, frame2) +
                             " has not been defined on the {}".format(self))
        return self._frame_changes[(frame1, frame2)]

    def changes_of_frame(self):
        r"""
        Return all the changes of local frames defined on ``self``.

        OUTPUT:

        - dictionary of vector bundle automorphisms representing
          the changes of frames, the keys being the pair of frames

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: auto_group = E.section_module().general_linear_group()
            sage: e_to_f = auto_group([[0,1],[1,0]]); e_to_f
            Automorphism of the Free module C^0(M;E) of sections on the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: f_in_e = auto_group([[0,1],[1,0]])
            sage: f = e.new_frame(f_in_e, 'f'); f
            Local frame (E|_M, (f_0,f_1))
            sage: E.changes_of_frame() # random
            {(Local frame (E|_M, (f_0,f_1)),
             Local frame (E|_M, (e_0,e_1))): Automorphism of the Free module
             C^0(M;E) of sections on the 3-dimensional topological manifold M
             with values in the real vector bundle E of rank 2,
             (Local frame (E|_M, (e_0,e_1)),
             Local frame (E|_M, (f_0,f_1))): Automorphism of the Free module
             C^0(M;E) of sections on the 3-dimensional topological manifold M
             with values in the real vector bundle E of rank 2}
        """
        return self._frame_changes.copy()

    def frames(self):
        r"""
        Return the list of local frames defined on ``self``.

        OUTPUT: list of local frames defined on ``self``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = E.local_frame('e', domain=V)
            sage: E.frames()
            [Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2))),
             Local frame (E|_V, (e_0,e_1))]
        """
        return list(self._frames)

    def coframes(self):
        r"""
        Return the list of coframes defined on ``self``.

        OUTPUT: list of coframes defined on ``self``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = E.local_frame('e', domain=V)
            sage: E.coframes()
            [Trivialization coframe (E|_U, ((phi_U^*e^1),(phi_U^*e^2))),
             Local coframe (E|_V, (e^0,e^1))]
        """
        return list(self._coframes)

    def default_frame(self):
        r"""
        Return the default frame of on ``self``.

        OUTPUT:

        - a local frame as an instance of
          :class:`~sage.manifolds.local_frame.LocalFrame`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: E.default_frame()
            Local frame (E|_M, (e_0,e_1))
        """
        return self._def_frame

    def set_default_frame(self, frame):
        r"""
        Set the default frame of ``self``.

        INPUT:

        - ``frame`` -- a local frame defined on ``self`` as an instance of
          :class:`~sage.manifolds.local_frame.LocalFrame`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: E.default_frame()
            Local frame (E|_M, (e_0,e_1))
            sage: f = E.local_frame('f')
            sage: E.set_default_frame(f)
            sage: E.default_frame()
            Local frame (E|_M, (f_0,f_1))
        """
        from sage.manifolds.local_frame import LocalFrame
        if not isinstance(frame, LocalFrame):
            raise TypeError("{} is not a local frame".format(frame))
        if not frame._domain.is_subset(self._base_space):
            raise ValueError("the frame must be defined on " +
                             "the {}".format(self))
        frame._fmodule.set_default_basis(frame)
        self._def_frame = frame

    def set_orientation(self, orientation):
        r"""
        Set the preferred orientation of ``self``.

        INPUT:

        - ``orientation`` -- a local frame or a list of local frames whose
          domains cover the base space

        .. WARNING::

            It is the user's responsibility that the orientation set here
            is indeed an orientation. There is no check going on in the
            background. See :meth:`orientation` for the definition of an
            orientation.

        EXAMPLES:

        Set an orientation on a vector bundle::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: f = E.local_frame('f'); f
            Local frame (E|_M, (f_0,f_1))
            sage: E.set_orientation(f)
            sage: E.orientation()
            [Local frame (E|_M, (f_0,f_1))]

        Set an orientation in the non-trivial case::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U); e
            Local frame (E|_U, (e_0,e_1))
            sage: f = E.local_frame('f', domain=V); f
            Local frame (E|_V, (f_0,f_1))
            sage: E.orientation()
            []
            sage: E.set_orientation([e, f])
            sage: E.orientation()
            [Local frame (E|_U, (e_0,e_1)),
             Local frame (E|_V, (f_0,f_1))]
        """
        from .local_frame import LocalFrame
        if isinstance(orientation, LocalFrame):
            orientation = [orientation]
        elif isinstance(orientation, (tuple, list)):
            orientation = list(orientation)
        else:
            raise TypeError("orientation must be a frame or a list/tuple of "
                            "frames")
        dom_union = None
        for frame in orientation:
            if frame not in self.frames():
                raise ValueError("{} must be a frame ".format(frame) +
                                 "defined on {}".format(self))
            dom = frame.domain()
            if dom_union is not None:
                dom_union = dom.union(dom_union)
            else:
                dom_union = dom
        base_space = self._base_space
        if dom_union != base_space:
            raise ValueError("the frames's domains must "
                             "cover {}".format(base_space))
        self._orientation = orientation

    def orientation(self):
        r"""
        Get the orientation of ``self`` if available.

        An *orientation* on a vector bundle is a choice of local frames whose

        1. union of domains cover the base space,
        2. changes of frames are pairwise orientation preserving, i.e. have
           positive determinant.

        A vector bundle endowed with an orientation is called *orientable*.

        The trivial case corresponds to ``self`` being trivial, i.e. ``self``
        can be covered by one frame. In that case, if no preferred
        orientation has been set before, one of those frames (usually the
        default frame) is set automatically to the preferred orientation and
        returned here.

        EXAMPLES:

        The trivial case is covered automatically::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: E.orientation()  # trivial case
            [Local frame (E|_M, (e_0,e_1))]

        The orientation can also be set by the user::

            sage: f = E.local_frame('f'); f
            Local frame (E|_M, (f_0,f_1))
            sage: E.set_orientation(f)
            sage: E.orientation()
            [Local frame (E|_M, (f_0,f_1))]

        In case of the non-trivial case, the orientation must be set manually,
        otherwise no orientation is returned::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U); e
            Local frame (E|_U, (e_0,e_1))
            sage: f = E.local_frame('f', domain=V); f
            Local frame (E|_V, (f_0,f_1))
            sage: E.orientation()
            []
            sage: E.set_orientation([e, f])
            sage: E.orientation()
            [Local frame (E|_U, (e_0,e_1)),
             Local frame (E|_V, (f_0,f_1))]
        """
        if not self._orientation:
            # Trivial case:
            if self.is_manifestly_trivial():
                # Try the default frame:
                def_frame = self._def_frame
                if def_frame is not None:
                    if def_frame._domain is self._base_space:
                        self._orientation = [def_frame]
                # Still no orientation? Choose arbitrary frame:
                if not self._orientation:
                    for frame in self.frames():
                        if frame._domain is self._base_space:
                            self._orientation = [frame]
                            break
        return list(self._orientation)

    def has_orientation(self):
        r"""
        Check whether ``self`` admits an obvious or by user set orientation.

        .. SEEALSO::

            Consult :meth:`orientation` for details about orientations.

        .. NOTE::

            Notice that if :meth:`has_orientation` returns ``False`` this does
            not necessarily mean that the vector bundle admits no orientation.
            It just means that the user has to set an orientation manually
            in that case, see :meth:`set_orientation`.

        EXAMPLES:

        The trivial case::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: E.has_orientation()  # trivial case
            True

        Non-trivial case::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: f = E.local_frame('f', domain=V)
            sage: E.has_orientation()
            False
            sage: E.set_orientation([e, f])
            sage: E.has_orientation()
            True
        """
        return bool(self.orientation())

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value `i_0` of the index;
          if none are provided, the value returned by
          :meth:`sage.manifolds.manifold.Manifold.start_index()` is assumed

        OUTPUT:

        - an iterable index, starting from `i_0` and ending at
          `i_0 + n - 1`, where `n` is the vector bundle's dimension

        EXAMPLES:

        Index range on a 4-dimensional vector bundle over a 5-dimensional
        manifold::

            sage: M = Manifold(5, 'M', structure='topological')
            sage: E = M.vector_bundle(4, 'E')
            sage: list(E.irange())
            [0, 1, 2, 3]
            sage: list(E.irange(2))
            [2, 3]

        Index range on a 4-dimensional vector bundle over a 5-dimensional
        manifold with starting index=1::

            sage: M = Manifold(5, 'M', structure='topological', start_index=1)
            sage: E = M.vector_bundle(4, 'E')
            sage: list(E.irange())
            [1, 2, 3, 4]
            sage: list(E.irange(2))
            [2, 3, 4]

        In general, one has always::

            sage: next(E.irange()) == M.start_index()
            True
        """
        si = self._base_space._sindex
        imax = self._rank + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1
