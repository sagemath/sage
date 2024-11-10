r"""
Trivializations

The class :class:`Trivialization` implements trivializations on vector bundles.
The corresponding transition maps between two trivializations are represented by
:class:`TransitionMap`.

AUTHORS:

- Michael Jung (2019) : initial version
"""

# ****************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.manifolds.local_frame import TrivializationFrame


class Trivialization(UniqueRepresentation, SageObject):
    r"""
    A local trivialization of a given vector bundle.

    Let `\pi:E \to M` be a vector bundle of rank `n` and class `C^k` over the
    field `K`
    (see :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle` or
    :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`).
    A *local trivialization* over an open subset `U \subset M` is a
    `C^k`-diffeomorphism `\varphi: \pi^{-1}(U) \to U \times K^n` such that
    `\pi \circ \varphi^{-1}=\mathrm{pr}_1` and `v \mapsto \varphi^{-1}(q,v)`
    is a linear isomorphism for any `q \in U`.

    .. NOTE::

        Notice that frames and trivializations are equivalent concepts (for
        further details see :class:`~sage.manifolds.local_frame.LocalFrame`).
        However, in order to facilitate applications and being consistent with
        the implementations of charts, trivializations are introduced
        separately.

    EXAMPLES:

    Local trivializations on a real rank 2 vector bundle over the 2-sphere::

        sage: S2 = Manifold(2, 'S^2', structure='top')
        sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
        sage: S2.declare_union(U,V)
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: W = U.intersection(V)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: E = S2.vector_bundle(2, 'E')
        sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
        ....:                          domain=U); phi_U
        Trivialization (phi_U, E|_U)
        sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
        ....:                          domain=V); phi_V
        Trivialization (phi_V, E|_V)
        sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]]); phi_U_to_phi_V
        Transition map from Trivialization (phi_U, E|_U) to Trivialization
         (phi_V, E|_V)

    The LaTeX output gives the following::

        sage: latex(phi_U)
        \varphi_U : E |_{U} \to U \times \Bold{R}^2
        sage: latex(phi_V)
        \varphi_V : E |_{V} \to V \times \Bold{R}^2

    The trivializations are part of the vector bundle atlas::

        sage: E.atlas()
        [Trivialization (phi_U, E|_U), Trivialization (phi_V, E|_V)]

    Each trivialization induces a local trivialization frame::

        sage: fU = phi_U.frame(); fU
        Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2)))
        sage: fV = phi_V.frame(); fV
        Trivialization frame (E|_V, ((phi_V^*e_1),(phi_V^*e_2)))

    and the transition map connects these two frames via a bundle automorphism::

        sage: aut = phi_U_to_phi_V.automorphism(); aut
        Automorphism phi_U^(-1)*phi_V of the Free module C^0(W;E) of sections on
         the Open subset W of the 2-dimensional topological manifold S^2 with
         values in the real vector bundle E of rank 2
        sage: aut.display(fU.restrict(W))
        phi_U^(-1)*phi_V = (phi_U^*e_1)⊗(phi_U^*e^2) + (phi_U^*e_2)⊗(phi_U^*e^1)
        sage: aut.display(fV.restrict(W))
        phi_U^(-1)*phi_V = (phi_V^*e_1)⊗(phi_V^*e^2) + (phi_V^*e_2)⊗(phi_V^*e^1)

    The automorphisms are listed in the frame changes of the vector bundle::

        sage: E.changes_of_frame() # random
        {(Local frame (E|_W, ((phi_U^*e_1),(phi_U^*e_2))),
         Local frame (E|_W, ((phi_V^*e_1),(phi_V^*e_2)))): Automorphism
         phi_U^(-1)*phi_V^(-1) of the Free module C^0(W;E) of sections on the
         Open subset W of the 2-dimensional topological manifold S^2 with values
         in the real vector bundle E of rank 2,
         (Local frame (E|_W, ((phi_V^*e_1),(phi_V^*e_2))),
         Local frame (E|_W, ((phi_U^*e_1),(phi_U^*e_2)))): Automorphism
         phi_U^(-1)*phi_V of the Free module C^0(W;E) of sections on the Open
         subset W of the 2-dimensional topological manifold S^2 with values in
         the real vector bundle E of rank 2}

    Let us check the components of ``fU`` with respect to the frame ``fV``::

        sage: fU[0].comp(fV.restrict(W))[:]
        [0, 1]
        sage: fU[1].comp(fV.restrict(W))[:]
        [1, 0]
    """

    def __init__(self, vector_bundle, name, domain, latex_name=None):
        r"""
        Construct a local trivialization of the vector bundle ``vector_bundle``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: TestSuite(phi).run()
        """
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._base_space = vector_bundle.base_space()
        self._vbundle = vector_bundle
        self._bdl_rank = vector_bundle.rank()
        self._base_field = vector_bundle.base_field()
        self._sindex = self._base_space.start_index()
        self._domain = domain
        # Add this trivialization to the atlas of the vector bundle:
        vector_bundle._atlas.append(self)
        self._frame = TrivializationFrame(self)
        self._coframe = self._frame._coframe

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: phi = E.trivialization('phi', domain=M)
            sage: phi._repr_()
            'Trivialization (phi, E|_M)'
        """
        desc = "Trivialization ("
        desc += self._name + ", "
        desc += "{}|_{})".format(self._vbundle._name, self._domain._name)
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: phi = E.trivialization('phi', domain=M, latex_name=r'\varphi')
            sage: phi._latex_()
            '\\varphi : E |_{M} \\to M \\times \\Bold{R}^1'
        """
        latex = self._latex_name + r' : '
        latex += r'{} |_{{{}}} \to {} \times {}^{}'.format(self._vbundle._latex_name,
                            self._domain._latex_(), self._domain._latex_(),
                            self._base_field._latex_(), self._bdl_rank)
        return latex

    def base_space(self):
        r"""
        Return the manifold on which the trivialization is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U)
            sage: phi.base_space()
            2-dimensional topological manifold M
        """
        return self._base_space

    def transition_map(self, other, transf, compute_inverse=True):
        r"""
        Return the transition map between ``self`` and ``other``.

        INPUT:

        - ``other`` -- the trivialization where the transition map from ``self``
          goes to
        - ``transf`` -- transformation of the transition map
        - ``intersection_name`` -- (default: ``None``) name to be given to the
          subset `U \cap V` if the latter differs from `U` or `V`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(U)
            sage: W = U.intersection(V)
            sage: XW = X.restrict(W)
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U.transition_map(phi_V, 1)
            Transition map from Trivialization (phi_U, E|_U) to Trivialization
             (phi_V, E|_V)
        """
        return TransitionMap(self, other, transf,
                             compute_inverse=compute_inverse)

    def vector_bundle(self):
        r"""
        Return the vector bundle on which the trivialization is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U)
            sage: phi.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
             2-dimensional topological manifold M
        """
        return self._vbundle

    def domain(self):
        r"""
        Return the domain on which the trivialization is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U)
            sage: phi.domain()
            Open subset U of the 2-dimensional topological manifold M
        """
        return self._domain

    def frame(self):
        r"""
        Return the standard frame induced by ``self``. If `\psi` is a
        trivialization then the corresponding frame can be obtained by the maps
        `p \mapsto \psi^{-1}(p,e_i)`, where `(e_1, \ldots, e_n)` is the standard
        basis of `K^n`. We briefly denote `(\psi^*e_i)` instead of
        `\psi^{-1}(\cdot,e_i)`.

        .. SEEALSO::

            :class:`~sage.manifolds.local_frame.LocalFrame`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: phi.frame()
            Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2)))
        """
        return self._frame

    def coframe(self):
        r"""
        Return the standard coframe induced by ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.local_frame.LocalCoFrame`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: phi.coframe()
            Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2)))
        """
        return self._frame._coframe

# *****************************************************************************


class TransitionMap(SageObject):
    r"""
    Transition map between two trivializations.

    Given a vector bundle `\pi : E \to M` of class `C^k` and rank `n` over the
    field `K`, and two trivializations
    `\varphi_U : \pi^{-1}(U) \to U \times K^n` and
    `\varphi_V : \pi^{-1}(V) \to V \times K^n`, the transition map from
    `\varphi_U` to `\varphi_V` is given by the composition

    .. MATH::

        \varphi_V \circ \varphi_U^{-1} : U \cap V \times K^n \to
        U \cap V \times K^n .

    This composition is of the form

    .. MATH::

        (p, v) \mapsto (p, g(p)v),

    where `p \mapsto g(p)` is a `C^k` family of invertible `n\times n` matrices.

    INPUT:

    - ``triv1`` -- trivialization 1
    - ``triv2`` -- trivialization 2
    - ``transf`` -- the transformation between both trivializations in form of a
      matrix of scalar fields (:class:`~sage.manifolds.scalarfield.ScalarField`)
      or coordinate functions (:class:`~sage.manifolds.chart_func.ChartFunction`),
      or a bundle automorphism
      (:class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`)
    - ``compute_inverse`` -- boolean (default: ``True``); determines whether the inverse
      shall be computed or not

    EXAMPLES:

    Transition map of two trivializations on a real rank 2 vector bundle of the
    2-sphere::

        sage: S2 = Manifold(2, 'S^2', structure='top')
        sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
        sage: S2.declare_union(U,V)
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: W = U.intersection(V)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: E = S2.vector_bundle(2, 'E')
        sage: phi_U = E.trivialization('phi_U', domain=U)
        sage: phi_V = E.trivialization('phi_V', domain=V)
        sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]])
        sage: phi_U_to_phi_V
        Transition map from Trivialization (phi_U, E|_U) to Trivialization
         (phi_V, E|_V)
    """
    def __init__(self, triv1, triv2, transf, compute_inverse=True):
        r"""
        Construct a transition map between two trivializations.

        TESTS::

            sage: S2 = Manifold(2, 'S^2', structure='top')
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: W = U.intersection(V)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: E = S2.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]])
            sage: TestSuite(phi_U_to_phi_V).run()
        """
        bs1 = triv1.base_space()
        bs2 = triv2.base_space()
        if bs1 is not bs2:
            raise ValueError("base spaces must coincide")
        self._base_space = bs1

        vb1 = triv1.vector_bundle()
        vb2 = triv2.vector_bundle()
        if vb1 is not vb2:
            raise ValueError("vector bundles must coincide")
        self._vbundle = vb1
        self._bdl_rank = self._vbundle.rank()

        dom1 = triv1.domain()
        dom2 = triv2.domain()
        dom = dom1.intersection(dom2)
        self._domain = dom
        self._frame1 = triv1._frame.restrict(dom)
        self._frame2 = triv2._frame.restrict(dom)
        self._triv1 = triv1
        self._triv2 = triv2
        self._inverse = None
        self._vbundle._transitions[(triv1, triv2)] = self
        self._name = triv2._name + "*" + triv1._name + "^(-1)"
        self._latex_name = triv2._latex_name + r'\circ ' + triv1._latex_name + \
                           r'^{-1}'
        ###
        # Define the automorphism
        auto_name = triv1._name + "^(-1)*" + triv2._name
        auto_lname = triv1._latex_name + r'^{-1} \circ ' + triv2._latex_name
        sec_module = self._vbundle.section_module(dom, force_free=True)
        auto_group = sec_module.general_linear_group()
        auto = auto_group(transf, basis=self._frame1, name=auto_name,
                          latex_name=auto_lname)
        self._automorphism = auto
        # Add this change of basis to the basis changes
        self._vbundle.set_change_of_frame(self._frame2, self._frame1, auto,
                                          compute_inverse=compute_inverse)
        if compute_inverse:
            self._inverse = type(self)(self._triv2, self._triv1, ~auto,
                                       compute_inverse=False)
            self._inverse._inverse = self

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(U)
            sage: W = U.intersection(V)
            sage: XW = X.restrict(W)
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: phi_U_to_phi_V._repr_()
            'Transition map from Trivialization (phi_U, E|_U) to Trivialization
             (phi_V, E|_V)'
            sage: repr(phi_U_to_phi_V)
            'Transition map from Trivialization (phi_U, E|_U) to Trivialization
             (phi_V, E|_V)'
            sage: phi_U_to_phi_V # indirect doctest
            Transition map from Trivialization (phi_U, E|_U) to Trivialization
             (phi_V, E|_V)
        """
        desc = "Transition map from {} to {}".format(self._triv1, self._triv2)
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(U)
            sage: W = U.intersection(V)
            sage: XW = X.restrict(W)
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: phi_U_to_phi_V._latex_()
            '\\varphi_V\\circ \\varphi_U^{-1}:U\\cap V\\times \\Bold{R}^2 \\to
             U\\cap V\\times \\Bold{R}^2'
            sage: latex(phi_U_to_phi_V)
            \varphi_V\circ \varphi_U^{-1}:U\cap V\times \Bold{R}^2 \to U\cap
             V\times \Bold{R}^2
        """
        vspace_lname = r'\times {}^{}'.format(self._triv1._base_field._latex_(),
                                              self._bdl_rank)
        latex = self._latex_name + r':'
        latex += self._domain._latex_name
        latex += vspace_lname + r' \to '
        latex += self._domain._latex_name
        latex += vspace_lname
        return latex

    def automorphism(self):
        r"""
        Return the automorphism connecting both trivializations.

        The family of matrices `p \mapsto g(p)` given by the transition map
        induce a bundle automorphism

        .. MATH::

            \varphi_U^{-1} \circ \varphi_V : \pi^{-1}(U \cap V) \to
            \pi^{-1}(U \cap V)

        correlating the local frames induced by the trivializations in the
        following way:

        .. MATH::

           (\varphi_U^{-1} \circ \varphi_V) (\varphi_V^*e_i) = \varphi_U^*e_i .

        Then, for each point `p \in M`, the matrix `g(p)` is the representation
        of the induced automorphism on the fiber `E_p=\pi^{-1}(p)` in the basis
        `\left( (\varphi_V^*e_i)(p)\right)_{i=1,\dots,n}`.

        EXAMPLES::

            sage: S2 = Manifold(2, 'S^2', structure='top')
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: W = U.intersection(V)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: E = S2.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U); phi_U
            Trivialization (phi_U, E|_U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V); phi_V
            Trivialization (phi_V, E|_V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]])
            sage: aut = phi_U_to_phi_V.automorphism(); aut
            Automorphism phi_U^(-1)*phi_V of the Free module C^0(W;E) of
             sections on the Open subset W of the 2-dimensional topological
             manifold S^2 with values in the real vector bundle E of rank 2
            sage: aut.display(phi_U.frame().restrict(W))
            phi_U^(-1)*phi_V = (phi_U^*e_1)⊗(phi_U^*e^2) +
             (phi_U^*e_2)⊗(phi_U^*e^1)
        """
        return self._automorphism

    def inverse(self):
        r"""
        Return the inverse transition map.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(U)
            sage: W = U.intersection(V)
            sage: XW = X.restrict(W)
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[1,1],[-1,1]],
            ....:                                       compute_inverse=False)
            sage: phi_V_to_phi_U = phi_U_to_phi_V.inverse(); phi_V_to_phi_U
            Transition map from Trivialization (phi_V, E|_V) to Trivialization (phi_U, E|_U)
            sage: phi_V_to_phi_U.automorphism() == phi_U_to_phi_V.automorphism().inverse()
            True
        """
        if self._inverse is None:
            self._vbundle.set_change_of_frame(self._frame1, self._frame2,
                                              ~self._automorphism)
            self._inverse = type(self)(self._triv2, self._triv1,
                                       ~self._automorphism,
                                       compute_inverse=False)
            self._inverse._inverse = self
        return self._inverse

    def det(self):
        r"""
        Return the determinant of ``self``.

        OUTPUT: an instance of :class:`~sage.manifolds.scalarfield.ScalarField`

        EXAMPLES::

            sage: S2 = Manifold(2, 'S^2', structure='top')
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: W = U.intersection(V)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: E = S2.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U); phi_U
            Trivialization (phi_U, E|_U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V); phi_V
            Trivialization (phi_V, E|_V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]])
            sage: det = phi_U_to_phi_V.det(); det
            Scalar field det(phi_U^(-1)*phi_V) on the Open subset W of the
             2-dimensional topological manifold S^2
            sage: det.display()
            det(phi_U^(-1)*phi_V): W → ℝ
                 (x, y) ↦ -1
                 (u, v) ↦ -1
        """
        aut = self._automorphism
        det = aut.det()
        name = "det({})".format(aut._name)
        latex_name = r'\det({})'.format(aut._latex_name)
        det.set_name(name=name, latex_name=latex_name)
        return det

    def matrix(self):
        r"""
        Return the matrix representation the transition map.

        EXAMPLES:

        Local trivializations on a real rank 2 vector bundle over the 2-sphere::

            sage: S2 = Manifold(2, 'S^2', structure='top')
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: W = U.intersection(V)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: E = S2.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U); phi_U
            Trivialization (phi_U, E|_U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V); phi_V
            Trivialization (phi_V, E|_V)

        The input is coerced into a bundle automorphism. From there, the matrix
        can be recovered::

            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[0,1],[1,0]])
            sage: matrix = phi_U_to_phi_V.matrix(); matrix
            [Scalar field zero on the Open subset W of the 2-dimensional
             topological manifold S^2      Scalar field 1 on the Open subset
             W of the 2-dimensional topological manifold S^2]
             [     Scalar field 1 on the Open subset W of the 2-dimensional
             topological manifold S^2 Scalar field zero on the Open subset W of
             the 2-dimensional topological manifold S^2]

        Let us check the matrix components::

            sage: matrix[0,0].display()
            zero: W → ℝ
                (x, y) ↦ 0
                (u, v) ↦ 0
            sage: matrix[0,1].display()
            1: W → ℝ
             (x, y) ↦ 1
             (u, v) ↦ 1
            sage: matrix[1,0].display()
            1: W → ℝ
             (x, y) ↦ 1
             (u, v) ↦ 1
            sage: matrix[1,1].display()
            zero: W → ℝ
                (x, y) ↦ 0
                (u, v) ↦ 0
        """
        return self._automorphism.matrix(self._frame1)

    def __eq__(self, other):
        r"""
        Equality operator.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(V)
            sage: XUV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[1,2],[3,-2]])
            sage: phi_U_to_phi_V == phi_U_to_phi_V
            True
            sage: phi_U_to_phi_V_1 = phi_U.transition_map(phi_V, [[1,2],[3,-2]])
            sage: phi_U_to_phi_V == phi_U_to_phi_V_1
            True
            sage: phi_U_to_phi_V_2 = phi_U.transition_map(phi_V, [[1,2],[2,-1]])
            sage: phi_U_to_phi_V == phi_U_to_phi_V_2
            False
            sage: psi_V = E.trivialization('psi_V', domain=V)
            sage: phi_U_to_psi_V = phi_U.transition_map(psi_V, [[1,2],[3,-2]])
            sage: phi_U_to_phi_V == phi_U_to_psi_V
            False
        """
        if other is self:
            return True
        if not isinstance(other, TransitionMap):
            return False
        return ((self._triv1 == other._triv1)
                and (self._triv2 == other._triv2)
                and (self._automorphism == other._automorphism))

    def __ne__(self, other):
        r"""
        Non-equality operator.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: XU = X.restrict(U); XV = X.restrict(V)
            sage: XUV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, [[1,2],[3,-2]])
            sage: phi_U_to_phi_V_1 = phi_U.transition_map(phi_V, [[1,2],[2,-1]])
            sage: phi_U_to_phi_V != phi_U_to_phi_V_1
            True
        """
        return not (self == other)
