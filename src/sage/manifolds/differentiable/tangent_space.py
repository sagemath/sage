r"""
Tangent Spaces

The class :class:`TangentSpace` implements tangent vector spaces to a
differentiable manifold.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- Chap. 3 of [Lee2013]_
"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING

from sage.manifolds.differentiable.tangent_vector import TangentVector
from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule

if TYPE_CHECKING:
    from sage.manifolds.point import ManifoldPoint


class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space to a differentiable manifold at a given point.

    Let `M` be a differentiable manifold of dimension `n` over a
    topological field `K` and `p \in M`. The tangent space `T_p M` is an
    `n`-dimensional vector space over `K` (without a distinguished basis).

    INPUT:

    - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
      point `p` at which the tangent space is defined

    EXAMPLES:

    Tangent space on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((-1,2), name='p')
        sage: Tp = M.tangent_space(p) ; Tp
        Tangent space at Point p on the 2-dimensional differentiable manifold M

    Tangent spaces are free modules of finite rank over
    :class:`~sage.symbolic.ring.SymbolicRing`
    (actually vector spaces of finite dimension over the manifold base
    field `K`, with `K=\RR` here)::

        sage: Tp.base_ring()
        Symbolic Ring
        sage: Tp.category()
        Category of finite dimensional vector spaces over Symbolic Ring
        sage: Tp.rank()
        2
        sage: dim(Tp)
        2

    The tangent space is automatically endowed with bases deduced from the
    vector frames around the point::

        sage: Tp.bases()
        [Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M]
        sage: M.frames()
        [Coordinate frame (M, (∂/∂x,∂/∂y))]

    At this stage, only one basis has been defined in the tangent space, but
    new bases can be added from vector frames on the manifold by means of the
    method :meth:`~sage.manifolds.differentiable.vectorframe.VectorFrame.at`,
    for instance, from the frame associated with some new coordinates::

        sage: c_uv.<u,v> = M.chart()
        sage: c_uv.frame().at(p)
        Basis (∂/∂u,∂/∂v) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M
        sage: Tp.bases()
        [Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M,
         Basis (∂/∂u,∂/∂v) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M]

    All the bases defined on ``Tp`` are on the same footing. Accordingly the
    tangent space is not in the category of modules with a distinguished
    basis::

        sage: Tp in ModulesWithBasis(SR)
        False

    It is simply in the category of modules::

        sage: Tp in Modules(SR)
        True

    Since the base ring is a field, it is actually in the category of
    vector spaces::

        sage: Tp in VectorSpaces(SR)
        True

    A typical element::

        sage: v = Tp.an_element() ; v
        Tangent vector at Point p on the
         2-dimensional differentiable manifold M
        sage: v.display()
        ∂/∂x + 2 ∂/∂y
        sage: v.parent()
        Tangent space at Point p on the
         2-dimensional differentiable manifold M

    The zero vector::

        sage: Tp.zero()
        Tangent vector zero at Point p on the
         2-dimensional differentiable manifold M
        sage: Tp.zero().display()
        zero = 0
        sage: Tp.zero().parent()
        Tangent space at Point p on the
         2-dimensional differentiable manifold M

    Tangent spaces are unique::

        sage: M.tangent_space(p) is Tp
        True
        sage: p1 = M.point((-1,2))
        sage: M.tangent_space(p1) is Tp
        True

    even if points are not::

        sage: p1 is p
        False

    Actually ``p1`` and ``p`` share the same tangent space because they
    compare equal::

        sage: p1 == p
        True

    The tangent-space uniqueness holds even if the points are created in
    different coordinate systems::

        sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y))
        sage: uv_to_xv = xy_to_uv.inverse()
        sage: p2 = M.point((1, -3), chart=c_uv, name='p_2')
        sage: p2 is p
        False
        sage: M.tangent_space(p2) is Tp
        True
        sage: p2 == p
        True

    An isomorphism of the tangent space with an inner product space with distinguished basis::

        sage: g = M.metric('g')
        sage: g[:] = ((1, 0), (0, 1))
        sage: Q_Tp_xy = g[c_xy.frame(),:](*p.coordinates(c_xy)); Q_Tp_xy
        [1 0]
        [0 1]
        sage: W_Tp_xy = VectorSpace(SR, 2, inner_product_matrix=Q_Tp_xy)
        sage: Tp.bases()[0]
        Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the
         2-dimensional differentiable manifold M
        sage: phi_Tp_xy = Tp.isomorphism_with_fixed_basis(Tp.bases()[0], codomain=W_Tp_xy)
        sage: phi_Tp_xy
        Generic morphism:
         From: Tangent space at Point p on the 2-dimensional differentiable manifold M
         To:   Ambient quadratic space of dimension 2 over Symbolic Ring
               Inner product matrix:
               [1 0]
               [0 1]

        sage: Q_Tp_uv = g[c_uv.frame(),:](*p.coordinates(c_uv)); Q_Tp_uv
        [1/2   0]
        [  0 1/2]
        sage: W_Tp_uv = VectorSpace(SR, 2, inner_product_matrix=Q_Tp_uv)
        sage: Tp.bases()[1]
        Basis (∂/∂u,∂/∂v) on the Tangent space at Point p on the
         2-dimensional differentiable manifold M
        sage: phi_Tp_uv = Tp.isomorphism_with_fixed_basis(Tp.bases()[1], codomain=W_Tp_uv)
        sage: phi_Tp_uv
        Generic morphism:
         From: Tangent space at Point p on the 2-dimensional differentiable manifold M
         To:   Ambient quadratic space of dimension 2 over Symbolic Ring
               Inner product matrix:
               [1/2   0]
               [  0 1/2]

        sage: t1, t2 = Tp.tensor((1,0)), Tp.tensor((1,0))
        sage: t1[:] = (8, 15)
        sage: t2[:] = (47, 11)
        sage: t1[Tp.bases()[0],:]
        [8, 15]
        sage: phi_Tp_xy(t1), phi_Tp_xy(t2)
        ((8, 15), (47, 11))
        sage: phi_Tp_xy(t1).inner_product(phi_Tp_xy(t2))
        541

        sage: Tp_xy_to_uv = M.change_of_frame(c_xy.frame(), c_uv.frame()).at(p); Tp_xy_to_uv
        Automorphism of the Tangent space at Point p on the
         2-dimensional differentiable manifold M
        sage: Tp.set_change_of_basis(Tp.bases()[0], Tp.bases()[1], Tp_xy_to_uv)
        sage: t1[Tp.bases()[1],:]
        [23, -7]
        sage: phi_Tp_uv(t1), phi_Tp_uv(t2)
        ((23, -7), (58, 36))
        sage: phi_Tp_uv(t1).inner_product(phi_Tp_uv(t2))
        541

    .. SEEALSO::

        :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
        for more documentation.
    """
    Element = TangentVector

    def __init__(self, point: ManifoldPoint, base_ring=None):
        r"""
        Construct the tangent space at a given point.

        TESTS::

            sage: from sage.manifolds.differentiable.tangent_space import TangentSpace
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = TangentSpace(p) ; Tp
            Tangent space at Point p on the 2-dimensional differentiable
             manifold M
            sage: TestSuite(Tp).run()
        """
        manif = point._manifold
        name = "T_{} {}".format(point._name, manif._name)
        latex_name = r"T_{%s}\,%s" % (point._latex_name, manif._latex_name)
        self._point = point
        self._manif = manif
        if base_ring is None:
            base_ring = SR
        FiniteRankFreeModule.__init__(self, base_ring, manif._dim, name=name,
                                      latex_name=latex_name,
                                      start_index=manif._sindex)
        # Initialization of bases of the tangent space from existing vector
        # frames around the point:

        # dictionary of bases of the tangent space derived from vector
        # frames around the point (keys: vector frames)
        self._frame_bases = {}
        for frame in point.parent()._top_frames:
            # the frame is used to construct a basis of the tangent space
            # only if it is a frame on the current manifold:
            if frame.destination_map().is_identity():
                if point in frame._domain:
                    coframe = frame.coframe()
                    basis = self.basis(frame._symbol,
                                       latex_symbol=frame._latex_symbol,
                                       indices=frame._indices,
                                       latex_indices=frame._latex_indices,
                                       symbol_dual=coframe._symbol,
                                       latex_symbol_dual=coframe._latex_symbol)
                    self._frame_bases[frame] = basis
        # The basis induced by the default frame of the manifold subset
        # in which the point has been created is declared the default
        # basis of self:
        def_frame = point.parent()._def_frame
        if def_frame in self._frame_bases:
            self._def_basis = self._frame_bases[def_frame]
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in point.parent()._frame_changes.items():
            if point in automorph.domain():
                frame1, frame2 = frame_pair[0], frame_pair[1]
                fr1, fr2 = None, None
                for frame in self._frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                for frame in self._frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr1 is not None and fr2 is not None:
                    basis1 = self._frame_bases[fr1]
                    basis2 = self._frame_bases[fr2]
                    auto = self.automorphism()
                    for frame, comp in automorph._components.items():
                        try:
                            basis = None
                            if frame is frame1:
                                basis = basis1
                            if frame is frame2:
                                basis = basis2
                            if basis is not None:
                                cauto = auto.add_comp(basis=basis)
                                for ind, val in comp._comp.items():
                                    cauto._comp[ind] = val(point)
                        except ValueError:
                            pass
                    self._basis_changes[(basis1, basis2)] = auto

    def construction(self):
        """
        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp.construction() is None
            True
        """
        return None

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp
            Tangent space at Point p on the
             2-dimensional differentiable manifold M
        """
        return "Tangent space at {}".format(self._point)

    def _an_element_(self):
        r"""
        Construct some (unnamed) vector in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp._an_element_()
            Tangent vector at Point p on the 2-dimensional differentiable
             manifold M
            sage: Tp._an_element_().display()
            ∂/∂x + 2 ∂/∂y
        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank+1)
        return resu

    def dimension(self):
        r"""
        Return the vector space dimension of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp.dimension()
            2

        A shortcut is ``dim()``::

            sage: Tp.dim()
            2

        One can also use the global function ``dim``::

            sage: dim(Tp)
            2
        """
        # The dimension is the rank of self as a free module:
        return self._rank

    dim = dimension

    def base_point(self):
        r"""
        Return the manifold point at which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp.base_point()
            Point p on the 2-dimensional differentiable manifold M
            sage: Tp.base_point() is p
            True
        """
        return self._point
