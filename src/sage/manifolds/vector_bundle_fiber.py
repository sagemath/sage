# sage_setup: distribution = sagemath-symbolics
r"""
Vector Bundle Fibers

The class :class:`VectorBundleFiber` implements fibers over a vector bundle.

AUTHORS:

- Michael Jung (2019): initial version
"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.vector_bundle_fiber_element import VectorBundleFiberElement

class VectorBundleFiber(FiniteRankFreeModule):
    r"""
    Fiber of a given vector bundle at a given point.

    Let `\pi: E \to M` be a vector bundle of rank `n` over the field `K`
    (see :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle`) and
    `p \in M`. The fiber `E_p` at `p` is defined via `E_p := \pi^{-1}(p)` and
    takes the structure of an `n`-dimensional vector space over the field `K`.

    INPUT:

    - ``vector_bundle`` -- :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle`;
      vector bundle `E` on which the fiber is defined
    - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
      point `p` at which the fiber is defined

    EXAMPLES:

    A vector bundle fiber in a trivial rank 2 vector bundle over a
    4-dimensional topological manifold::

        sage: M = Manifold(4, 'M', structure='top')
        sage: X.<x,y,z,t> = M.chart()
        sage: p = M((0,0,0,0), name='p')
        sage: E = M.vector_bundle(2, 'E')
        sage: e = E.local_frame('e')
        sage: Ep = E.fiber(p); Ep
        Fiber of E at Point p on the 4-dimensional topological manifold M

    Fibers are free modules of finite rank over
    :class:`~sage.symbolic.ring.SymbolicRing`
    (actually vector spaces of finite dimension over the vector bundle
    field `K`, here `K=\RR`)::

        sage: Ep.base_ring()
        Symbolic Ring
        sage: Ep.category()
        Category of finite dimensional vector spaces over Symbolic Ring
        sage: Ep.rank()
        2
        sage: dim(Ep)
        2

    The fiber is automatically endowed with bases deduced from the local frames
    around the point::

        sage: Ep.bases()
        [Basis (e_0,e_1) on the Fiber of E at Point p on the 4-dimensional
         topological manifold M]
        sage: E.frames()
        [Local frame (E|_M, (e_0,e_1))]

    At this stage, only one basis has been defined in the fiber, but new bases
    can be added from local frames on the vector bundle by means of the method
    :meth:`~sage.manifolds.local_frame.LocalFrame.at`::

        sage: aut = E.section_module().automorphism()
        sage: aut[:] = [[-1, x], [y, 2]]
        sage: f = e.new_frame(aut, 'f')
        sage: fp = f.at(p); fp
        Basis (f_0,f_1) on the Fiber of E at Point p on the 4-dimensional
         topological manifold M
        sage: Ep.bases()
        [Basis (e_0,e_1) on the Fiber of E at Point p on the 4-dimensional
         topological manifold M,
         Basis (f_0,f_1) on the Fiber of E at Point p on the 4-dimensional
         topological manifold M]

    The changes of bases are applied to the fibers::

        sage: f[1].display(e) # second component of frame f
        f_1 = x e_0 + 2 e_1
        sage: ep = e.at(p)
        sage: fp[1].display(ep) # second component of frame f at p
        f_1 = 2 e_1

    All the bases defined on ``Ep`` are on the same footing. Accordingly the
    fiber is not in the category of modules with a distinguished basis::

        sage: Ep in ModulesWithBasis(SR)
        False

    It is simply in the category of modules::

        sage: Ep in Modules(SR)
        True

    Since the base ring is a field, it is actually in the category of
    vector spaces::

        sage: Ep in VectorSpaces(SR)
        True

    A typical element::

        sage: v = Ep.an_element(); v
        Vector in the fiber of E at Point p on the 4-dimensional topological
         manifold M
        sage: v.display()
        e_0 + 2 e_1
        sage: v.parent()
        Fiber of E at Point p on the 4-dimensional topological manifold M

    The zero vector::

        sage: Ep.zero()
        Vector zero in the fiber of E at Point p on the 4-dimensional
         topological manifold M
        sage: Ep.zero().display()
        zero = 0
        sage: Ep.zero().parent()
        Fiber of E at Point p on the 4-dimensional topological manifold M

    Fibers are unique::

        sage: E.fiber(p) is Ep
        True
        sage: p1 = M.point((0,0,0,0))
        sage: E.fiber(p1) is Ep
        True

    even if points are different instances::

        sage: p1 is p
        False

    but ``p1`` and ``p`` share the same fiber because they compare equal::

        sage: p1 == p
        True

    .. SEEALSO::

        :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
        for more documentation.
    """
    Element = VectorBundleFiberElement

    def __init__(self, vector_bundle, point):
        r"""
        Construct a fiber of a vector bundle.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: Ep = E.fiber(p)
            sage: TestSuite(Ep).run()
        """
        if point._manifold is not vector_bundle._base_space:
            raise ValueError("Point must be an element "
                             "of {}".format(vector_bundle._manifold))
        name = "{}_{}".format(vector_bundle._name, point._name)
        latex_name = r'{}_{{{}}}'.format(vector_bundle._latex_name,
                                         point._latex_name)
        self._rank = vector_bundle._rank
        self._vbundle = vector_bundle
        self._point = point
        self._base_space = point._manifold
        FiniteRankFreeModule.__init__(self, SR, self._rank, name=name,
                                      latex_name=latex_name,
                                      start_index=self._base_space._sindex)
        ###
        # Construct basis
        self._frame_bases = {} # dictionary of bases of the vector bundle fiber
                        # derived from local frames around the point
                        # (keys: local frames)
        self._def_basis = None
        for frame in vector_bundle._frames:
            # the frame is used to construct a basis of the vector bundle fiber
            # only if it is a frame for the given point:
            if point in frame.domain():
                coframe = frame.coframe()
                basis = self.basis(frame._symbol,
                                   latex_symbol=frame._latex_symbol,
                                   indices=frame._indices,
                                   latex_indices=frame._latex_indices,
                                   symbol_dual=coframe._symbol,
                                   latex_symbol_dual=coframe._latex_symbol)
                self._frame_bases[frame] = basis
                if self._def_basis is None:
                    self._def_basis = basis # Declare the first basis as default
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in self._vbundle._frame_changes.items():
            if point in automorph._fmodule.domain():
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
                                cauto = auto.add_comp(basis)
                                for ind, val in comp._comp.items():
                                    cauto._comp[ind] = val(point)
                        except ValueError:
                            pass
                    self._basis_changes[(basis1, basis2)] = auto

    def construction(self):
        r"""
        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.fiber(p).construction() is None
            True
        """
        return None

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.fiber(p)._repr_()
            'Fiber of E at Point p on the 3-dimensional topological manifold M'
        """
        return "Fiber of {} at {}".format(self._vbundle._name,
                                          self._point)

    def dimension(self):
        r"""
        Return the vector space dimension of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: Ep = E.fiber(p)
            sage: Ep.dim()
            2
        """
        # The dimension is the rank of self as a free module:
        return self._rank

    dim = dimension

    def _an_element_(self):
        r"""
        Construct some (unnamed) vector in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: p = M.point((3,-2), name='p')
            sage: Ep = E.fiber(p)
            sage: Ep._an_element_()
            Vector in the fiber of E at Point p on the 2-dimensional topological
             manifold M
            sage: Ep._an_element_().display()
            e_0 + 2 e_1
        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank + 1)
        return resu

    def base_point(self):
        r"""
        Return the manifold point over which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: p = M.point((3,-2), name='p')
            sage: Ep = E.fiber(p)
            sage: Ep.base_point()
            Point p on the 2-dimensional topological manifold M
            sage: p is Ep.base_point()
            True
        """
        return self._point
