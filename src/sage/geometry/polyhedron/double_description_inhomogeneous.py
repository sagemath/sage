r"""
Double Description for Arbitrary Polyhedra

This module is part of the python backend for polyhedra. It uses the
double description method for cones
:mod:`~sage.geometry.polyhedron.double_description` to find minimal
H/V-representations of polyhedra. The latter works with cones
only. This is sufficient to treat general polyhedra by the following
construction: Any polyhedron can be embedded in one dimension higher
in the hyperplane `(1,*,\dots,*)`. The cone over the embedded
polyhedron will be called the *homogenized cone* in the
following. Conversely, intersecting the homogenized cone with the
hyperplane `x_0=1` gives you back the original polyhedron.

While slower than specialized C/C++ implementations, the
implementation is general and works with any field in Sage that allows
you to define polyhedra.

.. note::

    If you just want polyhedra over arbitrary fields then you should
    just use the
    :func:`~sage.geometry.polyhedron.constructor.Polyhedron`
    constructor.

EXAMPLES::

    sage: from sage.geometry.polyhedron.double_description_inhomogeneous \
    ....:     import Hrep2Vrep, Vrep2Hrep
    sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,4,3)], [])
    [-1/2|-1/2  1/2|]
    [   0| 2/3 -1/3|]

Note that the columns of the printed matrix are the vertices, rays,
and lines of the minimal V-representation. Dually, the rows of the
following are the inequalities and equations::

    sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
    [1 2 3]
    [2 4 3]
    [-----]
"""

# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.geometry.polyhedron.double_description import StandardAlgorithm as Algorithm

# Compare with PPL if the base ring is QQ. Can be left enabled since
# we don't use the Python fallback for polyhedra over QQ unless you
# construct one by hand.
VERIFY_RESULT = True


class PivotedInequalities(SageObject):

    def __init__(self, base_ring, dim):
        """
        Base class for inequalities that may contain linear subspaces

        INPUT:

        - ``base_ring`` -- a field.

        - ``dim`` -- integer. The ambient space dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous \
            ....:     import PivotedInequalities
            sage: piv = PivotedInequalities(QQ, 2)
            sage: piv._pivot_inequalities(matrix([(1,1,3), (5,5,7)]))
            [1 3]
            [5 7]
            sage: piv._pivots
            (0, 2)
            sage: piv._linear_subspace
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -1  0]
        """
        self.base_ring = base_ring
        self.dim = dim

    def _pivot_inequalities(self, A):
        """
        Pick pivots for inequalities.

        INPUT:

        - ``A`` -- matrix. The inequalities.

        OUTPUT:

        The matrix of pivot columns.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous \
            ....:     import PivotedInequalities
            sage: piv = PivotedInequalities(QQ, 2)
            sage: piv._pivot_inequalities(matrix([(1,1,3), (5,5,7)]))
            [1 3]
            [5 7]
        """
        self._linear_subspace = A.right_kernel()
        self._pivots = A.pivots()
        self._nonpivots = A.nonpivots()
        return A.matrix_from_columns(self._pivots)

    def _unpivot_ray(self, ray):
        """
        Undo the pivoting to go back to the original inequalities
        containing a linear subspace.

        INPUT:

        - ``ray`` -- ray in the pivoted coordinates.

        OUTPUT:

        Ray in the original coordinates.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous \
            ....:     import PivotedInequalities
            sage: piv = PivotedInequalities(QQ, 2)
            sage: piv._pivot_inequalities(matrix([(1,1,3), (5,5,7)]))
            [1 3]
            [5 7]
            sage: piv._unpivot_ray([1, 3])
            (1, 0, 3)
        """
        result = [self.base_ring.zero()] * (self.dim + 1)
        for r, i in zip(ray, self._pivots):
            result[i] = r
        return vector(self.base_ring, result)


class Hrep2Vrep(PivotedInequalities):

    def __init__(self, base_ring, dim, inequalities, equations):
        """
        Convert H-representation to a minimal V-representation.

        INPUT:

        - ``base_ring`` -- a field.

        - ``dim`` -- integer. The ambient space dimension.

        - ``inequalities`` -- list of inequalities. Each inequality
          is given as constant term, ``dim`` coefficients.

        - ``equations`` -- list of equations. Same notation as for
          inequalities.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,4,3)], [])
            [-1/2|-1/2  1/2|]
            [   0| 2/3 -1/3|]
            sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,-2,-3)], [])
            [   1 -1/2||   1]
            [   0    0||-2/3]
            sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,2,3)], [])
            [-1/2| 1/2|   1]
            [   0|   0|-2/3]
            sage: Hrep2Vrep(QQ, 2, [(8,7,-2), (1,-4,3), (4,-3,-1)], [])
            [ 1  0 -2||]
            [ 1  4 -3||]
            sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,4,3), (5,-1,-2)], [])
            [-19/5  -1/2| 2/33  1/11|]
            [ 22/5     0|-1/33 -2/33|]
            sage: Hrep2Vrep(QQ, 2, [(0,2,3), (0,4,3), (0,-1,-2)], [])
            [   0| 1/2  1/3|]
            [   0|-1/3 -1/6|]
            sage: Hrep2Vrep(QQ, 2, [], [(1,2,3), (7,8,9)])
            [-2||]
            [ 1||]
            sage: Hrep2Vrep(QQ, 2, [(1,0,0)], [])    # universe
            [0||1 0]
            [0||0 1]
            sage: Hrep2Vrep(QQ, 2, [(-1,0,0)], [])   # empty
            []
            sage: Hrep2Vrep(QQ, 2, [], [])   # universe
            [0||1 0]
            [0||0 1]
        """
        super(Hrep2Vrep, self).__init__(base_ring, dim)
        inequalities = [list(x) for x in inequalities]
        equations = [list(x) for x in equations]
        if not inequalities and not equations:
            # Adding a trivial inequality, so that the ambient dimension is passed to the algorithm.
            inequalities = [[self.base_ring.one()] + [self.base_ring.zero()] * self.dim]
        A = self._init_Vrep(inequalities, equations)
        DD = Algorithm(A).run()
        self._extract_Vrep(DD)
        if VERIFY_RESULT:
            self.verify(inequalities, equations)

    def _init_Vrep(self, inequalities, equations):
        """
        Split off the linear subspace from the inequalities and select pivots

        INPUT:

        - ``inequalities``, ``equations`` -- see :class:`Vrep2Hrep`.

        OUTPUT:

        The pivoted inequalities.

        TESTS::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H2V = Hrep2Vrep(QQ, 2, [], [])
            sage: H2V._init_Vrep([(1,0,3), (3,0,1)], [])
            [1 3]
            [3 1]
        """
        neg_eqns = [[-e for e in eqn] for eqn in equations]
        A = matrix(self.base_ring, equations + neg_eqns + inequalities)
        return self._pivot_inequalities(A)

    def _split_linear_subspace(self):
        r"""
        Split the linear subspace in a generator with `x_0\not=0` and the
        remaining generators with `x_0=0`.

        OUTPUT:

        Pair consisting of a line generator with its first coordinate
        scaled to one (if it exists, otherwise ``None``) and a list of
        remaining line generators whose first coordinate has been
        chosen to be zero.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H = Hrep2Vrep(QQ, 2, [(1,2,3)], [])
            sage: H._split_linear_subspace()
            ((1, 0, -1/3), [(0, 1, -2/3)])
            sage: H = Hrep2Vrep(QQ, 2, [(1,0,0)], [])
            sage: H._split_linear_subspace()
            (None, [(0, 1, 0), (0, 0, 1)])
        """
        lines = self._linear_subspace.basis_matrix().rows()
        L0 = []
        L1 = []
        zero = self.base_ring.zero()
        for l in lines:
            if l[0] == zero:
                L0.append(l)
            else:
                l = l / l[0]
                L1.append(l)
        if len(L1) == 0:
            return None, L0
        else:
            l1 = L1.pop()
            return l1, L0 + [l - l[0] * l1 for l in L1]

    def _extract_Vrep(self, DD):
        """
        Extract the V-representation from the extremal rays
        of the homogeneous cone.

        The V-representation is the intersection of the cone generated
        by the rays `R` and ``self._linear_subspace`` with the
        hyperplane `(1,*,*,...,*)`.

        INPUT:

        - ``DD`` -- a
          :class:`~sage.geometry.polyhedron.double_description.DoubleDescriptionPair`.

        TESTS::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H = Hrep2Vrep(QQ, 1, [(1,2)], [])
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: DD = StandardAlgorithm(matrix([[1,2], [3,5]])).run()
            sage: H._extract_Vrep(DD)
            sage: H.vertices
            [(-1/2)]
        """
        R = [self._unpivot_ray(_) for _ in DD.R]

        line1, L0 = self._split_linear_subspace()
        if line1:
            # we can shift all rays to have x_0 = 0, stay extremal
            L1 = [line1]
            R1 = []
            R0 = [r - r[0] * line1 for r in R]
        else:
            # have to really intersect with x_0 = 0
            L1 = []
            zero = self.base_ring.zero()
            R1 = [r / r[0] for r in R if r[0] > zero]
            DD0 = DD.first_coordinate_plane()
            R0 = [self._unpivot_ray(_) for _ in DD0.R]

        vertices = []
        one = self.base_ring.one()
        for v in R1 + L1:
            assert v[0] == one
            vertices.append(v[1:])
        self.vertices = vertices
        if len(vertices) > 0:
            self.rays = [r[1:] for r in R0]
            self.lines = [l[1:] for l in L0]
        else:
            # empty polyhedron
            self.rays = self.lines = []

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H = Hrep2Vrep(QQ, 1, [(1,2)], [])
            sage: H._repr_()
            '[-1/2| 1/2|]'
        """
        from sage.matrix.constructor import block_matrix

        def make_matrix(rows):
             return matrix(self.base_ring, len(rows), self.dim, rows).transpose()
        V = make_matrix(self.vertices)
        R = make_matrix(self.rays)
        L = make_matrix(self.lines)
        return str(block_matrix([[V, R, L]]))

    def verify(self, inequalities, equations):
        """
        Compare result to PPL if the base ring is QQ.

        This method is for debugging purposes and compares the
        computation with another backend if available.

        INPUT:

        - ``inequalities``, ``equations`` -- see :class:`Hrep2Vrep`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H = Hrep2Vrep(QQ, 1, [(1,2)], [])
            sage: H.verify([(1,2)], [])
        """
        from sage.rings.rational_field import QQ
        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.base_ring is not QQ:
            return
        P = Polyhedron(vertices=self.vertices, rays=self.rays, lines=self.lines,
                       base_ring=QQ, ambient_dim=self.dim, backend='ppl')
        Q = Polyhedron(ieqs=inequalities, eqns=equations,
                       base_ring=QQ, ambient_dim=self.dim, backend='ppl')
        if (P != Q) or \
           (len(self.vertices) != P.n_vertices()) or \
           (len(self.rays) != P.n_rays()) or \
           (len(self.lines) != P.n_lines()):
            print('incorrect!', end="")
            print(Q.Vrepresentation())
            print(P.Hrepresentation())


class Vrep2Hrep(PivotedInequalities):

    def __init__(self, base_ring, dim, vertices, rays, lines):
        """
        Convert V-representation to a minimal H-representation.

        INPUT:

        - ``base_ring`` -- a field.

        - ``dim`` -- integer. The ambient space dimension.

        - ``vertices`` -- list of vertices. Each vertex is given as
          list of ``dim`` coordinates.

        - ``rays`` -- list of rays. Each ray is given as
          list of ``dim`` coordinates, not all zero.

        - ``lines`` -- list of line generators. Each line is given as
          list of ``dim`` coordinates, not all zero.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
            sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
            [1 2 3]
            [2 4 3]
            [-----]

            sage: Vrep2Hrep(QQ, 2, [(1,0), (-1/2,0)], [], [(1,-2/3)])
            [ 1/3  2/3    1]
            [ 2/3 -2/3   -1]
            [--------------]

            sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [(1/2,0)], [(1,-2/3)])
            [1 2 3]
            [-----]

            sage: Vrep2Hrep(QQ, 2, [(1,1), (0,4), (-2,-3)], [], [])
            [ 8/13  7/13 -2/13]
            [ 1/13 -4/13  3/13]
            [ 4/13 -3/13 -1/13]
            [-----------------]

            sage: Vrep2Hrep(QQ, 2, [(-19/5,22/5), (-1/2,0)], [(2/33,-1/33), (1/11,-2/33)], [])
            [10/11 -2/11 -4/11]
            [ 66/5 132/5  99/5]
            [ 2/11  4/11  6/11]
            [-----------------]

            sage: Vrep2Hrep(QQ, 2, [(0,0)], [(1/2,-1/3), (1/3,-1/6)], [])
            [  0  -6 -12]
            [  0  12  18]
            [-----------]

            sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [], [(1,-2/3)])
            [-----]
            [1 2 3]

            sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [], [(1,-2/3), (1,0)])
            []
        """
        super(Vrep2Hrep, self).__init__(base_ring, dim)
        if rays or lines:
            assert len(vertices) > 0
        if not vertices and not rays and not lines:
            # The algorithm does not work, as the ambient dimension cannot be passed.
            # Manually setting a single equality in this case.
            one = self.base_ring.one()
            zero = self.base_ring.zero()
            self.equations = [[one] + [zero]*self.dim]
            self.inequalities = []
        else:
            A = self._init_Vrep(vertices, rays, lines)
            DD = Algorithm(A).run()
            self._extract_Hrep(DD)
        if VERIFY_RESULT:
            self.verify(vertices, rays, lines)

    def _init_Vrep(self, vertices, rays, lines):
        """
        Split off the linear subspace from the inequalities and select pivots.

        INPUT:

        - ``vertices``, ``rays``, ``lines`` -- see :class:`Vrep2Hrep`.

        OUTPUT:

        Matrix of pivoted inequalities for the dual homogenized cone.

        TESTS::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
            sage: V2H = Vrep2Hrep(QQ, 2, [(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
            sage: V2H._init_Vrep([(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
            [   1 -1/2    0]
            [   0 -1/2  2/3]
            [   0  1/2 -1/3]
        """
        one = self.base_ring.one()
        zero = self.base_ring.zero()
        homogeneous = \
            [[one] + list(v) for v in vertices] + \
            [[zero] + list(r) for r in rays] + \
            [[zero] + list(l) for l in lines] + \
            [[zero] + [-x for x in l] for l in lines]
        A = matrix(self.base_ring, homogeneous)
        return self._pivot_inequalities(A)

    def _extract_Hrep(self, DD):
        """
        Extract generators from the dual description of the homogenized cone.

        INPUT:

        - ``DD`` -- a
          :class:`~sage.geometry.polyhedron.double_description.DoubleDescriptionPair`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
            sage: V2H = Vrep2Hrep(QQ, 1, [(-1/2,), (2/3)], [], [])
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: DD = StandardAlgorithm(matrix([[1,2], [3,5]])).run()
            sage: V2H._extract_Hrep(DD)
        """
        zero = self.base_ring.zero()
        def is_trivial(ray):
            # trivial Hrep output 1 >= 0
            return ray[0] > zero and all(r == zero for r in ray[1:])
        ieqs = [self._unpivot_ray(_) for _ in DD.R]
        self.inequalities = [r for r in ieqs if not is_trivial(r)]
        self.equations = self._linear_subspace.matrix().rows()

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
            sage: V2H = Vrep2Hrep(QQ, 2, [(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
            sage: V2H._repr_()
            '[1 2 3]\n[2 4 3]\n[-----]'
        """
        from sage.matrix.constructor import block_matrix

        def make_matrix(cols):
            return matrix(self.base_ring, len(cols), self.dim + 1, cols)
        I = make_matrix(self.inequalities)
        E = make_matrix(self.equations)
        return str(block_matrix([[I], [E]]))

    def verify(self, vertices, rays, lines):
        """
        Compare result to PPL if the base ring is QQ.

        This method is for debugging purposes and compares the
        computation with another backend if available.

        INPUT:

        - ``vertices``, ``rays``, ``lines`` -- see :class:`Vrep2Hrep`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
            sage: vertices = [(-1/2,0)]
            sage: rays = [(-1/2,2/3), (1/2,-1/3)]
            sage: lines = []
            sage: V2H = Vrep2Hrep(QQ, 2, vertices, rays, lines)
            sage: V2H.verify(vertices, rays, lines)
        """
        from sage.rings.rational_field import QQ
        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.base_ring is not QQ:
            return
        P = Polyhedron(vertices=vertices, rays=rays, lines=lines,
                       base_ring=QQ, ambient_dim=self.dim)
        Q = Polyhedron(ieqs=self.inequalities, eqns=self.equations,
                       base_ring=QQ, ambient_dim=self.dim)
        if not P == Q:
            print('incorrect!', P, Q)
            print(Q.Vrepresentation())
            print(P.Hrepresentation())
