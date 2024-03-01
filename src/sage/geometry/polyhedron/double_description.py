r"""
Double Description Algorithm for Cones

This module implements the double description algorithm for extremal
vertex enumeration in a pointed cone following [FP1996]_. With a
little bit of preprocessing (see
:mod:`~sage.geometry.polyhedron.double_description_inhomogeneous`)
this defines a backend for polyhedral computations. But as far as this
module is concerned, *inequality* always means without a constant term
and the origin is always a point of the cone.

EXAMPLES::

    sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
    sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
    sage: alg = StandardAlgorithm(A);  alg
    Pointed cone with inequalities
    (1, 0, 1)
    (0, 1, 1)
    (-1, -1, 1)
    sage: DD, _ = alg.initial_pair();  DD
    Double description pair (A, R) defined by
        [ 1  0  1]        [ 2/3 -1/3 -1/3]
    A = [ 0  1  1],   R = [-1/3  2/3 -1/3]
        [-1 -1  1]        [ 1/3  1/3  1/3]

The implementation works over any exact field that is embedded in
`\RR`, for example::

    sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
    sage: A = matrix(AA, [(1,0,1), (0,1,1), (-AA(2).sqrt(),-AA(3).sqrt(),1),            # needs sage.rings.number_field
    ....:                 (-AA(3).sqrt(),-AA(2).sqrt(),1)])
    sage: alg = StandardAlgorithm(A)
    sage: alg.run().R                                                                   # needs sage.rings.number_field
    [(-0.4177376677004119?, 0.5822623322995881?, 0.4177376677004119?),
     (-0.2411809548974793?, -0.2411809548974793?, 0.2411809548974793?),
     (0.07665629029830300?, 0.07665629029830300?, 0.2411809548974793?),
     (0.5822623322995881?, -0.4177376677004119?, 0.4177376677004119?)]
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                     2015 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


#*****************************************************************************
# TODO
#
# The adjacency check should use caching and the "combinatorial
# criterion" instead of the "algebraic criterion", see [FP1996]
# for definition. Since coefficient arithmetic is relatively expensive
# we should avoid it as far as possible.
#
# Also, the variants of the double description algorithm described in
# [FP1996] should be implemented. The design of this module is
# such that variants of the basic algorithm should be easy to add as
# subclasses of DoubleDescriptionPair and Problem.
# *****************************************************************************


# Compare with PPL if the base ring is QQ. Can be left enabled since
# we don't use the Python fallback for polyhedra over QQ unless you
# construct one by hand.

VERIFY_RESULT = True

import itertools

from sage.misc.cachefunc import cached_method
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.matrix.matrix_space import MatrixSpace


def random_inequalities(d, n):
    """
    Random collections of inequalities for testing purposes.

    INPUT:

    - ``d`` -- integer. The dimension.

    - ``n``  -- integer. The number of random inequalities to generate.

    OUTPUT:

    A random set of inequalities as a :class:`StandardAlgorithm` instance.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.double_description import random_inequalities
        sage: P = random_inequalities(5, 10)
        sage: P.run().verify()
    """
    from sage.matrix.constructor import random_matrix
    while True:
        A = random_matrix(QQ, n, d)
        if A.rank() == min(n, d) and not any(a == 0 for a in A.rows()):
            break
    return StandardAlgorithm(A)


class DoubleDescriptionPair:

    def __init__(self, problem, A_rows, R_cols):
        r"""
        Base class for a double description pair `(A, R)`

        .. warning::

            You should use the :meth:`Problem.initial_pair` or
            :meth:`Problem.run` to generate double description pairs
            for a set of inequalities, and not generate
            ``DoubleDescriptionPair`` instances directly.

        INPUT:

        - ``problem`` -- instance of :class:`Problem`.

        - ``A_rows`` -- list of row vectors of the matrix `A`. These
          encode the inequalities.

        - ``R_cols`` -- list of column vectors of the matrix
          `R`. These encode the rays.

        TESTS::

            sage: from sage.geometry.polyhedron.double_description import \
            ....:     DoubleDescriptionPair, Problem
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: alg = Problem(A)
            sage: DoubleDescriptionPair(alg,
            ....:     [(1, 0, 1), (0, 1, 1), (-1, -1, 1)],
            ....:     [(2/3, -1/3, 1/3), (-1/3, 2/3, 1/3), (-1/3, -1/3, 1/3)])
            Double description pair (A, R) defined by
                [ 1  0  1]        [ 2/3 -1/3 -1/3]
            A = [ 0  1  1],   R = [-1/3  2/3 -1/3]
                [-1 -1  1]        [ 1/3  1/3  1/3]
        """
        self.problem = problem
        self.A = list(A_rows)
        self.R = list(R_cols)
        self.one = problem._field.one()
        self.zero = problem._field.zero()

        # a cache for scalar products (see the method zero_set)
        self.zero_set_cache = {}

    def _make_new(self, A_rows, R_cols):
        r"""
        Construct a new double description pair.

        INPUT:

        - ``A_rows`` -- list of row vectors of the matrix `A`. These
          encode the inequalities.

        - ``R_cols`` -- list of column vectors of the matrix
          `R`. These encode the rays.

        OUTPUT:

        A new double description pair of the same (sub)class of
        :class:`DoubleDescriptionProblem`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import \
            ....:     DoubleDescriptionPair, StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DDnew = DD._make_new(DD.A, DD.R);  DDnew
            Double description pair (A, R) defined by
                [ 1  0  1]        [ 2/3 -1/3 -1/3]
            A = [ 0  1  1],   R = [-1/3  2/3 -1/3]
                [-1 -1  1]        [ 1/3  1/3  1/3]
            sage: DDnew is DD
            False
            sage: DDnew.__class__ is DD.__class__
            True
        """
        return self.__class__(self.problem, A_rows, R_cols)

    def __repr__(self):
        r"""
        Return string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import \
            ....:     DoubleDescriptionPair, StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DD.__repr__()
            'Double description pair (A, R) defined by\n    [ 1  0  1]
             [ 2/3 -1/3 -1/3]\nA = [ 0  1  1],   R = [-1/3  2/3 -1/3]\n
             [-1 -1  1]        [ 1/3  1/3  1/3]'
        """
        from sage.typeset.ascii_art import ascii_art
        from sage.matrix.constructor import matrix
        s = ascii_art('Double description pair (A, R) defined by')
        A = ascii_art(matrix(self.A))
        A._baseline = (len(self.A) // 2)
        A = ascii_art('A = ') + A
        R = ascii_art(matrix(self.R).transpose())
        if len(self.R) > 0:
            R._baseline = (len(self.R[0]) // 2)
        else:
            R._baseline = 0
        R = ascii_art('R = ') + R
        return str(s * (A + ascii_art(',   ') + R))

    def inner_product_matrix(self):
        """
        Return the inner product matrix between the rows of `A`
        and the columns of `R`.

        OUTPUT:

        A matrix over the base ring. There is one row for each row of
        `A` and one column for each column of `R`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: alg = StandardAlgorithm(A)
            sage: DD, _ = alg.initial_pair()
            sage: DD.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        from sage.matrix.constructor import matrix
        return matrix(self.problem.base_ring(), [[a.inner_product(r) for r in self.R] for a in self.A])

    def cone(self):
        r"""
        Return the cone defined by `A`.

        This method is for debugging only. Assumes that the base ring
        is `\QQ`.

        OUTPUT:

        The cone defined by the inequalities as a
        :func:`~sage.geometry.polyhedron.constructor.Polyhedron`,
        using the PPL backend.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD.cone().Hrepresentation()
            (An inequality (-1, -1, 1) x + 0 >= 0,
             An inequality (0, 1, 1) x + 0 >= 0,
             An inequality (1, 0, 1) x + 0 >= 0)
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        assert self.problem.base_ring() == QQ    # required for PPL backend

        if not self.A:
            return Polyhedron(vertices=[[0] * self.problem.dim()], backend='ppl')
        else:
            ieqs = [[0] + list(a) for a in self.A]
            return Polyhedron(ieqs=ieqs, base_ring=self.problem.base_ring(), backend='ppl')

    def verify(self):
        r"""
        Validate the double description pair.

        This method used the PPL backend to check that the double
        description pair is valid. An assertion is triggered if it is
        not. Does nothing if the base ring is not `\QQ`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import \
            ....:     DoubleDescriptionPair, Problem
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: alg = Problem(A)
            sage: DD = DoubleDescriptionPair(alg,
            ....:     [(1, 0, 3), (0, 1, 1), (-1, -1, 1)],
            ....:     [(2/3, -1/3, 1/3), (-1/3, 2/3, 1/3), (-1/3, -1/3, 1/3)])
            sage: DD.verify()
            Traceback (most recent call last):
            ...
                assert A_cone == R_cone
            AssertionError
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.problem.base_ring() is not QQ:
            return
        A_cone = self.cone()
        R_cone = Polyhedron(vertices=[[self.zero] * self.problem.dim()], rays=self.R,
                            base_ring=self.problem.base_ring(), backend='ppl')
        assert A_cone == R_cone
        assert A_cone.n_inequalities() <= len(self.A)
        assert R_cone.n_rays() == len(self.R)

    def R_by_sign(self, a):
        """
        Classify the rays into those that are positive, zero, and negative on `a`.

        INPUT:

        - ``a`` -- vector. Coefficient vector of a homogeneous inequality.

        OUTPUT:

        A triple consisting of the rays (columns of `R`) that are
        positive, zero, and negative on `a`. In that order.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD.R_by_sign(vector([1,-1,0]))
            ([(2/3, -1/3, 1/3)], [(-1/3, -1/3, 1/3)], [(-1/3, 2/3, 1/3)])
            sage: DD.R_by_sign(vector([1,1,1]))
            ([(2/3, -1/3, 1/3), (-1/3, 2/3, 1/3)], [], [(-1/3, -1/3, 1/3)])
        """
        pos = []
        nul = []
        neg = []
        for r in self.R:
            sgn = a * r
            if sgn == self.zero:
                nul.append(r)
            elif sgn > self.zero:
                pos.append(r)
            else:
                neg.append(r)
        return pos, nul, neg

    def zero_set(self, ray):
        """
        Return the zero set (active set) `Z(r)`.

        INPUT:

        - ``ray`` -- a ray vector.

        OUTPUT:

        A set containing the inequality vectors that are zero on ``ray``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = Problem(A).initial_pair()
            sage: r = DD.R[0];  r
            (2/3, -1/3, 1/3)
            sage: DD.zero_set(r)
            {(-1, -1, 1), (0, 1, 1)}
        """
        if ray not in self.zero_set_cache:
            self.zero_set_cache[ray] = (0, set())
        n, t = self.zero_set_cache[ray]
        if n != len(self.A):
            t.update(self.A[i] for i in range(n,len(self.A)) if self.A[i].inner_product(ray) == self.zero)
            self.zero_set_cache[ray] = (len(self.A), t)
        return t

    def is_extremal(self, ray):
        """
        Test whether the ray is extremal.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DD.is_extremal(DD.R[0])
            True
        """
        from sage.matrix.constructor import matrix
        A_Zray = matrix(self.problem.base_ring(), list(self.zero_set(ray)))
        return A_Zray.rank() == self.problem.dim() - 1

    @cached_method
    def matrix_space(self, nrows, ncols):
        r"""
        Return a matrix space of size ``nrows`` and ``ncols`` over the base ring
        of ``self``.

        These matrix spaces are cached to avoid their creation in the very
        demanding :meth:`add_inequality` and more precisely :meth:`are_adjacent`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = Problem(A).initial_pair()
            sage: DD.matrix_space(2,2)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: DD.matrix_space(3,2)
            Full MatrixSpace of 3 by 2 dense matrices over Rational Field

            sage: # needs sage.rings.number_field
            sage: K.<sqrt2> = QuadraticField(2)
            sage: A = matrix([[1,sqrt2],[2,0]])
            sage: DD, _  = Problem(A).initial_pair()
            sage: DD.matrix_space(1,2)
            Full MatrixSpace of 1 by 2 dense matrices
             over Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
        """
        return MatrixSpace(self.problem.base_ring(), nrows, ncols)

    def are_adjacent(self, r1, r2):
        """
        Return whether the two rays are adjacent.

        INPUT:

        - ``r1``, ``r2`` -- two rays.

        OUTPUT:

        Boolean. Whether the two rays are adjacent.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DD.are_adjacent(DD.R[0], DD.R[1])
            True
            sage: DD.are_adjacent(DD.R[0], DD.R[2])
            True
            sage: DD.are_adjacent(DD.R[0], DD.R[3])
            False
        """
        Z = self.zero_set(r1).intersection(self.zero_set(r2))
        if not Z:
            return self.problem.dim() == 2
        Z = list(Z)

        # here we try to create a matrix as fast as possible
        # since the generic matrix constructor is very slow (trac #18231)
        A_Z12 = self.matrix_space(len(Z), len(Z[0])).matrix(Z, coerce=False)
        return A_Z12.rank() == self.problem.dim() - 2

    def dual(self):
        """
        Return the dual.

        OUTPUT:

        For the double description pair `(A, R)` this method returns
        the dual double description pair `(R^T, A^T)`

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD, _ = Problem(A).initial_pair()
            sage: DD
            Double description pair (A, R) defined by
                [ 0  1  0]        [0 1 0]
            A = [ 1  0  0],   R = [1 0 0]
                [ 0 -1  1]        [1 0 1]
            sage: DD.dual()
             Double description pair (A, R) defined by
                [0 1 1]        [ 0  1  0]
            A = [1 0 0],   R = [ 1  0 -1]
                [0 0 1]        [ 0  0  1]
        """
        return self._make_new(self.R, self.A)

    def first_coordinate_plane(self):
        """
        Restrict to the first coordinate plane.

        OUTPUT:

        A new double description pair with the constraint `x_0 = 0`
        added.

        EXAMPLES::

            sage: A = matrix([(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD
            Double description pair (A, R) defined by
            A = [ 1  1],   R = [ 1/2 -1/2]
                [-1  1]        [ 1/2  1/2]
            sage: DD.first_coordinate_plane()
            Double description pair (A, R) defined by
                [ 1  1]
            A = [-1  1],   R = [  0]
                [-1  0]        [1/2]
                [ 1  0]
        """
        R = self.problem.base_ring()
        d = self.problem.dim()
        a_neg = vector(R, [-self.one] + [self.zero] * (d - 1))
        a_pos = vector(R, [+self.one] + [self.zero] * (d - 1))
        new = self._make_new(self.A, self.R)
        new.add_inequality(a_neg)
        new.add_inequality(a_pos)
        return new


class Problem:

    pair_class = DoubleDescriptionPair

    def __init__(self, A):
        r"""
        Base class for implementations of the double description algorithm

        It does not make sense to instantiate the base class directly,
        it just provides helpers for implementations.

        INPUT:

        - ``A`` -- a matrix. The rows of the matrix are interpreted as
          homogeneous inequalities `A x \geq 0`. Must have maximal rank.

        TESTS::

            sage: A = matrix([(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A)
            Pointed cone with inequalities
            (1, 1)
            (-1, 1)
        """
        assert A.rank() == A.ncols()    # implementation assumes maximal rank
        if A.is_mutable():
            A = A.__copy__()
            A.set_immutable()
        self._A = A
        self._field = A.base_ring().fraction_field()

    @cached_method
    def A(self):
        """
        Return the rows of the defining matrix `A`.

        OUTPUT:

        The matrix `A` whose rows are the inequalities.

        EXAMPLES::

            sage: A = matrix([(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A).A()
            ((1, 1), (-1, 1))
        """
        rows = [a.change_ring(self._field) for a in self._A.rows()]
        for a in rows:
            a.set_immutable()
        return tuple(rows)

    def A_matrix(self):
        """
        Return the defining matrix `A`.

        OUTPUT:

        Matrix whose rows are the inequalities.

        EXAMPLES::

            sage: A = matrix([(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A).A_matrix()
            [ 1  1]
            [-1  1]
        """
        return self._A

    def base_ring(self):
        """
        Return the base field.

        OUTPUT:

        A field.

        EXAMPLES::

            sage: A = matrix(AA, [(1, 1), (-1, 1)])                                     # needs sage.rings.number_field
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A).base_ring()                                                # needs sage.rings.number_field
            Algebraic Real Field
        """
        return self._field

    @cached_method
    def dim(self):
        """
        Return the ambient space dimension.

        OUTPUT:

        Integer. The ambient space dimension of the cone.

        EXAMPLES::

            sage: A = matrix(QQ, [(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A).dim()
            2
        """
        return self._A.ncols()

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: A = matrix(QQ, [(1, 1), (-1, 1)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: Problem(A).__repr__()
            'Pointed cone with inequalities\n(1, 1)\n(-1, 1)'
        """
        return 'Pointed cone with inequalities\n' + '\n'.join(map(str, self.A()))

    def initial_pair(self):
        """
        Return an initial double description pair.

        Picks an initial set of rays by selecting a basis. This is
        probably the most efficient way to select the initial set.

        INPUT:

        - ``pair_class`` -- subclass of
          :class:`DoubleDescriptionPair`.

        OUTPUT:

        A pair consisting of a :class:`DoubleDescriptionPair` instance
        and the tuple of remaining unused inequalities.

        EXAMPLES::

            sage: A = matrix([(-1, 1), (-1, 2), (1/2, -1/2), (1/2, 2)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: DD, remaining = Problem(A).initial_pair()
            sage: DD.verify()
            sage: remaining
            [(1/2, -1/2), (1/2, 2)]
        """
        pivot_rows = self.A_matrix().pivot_rows()
        A0 = [self.A()[pivot] for pivot in pivot_rows]
        Ac = [self.A()[i] for i in range(len(self.A())) if i not in pivot_rows]
        from sage.matrix.constructor import identity_matrix, matrix
        I = identity_matrix(self.base_ring(), self.dim())
        R = matrix(self.base_ring(), A0).solve_right(I)
        return self.pair_class(self, A0, R.columns()), list(Ac)


class StandardDoubleDescriptionPair(DoubleDescriptionPair):
    """
    Double description pair for the "Standard Algorithm".

    See :class:`StandardAlgorithm`.

    TESTS::

        sage: A = matrix([(-1, 1, 0), (-1, 2, 1), (1/2, -1/2, -1)])
        sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
        sage: DD, _ = StandardAlgorithm(A).initial_pair()
    """

    def add_inequality(self, a):
        """
        Add the inequality ``a`` to the matrix `A` of the double description.

        INPUT:

        - ``a`` -- vector. An inequality.

        EXAMPLES::

            sage: A = matrix([(-1, 1, 0), (-1, 2, 1), (1/2, -1/2, -1)])
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD.add_inequality(vector([1,0,0]))
            sage: DD
            Double description pair (A, R) defined by
                [  -1    1    0]        [   1    1    0    0]
            A = [  -1    2    1],   R = [   1    1    1    1]
                [ 1/2 -1/2   -1]        [   0   -1 -1/2   -2]
                [   1    0    0]
        """
        R_pos, R_nul, R_neg = self.R_by_sign(a)
        if not R_neg:
            return
        R_new = []
        for rp, rn in itertools.product(R_pos, R_neg):
            if not self.are_adjacent(rp, rn):
                continue
            r = a.inner_product(rp) * rn - a.inner_product(rn) * rp
            r.set_immutable()
            R_new.append(r)
        self.R = R_pos + R_nul + R_new
        self.A.append(a)

class StandardAlgorithm(Problem):
    """
    Standard implementation of the double description algorithm

    See [FP1996]_ for the definition of the "Standard
    Algorithm".

    EXAMPLES::

        sage: A = matrix(QQ, [(1, 1), (-1, 1)])
        sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
        sage: DD = StandardAlgorithm(A).run()
        sage: DD.R    # the extremal rays
        [(1/2, 1/2), (-1/2, 1/2)]
    """
    pair_class = StandardDoubleDescriptionPair

    def run(self):
        """
        Run the Standard Algorithm.

        OUTPUT:

        A double description pair `(A, R)` of all inequalities as a
        :class:`DoubleDescriptionPair`.  By virtue of the double
        description algorithm, the columns of `R` are the extremal
        rays.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: StandardAlgorithm(A).run()
            Double description pair (A, R) defined by
                [ 0  1  0]        [0 0 1 1]
            A = [ 1  0  0],   R = [1 0 1 0]
                [ 0 -1  1]        [1 1 1 1]
                [-1  0  1]
        """
        DD, remaining = self.initial_pair()
        for a in remaining:
            DD.add_inequality(a)
            if VERIFY_RESULT:
                DD.verify()
        return DD
