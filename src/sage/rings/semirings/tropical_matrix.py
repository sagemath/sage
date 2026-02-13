r"""
Matrices over tropical semirings

AUTHORS:

- Xavier Caruso (2025-11): initial version
"""

# ****************************************************************************
#       Copyright (C) 2025 Xavier Caruso <xavier@caruso.ovh>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.rings.real_mpfr import RR
from sage.rings.infinity import infinity


class Matrix_tropical_dense(Matrix_generic_dense):
    r"""
    A class for dense matrices over a tropical semiring.

    EXAMPLES::

        sage: from sage.rings.semirings.tropical_matrix import Matrix_tropical_dense
        sage: T = TropicalSemiring(QQ)
        sage: M = matrix(T, [[1, 2], [3, 4]])
        sage: isinstance(M, Matrix_tropical_dense)
        True
    """
    def extremum_cycle_mean(self):
        r"""
        Return the extremal (that is, minimal if the addition is max
        and maximal is the addition is min) mean weight of this matrix
        It is also the smallest/largest eigenvalue of this matrix.

        ALGORITHM:

        We implement Karp's algorithm described in [But2010]_, Section 1.6.1.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: M = matrix(T, [[-2,  1, -3],
            ....:                [ 3,  0,  3],
            ....:                [ 5,  2,  1]])
            sage: M.extremum_cycle_mean()
            3

        ::

            sage: T = TropicalSemiring(QQ)
            sage: z = T.zero()
            sage: M = matrix(T, [[z, 1, 10, z],
            ....:                [z, z,  3, z],
            ....:                [z, z,  z, 2],
            ....:                [8, 0,  z, z]])
            sage: M.extremum_cycle_mean()
            5/3

        TESTS::

            sage: M = matrix(T, [])
            sage: M.extremum_cycle_mean()
            +infinity
            sage: M = matrix(T, [[]])
            sage: M.extremum_cycle_mean()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square

        ::

            sage: R.<x> = QQ[]
            sage: T = TropicalSemiring(R)
            sage: M = matrix(T, 2, 2, [1, x, x^2, x^3])
            sage: M.extremum_cycle_mean()
            Traceback (most recent call last):
            ...
            NotImplementedError: extremum cycle mean is only implemented for subrings of RR
        """
        T = self.base_ring()
        if not T._base.is_subring(RR):
            raise NotImplementedError("extremum cycle mean is only implemented for subrings of RR")
        n = self.ncols()
        if self.nrows() != n:
            raise TypeError("matrix must be square")
        if self.is_zero():
            return T.zero()
        v = matrix(1, n, n*[T.one()])
        vs = [v]
        for _ in range(n):
            v = v * self
            vs.append(v)
        w = [vs[n][0,j].lift() for j in range(n)]
        if T._use_min:
            f, fp = max, min
        else:
            f, fp = min, max
        ans = fp(f((w[j] - vs[k][0,j].lift()) / (n-k) for k in range(n))
                 for j in range(n) if w[j] is not infinity)
        return T(ans)

    def weak_transitive_closure(self):
        r"""
        Return the weak transitive closure of this matrix `M`,
        that is, by definition

        .. MATH::

            A \oplus A^2 \oplus A^3 \oplus A^4 \oplus \cdots

        or raise an error if this sum does not converge.

        ALGORITHM:

        We implement the Floyd-Warshall algorithm described in
        [But2010]_, Algorithm 1.6.21.

        .. SEEALSO::

            :meth:`strong_transitive_closure`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: z = T.zero()
            sage: M = matrix(T, [[z, 1, 10, z],
            ....:                [z, z,  3, z],
            ....:                [z, z,  z, 2],
            ....:                [8, 0,  z, z]])
            sage: M.weak_transitive_closure()
            [14  1  4  6]
            [13  5  3  5]
            [10  2  5  2]
            [ 8  0  3  5]

        We check that the minimal cycle mean of `M` is the largest
        value `a` such that `(-a) \otimes M` has a weak transitive
        closure::

            sage: M.extremum_cycle_mean()
            5/3
            sage: aM = T(-5/3) * M
            sage: aM.weak_transitive_closure()
            [22/3 -2/3  2/3    1]
            [   8    0  4/3  5/3]
            [20/3 -4/3    0  1/3]
            [19/3 -5/3 -1/3    0]
            sage: bM = T(-2) * M
            sage: bM.weak_transitive_closure()
            Traceback (most recent call last):
            ...
            ValueError: negative cycle exists

        TESTS::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: M = matrix(T, [[ 1,  2],
            ....:                [-2, -3]])
            sage: M.weak_transitive_closure()
            Traceback (most recent call last):
            ...
            ValueError: positive cycle exists
        """
        T = self.base_ring()
        if not T._base.is_subring(RR):
            raise NotImplementedError("extremum cycle mean is only implemented for subrings of RR")
        n = self.ncols()
        if self.nrows() != n:
            raise TypeError("matrix must be square")
        G = self.__copy__()
        for p in range(n):
            for i in range(n):
                if i == p:
                    continue
                for j in range(n):
                    if j == p:
                        continue
                    G[i,j] += G[i,p] * G[p,j]
                    if i == j:
                        if T._use_min and G[i,i].lift() < 0:
                            raise ValueError("negative cycle exists")
                        if not T._use_min and G[i,i].lift() > 0:
                            raise ValueError("positive cycle exists")
        return G

    def strong_transitive_closure(self):
        r"""
        Return the string transitive closure of this matrix `M`,
        that is, by definition

        .. MATH::

            I \oplus A \oplus A^2 \oplus A^3 \oplus A^4 \oplus \cdots

        or raise an error if this sum does not converge.

        ALGORITHM:

        We implement the Floyd-Warshall algorithm described in
        [But2010]_, Algorithm 1.6.21.

        .. SEEALSO::

            :meth:`weak_transitive_closure`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: M = matrix(T, [[-5, -2, -6],
            ....:                [ 0, -3,  0],
            ....:                [ 2, -1, -2]])
            sage: M.strong_transitive_closure()
            [ 0 -2 -2]
            [ 2  0  0]
            [ 2  0  0]

        ::

            sage: T = TropicalSemiring(QQ)
            sage: M = matrix(T, [[-5, -2, -6],
            ....:                [ 0, -3,  0],
            ....:                [ 2, -1, -2]])
            sage: M.strong_transitive_closure()
            Traceback (most recent call last):
            ...
            ValueError: negative cycle exists
        """
        return self.parent().identity_matrix() + self.weak_transitive_closure()

    kleene_star = strong_transitive_closure
