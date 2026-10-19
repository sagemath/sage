r"""
Determinantal Point Processes

This models determinant point processes (DPP) by using the sampling method
from Poulson [Pou2020]_.

AUTHORS:

- Travis Scrimshaw (2025-08): initial version
"""

# ****************************************************************************
#       Copyright (C) 2025 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.prandom import random
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ

class DeterminantalPointProcess(SageObject):
    r"""
    A determinantal point process (DPP).

    INPUT:

    - ``kernel`` -- function `X \times X \to [0, 1]`, where `X` is
      some set that will be used in the :meth:`sample`

    EXAMPLES::

        sage: DPP = DeterminantalPointProcess(lambda x,y: 0.02 * exp(-(y - x).norm()^2/0.018))
        sage: domain = [vector([RR(x/10), RR(y/10)]) for x in range(11) for y in range(11)]
        sage: X = DPP.sample(domain)
        sage: scatter_plot([domain[i] for i in X], xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1)
        Launched png viewer for Graphics object consisting of 1 graphics primitive
    """
    def __init__(self, kernel):
        self._kernel = kernel

    def sample(self, domain, algorithm=None):
        return self._unblocked_rl_nonhermitian(domain)[0]

    def kernel_matrix(self, domain):
        return matrix([[self._kernel(i, j) for j in domain] for i in domain])

    def _unblocked_rl_nonhermitian(self, domain):
        import numpy as np
        domain = tuple(domain)
        n = len(domain)
        sample = []
        K = np.array([[self._kernel(i, j) for j in domain] for i in domain])
        # TODO: Switch to numpy
        for j in range(n):
            if random() < K[j,j]:
                sample.append(j)
            else:
                K[j, j] -= 1
            K[j+1:, j] /= K[j, j]
            K[j+1:, j+1:] -= np.outer(K[j+1:, j], K[j, j+1:])
            #for i in range(j+1, n):
            #    K[i, j] /= val
            #temp = K[j+1:, j:j+1] * K[j:j+1, j+1:]
            #assert temp.nrows() == n - (j+1)
            #for i in range(j+1, n):
            #    for k in range(j+1, n):
            #        assert i-j-1 >= 0 and k-j-1 >= 0
            #        K[i, k] -= temp[i-j-1, k-j-1]
        return (sample, K)

    def check_kernel_matrix(self, K):
        if not K.is_square():
            return False
        n = K.nrows()
        if all(0 <= K[i,i] <= 1
               and min(K[i,i], 1-K[i,i]) > sum(abs(K[i,j]) for j in range(n) if j != i)
               for i in range(n)):
            return True
        half = QQ((1, 2))
        if (K - matrix.identity(n) * half).norm(p=2) <= half:
            return True
        if K.is_symmetric():
            return False
        return self._brunel_kernel_test(K)

    def _brunel_kernel_test(self, K, tol=0.000001):
        """
        Really slow fallback test.

        EXAMPLES::

            sage: DPP = DeterminantalPointProcess(lambda x,y: 0.02 * exp(-(y - x).norm()^2/0.018))
            sage: domain = [vector([RR(x/3), RR(y/3)]) for x in range(3) for y in range(3)]
            sage: K = DPP.kernel_matrix(domain)
            sage: DPP._brunel_kernel_test(K)
            True
        """
        if not K.is_square():
            return False
        # Run the test from Burnel (see Prop. 1 of [Pou2022]_)
        from itertools import combinations
        from copy import copy
        n = K.nrows()
        for N in range(1, n+1):
            for J in combinations(range(n), N):
                Kp = copy(K)
                for j in J:
                    Kp[j, j] -= 1
                if (-1)**N * Kp.det() < -tol:
                    return False, (J, (-1)**N * Kp.det())
        return True

class AztecDiamond(SageObject):
    """
    We follow *Asymptotic domino statistics in the Aztec diamond*
    by S. Chhita, K. Johansson, and B. Young :arXiv:`1212.5414`.
    """
    def __init__(self, n, a=QQ.one()):
        self._n = n
        self._a = a
        self._black_vertices = [(x1, x2) for x1 in range(0,2*n+1,2) for x2 in range(1,2*n,2)]
        self._white_vertices = [(x1, x2) for x1 in range(1,2*n,2) for x2 in range(0,2*n+1,2)]
        self._dpp = DeterminantalPointProcess(self.kernel)

    @cached_method
    def kasteleyn_matrix(self):
        from sage.rings.number_field.number_field import CyclotomicField
        Z4 = CyclotomicField(4)
        I = Z4.gen()
        n = self._n
        B = self._black_vertices
        W = self._white_vertices

        def entry(b, w):
            # e_1 = (1,1), e_2 = (-1,1)
            if w[0] == b[0] + 1 and w[1] == b[1] + 1: # w = b + (-1)^0 e_1
                return Z4((-1)**((b[0]+b[1]-1)//2))
            elif w[0] == b[0] - 1 and w[1] == b[1] - 1: # w = b + (-1)^1 e_1
                return Z4((-1)**(1+(b[0]+b[1]-1)//2))
            elif w[0] == b[0] + 1 and w[1] == b[1] - 1: # w = b 0 (-1)^0 e_2
                return Z4((-1)**((b[0]+b[1]-1)//2) * self._a * I)
            elif w[0] == b[0] - 1 and w[1] == b[1] + 1: # w = b - (-1)^1 e_2
                return Z4((-1)**(1+(b[0]+b[1]-1)//2) * self._a * I)
            return Z4.zero()

        return matrix([[entry(b, w) for w in W] for b in B])

    @cached_method
    def inverse_kasteleyn_matrix(self):
        return self.kasteleyn_matrix().inverse()

    @cached_method
    def domain(self):
        from itertools import product
        return tuple(product(self._black_vertices, self._white_vertices))

    def kernel(self, bw, wb):
        r1 = self._black_vertices.index(bw[0])
        c1 = self._white_vertices.index(bw[1])
        r2 = self._white_vertices.index(wb[1])
        c2 = self._black_vertices.index(bw[0])
        K = self.kasteleyn_matrix()
        Ki = self.inverse_kasteleyn_matrix()
        return K[r1,c1] * Ki[r2,c2]

    def kernel_matrix(self):
        return self._dpp.kernel_matrix(self.domain())

    def sample(self):
        X = self._dpp.sample(self.domain())
        D = self.domain()
        return [D[i] for i in X]

    def sample_plot(self):
        from sage.plot.line import line2d
        X = self.sample()
        n = self._n
        P  = [line2d([(0,y), (2*n-y,2*n)], color='blue') for y in range(1,2*n,2)]
        P += [line2d([(x,0), (2*n,2*n-x)], color='blue') for x in range(1,2*n,2)]
        P += [line2d([(0,y), (y,0)], color='blue') for y in range(1,2*n,2)]
        P += [line2d([(x,2*n), (2*n,x)], color='blue') for x in range(1,2*n,2)]
        for p1, p2 in X:
            P.append(line2d([p1, p2], color='red', thickness=4))
        P = sum(P)
        P.set_aspect_ratio(1)
        return P
