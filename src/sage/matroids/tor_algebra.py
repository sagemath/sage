r"""
Tor algebra of matroids

AUTHORS:

- Travis Scrimshaw (2024-12): initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2024 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule


class TorAlgebra(CombinatorialFreeModule):
    r"""
    The Tor algebra of a matroid.

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring

    REFERENCES:

    - [Binder2024]_

    EXAMPLES::

        sage: M = matroids.catalog.P8pp()
        sage: Tor = M.tor_algebra(QQ)
        sage: Tor
    """
    def __init__(self, R, M):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Tor = matroids.Wheel(3).tor_algebra(QQ)
            sage: TestSuite(Tor).run()
        """
        from itertools import combinations
        from sage.homology.chain_complex import ChainComplex
        from sage.matrix.constructor import matrix

        self._matroid = M
        E = tuple(M.groundset())
        ordering = {e: i for i, e in enumerate(E)}  # used in case E has incomparable elements
        n = E[-1]
        bases = M.bases()

        # get the squarefree basis elements
        # we can skip the full groundset (which necessarily has full rank)
        sf_basis = [X for k in range(len(E)) for X in combinations(E, k)
                    if not any(B.issubset(frozenset(X)) for B in bases)]
        lin_basis = []
        comp_max = []
        for i, X in enumerate(sf_basis):
            # we use the fact that the order of X is the same as in E
            X = frozenset(X)
            comp = [v for v in E if v not in X]
            comp_max.append(comp.pop())
            lin_basis.append(tuple(comp))

        # construct the Kozsul complex
        diffs = {}
        cur_basis = {(X, ()): i for i,X in enumerate(sf_basis)}
        for k in range(1, len(E)):
            print(cur_basis)
            prev_basis = cur_basis
            cur_basis = {}
            rows = []
            for X, L, comp in zip(sf_basis, lin_basis, comp_max):
                for Y in combinations(L, k):
                    cur_basis[X,Y] = len(cur_basis)
                    row = [0] * len(prev_basis)

                    for j in range(k):
                        # compute multipliciation by x_i from (x_i - x_{\ell_+})
                        assert Y[j] not in X
                        Yp = Y[:j] + Y[j+1:]
                        Xp = tuple(sorted(X + (Y[j],), key=ordering.__getitem__))
                        try:
                            row[prev_basis[Xp,Yp]] += (-1) ** j
                        except KeyError:
                            # Xp goes to 0, nothing more to do for this term
                            pass

                    # Compute the multiplication by all x_{\ell_+} from the entire
                    #   wedge product. We need to do this separately since
                    #   we need to combine the factors in different ways.
                    # If new_comp in Y, then this will be (-1)^deg (Y \ new_comp)
                    # Otherwise we use the following fact:
                    # d[(a1 - an) ^ (a2 - an) ^ ... ^ (ak - an)] = d[a1 ^ a2 ^ ... ^ ak]
                    # For example:
                    #   sage: E.<a,b,c,d,e> = ExteriorAlgebra(QQ)
                    #   sage: (b+e)*(c+e) - (a+e)*(c+e) + (a+e)*(b+e)
                    #   a*b - a*c + b*c
                    #   sage: (b+e)*(c+e)*(d+e) - (a+e)*(c+e)*(d+e) + (a+e)*(b+e)*(d+e) - (a+e)*(b+e)*(c+e)
                    #   -a*b*c + a*b*d - a*c*d + b*c*d
                    #   sage: e - (a+e)
                    #   -a
                    #   sage: (b+e)*e - (a+e)*e + (a+e)*(b+e)
                    #   a*b
                    #   sage: (b+e)*(c+e)*e - (a+e)*(c+e)*e + (a+e)*(b+e)*e - (a+e)*(b+e)*(c+e)
                    #   -a*b*c

                    Xp = tuple(sorted(X + (comp,), key=ordering.__getitem__))
                    try:
                        ind = sf_basis.index(Xp)
                    except ValueError:
                        # Xp goes to 0, nothing more to do
                        rows.append(row)
                        continue
                    new_comp = comp_max[ind]
                    assert ordering[new_comp] < ordering[comp], (new_comp, comp, X, Y)
                    if new_comp in Y:
                        assert Y[-1] == new_comp
                        Yp = Y[:-1]
                        row[prev_basis[Xp,Yp]] -= (-1) ** len(Yp)
                    else:
                        for j in range(k):
                            Yp = Y[:j] + Y[j+1:]
                            row[prev_basis[Xp,Yp]] -= (-1) ** j

                    rows.append(row)
            diffs[k] = matrix(R, rows).transpose()
        self._kozsul_complex = ChainComplex(diffs, degree_of_differential=-1)
        self._cohomology = self._kozsul_complex.homology()

        # finish the initialization
        C = GradedAlgebrasWithBasis(R)
        indices = [(deg, i) for deg, V in self._cohomology.items() for i in range(V.dimension())]

        CombinatorialFreeModule.__init__(self, R, indices, category=C)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M1 = matroids.catalog.Fano()
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch
            Chow ring of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
            over Rational Field
        """
        return "Tor algebra of {} over {}".format(self._matroid, self.base_ring())

    def _latex_(self):
        r"""
        Return the LaTeX output of the polynomial ring and Chow ring ideal.

        EXAMPLES::

            sage: M1 = matroids.Uniform(2,5)
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch._latex_()
            'A(\\begin{array}{l}\n\\text{\\texttt{U(2,{ }5):{ }Matroid{ }of{ }rank{ }2{ }on{ }5{ }elements{ }with{ }circuit{-}closures}}\\\\\n\\text{\\texttt{{\\char`\\{}2:{ }{\\char`\\{}{\\char`\\{}0,{ }1,{ }2,{ }3,{ }4{\\char`\\}}{\\char`\\}}{\\char`\\}}}}\n\\end{array})_{\\Bold{Q}}'
        """
        from sage.misc.latex import latex
        return r"\mathrm{{Tor}}({})_{{{}}}".format(latex(self._matroid), latex(self.base_ring()))

    def kozsul_complex(self):
        """
        Return the Kozsul complex of ``self``.
        """
        return self._kozsul_complex

    def matroid(self):
        r"""
        Return the matroid of ``self``.

        EXAMPLES::

            sage: ch = matroids.Uniform(3,6).chow_ring(QQ, True, 'fy')
            sage: ch.matroid()
            U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
            {3: {{0, 1, 2, 3, 4, 5}}}
        """
        return self._matroid

    def degree_on_basis(self, m):
        return m[0]

