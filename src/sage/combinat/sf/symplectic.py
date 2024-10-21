# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
"""
Symplectic Symmetric Functions

AUTHORS:

- Travis Scrimshaw (2013-11-10): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import sfa
import sage.libs.lrcalc.lrcalc as lrcalc
from sage.combinat.partition import Partitions
from sage.misc.cachefunc import cached_method


class SymmetricFunctionAlgebra_symplectic(sfa.SymmetricFunctionAlgebra_generic):
    r"""
    The symplectic symmetric function basis (or symplectic basis, to be short).

    The symplectic basis `\{ sp_{\lambda} \}` where `\lambda` is taken over
    all partitions is defined by the following change of basis with the
    Schur functions:

    .. MATH::

        s_{\lambda} = \sum_{\mu} \left( \sum_{\nu \in V} c^{\lambda}_{\mu\nu}
        \right) sp_{\mu}

    where `V` is the set of all partitions with even-height columns and
    `c^{\lambda}_{\mu\nu}` is the usual Littlewood-Richardson (LR)
    coefficients. By the properties of LR coefficients, this can be shown to
    be a upper unitriangular change of basis.

    .. NOTE::

        This is only a filtered basis, not a `\ZZ`-graded basis. However this
        does respect the induced `(\ZZ/2\ZZ)`-grading.

    INPUT:

    - ``Sym`` -- an instance of the ring of the symmetric functions

    REFERENCES:

    .. [ChariKleber2000] Vyjayanthi Chari and Michael Kleber.
       *Symmetric functions and representations of quantum affine algebras*.
       :arxiv:`math/0011161v1`

    .. [KoikeTerada1987] \K. Koike, I. Terada, *Young-diagrammatic methods for
       the representation theory of the classical groups of type Bn, Cn, Dn*.
       J. Algebra 107 (1987), no. 2, 466-511.

    .. [ShimozonoZabrocki2006] Mark Shimozono and Mike Zabrocki.
       *Deformed universal characters for classical and affine algebras*.
       Journal of Algebra, **299** (2006). :arxiv:`math/0404288`.

    EXAMPLES:

    Here are the first few symplectic symmetric functions, in various bases::

        sage: Sym = SymmetricFunctions(QQ)
        sage: sp = Sym.sp()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(sp([1]))
        p[1]
        sage: m(sp([1]))
        m[1]
        sage: e(sp([1]))
        e[1]
        sage: h(sp([1]))
        h[1]
        sage: s(sp([1]))
        s[1]

        sage: p(sp([2]))
        1/2*p[1, 1] + 1/2*p[2]
        sage: m(sp([2]))
        m[1, 1] + m[2]
        sage: e(sp([2]))
        e[1, 1] - e[2]
        sage: h(sp([2]))
        h[2]
        sage: s(sp([2]))
        s[2]

        sage: p(sp([3]))
        1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
        sage: m(sp([3]))
        m[1, 1, 1] + m[2, 1] + m[3]
        sage: e(sp([3]))
        e[1, 1, 1] - 2*e[2, 1] + e[3]
        sage: h(sp([3]))
        h[3]
        sage: s(sp([3]))
        s[3]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: sp = Sym.sp()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(sp([4]))
        m[1, 1, 1, 1] + m[2, 1, 1] + m[2, 2] + m[3, 1] + m[4]
        sage: e(sp([4]))
        e[1, 1, 1, 1] - 3*e[2, 1, 1] + e[2, 2] + 2*e[3, 1] - e[4]
        sage: h(sp([4]))
        h[4]
        sage: s(sp([4]))
        s[4]

    Some examples of conversions the other way::

        sage: sp(h[3])
        sp[3]
        sage: sp(e[3])
        sp[1] + sp[1, 1, 1]
        sage: sp(m[2,1])
        -sp[1] - 2*sp[1, 1, 1] + sp[2, 1]
        sage: sp(p[3])
        sp[1, 1, 1] - sp[2, 1] + sp[3]

    Some multiplication::

        sage: sp([2]) * sp([1,1])
        sp[1, 1] + sp[2] + sp[2, 1, 1] + sp[3, 1]
        sage: sp([2,1,1]) * sp([2])
        sp[1, 1] + sp[1, 1, 1, 1] + 2*sp[2, 1, 1] + sp[2, 2] + sp[2, 2, 1, 1]
         + sp[3, 1] + sp[3, 1, 1, 1] + sp[3, 2, 1] + sp[4, 1, 1]
        sage: sp([1,1]) * sp([2,1])
        sp[1] + sp[1, 1, 1] + 2*sp[2, 1] + sp[2, 1, 1, 1] + sp[2, 2, 1]
         + sp[3] + sp[3, 1, 1] + sp[3, 2]

    Examples of the Hopf algebra structure::

        sage: sp([1]).antipode()
        -sp[1]
        sage: sp([2]).antipode()
        sp[] + sp[1, 1]
        sage: sp([1]).coproduct()
        sp[] # sp[1] + sp[1] # sp[]
        sage: sp([2]).coproduct()
        sp[] # sp[2] + sp[1] # sp[1] + sp[2] # sp[]
        sage: sp([1]).counit()
        0
        sage: sp.one().counit()
        1
    """

    def __init__(self, Sym):
        """
        Initialize ``self``.

        TESTS::

            sage: sp = SymmetricFunctions(QQ).sp()
            sage: TestSuite(sp).run()
        """
        sfa.SymmetricFunctionAlgebra_generic.__init__(self, Sym, "symplectic",
                                                      'sp', graded=False)

        # We make a strong reference since we use it for our computations
        #   and so we can define the coercion below (only codomains have
        #   strong references)
        self._s = Sym.schur()

        # Setup the coercions
        M = self._s.module_morphism(self._s_to_sp_on_basis, codomain=self,
                                    triangular='upper', unitriangular=True)
        M.register_as_coercion()
        Mi = self.module_morphism(self._sp_to_s_on_basis, codomain=self._s,
                                  triangular='upper', unitriangular=True)
        Mi.register_as_coercion()

    @cached_method
    def _sp_to_s_on_basis(self, lam):
        r"""
        Return the symplectic symmetric function ``sp[lam]`` expanded in
        the Schur basis, where ``lam`` is a partition.

        TESTS:

        Check that this is the inverse::

            sage: sp = SymmetricFunctions(QQ).sp()
            sage: s = SymmetricFunctions(QQ).s()
            sage: all(sp(s(sp[la])) == sp[la] for i in range(5) for la in Partitions(i))
            True
            sage: all(s(sp(s[la])) == s[la] for i in range(5) for la in Partitions(i))
            True
        """
        R = self.base_ring()
        n = sum(lam)
        return self._s._from_dict({ mu: R.sum( (-1)**j * lrcalc.lrcoef_unsafe(lam, mu, nu)
                                               for nu in Partitions(2*j)
                                                   if all(nu.leg_length(i,i) == nu.arm_length(i,i)+1
                                                          for i in range(nu.frobenius_rank()))
                                             )
                                    for j in range(n//2+1) # // 2 for horizontal dominoes
                                    for mu in Partitions(n-2*j) })

    @cached_method
    def _s_to_sp_on_basis(self, lam):
        r"""
        Return the Schur symmetric function ``s[lam]`` expanded in
        the symplectic basis, where ``lam`` is a partition.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT: the expansion of ``s[lam]`` in the symplectic basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: sp = Sym.symplectic()
            sage: sp._s_to_sp_on_basis(Partition([]))
            sp[]
            sage: sp._s_to_sp_on_basis(Partition([4,2,1]))
            sp[2, 1] + sp[3] + sp[3, 1, 1] + sp[3, 2] + sp[4, 1] + sp[4, 2, 1]
            sage: s(sp._s_to_sp_on_basis(Partition([3,1]))) == s[3,1]
            True
        """
        R = self.base_ring()
        n = sum(lam)
        return self._from_dict({ mu: R.sum( lrcalc.lrcoef_unsafe(lam, mu, sum([[x,x] for x in nu], []))
                                            for nu in Partitions(j) )
                                 for j in range(n//2+1) # // 2 for vertical dominoes
                                 for mu in Partitions(n-2*j) })
