# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
"""
Orthogonal Symmetric Functions

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


class SymmetricFunctionAlgebra_orthogonal(sfa.SymmetricFunctionAlgebra_generic):
    r"""
    The orthogonal symmetric function basis (or orthogonal basis, to be short).

    The orthogonal basis `\{ o_{\lambda} \}` where `\lambda` is taken over
    all partitions is defined by the following change of basis with the
    Schur functions:

    .. MATH::

        s_{\lambda} = \sum_{\mu} \left( \sum_{\nu \in H} c^{\lambda}_{\mu\nu}
        \right) o_{\mu}

    where `H` is the set of all partitions with even-width rows and
    `c^{\lambda}_{\mu\nu}` is the usual Littlewood-Richardson (LR)
    coefficients. By the properties of LR coefficients, this can be shown to
    be a upper unitriangular change of basis.

    .. NOTE::

        This is only a filtered basis, not a `\ZZ`-graded basis. However this
        does respect the induced `(\ZZ/2\ZZ)`-grading.

    INPUT:

    - ``Sym`` -- an instance of the ring of the symmetric functions

    REFERENCES:

    - [ChariKleber2000]_
    - [KoikeTerada1987]_
    - [ShimozonoZabrocki2006]_

    EXAMPLES:

    Here are the first few orthogonal symmetric functions, in various bases::

        sage: Sym = SymmetricFunctions(QQ)
        sage: o = Sym.o()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: p = Sym.p()
        sage: s = Sym.s()
        sage: m = Sym.m()

        sage: p(o([1]))
        p[1]
        sage: m(o([1]))
        m[1]
        sage: e(o([1]))
        e[1]
        sage: h(o([1]))
        h[1]
        sage: s(o([1]))
        s[1]

        sage: p(o([2]))
        -p[] + 1/2*p[1, 1] + 1/2*p[2]
        sage: m(o([2]))
        -m[] + m[1, 1] + m[2]
        sage: e(o([2]))
        -e[] + e[1, 1] - e[2]
        sage: h(o([2]))
        -h[] + h[2]
        sage: s(o([2]))
        -s[] + s[2]

        sage: p(o([3]))
        -p[1] + 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
        sage: m(o([3]))
        -m[1] + m[1, 1, 1] + m[2, 1] + m[3]
        sage: e(o([3]))
        -e[1] + e[1, 1, 1] - 2*e[2, 1] + e[3]
        sage: h(o([3]))
        -h[1] + h[3]
        sage: s(o([3]))
        -s[1] + s[3]

        sage: Sym = SymmetricFunctions(ZZ)
        sage: o = Sym.o()
        sage: e = Sym.e()
        sage: h = Sym.h()
        sage: s = Sym.s()
        sage: m = Sym.m()
        sage: p = Sym.p()
        sage: m(o([4]))
        -m[1, 1] + m[1, 1, 1, 1] - m[2] + m[2, 1, 1] + m[2, 2] + m[3, 1] + m[4]
        sage: e(o([4]))
        -e[1, 1] + e[1, 1, 1, 1] + e[2] - 3*e[2, 1, 1] + e[2, 2] + 2*e[3, 1] - e[4]
        sage: h(o([4]))
        -h[2] + h[4]
        sage: s(o([4]))
        -s[2] + s[4]

    Some examples of conversions the other way::

        sage: o(h[3])
        o[1] + o[3]
        sage: o(e[3])
        o[1, 1, 1]
        sage: o(m[2,1])
        o[1] - 2*o[1, 1, 1] + o[2, 1]
        sage: o(p[3])
        o[1, 1, 1] - o[2, 1] + o[3]

    Some multiplication::

        sage: o([2]) * o([1,1])
        o[1, 1] + o[2] + o[2, 1, 1] + o[3, 1]
        sage: o([2,1,1]) * o([2])
        o[1, 1] + o[1, 1, 1, 1] + 2*o[2, 1, 1] + o[2, 2] + o[2, 2, 1, 1]
         + o[3, 1] + o[3, 1, 1, 1] + o[3, 2, 1] + o[4, 1, 1]
        sage: o([1,1]) * o([2,1])
        o[1] + o[1, 1, 1] + 2*o[2, 1] + o[2, 1, 1, 1] + o[2, 2, 1]
         + o[3] + o[3, 1, 1] + o[3, 2]

    Examples of the Hopf algebra structure::

        sage: o([1]).antipode()
        -o[1]
        sage: o([2]).antipode()
        -o[] + o[1, 1]
        sage: o([1]).coproduct()
        o[] # o[1] + o[1] # o[]
        sage: o([2]).coproduct()
        o[] # o[] + o[] # o[2] + o[1] # o[1] + o[2] # o[]
        sage: o([1]).counit()
        0
        sage: o.one().counit()
        1
    """

    def __init__(self, Sym):
        """
        Initialize ``self``.

        TESTS::

            sage: o = SymmetricFunctions(QQ).o()
            sage: TestSuite(o).run()
        """
        sfa.SymmetricFunctionAlgebra_generic.__init__(self, Sym, "orthogonal",
                                                      'o', graded=False)

        # We make a strong reference since we use it for our computations
        #   and so we can define the coercion below (only codomains have
        #   strong references)
        self._s = Sym.schur()

        # Setup the coercions
        M = self._s.module_morphism(self._s_to_o_on_basis, codomain=self,
                                    triangular='upper', unitriangular=True)
        M.register_as_coercion()
        Mi = self.module_morphism(self._o_to_s_on_basis, codomain=self._s,
                                  triangular='upper', unitriangular=True)
        Mi.register_as_coercion()

    @cached_method
    def _o_to_s_on_basis(self, lam):
        r"""
        Return the orthogonal symmetric function ``o[lam]`` expanded in
        the Schur basis, where ``lam`` is a partition.

        TESTS:

        Check that this is the inverse::

            sage: o = SymmetricFunctions(QQ).o()
            sage: s = SymmetricFunctions(QQ).s()
            sage: all(o(s(o[la])) == o[la] for i in range(5) for la in Partitions(i))
            True
            sage: all(s(o(s[la])) == s[la] for i in range(5) for la in Partitions(i))
            True
        """
        R = self.base_ring()
        n = sum(lam)
        return self._s._from_dict({ mu: R.sum( (-1)**j * lrcalc.lrcoef_unsafe(lam, mu, nu)
                                               for nu in Partitions(2*j)
                                                   if all(nu.arm_length(i,i) == nu.leg_length(i,i)+1
                                                          for i in range(nu.frobenius_rank()))
                                             )
                                    for j in range(n//2+1) # // 2 for horizontal dominoes
                                    for mu in Partitions(n-2*j) })

    @cached_method
    def _s_to_o_on_basis(self, lam):
        r"""
        Return the Schur symmetric function ``s[lam]`` expanded in
        the orthogonal basis, where ``lam`` is a partition.

        INPUT:

        - ``lam`` -- a partition

        OUTPUT: the expansion of ``s[lam]`` in the orthogonal basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: o = Sym.orthogonal()
            sage: o._s_to_o_on_basis(Partition([]))
            o[]
            sage: o._s_to_o_on_basis(Partition([4,2,1]))
            o[1] + 2*o[2, 1] + o[2, 2, 1] + o[3]
             + o[3, 1, 1] + o[3, 2] + o[4, 1] + o[4, 2, 1]
            sage: s(o._s_to_o_on_basis(Partition([3,1]))) == s[3,1]
            True
        """
        R = self.base_ring()
        n = sum(lam)
        return self._from_dict({ mu: R.sum( lrcalc.lrcoef_unsafe(lam, mu, [2*x for x in nu])
                                            for nu in Partitions(j) )
                                 for j in range(n//2+1) # // 2 for horizontal dominoes
                                 for mu in Partitions(n-2*j) })
