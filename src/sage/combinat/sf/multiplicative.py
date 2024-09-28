# sage.doctest: needs sage.combinat sage.modules
r"""
Multiplicative symmetric functions

A realization `h` of the ring of symmetric functions is multiplicative if for
a partition `\lambda = (\lambda_1,\lambda_2,\ldots)` we have
`h_\lambda = h_{\lambda_1} h_{\lambda_2} \cdots`.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from . import classical
import sage.combinat.partition


class SymmetricFunctionAlgebra_multiplicative(classical.SymmetricFunctionAlgebra_classical):
    r"""
    The class of multiplicative bases of the ring of symmetric functions.

    A realization `q` of the ring of symmetric functions is multiplicative if
    for a partition `\lambda = (\lambda_1,\lambda_2,\ldots)` we have
    `q_\lambda = q_{\lambda_1} q_{\lambda_2} \cdots` (with `q_0` meaning `1`).

    Examples of multiplicative realizations are the elementary symmetric basis,
    the complete homogeneous basis, the powersum basis (if the base ring is a
    `\QQ`-algebra), and the Witt basis (but not the Schur basis or the
    monomial basis).
    """

    def product_on_basis(self, left, right):
        """
        Return the product of ``left`` and ``right``.

        INPUT:

        - ``left``, ``right`` -- partitions

        OUTPUT: an element of ``self``

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e([2,1])^2  # indirect doctest
            e[2, 2, 1, 1]

        ::

            sage: h = SymmetricFunctions(QQ).h()
            sage: h([2,1])^2
            h[2, 2, 1, 1]

        ::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p([2,1])^2
            p[2, 2, 1, 1]

        ::

            sage: QQx.<x> = QQ[]
            sage: p = SymmetricFunctions(QQx).p()
            sage: (x*p([2]))^2
            x^2*p[2, 2]

            sage: TestSuite(p).run() # to silence sage -coverage
        """
        m = list(left) + list(right)
        m.sort(reverse=True)
        return self.monomial(sage.combinat.partition.Partition(m))

    def coproduct_on_basis(self, mu):
        r"""
        Return the coproduct on a basis element for multiplicative bases.

        INPUT:

        - ``mu`` -- a partition

        OUTPUT:

        - the image of ``self[mu]`` under comultiplication; this is an
          element of the tensor square of ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.powersum()
            sage: p.coproduct_on_basis([2,1])
            p[] # p[2, 1] + p[1] # p[2] + p[2] # p[1] + p[2, 1] # p[]

            sage: e = Sym.elementary()
            sage: e.coproduct_on_basis([3,1])
            e[] # e[3, 1] + e[1] # e[2, 1] + e[1] # e[3] + e[1, 1] # e[2] + e[2] # e[1, 1] + e[2, 1] # e[1] + e[3] # e[1] + e[3, 1] # e[]

            sage: h = Sym.homogeneous()
            sage: h.coproduct_on_basis([3,1])
            h[] # h[3, 1] + h[1] # h[2, 1] + h[1] # h[3] + h[1, 1] # h[2] + h[2] # h[1, 1] + h[2, 1] # h[1] + h[3] # h[1] + h[3, 1] # h[]
        """
        T = self.tensor_square()
        return T.prod(self.coproduct_on_generators(p) for p in mu)
