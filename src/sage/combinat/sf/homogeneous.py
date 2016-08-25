"""
Homogeneous symmetric functions

By this we mean the basis formed of the complete homogeneous
symmetric functions `h_\lambda`, not an arbitrary graded basis.
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

####################################
#                                  #
# Homogeneous Symmetric Functions  #
#                                  #
####################################
from . import multiplicative, classical
from sage.combinat.partition import Partition

class SymmetricFunctionAlgebra_homogeneous(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class of methods specific to the homogeneous basis of
        symmetric functions.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: h = SymmetricFunctions(QQ).e()
            sage: h == loads(dumps(h))
            True
            sage: TestSuite(h).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(h).run(elements = [h[1,1]+h[2], h[1]+2*h[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "homogeneous", 'h')

    def _dual_basis_default(self):
        r"""
        Returns the dual basis to ``self``.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``scalar`` -- optional input which specifies a function ``zee``
          on partitions. The function ``zee`` determines the scalar
          product on the power sum basis with normalization
          `\langle p_\mu, p_\mu \rangle = \mathrm{zee}(mu)`.
          (default: uses standard ``zee`` function)
        - ``scalar_name`` -- specifies the name of the scalar function
          (optional)
        - ``prefix`` -- optional input, specifies the prefix to be
          used to display the basis.

        OUTPUT:

        The dual basis of the homogeneous basis with respect to the
        standard scalar product (the monomial basis).  If a function
        ``zee`` is specified, the dual basis is with respect to the
        modified scalar product.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: h = SymmetricFunctions(QQ).h()
            sage: h.dual_basis() == m
            True

            sage: zee = lambda x : 2
            sage: hh = h.dual_basis(zee); hh
            Dual basis to Symmetric Functions over Rational Field in the homogeneous basis
            sage: hh[2,1].scalar(h[2,1])
            1
            sage: hh[2,2].scalar(h[2,2])
            4

        TESTS::

            sage: h._dual_basis_default() is h.dual_basis()
            True
        """
        return self.realization_of().m()

    def coproduct_on_generators(self, i):
        r"""
        Returns the coproduct on `h_i`.

        INPUT:

        - ``self`` -- a homogeneous basis of symmetric functions
        - ``i`` -- a nonnegative integer

        OUTPUT:

        - the sum `\sum_{r=0}^i h_r \otimes h_{i-r}`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.homogeneous()
            sage: h.coproduct_on_generators(2)
            h[] # h[2] + h[1] # h[1] + h[2] # h[]
            sage: h.coproduct_on_generators(0)
            h[] # h[]
        """
        def P(i): return Partition([i]) if i else Partition([])
        T = self.tensor_square()
        return T.sum_of_monomials( (P(j), P(i-j)) for j in range(i+1) )


    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        def omega(self):
            r"""
            Return the image of ``self`` under the omega automorphism.

            The *omega automorphism* is defined to be the unique algebra
            endomorphism `\omega` of the ring of symmetric functions that
            satisfies `\omega(e_k) = h_k` for all positive integers `k`
            (where `e_k` stands for the `k`-th elementary symmetric
            function, and `h_k` stands for the `k`-th complete homogeneous
            symmetric function). It furthermore is a Hopf algebra
            endomorphism and an involution, and it is also known as the
            *omega involution*. It sends the power-sum symmetric function
            `p_k` to `(-1)^{k-1} p_k` for every positive integer `k`.

            The images of some bases under the omega automorphism are given by

            .. MATH::

                \omega(e_{\lambda}) = h_{\lambda}, \qquad
                \omega(h_{\lambda}) = e_{\lambda}, \qquad
                \omega(p_{\lambda}) = (-1)^{|\lambda| - \ell(\lambda)}
                p_{\lambda}, \qquad
                \omega(s_{\lambda}) = s_{\lambda^{\prime}},

            where `\lambda` is any partition, where `\ell(\lambda)` denotes
            the length (:meth:`~sage.combinat.partition.Partition.length`)
            of the partition `\lambda`, where `\lambda^{\prime}` denotes the
            conjugate partition
            (:meth:`~sage.combinat.partition.Partition.conjugate`) of
            `\lambda`, and where the usual notations for bases are used
            (`e` = elementary, `h` = complete homogeneous, `p` = powersum,
            `s` = Schur).

            :meth:`omega_involution()` is a synonym for the :meth:`omega()`
            method.

            OUTPUT:

            - the image of ``self`` under the omega automorphism

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: a = h([2,1]); a
                h[2, 1]
                sage: a.omega()
                h[1, 1, 1] - h[2, 1]
                sage: e = SymmetricFunctions(QQ).e()
                sage: e(h([2,1]).omega())
                e[2, 1]
            """
            e = self.parent().realization_of().e()
            return self.parent()(e._from_element(self))

        omega_involution = omega

        def expand(self, n, alphabet='x'):
            """
            Expand the symmetric function ``self`` as a symmetric polynomial
            in ``n`` variables.

            INPUT:

            - ``n`` -- a nonnegative integer

            - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

            OUTPUT:

            A monomial expansion of ``self`` in the `n` variables
            labelled by ``alphabet``.

            EXAMPLES::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h([3]).expand(2)
                x0^3 + x0^2*x1 + x0*x1^2 + x1^3
                sage: h([1,1,1]).expand(2)
                x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
                sage: h([2,1]).expand(3)
                x0^3 + 2*x0^2*x1 + 2*x0*x1^2 + x1^3 + 2*x0^2*x2 + 3*x0*x1*x2 + 2*x1^2*x2 + 2*x0*x2^2 + 2*x1*x2^2 + x2^3
                sage: h([3]).expand(2,alphabet='y')
                y0^3 + y0^2*y1 + y0*y1^2 + y1^3
                sage: h([3]).expand(2,alphabet='x,y')
                x^3 + x^2*y + x*y^2 + y^3
                sage: h([3]).expand(3,alphabet='x,y,z')
                x^3 + x^2*y + x*y^2 + y^3 + x^2*z + x*y*z + y^2*z + x*z^2 + y*z^2 + z^3
                sage: (h([]) + 2*h([1])).expand(3)
                2*x0 + 2*x1 + 2*x2 + 1
                sage: h([1]).expand(0)
                0
                sage: (3*h([])).expand(0)
                3
            """
            if n == 0:   # Symmetrica crashes otherwise...
                return self.counit()
            condition = lambda part: False
            return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.homogeneous', 'SymmetricFunctionAlgebraElement_homogeneous',  SymmetricFunctionAlgebra_homogeneous.Element)
