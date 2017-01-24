"""
Hochschild Complexes
"""

########################################################################
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import CategoryObject
from sage.categories.category_types import ChainComplexes
from sage.categories.tensor import tensor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.homology.chain_complex import ChainComplex

class HochschildComplex(UniqueRepresentation, CategoryObject):
    r"""
    The Hochschild complex.

    Let `A` be an algebra over a commutative ring `R` such
    that `A` a projective `R`-module, and `M` an `A`-bimodule.
    The *Hochschild complex* is the chain complex given by

    .. MATH::

        C_n(A, M) := M \otimes A^{\otimes n}

    with the boundary operators given as follows. For fixed `n`, define
    the face maps

    .. MATH::

        f_{n,i}(m \otimes a_1 \otimes \cdots \otimes a_n) = \begin{cases}
        m a_1 \otimes \cdots \otimes a_n & \text{if } i = 0, \\
        a_n m \otimes a_1 \otimes \cdots \otimes a_{n-1}
        & \text{if } i = n, \\
        m \otimes a_1 \otimes \cdots \otimes a_i a_{i+1} \otimes
        \cdots \otimes a_n & \text{otherwise.}
        \end{cases}

    We define the boundary operators as

    .. MATH::

        d_n = \sum_{i=0}^n (-1)^i f_{n,i}.

    The *Hochschild homology* of `A` is the homology of this complex.
    Alternatively, the Hochschild homology can be described by
    `HH_n(A, M) = \operatorname{Tor}_n^{A^e}(A, M)`, where
    `A^e = A \otimes A^o` (`A^o` is the opposite algebra of `A`)
    is the enveloping algebra of `A`.

    *Hochschild cohomology* is the homology of the dual complex and
    can be described by `HH^n(A, M) = \operatorname{Ext}^n_{A^e}(A, M)`.

    Another perspective on Hochschild homology is that `f_{n,i}`
    make the family `C_n(A, M)` a simplicial object in the
    category of `R`-modules, and the degeneracy maps are

    .. MATH::

        s_i(a_0 \otimes \cdots \otimes a_n) =
        a_0 \otimes \cdots \otimes a_i \otimes 1 \otimes a_{i+1}
        \otimes \cdots \otimes a_n

    The Hochschild homology can also be constructed as the homology
    of this simplicial module.

    REFERENCES:

    - :wikipedia:`Hochschild_homology`
    - https://ncatlab.org/nlab/show/Hochschild+cohomology
    - [Red2001]_
    """
    def __init__(self, A, M):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)

            sage: H.category()
            Category of chain complexes over Rational Field
            sage: H in ChainComplexes(QQ)
            True

        Some methods required by the category are not implemented::

            sage: TestSuite(H).run()  # known bug (#21386)
        """
        self._A = A
        self._M = M
        CategoryObject.__init__(self, base=A.base_ring(),
                                category=ChainComplexes(A.base_ring()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: T.rename("Trivial representation of SGA")
            sage: SGA.hochschild_complex(T)
            Hochschild complex of Symmetric group algebra of order 3 over Rational Field
             with coefficients in Trivial representation of SGA
            sage: T.rename()  # reset the name
        """
        return "Hochschild complex of {} with coefficients in {}".format(self._A, self._M)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: latex(H)
            C_{\bullet}\left(..., ...\right)
        """
        from sage.misc.latex import latex
        return "C_{{\\bullet}}\\left({}, {}\\right)".format(latex(self._A), latex(self._M))

    def algebra(self):
        """
        Return the defining algebra of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.algebra()
            Symmetric group algebra of order 3 over Rational Field
        """
        return self._A

    def coefficients(self):
        """
        Return the coefficients of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.coefficients()
            Trivial representation of Standard permutations of 3 over Rational Field
        """
        return self._M

    def free_module(self, d):
        """
        Return the free module in degree ``d``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.free_module(0)
            Trivial representation of Standard permutations of 3 over Rational Field
            sage: H.free_module(1)
            Trivial representation of Standard permutations of 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
            sage: H.free_module(2)
            Trivial representation of Standard permutations of 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
        """
        if d < 0:
            raise ValueError("only defined for non-negative degree")
        return tensor([self._M] + [self._A]*d)

    @cached_method
    def trivial_module(self):
        """
        Return the trivial module of ``self``.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: H.trivial_module()
            Free module generated by {} over Rational Field
        """
        return CombinatorialFreeModule(self._A.base_ring(), [])

    def boundary(self, d):
        """
        Return the boundary operator in degree ``d``.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: d1 = H.boundary(1)
            sage: z = d1.domain().an_element(); z
            2*1 # 1 + 2*1 # x + 3*1 # y
            sage: d1(z)
            0
            sage: d1.matrix()
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  2  0  0 -2  0  0  0  0  0  0]

            sage: s = SymmetricFunctions(QQ).s()
            sage: H = s.hochschild_complex(s)
            sage: d1 = H.boundary(1)
            sage: x = d1.domain().an_element(); x
            2*s[] # s[] + 2*s[] # s[1] + 3*s[] # s[2]
            sage: d1(x)
            0
            sage: y = tensor([s.an_element(), s.an_element()])
            sage: d1(y)
            0
            sage: z = tensor([s[2,1] + s[3], s.an_element()])
            sage: d1(z)
            0

        TESTS::

            sage: def test_complex(H, n):
            ....:     phi = H.boundary(n)
            ....:     psi = H.boundary(n+1)
            ....:     comp = phi * psi
            ....:     zero = H.free_module(n-1).zero()
            ....:     return all(comp(b) == zero for b in H.free_module(n+1).basis())

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: H = SGA.hochschild_complex(SGA)
            sage: test_complex(H, 1)
            True
            sage: test_complex(H, 2)
            True
            sage: test_complex(H, 3) # long time
            True

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: test_complex(H, 1)
            True
            sage: test_complex(H, 2)
            True
            sage: test_complex(H, 3)
            True
        """
        R = self._A.base_ring()
        one = R.one()
        if d == 0:
            t = self.trivial_module()
            zero = t.zero()
            return self.free_module(0).module_morphism(lambda x: zero, codomain=t)
        Fd = self.free_module(d-1)
        Fd1 = self.free_module(d)
        mone = -one
        def on_basis(k):
            p = self._M.monomial(k[0]) * self._A.monomial(k[1])
            ret = Fd._from_dict({(m,) + k[2:]: c for m,c in p}, remove_zeros=False)
            for i in range(1, d):
                p = self._A.monomial(k[i]) * self._A.monomial(k[i+1])
                ret += mone**i * Fd._from_dict({k[:i] + (m,) + k[i+2:]: c
                                                   for m,c in p}, remove_zeros=False)
            p = self._A.monomial(k[-1]) * self._M.monomial(k[0])
            ret += mone**d * Fd._from_dict({(m,) + k[1:-1]: c for m,c in p},
                                           remove_zeros=False)
            return ret
        return Fd1.module_morphism(on_basis, codomain=Fd)

    def coboundary(self, d):
        """
        Return the coboundary morphism of degree ``d``.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: del1 = H.coboundary(1)
            sage: z = del1.domain().an_element(); z
            2 + 2*x + 3*y
            sage: del1(z)
            0
            sage: del1.matrix()
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  2]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0 -2]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]

        TESTS::

            sage: def test_complex(H, n):
            ....:     phi = H.coboundary(n)
            ....:     psi = H.coboundary(n+1)
            ....:     comp = psi * phi
            ....:     zero = H.free_module(n+1).zero()
            ....:     return all(comp(b) == zero for b in H.free_module(n-1).basis())

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: H = SGA.hochschild_complex(SGA)
            sage: test_complex(H, 1)
            True
            sage: test_complex(H, 2)
            True

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: test_complex(H, 1)
            True
            sage: test_complex(H, 2)
            True
            sage: test_complex(H, 3)
            True
        """
        if self._A.category() is not self._A.category().FiniteDimensional():
            raise NotImplementedError("the algebra must be finite dimensional")
        bdry = self.boundary(d)
        dom = bdry.domain()
        cod = bdry.codomain()
        return cod.module_morphism(matrix=bdry.matrix().transpose(), codomain=dom)

    def homology(self, d):
        """
        Return the ``d``-th homology group.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: H.homology(0)
            Vector space of dimension 3 over Rational Field
            sage: H.homology(1)
            Vector space of dimension 4 over Rational Field
            sage: H.homology(2)
            Vector space of dimension 6 over Rational Field

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.homology(0)
            Vector space of dimension 1 over Rational Field
            sage: H.homology(1)
            Vector space of dimension 0 over Rational Field
            sage: H.homology(2)
            Vector space of dimension 0 over Rational Field

        When working over general rings (except `\ZZ`) and we can
        construct a unitriangular basis for the image quotient,
        we fallback to a slower implementation using (combinatorial)
        free modules::

            sage: R.<x,y> = QQ[]
            sage: SGA = SymmetricGroupAlgebra(R, 2)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.homology(1)
            Free module generated by {} over Multivariate Polynomial Ring in x, y over Rational Field
        """
        if self._A.category() is not self._A.category().FiniteDimensional():
            raise NotImplementedError("the algebra must be finite dimensional")

        maps = {d: self.boundary(d).matrix(), d+1: self.boundary(d+1).matrix()}
        C = ChainComplex(maps, degree_of_differential=-1)
        try:
            return C.homology(d)
        except NotImplementedError:
            pass
        # Fallback if we are not working over a field or \ZZ
        bdry = self.boundary(d)
        bdry1 = self.boundary(d+1)
        ker = bdry.kernel()
        im_retract = ker.submodule([ker.retract(b) for b in bdry1.image_basis()],
                                   unitriangular=True)
        return ker.quotient_module(im_retract)

    def cohomology(self, d):
        """
        Return the ``d``-th cohomology group.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: H.cohomology(0)
            Vector space of dimension 3 over Rational Field
            sage: H.cohomology(1)
            Vector space of dimension 4 over Rational Field
            sage: H.cohomology(2)
            Vector space of dimension 6 over Rational Field

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.cohomology(0)
            Vector space of dimension 1 over Rational Field
            sage: H.cohomology(1)
            Vector space of dimension 0 over Rational Field
            sage: H.cohomology(2)
            Vector space of dimension 0 over Rational Field

        When working over general rings (except `\ZZ`) and we can
        construct a unitriangular basis for the image quotient,
        we fallback to a slower implementation using (combinatorial)
        free modules::

            sage: R.<x,y> = QQ[]
            sage: SGA = SymmetricGroupAlgebra(R, 2)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.cohomology(1)
            Free module generated by {} over Multivariate Polynomial Ring in x, y over Rational Field
        """
        if self._A.category() is not self._A.category().FiniteDimensional():
            raise NotImplementedError("the algebra must be finite dimensional")

        maps = {d+1: self.coboundary(d+1).matrix(), d: self.coboundary(d).matrix()}
        C = ChainComplex(maps, degree_of_differential=1)
        try:
            return C.homology(d+1)
        except NotImplementedError:
            pass
        # Fallback if we are not working over a field or \ZZ
        cb = self.coboundary(d)
        cb1 = self.coboundary(d+1)
        ker = cb1.kernel()
        im_retract = ker.submodule([ker.retract(b) for b in cb.image_basis()],
                                   unitriangular=True)
        return ker.quotient_module(im_retract)

