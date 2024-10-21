# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat        (because all doctests use FreeAlgebra, SymmetricGroupAlgebra, etc.)
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
#                  https://www.gnu.org/licenses/
########################################################################

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import ModuleElement, parent
from sage.structure.richcmp import richcmp
from sage.categories.chain_complexes import ChainComplexes
from sage.categories.tensor import tensor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.homology.chain_complex import ChainComplex, Chain_class


class HochschildComplex(UniqueRepresentation, Parent):
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

            sage: TestSuite(H).run()
        """
        self._A = A
        self._M = M
        Parent.__init__(self, base=A.base_ring(),
                        category=ChainComplexes(A.base_ring()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: T.rename('Trivial representation of SGA')
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

    def module(self, d):
        """
        Return the module in degree ``d``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H.module(0)
            Trivial representation of Standard permutations of 3 over Rational Field
            sage: H.module(1)
            Trivial representation of Standard permutations of 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
            sage: H.module(2)
            Trivial representation of Standard permutations of 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
             # Symmetric group algebra of order 3 over Rational Field
        """
        if d < 0:
            raise ValueError("only defined for nonnegative degree")
        return tensor([self._M] + [self._A] * d)

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
            ....:     zero = H.module(n-1).zero()
            ....:     return all(comp(b) == zero for b in H.module(n+1).basis())

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
            return self.module(0).module_morphism(lambda x: zero, codomain=t)
        Fd = self.module(d-1)
        Fd1 = self.module(d)
        mone = -one

        def on_basis(k):
            p = self._M.monomial(k[0]) * self._A.monomial(k[1])
            ret = Fd._from_dict({(m,) + k[2:]: c for m, c in p}, remove_zeros=False)
            for i in range(1, d):
                p = self._A.monomial(k[i]) * self._A.monomial(k[i+1])
                ret += mone**i * Fd._from_dict({k[:i] + (m,) + k[i+2:]: c
                                                for m, c in p}, remove_zeros=False)
            p = self._A.monomial(k[-1]) * self._M.monomial(k[0])
            ret += mone**d * Fd._from_dict({(m,) + k[1:-1]: c for m, c in p},
                                           remove_zeros=False)
            return ret
        return Fd1.module_morphism(on_basis, codomain=Fd)

    differential = boundary

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
            ....:     zero = H.module(n+1).zero()
            ....:     return all(comp(b) == zero for b in H.module(n-1).basis())

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
        r"""
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
        r"""
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

    def _element_constructor_(self, vectors):
        """
        Construct an element of ``self`` from ``vectors``.

        TESTS::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: H = E.hochschild_complex(E)
            sage: H(0)
            Trivial chain
            sage: H(2)
            Chain(0: 2)
            sage: H(x+2*y)
            Chain(0: x + 2*y)
            sage: H({0: H.module(0).an_element()})
            Chain(0: 2 + 2*x + 3*y)
            sage: H({2: H.module(2).an_element()})
            Chain(2: 2*1 # 1 # 1 + 2*1 # 1 # x + 3*1 # 1 # y)
            sage: H({0:x-y, 2: H.module(2).an_element()})
            Chain with 2 nonzero terms over Rational Field
            sage: H([2])
            Traceback (most recent call last):
            ...
            ValueError: cannot construct an element from [2]
        """
        if not vectors:  # special case: the zero chain
            return self.element_class(self, {})
        # special case: an element of the defining module
        if self._M.has_coerce_map_from(parent(vectors)):
            vectors = self._M(vectors)
        if parent(vectors) is self._M:
            mc = vectors.monomial_coefficients(copy=False)
            vec = self.module(0)._from_dict({(k,): mc[k] for k in mc})
            return self.element_class(self, {0: vec})
        if isinstance(vectors, (Chain_class, self.element_class)):
            vectors = vectors._vec
        data = {}
        if not isinstance(vectors, dict):
            raise ValueError("cannot construct an element from {}".format(vectors))
        # Special handling for the 0 free module
        # FIXME: Allow coercions between the 0 free module and the defining module
        if 0 in vectors:
            vec = vectors.pop(0)
            if parent(vec) is self._M:
                mc = vec.monomial_coefficients(copy=False)
                data[0] = self.module(0)._from_dict({(k,): mc[k] for k in mc})
            else:
                data[0] = self.module(0)(vec)
        for degree in vectors:
            vec = self.module(degree)(vectors[degree])
            if not vec:
                continue
            data[degree] = vec
        return self.element_class(self, data)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(ZZ)
            sage: H = F.hochschild_complex(F)
            sage: v = H.an_element()
            sage: [v.vector(i) for i in range(6)]
            [2*F[1] + 2*F[x] + 3*F[y],
             2*F[1] # F[1] + 2*F[1] # F[x] + 3*F[1] # F[y],
             2*F[1] # F[1] # F[1] + 2*F[1] # F[1] # F[x] + 3*F[1] # F[1] # F[y],
             2*F[1] # F[1] # F[1] # F[1] + 2*F[1] # F[1] # F[1] # F[x]
                + 3*F[1] # F[1] # F[1] # F[y],
             0,
             0]
        """
        return self.element_class(self, {d: self.module(d).an_element()
                                         for d in range(4)})

    class Element(ModuleElement):
        """
        A chain of the Hochschild complex.

        INPUT:

        Can be one of the following:

        - A dictionary whose keys are the degree and whose `d`-th
          value is an element in the degree `d` module.
        - An element in the coefficient module `M`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: T = SGA.trivial_representation()
            sage: H = SGA.hochschild_complex(T)
            sage: H(T.an_element())
            Chain(0: 2*B['v'])
            sage: H({0: T.an_element()})
            Chain(0: 2*B['v'])
            sage: H({1: H.module(1).an_element()})
            Chain(1: 2*B['v'] # [1, 2, 3] + 2*B['v'] # [1, 3, 2] + 3*B['v'] # [2, 1, 3])
            sage: H({0: H.module(0).an_element(), 3: H.module(3).an_element()})
            Chain with 2 nonzero terms over Rational Field

            sage: F.<x,y> = FreeAlgebra(ZZ)
            sage: H = F.hochschild_complex(F)
            sage: H(x + 2*y^2)
            Chain(0: F[x] + 2*F[y^2])
            sage: H({0: x*y - x})
            Chain(0: -F[x] + F[x*y])
            sage: H(2)
            Chain(0: 2*F[1])
            sage: H({0: x-y, 2: H.module(2).basis().an_element()})
            Chain with 2 nonzero terms over Integer Ring
        """
        def __init__(self, parent, vectors):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x-y, 2: H.module(2).basis().an_element()})
                sage: TestSuite(a).run()
            """
            self._vec = vectors
            ModuleElement.__init__(self, parent)

        def vector(self, degree):
            """
            Return the free module element in ``degree``.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x-y, 2: H.module(2).basis().an_element()})
                sage: [a.vector(i) for i in range(3)]
                [F[x] - F[y], 0, F[1] # F[1] # F[1]]
            """
            try:
                return self._vec[degree]
            except KeyError:
                return self.parent().module(degree).zero()

        def _repr_(self):
            """
            Print representation.

            EXAMPLES::

                sage: E.<x,y> = ExteriorAlgebra(QQ)
                sage: H = E.hochschild_complex(E)
                sage: H(0)
                Trivial chain
                sage: H(x+2*y)
                Chain(0: x + 2*y)
                sage: H({2: H.module(2).an_element()})
                Chain(2: 2*1 # 1 # 1 + 2*1 # 1 # x + 3*1 # 1 # y)
                sage: H({0:x-y, 2: H.module(2).an_element()})
                Chain with 2 nonzero terms over Rational Field
            """
            n = len(self._vec)
            if n == 0:
                return 'Trivial chain'

            if n == 1:
                (deg, vec), = self._vec.items()
                return f'Chain({deg}: {vec})'

            return f'Chain with {n} nonzero terms over {self.parent().base_ring()}'

        def _ascii_art_(self):
            """
            Return an ascii art representation.

            Note that arrows go to the left so that composition of
            differentials is the usual matrix multiplication.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x - y,
                ....:        1: H.module(1).basis().an_element(),
                ....:        2: H.module(2).basis().an_element()})
                sage: ascii_art(a)
                   d_0           d_1         d_2             d_3
                0 <---- F  - F  <---- 1 # 1 <---- 1 # 1 # 1 <---- 0
                         x    y
            """
            from sage.typeset.ascii_art import AsciiArt, ascii_art

            if not self._vec:   # 0 chain
                return AsciiArt(['0'])

            def arrow_art(d):
                d_str = ['  d_{0}  '.format(d)]
                arrow = ' <' + '-'*(len(d_str[0])-3) + ' '
                d_str.append(arrow)
                return AsciiArt(d_str, baseline=0)

            result = AsciiArt(['0'])
            max_deg = max(self._vec)
            for deg in range(min(self._vec), max_deg+1):
                A = ascii_art(self.vector(deg))
                A._baseline = A.height() // 2
                result += arrow_art(deg) + A
            return result + arrow_art(max_deg+1) + AsciiArt(['0'])

        def _add_(self, other):
            """
            Module addition.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x - y,
                ....:        1: H.module(1).basis().an_element(),
                ....:        2: H.module(2).basis().an_element()})
                sage: [a.vector(i) for i in range(3)]
                [F[x] - F[y], F[1] # F[1], F[1] # F[1] # F[1]]
                sage: [H.an_element().vector(i) for i in range(3)]
                [2*F[1] + 2*F[x] + 3*F[y],
                 2*F[1] # F[1] + 2*F[1] # F[x] + 3*F[1] # F[y],
                 2*F[1] # F[1] # F[1] + 2*F[1] # F[1] # F[x] + 3*F[1] # F[1] # F[y]]

                sage: v = a + H.an_element()
                sage: [v.vector(i) for i in range(3)]
                [2*F[1] + 3*F[x] + 2*F[y],
                 3*F[1] # F[1] + 2*F[1] # F[x] + 3*F[1] # F[y],
                 3*F[1] # F[1] # F[1] + 2*F[1] # F[1] # F[x] + 3*F[1] # F[1] # F[y]]
            """
            vectors = dict(self._vec)  # Make a (shallow) copy
            for d in other._vec:
                if d in vectors:
                    vectors[d] += other._vec[d]
                    if not vectors[d]:
                        del vectors[d]
                else:
                    vectors[d] = other._vec
            parent = self.parent()
            return parent.element_class(parent, vectors)

        def _lmul_(self, scalar):
            """
            Scalar multiplication.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x - y,
                ....:        1: H.module(1).basis().an_element(),
                ....:        2: H.module(2).basis().an_element()})
                sage: v = 3*a
                sage: [v.vector(i) for i in range(3)]
                [3*F[x] - 3*F[y], 3*F[1] # F[1], 3*F[1] # F[1] # F[1]]
            """
            if scalar == 0:
                return self.zero()
            vectors = {}
            for d in self._vec:
                vec = scalar * self._vec[d]
                if vec:
                    vectors[d] = vec
            return self.__class__(self.parent(), vectors)

        def _richcmp_(self, other, op):
            """
            Rich comparison of ``self`` to ``other``.

            EXAMPLES::

                sage: F.<x,y> = FreeAlgebra(ZZ)
                sage: H = F.hochschild_complex(F)
                sage: a = H({0: x - y,
                ....:        1: H.module(1).basis().an_element(),
                ....:        2: H.module(2).basis().an_element()})
                sage: a == 3*a
                False
                sage: a + a == 2*a
                True
                sage: a == H.zero()
                False

                sage: a != 3*a
                True
                sage: a + a != 2*a
                False
                sage: a != H.zero()
                True
            """
            return richcmp(self._vec, other._vec, op)
