"""
Classical Lie Algebras

These are the Lie algebras corresponding to types `A_n`, `B_n`, `C_n`,
and `D_n`. We also include support for the exceptional types
`E_{6,7,8}`, `F_4`, and `G_2` in the Chevalley basis, and we
give the matrix representation given in [HRT2000]_.

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
- Sebastian Oehms  (2018-03-18): matrix method of the element class
  of ClassicalMatrixLieAlgebra added
- Travis Scrimshaw (2019-07-09): Implemented compact real form
"""

# ****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from collections import OrderedDict

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.indexed_generators import IndexedGenerators
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.triangular_kac_moody_algebras import TriangularKacMoodyAlgebras

from sage.algebras.lie_algebras.lie_algebra import MatrixLieAlgebraFromAssociative, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
from sage.matrix.matrix_space import MatrixSpace
from sage.sets.family import Family
from sage.modules.free_module import FreeModule


class ClassicalMatrixLieAlgebra(MatrixLieAlgebraFromAssociative):
    """
    A classical Lie algebra represented using matrices.

    This means a classical Lie algebra given as a Lie
    algebra of matrices, with commutator as Lie bracket.

    INPUT:

    - ``R`` -- the base ring
    - ``ct`` -- the finite Cartan type

    EXAMPLES::

        sage: lie_algebras.ClassicalMatrix(QQ, ['A', 4])
        Special linear Lie algebra of rank 5 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, CartanType(['B',4]))
        Special orthogonal Lie algebra of rank 9 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, 'C4')
        Symplectic Lie algebra of rank 8 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, cartan_type=['D',4])
        Special orthogonal Lie algebra of rank 8 over Rational Field
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Return the correct parent based on input.

        EXAMPLES::

            sage: lie_algebras.ClassicalMatrix(QQ, ['A', 4])
            Special linear Lie algebra of rank 5 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, CartanType(['B',4]))
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, 'C4')
            Symplectic Lie algebra of rank 8 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, cartan_type=['D',4])
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        if isinstance(cartan_type, (CartanMatrix, DynkinDiagram_class)):
            cartan_type = cartan_type.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)

        if not cartan_type.is_finite():
            raise ValueError("only for finite types")

        if cartan_type.type() == 'A':
            return sl(R, cartan_type.rank() + 1)
        if cartan_type.type() == 'B':
            return so(R, 2*cartan_type.rank() + 1)
        if cartan_type.type() == 'C':
            return sp(R, 2*cartan_type.rank())
        if cartan_type.type() == 'D':
            return so(R, 2*cartan_type.rank())
        if cartan_type.type() == 'E':
            if cartan_type.rank() == 6:
                return e6(R)
            if cartan_type.rank() == 7:
                return e7(R)
            if cartan_type.rank() == 8:
                return e8(R)
        if cartan_type.type() == 'F' and cartan_type.rank() == 4:
            return f4(R)
        if cartan_type.type() == 'G' and cartan_type.rank() == 2:
            return g2(R)
        raise ValueError("invalid Cartan type")

    def __init__(self, R, ct, e, f, h, sparse=True):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- the base ring
        - ``ct`` -- the Cartan type
        - ``e`` -- the `e` generators
        - ``f`` -- the `f` generators
        - ``h`` -- the `h` generators
        - ``sparse`` -- boolean (default: ``True``); use the sparse vectors
          for the basis computation

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: TestSuite(g).run()

        TESTS:

        Check that :issue:`23266` is fixed::

            sage: sl2 = lie_algebras.sl(QQ, 2, 'matrix')
            sage: isinstance(sl2.indices(), FiniteEnumeratedSet)
            True

        Check that elements are hashable (see :issue:`28961`)::

            sage: sl2 = lie_algebras.sl(QQ, 2, 'matrix')
            sage: e,f,h = list(sl2.basis())
            sage: len(set([e, e+f]))
            2
        """
        n = len(e)
        I = ct.index_set()
        names = ['e%s' % i for i in I]
        names += ['f%s' % i for i in I]
        names += ['h%s' % i for i in I]
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(names)
        MatrixLieAlgebraFromAssociative.__init__(self, e[0].parent(),
                                                 gens=tuple(e + f + h),
                                                 names=tuple(names),
                                                 index_set=index_set,
                                                 category=category)
        self._cartan_type = ct
        self._sparse = sparse

        gens = tuple(self.gens())
        self._e = Family({i: gens[c] for c, i in enumerate(I)})
        self._f = Family({i: gens[n+c] for c, i in enumerate(I)})
        self._h = Family({i: gens[2*n+c] for c, i in enumerate(I)})

    def e(self, i):
        r"""
        Return the generator `e_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.e(2)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._e[i]

    def f(self, i):
        r"""
        Return the generator `f_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.f(2)
            [0 0 0]
            [0 0 0]
            [0 1 0]
        """
        return self._f[i]

    def h(self, i):
        """
        Return the generator `h_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.h(2)
            [ 0  0  0]
            [ 0  1  0]
            [ 0  0 -1]
        """
        return self._h[i]

    @cached_method
    def index_set(self):
        """
        Return the index_set of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.index_set()
            (1, 2)
        """
        return self._cartan_type.index_set()

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def epsilon(self, i, h):
        r"""
        Return the action of the functional
        `\varepsilon_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.epsilon(1, g.h(1))
            1
            sage: g.epsilon(2, g.h(1))
            -1
            sage: g.epsilon(3, g.h(1))
            0
        """
        return h[i-1,i-1]

    # Do we want this to be optional or required?
    # There probably is a generic implementation we can do.
    @abstract_method(optional=True)
    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.simple_root(1, g.h(1))
            2
            sage: g.simple_root(1, g.h(2))
            -1
        """

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.
        If ``pos`` is ``True``, then returns `e_{\theta}`, otherwise it
        returns `f_{\theta}`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.highest_root_basis_elt()
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        RL = self._cartan_type.root_system().root_lattice()
        coroots = RL.simple_coroots()
        theta = RL.highest_root()
        i,w = theta.to_simple_root(True)
        r = RL.simple_root(i)
        if pos:
            gens = self._e
        else:
            gens = self._f
        cur = gens[i]
        for j in reversed(w):
            for _ in range(-r.scalar(coroots[j])):
                cur = self.bracket(gens[j], cur)
            r = r.reflection(coroots[j], True)
        return cur

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: M = LieAlgebra(ZZ, cartan_type=['A',2], representation='matrix')
            sage: list(M.basis())
            [
            [ 1  0  0]  [0 1 0]  [0 0 1]  [0 0 0]  [ 0  0  0]  [0 0 0]  [0 0 0]
            [ 0  0  0]  [0 0 0]  [0 0 0]  [1 0 0]  [ 0  1  0]  [0 0 1]  [0 0 0]
            [ 0  0 -1], [0 0 0], [0 0 0], [0 0 0], [ 0  0 -1], [0 0 0], [1 0 0],
            <BLANKLINE>
            [0 0 0]
            [0 0 0]
            [0 1 0]
            ]

        Sparse version::

            sage: e6 = LieAlgebra(QQ, cartan_type=['E',6], representation='matrix')
            sage: len(e6.basis())  # long time
            78
        """
        # This is a fairly generic method of constructing a basis
        from sage.matrix.constructor import matrix

        R = self.base_ring()
        basis_pivots = set()
        gens = list(self.lie_algebra_generators())
        added = gens
        m = self._assoc.ncols()
        adim = self._assoc.dimension()
        cur_mat = matrix(R, 0, adim, sparse=self._sparse)

        # Helper functions for sparse matrices
        def set_row(mat, row, val):
            for k, v in val.dict().items():
                a, b = k
                mat[row, a*m+b] = v

        def build_assoc(row):
            ret = {}
            for i, v in row.dict().items():
                ret[i//m, i % m] = v
            return self._assoc(ret)

        while added:
            if self._sparse:
                mat = {}
                count = 0
                for x in added:
                    set_row(mat, count, x.value)
                    count += 1
                    for y in gens:
                        ret = x.bracket(y)
                        if ret:
                            set_row(mat, count, ret.value)
                            count += 1
                mat = matrix(R, count, adim, mat, sparse=True)
            else:
                mat = []
                for x in added:
                    mat.append(x.value.list())
                    for y in gens:
                        ret = x.bracket(y)
                        if ret:
                            mat.append(ret.value.list())
                mat = matrix(R, mat)
            cur_mat = cur_mat.stack(mat)
            cur_mat.echelonize()
            pivots = cur_mat.pivots()
            added = []
            if len(pivots) != len(basis_pivots):
                for i,p in enumerate(pivots):
                    if p in basis_pivots:
                        continue
                    basis_pivots.add(p)
                    if self._sparse:
                        added.append(self.element_class(self, build_assoc(cur_mat[i])))
                    else:
                        added.append(self.element_class(self, self._assoc(cur_mat[i].list())))
                cur_mat = cur_mat.submatrix(nrows=len(pivots))
        if self._sparse:
            basis = [self.element_class(self, build_assoc(cur_mat[i]))
                     for i in range(cur_mat.rank())]
        else:
            basis = [self.element_class(self, self._assoc(cur_mat[i].list()))
                     for i in range(cur_mat.rank())]
        return Family(basis)

    def affine(self, kac_moody=True):
        """
        Return the affine (Kac-Moody) Lie algebra of ``self``.

        EXAMPLES::

            sage: so5 = lie_algebras.so(QQ, 5, 'matrix')
            sage: so5
            Special orthogonal Lie algebra of rank 5 over Rational Field
            sage: so5.affine()
            Affine Special orthogonal Kac-Moody algebra of rank 5 over Rational Field
            sage: so5.affine(False)
            Affine Special orthogonal Lie algebra of rank 5 over Rational Field
        """
        from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
        return AffineLieAlgebra(self, kac_moody=kac_moody)


class gl(MatrixLieAlgebraFromAssociative):
    r"""
    The matrix Lie algebra `\mathfrak{gl}_n`.

    The Lie algebra `\mathfrak{gl}_n` which consists of all `n \times n`
    matrices.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: TestSuite(g).run()

        TESTS:

        Check that :issue:`23266` is fixed::

            sage: gl2 = lie_algebras.gl(QQ, 2)
            sage: isinstance(gl2.basis().keys(), FiniteEnumeratedSet)
            True
            sage: Ugl2 = gl2.pbw_basis()
            sage: prod(Ugl2.gens())
            PBW['E_0_0']*PBW['E_0_1']*PBW['E_1_0']*PBW['E_1_1']
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        names = []
        gens = []
        for i in range(n):
            for j in range(n):
                names.append('E_{}_{}'.format(i, j))
                mat = MS({(i, j): one})
                mat.set_immutable()
                gens.append(mat)
        self._n = n
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(names)
        MatrixLieAlgebraFromAssociative.__init__(self, MS, tuple(gens),
                                                 names=tuple(names),
                                                 index_set=index_set,
                                                 category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.gl(QQ, 4)
            General linear Lie algebra of rank 4 over Rational Field
        """
        return "General linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{gl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy) - 2 \mathrm{tr}(x)
            \mathrm{tr}(y).

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: x = g.an_element()
            sage: y = g.gens()[1]
            sage: g.killing_form(x, y)
            8
        """
        return (2 * self._n * (x.value * y.value).trace()
                - 2 * x.value.trace() * y.value.trace())

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 2)
            sage: tuple(g.basis())
            (
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            )
        """
        G = self.gens()
        return Family(self._indices, lambda i: G[self._indices.index(i)])

    def monomial(self, i):
        r"""
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: gl4 = lie_algebras.gl(QQ, 4)
            sage: gl4.monomial('E_2_1')
            [0 0 0 0]
            [0 0 0 0]
            [0 1 0 0]
            [0 0 0 0]
            sage: gl4.monomial((2,1))
            [0 0 0 0]
            [0 0 0 0]
            [0 1 0 0]
            [0 0 0 0]
        """
        if isinstance(i, tuple):
            return self.basis()['E_{}_{}'.format(*i)]
        return self.basis()[i]

    class Element(MatrixLieAlgebraFromAssociative.Element):
        def monomial_coefficients(self, copy=True):
            r"""
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: gl4 = lie_algebras.gl(QQ, 4)
                sage: x = gl4.monomial('E_2_1') + 3*gl4.monomial('E_0_3')
                sage: x.monomial_coefficients()
                {'E_0_3': 3, 'E_2_1': 1}
            """
            d = {}
            for k in self.value.nonzero_positions():
                d['E_{}_{}'.format(*k)] = self.value[k]
            return d


class sl(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{sl}_n`.

    The Lie algebra `\mathfrak{sl}_n`, which consists of all `n \times n`
    matrices with trace 0. This is the Lie algebra of type `A_{n-1}`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        e = [MS({(i,i+1):one}) for i in range(n-1)]
        f = [MS({(i+1,i):one}) for i in range(n-1)]
        h = [MS({(i,i):one, (i+1,i+1):-one}) for i in range(n-1)]
        self._n = n
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['A', n-1]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sl(QQ, 5, representation='matrix')
            Special linear Lie algebra of rank 5 over Rational Field
        """
        return "Special linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            10
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        return h[i,i] - h[i+1,i+1]


class so(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{so}_n`.

    The Lie algebra `\mathfrak{so}_n`, which is isomorphic to the
    Lie algebra of all anti-symmetric `n \times n` matrices.
    The implementation here uses a different bilinear form and follows
    the description in Chapter 8 of [HK2002]_. More precisely, this
    is the set of matrices:

    .. MATH::

        \begin{pmatrix}
        A & B \\
        C & D
        \end{pmatrix}

    such that `A^t = -D`, `B^t = -B`, `C^t = -C` for `n` even and

    .. MATH::

        \begin{pmatrix}
        A & B & a \\
        C & D & b \\
        c & d & 0
        \end{pmatrix}

    such that `A^t = -D`, `B^t = -B`, `C^t = -C`, `a^t = -d`,
    and `b^t = -c` for `n` odd.

    This is the Lie algebra of type `B_{(n-1)/2}` or `D_{n/2}` if `n`
    is odd or even respectively.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n)
        one = R.one()
        self._n = n
        if n % 2 == 0:  # Even
            m = n // 2 - 1  # -1 for indexing
            n -= 1
            e = [MS({(m-1, n): one, (m, n-1): -one})]
            f = [MS({(n, m-1): one, (n-1, m): -one})]
            h = [MS({(m-1, m-1): one, (m, m): one, (n-1, n-1): -one, (n, n): -one})]
            m += 1
            ct = CartanType(['D', m])
        else:  # Odd
            m = (n-1) // 2 - 1  # -1 for indexing
            n -= 1
            e = [MS({(m, n): 2, (n, n-1): -2})]
            f = [MS({(n, m): one, (n-1, n): -one})]
            h = [MS({(m, m): 2, (n-1, n-1): -2})]
            m += 1
            ct = CartanType(['B', m])
        e = [MS({(i, i+1): one, (m+i+1, m+i): -one}) for i in range(m-1)] + e
        f = [MS({(i+1, i): one, (m+i, m+i+1): -one}) for i in range(m-1)] + f
        h = [MS({(i, i): one, (i+1, i+1): -one, (m+i, m+i): -one, (m+i+1, m+i+1): one})
             for i in range(m-1)] + h
        ClassicalMatrixLieAlgebra.__init__(self, R, ct, e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['B', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: LieAlgebra(QQ, cartan_type=['D', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        return "Special orthogonal Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{so}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (n - 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            12
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            14
        """
        return (self._n - 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES:

        The even or type `D` case::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1 -1]
            [ 0 -1  2  0]
            [ 0 -1  0  2]

        The odd or type `B` case::

            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -2  2]
        """
        i = self.index_set().index(i)
        if i == len(self.index_set()) - 1:
            if self._n % 2 == 0:
                return h[i-1, i-1] + h[i, i]
            # otherwise we are odd
            return h[i, i]
        return h[i, i] - h[i+1, i+1]


class sp(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{sp}_n`.

    The Lie algebra `\mathfrak{sp}_{2k}`, which consists of all
    `2k \times 2k` matrices `X` that satisfy the equation:

    .. MATH::

        X^T M - M X = 0

    where

    .. MATH::

        M = \begin{pmatrix}
        0 & I_k \\
        -I_k & 0
        \end{pmatrix}.

    This is the Lie algebra of type `C_k`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        self._n = n
        n = n // 2
        e = [MS({(i,i+1):one, (n+i+1,n+i):-one}) for i in range(n-1)]
        e.append(MS({(n-1,2*n-1):one})) # -1 for indexing
        f = [MS({(i+1,i):one, (n+i,n+i+1):-one}) for i in range(n-1)]
        f.append(MS({(2*n-1,n-1):one})) # -1 for indexing
        h = [MS({(i,i):one, (i+1,i+1):-one, (n+i,n+i):-one, (n+i+1,n+i+1):one}) for i in range(n-1)]
        h.append(MS({(n-1,n-1):one, (2*n-1,2*n-1):-one})) # -1 for indexing
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['C', n]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sp(QQ, 8, representation='matrix')
            Symplectic Lie algebra of rank 8 over Rational Field
        """
        return "Symplectic Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sp}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (2n + 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            36
        """
        return (2 * self._n + 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        if i == self._n / 2 - 1:
            return 2*h[i,i]
        return h[i,i] - h[i+1,i+1]


class ExceptionalMatrixLieAlgebra(ClassicalMatrixLieAlgebra):
    """
    A matrix Lie algebra of exceptional type.
    """
    def __init__(self, R, cartan_type, e, f, h=None, sparse=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E',6], representation='matrix')
            sage: all(g.h(i) == g.e(i).bracket(g.f(i)) for i in range(1,7))
            True
        """
        if h is None:
            h = [e[i] * f[i] - f[i] * e[i] for i in range(len(e))]
        ClassicalMatrixLieAlgebra.__init__(self, R, cartan_type, e, f, h, sparse=sparse)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['G',2], representation='matrix')
            Simple matrix Lie algebra of type ['G', 2] over Rational Field
        """
        return "Simple matrix Lie algebra of type {} over {}".format(self.cartan_type(), self.base_ring())


class e6(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{e}_6`.

    The simple Lie algebra `\mathfrak{e}_6` of type `E_6`. The matrix
    representation is given following [HRT2000]_.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E',6], representation='matrix')
            sage: TestSuite(g).run()  # long time
        """
        MS = MatrixSpace(R, 27, sparse=True)
        one = R.one()
        coords = [[(0,1), (10,12), (13,15), (16,17), (18,19), (20,21)],
                  [(3,4), (5,6), (7,9), (18,20), (19,21), (22,23)],
                  [(1,2), (8,10), (11,13), (14,16), (19,22), (21,23)],
                  [(2,3), (6,8), (9,11), (16,18), (17,19), (23,24)],
                  [(3,5), (4,6), (11,14), (13,16), (15,17), (24,25)],
                  [(5,7), (6,9), (8,11), (10,13), (12,15), (25,26)]]
        e = [MS({c: one for c in coord}) for coord in coords]
        f = [MS({(c[1],c[0]): one for c in coord}) for coord in coords]
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['E', 6]), e, f)


class e7(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{e}_7`.

    The simple Lie algebra `\mathfrak{e}_7` of type `E_7`. The matrix
    representation is given following [HRT2000]_.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 7], representation='matrix')
            sage: g
            Simple matrix Lie algebra of type ['E', 7] over Rational Field

            sage: len(g.basis())  # long time
            133
            sage: TestSuite(g).run()  # long time
        """
        MS = MatrixSpace(R, 56, sparse=True)
        one = R.one()
        coords = [[(6,7), (8,9), (10,11), (12,14), (15,17), (18,21), (34,37), (38,40), (41,43), (44,45), (46,47), (48,49)],
                  [(4,5), (6,8), (7,9), (19,22), (23,25), (26,28), (27,29), (30,32), (33,36), (46,48), (47,49), (50,51)],
                  [(4,6), (5,8), (11,13), (14,16), (17,20), (21,24), (31,34), (35,38), (39,41), (42,44), (47,50), (49,51)],
                  [(3,4), (8,10), (9,11), (16,19), (20,23), (24,27), (28,31), (32,35), (36,39), (44,46), (45,47), (51,52)],
                  [(2,3), (10,12), (11,14), (13,16), (23,26), (25,28), (27,30), (29,32), (39,42), (41,44), (43,45), (52,53)],
                  [(1,2), (12,15), (14,17), (16,20), (19,23), (22,25), (30,33), (32,36), (35,39), (38,41), (40,43), (53,54)],
                  [(0,1), (15,18), (17,21), (20,24), (23,27), (25,29), (26,30), (28,32), (31,35), (34,38), (37,40), (54,55)]]
        e = [MS({c: one for c in coord}) for coord in coords]
        f = [MS({(c[1], c[0]): one for c in coord}) for coord in coords]
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['E', 7]), e, f)


class e8(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{e}_8`.

    The simple Lie algebra `\mathfrak{e}_8` of type `E_8` built from the
    adjoint representation in the Chevalley basis.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 8], representation='matrix')  # long time
            sage: g  # long time
            Simple matrix Lie algebra of type ['E', 8] over Rational Field

        We skip the not implemented methods test as it takes too much time::

            sage: TestSuite(g).run(skip='_test_not_implemented_methods')  # long time
        """
        ct = CartanType(['E', 8])
        g = LieAlgebraChevalleyBasis(R, ct)
        e = [ge.adjoint_matrix(sparse=True) for ge in g.e()]
        f = [gf.adjoint_matrix(sparse=True) for gf in g.f()]
        ExceptionalMatrixLieAlgebra.__init__(self, R, ct, e, f)

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 8], representation='matrix')  # long time
            sage: len(g.basis())  # long time
            248
        """
        g = LieAlgebraChevalleyBasis(self.base_ring(), self.cartan_type())
        return Family([ge.adjoint_matrix(sparse=True) for ge in g.basis()])


class f4(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{f}_4`.

    The simple Lie algebra `\mathfrak{f}_f` of type `F_4`. The matrix
    representation is given following [HRT2000]_ but indexed in the
    reversed order (i.e., interchange 1 with 4 and 2 with 3).
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['F',4], representation='matrix')
            sage: TestSuite(g).run()  # long time
        """
        MS = MatrixSpace(R, 26, sparse=True)
        one = R.one()

        coords = [[(0,1), (5,7), (6,9), (8,11), (10,12), (10,13), (12,14),
                   (15,16), (17,18), (19,20), (24,25)],
                  [(1,2), (3,5), (4,6), (8,10), (11,12), (11,13), (13,15),
                   (14,16), (18,21), (20,22), (23,24)],
                  [(2,3), (6,8), (9,11), (15,17), (16,18), (22,23)],
                  [(3,4), (5,6), (7,9), (17,19), (18,20), (21,22)]]
        e = [MS({c: one for c in coord}) for coord in coords]
        # Double (10, 12) in e1 and (11,13) in e2
        e[0][10,12] = 2*one
        e[1][11,13] = 2*one

        coords = [[(1,0), (7,5), (9,6), (11,8), (12,10), (14,12), (14,13),
                   (16,15), (18,17), (20,19), (25,24)],
                  [(2,1), (5,3), (6,4), (10,8), (13,11), (15,12), (15,13),
                   (16,14), (21,18), (22,20), (24,23)],
                  [(3,2), (8,6), (11,9), (17,15), (18,16), (23,22)],
                  [(4,3), (6,5), (9,7), (19,17), (20,18), (22,21)]]
        f = [MS({c: one for c in coord}) for coord in coords]
        # Double (14, 12) in f1 and (15,13) in f2
        f[0][14,12] = 2*one
        f[1][15,13] = 2*one

        # Our Cartan matrix convention is dual to that of [HRT2000]_
        e.reverse()
        f.reverse()
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['F', 4]), e, f)


class g2(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{g}_2`.

    The simple Lie algebra `\mathfrak{g}_2` of type `G_2`. The matrix
    representation is given following [HRT2000]_.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G',2], representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, 7, sparse=True)
        one = R.one()
        e = [MS({(0,1): one, (2,3): 2*one, (3,4): one, (5,6): one}),
             MS({(1,2): one, (4,5): one})]
        f = [MS({(1,0): one, (3,2): one, (4,3): 2*one, (6,5): one}),
             MS({(2,1): one, (5,4): one})]
        h = [MS({(0,0): one, (1,1): -one, (2,2): 2*one, (4,4): -2*one, (5,5): one, (6,6): -one}),
             MS({(1,1): one, (2,2): -one, (4,4): one, (5,5): -one})]
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['G', 2]), e, f, h)


#######################################
# Compact real form

class MatrixCompactRealForm(FinitelyGeneratedLieAlgebra):
    r"""
    The compact real form of a matrix Lie algebra.

    Let `L` be a classical (i.e., type `ABCD`) Lie algebra over `\RR`
    given as matrices that is invariant under matrix transpose (i.e.,
    `X^T \in L` for all `X \in L`). Then we can perform the
    *Cartan decomposition* of `L` by `L = K \oplus S`, where `K`
    (resp. `S`) is the set of skew-symmetric (resp. symmetric) matrices
    in `L`. Then the Lie algebra `U = K \oplus i S` is an `\RR`-subspace
    of the complexification of `L` that is closed under commutators and
    has skew-hermitian matrices. Hence, the Killing form is negative
    definitive (i.e., `U` is a compact Lie algebra), and thus `U` is
    the complex real form of the complexification of `L`.

    EXAMPLES::

        sage: U = LieAlgebra(QQ, cartan_type=['A',1], representation="compact real")
        sage: list(U.basis())
        [
        [ 0  1]  [ i  0]  [0 i]
        [-1  0], [ 0 -i], [i 0]
        ]
        sage: U.killing_form_matrix()
        [-8  0  0]
        [ 0 -8  0]
        [ 0  0 -8]

    Computations are only (currently) possible if this is defined
    over a field::

        sage: U = LieAlgebra(ZZ, cartan_type=['A',1], representation="compact real")
        sage: list(U.basis())
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    def __init__(self, R, cartan_type):
        """
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation="compact real")
            sage: TestSuite(L).run()
        """
        if not cartan_type.is_finite():
            raise ValueError("the Cartan type must be finite type")

        self._classical = ClassicalMatrixLieAlgebra(R, cartan_type)
        self._MS = self._classical._assoc
        dim = self._classical.dimension()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(range(dim))
        names = tuple(['CR%s' % s for s in range(dim)])
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names=names,
                                             index_set=index_set,
                                             category=category)

    @cached_method
    def basis(self):
        """
        Compute a basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['B',2], representation="compact real")
            sage: list(L.basis())
            [
            [ 0  1  0  0  0]  [ 0  0  0  1  0]  [ 0  0  0  0  1]  [ 0  0  0  0  0]
            [-1  0  0  0  0]  [ 0  0 -1  0  0]  [ 0  0  0  0  0]  [ 0  0  0  0  1]
            [ 0  0  0  1  0]  [ 0  1  0  0  0]  [ 0  0  0  0  1]  [ 0  0  0  0  0]
            [ 0  0 -1  0  0]  [-1  0  0  0  0]  [ 0  0  0  0  0]  [ 0  0  0  0  1]
            [ 0  0  0  0  0], [ 0  0  0  0  0], [-1  0 -1  0  0], [ 0 -1  0 -1  0],
            <BLANKLINE>
            [ i  0  0  0  0]  [ 0  i  0  0  0]  [ 0  0  0  i  0]  [ 0  0  0  0  i]
            [ 0  0  0  0  0]  [ i  0  0  0  0]  [ 0  0 -i  0  0]  [ 0  0  0  0  0]
            [ 0  0 -i  0  0]  [ 0  0  0 -i  0]  [ 0 -i  0  0  0]  [ 0  0  0  0 -i]
            [ 0  0  0  0  0]  [ 0  0 -i  0  0]  [ i  0  0  0  0]  [ 0  0  0  0  0]
            [ 0  0  0  0  0], [ 0  0  0  0  0], [ 0  0  0  0  0], [ i  0 -i  0  0],
            <BLANKLINE>
            [ 0  0  0  0  0]  [ 0  0  0  0  0]
            [ 0  i  0  0  0]  [ 0  0  0  0  i]
            [ 0  0  0  0  0]  [ 0  0  0  0  0]
            [ 0  0  0 -i  0]  [ 0  0  0  0 -i]
            [ 0  0  0  0  0], [ 0  i  0 -i  0]
            ]
        """
        from sage.matrix.constructor import matrix
        zero = self._MS.zero()
        basis = self._classical.basis()
        R = self.base_ring()
        mat = matrix(R, [((b.value - b.value.transpose()) / 2).list() for b in basis],
                     sparse=self._MS.is_sparse())
        mat.echelonize()
        ret = [self.element_class(self, self._MS(mat[i].list()), zero)
               for i in range(mat.rank())]
        mat = matrix(R, [((b.value + b.value.transpose()) / 2).list() for b in basis],
                     sparse=self._MS.is_sparse())
        mat.echelonize()
        ret += [self.element_class(self, zero, self._MS(mat[i].list()))
                for i in range(mat.rank())]
        return Family(ret)

    @cached_method
    def zero(self):
        """
        Return the element `0`.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['D',4], representation="compact real")
            sage: L.zero()
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
        """
        return self.element_class(self, self._MS.zero(), self._MS.zero())

    def monomial(self, i):
        """
        Return the monomial indexed by ``i``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',3], representation="compact real")
            sage: L.monomial(0)
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  0]
            [ 0  0  0  0]
        """
        return self.basis()[i]

    def term(self, i, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['C',3], representation="compact real")
            sage: L.term(4, 7/2)
            [   0    0    0    0    0  7/2]
            [   0    0    0    0    0    0]
            [   0    0    0  7/2    0    0]
            [   0    0 -7/2    0    0    0]
            [   0    0    0    0    0    0]
            [-7/2    0    0    0    0    0]
        """
        if c is None:
            c = self.base_ring().one()
        else:
            c = self.base_ring()(c)
        return c * self.basis()[i]

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',1], representation="compact real")
            sage: L._repr_option("element_ascii_art")
            True
        """
        if key == "element_ascii_art":
            return True
        return FinitelyGeneratedLieAlgebra._repr_option(self, key)

    class Element(Element):
        """
        An element of a matrix Lie algebra in its compact real form.
        """
        def __init__(self, parent, real, imag):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['D',4], representation="compact real")
                sage: TestSuite(L.an_element()).run()
            """
            Element.__init__(self, parent)
            self._real = real
            self._imag = imag
            self._real.set_immutable()
            self._imag.set_immutable()
            self._mc = None

        def _combined_matrix(self):
            r"""
            Return a single matrix representative of ``self``.

            .. NOTE::

                The resulting base ring is `R[i]`, where `R` is the
                base ring of the Lie algebra.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation="compact real")
                sage: x = L.sum((i+1)/7*b for i,b in enumerate(L.basis()))
                sage: M = x._combined_matrix()
                sage: M
                [      4/7*i 5/7*i + 1/7 6/7*i + 2/7]
                [5/7*i - 1/7           i 8/7*i + 3/7]
                [6/7*i - 2/7 8/7*i - 3/7     -11/7*i]
                sage: M.parent()
                Full MatrixSpace of 3 by 3 sparse matrices over
                 Univariate Polynomial Ring in i over Rational Field
            """
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            MS = self.parent()._MS
            R = PolynomialRing(MS.base_ring(), 'i')
            return self._real + R.gen() * self._imag

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['B',2], representation="compact real")
                sage: L.sum((i+1)/7*b for i,b in enumerate(L.basis()))
                [        5/7*i   6/7*i + 1/7             0       i + 2/7   8/7*i + 3/7]
                [  6/7*i - 1/7         9/7*i      -i - 2/7             0  10/7*i + 4/7]
                [            0      -i + 2/7        -5/7*i  -6/7*i + 1/7  -8/7*i + 3/7]
                [      i - 2/7             0  -6/7*i - 1/7        -9/7*i -10/7*i + 4/7]
                [  8/7*i - 3/7  10/7*i - 4/7  -8/7*i - 3/7 -10/7*i - 4/7             0]
            """
            return repr(self._combined_matrix())

        def _latex_(self):
            r"""
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['B',2], representation="compact real")
                sage: x = L.sum((i+1)/7*b for i,b in enumerate(L.basis()))
                sage: latex(x)
                \left(\begin{array}{rrrrr}
                \frac{5}{7} i & \frac{6}{7} i + \frac{1}{7} & 0 & i + \frac{2}{7} & \frac{8}{7} i + \frac{3}{7} \\
                \frac{6}{7} i - \frac{1}{7} & \frac{9}{7} i & -i - \frac{2}{7} & 0 & \frac{10}{7} i + \frac{4}{7} \\
                0 & -i + \frac{2}{7} & -\frac{5}{7} i & -\frac{6}{7} i + \frac{1}{7} & -\frac{8}{7} i + \frac{3}{7} \\
                i - \frac{2}{7} & 0 & -\frac{6}{7} i - \frac{1}{7} & -\frac{9}{7} i & -\frac{10}{7} i + \frac{4}{7} \\
                \frac{8}{7} i - \frac{3}{7} & \frac{10}{7} i - \frac{4}{7} & -\frac{8}{7} i - \frac{3}{7} & -\frac{10}{7} i - \frac{4}{7} & 0
                \end{array}\right)
            """
            from sage.misc.latex import latex
            return latex(self._combined_matrix())

        def _ascii_art_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation="compact real")
                sage: x = L.sum((i+1)/7*b for i,b in enumerate(L.basis()))
                sage: ascii_art(x)
                [      4/7*i 5/7*i + 1/7 6/7*i + 2/7]
                [5/7*i - 1/7           i 8/7*i + 3/7]
                [6/7*i - 2/7 8/7*i - 3/7     -11/7*i]
            """
            from sage.typeset.ascii_art import ascii_art
            return ascii_art(self._combined_matrix())

        def _unicode_art_(self):
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation="compact real")
                sage: x = L.sum((i+1)/7*b for i,b in enumerate(L.basis()))
                sage: unicode_art(x)
                ⎛      4/7*i 5/7*i + 1/7 6/7*i + 2/7⎞
                ⎜5/7*i - 1/7           i 8/7*i + 3/7⎟
                ⎝6/7*i - 2/7 8/7*i - 3/7     -11/7*i⎠
            """
            from sage.typeset.unicode_art import unicode_art
            return unicode_art(self._combined_matrix())

        def __bool__(self) -> bool:
            r"""
            Return if ``self`` is nonzero.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['C',3], representation="compact real")
                sage: all(b for b in L.basis() if b != 0)
                True
                sage: bool(L.zero())
                False
            """
            return bool(self._real) or bool(self._imag)

        def __hash__(self):
            r"""
            Return the hash of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',2], representation="compact real")
                sage: x = L.an_element()
                sage: hash(x) == hash((x._real, x._imag))
                True
            """
            return hash((self._real, self._imag))

        def _richcmp_(self, other, op):
            r"""
            Return the richcmp of ``self`` and ``other`` by ``op``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',1], representation="compact real")
                sage: sorted(L.basis())
                [
                [0 i]  [ i  0]  [ 0  1]
                [i 0], [ 0 -i], [-1  0]
                ]
            """
            return richcmp((self._real, self._imag), (other._real, other._imag), op)

        def _add_(self, other):
            r"""
            Add ``self`` and ``other``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['C',2], representation="compact real")
                sage: B = L.basis()
                sage: B[0] + B[6]
                [ 0  1  i  0]
                [-1  0  0  0]
                [ i  0  0  1]
                [ 0  0 -1  0]
                sage: L.sum(B)
                [     i  i + 1  i + 1  i + 1]
                [ i - 1      i  i + 1  i + 1]
                [ i - 1  i - 1     -i -i + 1]
                [ i - 1  i - 1 -i - 1     -i]
            """
            P = self.parent()
            return P.element_class(P, self._real + other._real,
                                   self._imag + other._imag)

        def _sub_(self, other):
            r"""
            Subtract ``self`` and ``other``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['C',2], representation="compact real")
                sage: B = L.basis()
                sage: B[0] - B[6]
                [ 0  1 -i  0]
                [-1  0  0  0]
                [-i  0  0  1]
                [ 0  0 -1  0]
                sage: all(b - b == L.zero() for b in B)
                True
            """
            P = self.parent()
            return P.element_class(P, self._real - other._real,
                                   self._imag - other._imag)

        def _neg_(self):
            r"""
            Negate ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['C',2], representation="compact real")
                sage: B = L.basis()
                sage: -(B[0] + B[6])
                [ 0 -1 -i  0]
                [ 1  0  0  0]
                [-i  0  0 -1]
                [ 0  0  1  0]
                sage: all(-(-b) == b for b in B)
                True
            """
            P = self.parent()
            return P.element_class(P, -self._real, -self._imag)

        def _bracket_(self, other):
            r"""
            Return the Lie bracket of ``self`` and ``other``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['A',1], representation="compact real")
                sage: B = L.basis()
                sage: list(B)
                [
                [ 0  1]  [ i  0]  [0 i]
                [-1  0], [ 0 -i], [i 0]
                ]
                sage: [b._bracket_(bp) for b in B for bp in B]
                [
                [0 0]  [   0 -2*i]  [ 2*i    0]  [  0 2*i]  [0 0]  [ 0 -2]
                [0 0], [-2*i    0], [   0 -2*i], [2*i   0], [0 0], [ 2  0],
                <BLANKLINE>
                [-2*i    0]  [ 0  2]  [0 0]
                [   0  2*i], [-2  0], [0 0]
                ]
            """
            A, B = self._real, self._imag
            X, Y = other._real, other._imag
            P = self.parent()
            return P.element_class(P, A*X - X*A - B*Y + Y*B,
                                   A*Y - Y*A + B*X - X*B)

        def _acted_upon_(self, x, self_on_left):
            r"""
            Return the action of ``x`` on ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['D',4], representation="compact real")
                sage: B = L.basis()
                sage: (3/5) * B[21]
                [     0      0      0      0      0      0      0      0]
                [     0      0      0  3/5*i      0      0      0      0]
                [     0      0      0      0      0      0      0      0]
                [     0  3/5*i      0      0      0      0      0      0]
                [     0      0      0      0      0      0      0      0]
                [     0      0      0      0      0      0      0 -3/5*i]
                [     0      0      0      0      0      0      0      0]
                [     0      0      0      0      0 -3/5*i      0      0]
                sage: B[7] * 7
                [ 0  0  0  0  0  0  0  0]
                [ 0  0  0  7  0  0  0  0]
                [ 0  0  0  0  0  0  0  0]
                [ 0 -7  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  7]
                [ 0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0 -7  0  0]
            """
            P = self.parent()
            return P.element_class(P, x*self._real, x*self._imag)

        def monomial_coefficients(self, copy=False):
            """
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, cartan_type=['C',3], representation="compact real")
                sage: B = L.basis()
                sage: x = L.sum(i*B[i] for i in range(len(B)))
                sage: x.monomial_coefficients() == {i: i for i in range(1,len(B))}
                True
            """
            if self._mc is None:
                P = self.parent()
                B = [b._real.list() + b._imag.list() for b in P.basis()]
                B.append(self._real.list() + self._imag.list())
                R = self.base_ring()
                F = FreeModule(R, len(B[0]))
                dep = list(F.linear_dependence([F(b) for b in B])[0])
                last = dep.pop()
                self._mc = {i: R(-val / last) for i,val in enumerate(dep) if val != 0}
            if copy:
                return dict(self._mc)
            return self._mc


#######################################
# Chevalley Basis

class LieAlgebraChevalleyBasis(LieAlgebraWithStructureCoefficients):
    r"""
    A simple finite dimensional Lie algebra in the Chevalley basis.

    Let `L` be a simple (complex) Lie algebra with roots `\Phi`, then the
    Chevalley basis is given by `e_{\alpha}` for all `\alpha \in \Phi` and
    `h_{\alpha_i} := h_i` where `\alpha_i` is a simple root subject. These
    generators are subject to the relations:

    .. MATH::

        \begin{aligned}
        [h_i, h_j] & = 0,
        \\ [h_i, e_{\beta}] & = A_{\alpha_i, \beta} e_{\beta},
        \\ [e_{\beta}, e_{-\beta}] & = \sum_i A_{\beta, \alpha_i} h_i,
        \\ [e_{\beta}, e_{\gamma}] & = \begin{cases}
        N_{\beta,\gamma} e_{\beta + \gamma} & \beta + \gamma \in \Phi, \\
        0 & \text{otherwise,} \end{cases}
        \end{aligned}

    where `A_{\alpha, \beta} = \frac{2 (\alpha, \beta)}{(\alpha, \alpha)}`
    and `N_{\alpha, \beta}` is the maximum such that
    `\alpha - N_{\alpha, \beta} \beta \in \Phi`.

    For computing the signs of the coefficients, see Section 3 of [CMT2003]_.

    .. SEEALSO::

        For simply-laced types, an alternative construction using an asymmetry
        function is given by :class:`LieAlgebraChevalleyBasis_simply_laced`.
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type, epsilon=None):
        """
        Normalize ``self`` to ensure a unique representation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L2 = LieAlgebra(QQ, cartan_type=CartanType(['A', 2]))
            sage: L3 = LieAlgebra(QQ, cartan_type=CartanMatrix(['A', 2]))
            sage: L1 is L2 and L2 is L3
            True


            sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(1,2)])
            sage: type(L)
            <class '...LieAlgebraChevalleyBasis_simply_laced_with_category'>
            sage: L = LieAlgebra(QQ, cartan_type=['A', 1], epsilon=[])
            sage: type(L)
            <class '...LieAlgebraChevalleyBasis_simply_laced_with_category'>

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(2,3)])
            Traceback (most recent call last):
            ...
            ValueError: not a valid Dynkin orientation
            sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(1,2), (2,1)])
            Traceback (most recent call last):
            ...
            ValueError: not a valid Dynkin orientation
            sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(1,2), (1,1)])
            Traceback (most recent call last):
            ...
            ValueError: not a valid Dynkin orientation
            sage: L = LieAlgebra(QQ, cartan_type=['A', 1], epsilon=[(1,1)])
            Traceback (most recent call last):
            ...
            ValueError: not a valid Dynkin orientation
        """
        if isinstance(cartan_type, (CartanMatrix, DynkinDiagram_class)):
            cartan_type = cartan_type.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
        if epsilon is not None:
            if not cartan_type.is_simply_laced():
                raise ValueError("the Cartan type must be simply-laced with an asymmetry function")
            epsilon = frozenset([tuple(p) for p in epsilon])
            if cartan_type.rank() == 1:
                if epsilon:
                    raise ValueError("not a valid Dynkin orientation")
            else:
                from sage.graphs.graph import Graph
                G = Graph(epsilon, multiedges=True, loops=True, format='list_of_edges')
                if (G.has_multiple_edges() or G.has_loops()
                    or cartan_type.dynkin_diagram().to_undirected() != G.to_simple()):
                    raise ValueError("not a valid Dynkin orientation")
            return LieAlgebraChevalleyBasis_simply_laced(R, cartan_type, epsilon)
        return super().__classcall__(cls, R, cartan_type)

    def __init__(self, R, cartan_type):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2])
            sage: TestSuite(L).run()  # long time
        """
        self._cartan_type = cartan_type
        self._Q = cartan_type.root_system().root_lattice()
        p_roots = list(self._Q.positive_roots_by_height())
        n_roots = [-x for x in p_roots]
        self._p_roots_index = OrderedDict((al, i)
                                          for i, al in enumerate(p_roots))

        alphacheck = self._Q.simple_coroots()
        # We pass p_roots and n_roots so we don't have to reconstruct them
        s_coeffs = self._construct_struct_coeffs(R, p_roots)

        # Make sure a < b for all (a, b) in the coefficients and flip if necessary
        for k in list(s_coeffs):
            a, b = k
            if self._basis_key(a) > self._basis_key(b):
                s_coeffs[(b, a)] = [(index, -v) for index, v in s_coeffs[k].items()]
                del s_coeffs[k]
            else:
                s_coeffs[k] = s_coeffs[k].items()

        I = self._cartan_type.index_set()
        names = ['e{}'.format(i) for i in I]
        names += ['f{}'.format(i) for i in I]
        names += ['h{}'.format(i) for i in I]
        category = TriangularKacMoodyAlgebras(R).FiniteDimensional()
        index_set = p_roots + list(alphacheck) + n_roots
        self._cartan_indices = range(len(p_roots), len(p_roots) + len(alphacheck))
        names = tuple(names)
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(index_set)
        LieAlgebraWithStructureCoefficients.__init__(self, R, s_coeffs, names, index_set,
                                                     category, prefix='E', bracket='[',
                                                     sorting_key=self._basis_key)

    def _construct_struct_coeffs(self, R, p_roots):
        """
        Construct the structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',1], epsilon=[])
            sage: sorted(L._construct_struct_coeffs(QQ, list(L._Q.positive_roots())).items(), key=str)
            [((alpha[1], -alpha[1]), {alphacheck[1]: 1}),
             ((alpha[1], alphacheck[1]), {alpha[1]: -2}),
             ((alphacheck[1], -alpha[1]), {-alpha[1]: -2})]

        TESTS:

        Check that we can construct Lie algebras over positive characteristic
        fields (:issue:`37773`)::

            sage: sl3 = LieAlgebra(GF(3), cartan_type=['A',2])
            sage: sl3.center().basis()
            Family (2*h1 + h2,)
            sage: sl4 = lie_algebras.sl(GF(3), 4)
            sage: sl4.center().dimension()
            0
            sage: sl4.is_nilpotent()
            False
            sage: sl4.lower_central_series()
            (Lie algebra of ['A', 3] in the Chevalley basis,)
            sage: sl4.is_solvable()
            False
            sage: sl4.is_semisimple()
            True
            sage: sl4.killing_form_matrix().det()
            2
            sage: sl5 = lie_algebras.sl(GF(3), 5)
            sage: sl5.killing_form_matrix().det()
            2

        This also includes characteristic 2::

            sage: sl4 = LieAlgebra(GF(2), cartan_type=['A',3])
            sage: sl4.center().basis()
            Family (h1 + h3,)
            sage: sp6 = LieAlgebra(GF(2), cartan_type=['C',3])
            sage: sp6.killing_form_matrix().det()
            0
            sage: F4 = LieAlgebra(GF(2), cartan_type=['F',4])
            sage: F4.killing_form_matrix().det()  # long time
            0
            sage: G2 = LieAlgebra(GF(2), cartan_type=['G',2])
            sage: G2.killing_form_matrix().det()
            0
        """
        alphacheck = self._Q.simple_coroots()
        roots = frozenset(self._Q.roots())
        # We do everything initially over QQ and then convert to R at the end
        #   since this is a ZZ-basis.
        from sage.rings.rational_field import QQ
        one = QQ.one()

        # Determine the signs for the structure coefficients from the root system
        # We first create the special roots
        sp_sign = {}
        for i,a in enumerate(p_roots):
            for b in p_roots[i+1:]:
                if a + b not in p_roots:
                    continue

                # Compute the sign for the extra special pair
                x, y = (a + b).extraspecial_pair()

                if (x, y) == (a, b): # If it already is an extra special pair
                    if (x, y) not in sp_sign:
                        # This swap is so the structure coefficients match with GAP
                        if (sum(x.coefficients()) == sum(y.coefficients())
                            and str(x) > str(y)):
                            y,x = x,y
                        sp_sign[(x, y)] = -one
                        sp_sign[(y, x)] = one
                    continue

                if b - x in roots:
                    t1 = ((b-x).norm_squared() / b.norm_squared()
                          * sp_sign[(x, b-x)] * sp_sign[(a, y-a)])
                else:
                    t1 = 0
                if a - x in roots:
                    t2 = ((a-x).norm_squared() / a.norm_squared()
                          * sp_sign[(x, a-x)] * sp_sign[(b, y-b)])
                else:
                    t2 = 0

                if t1 - t2 > 0:
                    sp_sign[(a,b)] = -one
                elif t2 - t1 > 0:
                    sp_sign[(a,b)] = one
                sp_sign[(b,a)] = -sp_sign[(a,b)]

        # Function to construct the structure coefficients (up to sign)
        def e_coeff(r, s):
            p = 1
            while r - p*s in roots:
                p += 1
            return p

        # Now we can compute all necessary structure coefficients
        s_coeffs = {}
        for i,r in enumerate(p_roots):
            # [e_r, h_i] and [h_i, f_r]
            for ac in alphacheck:
                c = R(r.scalar(ac))
                if c == 0:
                    continue
                s_coeffs[(r, ac)] = {r: -c}
                s_coeffs[(ac, -r)] = {-r: -c}

            # [e_r, f_r]
            s_coeffs[(r, -r)] = {alphacheck[j]: Rc
                                 for j, c in r.associated_coroot() if (Rc := R(c))}

            # [e_r, e_s] and [e_r, f_s] with r != +/-s
            # We assume s is positive, as otherwise we negate
            #   both r and s and the resulting coefficient
            for j, s in enumerate(p_roots[i+1:], start=i+1):
                #j += i + 1  # Offset
                # Since h(s) >= h(r), we have s - r > 0 when s - r is a root
                # [f_r, e_s]
                if s - r in p_roots:
                    c = e_coeff(r, -s)
                    a, b = s-r, r
                    if self._p_roots_index[a] > self._p_roots_index[b]: # Note a != b
                        c *= -sp_sign[(b, a)]
                    else:
                        c *= sp_sign[(a, b)]
                    c = R(c)
                    if c:
                        s_coeffs[(-r, s)] = {a: -c}
                        s_coeffs[(r, -s)] = {-a: c}

                # [e_r, e_s]
                a = r + s
                if a in p_roots:
                    # (r, s) is a special pair
                    c = R(e_coeff(r, s) * sp_sign[(r, s)])
                    if c:
                        s_coeffs[(r, s)] = {a: c}
                        s_coeffs[(-r, -s)] = {-a: -c}

        return s_coeffs

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['A', 2])
            Lie algebra of ['A', 2] in the Chevalley basis
        """
        return "Lie algebra of {} in the Chevalley basis".format(self._cartan_type)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: latex(g)
            \mathfrak{g}(A_{2})_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return r"\mathfrak{{g}}({})_{{{}}}".format(latex(self._cartan_type), latex(self.base_ring()))

    def _test_structure_coeffs(self, **options):
        """
        Check the structure coefficients against the GAP implementation.

        EXAMPLES::

            sage: L = LieAlgebra(ZZ, cartan_type=['G',2])
            sage: L._test_structure_coeffs()
        """
        tester = self._tester(**options)
        ct = self.cartan_type()

        # Setup the GAP objects
        from sage.libs.gap.libgap import libgap
        L = libgap.SimpleLieAlgebra(ct.letter, ct.n, libgap(self.base_ring()))
        pos_B, neg_B, _ = libgap.ChevalleyBasis(L)
        gap_p_roots = libgap.PositiveRoots(libgap.RootSystem(L)).sage()
        # E, F, H = libgap.CanonicalGenerators(L)

        # Setup the conversion between the Sage roots and GAP roots.
        #   The GAP roots are given in terms of the weight lattice.
        p_roots = list(self._Q.positive_roots_by_height())
        WL = ct.root_system().weight_lattice()
        La = WL.fundamental_weights()
        convert = {WL(root): root for root in p_roots}
        index = {convert[sum(c*La[j+1] for j,c in enumerate(rt))]: i
                 for i, rt in enumerate(gap_p_roots)}

        # Run the check
        basis = self.basis()
        roots = frozenset(p_roots)
        for i,x in enumerate(p_roots):
            for y in p_roots[i+1:]:
                if x + y in roots:
                    c = basis[x].bracket(basis[y]).leading_coefficient()
                    a, b = (x + y).extraspecial_pair()
                    if (x, y) == (a, b): # If it already is an extra special pair
                        tester.assertEqual(pos_B[index[x]] * pos_B[index[y]],
                                           c * pos_B[index[x+y]],
                                           "extra special pair differ for [{}, {}]".format(x, y))
                    else:
                        tester.assertEqual(pos_B[index[x]] * pos_B[index[y]],
                                           c * pos_B[index[x+y]],
                                           "incorrect structure coefficient for [{}, {}]".format(x, y))
                if x - y in roots: # This must be a negative root if it is a root
                    c = basis[x].bracket(basis[-y]).leading_coefficient()
                    tester.assertEqual(pos_B[index[x]] * neg_B[index[y]],
                                       c * neg_B[index[x-y]],
                                       "incorrect structure coefficient for [{}, {}]".format(x, y))

    def _repr_generator(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: K = L.basis().keys()
            sage: L._repr_generator(K[0])
            'E[alpha[2]]'
            sage: L._repr_generator(K[4])
            'h2'
        """
        if m in self._Q.simple_coroots():
            return "h{}".format(m.support()[0])
        return IndexedGenerators._repr_generator(self, m)

    def _latex_generator(self, m):
        r"""
        Return a latex representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: K = L.basis().keys()
            sage: L._latex_generator(K[0])
            'E_{\\alpha_{2}}'
            sage: L._latex_generator(K[4])
            'h_{2}'
        """
        if m in self._Q.simple_coroots():
            return "h_{{{}}}".format(m.support()[0])
        return IndexedGenerators._latex_generator(self, m)

    def _basis_key(self, x):
        """
        Compare two basis element indices. We order the basis elements by
        positive roots, coroots, and negative roots and then according to
        height.

        OUTPUT:

        If ``x == y``, return 0. If ``x < y``, return -1. Else return 1.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: K = L.basis().keys()
            sage: S = sorted(K, key=L._basis_key); S
            [alpha[2],
             alpha[1],
             alpha[1] + alpha[2],
             alpha[1] + 2*alpha[2],
             alphacheck[1],
             alphacheck[2],
             -alpha[2],
             -alpha[1],
             -alpha[1] - alpha[2],
             -alpha[1] - 2*alpha[2]]
            sage: S == K
            False
        """
        if x in self._p_roots_index:
            return self._p_roots_index[x]
        if -x in self._p_roots_index:
            return (len(self._p_roots_index)
                    + self._cartan_type.rank()
                    + self._p_roots_index[-x])
        alphacheck = list(self._Q.simple_coroots())
        try:
            return len(self._p_roots_index) + alphacheck.index(x)
        except ValueError:
            raise KeyError(x)

    def degree_on_basis(self, m):
        """
        Return the degree of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: [L.degree_on_basis(m) for m in L.basis().keys()]
            [alpha[2], alpha[1], alpha[1] + alpha[2],
             2*alpha[1] + alpha[2], 3*alpha[1] + alpha[2],
             3*alpha[1] + 2*alpha[2],
             0, 0,
             -alpha[2], -alpha[1], -alpha[1] - alpha[2],
             -2*alpha[1] - alpha[2], -3*alpha[1] - alpha[2],
             -3*alpha[1] - 2*alpha[2]]
        """
        if m.parent() is self._Q:
            return m
        return self._Q.zero()

    def _negative_half_index_set(self):
        """
        Return the index set of the basis for the negative half of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 4)
            sage: L._negative_half_index_set()
            [-alpha[2], -alpha[1], -alpha[1] - alpha[2],
             -2*alpha[1] - alpha[2]]
        """
        return [-x for x in self._p_roots_index]

    def _weight_action(self, m, wt):
        """
        Return the action of the basis element indexed by ``m`` on ``wt``.

        INPUT:

        - ``m`` -- an index of a basis element of the Cartan subalgebra
        - ``wt`` -- a weight

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: mu = La[1] - 3/5*La[2]
            sage: ac = L.cartan_type().root_system().coroot_lattice().simple_roots()
            sage: L._weight_action(ac[1], mu)
            1
            sage: L._weight_action(ac[2], mu)
            -3/5
            sage: L._weight_action(ac[3], mu)
            0
        """
        ac = self._cartan_type.root_system().root_lattice().simple_coroots()
        aci = dict(ac.inverse_family())
        if m not in aci:
            raise ValueError("not an element in the Cartan subalgebra")
        R = self.base_ring()
        # This is a little ugly way to make sure we have the correct
        #   coroots. However, it does work as :meth:`scalar` is not smart
        #   enough in the ambient space to correctly convert things to do
        #   the scalar product.
        alc = wt.parent().simple_coroots()
        return R(wt.scalar(alc[aci[m]]))

    def affine(self, kac_moody=True):
        r"""
        Return the affine Lie algebra of ``self``.

        EXAMPLES::

            sage: sp6 = lie_algebras.sp(QQ, 6)
            sage: sp6
            Lie algebra of ['C', 3] in the Chevalley basis
            sage: sp6.affine()
            Affine Kac-Moody algebra of ['C', 3] in the Chevalley basis

            sage: L = LieAlgebra(QQ, cartan_type=['A',3], epsilon=[(1,2),(3,2)])
            sage: L.affine(False)
            Affine Lie algebra of ['A', 3] in the Chevalley basis
            sage: L.affine(True)
            Affine Kac-Moody algebra of ['A', 3] in the Chevalley basis
        """
        from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
        return AffineLieAlgebra(self, kac_moody=kac_moody)

    # Useful in creating the UEA
    @cached_method
    def indices_to_positive_roots_map(self):
        r"""
        Return the map from indices to positive roots.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.indices_to_positive_roots_map()
            {1: alpha[1], 2: alpha[2], 3: alpha[1] + alpha[2]}
        """
        return {i+1: r for i, r in enumerate(self._Q.positive_roots())}

    @cached_method
    def lie_algebra_generators(self, str_keys=False):
        r"""
        Return the Chevalley Lie algebra generators of ``self``.

        INPUT:

        - ``str_keys`` -- boolean (default: ``False``); set to ``True`` to have
          the indices indexed by strings instead of simple (co)roots

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: L.lie_algebra_generators()
            Finite family {alpha[1]: E[alpha[1]], -alpha[1]: E[-alpha[1]], alphacheck[1]: h1}
            sage: L.lie_algebra_generators(True)
            Finite family {'e1': E[alpha[1]], 'f1': E[-alpha[1]], 'h1': h1}
        """
        index_set = self._cartan_type.index_set()
        alpha = self._Q.simple_roots()
        alphacheck = self._Q.simple_coroots()
        B = self.basis()
        ret = {}

        if str_keys:
            for i in index_set:
                al = alpha[i]
                ret['e{}'.format(i)] = B[al]
                ret['f{}'.format(i)] = B[-al]
                ret['h{}'.format(i)] = B[alphacheck[i]]
            keys = (['e{}'.format(i) for i in index_set]
                    + ['f{}'.format(i) for i in index_set]
                    + ['h{}'.format(i) for i in index_set])
        else:
            for i in index_set:
                al = alpha[i]
                ret[al] = B[al]
                ret[-al] = B[-al]
                ret[alphacheck[i]] = B[alphacheck[i]]
            keys = ([alpha[i] for i in index_set]
                    + [-alpha[i] for i in index_set]
                    + [alphacheck[i] for i in index_set])

        return Family(keys, ret.__getitem__)

    @cached_method
    def _part_generators(self, positive=False):
        r"""
        Return the Lie algebra generators for the positive or
        negative half of ``self``.

        INPUT:

        - ``positive`` -- boolean (default: ``False``); if ``True``
          then return positive part generators, otherwise the return
          the negative part generators

        OUTPUT:

        A :func:`~sage.sets.family.Family` whose keys are the
        index set of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: list(L._part_generators(True))
            [E[alpha[1]],
             E[alpha[2]]]
        """
        I = self._cartan_type.index_set()
        al = self._Q.simple_roots()
        G = self.lie_algebra_generators()
        if positive:
            d = {i: G[al[i]] for i in I}
        else:
            d = {i: G[-al[i]] for i in I}
        return Family(I, d.__getitem__)

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self`` in the order of `e_i`, `f_i`,
        and `h_i`.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.gens()
            (E[alpha[1]], E[alpha[2]], E[-alpha[1]], E[-alpha[2]], h1, h2)
        """
        index_set = self._cartan_type.index_set()
        alpha = self._Q.simple_roots()
        alphacheck = self._Q.simple_coroots()
        B = self.basis()

        ret = [B[alpha[i]] for i in index_set]
        ret.extend(B[-alpha[i]] for i in index_set)
        ret.extend(B[alphacheck[i]] for i in index_set)
        return tuple(ret)

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.

        INPUT:

        - ``pos`` -- boolean (default: ``True``); if ``True``, then return
          `e_{\theta}`, otherwise return `f_{\theta}`

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.highest_root_basis_elt()
            E[alpha[1] + alpha[2]]
            sage: L.highest_root_basis_elt(False)
            E[-alpha[1] - alpha[2]]
        """
        theta = self._Q.highest_root()
        B = self.basis()
        if pos:
            return B[theta]
        return B[-theta]

    @cached_method
    def killing_form_matrix(self):
        r"""
        Return the matrix of the Killing form of ``self``.

        The rows and the columns of this matrix are indexed by the
        elements of the basis of ``self`` (in the order provided by
        :meth:`basis`).

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: g.killing_form_matrix()
            [ 0  0  0  0  0  6  0  0]
            [ 0  0  0  0  0  0  6  0]
            [ 0  0  0  0  0  0  0  6]
            [ 0  0  0 12 -6  0  0  0]
            [ 0  0  0 -6 12  0  0  0]
            [ 6  0  0  0  0  0  0  0]
            [ 0  6  0  0  0  0  0  0]
            [ 0  0  6  0  0  0  0  0]
        """
        B = self.basis()
        Q = self._Q
        from sage.matrix.constructor import matrix
        ret = matrix.zero(self.base_ring(), self._M.rank())
        keys = list(B.keys())
        for i, a in enumerate(keys):
            for j, b in enumerate(keys[i:], start=i):
                if ((a in Q) != (b in Q)) or (a in Q and a + b):
                    continue
                # either a and b are both coroots or a + b == 0
                temp = self.killing_matrix(B[a], B[b]).trace()
                ret[i, j] = temp
                ret[j, i] = temp
        ret.set_immutable()
        return ret

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``, where ``x``
        and ``y`` are two elements of ``self``.

        The Killing form is defined as

        .. MATH::

            \langle x \mid y \rangle
            = \operatorname{tr}\left( \operatorname{ad}_x
            \circ \operatorname{ad}_y \right).

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.killing_form(L.an_element(), L.an_element())
            36
            sage: B = L.basis()
            sage: matrix([[L.killing_form(a, b) for a in B] for b in B])
            [ 0  0  0  0  0  6  0  0]
            [ 0  0  0  0  0  0  6  0]
            [ 0  0  0  0  0  0  0  6]
            [ 0  0  0 12 -6  0  0  0]
            [ 0  0  0 -6 12  0  0  0]
            [ 6  0  0  0  0  0  0  0]
            [ 0  6  0  0  0  0  0  0]
            [ 0  0  6  0  0  0  0  0]
        """
        return x.value * self.killing_form_matrix() * y.value


class LieAlgebraChevalleyBasis_simply_laced(LieAlgebraChevalleyBasis):
    r"""
    A finite dimensional simply-laced Lie algebra in the Chevalley basis
    with structure coefficients given by an orientation of the Dynkin
    diagram.

    We follow Chapter 7.7 of [Ka1990]_, where the structure coefficients
    are given by an :meth:`asymmetry function <asymmetry_function>` defined by
    `\varepsilon(\alpha_i, \alpha_j) = -1` if there is an arrow `i \to j` in
    the Dynkin quiver (an orientation of the Dynkin diagram). However we twist
    `E_{\alpha}` by `\mathrm{sign}(\alpha)` so that `F_i = E_{-\alpha_i}`
    rather than its negative.

    EXAMPLES::

        sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(2, 1)])
        sage: L.e(1).bracket(L.e(2))
        E[alpha[1] + alpha[2]]

        sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(1, 2)])
        sage: L.e(1).bracket(L.e(2))
        -E[alpha[1] + alpha[2]]
    """
    def __init__(self, R, cartan_type, epsilon):
        """
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2], epsilon=[(2,1)])
            sage: TestSuite(L).run(elements=list(L.basis()))
            sage: L = LieAlgebra(QQ, cartan_type=['D', 4], epsilon=[(2,1), (3,2), (4,2)])
            sage: TestSuite(L).run(elements=list(L.basis()))  # long time
        """
        self._epsilon = epsilon
        super().__init__(R, cartan_type)

    def _construct_struct_coeffs(self, R, p_roots):
        r"""
        Construct the structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',1], epsilon=[])
            sage: sorted(L._construct_struct_coeffs(QQ, list(L._Q.positive_roots())).items(), key=str)
            [((alpha[1], -alpha[1]), {alphacheck[1]: 1}),
             ((alpha[1], alphacheck[1]), {alpha[1]: -2}),
             ((alphacheck[1], -alpha[1]), {-alpha[1]: -2})]
        """
        p_roots_set = set(p_roots)
        alphacheck = self._Q.simple_coroots()

        s_coeffs = {}
        for i, r in enumerate(p_roots):
            # [e_r, h_i] and [h_i, f_r]
            for ac in alphacheck:
                c = r.scalar(ac)
                if c == 0:
                    continue
                s_coeffs[(r, ac)] = {r: -c}
                s_coeffs[(ac, -r)] = {-r: -c}

            # [e_r, f_r]
            s_coeffs[(r, -r)] = {alphacheck[j]: c
                                 for j, c in r.associated_coroot()}

            # [e_r, e_s] and [e_r, f_s] with r != +/-s
            # We assume s is positive, as otherwise we negate
            #   both r and s and the resulting coefficient
            for j, s in enumerate(p_roots[i+1:], start=i+1):
                if r + s in p_roots_set:
                    coeff = R.prod((-1)**(ca*cb) if (ii, jj) in self._epsilon or ii == jj else 1
                                   for ii, ca in r._monomial_coefficients.items()
                                   for jj, cb in s._monomial_coefficients.items())
                    s_coeffs[r, s] = {r+s: coeff}
                    ht = sum(r.coefficients()) + sum(s.coefficients())
                    s_coeffs[-r, -s] = {-r-s: -coeff}
                if r - s in p_roots_set or s - r in p_roots_set:
                    coeff = R.prod((-1)**(ca*cb) if (ii, jj) in self._epsilon or ii == jj else 1
                                   for ii, ca in r._monomial_coefficients.items()
                                   for jj, cb in s._monomial_coefficients.items())
                    if r - s in p_roots_set:
                        s_coeffs[r, -s] = {r-s: -coeff}
                        s_coeffs[-r, s] = {s-r: coeff}
                    else:
                        s_coeffs[r, -s] = {r-s: coeff}
                        s_coeffs[-r, s] = {s-r: -coeff}

        return s_coeffs

    def asymmetry_function(self):
        r"""
        Return the asymmetry function of ``self``.

        An *asymmetry function* is a function `\varepsilon : Q \times Q
        \to \{1, -1\}` that satisfies the following properties:

        1. `\varepsilon(\alpha, \alpha) = (-1)^{(\alpha|\alpha)/2}`
        2. bimultiplicativity `\varepsilon(alpha + \alpha', \beta) =
           \varepsilon(\alpha, \beta) \varepsilon(\alpha', \beta)` and
           `\varepsilon(alpha, \beta + \beta') =
           \varepsilon(\alpha, \beta) \varepsilon(\alpha', \beta)`,

        where `(\alpha | \beta)` is the symmetric bilinear form on `Q` given
        by the Cartan matrix. Some consequences of these properties are that
        `\varepsilon(\alpha, 0) = \varepsilon(0, \beta) = 1` and
        `varepsilon(\alpha, \beta) \varepsilon(\beta, \alpha) =
        (-1)^{(\alpha|\beta)}`.

        OUTPUT:

        The asymmetry function as a ``dict`` consisting of pairs of all of
        the roots of `Q` and `0`.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2], epsilon=[(2,1)])
            sage: ep = L.asymmetry_function()
            sage: al = L.cartan_type().root_system().root_lattice().simple_roots()
            sage: ep[al[1], al[2]]
            1
            sage: ep[al[2],al[1]]
            -1

            sage: L = LieAlgebra(QQ, cartan_type=['A',2], epsilon=[(1,2)])
            sage: ep = L.asymmetry_function()
            sage: al = L.cartan_type().root_system().root_lattice().simple_roots()
            sage: ep[al[1], al[2]]
            -1
            sage: ep[al[2],al[1]]
            1
        """
        roots = set(self._Q.roots())
        al = self._Q.simple_roots()

        ep = {(r, r): (-1)**(r.scalar(r.associated_coroot()) // 2) for r in roots}
        next_level = set()
        for i in self._Q.index_set():
            # ep[i,0] = ep[i,-j] * ep[i,j] = ep[i,0]^2
            # => ep[i,0] = 1 = ep[i,-j] * ep[i,j]
            # => ep[-i,0] = 1 = ep[-i,-j] * ep[-i,j]
            ep[al[i], self._Q.zero()] = 1
            ep[-al[i], self._Q.zero()] = 1
            for j in self._Q.index_set():
                if i == j or (i, j) in self._epsilon:
                    ep[al[i], al[j]] = -1
                else:
                    ep[al[i], al[j]] = 1
                next_level.add((al[i], al[j]))
                ep[al[i], -al[j]] = ep[al[i], al[j]]
                next_level.add((al[i], -al[j]))
                ep[-al[i], al[j]] = ep[al[i], al[j]]
                next_level.add((-al[i], al[j]))
                ep[-al[i], -al[j]] = ep[al[i], al[j]]
                next_level.add((-al[i], -al[j]))

        while next_level:
            cur = next_level
            next_level = set()
            for r, s in cur:
                for a in al:
                    prev = r
                    temp = r + a
                    while temp in roots or temp == 0:
                        if (temp, s) in ep:
                            break
                        next_level.add((temp, s))
                        ep[temp, s] = ep[prev, s] * ep[a, s]
                        prev = temp
                        temp += a

                    prev = r
                    temp = r - a
                    while temp in roots or temp == 0:
                        if (temp, s) in ep:
                            break
                        next_level.add((temp, s))
                        ep[temp, s] = ep[prev, s] * ep[-a, s]
                        prev = temp
                        temp -= a

                    prev = s
                    temp = s + a
                    while temp in roots or temp == 0:
                        if (r, temp) in ep:
                            break
                        next_level.add((r, temp))
                        ep[r, temp] = ep[r, prev] * ep[r, a]
                        prev = temp
                        temp += a

                    prev = s
                    temp = s - a
                    while temp in roots or temp == 0:
                        if (r, temp) in ep:
                            break
                        next_level.add((r, temp))
                        ep[r, temp] = ep[r, prev] * ep[r, -a]
                        prev = temp
                        temp -= a

        return ep

    def _test_structure_coeffs(self, **options):
        r"""
        Check the structure coefficients using the :meth:`asymmetry_function`.

        EXAMPLES::

            sage: L = LieAlgebra(ZZ, cartan_type=['A',4], epsilon=[(1,2), (3,2), (3,4)])
            sage: L._test_structure_coeffs()
        """
        tester = self._tester(**options)

        ep = self.asymmetry_function()

        roots = set(self._Q.roots())
        p_roots = set(self._Q.positive_roots())
        B = self.basis()
        for r in roots:
            for s in roots:
                if r + s not in roots:
                    continue
                x = B[r].bracket(B[s])
                tester.assertEqual(list(x.support()), [r+s], f"[{r}, {s}] = {x} is not a root vector")
                sign = 1 if (r in p_roots) == (s in p_roots) else -1
                if (r + s) not in p_roots:
                    sign = -sign
                tester.assertEqual(x[r+s], sign * ep[r, s], f"[{r}, {s}] = {x[r+s]} != {sign*ep[r,s]}")
