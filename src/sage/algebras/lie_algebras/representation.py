r"""
Representations of Lie algebras

AUTHORS:

- Travis Scrimshaw (2023-08-31): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family, AbstractFamily
from sage.matrix.constructor import matrix
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.modules import Modules


class Representation_abstract:
    r"""
    Mixin class for (left) representations of Lie algebras.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra
    """
    def __init__(self, lie_algebra):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: R = L.trivial_representation()
            sage: TestSuite(R).run()
        """
        self._lie_algebra = lie_algebra

    def lie_algebra(self):
        r"""
        Return the Lie algebra whose representation ``self`` is.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 4)
            sage: R = L.trivial_representation()
            sage: R.lie_algebra() is L
            True
        """
        return self._lie_algebra

    def side(self):
        r"""
        Return that ``self`` is a left representation.

        OUTPUT: the string ``'left'``

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 4)
            sage: R = L.trivial_representation()
            sage: R.side()
            'left'
        """
        return 'left'

    def _test_representation(self, **options):
        r"""
        Check (on some elements) that ``self`` is a representation of the
        given Lie algebra using the basis of the Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: f = {b: b.adjoint_matrix() for b in L.basis()}
            sage: R = L.representation(f)
            sage: R._test_representation()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        elts = self._lie_algebra.basis()
        if elts.cardinality() == float('inf'):
            elts = list(elts.some_elements())
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(elts, 2, tester._max_runs):
            for v in S:
                tester.assertEqual(x.bracket(y) * v, x * (y * v) - y * (x * v))

    def representation_matrix(self, elt):
        """
        Return the matrix for the action of ``elt`` on ``self``.

        EXAMPLES::

            sage: H1 = lie_algebras.Heisenberg(QQ, 1)
            sage: F = H1.faithful_representation(algorithm='minimal')
            sage: P1 = F.representation_matrix(H1.gen(0)); P1
            [0 0 0]
            [0 0 0]
            [1 0 0]
            sage: Q1 = F.representation_matrix(H1.gen(1)); Q1
            [ 0  0  0]
            [ 0  0 -1]
            [ 0  0  0]
            sage: Z = P1.commutator(Q1); Z
            [0 0 0]
            [1 0 0]
            [0 0 0]
            sage: P1.commutator(Z) == Q1.commutator(Z) == 0
            True
            sage: (H1.gen(0) * F.an_element()).to_vector()
            (0, 0, 2)
            sage: P1 * F.an_element().to_vector()
            (0, 0, 2)
            sage: (H1.gen(1) * F.an_element()).to_vector()
            (0, -3, 0)
            sage: Q1 * F.an_element().to_vector()
            (0, -3, 0)
            sage: (H1.basis()['z'] * F.an_element()).to_vector()
            (0, 2, 0)
            sage: Z * F.an_element().to_vector()
            (0, 2, 0)
        """
        B = self.basis()
        return matrix([(elt * B[k]).to_vector() for k in self.get_order()]).transpose()


class RepresentationByMorphism(CombinatorialFreeModule, Representation_abstract):
    r"""
    Representation of a Lie algebra defined by a Lie algebra morphism.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra
    - ``f`` -- the Lie algebra morphism defining the action of the basis
      elements of ``lie_algebra``
    - ``index_set`` -- (optional) the index set of the module basis
    - ``on_basis`` -- boolean (default: ``False``); the function `f` defines a
      map from the basis elements or from a generic element of ``lie_algebra``

    If `f` is encoded as a ``dict`` or ``Family``, then the keys must
    be indices of the basis of ``lie_algebra`` and the values being the
    corresponding matrix defining the action. This sets ``on_basis=True``.

    EXAMPLES::

        sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
        sage: f = {x: Matrix([[1,0],[0,0]]), y: Matrix([[0,1],[0,0]])}
        sage: L.representation(f)
        Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
               [1 0]
        x |--> [0 0]
               [0 1]
        y |--> [0 0]

    We construct the direct sum of two copies of the trivial representation
    for an infinite dimensional Lie algebra::

        sage: L = lie_algebras.Affine(QQ, ['E',6,1])
        sage: R = L.representation(lambda b: matrix.zero(QQ, 2), index_set=['a','b'])
        sage: x = L.an_element()
        sage: v = R.an_element(); v
        2*R['a'] + 2*R['b']
        sage: x * v
        0

    We construct a finite dimensional representation of the affline Lie algebra
    of type `A_2^{(1)}`::

        sage: L = lie_algebras.Affine(QQ, ['A',2,1]).derived_subalgebra()
        sage: Phi_plus = list(RootSystem(['A',2]).root_lattice().positive_roots())
        sage: def aff_action(key):
        ....:     mat = matrix.zero(QQ, 3)
        ....:     if key == 'c':  # central element
        ....:         return mat
        ....:     b, ell = key
        ....:     if b in Phi_plus:  # positive root
        ....:         ind = tuple(sorted(b.to_ambient().support()))
        ....:         mat[ind] = 1
        ....:         if ind[0] + 1 != ind[1]:
        ....:             mat[ind] = -1
        ....:     elif -b in Phi_plus:  # negative root
        ....:         ind = tuple(sorted(b.to_ambient().support(), reverse=True))
        ....:         mat[ind] = 1
        ....:         if ind[0] - 1 != ind[1]:
        ....:             mat[ind] = -1
        ....:     else:  # must be in the Cartan
        ....:         i = b.leading_support()
        ....:         mat[i,i] = -1
        ....:         mat[i-1,i-1] = 1
        ....:     return mat
        sage: F = Family(L.basis(), aff_action, name="lifted natural repr")
        sage: R = L.representation(index_set=range(1,4), on_basis=F)
        sage: x = L.an_element(); x
        (E[alpha[2]] + E[alpha[1]] + h1 + h2 + E[-alpha[2]] + E[-alpha[1]])#t^0
         + (E[-alpha[1] - alpha[2]])#t^1 + (E[alpha[1] + alpha[2]])#t^-1 + c
        sage: v = R.an_element(); v
        2*R[1] + 2*R[2] + 3*R[3]
        sage: x * v
        R[1] + 5*R[2] - 3*R[3]
        sage: R._test_representation()  # verify that it is a representation
    """
    @staticmethod
    def __classcall_private__(cls, lie_algebra, f=None, index_set=None, on_basis=False, **kwargs):
        r"""
        Normalize inpute to ensure a unique representation.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f1 = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R1 = L.representation(f1)
            sage: f2 = Family({x: Matrix([[1,0],[0,0]]), y: Matrix(QQ, [[0,1],[0,0]])})
            sage: R2 = L.representation(f2)
            sage: R1 is R2
            True

        TESTS::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0]]), 'y': Matrix([[0,1]])}
            sage: L.representation(f)
            Traceback (most recent call last):
            ...
            ValueError: all matrices must be square

            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0]])}
            sage: L.representation(f)
            Traceback (most recent call last):
            ...
            ValueError: all matrices must be square of size 2

            sage: L.representation(index_set=[1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: either 'f' or 'on_basis' must be specified
            sage: L.representation(on_basis=lambda x: QQ.zero())
            Traceback (most recent call last):
            ...
            ValueError: the index set needs to be specified
        """
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        base = lie_algebra.base_ring()
        C = Modules(base).WithBasis().FiniteDimensional()
        C = C.or_subcategory(kwargs.pop('category', C))
        B = lie_algebra.basis()
        if not isinstance(on_basis, bool):
            f = on_basis
            on_basis = True
        if isinstance(f, AbstractFamily):
            if f.cardinality() < float('inf'):
                f = dict(f)
            on_basis = True
        if isinstance(f, dict):
            data = {}
            dim = None
            for k, mat in f.items():
                if k in B:
                    k = k.leading_support()
                if not mat.is_square():
                    raise ValueError("all matrices must be square")
                if dim is None:
                    dim = mat.nrows()
                elif mat.nrows() != dim or mat.ncols() != dim:
                    raise ValueError("all matrices must be square of size {}".format(dim))
                data[k] = mat.change_ring(base)
                data[k].set_immutable()

            if index_set is None:
                index_set = FiniteEnumeratedSet(range(dim))
            f = Family(data)
            on_basis = True

        if f is None:
            raise ValueError("either 'f' or 'on_basis' must be specified")
        if index_set is None:
            raise ValueError("the index set needs to be specified")

        index_set = FiniteEnumeratedSet(index_set)

        return super(cls, RepresentationByMorphism).__classcall__(cls, lie_algebra,
             f, index_set, on_basis, category=C, **kwargs)

    def __init__(self, lie_algebra, f, index_set, on_basis, category, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: R = L.representation(f)
            sage: TestSuite(R).run()
        """
        if on_basis:
            self._family = f
            self._f = f.__getitem__
        else:
            self._f = f
        prefix = kwargs.pop("prefix", 'R')
        self._on_basis = on_basis

        Representation_abstract.__init__(self, lie_algebra)
        CombinatorialFreeModule.__init__(self, lie_algebra.base_ring(), index_set,
                                         category=category, prefix=prefix, **kwargs)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
            sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
            sage: L.representation(f)
            Representation of Lie algebra on 2 generators (x, y) over Rational Field defined by:
                   [1 0]
            x |--> [0 0]
                   [0 1]
            y |--> [0 0]

            sage: L = lie_algebras.Affine(QQ, ['E',6,1])
            sage: F = Family(L.basis(), lambda b: matrix.zero(QQ, 2), name="zero map")
            sage: L.representation(F, index_set=['a','b'], on_basis=True)
            Representation of Affine Kac-Moody algebra of ['E', 6] in the Chevalley basis defined by:
            Lazy family (zero map(i))_{i in Lazy family...}

            sage: L.representation(lambda b: matrix.zero(QQ, 2), index_set=['a','b'])
            Representation of Affine Kac-Moody algebra of ['E', 6] in the Chevalley basis defined by:
            <function <lambda> at 0x...>
        """
        ret = "Representation of {} defined by:".format(self._lie_algebra)
        from sage.typeset.ascii_art import ascii_art
        if self._on_basis:
            B = self._lie_algebra.basis()
            if B.cardinality() < float('inf'):
                for k in B.keys():
                    ret += '\n' + repr(ascii_art(B[k], self._f(k), sep=" |--> ", sep_baseline=0))
            else:
                ret += '\n' + repr(self._family)
        else:
            ret += '\n' + repr(self._f)
        return ret

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'y':1}})
                sage: f = {'x': Matrix([[1,0],[0,0]]), 'y': Matrix([[0,1],[0,0]])}
                sage: R = L.representation(f)
                sage: v = R.an_element(); v
                2*R[0] + 2*R[1]
                sage: x * v
                2*R[0]
                sage: y * v
                2*R[0]
                sage: (2*x + 5*y) * v
                14*R[0]
                sage: v * x
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...

                sage: v = sum((i+4) * b for i, b in enumerate(R.basis())); v
                4*R[0] + 5*R[1]
                sage: (1/3*x - 5*y) * v
                -71/3*R[0]

                sage: L = lie_algebras.Affine(QQ, ['E',6,1])
                sage: F = Family(L.basis(), lambda b: matrix.zero(QQ, 2), name="zero map")
                sage: R = L.representation(F, index_set=['a','b'], on_basis=True)
                sage: R.an_element()
                2*R['a'] + 2*R['b']
                sage: L.an_element() * R.an_element()
                0
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                if not self:  # we are (already) the zero vector
                    return self
                scalar = P._lie_algebra(scalar)
                if not scalar:  # we are acting by zero
                    return P.zero()
                if P._on_basis:
                    mat = sum(c * P._f(k) for k, c in scalar.monomial_coefficients(copy=False).items())
                else:
                    mat = P._f(scalar)
                return P.from_vector(mat * self.to_vector())

            return super()._acted_upon_(scalar, self_on_left)


class TrivialRepresentation(CombinatorialFreeModule, Representation_abstract):
    r"""
    The trivial representation of a Lie algebra.

    The trivial representation of a Lie algebra `L` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `L`
    acts by zero.

    INPUT:

    - ``lie_algebra`` -- a Lie algebra

    REFERENCES:

    - :wikipedia:`Trivial_representation`
    """
    def __init__(self, lie_algebra, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: R = L.trivial_representation()
            sage: TestSuite(R).run()
        """
        R = lie_algebra.base_ring()
        cat = Modules(R).WithBasis().FiniteDimensional()
        Representation_abstract.__init__(self, lie_algebra)
        CombinatorialFreeModule.__init__(self, R, ['v'], prefix='T', category=cat, **kwargs)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.VirasoroAlgebra(QQ)
            sage: L.trivial_representation()
            Trivial representation of The Virasoro algebra over Rational Field
        """
        return "Trivial representation of {}".format(self._lie_algebra)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: L = lie_algebras.VirasoroAlgebra(QQ)
                sage: R = L.trivial_representation()
                sage: L.an_element() * R.an_element()
                0
                sage: R.an_element() * L.an_element()
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...
                sage: 3 / 5 * R.an_element()
                6/5*T['v']
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                return P.zero()
            return super()._acted_upon_(scalar, self_on_left)


class FaithfulRepresentationNilpotentPBW(CombinatorialFreeModule, Representation_abstract):
    r"""
    Return a faithful reprensetation of a nilpotent Lie algebra
    constructed using the PBW basis.

    Let `L` be a `k`-step nilpotent Lie algebra. Define a weight function
    on elements in `L` by the lower central series of `L`. Then a faithful
    representation of `L` is `U(L) / U(L)^{k+1}`, where `U(L)^{k+1}`
    is the (twosided) ideal of `U(L)` generated by all monomials
    of weight at least `k + 1`.

    We can also expand the ideal keeping the property that `I \cap Z(L) = 0`.
    The resulting quotient `U(L) / I` remains faithful and is a minimal
    faithful representation of `L` in the sense that it has no faithful
    submodules or quotients. (Note: this is not necessarily the smallest
    dimensional faithful representation of `L`.)

    We consider an example of the rank 2 Heisenberg Lie algebra,
    but with a non-standard basis given by `a = p_1 + z`, `b = q_1`,
    and `c = q_1 + z`::

        sage: scoeffs = {('a','b'): {'b':-1, 'c':1}, ('a','c'): {'b':-1, 'c':1}}
        sage: L.<a,b,c> = LieAlgebra(QQ, scoeffs)
        sage: TestSuite(L).run(elements=list(L.basis()))
        sage: L.is_nilpotent()
        True
        sage: L.derived_series()
        (Lie algebra on 3 generators (a, b, c) over Rational Field,
         Ideal (b - c) of Lie algebra on 3 generators (a, b, c) over Rational Field,
         Ideal () of Lie algebra on 3 generators (a, b, c) over Rational Field)
        sage: F = L.faithful_representation()
        sage: L.an_element() * F.an_element()
        2*F[1, 0, 0] + 8*F[1, 1, 0] + 3*F[2, 0, 0] + 4*F[0, 1, 0]
         + 4*F[0, 2, 0] + 4*F[0, 0, 1]

        sage: MF = L.faithful_representation(algorithm='minimal')
        sage: MF.dimension()
        3
        sage: [MF.representation_matrix(be) for be in L.basis()]
        [
        [0 0 0]  [ 0  0  0]  [ 0  0  0]
        [0 0 0]  [ 0  0 -1]  [ 1  0 -1]
        [1 0 0], [ 0  0  0], [ 0  0  0]
        ]

    An example with ``minimal=True`` for `H_2 \oplus A_1`, where `A_1` is
    a `1`-dimensional Abelian Lie algebra::

        sage: scoeffs = {('a','b'): {'b':-1, 'c':1}, ('a','c'): {'b':-1, 'c':1}}
        sage: L.<a,b,c,d> = LieAlgebra(QQ, scoeffs)
        sage: F = L.faithful_representation(); F
        Faithful 11 dimensional representation of Lie algebra on 4
         generators (a, b, c, d) over Rational Field
        sage: MF = L.faithful_representation(algorithm='minimal'); MF
        Minimal faithful representation of Lie algebra on 4
         generators (a, b, c, d) over Rational Field
        sage: MF.dimension()
        4

    INPUT:

    - ``minimal`` -- boolean (default: ``False``); whether to construct
      the minimal basis or not

    REFERENCES:

    - [BEdG2009]_
    """
    def __init__(self, L, minimal=False):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: H2 = lie_algebras.Heisenberg(QQ, 2)
            sage: F = H2.faithful_representation()
            sage: TestSuite(F).run(elements=list(F.basis()))
            sage: MF = H2.faithful_representation(algorithm='minimal')
            sage: TestSuite(MF).run(elements=list(MF.basis()))

            sage: sc = {('a','b'): {'b':-1, 'c':1}, ('a','c'): {'b':-1, 'c':1}}
            sage: L.<a,b,c> = LieAlgebra(QQ, sc)
            sage: F = L.faithful_representation()
            sage: TestSuite(F).run(elements=list(F.basis()))
            sage: MF = L.faithful_representation(algorithm='minimal')
            sage: TestSuite(MF).run(elements=list(MF.basis()))
        """
        LCS = L.lower_central_series()
        if LCS[-1].dimension() != 0:
            raise ValueError("the Lie algebra must be nilpotent")
        # construct an appropriate basis of L
        basis_by_deg = {}
        self._step = len(LCS) - 1
        self._minimal = minimal
        if self._minimal:
            Z = L.center()
            ZB = [L(b) for b in Z.basis()]
            prev = LCS[-1]
            for D in reversed(LCS[:-1]):
                cur = []
                for ind in range(len(ZB) - 1, -1, -1):
                    z = ZB[ind]
                    if z in D:
                        ZB.pop(ind)
                        cur.append(z)
                k = self._step - len(basis_by_deg)
                basis_by_deg[k] = cur
                temp = [bred for b in D.basis() if (bred := Z.reduce(prev.reduce(L(b))))]
                basis_by_deg[k].extend(L.echelon_form(temp))
                prev = D
        else:
            prev = LCS[-1]
            for D in reversed(LCS[:-1]):
                temp = [L(bred) for b in D.basis() if (bred := prev.reduce(L(b)))]
                basis_by_deg[self._step - len(basis_by_deg)] = L.echelon_form(temp)
                prev = D

        L_basis = sum((basis_by_deg[deg] for deg in sorted(basis_by_deg)), [])

        if all(len(b.support()) == 1 for b in L_basis):
            self._Lp = L
        else:
            cob = matrix([b._vector_() for b in L_basis]).transpose()
            self._invcob = cob.inverse()
            scoeffs = {}
            for i, b in enumerate(L_basis):
                for j, bp in enumerate(L_basis[i+1:], start=i + 1):
                    scoeffs[i, j] = (self._invcob * b.bracket(bp)._vector_()).dict()
            index_set = tuple(range(L.dimension()))
            from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
            self._Lp = LieAlgebra(L.base_ring(), scoeffs, index_set=index_set)

        self._pbw = self._Lp.pbw_basis()
        self._degrees = tuple(sum(([deg] * len(B) for deg, B in sorted(basis_by_deg.items())), []))

        from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
        from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
        indices = DisjointUnionEnumeratedSets([WeightedIntegerVectors(n, self._degrees)
                                               for n in range(self._step+1)])

        if self._minimal:
            X = {tuple(index) for index in indices}
            monoid = self._pbw._indices
            I = monoid._indices
            one = L.base_ring().one()
            pbw_gens = self._pbw.algebra_generators()
            ZB = frozenset([L(b) for b in Z.basis()])
            Zind = [i for i, b in enumerate(L_basis) if b in ZB]
            Ztup = set()
            for i in Zind:
                vec = [0] * L.dimension()
                vec[i] = 1
                Ztup.add(tuple(vec))

            def as_exp(s):
                sm = s._monomial
                return tuple([sm[i] if i in sm else 0 for i in I])

            def test_ideal(m, X):
                elt = self._pbw.element_class(self._pbw, {monoid(list(zip(I, m))): one})
                for g in pbw_gens:
                    gelt = g * elt
                    if any(as_exp(s) in X for s in gelt.support()):
                        return False
                return True

            to_remove = {None}
            while to_remove:
                X -= to_remove
                to_remove = set()
                for m in X:
                    m = tuple(m)
                    if m in Ztup or not test_ideal(m, X):
                        continue
                    to_remove.add(m)
            indices = sorted(X)

        Representation_abstract.__init__(self, L)
        CombinatorialFreeModule.__init__(self, L.base_ring(), indices, prefix='F', bracket=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: H2 = lie_algebras.Heisenberg(QQ, 2)
            sage: H2.faithful_representation()
            Faithful 16 dimensional representation of Heisenberg algebra
             of rank 2 over Rational Field
        """
        if self._minimal:
            return "Minimal faithful representation of {}".format(self._lie_algebra)
        return "Faithful {} dimensional representation of {}".format(self.dimension(), self._lie_algebra)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: H2 = lie_algebras.Heisenberg(QQ, 2)
            sage: latex(H2.faithful_representation())
            U(\text{\texttt{Heisenberg...}}) / U(\text{\texttt{Heisenberg...}})^{3}
        """
        from sage.misc.latex import latex
        g = latex(self._lie_algebra)
        ret = "U({0}) / U({0})^{{{1}}}".format(g, self._step + 1)
        if self._minimal:
            return "\\min " + ret
        return ret

    def _project(self, elt):
        r"""
        The projection to ``self`` from the PBW basis.

        EXAMPLES::

            sage: sc = {('a','b'): {'b':-1, 'c':1}, ('a','c'): {'b':-1, 'c':1}}
            sage: L.<a,b,c> = LieAlgebra(QQ, sc)
            sage: F = L.faithful_representation()
            sage: elt = F._to_pbw(a + b + c)^2; elt
            PBW[0]^2 + 4*PBW[0]*PBW[1] - 2*PBW[0]*PBW[2] + 4*PBW[1]^2
             - 4*PBW[1]*PBW[2] + PBW[2]^2 + 2*PBW[2]
            sage: F._project(elt)
            2*F[0, 0, 1] + 4*F[0, 2, 0] + 4*F[1, 1, 0] + F[2, 0, 0]
            sage: F._project(F._to_pbw(a + b + c)^3)
            0
        """
        ret = {}
        I = self._pbw._indices._indices
        if self._minimal:
            for m, c in elt._monomial_coefficients.items():
                mm = m._monomial
                vec = tuple([mm[i] if i in mm else 0 for i in I])
                if vec in self._indices:
                    ret[self._indices(vec)] = c
        else:
            for m, c in elt._monomial_coefficients.items():
                mm = m._monomial
                vec = [mm[i] if i in mm else 0 for i in I]
                if sum(e * d for e, d in zip(vec, self._degrees)) <= self._step:
                    ret[self._indices(vec)] = c
        return self.element_class(self, ret)

    def _to_pbw(self, elt):
        """
        Return the PBW element corresponding to ``elt``.

        EXAMPLES::

            sage: sc = {('a','b'): {'b':-1, 'c':1}, ('a','c'): {'b':-1, 'c':1}}
            sage: L.<a,b,c> = LieAlgebra(QQ, sc)
            sage: F = L.faithful_representation()
            sage: F._to_pbw(a)
            PBW[0]
            sage: F._to_pbw(b)
            PBW[1]
            sage: F._to_pbw(c)
            PBW[1] - PBW[2]

            sage: H2 = lie_algebras.Heisenberg(QQ, 2)
            sage: F = H2.faithful_representation()
            sage: F._to_pbw(sum(H2.basis()))
            PBW['p1'] + PBW['p2'] + PBW['q1'] + PBW['q2'] + PBW['z']
        """
        if self._Lp is self._lie_algebra:
            return self._pbw(elt)
        return self._pbw(self._Lp.from_vector(self._invcob * elt._vector_()))

    class Element(CombinatorialFreeModule.Element):
        def _lift_pbw(self):
            """
            Return ``self`` as an element of the PBW basis.

            EXAMPLES::

                sage: H2 = lie_algebras.Heisenberg(QQ, 2)
                sage: F = H2.faithful_representation()
                sage: F.an_element()._lift_pbw()
                3*PBW['q1'] + 2*PBW['q2'] + 2
            """
            P = self.parent()
            monoid = P._pbw._indices
            I = monoid._indices
            return P._pbw.element_class(P._pbw, {monoid(list(zip(I, m))): coeff
                                                 for m, coeff in self._monomial_coefficients.items()})

        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: H2 = lie_algebras.Heisenberg(QQ, 2)
                sage: F = H2.faithful_representation()
                sage: H2.an_element()
                p1
                sage: F.an_element()
                2*F[0, 0, 0, 0, 0] + 2*F[0, 0, 0, 1, 0] + 3*F[0, 0, 1, 0, 0]
                sage: H2.an_element() * F.an_element()
                2*F[1, 0, 0, 0, 0] + 2*F[1, 0, 0, 1, 0] + 3*F[1, 0, 1, 0, 0]
                sage: 5 * F.an_element()
                10*F[0, 0, 0, 0, 0] + 10*F[0, 0, 0, 1, 0] + 15*F[0, 0, 1, 0, 0]
            """
            P = self.parent()
            if scalar in P._lie_algebra:
                if self_on_left:
                    return None
                if not self:  # we are (already) the zero vector
                    return self
                scalar = P._lie_algebra(scalar)
                return P._project(P._to_pbw(scalar) * self._lift_pbw())

            return super()._acted_upon_(scalar, self_on_left)


class FaithfulRepresentationPBWPosChar(CombinatorialFreeModule, Representation_abstract):
    r"""
    A faithful representation of a finite dimensional Lie algebra
    in positive characteristic.

    .. WARNING::

        This is often a very large dimensional representation relative
        to the dimension of the Lie algebra.

    ALGORITHM:

    We implement the algorithm given in [deG2000] Section 6.6. Let `L`
    be a finite dimensional Lie algebra over a ring of characteristic `p`
    with basis `(b_1, \ldots, b_n)`. We compute (monic) `p`-polynomials
    `f_i` such that `A = \mathrm{ad}(b_i)` (the adjoint action of `b_i`)
    solves `f_i(A) = 0` by using minimal polynomial of `A`. The
    `(f_1, \ldots, f_n)` is a Gröbner basis for an ideal `I` of the
    universal enveloping algebra `U(L)` such that the quotient `U(L) / I`
    is a faithful representation of `L`.

    EXAMPLES::

        sage: sl2 = LieAlgebra(GF(3), cartan_type=['A',1])
        sage: F = sl2.faithful_representation()
        sage: F
        Faithful representation with p-multiplicities (1, 3, 1) of Lie algebra
         of ['A', 1] in the Chevalley basis
        sage: F.dimension()
        243
    """
    def __init__(self, L):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: sl2 = LieAlgebra(GF(3), cartan_type=['A',1])
            sage: F = sl2.faithful_representation()
            sage: TestSuite(F).run()
        """
        R = L.base_ring()
        self._p = R.characteristic()
        if self._p == 0:
            raise ValueError("the Lie algebra must be over a ring of positive characteristic")

        self._pbw = L.pbw_basis()
        self._key_order = tuple(self._pbw.algebra_generators().keys())

        # calculate the Gröbner basis and p-exponents
        gb = []
        p_exp = []
        B = L.basis()
        for k in self._key_order:
            b = B[k]
            ad = b.adjoint_matrix()
            g = ad.minpoly()
            d = g.degree()
            # TODO: Use the sparse polynomial ring?
            x = g.parent().gen()
            r = [x**(self._p**i) % g for i in range(d+1)]
            deg = max(ri.degree() for ri in r)
            mat = matrix(R, [[ri[j] for ri in r] for j in range(deg+1)])
            la = mat.right_kernel_matrix()[0]
            if la:
                mongen = self._pbw._indices.monoid_generators()[k]
                gb.append(self._pbw._from_dict({mongen ** (self._p ** i): val
                                                for i, val in enumerate(la) if val},
                                               remove_zeros=False))
                p_exp.append(max(la.support()))

        self._groebner_basis = gb
        self._p_exp = tuple(p_exp)
        self._degrees = [self._p ** m for m in self._p_exp]

        from sage.groups.abelian_gps.abelian_group import AbelianGroup
        indices = AbelianGroup(self._degrees)

        Representation_abstract.__init__(self, L)
        CombinatorialFreeModule.__init__(self, R, indices, prefix='', bracket=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: sl3 = LieAlgebra(GF(3), cartan_type=['A',2])
            sage: sl3.faithful_representation()
            Faithful representation with p-multiplicities (1, 1, 1, 3, 3, 1, 1, 1)
             of Lie algebra of ['A', 2] in the Chevalley basis
        """
        return "Faithful representation with p-multiplicities {} of {}".format(self.p_exponents(), self._lie_algebra)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: sl2 = LieAlgebra(GF(3), cartan_type=['A',1])
            sage: latex(sl2.faithful_representation())
            U(\mathfrak{g}(A_{1})_{\Bold{F}_{3}}) / \langle PBW_{\alpha_{1}}^{3},
             2 PBW_{\alpha^\vee_{1}}^{27} + PBW_{\alpha^\vee_{1}},
             PBW_{-\alpha_{1}}^{3} \rangle
        """
        from sage.misc.latex import latex
        g = latex(self._lie_algebra)
        data = ', '.join(latex(f) for f in self._groebner_basis)
        return "U({}) / \\langle {} \\rangle".format(g, data)

    @cached_method
    def p_exponents(self):
        """
        Return the `p`-exponents of ``self``.

        Let `p` be the characteristic of the base ring of ``self``.
        The `p`-*exponents* are the exponents `m_i` such that the `i`-th
        `p`-polynomial `f_i` is of degree `p^{m_i}`.

        EXAMPLES::

            sage: sp4 = LieAlgebra(GF(3), cartan_type=['C',2])
            sage: F = sp4.faithful_representation()
            sage: F.p_exponents()
            (1, 1, 1, 1, 3, 3, 1, 1, 1, 1)
        """
        return self._p_exp

    def groebner_basis(self):
        """
        Return the defining Gröbner basis of ``self``.

        EXAMPLES::

            sage: sp4 = LieAlgebra(GF(3), cartan_type=['C',2])
            sage: F = sp4.faithful_representation()
            sage: F.groebner_basis()
            [PBW[alpha[2]]^3,
             PBW[alpha[1]]^3,
             PBW[alpha[1] + alpha[2]]^3,
             PBW[2*alpha[1] + alpha[2]]^3,
             2*PBW[alphacheck[1]]^27 + PBW[alphacheck[1]],
             2*PBW[alphacheck[2]]^27 + PBW[alphacheck[2]],
             PBW[-alpha[2]]^3,
             PBW[-alpha[1]]^3,
             PBW[-alpha[1] - alpha[2]]^3,
             PBW[-2*alpha[1] - alpha[2]]^3]
        """
        return self._groebner_basis

    def _project(self, x):
        r"""
        The projection to ``self`` from the PBW basis.

        EXAMPLES::

            sage: sl2 = LieAlgebra(GF(3), cartan_type=['A',1])
            sage: F = sl2.faithful_representation()
            sage: PBW = F._pbw
            sage: elt = PBW.an_element(); elt
            PBW[alpha[1]]^2*PBW[alphacheck[1]]^2*PBW[-alpha[1]]^3
             + 2*PBW[alpha[1]] + 1
            sage: F._project(elt)
            1 + 2*f0
            sage: F._project(elt^2)
            1 + f0 + f0^2
            sage: F._project(elt^3)
            1

            sage: elt = PBW(sum(sl2.basis())); elt
            PBW[alpha[1]] + PBW[alphacheck[1]] + PBW[-alpha[1]]
            sage: F._project(elt)
            f2 + f1 + f0
            sage: F._project(elt^2)
            2*f2 + f2^2 + 2*f1 + 2*f1*f2 + f1^2 + 2*f0 + 2*f0*f2 + 2*f0*f1 + f0^2
            sage: F._project(elt^3)
            2*f2 + f1 + f1^3 + 2*f0
            sage: F._project(elt^4)
            f2 + 2*f2^2 + f1 + f1^2 + f1^3*f2 + f1^4 + f0 + f0*f2 + f0*f1^3 + 2*f0^2
        """
        reduction = True
        while reduction:
            reduction = False
            mc = x._monomial_coefficients
            for m, c in mc.items():
                d = m.dict()
                for k, e, g in zip(self._key_order, self._degrees, self._groebner_basis):
                    if k not in d:
                        continue
                    if d[k] >= e:
                        d[k] -= e
                        x -= self._pbw.monomial(self._pbw._indices(d)) * g
                        reduction = True
                        break
        data = {}
        for m, c in x._monomial_coefficients.items():
            d = m.dict()
            data[self._indices([d.get(k, 0) for k in self._key_order])] = c
        return self.element_class(self, data)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: sl2 = LieAlgebra(GF(3), cartan_type=['A',1])
                sage: F = sl2.faithful_representation()
                sage: v = F.an_element(); v
                1 + 2*f2 + f0*f1*f2
                sage: sl2.an_element() * v
                f2 + 2*f2^2 + f1 + 2*f1*f2 + 2*f1^2*f2 + f0 + 2*f0*f2 + 2*f0*f2^2
                 + 2*f0*f1*f2 + f0*f1*f2^2 + f0*f1^2*f2 + f0^2*f1*f2
                sage: sl2.pbw_basis().an_element() * v
                1 + 2*f2 + 2*f0 + f0*f2 + f0*f1*f2 + 2*f0^2*f1*f2
                sage: 5 * v
                2 + f2 + 2*f0*f1*f2
                sage: v * 5
                2 + f2 + 2*f0*f1*f2
                sage: v._acted_upon_(sl2.an_element(), True) is None
                True
            """
            P = self.parent()
            if scalar in P._lie_algebra or scalar in P._pbw:
                if self_on_left:
                    return None
                if not self:  # we are (already) the zero vector
                    return self
                scalar = P._pbw(scalar)
                monoid = P._pbw._indices
                I = P._key_order
                lift = P._pbw.element_class(P._pbw, {monoid(list(zip(I, m.exponents()))): coeff
                                                     for m, coeff in self._monomial_coefficients.items()})
                return P._project(scalar * lift)

            return super()._acted_upon_(scalar, self_on_left)
