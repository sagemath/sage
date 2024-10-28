# sage.doctest: needs sage.combinat sage.modules
r"""
Quantum Matrix Coordinate Algebras

AUTHORS:

- Travis Scrimshaw (01-2016): initial version
"""

##############################################################################
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.family import Family
from sage.categories.algebras import Algebras
from sage.categories.bialgebras import Bialgebras
from sage.categories.hopf_algebras import HopfAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import ZZ


class QuantumMatrixCoordinateAlgebra_abstract(CombinatorialFreeModule):
    """
    Abstract base class for quantum coordinate algebras of a set
    of matrices.
    """
    @staticmethod
    def __classcall__(cls, q=None, bar=None, R=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: O1 = algebras.QuantumMatrixCoordinate(4)
            sage: O2 = algebras.QuantumMatrixCoordinate(4, 4, q=q)
            sage: O3 = algebras.QuantumMatrixCoordinate(4, R=ZZ)
            sage: O4 = algebras.QuantumMatrixCoordinate(4, R=R, q=q)
            sage: O1 is O2 and O2 is O3 and O3 is O4
            True
            sage: O5 = algebras.QuantumMatrixCoordinate(4, R=QQ)
            sage: O1 is O5
            False
        """
        if R is None:
            R = ZZ
        else:
            if q is not None:
                q = R(q)
        if q is None:
            q = LaurentPolynomialRing(R, 'q').gen()
        return super().__classcall__(cls,
                                     q=q, bar=bar, R=q.parent(), **kwds)

    def __init__(self, gp_indices, n, q, bar, R, category, indices_key=None):
        """
        Initialize ``self``.

        TESTS::

            sage: O = algebras.QuantumMatrixCoordinate(3, 2)
            sage: TestSuite(O).run()
        """
        self._n = n
        self._q = q
        if bar is None:
            def bar(x):
                return x.subs(q=~self._q)
        self._bar = bar
        if indices_key is None:
            indices = IndexedFreeAbelianMonoid(gp_indices)
        else:
            indices = IndexedFreeAbelianMonoid(gp_indices, sorting_key=indices_key)
        CombinatorialFreeModule.__init__(self, R, indices, category=category)

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: I = O.indices()
            sage: x = I.an_element(); x
            F[(1, 1)]^2*F[(1, 2)]^2*F[(1, 3)]^3
            sage: O._repr_term(x)
            'x[1,1]^2*x[1,2]^2*x[1,3]^3'
            sage: O._repr_term(I.one())
            '1'
            sage: O.q() * O.one()
            q
        """
        S = m._sorted_items()
        if not S:
            return '1'

        def exp(e):
            return '^{}'.format(e) if e > 1 else ''
        return '*'.join(('x[{},{}]'.format(*k) if k != 'c' else 'c') + exp(e)
                        for k, e in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a latex representation of the term indexed by ``m``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: I = O.indices()
            sage: x = I.an_element(); x
            F[(1, 1)]^2*F[(1, 2)]^2*F[(1, 3)]^3
            sage: O._latex_term(x)
            'x_{1,1}^{2} x_{1,2}^{2} x_{1,3}^{3}'
            sage: O._latex_term(I.one())
            '1'
            sage: latex(O.q() * O.one())
            q
        """
        S = m._sorted_items()
        if not S:
            return '1'

        def exp(e):
            return '^{{{}}}'.format(e) if e > 1 else ''
        return ' '.join(('x_{{{},{}}}'.format(*k) if k != 'c' else 'c') + exp(e)
                        for k, e in m._sorted_items())

    def n(self):
        """
        Return the value `n`.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O.n()
            4
            sage: O = algebras.QuantumMatrixCoordinate(4, 6)
            sage: O.n()
            6
        """
        return self._n

    def q(self):
        """
        Return the variable ``q``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O.q()
            q
            sage: O.q().parent()
            Univariate Laurent Polynomial Ring in q over Integer Ring
            sage: O.q().parent() is O.base_ring()
            True
        """
        return self._q

    @cached_method
    def one_basis(self):
        """
        Return the basis element indexing `1`.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O.one_basis()
            1
            sage: O.one()
            1

        TESTS::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O.one_basis() == O.indices().one()
            True
        """
        return self._indices.one()

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self`` as a tuple.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(3)
            sage: O.gens()
            (x[1,1], x[1,2], x[1,3],
             x[2,1], x[2,2], x[2,3],
             x[3,1], x[3,2], x[3,3])
        """
        return tuple(self.algebra_generators())

    @cached_method
    def quantum_determinant(self):
        r"""
        Return the quantum determinant of ``self``.

        The quantum determinant is defined by

        .. MATH::

            \det_q = \sum_{\sigma \in S_n} (-q)^{\ell(\sigma)}
            x_{1, \sigma(1)} x_{2, \sigma(2)} \cdots x_{n, \sigma(n)}.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(2)
            sage: O.quantum_determinant()
            x[1,1]*x[2,2] - q*x[1,2]*x[2,1]

        We verify that the quantum determinant is central::

            sage: for n in range(2,5):
            ....:     O = algebras.QuantumMatrixCoordinate(n)
            ....:     qdet = O.quantum_determinant()
            ....:     assert all(g * qdet == qdet * g for g in O.algebra_generators())

        We also verify that it is group-like::

            sage: for n in range(2,4):
            ....:     O = algebras.QuantumMatrixCoordinate(n)
            ....:     qdet = O.quantum_determinant()
            ....:     assert qdet.coproduct() == tensor([qdet, qdet])
        """
        if hasattr(self, '_m') and self._m != self._n:
            raise ValueError("undefined for non-square quantum matrices")
        from sage.combinat.permutation import Permutations
        q = self._q
        return self._from_dict({self._indices({(i, p(i)): 1 for i in range(1, self._n + 1)}):
                               (-q) ** p.length() for p in Permutations(self._n)})

    def product_on_basis(self, a, b):
        """
        Return the product of basis elements indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: x = O.algebra_generators()
            sage: b = x[1,4] * x[2,1] * x[3,4]  # indirect doctest
            sage: b * (b * b) == (b * b) * b
            True
            sage: p = prod(list(O.algebra_generators())[:10])
            sage: p * (p * p) == (p * p) * p  # long time
            True
            sage: x = O.an_element()
            sage: y = x^2 + x[4,4] * x[3,3] * x[1,2]
            sage: z = x[2,2] * x[1,4] * x[3,4] * x[1,1]
            sage: x * (y * z) == (x * y) * z
            True
        """
        al = a._sorted_items()
        bl = b._sorted_items()
        # Check for multiplication by 1
        if not al:
            return self.monomial(b)
        if not bl:
            return self.monomial(a)
        if al[-1][0] < bl[0][0]:  # Already in order
            return self.monomial(a * b)
        G = self._indices.monoid_generators()
        one = self.base_ring().one()
        q = self._q
        qi = q ** -1
        monomial = b
        coeff = one
        for pos in range(len(al) - 1, -1, -1):
            ax, ae = al[pos]
            for bx, be in bl:
                if ax[0] < bx[0]:
                    # In order, so nothing more to do
                    break
                elif ax[0] == bx[0]:
                    if ax[1] > bx[1]:
                        # x_{it} x_{ij} = q^{-1} x_{ij} x_{it} if t < j
                        coeff *= qi ** (ae * be)
                    else:
                        # In order, so nothing more to do
                        break
                elif ax[1] == bx[1]:
                    # x_{sj} x_{ij} = q^{-1} x_{ij} x_{sj} if s > i
                    coeff *= qi ** (ae * be)
                elif ax[1] > bx[1]:  # By this point, we must have ax[0] > bx[0]
                    # x_{st} x_{ij} = x_{ij} x_{st} + (q^-1 - q) x_{it} x_{sj}
                    # if s > i, t > j

                    # By Lemma 2.7 (with fixed typo) in H. Zhang and R.B. Zhang:
                    # x_{st} x_{ij}^k = x_{ij}^k x_{st}
                    #                  + (q^{1-2k} - q) x_{ij}^{k-1} x_{it} x_{sj}
                    m1 = G[bx] ** be * G[ax]
                    m2 = G[bx] ** (be - 1) * G[(bx[0], ax[1])] * G[(ax[0], bx[1])]
                    ret = self._from_dict({m1: one, m2: (q ** (1 - 2 * be) - q)})
                    ml = monomial._sorted_items()
                    index = ml.index((bx, be))
                    a_key = self._indices(dict(al[:pos]))
                    bp_key = self._indices(dict(ml[:index])) * G[ax] ** (ae - 1)
                    return (self.monomial(a_key) *
                            self.monomial(bp_key) *
                            ret *
                            self.term(self._indices(dict(ml[index + 1:])),
                                      coeff))

                # Otherwise ax[1] > bx[1], but for this case they commute:
                # x_{st} x_{ij} = x_{ij} x_{st} if s > i, t < j
                # So there is nothing to do to coeff
            monomial *= G[ax] ** ae
        return self.term(monomial, coeff)

    @cached_method
    def _bar_on_basis(self, x):
        """
        Return the bar involution on the basis element indexed by ``x``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O._bar_on_basis(O._indices.an_element())
            (q^-16)*x[1,1]^2*x[1,2]^2*x[1,3]^3
        """
        ret = self.one()
        for k, e in reversed(x._sorted_items()):
            ret *= self.monomial(self._indices({k: e}))
        return ret

    def counit_on_basis(self, x):
        r"""
        Return the counit on the basis element indexed by ``x``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: G = O.algebra_generators()
            sage: I = [1,2,3,4]
            sage: matrix([[G[i,j].counit() for i in I] for j in I])  # indirect doctest
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        if all(t == 'c' or t[0] == t[1] for t, e in x._sorted_items()):
            return self.base_ring().one()
        else:
            return self.base_ring().zero()

    class Element(CombinatorialFreeModule.Element):
        """
        An element of a quantum matrix coordinate algebra.
        """
        def bar(self):
            r"""
            Return the image of ``self`` under the bar involution.

            The bar involution is the `\QQ`-algebra anti-automorphism
            defined by `x_{ij} \mapsto x_{ji}` and `q \mapsto q^{-1}`.

            EXAMPLES::

                sage: O = algebras.QuantumMatrixCoordinate(4)
                sage: x = O.an_element()
                sage: x.bar()
                1 + 2*x[1,1] + (q^-16)*x[1,1]^2*x[1,2]^2*x[1,3]^3 + 3*x[1,2]
                sage: x = O.an_element() * O.algebra_generators()[2,4]; x
                x[1,1]^2*x[1,2]^2*x[1,3]^3*x[2,4] + 2*x[1,1]*x[2,4]
                 + 3*x[1,2]*x[2,4] + x[2,4]
                sage: xb = x.bar(); xb
                (q^-16)*x[1,1]^2*x[1,2]^2*x[1,3]^3*x[2,4]
                 + (q^-21-q^-15)*x[1,1]^2*x[1,2]^2*x[1,3]^2*x[1,4]*x[2,3]
                 + (q^-22-q^-18)*x[1,1]^2*x[1,2]*x[1,3]^3*x[1,4]*x[2,2]
                 + (q^-24-q^-20)*x[1,1]*x[1,2]^2*x[1,3]^3*x[1,4]*x[2,1]
                 + 2*x[1,1]*x[2,4] + 3*x[1,2]*x[2,4]
                 + (2*q^-1-2*q)*x[1,4]*x[2,1]
                 + (3*q^-1-3*q)*x[1,4]*x[2,2] + x[2,4]
                sage: xb.bar() == x
                True
            """
            P = self.parent()
            return P.sum(P._bar(c) * P._bar_on_basis(m) for m, c in self)


class QuantumMatrixCoordinateAlgebra(QuantumMatrixCoordinateAlgebra_abstract):
    r"""
    A quantum matrix coordinate algebra.

    Let `R` be a commutative ring. The quantum matrix coordinate algebra
    of `M(m, n)` is the associative algebra over `R[q, q^{-1}]`
    generated by `x_{ij}`, for `i = 1, 2, \ldots, m`, `j = 1, 2, \ldots, n`,
    and subject to the following relations:

    .. MATH::

        \begin{array}{ll}
        x_{it} x_{ij} = q^{-1} x_{ij} x_{it}   & \text{if } j < t, \\
        x_{sj} x_{ij} = q^{-1} x_{ij} x_{sj}   & \text{if } i < s, \\
        x_{st} x_{ij} = x_{ij} x_{st}          & \text{if } i < s, j > t, \\
        x_{st} x_{ij} = x_{ij} x_{st} + (q^{-1} - q) x_{it} x_{sj}
                                               & \text{if } i < s, j < t. \\
        \end{array}

    The quantum matrix coordinate algebra is denoted by
    `\mathcal{O}_q(M(m, n))`. For `m = n`, it is also a bialgebra given by

    .. MATH::

        \Delta(x_{ij}) = \sum_{k=1}^n x_{ik} \otimes x_{kj},
        \varepsilon(x_{ij}) = \delta_{ij}.

    Moreover, there is a central group-like element called the
    *quantum determinant* that is defined by

    .. MATH::

        \det_q = \sum_{\sigma \in S_n} (-q)^{\ell(\sigma)}
        x_{1,\sigma(1)} x_{2,\sigma(2)} \cdots x_{n,\sigma(n)}.

    The quantum matrix coordinate algebra also has natural inclusions
    when restricting to submatrices. That is, let
    `I \subseteq \{1, 2, \ldots, m\}` and `J \subseteq \{1, 2, \ldots, n\}`.
    Then the subalgebra generated by `\{ x_{ij} \mid i \in I, j \in J \}`
    is naturally isomorphic to `\mathcal{O}_q(M(|I|, |J|))`.

    .. NOTE::

        The `q` considered here is `q^2` in some references, e.g., [ZZ2005]_.

    INPUT:

    - ``m`` -- the integer `m`
    - ``n`` -- the integer `n`
    - ``R`` -- (optional) the ring `R` if `q` is not specified
      (the default is `\ZZ`); otherwise the ring containing `q`
    - ``q`` -- (optional) the variable `q`; the default is
      `q \in R[q, q^{-1}]`
    - ``bar`` -- (optional) the involution on the base ring; the
      default is `q \mapsto q^{-1}`

    EXAMPLES:

    We construct `\mathcal{O}_q(M(2,3))` and the variables::

        sage: O = algebras.QuantumMatrixCoordinate(2,3)
        sage: O.inject_variables()
        Defining x11, x12, x13, x21, x22, x23

    We do some basic computations::

        sage: x21 * x11
        (q^-1)*x[1,1]*x[2,1]
        sage: x23 * x12 * x11
        (q^-1)*x[1,1]*x[1,2]*x[2,3] + (q^-2-1)*x[1,1]*x[1,3]*x[2,2]
         + (q^-3-q^-1)*x[1,2]*x[1,3]*x[2,1]

    We construct the maximal quantum minors::

        sage: q = O.q()
        sage: qm12 = x11*x22 - q*x12*x21
        sage: qm13 = x11*x23 - q*x13*x21
        sage: qm23 = x12*x23 - q*x13*x22

    However, unlike for the quantum determinant, they are not central::

        sage: all(qm12 * g == g * qm12 for g in O.algebra_generators())
        False
        sage: all(qm13 * g == g * qm13 for g in O.algebra_generators())
        False
        sage: all(qm23 * g == g * qm23 for g in O.algebra_generators())
        False

    REFERENCES:

    - [FRT1990]_
    - [ZZ2005]_
    """
    @staticmethod
    def __classcall_private__(cls, m, n=None, q=None, bar=None, R=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: O1 = algebras.QuantumMatrixCoordinate(4)
            sage: O2 = algebras.QuantumMatrixCoordinate(4, 4, q=q)
            sage: O3 = algebras.QuantumMatrixCoordinate(4, R=ZZ)
            sage: O4 = algebras.QuantumMatrixCoordinate(4, R=R, q=q)
            sage: O1 is O2 and O2 is O3 and O3 is O4
            True
            sage: O5 = algebras.QuantumMatrixCoordinate(4, R=QQ)
            sage: O1 is O5
            False
        """
        if n is None:
            n = m
        return super().__classcall__(cls, m=m, n=n,
                                     q=q, bar=bar,
                                     R=R)

    def __init__(self, m, n, q, bar, R):
        """
        Initialize ``self``.

        TESTS::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: TestSuite(O).run()

            sage: O = algebras.QuantumMatrixCoordinate(10)
            sage: O.variable_names()
            ('x0101', ..., 'x1010')
            sage: O = algebras.QuantumMatrixCoordinate(11,3)
            sage: O.variable_names()
            ('x011', ..., 'x113')
            sage: O = algebras.QuantumMatrixCoordinate(3,11)
            sage: O.variable_names()
            ('x101', ..., 'x311')
        """
        gp_indices = [(i, j) for i in range(1, m + 1) for j in range(1, n + 1)]

        if m == n:
            cat = Bialgebras(R.category()).WithBasis()
        else:
            cat = Algebras(R.category()).WithBasis()

        self._m = m
        QuantumMatrixCoordinateAlgebra_abstract.__init__(self, gp_indices, n, q, bar, R, cat)
        # Set the names
        mb = len(str(m))
        nb = len(str(n))
        base = 'x{{:0>{}}}{{:0>{}}}'.format(mb,nb)
        names = [base.format(*k) for k in gp_indices]
        self._assign_names(names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.QuantumMatrixCoordinate(4)
            Quantized coordinate algebra of M(4, 4) with q=q over
             Univariate Laurent Polynomial Ring in q over Integer Ring

            sage: algebras.QuantumMatrixCoordinate(4, 2)
            Quantized coordinate algebra of M(4, 2) with q=q over
             Univariate Laurent Polynomial Ring in q over Integer Ring
        """
        txt = "Quantized coordinate algebra of M({}, {}) with q={} over {}"
        return txt.format(self._m, self._n, self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: latex(O)
            \mathcal{O}_{q}(M(4, 4))
        """
        return "\\mathcal{O}_{%s}(M(%s, %s))" % (self._q, self._m, self._n)

    def m(self):
        """
        Return the value `m`.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4, 6)
            sage: O.m()
            4
            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: O.m()
            4
        """
        return self._m

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(2)
            sage: O.algebra_generators()
            Finite family {(1, 1): x[1,1], (1, 2): x[1,2], (2, 1): x[2,1], (2, 2): x[2,2]}
        """
        l = [(i, j) for i in range(1, self._m + 1)
             for j in range(1, self._n + 1)]
        G = self._indices.monoid_generators()
        one = self.base_ring().one()
        return Family(l, lambda x: self.element_class(self, {G[x]: one}))

    def coproduct_on_basis(self, x):
        r"""
        Return the coproduct on the basis element indexed by ``x``.

        EXAMPLES::

            sage: O = algebras.QuantumMatrixCoordinate(4)
            sage: x24 = O.algebra_generators()[2,4]
            sage: O.coproduct_on_basis(x24.leading_support())
            x[2,1] # x[1,4] + x[2,2] # x[2,4] + x[2,3] # x[3,4] + x[2,4] # x[4,4]

        TESTS:

        We check that it is an algebra morphism::

            sage: O = algebras.QuantumMatrixCoordinate(3)
            sage: G = O.algebra_generators()
            sage: all(x.coproduct() * y.coproduct() == (x * y).coproduct()
            ....:     for x in G for y in G)
            True
        """
        if self._m != self._n:
            raise ValueError("undefined for non-square quantum matrices")
        T = self.tensor_square()
        I = self._indices.monoid_generators()
        return T.prod(T.sum_of_monomials((I[t[0], k], I[k, t[1]])
                                         for k in range(1, self._n + 1)) ** e
                      for t, e in x._sorted_items())


class QuantumGL(QuantumMatrixCoordinateAlgebra_abstract):
    r"""
    Quantum coordinate algebra of `GL(n)`.

    The quantum coordinate algebra of `GL(n)`, or quantum `GL(n)`
    for short and denoted by `\mathcal{O}_q(GL(n))`, is the quantum
    coordinate algebra of `M_R(n, n)` with the addition of the
    additional central group-like element `c` which satisfies
    `c d = d c = 1`, where `d` is the quantum determinant.

    Quantum `GL(n)` is a Hopf algebra where `\varepsilon(c) = 1`
    and the antipode `S` is given by the (quantum) matrix inverse.
    That is to say, we have `S(c) = c^-1 = d` and

    .. MATH::

        S(x_{ij}) = c * (-q)^{i-j} * \tilde{t}_{ji},

    where we have the quantum minor

    .. MATH::

        \tilde{t}_{ij} = \sum_{\sigma} (-q)^{\ell(\sigma)}
        x_{1, \sigma(1)} \cdots x_{i-1, \sigma(i-1)} x_{i+1, \sigma(i+1)}
        \cdots x_{n, \sigma(n)}

    with the sum over permutations `\sigma \colon \{1, \ldots, i-1, i+1,
    \ldots n\} \to \{1, \ldots, j-1, j+1, \ldots, n\}`.

    .. SEEALSO::

        :class:`QuantumMatrixCoordinateAlgebra`

    INPUT:

    - ``n`` -- the integer `n`
    - ``R`` -- (optional) the ring `R` if `q` is not specified
      (the default is `\ZZ`); otherwise the ring containing `q`
    - ``q`` -- (optional) the variable `q`; the default is
      `q \in R[q, q^{-1}]`
    - ``bar`` -- (optional) the involution on the base ring; the
      default is `q \mapsto q^{-1}`

    EXAMPLES:

    We construct `\mathcal{O}_q(GL(3))` and the variables::

        sage: O = algebras.QuantumGL(3)
        sage: O.inject_variables()
        Defining x11, x12, x13, x21, x22, x23, x31, x32, x33, c

    We do some basic computations::

        sage: x33 * x12
        x[1,2]*x[3,3] + (q^-1-q)*x[1,3]*x[3,2]
        sage: x23 * x12 * x11
        (q^-1)*x[1,1]*x[1,2]*x[2,3] + (q^-2-1)*x[1,1]*x[1,3]*x[2,2]
         + (q^-3-q^-1)*x[1,2]*x[1,3]*x[2,1]
        sage: c * O.quantum_determinant()
        1

    We verify the quantum determinant is in the center and is group-like::

        sage: qdet = O.quantum_determinant()
        sage: all(qdet * g == g * qdet for g in O.algebra_generators())
        True
        sage: qdet.coproduct() == tensor([qdet, qdet])
        True

    We check that the inverse of the quantum determinant is also in
    the center and group-like::

        sage: all(c * g == g * c for g in O.algebra_generators())
        True
        sage: c.coproduct() == tensor([c, c])
        True

    Moreover, the antipode interchanges the quantum determinant and
    its inverse::

        sage: c.antipode() == qdet
        True
        sage: qdet.antipode() == c
        True

    REFERENCES:

    - [DD1991]_
    - [Kar1993]_
    """
    @staticmethod
    def __classcall_private__(cls, n, q=None, bar=None, R=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: O1 = algebras.QuantumGL(4)
            sage: O2 = algebras.QuantumGL(4, R=ZZ)
            sage: O3 = algebras.QuantumGL(4, R=R, q=q)
            sage: O1 is O2 and O2 is O3
            True
            sage: O4 = algebras.QuantumGL(4, R=QQ)
            sage: O1 is O4
            False
        """
        return super().__classcall__(cls, n=n, q=q, bar=bar, R=R)

    def __init__(self, n, q, bar, R):
        """
        Initialize ``self``.

        TESTS::

            sage: O = algebras.QuantumGL(2)
            sage: elts = list(O.algebra_generators())
            sage: elts += [O.quantum_determinant(), O.an_element()]
            sage: TestSuite(O).run(elements=elts) # long time
        """
        # Set the names
        gp_indices = [(i, j) for i in range(1, n + 1) for j in range(1, n + 1)]
        gp_indices.append('c')
        cat = HopfAlgebras(R.category()).WithBasis()
        QuantumMatrixCoordinateAlgebra_abstract.__init__(self, gp_indices, n, q,
                                                         bar, R, cat,
                                                         indices_key=_generator_key)
        names = ['x{}{}'.format(*k) for k in gp_indices[:-1]]
        names.append('c')
        self._assign_names(names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.QuantumGL(4)
            Quantized coordinate algebra of GL(4) with q=q over
             Univariate Laurent Polynomial Ring in q over Integer Ring
        """
        txt = "Quantized coordinate algebra of GL({}) with q={} over {}"
        return txt.format(self._n, self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(4)
            sage: latex(O)
            \mathcal{O}_{q}(GL(4))
        """
        return "\\mathcal{O}_{%s}(GL(%s))" % (self._q, self._n)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(2)
            sage: O.algebra_generators()
            Finite family {(1, 1): x[1,1], (1, 2): x[1,2], (2, 1): x[2,1], (2, 2): x[2,2], 'c': c}
        """
        l = [(i, j) for i in range(1, self._n + 1)
             for j in range(1, self._n + 1)]
        l.append('c')
        G = self._indices.monoid_generators()
        one = self.base_ring().one()
        return Family(l, lambda x: self.element_class(self, {G[x]: one}))

    @lazy_attribute
    def _qdet_cancel_monomial(self):
        """
        Return the trailing monomial of the quantum determinant.

        EXAMPLES::

            sage: O = algebras.QuantumGL(2)
            sage: O._qdet_cancel_monomial
            F[(1, 1)]*F[(2, 2)]
        """
        I = self._indices
        gens = I.monoid_generators()
        return I.prod(gens[i, i] for i in range(1, self._n + 1))

    @lazy_attribute
    def _qdet_remaining(self):
        r"""
        Return the remaining terms when cancelling the leading term.

        Consider `d = m + L`, where `m` is the leading term of the
        quantum determinant `d`. Then we have `c d = cm + cL = 1`,
        which we rewrite as `cm = 1 - cL`. This lazy attribute
        is `1 - cL`.

        EXAMPLES::

            sage: O = algebras.QuantumGL(2)
            sage: O._qdet_remaining
            1 + q*c*x[1,2]*x[2,1]
        """
        temp = self.monomial(self._qdet_cancel_monomial) - self.quantum_determinant()
        c = self._indices.monoid_generators()['c']
        ret = {c * mon: coeff for mon, coeff in temp}
        return self._from_dict(ret, remove_zeros=False) + self.one()

    def product_on_basis(self, a, b):
        r"""
        Return the product of basis elements indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(2)
            sage: I = O.indices().monoid_generators()
            sage: O.product_on_basis(I[1,1], I[2,2])
            x[1,1]*x[2,2]
            sage: O.product_on_basis(I[2,2], I[1,1])
            x[1,1]*x[2,2] + (q^-1-q)*x[1,2]*x[2,1]

        TESTS::

            sage: x11,x12,x21,x22,c = O.algebra_generators()
            sage: x11 * x22
            x[1,1]*x[2,2]
            sage: x22 * x12
            (q^-1)*x[1,2]*x[2,2]
            sage: x22 * x11
            x[1,1]*x[2,2] + (q^-1-q)*x[1,2]*x[2,1]
            sage: c * (x11 * O.quantum_determinant())
            x[1,1]
        """
        I = self._indices
        c_exp = 0
        if 'c' in a._monomial:
            da = dict(a._monomial)  # Make a copy
            c_exp += da.pop('c')
            a = I(da)
        if 'c' in b._monomial:
            db = dict(b._monomial)  # Make a copy
            c_exp += db.pop('c')
            b = I(db)
        # a and b contain no powers of c
        p = super().product_on_basis(a, b)
        if c_exp == 0:
            return p
        c = self._indices.monoid_generators()['c']
        ret = {}
        other = self.zero()
        for mon, coeff in p:
            try:
                # Given that cz = R and we have a monomial ab, we need to
                #   rewrite zx in terms of ab plus lower order terms L:
                # zx = X * ab + L
                # c * zx = R * x = c * X * ab + c * L
                # c * ab = (R * x - c * L) / X
                rem = self.monomial(mon // self._qdet_cancel_monomial)
                L = self.monomial(self._qdet_cancel_monomial) * rem
                co = L[mon]
                del L._monomial_coefficients[mon]
                temp = self.term(c ** (c_exp - 1), coeff) * self._qdet_remaining * rem
                if L != self.zero():
                    temp -= self.term(c ** c_exp, coeff) * L
                for k in temp._monomial_coefficients:
                    temp._monomial_coefficients[k] //= co
                other += temp
            except ValueError:  # We cannot cancel, so we just add on the correct power of c
                ret[c ** c_exp * mon] = coeff
        return self._from_dict(ret, remove_zeros=False) + other

    @cached_method
    def _antipode_on_generator(self, i, j):
        """
        Return the antipode on the generator indexed by ``(i, j)``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(2)
            sage: [[O._antipode_on_generator(i, j) for i in [1,2]] for j in [1,2]]
            [[c*x[2,2], -q*c*x[2,1]],
             [-(q^-1)*c*x[1,2], c*x[1,1]]]
        """
        from sage.combinat.permutation import Permutations
        q = self._q
        I = list(range(1, j)) + list(range(j + 1, self._n + 1))

        def lift(p):
            return [val if val < i else val + 1 for val in p]
        gens = self.algebra_generators()
        t_tilde = self.sum((-q) ** p.length() * gens['c'] *
                           self.prod(gens[I[k], val]
                                     for k, val in enumerate(lift(p)))
                           for p in Permutations(self._n - 1))
        return (-q) ** (i - j) * t_tilde

    def antipode_on_basis(self, x):
        r"""
        Return the antipode of the basis element indexed by ``x``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(3)
            sage: x = O.indices().monoid_generators()
            sage: O.antipode_on_basis(x[1,2])
            -(q^-1)*c*x[1,2]*x[3,3] + c*x[1,3]*x[3,2]
            sage: O.antipode_on_basis(x[2,2])
            c*x[1,1]*x[3,3] - q*c*x[1,3]*x[3,1]
            sage: O.antipode_on_basis(x['c']) == O.quantum_determinant()
            True
        """
        ret = self.one()
        for k, e in reversed(x._sorted_items()):
            if k == 'c':
                ret *= self.quantum_determinant() ** e
            else:
                ret *= self._antipode_on_generator(*k) ** e
        return ret

    def coproduct_on_basis(self, x):
        r"""
        Return the coproduct on the basis element indexed by ``x``.

        EXAMPLES::

            sage: O = algebras.QuantumGL(3)
            sage: x = O.indices().monoid_generators()
            sage: O.coproduct_on_basis(x[1,2])
            x[1,1] # x[1,2] + x[1,2] # x[2,2] + x[1,3] # x[3,2]
            sage: O.coproduct_on_basis(x[2,2])
            x[2,1] # x[1,2] + x[2,2] # x[2,2] + x[2,3] # x[3,2]
            sage: O.coproduct_on_basis(x['c'])
            c # c
        """
        T = self.tensor_square()
        I = self._indices.monoid_generators()
        return T.prod(T.sum_of_monomials((I[t[0], k], I[k, t[1]])
                                         for k in range(1, self._n + 1)) ** e
                      if t != 'c' else T.monomial((I['c'], I['c'])) ** e
                      for t, e in x._sorted_items())


def _generator_key(t):
    """
    Helper function to make ``'c'`` less that all other indices for
    sorting the monomials in :class:`QuantumGL`.

    INPUT:

    - ``t`` -- tuple (index, exponent)

    OUTPUT: a tuple made from the index only

    EXAMPLES::

        sage: from sage.algebras.quantum_matrix_coordinate_algebra import _generator_key as k
        sage: k(((1,2),1)) < k(('c',1))
        False
        sage: k(((1,2),1)) < k(((1,3),1))
        True
        sage: k(((1,2),1)) < k(((3,1),1))
        True
        sage: k(('c',2)) < k(((1,1),1))
        True
    """
    t = t[0]
    if isinstance(t, tuple):
        return t
    return ()
