"""
Dual Symmetric Functions in Non-Commuting Variables

AUTHORS:

- Travis Scrimshaw (08-04-2013): Initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras

from sage.combinat.ncsym.bases import NCSymDualBases, NCSymBasis_abstract
from sage.combinat.partition import Partition
from sage.combinat.set_partition import SetPartitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.subset import Subsets
from sage.functions.other import factorial
from sage.sets.set import Set

class SymmetricFunctionsNonCommutingVariablesDual(UniqueRepresentation, Parent):
    r"""
    The Hopf dual to the symmetric functions in non-commuting variables.

    See Section 2.3 of [BZ05]_ for a study.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(SymmetricFunctionsNonCommutingVariables(QQ).dual()).run()
        """
        self._base = R # Won't be needed once CategoryObject won't override base_ring
        category = GradedHopfAlgebras(R)  # TODO: .Commutative()
        Parent.__init__(self, category=category.WithRealizations())

        # Bases
        w = self.w()

        # Embedding of Sym in the homogeneous bases into DNCSym in the w basis
        Sym = SymmetricFunctions(self.base_ring())
        Sym_h_to_w = Sym.h().module_morphism(w.sum_of_partitions,
                                             triangular='lower',
                                             inverse_on_support=w._set_par_to_par,
                                             codomain=w, category=category)
        Sym_h_to_w.register_as_coercion()
        self.to_symmetric_function = Sym_h_to_w.section()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctionsNonCommutingVariables(ZZ).dual()
            Dual symmetric functions in non-commuting variables over the Integer Ring
        """
        return "Dual symmetric functions in non-commuting variables over the %s"%self.base_ring()

    def a_realization(self):
        r"""
        Return the realization of the `\mathbf{w}` basis of ``self``.

        EXAMPLES::

            sage: SymmetricFunctionsNonCommutingVariables(QQ).dual().a_realization()
            Dual symmetric functions in non-commuting variables over the Rational Field in the w basis
        """
        return self.w()

    _shorthands = tuple(['w'])

    def dual(self):
        r"""
        Return the dual Hopf algebra of the dual symmetric functions in
        non-commuting variables.

        EXAMPLES::

            sage: NCSymD = SymmetricFunctionsNonCommutingVariables(QQ).dual()
            sage: NCSymD.dual()
            Symmetric functions in non-commuting variables over the Rational Field
        """
        from sage.combinat.ncsym.ncsym import SymmetricFunctionsNonCommutingVariables
        return SymmetricFunctionsNonCommutingVariables(self.base_ring())

    class w(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the `\mathbf{w}` basis.

        EXAMPLES::

            sage: NCSymD = SymmetricFunctionsNonCommutingVariables(QQ).dual()
            sage: w = NCSymD.w()

        We have the embedding `\chi^*` of `Sym` into `NCSym^*` available as
        a coercion::

            sage: h = SymmetricFunctions(QQ).h()
            sage: w(h[2,1])
            w{{1}, {2, 3}} + w{{1, 2}, {3}} + w{{1, 3}, {2}}

        Similarly we can pull back when we are in the image of `\chi^*`::

            sage: elt = 3*(w[[1],[2,3]] + w[[1,2],[3]] + w[[1,3],[2]])
            sage: h(elt)
            3*h[2, 1]
        """
        def __init__(self, NCSymD):
            """
            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: TestSuite(w).run()
            """
            def lt_set_part(A, B):
                A = sorted(map(sorted, A))
                B = sorted(map(sorted, B))
                for i in range(len(A)):
                    if A[i] > B[i]:
                        return 1
                    elif A[i] < B[i]:
                        return -1
                return 0
            CombinatorialFreeModule.__init__(self, NCSymD.base_ring(), SetPartitions(),
                                             prefix='w', bracket=False,
                                             monomial_cmp=lt_set_part,
                                             category=NCSymDualBases(NCSymD))

        @lazy_attribute
        def to_symmetric_function(self):
            r"""
            The preimage of `\chi^*` in the `\mathbf{w}` basis.

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.to_symmetric_function
                Generic morphism:
                  From: Dual symmetric functions in non-commuting variables over the Rational Field in the w basis
                  To:   Symmetric Functions over Rational Field in the homogeneous basis
            """
            return self.realization_of().to_symmetric_function

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{w}` basis.

            The dual basis to the `\mathbf{w}` basis is the monomial basis
            of the symmetric functions in non-commuting variables.

            OUTPUT:

            - the monomial basis of the symmetric functions in non-commuting variables

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.dual_basis()
                Symmetric functions in non-commuting variables over the Rational Field in the monomial basis
            """
            return self.realization_of().dual().m()

        def product_on_basis(self, A, B):
            r"""
            The product on `\mathbf{w}` basis elements.

            The product on the `\mathbf{w}` is the dual to the coproduct on the
            `\mathbf{m}` basis.  On the basis `\mathbf{w}` it is defined as

            .. MATH::

                \mathbf{w}_A \mathbf{w}_B = \sum_{S \subseteq [n]}
                \mathbf{w}_{A\uparrow_S \cup B\uparrow_{S^c}}

            where the sum is over all possible subsets `S` of `[n]` such that
            `|S| = |A|` with a term indexed the union of `A \uparrow_S` and
            `B \uparrow_{S^c}`. The notation `A \uparrow_S` represents the
            unique set partition of the set `S` such that the standardization
            is `A`.  This product is commutative.

            INPUT:

            - ``A``, ``B`` -- set partitions

            OUTPUT:

            - an element of the `\mathbf{w}` basis

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: A = SetPartition([[1], [2,3]])
                sage: B = SetPartition([[1, 2, 3]])
                sage: w.product_on_basis(A, B)
                w{{1}, {2, 3}, {4, 5, 6}} + w{{1}, {2, 3, 4}, {5, 6}}
                 + w{{1}, {2, 3, 5}, {4, 6}} + w{{1}, {2, 3, 6}, {4, 5}}
                 + w{{1}, {2, 4}, {3, 5, 6}} + w{{1}, {2, 4, 5}, {3, 6}}
                 + w{{1}, {2, 4, 6}, {3, 5}} + w{{1}, {2, 5}, {3, 4, 6}}
                 + w{{1}, {2, 5, 6}, {3, 4}} + w{{1}, {2, 6}, {3, 4, 5}}
                 + w{{1, 2, 3}, {4}, {5, 6}} + w{{1, 2, 4}, {3}, {5, 6}}
                 + w{{1, 2, 5}, {3}, {4, 6}} + w{{1, 2, 6}, {3}, {4, 5}}
                 + w{{1, 3, 4}, {2}, {5, 6}} + w{{1, 3, 5}, {2}, {4, 6}}
                 + w{{1, 3, 6}, {2}, {4, 5}} + w{{1, 4, 5}, {2}, {3, 6}}
                 + w{{1, 4, 6}, {2}, {3, 5}} + w{{1, 5, 6}, {2}, {3, 4}}
                sage: B = SetPartition([[1], [2]])
                sage: w.product_on_basis(A, B)
                3*w{{1}, {2}, {3}, {4, 5}} + 2*w{{1}, {2}, {3, 4}, {5}}
                 + 2*w{{1}, {2}, {3, 5}, {4}} + w{{1}, {2, 3}, {4}, {5}}
                 + w{{1}, {2, 4}, {3}, {5}} + w{{1}, {2, 5}, {3}, {4}}
                sage: w.product_on_basis(A, SetPartition([]))
                w{{1}, {2, 3}}
            """
            if len(A) == 0:
                return self.monomial(B)
            if len(B) == 0:
                return self.monomial(A)

            P = SetPartitions()
            n = A.size()
            k = B.size()
            def unions(s):
                a = sorted(s)
                b = sorted(Set(range(1, n+k+1)).difference(s))
                # -1 for indexing
                ret = [[a[i-1] for i in sorted(part)] for part in A]
                ret += [[b[i-1] for i in sorted(part)] for part in B]
                return P(ret)
            return self.sum_of_terms([(unions(s), 1)
                    for s in Subsets(n+k, n)])

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct of a `\mathbf{w}` basis element.

            The coproduct on the basis element `\mathbf{w}_A` is the sum over
            tensor product terms `\mathbf{w}_B \otimes \mathbf{w}_C` where
            `B` is the restriction of `A` to `\{1,2,\ldots,k\}` and `C` is
            the restriction of `A` to `\{k+1, k+2, \ldots, n\}`.

            INPUT:

            - ``A`` -- a set partition

            OUTPUT:

            - The coproduct applied to the `\mathbf{w}` dual symmetric function
              in non-commuting variables indexed by ``A`` expressed in the
              `\mathbf{w}` basis.

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w[[1], [2,3]].coproduct()
                w{} # w{{1}, {2, 3}} + w{{1}} # w{{1, 2}}
                 + w{{1}, {2}} # w{{1}} + w{{1}, {2, 3}} # w{}
                sage: w.coproduct_on_basis(SetPartition([]))
                w{} # w{}
            """
            n = A.size()
            return self.tensor_square().sum_of_terms([
                (( A.restriction(range(1, i+1)).standardization(),
                   A.restriction(range(i+1, n+1)).standardization() ), 1)
                for i in range(n+1)], distinct=True)

        def antipode_on_basis(self, A):
            r"""
            Return the antipode applied to the basis element indexed by ``A``.

            INPUT:

            - ``A`` -- a set partition

            OUTPUT:

            - an element in the basis ``self``

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.antipode_on_basis(SetPartition([[1],[2,3]]))
                -3*w{{1}, {2}, {3}} + w{{1, 2}, {3}} + w{{1, 3}, {2}}
                sage: F = w[[1,3],[5],[2,4]].coproduct()
                sage: F.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            if A.size() == 0:
                return self.one()
            if A.size() == 1:
                return -self(A)
            cpr = self.coproduct_on_basis(A)
            return -sum( c*self.monomial(B1)*self.antipode_on_basis(B2)
                         for ((B1,B2),c) in cpr if B2 != A )

        def duality_pairing(self, x, y):
            r"""
            Compute the pairing between an element of ``self`` and an
            element of the dual.

            INPUT:

            - ``x`` -- an element of the dual of symmetric functions in
              non-commuting variables
            - ``y`` -- an element of the symmetric functions in non-commuting
              variables

            OUTPUT:

            - an element of the base ring of ``self``

            EXAMPLES::

                sage: DNCSym = SymmetricFunctionsNonCommutingVariablesDual(QQ)
                sage: w = DNCSym.w()
                sage: m = w.dual_basis()
                sage: matrix([[w(A).duality_pairing(m(B)) for A in SetPartitions(3)] for B in SetPartitions(3)])
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]
                sage: (w[[1,2],[3]] + 3*w[[1,3],[2]]).duality_pairing(2*m[[1,3],[2]] + m[[1,2,3]] + 2*m[[1,2],[3]])
                8
                sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                sage: matrix([[w(A).duality_pairing(h(B)) for A in SetPartitions(3)] for B in SetPartitions(3)])
                [6 2 2 2 1]
                [2 2 1 1 1]
                [2 1 2 1 1]
                [2 1 1 2 1]
                [1 1 1 1 1]
                sage: (2*w[[1,3],[2]] + w[[1,2,3]] + 2*w[[1,2],[3]]).duality_pairing(h[[1,2],[3]] + 3*h[[1,3],[2]])
                32
            """
            x = self(x)
            y = self.dual_basis()(y)
            return sum(coeff * y[I] for (I, coeff) in x)

        def sum_of_partitions(self, la):
            """
            Return the sum over all sets partitions whose shape is ``la``,
            scaled by `\prod_i m_i!` where `m_i` is the multiplicity
            of `i` in ``la``.

            INPUT:

            - ``la`` -- an integer partition

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.sum_of_partitions([2,1,1])
                2*w{{1}, {2}, {3, 4}} + 2*w{{1}, {2, 3}, {4}} + 2*w{{1}, {2, 4}, {3}}
                 + 2*w{{1, 2}, {3}, {4}} + 2*w{{1, 3}, {2}, {4}} + 2*w{{1, 4}, {2}, {3}}
            """
            la = Partition(la)
            c = prod(map(factorial, la.to_exp()))
            P = SetPartitions()
            return self.sum_of_terms([(P(m), c) for m in SetPartitions(sum(la), la)], distinct=True)

        def _set_par_to_par(self, A):
            r"""
            Return the the shape of ``A`` if ``A`` is the canonical standard
            set partition `A_1 | A_2 | \cdots | A_k` where `|` is the pipe
            operation (see
            :meth:~sage.combinat.set_partition.SetPartition.pipe()` )
            and `A_i = [\lambda_i]` where `\lambda_1 \leq \lambda_2 \leq
            \cdots \leq \lambda_k`. Otherwise, return ``None``.

            This is the trailing term of `h_{\lambda}` mapped by `\chi` to
            the `\mathbf{w}` basis and is used by the coercion framework to
            constuct the preimage `\chi^{-1}`.

            INPUT:

            - ``A`` -- a set partition

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w._set_par_to_par(SetPartition([[1], [2], [3,4,5]]))
                [3, 1, 1]
                sage: w._set_par_to_par(SetPartition([[1,2,3],[4],[5]]))
                sage: w._set_par_to_par(SetPartition([[1],[2,3,4],[5]]))
                sage: w._set_par_to_par(SetPartition([[1],[2,3,5],[4]]))

            TESTS:

            This is used in the coercion between `\mathbf{w}` and the
            homogeneous symmetric functions. ::

                sage: w = SymmetricFunctionsNonCommutingVariablesDual(QQ).w()
                sage: h = SymmetricFunctions(QQ).h()
                sage: h(w[[1,3],[2]])
                Traceback (most recent call last):
                ...
                ValueError: w{{1, 3}, {2}} is not in the image of Generic morphism:
                  From: Symmetric Functions over Rational Field in the homogeneous basis
                  To:   Dual symmetric functions in non-commuting variables over the Rational Field in the w basis
                sage: h(w(h[2,1])) == w(h[2,1]).to_symmetric_function()
                True
            """
            cur = 1
            prev_len = 0
            for p in A:
                if prev_len > len(p) or list(p) != range(cur, cur+len(p)):
                    return None
                prev_len = len(p)
                cur += len(p)
            return A.shape()

        class Element(CombinatorialFreeModule.Element):
            r"""
            An element in the `\mathbf{w}` basis.
            """
            def expand(self, n, letter='x'):
                r"""
                Expand ``self`` written in the `\mathbf{w}` basis in `n^2`
                commuting variables which satisfy the relation
                `x_{ij} x_{ik} = 0` for all `i`, `j`, and `k`.

                The expansion of an element of the `\mathbf{w}` basis is
                given by equations (26) and (55) in [HNT06]_.

                INPUT:

                - ``n`` -- an integer
                - ``letter`` -- (default: ``'x'``) a string

                OUTPUT:

                - The symmetric function of ``self`` expressed in the ``n*n``
                  non-commuting variables described by ``letter``.

                REFERENCES:

                .. [HNT06] F. Hivert, J.-C. Novelli, J.-Y. Thibon.
                   *Commutative combinatorial Hopf algebras*. (2006).
                   :arxiv:`0605262v1`.

                EXAMPLES::

                    sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                    sage: w[[1,3],[2]].expand(4)
                    x02*x11*x20 + x03*x11*x30 + x03*x22*x30 + x13*x22*x31

                One can use a different set of variable by using the
                optional argument ``letter``::

                    sage: w[[1,3],[2]].expand(3, letter='y')
                    y02*y11*y20
                """
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                from sage.combinat.permutation import Permutations

                m = self.parent()
                names = ['{}{}{}'.format(letter, i, j) for i in range(n) for j in range(n)]
                R = PolynomialRing(m.base_ring(), n*n, names)
                x = [[R.gens()[i*n+j] for j in range(n)] for i in range(n)]
                I = R.ideal([x[i][j]*x[i][k] for j in range(n) for k in range(n) for i in range(n)])
                Q = R.quotient(I, names)
                x = [[Q.gens()[i*n+j] for j in range(n)] for i in range(n)]
                P = SetPartitions()

                def on_basis(A):
                    k = A.size()
                    ret = R.zero()
                    if n < k:
                        return ret

                    for p in Permutations(k):
                        if P(p.to_cycles()) == A:
                            # -1 for indexing
                            ret += R.sum(prod(x[I[i]][I[p[i]-1]] for i in range(k))
                                         for I in Subsets(range(n), k))
                    return ret

                return m._apply_module_morphism(self, on_basis, codomain=R)

            def is_symmetric(self):
                r"""
                Determine if a `NCSym^*` function, expressed in the
                `\mathbf{w}` basis, is symmetric.

                A function `f` in the `\mathbf{w}` basis is a symmetric
                function if it is in the image of `\chi^*`. That is to say we
                have

                .. MATH::

                    f = \sum_{\lambda} c_{\lambda} \prod_i m_i(\lambda)!
                    \sum_{\lambda(A) = \lambda} \mathbf{w}_A

                where the second sum is over all set partitions `A` whose
                shape `\lambda(A)` is equal to `\lambda` and `m_i(\mu)` is
                the multiplicity of `i` in the partition `\mu`.

                OUTPUT:

                - ``True`` if `\lambda(A)=\lambda(B)` implies the coefficients of
                  `\mathbf{w}_A` and `\mathbf{w}_B` are equal, ``False`` otherwise

                EXAMPLES::

                    sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                    sage: elt = w.sum_of_partitions([2,1,1])
                    sage: elt.is_symmetric()
                    True
                    sage: elt -= 3*w.sum_of_partitions([1,1])
                    sage: elt.is_symmetric()
                    True
                    sage: w = SymmetricFunctionsNonCommutingVariables(ZZ).dual().w()
                    sage: elt = w.sum_of_partitions([2,1,1]) / 2
                    sage: elt.is_symmetric()
                    False
                    sage: elt = w[[1,3],[2]]
                    sage: elt.is_symmetric()
                    False
                    sage: elt = w[[1],[2,3]] + w[[1,2],[3]] + 2*w[[1,3],[2]]
                    sage: elt.is_symmetric()
                    False
                """
                d = {}
                R = self.base_ring()
                for A, coeff in self:
                    la = A.shape()
                    exp = prod(map(factorial, la.to_exp()))
                    if la not in d:
                        if coeff / exp not in R:
                            return False
                        d[la] = [coeff, 1]
                    else:
                        if d[la][0] != coeff:
                            return False
                        d[la][1] += 1
                # Make sure we've seen each set partition of the shape
                return all(d[la][1] == SetPartitions(la.size(), la).cardinality() for la in d)

            def to_symmetric_function(self):
                r"""
                Take a function in the `\mathbf{w}` basis, and return its
                symmetric realization, when possible, expressed in the
                homogeneous basis of symmetric functions.

                OUTPUT:

                - If ``self`` is a symmetric function, then the expansion
                  in the homogeneous basis of the symmetric functions is returned.
                  Otherwise an error is raised.

                EXAMPLES::

                    sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                    sage: elt = w[[1],[2,3]] + w[[1,2],[3]] + w[[1,3],[2]]
                    sage: elt.to_symmetric_function()
                    h[2, 1]
                    sage: elt = w.sum_of_partitions([2,1,1]) / 2
                    sage: elt.to_symmetric_function()
                    1/2*h[2, 1, 1]

                TESTS::

                    sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                    sage: w(0).to_symmetric_function()
                    0
                    sage: w([]).to_symmetric_function()
                    h[]
                    sage: (2*w([])).to_symmetric_function()
                    2*h[]
                """
                if not self.is_symmetric():
                    raise ValueError("not a symmetric function")
                h = SymmetricFunctions(self.parent().base_ring()).homogeneous()
                d = {A.shape(): c for A,c in self}
                return h.sum_of_terms([( AA, cc / prod(map(factorial, AA.to_exp())) )
                                        for AA,cc in d.items()], distinct=True)
