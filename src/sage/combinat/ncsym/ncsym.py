# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
"""
Symmetric Functions in Non-Commuting Variables

AUTHORS:

- Travis Scrimshaw (08-04-2013): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import repeat
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.rings import Rings
from sage.categories.fields import Fields

from sage.arith.misc import factorial
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ncsym.bases import NCSymBases, MultiplicativeNCSymBases, NCSymBasis_abstract
from sage.combinat.set_partition import SetPartitions
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.posets.posets import Poset
from sage.combinat.sf.sf import SymmetricFunctions
from sage.matrix.matrix_space import MatrixSpace
from sage.sets.set import Set
from sage.rings.integer_ring import ZZ
from functools import reduce


def matchings(A, B):
    """
    Iterate through all matchings of the sets `A` and `B`.

    EXAMPLES::

        sage: from sage.combinat.ncsym.ncsym import matchings
        sage: list(matchings([1, 2, 3], [-1, -2]))
        [[[1], [2], [3], [-1], [-2]],
         [[1], [2], [3, -1], [-2]],
         [[1], [2], [3, -2], [-1]],
         [[1], [2, -1], [3], [-2]],
         [[1], [2, -1], [3, -2]],
         [[1], [2, -2], [3], [-1]],
         [[1], [2, -2], [3, -1]],
         [[1, -1], [2], [3], [-2]],
         [[1, -1], [2], [3, -2]],
         [[1, -1], [2, -2], [3]],
         [[1, -2], [2], [3], [-1]],
         [[1, -2], [2], [3, -1]],
         [[1, -2], [2, -1], [3]]]
    """
    lst_A = list(A)
    lst_B = list(B)
    # Handle corner cases
    if not lst_A:
        if not lst_B:
            yield []
        else:
            yield [[b] for b in lst_B]
        return
    if not lst_B:
        yield [[a] for a in lst_A]
        return

    rem_A = lst_A[:]
    a = rem_A.pop(0)
    for m in matchings(rem_A, lst_B):
        yield [[a]] + m
    for i in range(len(lst_B)):
        rem_B = lst_B[:]
        b = rem_B.pop(i)
        for m in matchings(rem_A, rem_B):
            yield [[a, b]] + m


def nesting(la, nu):
    r"""
    Return the nesting number of ``la`` inside of ``nu``.

    If we consider a set partition `A` as a set of arcs `i - j` where `i`
    and `j` are in the same part of `A`. Define

    .. MATH::

        \operatorname{nst}_{\lambda}^{\nu} = \#\{ i < j < k < l \mid
        i - l \in \nu, j - k \in \lambda \},

    and this corresponds to the number of arcs of `\lambda` strictly
    contained inside of `\nu`.

    EXAMPLES::

        sage: from sage.combinat.ncsym.ncsym import nesting
        sage: nu = SetPartition([[1,4], [2], [3]])
        sage: mu = SetPartition([[1,4], [2,3]])
        sage: nesting(set(mu).difference(nu), nu)
        1

        sage: A = SetPartition([[1], [2,5], [3,4]])
        sage: B = SetPartition([[1,3,4], [2,5]])
        sage: nesting(A, B)
        1
        sage: nesting(B, A)
        1

    ::

        sage: lst = list(SetPartitions(4))
        sage: d = {}
        sage: for i, nu in enumerate(lst):
        ....:    for mu in nu.coarsenings():
        ....:        if set(nu.arcs()).issubset(mu.arcs()):
        ....:            d[i, lst.index(mu)] = nesting(set(mu).difference(nu), nu)
        sage: matrix(d)
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    """
    arcs = []
    for p in nu:
        p = sorted(p)
        arcs += [(p[i], p[i+1]) for i in range(len(p)-1)]
    nst = 0
    for p in la:
        p = sorted(p)
        for a in arcs:
            if p[-1] < a[0]:
                continue
            for i in range(len(p)-1):
                if a[1] <= p[i+1]:
                    break
                if a[0] < p[i]:
                    nst += 1
    return nst


class SymmetricFunctionsNonCommutingVariables(UniqueRepresentation, Parent):
    r"""
    Symmetric functions in non-commutative variables.

    The ring of symmetric functions in non-commutative variables,
    which is not to be confused with the :class:`non-commutative symmetric
    functions<NonCommutativeSymmetricFunctions>`, is the ring of all
    bounded-degree noncommutative power series in countably many
    indeterminates (i.e., elements in
    `R \langle \langle x_1, x_2, x_3, \ldots \rangle \rangle` of bounded
    degree) which are invariant with respect to the action of the
    symmetric group `S_{\infty}` on the indices of the indeterminates.
    It can be regarded as a direct limit over all `n \to \infty` of rings
    of `S_n`-invariant polynomials in `n` non-commuting variables
    (that is, `S_n`-invariant elements of `R\langle x_1, x_2, \ldots, x_n \rangle`).

    This ring is implemented as a Hopf algebra whose basis elements are
    indexed by set partitions.

    Let `A = \{A_1, A_2, \ldots, A_r\}` be a set partition of the integers
    `[k] := \{ 1, 2, \ldots, k \}`.  This partition `A` determines an
    equivalence relation `\sim_A` on `[k]`, which has `c \sim_A d` if and
    only if `c` and `d` are in the same part `A_j` of `A`.
    The monomial basis element `\mathbf{m}_A` indexed by `A` is the sum of
    monomials `x_{i_1} x_{i_2} \cdots x_{i_k}` such that `i_c = i_d` if
    and only if `c \sim_A d`.

    The `k`-th graded component of the ring of symmetric functions in
    non-commutative variables has its dimension equal to the number of
    set partitions of `[k]`. (If we work, instead, with finitely many --
    say, `n` -- variables, then its dimension is equal to the number of
    set partitions of `[k]` where the number of parts is at most `n`.)

    .. NOTE::

        All set partitions are considered standard (i.e., set partitions
        of `[n]` for some `n`) unless otherwise stated.

    REFERENCES:

    .. [BZ05] \N. Bergeron, M. Zabrocki. *The Hopf algebra of symmetric
       functions and quasisymmetric functions in non-commutative variables
       are free and cofree*. (2005). :arxiv:`math/0509265v3`.

    .. [BHRZ06] \N. Bergeron, C. Hohlweg, M. Rosas, M. Zabrocki.
       *Grothendieck bialgebras, partition lattices, and symmetric
       functions in noncommutative variables*. Electronic Journal of
       Combinatorics. **13** (2006).

    .. [RS06] \M. Rosas, B. Sagan. *Symmetric functions in noncommuting
       variables*. Trans. Amer. Math. Soc. **358** (2006). no. 1, 215-232.
       :arxiv:`math/0208168`.

    .. [BRRZ08] \N. Bergeron, C. Reutenauer, M. Rosas, M. Zabrocki.
       *Invariants and coinvariants of the symmetric group in noncommuting
       variables*. Canad. J. Math. **60** (2008). 266-296.
       :arxiv:`math/0502082`

    .. [BT13] \N. Bergeron, N. Thiem. *A supercharacter table decomposition
       via power-sum symmetric functions*. Int. J. Algebra Comput. **23**,
       763 (2013). :doi:`10.1142/S0218196713400171`. :arxiv:`1112.4901`.

    EXAMPLES:

    We begin by first creating the ring of `NCSym` and the bases that are
    analogues of the usual symmetric functions::

        sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
        sage: m = NCSym.m()
        sage: e = NCSym.e()
        sage: h = NCSym.h()
        sage: p = NCSym.p()
        sage: m
        Symmetric functions in non-commuting variables over the Rational Field in the monomial basis

    The basis is indexed by set partitions, so we create a few elements and
    convert them between these bases::

        sage: elt = m(SetPartition([[1,3],[2]])) - 2*m(SetPartition([[1],[2]])); elt
        -2*m{{1}, {2}} + m{{1, 3}, {2}}
        sage: e(elt)
        1/2*e{{1}, {2, 3}} - 2*e{{1, 2}} + 1/2*e{{1, 2}, {3}} - 1/2*e{{1, 2, 3}} - 1/2*e{{1, 3}, {2}}
        sage: h(elt)
        -4*h{{1}, {2}} - 2*h{{1}, {2}, {3}} + 1/2*h{{1}, {2, 3}} + 2*h{{1, 2}}
         + 1/2*h{{1, 2}, {3}} - 1/2*h{{1, 2, 3}} + 3/2*h{{1, 3}, {2}}
        sage: p(elt)
        -2*p{{1}, {2}} + 2*p{{1, 2}} - p{{1, 2, 3}} + p{{1, 3}, {2}}
        sage: m(p(elt))
        -2*m{{1}, {2}} + m{{1, 3}, {2}}

        sage: elt = p(SetPartition([[1,3],[2]])) - 4*p(SetPartition([[1],[2]])) + 2; elt
        2*p{} - 4*p{{1}, {2}} + p{{1, 3}, {2}}
        sage: e(elt)
        2*e{} - 4*e{{1}, {2}} + e{{1}, {2}, {3}} - e{{1, 3}, {2}}
        sage: m(elt)
        2*m{} - 4*m{{1}, {2}} - 4*m{{1, 2}} + m{{1, 2, 3}} + m{{1, 3}, {2}}
        sage: h(elt)
        2*h{} - 4*h{{1}, {2}} - h{{1}, {2}, {3}} + h{{1, 3}, {2}}
        sage: p(m(elt))
        2*p{} - 4*p{{1}, {2}} + p{{1, 3}, {2}}

    There is also a shorthand for creating elements. We note that we must use
    ``p[[]]`` to create the empty set partition due to python's syntax. ::

        sage: eltm = m[[1,3],[2]] - 3*m[[1],[2]]; eltm
        -3*m{{1}, {2}} + m{{1, 3}, {2}}
        sage: elte = e[[1,3],[2]]; elte
        e{{1, 3}, {2}}
        sage: elth = h[[1,3],[2,4]]; elth
        h{{1, 3}, {2, 4}}
        sage: eltp = p[[1,3],[2,4]] + 2*p[[1]] - 4*p[[]]; eltp
        -4*p{} + 2*p{{1}} + p{{1, 3}, {2, 4}}

    There is also a natural projection to the usual symmetric functions by
    letting the variables commute.  This projection map preserves the product
    and coproduct structure.  We check that Theorem 2.1 of [RS06]_ holds::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sm = Sym.m()
        sage: Se = Sym.e()
        sage: Sh = Sym.h()
        sage: Sp = Sym.p()
        sage: eltm.to_symmetric_function()
        -6*m[1, 1] + m[2, 1]
        sage: Sm(p(eltm).to_symmetric_function())
        -6*m[1, 1] + m[2, 1]
        sage: elte.to_symmetric_function()
        2*e[2, 1]
        sage: Se(h(elte).to_symmetric_function())
        2*e[2, 1]
        sage: elth.to_symmetric_function()
        4*h[2, 2]
        sage: Sh(m(elth).to_symmetric_function())
        4*h[2, 2]
        sage: eltp.to_symmetric_function()
        -4*p[] + 2*p[1] + p[2, 2]
        sage: Sp(e(eltp).to_symmetric_function())
        -4*p[] + 2*p[1] + p[2, 2]
    """

    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: NCSym1 = SymmetricFunctionsNonCommutingVariables(FiniteField(23))
            sage: NCSym2 = SymmetricFunctionsNonCommutingVariables(Integers(23))
            sage: TestSuite(SymmetricFunctionsNonCommutingVariables(QQ)).run()
        """
        # change the line below to assert R in Rings() once MRO issues from #15536, #15475 are resolved
        assert R in Fields() or R in Rings()  # side effect of this statement assures MRO exists for R
        self._base = R  # Won't be needed once CategoryObject won't override base_ring
        category = GradedHopfAlgebras(R).Cocommutative()
        Parent.__init__(self, category=category.WithRealizations())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctionsNonCommutingVariables(ZZ)
            Symmetric functions in non-commuting variables over the Integer Ring
        """
        return "Symmetric functions in non-commuting variables over the %s" % self.base_ring()

    def a_realization(self):
        r"""
        Return the realization of the powersum basis of ``self``.

        OUTPUT: the powersum basis of symmetric functions in non-commuting variables

        EXAMPLES::

            sage: SymmetricFunctionsNonCommutingVariables(QQ).a_realization()
            Symmetric functions in non-commuting variables over the Rational Field in the powersum basis
        """
        return self.powersum()

    _shorthands = tuple(['chi', 'cp', 'm', 'e', 'h', 'p', 'rho', 'x'])

    def dual(self):
        r"""
        Return the dual Hopf algebra of the symmetric functions in
        non-commuting variables.

        EXAMPLES::

            sage: SymmetricFunctionsNonCommutingVariables(QQ).dual()
            Dual symmetric functions in non-commuting variables over the Rational Field
        """
        from sage.combinat.ncsym.dual import SymmetricFunctionsNonCommutingVariablesDual
        return SymmetricFunctionsNonCommutingVariablesDual(self.base_ring())

    class monomial(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the monomial basis.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: m = NCSym.m()
            sage: m[[1,3],[2]]*m[[1,2]]
            m{{1, 3}, {2}, {4, 5}} + m{{1, 3}, {2, 4, 5}} + m{{1, 3, 4, 5}, {2}}
            sage: m[[1,3],[2]].coproduct()
            m{} # m{{1, 3}, {2}} + m{{1}} # m{{1, 2}} + m{{1, 2}} # m{{1}} + m{{1,
             3}, {2}} # m{}
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.m()).run()
            """
            R = NCSym.base_ring()
            category = GradedHopfAlgebras(R).Cocommutative()
            category &= NCSymBases(NCSym)

            CombinatorialFreeModule.__init__(self, R, SetPartitions(),
                                             prefix='m', bracket=False,
                                             category=category)

        @cached_method
        def _m_to_p_on_basis(self, A):
            r"""
            Return `\mathbf{m}_A` in terms of the powersum basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the powersum basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: all(m(m._m_to_p_on_basis(A)) == m[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            def lt(s, t):
                if s == t:
                    return False
                for p in s:
                    if len([z for z in t if z.intersection(p)]) != 1:
                        return False
                return True

            p = self.realization_of().p()
            P = Poset((A.coarsenings(), lt))
            R = self.base_ring()
            return p._from_dict({B: R(P.moebius_function(A, B)) for B in P})

        @cached_method
        def _m_to_cp_on_basis(self, A):
            r"""
            Return `\mathbf{m}_A` in terms of the `\mathbf{cp}` basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{cp}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: all(m(m._m_to_cp_on_basis(A)) == m[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            cp = self.realization_of().cp()
            arcs = set(A.arcs())
            R = self.base_ring()
            return cp._from_dict({B: R((-1)**len(set(B.arcs()).difference(A.arcs())))
                                  for B in A.coarsenings() if arcs.issubset(B.arcs())},
                                 remove_zeros=False)

        def from_symmetric_function(self, f):
            r"""
            Return the image of the symmetric function ``f`` in ``self``.

            This is performed by converting to the monomial basis and
            extending the method :meth:`sum_of_partitions` linearly.  This is a
            linear map from the symmetric functions to the symmetric functions
            in non-commuting variables that does not preserve the product or
            coproduct structure of the Hopf algebra.

            .. SEEALSO:: :meth:`~Element.to_symmetric_function`

            INPUT:

            - ``f`` -- an element of the symmetric functions

            OUTPUT: an element of the `\mathbf{m}` basis

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: mon = SymmetricFunctions(QQ).m()
                sage: elt = m.from_symmetric_function(mon[2,1,1]); elt
                1/12*m{{1}, {2}, {3, 4}} + 1/12*m{{1}, {2, 3}, {4}} + 1/12*m{{1}, {2, 4}, {3}}
                 + 1/12*m{{1, 2}, {3}, {4}} + 1/12*m{{1, 3}, {2}, {4}} + 1/12*m{{1, 4}, {2}, {3}}
                sage: elt.to_symmetric_function()
                m[2, 1, 1]
                sage: e = SymmetricFunctionsNonCommutingVariables(QQ).e()
                sage: elm = SymmetricFunctions(QQ).e()
                sage: e(m.from_symmetric_function(elm[4]))
                1/24*e{{1, 2, 3, 4}}
                sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                sage: hom = SymmetricFunctions(QQ).h()
                sage: h(m.from_symmetric_function(hom[4]))
                1/24*h{{1, 2, 3, 4}}
                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).p()
                sage: pow = SymmetricFunctions(QQ).p()
                sage: p(m.from_symmetric_function(pow[4]))
                p{{1, 2, 3, 4}}
                sage: p(m.from_symmetric_function(pow[2,1]))
                1/3*p{{1}, {2, 3}} + 1/3*p{{1, 2}, {3}} + 1/3*p{{1, 3}, {2}}
                sage: p([[1,2]])*p([[1]])
                p{{1, 2}, {3}}

            Check that `\chi \circ \widetilde{\chi}` is the identity on `Sym`::

                sage: all(m.from_symmetric_function(pow(la)).to_symmetric_function() == pow(la)
                ....:     for la in Partitions(4))
                True
            """
            m = SymmetricFunctions(self.base_ring()).m()
            return self.sum([c * self.sum_of_partitions(i) for i,c in m(f)])

        def dual_basis(self):
            r"""
            Return the dual basis to the monomial basis.

            OUTPUT: the `\mathbf{w}` basis of the dual Hopf algebra

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m.dual_basis()
                Dual symmetric functions in non-commuting variables over the Rational Field in the w basis
            """
            return self.realization_of().dual().w()

        def duality_pairing(self, x, y):
            r"""
            Compute the pairing between an element of ``self`` and an element
            of the dual.

            INPUT:

            - ``x`` -- an element of symmetric functions in non-commuting
              variables
            - ``y`` -- an element of the dual of symmetric functions in
              non-commuting variables

            OUTPUT: an element of the base ring of ``self``

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: w = m.dual_basis()
                sage: matrix([[m(A).duality_pairing(w(B)) for A in SetPartitions(3)] for B in SetPartitions(3)])
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]
                sage: (m[[1,2],[3]] + 3*m[[1,3],[2]]).duality_pairing(2*w[[1,3],[2]] + w[[1,2,3]] + 2*w[[1,2],[3]])
                8
            """
            x = self(x)
            y = self.dual_basis()(y)
            return sum(coeff * y[I] for (I, coeff) in x)

        def product_on_basis(self, A, B):
            r"""
            The product on monomial basis elements.

            The product of the basis elements indexed by two set partitions `A`
            and `B` is the sum of the basis elements indexed by set partitions
            `C` such that `C \wedge ([n] | [k]) = A | B` where `n = |A|`
            and `k = |B|`. Here `A \wedge B` is the infimum of `A` and `B`
            and `A | B` is the
            :meth:`SetPartition.pipe` operation.
            Equivalently we can describe all `C` as matchings between the
            parts of `A` and `B` where if `a \in A` is matched
            with `b \in B`, we take `a \cup b` instead of `a` and `b` in `C`.

            INPUT:

            - ``A``, ``B`` -- set partitions

            OUTPUT: an element of the `\mathbf{m}` basis

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                sage: A = SetPartition([[1], [2,3]])
                sage: B = SetPartition([[1], [3], [2,4]])
                sage: m.product_on_basis(A, B)
                m{{1}, {2, 3}, {4}, {5, 7}, {6}} + m{{1}, {2, 3, 4}, {5, 7}, {6}}
                 + m{{1}, {2, 3, 5, 7}, {4}, {6}} + m{{1}, {2, 3, 6}, {4}, {5, 7}}
                 + m{{1, 4}, {2, 3}, {5, 7}, {6}} + m{{1, 4}, {2, 3, 5, 7}, {6}}
                 + m{{1, 4}, {2, 3, 6}, {5, 7}} + m{{1, 5, 7}, {2, 3}, {4}, {6}}
                 + m{{1, 5, 7}, {2, 3, 4}, {6}} + m{{1, 5, 7}, {2, 3, 6}, {4}}
                 + m{{1, 6}, {2, 3}, {4}, {5, 7}} + m{{1, 6}, {2, 3, 4}, {5, 7}}
                 + m{{1, 6}, {2, 3, 5, 7}, {4}}
                sage: B = SetPartition([[1], [2]])
                sage: m.product_on_basis(A, B)
                m{{1}, {2, 3}, {4}, {5}} + m{{1}, {2, 3, 4}, {5}}
                 + m{{1}, {2, 3, 5}, {4}} + m{{1, 4}, {2, 3}, {5}} + m{{1, 4}, {2, 3, 5}}
                 + m{{1, 5}, {2, 3}, {4}} + m{{1, 5}, {2, 3, 4}}
                sage: m.product_on_basis(A, SetPartition([]))
                m{{1}, {2, 3}}

            TESTS:

            We check that we get all of the correct set partitions::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                sage: A = SetPartition([[1], [2,3]])
                sage: B = SetPartition([[1], [2]])
                sage: S = SetPartition([[1,2,3], [4,5]])
                sage: AB = SetPartition([[1], [2,3], [4], [5]])
                sage: L = sorted(filter(lambda x: S.inf(x) == AB, SetPartitions(5)), key=str)
                sage: list(map(list, L)) == list(map(list, sorted(m.product_on_basis(A, B).support(), key=str)))
                True
            """
            if not A:
                return self.monomial(B)
            if not B:
                return self.monomial(A)

            P = SetPartitions()
            n = A.size()
            B = [Set([y+n for y in b]) for b in B] # Shift B by n
            unions = lambda m: [reduce(lambda a,b: a.union(b), x) for x in m]
            one = self.base_ring().one()
            return self._from_dict({P(unions(m)): one for m in matchings(A, B)},
                                   remove_zeros=False)

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct of a monomial basis element.

            INPUT:

            - ``A`` -- set partition

            OUTPUT:

            - The coproduct applied to the monomial symmetric function in
              non-commuting variables indexed by ``A`` expressed in the
              monomial basis.

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                sage: m[[1, 3], [2]].coproduct()
                m{} # m{{1, 3}, {2}} + m{{1}} # m{{1, 2}} + m{{1, 2}} # m{{1}} + m{{1, 3}, {2}} # m{}
                sage: m.coproduct_on_basis(SetPartition([]))
                m{} # m{}
                sage: m.coproduct_on_basis(SetPartition([[1,2,3]]))
                m{} # m{{1, 2, 3}} + m{{1, 2, 3}} # m{}
                sage: m[[1,5],[2,4],[3,7],[6]].coproduct()
                m{} # m{{1, 5}, {2, 4}, {3, 7}, {6}} + m{{1}} # m{{1, 5}, {2, 4}, {3, 6}}
                 + 2*m{{1, 2}} # m{{1, 3}, {2, 5}, {4}} + m{{1, 2}} # m{{1, 4}, {2, 3}, {5}}
                 + 2*m{{1, 2}, {3}} # m{{1, 3}, {2, 4}} + m{{1, 3}, {2}} # m{{1, 4}, {2, 3}}
                 + 2*m{{1, 3}, {2, 4}} # m{{1, 2}, {3}} + 2*m{{1, 3}, {2, 5}, {4}} # m{{1, 2}}
                 + m{{1, 4}, {2, 3}} # m{{1, 3}, {2}} + m{{1, 4}, {2, 3}, {5}} # m{{1, 2}}
                 + m{{1, 5}, {2, 4}, {3, 6}} # m{{1}} + m{{1, 5}, {2, 4}, {3, 7}, {6}} # m{}
            """
            P = SetPartitions()
            # Handle corner cases
            if not A:
                return self.tensor_square().monomial(( P([]), P([]) ))
            if len(A) == 1:
                return self.tensor_square().sum_of_monomials([(P([]), A), (A, P([]))])

            ell_set = list(range(1, len(A) + 1))  # +1 for indexing
            L = [[[], ell_set]] + list(SetPartitions(ell_set, 2))

            def to_basis(S):
                if not S:
                    return P([])
                sub_parts = [list(A[i-1]) for i in S] # -1 for indexing
                mins = [min(p) for p in sub_parts]
                over_max = max([max(p) for p in sub_parts]) + 1
                ret = [[] for _ in repeat(None, len(S))]
                cur = 1
                while min(mins) != over_max:
                    m = min(mins)
                    i = mins.index(m)
                    ret[i].append(cur)
                    cur += 1
                    sub_parts[i].pop(sub_parts[i].index(m))
                    if sub_parts[i]:
                        mins[i] = min(sub_parts[i])
                    else:
                        mins[i] = over_max
                return P(ret)
            L1 = [(to_basis(S), to_basis(C)) for S,C in L]
            L2 = [(M, N) for N,M in L1]
            return self.tensor_square().sum_of_monomials(L1 + L2)

        def internal_coproduct_on_basis(self, A):
            r"""
            Return the internal coproduct of a monomial basis element.

            The internal coproduct is defined by

            .. MATH::

                \Delta^{\odot}(\mathbf{m}_A) = \sum_{B \wedge C = A}
                \mathbf{m}_B \otimes \mathbf{m}_C

            where we sum over all pairs of set partitions `B` and `C`
            whose infimum is `A`.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the tensor square of the `\mathbf{m}` basis

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                sage: m.internal_coproduct_on_basis(SetPartition([[1,3],[2]]))
                m{{1, 2, 3}} # m{{1, 3}, {2}} + m{{1, 3}, {2}} # m{{1, 2, 3}} + m{{1, 3}, {2}} # m{{1, 3}, {2}}
            """
            P = SetPartitions()
            SP = SetPartitions(A.size())
            ret = [[A,A]]
            for i, B in enumerate(SP):
                for C in SP[i+1:]:
                    if B.inf(C) == A:
                        B_std = P(list(B.standardization()))
                        C_std = P(list(C.standardization()))
                        ret.append([B_std, C_std])
                        ret.append([C_std, B_std])
            return self.tensor_square().sum_of_monomials((B, C) for B,C in ret)

        def sum_of_partitions(self, la):
            r"""
            Return the sum over all set partitions whose shape is ``la``
            with a fixed coefficient `C` defined below.

            Fix a partition `\lambda`, we define
            `\lambda! := \prod_i \lambda_i!` and `\lambda^! := \prod_i m_i!`.
            Recall that  `|\lambda| = \sum_i \lambda_i` and `m_i` is the
            number of parts of length `i` of `\lambda`. Thus we defined the
            coefficient as

            .. MATH::

                C := \frac{\lambda! \lambda^!}{|\lambda|!}.

            Hence we can define a lift `\widetilde{\chi}` from `Sym`
            to `NCSym` by

            .. MATH::

                m_{\lambda} \mapsto C \sum_A \mathbf{m}_A

            where the sum is over all set partitions whose shape
            is `\lambda`.

            INPUT:

            - ``la`` -- integer partition

            OUTPUT: an element of the `\mathbf{m}` basis

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m.sum_of_partitions(Partition([2,1,1]))
                1/12*m{{1}, {2}, {3, 4}} + 1/12*m{{1}, {2, 3}, {4}} + 1/12*m{{1}, {2, 4}, {3}}
                 + 1/12*m{{1, 2}, {3}, {4}} + 1/12*m{{1, 3}, {2}, {4}} + 1/12*m{{1, 4}, {2}, {3}}

            TESTS:

            Check that `\chi \circ \widetilde{\chi}` is the identity on `Sym`::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: mon = SymmetricFunctions(QQ).monomial()
                sage: all(m.from_symmetric_function(mon[la]).to_symmetric_function() == mon[la]
                ....:     for i in range(6) for la in Partitions(i))
                True
            """
            from sage.combinat.partition import Partition
            la = Partition(la) # Make sure it is a partition
            R = self.base_ring()
            P = SetPartitions()
            c = R( prod(factorial(i) for i in la) / ZZ(factorial(la.size())) )
            return self._from_dict({P(m): c for m in SetPartitions(sum(la), la)},
                                   remove_zeros=False)

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the monomial basis of `NCSym`.
            """

            def expand(self, n, alphabet='x'):
                r"""
                Expand ``self`` written in the monomial basis in `n`
                non-commuting variables.

                INPUT:

                - ``n`` -- integer
                - ``alphabet`` -- string (default: ``'x'``)

                OUTPUT:

                - The symmetric function of ``self`` expressed in the ``n``
                  non-commuting variables described by ``alphabet``.

                EXAMPLES::

                    sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                    sage: m[[1,3],[2]].expand(4)
                    x0*x1*x0 + x0*x2*x0 + x0*x3*x0 + x1*x0*x1 + x1*x2*x1 + x1*x3*x1
                     + x2*x0*x2 + x2*x1*x2 + x2*x3*x2 + x3*x0*x3 + x3*x1*x3 + x3*x2*x3

                One can use a different set of variables by using the
                optional argument ``alphabet``::

                    sage: m[[1],[2,3]].expand(3,alphabet='y')
                    y0*y1^2 + y0*y2^2 + y1*y0^2 + y1*y2^2 + y2*y0^2 + y2*y1^2
                """
                from sage.algebras.free_algebra import FreeAlgebra
                from sage.combinat.permutation import Permutations
                m = self.parent()
                F = FreeAlgebra(m.base_ring(), n, alphabet)

                x = F.gens()

                def on_basis(A):
                    basic_term = [0] * A.size()
                    for index, part in enumerate(A):
                        for i in part:
                            basic_term[i-1] = index # -1 for indexing
                    return sum( prod(x[p[i]-1] for i in basic_term) # -1 for indexing
                                for p in Permutations(n, len(A)) )
                return m._apply_module_morphism(self, on_basis, codomain=F)

            def to_symmetric_function(self):
                r"""
                The projection of ``self`` to the symmetric functions.

                Take a symmetric function in non-commuting variables
                expressed in the `\mathbf{m}` basis, and return the projection of
                expressed in the monomial basis of symmetric functions.

                The map `\chi \colon NCSym \to Sym` is defined by

                .. MATH::

                    \mathbf{m}_A \mapsto
                    m_{\lambda(A)} \prod_i n_i(\lambda(A))!

                where `\lambda(A)` is the partition associated with `A` by
                taking the sizes of the parts and `n_i(\mu)` is the
                multiplicity of `i` in `\mu`.

                OUTPUT: an element of the symmetric functions in the monomial basis

                EXAMPLES::

                    sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                    sage: m[[1,3],[2]].to_symmetric_function()
                    m[2, 1]
                    sage: m[[1],[3],[2]].to_symmetric_function()
                    6*m[1, 1, 1]
                """
                m = SymmetricFunctions(self.parent().base_ring()).monomial()
                c = lambda la: prod(factorial(i) for i in la.to_exp())
                return m.sum_of_terms((i.shape(), coeff*c(i.shape()))
                                      for (i, coeff) in self)

    m = monomial

    class elementary(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the elementary basis.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: e = NCSym.e()
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.e()).run()
            """
            CombinatorialFreeModule.__init__(self, NCSym.base_ring(), SetPartitions(),
                                             prefix='e', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            ## Register coercions
            # monomials
            m = NCSym.m()
            self.module_morphism(self._e_to_m_on_basis, codomain=m).register_as_coercion()
            # powersum
            # NOTE: Keep this ahead of creating the homogeneous basis to
            #   get the coercion path m -> p -> e
            p = NCSym.p()
            self.module_morphism(self._e_to_p_on_basis, codomain=p,
                                 triangular='upper').register_as_coercion()
            p.module_morphism(p._p_to_e_on_basis, codomain=self,
                              triangular='upper').register_as_coercion()
            # homogeneous
            h = NCSym.h()
            self.module_morphism(self._e_to_h_on_basis, codomain=h,
                                 triangular='upper').register_as_coercion()
            h.module_morphism(h._h_to_e_on_basis, codomain=self,
                              triangular='upper').register_as_coercion()

        @cached_method
        def _e_to_m_on_basis(self, A):
            r"""
            Return `\mathbf{e}_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: e = NCSym.e()
                sage: all(e(e._e_to_m_on_basis(A)) == e[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            m = self.realization_of().m()
            n = A.size()
            P = SetPartitions(n)
            min_elt = P([[i] for i in range(1, n+1)])
            one = self.base_ring().one()
            return m._from_dict({B: one for B in P if A.inf(B) == min_elt},
                                remove_zeros=False)

        @cached_method
        def _e_to_h_on_basis(self, A):
            r"""
            Return `\mathbf{e}_A` in terms of the homogeneous basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{h}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: e = NCSym.e()
                sage: all(e(e._e_to_h_on_basis(A)) == e[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            h = self.realization_of().h()
            sign = lambda B: (-1)**(B.size() - len(B))
            coeff = lambda B: sign(B) * prod(factorial(sum( 1 for part in B if part.issubset(big) )) for big in A)
            R = self.base_ring()
            return h._from_dict({B: R(coeff(B)) for B in A.refinements()},
                                remove_zeros=False)

        @cached_method
        def _e_to_p_on_basis(self, A):
            r"""
            Return `\mathbf{e}_A` in terms of the powersum basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{p}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: e = NCSym.e()
                sage: all(e(e._e_to_p_on_basis(A)) == e[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            p = self.realization_of().p()
            coeff = lambda B: prod([(-1)**(i-1) * factorial(i-1) for i in B.shape()])
            R = self.base_ring()
            return p._from_dict({B: R(coeff(B)) for B in A.refinements()},
                                remove_zeros=False)

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the elementary basis of `NCSym`.
            """

            def omega(self):
                r"""
                Return the involution `\omega` applied to ``self``.

                The involution `\omega` on `NCSym` is defined by
                `\omega(\mathbf{e}_A) = \mathbf{h}_A`.

                OUTPUT: an element in the basis ``self``

                EXAMPLES::

                    sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                    sage: e = NCSym.e()
                    sage: h = NCSym.h()
                    sage: elt = e[[1,3],[2]].omega(); elt
                    2*e{{1}, {2}, {3}} - e{{1, 3}, {2}}
                    sage: elt.omega()
                    e{{1, 3}, {2}}
                    sage: h(elt)
                    h{{1, 3}, {2}}
                """
                P = self.parent()
                h = P.realization_of().h()
                return P(h.sum_of_terms(self))

            def to_symmetric_function(self):
                r"""
                The projection of ``self`` to the symmetric functions.

                Take a symmetric function in non-commuting variables
                expressed in the `\mathbf{e}` basis, and return the projection of
                expressed in the elementary basis of symmetric functions.

                The map `\chi \colon NCSym \to Sym` is given by

                .. MATH::

                    \mathbf{e}_A \mapsto
                    e_{\lambda(A)} \prod_i \lambda(A)_i!

                where `\lambda(A)` is the partition associated with `A` by
                taking the sizes of the parts.

                OUTPUT: an element of the symmetric functions in the elementary basis

                EXAMPLES::

                    sage: e = SymmetricFunctionsNonCommutingVariables(QQ).e()
                    sage: e[[1,3],[2]].to_symmetric_function()
                    2*e[2, 1]
                    sage: e[[1],[3],[2]].to_symmetric_function()
                    e[1, 1, 1]
                """
                e = SymmetricFunctions(self.parent().base_ring()).e()
                c = lambda la: prod(factorial(i) for i in la)
                return e.sum_of_terms((i.shape(), coeff*c(i.shape()))
                                      for (i, coeff) in self)

    e = elementary

    class homogeneous(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the homogeneous basis.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: h = NCSym.h()
            sage: h[[1,3],[2,4]]*h[[1,2,3]]
            h{{1, 3}, {2, 4}, {5, 6, 7}}
            sage: h[[1,2]].coproduct()
            h{} # h{{1, 2}} + 2*h{{1}} # h{{1}} + h{{1, 2}} # h{}
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.h()).run()
            """
            CombinatorialFreeModule.__init__(self, NCSym.base_ring(), SetPartitions(),
                                             prefix='h', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            # Register coercions
            m = NCSym.m()
            self.module_morphism(self._h_to_m_on_basis, codomain=m).register_as_coercion()
            p = NCSym.p()
            self.module_morphism(self._h_to_p_on_basis, codomain=p).register_as_coercion()
            p.module_morphism(p._p_to_h_on_basis, codomain=self).register_as_coercion()

        @cached_method
        def _h_to_m_on_basis(self, A):
            r"""
            Return `\mathbf{h}_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: h = NCSym.h()
                sage: all(h(h._h_to_m_on_basis(A)) == h[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            P = SetPartitions()
            m = self.realization_of().m()
            coeff = lambda B: prod(factorial(i) for i in B.shape())
            R = self.base_ring()
            return m._from_dict({P(B): R( coeff(A.inf(B)) )
                                 for B in SetPartitions(A.size())}, remove_zeros=False)

        @cached_method
        def _h_to_e_on_basis(self, A):
            r"""
            Return `\mathbf{h}_A` in terms of the elementary basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{e}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: h = NCSym.h()
                sage: all(h(h._h_to_e_on_basis(A)) == h[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            e = self.realization_of().e()
            sign = lambda B: (-1)**(B.size() - len(B))
            coeff = lambda B: (sign(B) * prod(factorial(sum( 1 for part in B if part.issubset(big) ))
                                              for big in A))
            R = self.base_ring()
            return e._from_dict({B: R(coeff(B)) for B in A.refinements()},
                                remove_zeros=False)

        @cached_method
        def _h_to_p_on_basis(self, A):
            r"""
            Return `\mathbf{h}_A` in terms of the powersum basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{p}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: h = NCSym.h()
                sage: all(h(h._h_to_p_on_basis(A)) == h[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            p = self.realization_of().p()
            coeff = lambda B: abs( prod([(-1)**(i-1) * factorial(i-1) for i in B.shape()]) )
            R = self.base_ring()
            return p._from_dict({B: R(coeff(B)) for B in A.refinements()},
                                remove_zeros=False)

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the homogeneous basis of `NCSym`.
            """

            def omega(self):
                r"""
                Return the involution `\omega` applied to ``self``.

                The involution `\omega` on `NCSym` is defined by
                `\omega(\mathbf{h}_A) = \mathbf{e}_A`.

                OUTPUT: an element in the basis ``self``

                EXAMPLES::

                    sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                    sage: h = NCSym.h()
                    sage: e = NCSym.e()
                    sage: elt = h[[1,3],[2]].omega(); elt
                    2*h{{1}, {2}, {3}} - h{{1, 3}, {2}}
                    sage: elt.omega()
                    h{{1, 3}, {2}}
                    sage: e(elt)
                    e{{1, 3}, {2}}
                """
                P = self.parent()
                e = self.parent().realization_of().e()
                return P(e.sum_of_terms(self))

            def to_symmetric_function(self):
                r"""
                The projection of ``self`` to the symmetric functions.

                Take a symmetric function in non-commuting variables
                expressed in the `\mathbf{h}` basis, and return the projection of
                expressed in the complete basis of symmetric functions.

                The map `\chi \colon NCSym \to Sym` is given by

                .. MATH::

                    \mathbf{h}_A \mapsto
                    h_{\lambda(A)} \prod_i \lambda(A)_i!

                where `\lambda(A)` is the partition associated with `A` by
                taking the sizes of the parts.

                OUTPUT: an element of the symmetric functions in the complete basis

                EXAMPLES::

                    sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                    sage: h[[1,3],[2]].to_symmetric_function()
                    2*h[2, 1]
                    sage: h[[1],[3],[2]].to_symmetric_function()
                    h[1, 1, 1]
                """
                h = SymmetricFunctions(self.parent().base_ring()).h()
                c = lambda la: prod(factorial(i) for i in la)
                return h.sum_of_terms((i.shape(), coeff*c(i.shape()))
                                      for (i, coeff) in self)

    h = homogeneous

    class powersum(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the powersum basis.

        The powersum basis is given by

        .. MATH::

            \mathbf{p}_A = \sum_{A \leq B} \mathbf{m}_B,

        where we sum over all coarsenings of the set partition `A`. If we
        allow our variables to commute, then `\mathbf{p}_A` goes to the
        usual powersum symmetric function `p_{\lambda}` whose (integer)
        partition `\lambda` is the shape of `A`.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: p = NCSym.p()

            sage: x = p.an_element()**2; x
            4*p{} + 8*p{{1}} + 4*p{{1}, {2}} + 6*p{{1}, {2, 3}}
             + 12*p{{1, 2}} + 6*p{{1, 2}, {3}} + 9*p{{1, 2}, {3, 4}}
            sage: x.to_symmetric_function()
            4*p[] + 8*p[1] + 4*p[1, 1] + 12*p[2] + 12*p[2, 1] + 9*p[2, 2]
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.p()).run()
            """
            CombinatorialFreeModule.__init__(self, NCSym.base_ring(), SetPartitions(),
                                             prefix='p', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            # Register coercions
            m = NCSym.m()
            self.module_morphism(self._p_to_m_on_basis, codomain=m,
                                 unitriangular='lower').register_as_coercion()
            m.module_morphism(m._m_to_p_on_basis, codomain=self,
                              unitriangular='lower').register_as_coercion()
            x = NCSym.x()
            self.module_morphism(self._p_to_x_on_basis, codomain=x,
                                 unitriangular='upper').register_as_coercion()
            x.module_morphism(x._x_to_p_on_basis, codomain=self,
                              unitriangular='upper').register_as_coercion()

        @cached_method
        def _p_to_m_on_basis(self, A):
            r"""
            Return `\mathbf{p}_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: p = NCSym.p()
                sage: all(p(p._p_to_m_on_basis(A)) == p[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            m = self.realization_of().m()
            one = self.base_ring().one()
            return m._from_dict({B: one for B in A.coarsenings()}, remove_zeros=False)

        @cached_method
        def _p_to_e_on_basis(self, A):
            r"""
            Return `\mathbf{p}_A` in terms of the elementary basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{e}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: p = NCSym.p()
                sage: all(p(p._p_to_e_on_basis(A)) == p[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            e = self.realization_of().e()
            P_refine = Poset((A.refinements(), A.parent().lt))
            c = prod((-1)**(i-1) * factorial(i-1) for i in A.shape())
            R = self.base_ring()
            return e._from_dict({B: R(P_refine.moebius_function(B, A) / ZZ(c))
                                 for B in P_refine}, remove_zeros=False)

        @cached_method
        def _p_to_h_on_basis(self, A):
            r"""
            Return `\mathbf{p}_A` in terms of the homogeneous basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{h}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: p = NCSym.p()
                sage: all(p(p._p_to_h_on_basis(A)) == p[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            h = self.realization_of().h()
            P_refine = Poset((A.refinements(), A.parent().lt))
            c = abs(prod((-1)**(i-1) * factorial(i-1) for i in A.shape()))
            R = self.base_ring()
            return h._from_dict({B: R(P_refine.moebius_function(B, A) / ZZ(c))
                                 for B in P_refine}, remove_zeros=False)

        @cached_method
        def _p_to_x_on_basis(self, A):
            r"""
            Return `\mathbf{p}_A` in terms of the `\mathbf{x}` basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{x}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: p = NCSym.p()
                sage: all(p(p._p_to_x_on_basis(A)) == p[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            x = self.realization_of().x()
            one = self.base_ring().one()
            return x._from_dict({B: one for B in A.refinements()}, remove_zeros=False)

        # Note that this is the same as the monomial coproduct_on_basis
        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct of a monomial basis element.

            INPUT:

            - ``A`` -- set partition

            OUTPUT:

            - The coproduct applied to the monomial symmetric function in
              non-commuting variables indexed by ``A`` expressed in the
              monomial basis.

            EXAMPLES::

                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).powersum()
                sage: p[[1, 3], [2]].coproduct()
                p{} # p{{1, 3}, {2}} + p{{1}} # p{{1, 2}} + p{{1, 2}} # p{{1}} + p{{1, 3}, {2}} # p{}
                sage: p.coproduct_on_basis(SetPartition([[1]]))
                p{} # p{{1}} + p{{1}} # p{}
                sage: p.coproduct_on_basis(SetPartition([]))
                p{} # p{}
            """
            P = SetPartitions()
            # Handle corner cases
            if not A:
                return self.tensor_square().monomial(( P([]), P([]) ))
            if len(A) == 1:
                return self.tensor_square().sum_of_monomials([(P([]), A), (A, P([]))])

            ell_set = list(range(1, len(A) + 1))  # +1 for indexing
            L = [[[], ell_set]] + list(SetPartitions(ell_set, 2))

            def to_basis(S):
                if not S:
                    return P([])
                sub_parts = [list(A[i-1]) for i in S] # -1 for indexing
                mins = [min(p) for p in sub_parts]
                over_max = max([max(p) for p in sub_parts]) + 1
                ret = [[] for _ in repeat(None, len(S))]
                cur = 1
                while min(mins) != over_max:
                    m = min(mins)
                    i = mins.index(m)
                    ret[i].append(cur)
                    cur += 1
                    sub_parts[i].pop(sub_parts[i].index(m))
                    if sub_parts[i]:
                        mins[i] = min(sub_parts[i])
                    else:
                        mins[i] = over_max
                return P(ret)
            L1 = [(to_basis(S), to_basis(C)) for S,C in L]
            L2 = [(M, N) for N,M in L1]
            return self.tensor_square().sum_of_monomials(L1 + L2)

        def internal_coproduct_on_basis(self, A):
            r"""
            Return the internal coproduct of a powersum basis element.

            The internal coproduct is defined by

            .. MATH::

                \Delta^{\odot}(\mathbf{p}_A) = \mathbf{p}_A \otimes
                \mathbf{p}_A

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the tensor square of ``self``

            EXAMPLES::

                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).powersum()
                sage: p.internal_coproduct_on_basis(SetPartition([[1,3],[2]]))
                p{{1, 3}, {2}} # p{{1, 3}, {2}}
            """
            return self.tensor_square().monomial((A, A))

        def antipode_on_basis(self, A):
            r"""
            Return the result of the antipode applied to a powersum basis element.

            Let `A` be a set partition. The antipode given in [LM2011]_ is

            .. MATH::

                S(\mathbf{p}_A) = \sum_{\gamma} (-1)^{\ell(\gamma)}
                \mathbf{p}_{\gamma[A]}

            where we sum over all ordered set partitions (i.e. set
            compositions) of `[\ell(A)]` and

            .. MATH::

                \gamma[A] = A_{\gamma_1}^{\downarrow} | \cdots |
                A_{\gamma_{\ell(A)}}^{\downarrow}

            is the action of `\gamma` on `A` defined in
            :meth:`SetPartition.ordered_set_partition_action()`.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element in the basis ``self``

            EXAMPLES::

                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).powersum()
                sage: p.antipode_on_basis(SetPartition([[1], [2,3]]))
                p{{1, 2}, {3}}
                sage: p.antipode_on_basis(SetPartition([]))
                p{}
                sage: F = p[[1,3],[5],[2,4]].coproduct()
                sage: F.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            P = SetPartitions()

            def action(gamma):
                cur = 1
                ret = []
                for S in gamma:
                    sub_parts = [list(A[i - 1]) for i in S]  # -1 for indexing
                    mins = [min(p) for p in sub_parts]
                    over_max = max([max(p) for p in sub_parts]) + 1
                    temp = [[] for _ in repeat(None, len(S))]
                    while min(mins) != over_max:
                        m = min(mins)
                        i = mins.index(m)
                        temp[i].append(cur)
                        cur += 1
                        sub_parts[i].pop(sub_parts[i].index(m))
                        if sub_parts[i]:
                            mins[i] = min(sub_parts[i])
                        else:
                            mins[i] = over_max
                    ret += temp
                return P(ret)
            return self.sum_of_terms( (A.ordered_set_partition_action(gamma), (-1)**len(gamma))
                                      for gamma in OrderedSetPartitions(len(A)) )

        def primitive(self, A, i=1):
            r"""
            Return the primitive associated to ``A`` in ``self``.

            Fix some `i \in S`. Let `A` be an atomic set partition of `S`,
            then the primitive `p(A)` given in [LM2011]_ is

            .. MATH::

                p(A) = \sum_{\gamma} (-1)^{\ell(\gamma)-1}
                \mathbf{p}_{\gamma[A]}

            where we sum over all ordered set partitions of `[\ell(A)]` such
            that `i \in \gamma_1` and `\gamma[A]` is the action of `\gamma`
            on `A` defined in
            :meth:`SetPartition.ordered_set_partition_action()`.
            If `A` is not atomic, then `p(A) = 0`.

            .. SEEALSO:: :meth:`SetPartition.is_atomic`

            INPUT:

            - ``A`` -- set partition
            - ``i`` -- (default: 1) index in the base set for ``A`` specifying
              which set of primitives this belongs to

            OUTPUT: an element in the basis ``self``

            EXAMPLES::

                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).powersum()
                sage: elt = p.primitive(SetPartition([[1,3], [2]])); elt
                -p{{1, 2}, {3}} + p{{1, 3}, {2}}
                sage: elt.coproduct()
                -p{} # p{{1, 2}, {3}} + p{} # p{{1, 3}, {2}} - p{{1, 2}, {3}} # p{} + p{{1, 3}, {2}} # p{}
                sage: p.primitive(SetPartition([[1], [2,3]]))
                0
                sage: p.primitive(SetPartition([]))
                p{}
            """
            if not A:
                return self.one()
            A = SetPartitions()(A) # Make sure it's a set partition
            if not A.is_atomic():
                return self.zero()
            return self.sum_of_terms( (A.ordered_set_partition_action(gamma), (-1)**(len(gamma)-1))
                                      for gamma in OrderedSetPartitions(len(A)) if i in gamma[0] )

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the powersum basis of `NCSym`.
            """

            def to_symmetric_function(self):
                r"""
                The projection of ``self`` to the symmetric functions.

                Take a symmetric function in non-commuting variables
                expressed in the `\mathbf{p}` basis, and return the projection of
                expressed in the powersum basis of symmetric functions.

                The map `\chi \colon NCSym \to Sym` is given by

                .. MATH::

                    \mathbf{p}_A \mapsto p_{\lambda(A)}

                where `\lambda(A)` is the partition associated with `A` by
                taking the sizes of the parts.

                OUTPUT: an element of symmetric functions in the power sum basis

                EXAMPLES::

                    sage: p = SymmetricFunctionsNonCommutingVariables(QQ).p()
                    sage: p[[1,3],[2]].to_symmetric_function()
                    p[2, 1]
                    sage: p[[1],[3],[2]].to_symmetric_function()
                    p[1, 1, 1]
                """
                p = SymmetricFunctions(self.parent().base_ring()).p()
                return p.sum_of_terms((i.shape(), coeff) for (i, coeff) in self)

    p = powersum

    class coarse_powersum(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the `\mathbf{cp}` basis.

        This basis was defined in [BZ05]_ as

        .. MATH::

            \mathbf{cp}_A = \sum_{A \leq_* B} \mathbf{m}_B,

        where we sum over all strict coarsenings of the set partition `A`.
        An alternative description of this basis was given in [BT13]_ as

        .. MATH::

            \mathbf{cp}_A = \sum_{A \subseteq B} \mathbf{m}_B,

        where we sum over all set partitions whose arcs are a subset of
        the arcs of the set partition `A`.

        .. NOTE::

            In [BZ05]_, this basis was denoted by `\mathbf{q}`. In [BT13]_,
            this basis was called the powersum basis and denoted by `p`.
            However it is a coarser basis than the usual powersum basis in
            the sense that it does not yield the usual powersum basis
            of the symmetric function under the natural map of letting
            the variables commute.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: cp = NCSym.cp()
            sage: cp[[1,3],[2,4]]*cp[[1,2,3]]
            cp{{1, 3}, {2, 4}, {5, 6, 7}}
            sage: cp[[1,2],[3]].internal_coproduct()
            cp{{1, 2}, {3}} # cp{{1, 2}, {3}}
            sage: ps = SymmetricFunctions(NCSym.base_ring()).p()
            sage: ps(cp[[1,3],[2]].to_symmetric_function())
            p[2, 1] - p[3]
            sage: ps(cp[[1,2],[3]].to_symmetric_function())
            p[2, 1]
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.cp()).run()
            """
            CombinatorialFreeModule.__init__(self, NCSym.base_ring(), SetPartitions(),
                                             prefix='cp', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            # Register coercions
            m = NCSym.m()
            self.module_morphism(self._cp_to_m_on_basis, codomain=m,
                                 unitriangular='lower').register_as_coercion()
            m.module_morphism(m._m_to_cp_on_basis, codomain=self,
                              unitriangular='lower').register_as_coercion()

        @cached_method
        def _cp_to_m_on_basis(self, A):
            r"""
            Return `\mathbf{cp}_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: cp = NCSym.cp()
                sage: all(cp(cp._cp_to_m_on_basis(A)) == cp[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            m = self.realization_of().m()
            one = self.base_ring().one()
            return m._from_dict({B: one for B in A.strict_coarsenings()},
                                remove_zeros=False)

    cp = coarse_powersum

    class x_basis(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the `\mathbf{x}` basis.

        This basis is defined in [BHRZ06]_ by the formula:

        .. MATH::

            \mathbf{x}_A = \sum_{B \leq A} \mu(B, A) \mathbf{p}_B

        and has the following properties:

        .. MATH::

            \mathbf{x}_A \mathbf{x}_B = \mathbf{x}_{A|B}, \quad \quad
            \Delta^{\odot}(\mathbf{x}_C) = \sum_{A \vee B = C} \mathbf{x}_A
            \otimes \mathbf{x}_B.

        EXAMPLES::

            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: x = NCSym.x()
            sage: x[[1,3],[2,4]]*x[[1,2,3]]
            x{{1, 3}, {2, 4}, {5, 6, 7}}
            sage: x[[1,2],[3]].internal_coproduct()
            x{{1}, {2}, {3}} # x{{1, 2}, {3}} + x{{1, 2}, {3}} # x{{1}, {2}, {3}} +
             x{{1, 2}, {3}} # x{{1, 2}, {3}}
        """

        def __init__(self, NCSym):
            """
            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: TestSuite(NCSym.x()).run()
            """
            CombinatorialFreeModule.__init__(self, NCSym.base_ring(), SetPartitions(),
                                             prefix='x', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))

        @cached_method
        def _x_to_p_on_basis(self, A):
            r"""
            Return `\mathbf{x}_A` in terms of the powersum basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{p}` basis

            TESTS::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: x = NCSym.x()
                sage: all(x(x._x_to_p_on_basis(A)) == x[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            def lt(s, t):
                if s == t:
                    return False
                for p in s:
                    if len([z for z in t if z.intersection(p)]) != 1:
                        return False
                return True

            p = self.realization_of().p()
            P_refine = Poset((A.refinements(), lt))
            R = self.base_ring()
            return p._from_dict({B: R(P_refine.moebius_function(B, A))
                                 for B in P_refine})

    x = x_basis

    class deformed_coarse_powersum(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the `\rho` basis.

        This basis was defined in [BT13]_ as a `q`-deformation of the
        `\mathbf{cp}` basis:

        .. MATH::

            \rho_A = \sum_{A \subseteq B}
            \frac{1}{q^{\operatorname{nst}_{B-A}^A}} \mathbf{m}_B,

        where we sum over all set partitions whose arcs are a subset of
        the arcs of the set partition `A`.

        INPUT:

        - ``q`` -- (default: ``2``) the parameter `q`

        EXAMPLES::

            sage: R = QQ['q'].fraction_field()
            sage: q = R.gen()
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
            sage: rho = NCSym.rho(q)

        We construct Example 3.1 in [BT13]_::

            sage: rnode = lambda A: sorted([a[1] for a in A.arcs()], reverse=True)
            sage: dimv = lambda A: sorted([a[1]-a[0] for a in A.arcs()], reverse=True)
            sage: lst = list(SetPartitions(4))
            sage: S = sorted(lst, key=lambda A: (dimv(A), rnode(A)))
            sage: m = NCSym.m()
            sage: matrix([[m(rho[A])[B] for B in S] for A in S])
            [  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]
            [  0   1   0   0   1   1   0   1   0   0   1   0   0   0   0]
            [  0   0   1   0   1   0   1   1   0   0   0   0   0   0   1]
            [  0   0   0   1   0   1   1   1   0   0   0   1   0   0   0]
            [  0   0   0   0   1   0   0   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   1   0   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   1   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   1   0   0   1   1   0   0]
            [  0   0   0   0   0   0   0   0   0   1   1   0   1   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1 1/q]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
        """

        def __init__(self, NCSym, q=2):
            """
            EXAMPLES::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: TestSuite(NCSym.rho(q)).run()
            """
            R = NCSym.base_ring()
            self._q = R(q)
            CombinatorialFreeModule.__init__(self, R, SetPartitions(),
                                             prefix='rho', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            # Register coercions
            m = NCSym.m()
            self.module_morphism(self._rho_to_m_on_basis, codomain=m).register_as_coercion()
            m.module_morphism(self._m_to_rho_on_basis, codomain=self).register_as_coercion()

        def q(self):
            """
            Return the deformation parameter `q` of ``self``.

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: rho = NCSym.rho(5)
                sage: rho.q()
                5

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: rho = NCSym.rho(q)
                sage: rho.q() == q
                True
            """
            return self._q

        @cached_method
        def _rho_to_m_on_basis(self, A):
            r"""
            Return `\rho_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: rho = NCSym.rho(q)
                sage: all(rho(rho._rho_to_m_on_basis(A)) == rho[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            m = self.realization_of().m()
            arcs = set(A.arcs())
            return m._from_dict({B: self._q**-nesting(set(B).difference(A), A)
                                 for B in A.coarsenings() if arcs.issubset(B.arcs())},
                                remove_zeros=False)

        @cached_method
        def _m_to_rho_on_basis(self, A):
            r"""
            Return `\mathbf{m}_A` in terms of the `\rho` basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\rho` basis

            TESTS::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: rho = NCSym.rho(q)
                sage: m = NCSym.m()
                sage: all(m(rho._m_to_rho_on_basis(A)) == m[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            coeff = lambda A,B: ((-1)**len(set(B.arcs()).difference(A.arcs()))
                                 / self._q**nesting(set(B).difference(A), B))
            arcs = set(A.arcs())
            return self._from_dict({B: coeff(A,B) for B in A.coarsenings()
                                    if arcs.issubset(B.arcs())},
                                   remove_zeros=False)

    rho = deformed_coarse_powersum

    class supercharacter(NCSymBasis_abstract):
        r"""
        The Hopf algebra of symmetric functions in non-commuting variables
        in the supercharacter `\chi` basis.

        This basis was defined in [BT13]_ as a `q`-deformation of the
        supercharacter basis.

        .. MATH::

            \chi_A = \sum_B \chi_A(B) \mathbf{m}_B,

        where we sum over all set partitions `A` and `\chi_A(B)` is the
        evaluation of the supercharacter `\chi_A` on the superclass `\mu_B`.

        .. NOTE::

            The supercharacters considered in [BT13]_ are coarser than
            those considered by Aguiar et. al.

        INPUT:

        - ``q`` -- (default: ``2``) the parameter `q`

        EXAMPLES::

            sage: R = QQ['q'].fraction_field()
            sage: q = R.gen()
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
            sage: chi = NCSym.chi(q)
            sage: chi[[1,3],[2]]*chi[[1,2]]
            chi{{1, 3}, {2}, {4, 5}}
            sage: chi[[1,3],[2]].coproduct()
            chi{} # chi{{1, 3}, {2}} + (2*q-2)*chi{{1}} # chi{{1}, {2}} +
             (3*q-2)*chi{{1}} # chi{{1, 2}} + (2*q-2)*chi{{1}, {2}} # chi{{1}} +
             (3*q-2)*chi{{1, 2}} # chi{{1}} + chi{{1, 3}, {2}} # chi{}
            sage: chi2 = NCSym.chi()
            sage: chi(chi2[[1,2],[3]])
            ((-q+2)/q)*chi{{1}, {2}, {3}} + 2/q*chi{{1, 2}, {3}}
            sage: chi2
            Symmetric functions in non-commuting variables over the Fraction Field
             of Univariate Polynomial Ring in q over Rational Field in the
             supercharacter basis with parameter q=2
        """

        def __init__(self, NCSym, q=2):
            """
            EXAMPLES::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: TestSuite(NCSym.chi(q)).run()
            """
            R = NCSym.base_ring()
            self._q = R(q)
            CombinatorialFreeModule.__init__(self, R, SetPartitions(),
                                             prefix='chi', bracket=False,
                                             category=MultiplicativeNCSymBases(NCSym))
            # Register coercions
            m = NCSym.m()
            self.module_morphism(self._chi_to_m_on_basis, codomain=m).register_as_coercion()
            m.module_morphism(self._m_to_chi_on_basis, codomain=self).register_as_coercion()

        def q(self):
            """
            Return the deformation parameter `q` of ``self``.

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: chi = NCSym.chi(5)
                sage: chi.q()
                5

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: chi = NCSym.chi(q)
                sage: chi.q() == q
                True
            """
            return self._q

        @cached_method
        def _chi_to_m_on_basis(self, A):
            r"""
            Return `\chi_A` in terms of the monomial basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\mathbf{m}` basis

            TESTS::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: chi = NCSym.chi(q)
                sage: all(chi(chi._chi_to_m_on_basis(A)) == chi[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            m = self.realization_of().m()
            q = self._q
            arcs = set(A.arcs())
            ret = {}
            for B in SetPartitions(A.size()):
                Barcs = B.arcs()
                if any((a[0] == b[0] and b[1] < a[1])
                       or (b[0] > a[0] and a[1] == b[1])
                       for a in arcs for b in Barcs):
                    continue
                ret[B] = ((-1)**len(arcs.intersection(Barcs))
                          * (q - 1)**(len(arcs) - len(arcs.intersection(Barcs)))
                          * q**(sum(a[1] - a[0] for a in arcs) - len(arcs))
                          / q**nesting(B, A))
            return m._from_dict(ret, remove_zeros=False)

        @cached_method
        def _graded_inverse_matrix(self, n):
            r"""
            Return the inverse of the transition matrix of the ``n``-th
            graded part from the `\chi` basis to the monomial basis.

            EXAMPLES::

                sage: R = QQ['q'].fraction_field(); q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: chi = NCSym.chi(q); m = NCSym.m()
                sage: lst = list(SetPartitions(2))
                sage: m = matrix([[m(chi[A])[B] for A in lst] for B in lst]); m
                [   -1     1]
                [q - 1     1]
                sage: chi._graded_inverse_matrix(2)
                [     -1/q       1/q]
                [(q - 1)/q       1/q]
                sage: chi._graded_inverse_matrix(2) * m
                [1 0]
                [0 1]
            """
            lst = SetPartitions(n)
            MS = MatrixSpace(self.base_ring(), lst.cardinality())
            m = self.realization_of().m()
            m = MS([[m(self[A])[B] for A in lst] for B in lst])
            return ~m

        @cached_method
        def _m_to_chi_on_basis(self, A):
            r"""
            Return `\mathbf{m}_A` in terms of the `\chi` basis.

            INPUT:

            - ``A`` -- set partition

            OUTPUT: an element of the `\chi` basis

            TESTS::

                sage: R = QQ['q'].fraction_field()
                sage: q = R.gen()
                sage: NCSym = SymmetricFunctionsNonCommutingVariables(R)
                sage: chi = NCSym.chi(q)
                sage: m = NCSym.m()
                sage: all(m(chi._m_to_chi_on_basis(A)) == m[A] for i in range(5)
                ....:     for A in SetPartitions(i))
                True
            """
            n = A.size()
            lst = list(SetPartitions(n))
            m = self._graded_inverse_matrix(n)
            i = lst.index(A)
            return self._from_dict({B: m[j,i] for j,B in enumerate(lst)})

    chi = supercharacter
