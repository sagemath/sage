# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Partition Shifting Algebras

This module contains families of operators that act on partitions or,
more generally, integer sequences. In particular, this includes Young's
raising operators `R_{ij}`, which act on integer sequences by adding `1`
to the `i`-th entry and subtracting `1` to the `j`-th entry. A special
case is acting on partitions.

AUTHORS:

- Matthew Lancellotti, George H. Seelinger (2018): Initial version
"""
# ****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.sets_cat import Sets
from sage.categories.algebras import Algebras
from sage.structure.parent import Parent
from sage.combinat.composition import Composition
from sage.combinat.partition import _Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.fast_methods import Singleton
from sage.misc.cachefunc import cached_method
from sage.rings.rational_field import QQ
from sage.rings.semirings.non_negative_integer_semiring import NonNegativeIntegerSemiring
from sage.rings.integer_ring import ZZ


class ShiftingSequenceSpace(Singleton, Parent):
    r"""
    A helper for :class:`ShiftingOperatorAlgebra` that contains all
    tuples with entries in `\ZZ` of finite support with no trailing `0`'s.

    EXAMPLES::

        sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
        sage: S = ShiftingSequenceSpace()
        sage: (1, -1) in S
        True
        sage: (1, -1, 0, 9) in S
        True
        sage: [1, -1] in S
        False
        sage: (0.5, 1) in S
        False
    """

    def __init__(self):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
            sage: S = ShiftingSequenceSpace()
        """
        Parent.__init__(self, facade=(tuple,),
                        category=Sets().Infinite().Facade())

    def __contains__(self, seq):
        r"""
        Return ``True`` if and only if ``seq`` is a valid shifting sequence.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
            sage: S = ShiftingSequenceSpace()
            sage: (1, -1) in S
            True
            sage: (1, -1, 0, 9) in S
            True
            sage: (1, -1, 0) in S
            False
            sage: () in S
            True
            sage: [1, -1] in S
            False
            sage: (0.5, 1) in S
            False
        """
        return (isinstance(seq, tuple) and all(i in ZZ for i in seq)
                and (not seq or seq[-1]))

    def check(self, seq):
        r"""
        Verify that ``seq`` is a valid shifting sequence.

        If it is not, raise a :exc:`ValueError`.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
            sage: S = ShiftingSequenceSpace()
            sage: S.check((1, -1))
            sage: S.check((1, -1, 0, 9))
            sage: S.check([1, -1])
            Traceback (most recent call last):
            ...
            ValueError: invalid index [1, -1]
            sage: S.check((0.5, 1))
            Traceback (most recent call last):
            ...
            ValueError: invalid index (0.500000000000000, 1)
        """
        if not self.__contains__(seq):
            raise ValueError('invalid index {}'.format(seq))


class ShiftingOperatorAlgebra(CombinatorialFreeModule):
    r"""
    An algebra of shifting operators.

    Let `R` be a commutative ring. The algebra of shifting operators
    is isomorphic as an `R`-algebra to the Laurent polynomial ring
    `R[x_1^\pm, x_2^\pm, x_3^\pm, \ldots]`. Moreover, the monomials of
    the shifting operator algebra act on any integer sequence `\lambda
    = (\lambda_1, \lambda_2, \ldots, \lambda_{\ell})` as follows. Let
    `S` be our algebra of shifting operators. Then, for any monomial
    `s = x_1^{a_1}x_2^{a_2} \cdots x_r^{a_r} \in S` where `a_i \in
    \ZZ` and `r \geq \ell`, we get that `s.\lambda = (\lambda_1
    + a_1, \lambda_2 + a_2,\ldots,\lambda_r+a_r)` where we pad
    `\lambda` with `r-\ell` zeros. In particular, we can recover
    Young's raising operator, `R_{ij}`, for `i < j`, acting on
    partitions by having `\frac{x_i}{x_j}` act on a partition
    `\lambda`.

    One can extend the action of these shifting operators to a basis
    of symmetric functions, but at the expense of no longer actually
    having a well-defined operator. Formally, to extend the action of
    the shifting operators on a symmetric function basis `B =
    \{b_{\lambda}\}_{\lambda}`, we define an `R`-module homomorphism
    `\phi : R[x_1^\pm, x_2^\pm, \ldots] \to B`. Then we compute
    `x_1^{a_1} \cdots x_r^{a_r}.b_\lambda` by first computing
    `(x_1^{a_1} \cdots x_r^{a_r})x_1^{\lambda_1} \cdots
    x_\ell^{\lambda_\ell}` and then applying `\phi` to the result. For
    examples of what this looks like with specific bases, see below.

    This implementation is consistent with how many references work
    formally with raising operators. For instance, see exposition
    surrounding [BMPS2018]_ Equation (4.1).

    We follow the following convention for creating elements: ``S(1,
    0, -1, 2)`` is the shifting operator that raises the first part by
    `1`, lowers the third part by `1`, and raises the fourth part by
    `2`.

    In addition to acting on partitions (or any integer sequence), the
    shifting operators can also act on symmetric functions in a basis
    `B` when a conversion to `B` has been registered, preferably using
    :meth:`build_and_register_conversion`.

    For a definition of raising operators, see [BMPS2018]_ Definition
    2.1. See :meth:`ij` to create operators using the notation in
    [BMPS2018]_.

    INPUT:

    - ``base_ring`` -- (default: ``QQ['t']``) the base ring

    - ``prefix`` -- (default: ``'S'``) the label for the shifting operators

    EXAMPLES::

        sage: S = ShiftingOperatorAlgebra()

        sage: elm = S[1, -1, 2]; elm
        S(1, -1, 2)
        sage: elm([5, 4])
        [([6, 3, 2], 1)]

    The shifting operator monomials can act on a complete homogeneous symmetric
    function or a Schur function::

        sage: s = SymmetricFunctions(QQ['t']).s()
        sage: h = SymmetricFunctions(QQ['t']).h()

        sage: elm(s[5, 4])
        s[6, 3, 2]
        sage: elm(h[5, 4])
        h[6, 3, 2]

        sage: S[1, -1](s[5, 4])
        s[6, 3]
        sage: S[1, -1](h[5, 4])
        h[6, 3]

    In fact, we can extend this action by linearity::

        sage: elm = (1 - S[1,-1]) * (1 - S[4])
        sage: elm == S([]) - S([1, -1]) - S([4]) + S([5, -1])
        True
        sage: elm(s[2, 2, 1])
        s[2, 2, 1] - s[3, 1, 1] - s[6, 2, 1] + s[7, 1, 1]

        sage: elm = (1 - S[1,-1]) * (1 - S[0,1,-1])
        sage: elm == 1 - S[0,1,-1] - S[1,-1] + S[1,0,-1]
        True
        sage: elm(s[2, 2, 1])
        s[2, 2, 1] - s[3, 1, 1] + s[3, 2]

    The algebra also comes equipped with homomorphisms to various
    symmetric function bases; these homomorphisms are how the action of
    ``S`` on the specific symmetric function bases is implemented::

        sage: elm = S([3,1,2]); elm
        S(3, 1, 2)
        sage: h(elm)
        h[3, 2, 1]
        sage: s(elm)
        0

    However, not all homomorphisms are equivalent, so the action is basis
    dependent::

        sage: elm = S([3,2,1]); elm
        S(3, 2, 1)
        sage: h(elm)
        h[3, 2, 1]
        sage: s(elm)
        s[3, 2, 1]
        sage: s(elm) == s(h(elm))
        False

    We can also use raising operators to implement the Jacobi-Trudi identity::

        sage: op = (1-S[(1,-1)]) * (1-S[(1,0,-1)]) * (1-S[(0,1,-1)])
        sage: s(op(h[3,2,1]))
        s[3, 2, 1]
    """

    def __init__(self, base_ring=QQ['t'], prefix='S'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ['t'])
            sage: TestSuite(S).run()
        """
        indices = ShiftingSequenceSpace()
        cat = Algebras(base_ring).WithBasis()
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         prefix=prefix,
                                         bracket=False, category=cat)

        # Setup default conversions
        sym = SymmetricFunctions(base_ring)
        self._sym_h = sym.h()
        self._sym_s = sym.s()
        self._sym_h.register_conversion(self.module_morphism(self._supp_to_h, codomain=self._sym_h))
        self._sym_s.register_conversion(self.module_morphism(self._supp_to_s, codomain=self._sym_s))

    def _repr_(self):
        r"""
        Return a string describing ``self``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S
            Shifting Operator Algebra over Univariate Polynomial Ring in t
            over Rational Field
        """
        return "Shifting Operator Algebra over {}".format(self.base_ring())

    def __getitem__(self, seq):
        r"""
        Return the shifting operator whose index is ``seq``.

        This method is only for basis indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S[1, 1, -9]
            S(1, 1, -9)
            sage: S[(1, 1, -9)]
            S(1, 1, -9)
            sage: S[5]
            S(5,)
        """
        if not isinstance(seq, tuple):
            if seq in ZZ:
                seq = (seq,)
            else:
                seq = tuple(seq)
        return self._element_constructor_(seq)

    def _prepare_seq(self, seq):
        r"""
        Standardize the input ``seq`` to be a basis element of ``self``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S._prepare_seq([0, 2, 0])
            (0, 2)
            sage: S._prepare_seq((0, 0, 0, 0))
            ()
            sage: S._prepare_seq(Partition([5,2,2,1]))
            (5, 2, 2, 1)
        """
        seq = tuple(seq)
        index = len(seq) - 1
        while index >= 0 and seq[index] == 0:
            index -= 1
        seq = seq[:index + 1]
        self._indices.check(seq)
        return seq

    def _element_constructor_(self, seq):
        r"""
        Return the shifting operator whose index is ``seq``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S._element_constructor_([1, 1, -9])
            S(1, 1, -9)
            sage: S[1, 1, -9]
            S(1, 1, -9)
            sage: S[0, 1, 0]
            S(0, 1)
        """
        if seq in self.base_ring():
            return self.term(self.one_basis(), self.base_ring()(seq))
        return self.monomial(self._prepare_seq(seq))

    def product_on_basis(self, x, y):
        r"""
        Return the product of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S.product_on_basis((0, 5, 2), (3, 2, -2, 5))
            S(3, 7, 0, 5)
            sage: S.product_on_basis((1, -2, 0, 3, -6), (-1, 2, 2))
            S(0, 0, 2, 3, -6)
            sage: S.product_on_basis((1, -2, -2), (-1, 2, 2))
            S()
        """
        # Make x have the longer length
        if len(x) < len(y):
            x, y = y, x
        x = list(x)  # Make a mutable copy
        for i, val in enumerate(y):
            x[i] += val
        # strip trailing 0's
        index = len(x) - 1
        while index >= 0 and x[index] == 0:
            index -= 1
        return self.monomial(tuple(x[:index + 1]))

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element for `1`.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S.one_basis()
            ()
        """
        return ()

    def _supp_to_h(self, supp):
        r"""
        This is a helper function that is not meant to be called directly.

        Given the support of an element
        `x_1^{\gamma_1} x_2^{\gamma_2} \cdots x_\ell^{\gamma_\ell}` in the
        ``ShiftingOperatorAlgebra`` and a
        symmetric function algebra basis `b` generated by
        `\{b_1, b_2, b_3,\ldots\}`, return the element
        `b_{\gamma_1} b_{\gamma_2} \cdots b_{\gamma_\ell}` where `b_0 = 1` and
        `b_{-n} = 0` for all positive
        integers `n`. The canonical example for `b` in this case would be `h`.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: h = SymmetricFunctions(QQ).h()
            sage: S._supp_to_h(S([3,2,1]).support_of_term())
            h[3, 2, 1]
            sage: S._supp_to_h(S([2,3,1]).support_of_term())
            h[3, 2, 1]
            sage: S._supp_to_h(S([2,3,-1]).support_of_term()) == h.zero()
            True
            sage: S._supp_to_h(S([2,3,0]).support_of_term())
            h[3, 2]
        """
        gamma = sorted(supp, reverse=True)
        if gamma in _Partitions:
            return self._sym_h(gamma)
        else:
            return self._sym_h.zero()

    def _supp_to_s(self, gamma):
        r"""
        This is a helper function that is not meant to be called directly.

        Given the support of an element `x_1^{\gamma_1} x_2^{\gamma_2} \cdots
        x_\ell^{\gamma_\ell}` in the ``ShiftingOperatorAlgebra``, return
        the appropriate `s_\gamma` in the Schur basis using
        "Schur function straightening" in [BMPS2018]_ Proposition 4.1.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: s = SymmetricFunctions(QQ).s()
            sage: S._supp_to_s(S([3,2,1]).support_of_term())
            s[3, 2, 1]
            sage: S._supp_to_s(S([2,3,1]).support_of_term())
            0
            sage: S._supp_to_s(S([2,4,-1,1]).support_of_term())
            s[3, 3]
            sage: S._supp_to_s(S([3,2,0]).support_of_term())
            s[3, 2]
        """
        def number_of_noninversions(lis):
            return sum(1 for i, val in enumerate(lis)
                       for j in range(i + 1, len(lis))
                       if val < lis[j])  # i < j is already enforced

        rho = list(range(len(gamma) - 1, -1, -1))
        combined = [g + r for g, r in zip(gamma, rho)]
        if len(set(combined)) == len(combined) and all(e >= 0
                                                       for e in combined):
            sign = (-1) ** number_of_noninversions(combined)
            sort_combined = sorted(combined, reverse=True)
            new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
            return sign * self._sym_s(_Partitions(new_gamma))
        else:
            return self._sym_s.zero()

    def build_and_register_conversion(self, support_map, codomain):
        r"""
        Build a module homomorphism from a map sending integer sequences to
        ``codomain`` and registers the result into Sage's conversion model.

        The intended use is to define a morphism from
        ``self`` to a basis `B` of symmetric functions that will be used by
        :class:`ShiftingOperatorAlgebra` to define the action of the
        operators on `B`.

        .. NOTE::

            The actions on the complete homogeneous symmetric functions and
            on the Schur functions by morphisms are already registered.

        .. WARNING::

            Because :class:`ShiftingOperatorAlgebra` inherits from
            :class:`UniqueRepresentation`, once you register a conversion, this
            will apply to all instances of :class:`ShiftingOperatorAlgebra`
            over the same base ring with the same prefix.

        INPUT:

        - ``support_map`` -- a map from integer sequences to ``codomain``

        - ``codomain`` -- the codomain of ``support_map``, usually a basis
          of symmetric functions

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: p = sym.p()
            sage: zero_map = lambda part: p.zero()
            sage: S.build_and_register_conversion(zero_map, p)
            sage: p(2*S([1,0,-1]) + S([2,1,0]) - 3*S([0,1,3]))
            0
            sage: op = S((1, -1))
            sage: op(2*p[4,3] + 5*p[2,2] + 7*p[2]) == p.zero()
            True

        For a more illustrative example, we can implement a simple
        (but not mathematically justified!) conversion on the monomial basis::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: m = sym.m()
            sage: def supp_map(gamma):
            ....:     gsort = sorted(gamma, reverse=True)
            ....:     return m(gsort) if gsort in Partitions() else m.zero()
            sage: S.build_and_register_conversion(supp_map, m)
            sage: op = S.ij(0, 1)
            sage: op(2*m[4,3] + 5*m[2,2] + 7*m[2]) == 2*m[5, 2] + 5*m[3, 1]
            True
        """
        module_morphism = self.module_morphism(support_map,
                                               codomain=codomain)
        codomain.register_conversion(module_morphism)

    def ij(self, i, j):
        r"""
        Return the raising operator `R_{ij}` as notated in [BMPS2018]_ Definition 2.1.

        Shorthand element constructor that allows you to create raising
        operators using the familiar `R_{ij}` notation found in
        [BMPS2018]_ Definition 2.1, with the exception that indices
        here are 0-based, not 1-based.

        EXAMPLES:

        Create the raising operator which raises part 0 and lowers part 2
        (indices are 0-based)::

            sage: R = ShiftingOperatorAlgebra()
            sage: R.ij(0, 2)
            S(1, 0, -1)
        """
        if i not in NonNegativeIntegerSemiring():
            raise ValueError('i (={}) must be a natural number'.format(i))
        if j not in NonNegativeIntegerSemiring():
            raise ValueError('j (={}) must be a natural number'.format(j))
        if not i < j:
            raise ValueError('index j (={j}) must be greater than index i (={i})'.format(i=i, j=j))
        seq = [0] * (max(i, j) + 1)
        seq[i] = 1
        seq[j] = -1
        return self._element_constructor_(seq)

    class Element(CombinatorialFreeModule.Element):
        r"""
        An element of a :class:`ShiftingOperatorAlgebra`.
        """

        def __call__(self, operand):
            r"""
            Call method for shifting sequence operators to act on objects.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra(QQ)
                sage: op = S([1,1]) + 2*S([0,1,0,1])
                sage: sorted(op([1,1,1,1])) == sorted([([2,2,1,1], 1), ([1,2,1,2], 2)])
                True
                sage: sorted(op(Partition([1,1,1,1]))) == sorted([([2,2,1,1], 1), ([1,2,1,2], 2)])
                True
                sage: sym = SymmetricFunctions(QQ)
                sage: h = sym.h()
                sage: op(h[1,1,1,1])
                3*h[2, 2, 1, 1]
                sage: s = sym.s()
                sage: op(s[1,1,1,1])
                s[2, 2, 1, 1]
                sage: e = sym.e()
                sage: sorted(op(e[1,1,1,1])) == sorted([([2,2,1,1], 1), ([1,2,1,2], 2)])
                True
            """
            P = self.parent()
            if isinstance(operand, (list, tuple, Composition, Partition)):
                def add_lists(x, y):
                    # Make x have the longer length
                    if len(x) < len(y):
                        x, y = y, x
                    x = list(x)  # Make a mutable copy
                    for i, val in enumerate(y):
                        x[i] += val
                    return x
                return [(add_lists(index, operand), coeff)
                        for index, coeff in self]

            R = self.base_ring()
            lift_operand = P._from_dict({P._prepare_seq(p): R(c)
                                         for p, c in operand}, coerce=False)
            result = self * lift_operand
            operand_parent = operand.parent()
            try:
                return operand_parent(result)
            except TypeError:
                return [(list(index), coeff) for index, coeff in result]
