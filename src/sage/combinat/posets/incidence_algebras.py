# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.modules
r"""
Incidence Algebras
"""
# ****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.algebras import Algebras
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.matrix_space import MatrixSpace

from copy import copy


class IncidenceAlgebra(CombinatorialFreeModule):
    r"""
    The incidence algebra of a poset.

    Let `P` be a poset and `R` be a commutative unital associative ring.
    The *incidence algebra* `I_P` is the algebra of functions
    `\alpha \colon P \times P \to R` such that `\alpha(x, y) = 0`
    if `x \not\leq y` where multiplication is given by convolution:

    .. MATH::

        (\alpha \ast \beta)(x, y) = \sum_{x \leq k \leq y}
        \alpha(x, k) \beta(k, y).

    This has a natural basis given by indicator functions for the
    interval `[a, b]`, i.e. `X_{a,b}(x,y) = \delta_{ax} \delta_{by}`.
    The incidence algebra is a unital algebra with the identity given
    by the Kronecker delta `\delta(x, y) = \delta_{xy}`. The Möbius
    function of `P` is another element of `I_p` whose inverse is the
    `\zeta` function of the poset (so `\zeta(x, y) = 1` for
    every interval `[x, y]`).

    .. TODO::

        Implement the incidence coalgebra.

    REFERENCES:

    - :wikipedia:`Incidence_algebra`
    """
    def __init__(self, R, P, prefix='I') -> None:
        """
        Initialize ``self``.

        TESTS::

            sage: P = posets.BooleanLattice(3)
            sage: I = P.incidence_algebra(QQ)
            sage: TestSuite(I).run()  # long time
        """
        cat = Algebras(R).WithBasis()
        if P in FiniteEnumeratedSets():
            cat = cat.FiniteDimensional()
        self._poset = P
        CombinatorialFreeModule.__init__(self, R, map(tuple, P.relations()),
                                         prefix=prefix, category=cat)

    def _repr_term(self, A):
        """
        Return a string representation of the term labeled by ``A``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I._repr_term((4, 12))
            'I[4, 12]'
        """
        return self.prefix() + str(list(A))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: P.incidence_algebra(QQ)
            Incidence algebra of Finite lattice containing 16 elements
             over Rational Field
        """
        return f"Incidence algebra of {self._poset} over {self.base_ring()}"

    def _coerce_map_from_(self, R):
        """
        Return the coerce map from ``R`` into ``self`` if it exists
        or ``False`` otherwise.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: R = I.reduced_subalgebra()
            sage: I.has_coerce_map_from(R)
            True
            sage: I.has_coerce_map_from(QQ)
            True
            sage: Pp = posets.BooleanLattice(3)
            sage: Rp = Pp.incidence_algebra(QQ).reduced_subalgebra()
            sage: I.has_coerce_map_from(Rp)
            False
        """
        if isinstance(R, ReducedIncidenceAlgebra) and R._ambient is self:
            return copy(R.lift)
        return super()._coerce_map_from_(R)

    def reduced_subalgebra(self, prefix='R'):
        """
        Return the reduced incidence subalgebra.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.reduced_subalgebra()
            Reduced incidence algebra of Finite lattice containing 16 elements
             over Rational Field
        """
        return ReducedIncidenceAlgebra(self, prefix)

    def poset(self):
        """
        Return the defining poset of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.poset()
            Finite lattice containing 16 elements
            sage: I.poset() == P
            True
        """
        return self._poset

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(1)
            sage: I = P.incidence_algebra(QQ)
            sage: Ielts = I.some_elements(); Ielts # random
            [2*I[0, 0] + 2*I[0, 1] + 3*I[1, 1],
             I[0, 0] - I[0, 1] + I[1, 1],
             I[0, 0] + I[0, 1] + I[1, 1]]
            sage: [a in I for a in Ielts]
            [True, True, True]
        """
        return [self.an_element(), self.moebius(), self.zeta()]

    def product_on_basis(self, A, B):
        r"""
        Return the product of basis elements indexed by ``A`` and ``B``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.product_on_basis((1, 3), (3, 11))
            I[1, 11]
            sage: I.product_on_basis((1, 3), (2, 2))
            0
        """
        if A[1] == B[0]:
            return self.monomial((A[0], B[1]))
        return self.zero()

    @cached_method
    def one(self):
        r"""
        Return the element `1` in ``self`` (which is the Kronecker
        delta `\delta(x, y)`).

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.one()
            I[0, 0] + I[1, 1] + I[2, 2] + I[3, 3] + I[4, 4] + I[5, 5]
             + I[6, 6] + I[7, 7] + I[8, 8] + I[9, 9] + I[10, 10]
             + I[11, 11] + I[12, 12] + I[13, 13] + I[14, 14] + I[15, 15]
        """
        return self.sum_of_monomials((x, x) for x in self._poset)

    delta = one

    @cached_method
    def zeta(self):
        r"""
        Return the `\zeta` function in ``self``.

        The `\zeta` function on a poset `P` is given by

        .. MATH::

            \zeta(x, y) = \begin{cases}
            1 & x \leq y, \\
            0 & x \not\leq y.
            \end{cases}

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.zeta() * I.moebius() == I.one()
            True
        """
        return self.sum(self.basis())

    @cached_method
    def moebius(self):
        """
        Return the Möbius function of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(2)
            sage: I = P.incidence_algebra(QQ)
            sage: I.moebius()
            I[0, 0] - I[0, 1] - I[0, 2] + I[0, 3] + I[1, 1]
             - I[1, 3] + I[2, 2] - I[2, 3] + I[3, 3]
        """
        mu = self._poset.moebius_function
        R = self.base_ring()
        return self.sum_of_terms((A, R(mu(*A))) for A in self.basis().keys())

    def __getitem__(self, A):
        """
        Return the basis element indexed by ``A``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I[2, 6]
            I[2, 6]
            sage: I[(2, 6)]
            I[2, 6]
            sage: I[[2, 6]]
            I[2, 6]
            sage: I[2]
            I[2, 2]

            sage: I[2, 5]
            Traceback (most recent call last):
            ...
            ValueError: not an interval

            sage: I[-1]
            Traceback (most recent call last):
            ...
            ValueError: not an element of the poset
        """
        if not isinstance(A, (list, tuple)):
            if A not in self._poset.list():
                raise ValueError("not an element of the poset")
            A = (A, A)
        else:
            A = tuple(A)
            if len(A) != 2 or A not in self.basis().keys():
                raise ValueError("not an interval")
        return self.monomial(A)

    @lazy_attribute
    def _linear_extension(self):
        """
        Return a fixed linear extension of the defining poset of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I._linear_extension
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
        """
        return tuple(self._poset.linear_extension())

    class Element(CombinatorialFreeModule.Element):
        """
        An element of an incidence algebra.
        """
        def __call__(self, x, y):
            """
            Return ``self(x, y)``.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.moebius()
                sage: mu(0, 12)
                1
                sage: mu(0, 7)
                -1
                sage: mu(15, 0)
                0
            """
            P = self.parent()._poset
            x = P(x)
            y = P(y)
            return self[x, y]

        @cached_method
        def to_matrix(self):
            r"""
            Return ``self`` as a matrix.

            We define a matrix `M_{xy} = \alpha(x, y)` for some element
            `\alpha \in I_P` in the incidence algebra `I_P` and we order
            the elements `x,y \in P` by some linear extension of `P`. This
            defines an algebra (iso)morphism; in particular, multiplication
            in the incidence algebra goes to matrix multiplication.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: I.moebius().to_matrix()
                [ 1 -1 -1  1]
                [ 0  1  0 -1]
                [ 0  0  1 -1]
                [ 0  0  0  1]
                sage: I.zeta().to_matrix()
                [1 1 1 1]
                [0 1 0 1]
                [0 0 1 1]
                [0 0 0 1]

            TESTS:

            We check that this is an algebra (iso)morphism::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.moebius()
                sage: (mu*mu).to_matrix() == mu.to_matrix() * mu.to_matrix()
                True
            """
            P = self.parent()
            MS = MatrixSpace(P.base_ring(), P._poset.cardinality(), sparse=True)
            L = P._linear_extension
            M = copy(MS.zero())
            for i, c in self:
                M[L.index(i[0]), L.index(i[1])] = c
            M.set_immutable()
            return M

        def is_unit(self):
            """
            Return if ``self`` is a unit.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.moebius()
                sage: mu.is_unit()
                True
                sage: zeta = I.zeta()
                sage: zeta.is_unit()
                True
                sage: x = mu - I.zeta() + I[2,2]
                sage: x.is_unit()
                False
                sage: y = I.moebius() + I.zeta()
                sage: y.is_unit()
                True

            This depends on the base ring::

                sage: I = P.incidence_algebra(ZZ)
                sage: y = I.moebius() + I.zeta()
                sage: y.is_unit()
                False
            """
            return all(self[x, x].is_unit() for x in self.parent()._poset)

        def __invert__(self):
            """
            Return the inverse of ``self`` if it exists.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.moebius()
                sage: ~mu
                I[0, 0] + I[0, 1] + I[0, 2] + I[0, 3] + I[1, 1]
                 + I[1, 3] + I[2, 2] + I[2, 3] + I[3, 3]
                sage: x = mu - I.zeta() + I[2,2]
                sage: ~x
                Traceback (most recent call last):
                ...
                ValueError: element is not invertible

            TESTS::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.moebius()
                sage: ~mu == I.zeta()
                True
                sage: ~I.one() == I.one()
                True
            """
            M = self.to_matrix()
            if not M.is_invertible():
                raise ValueError("element is not invertible")
            inv = ~M
            L = self.parent()._linear_extension
            return self.parent().sum_of_terms(
                ((L[i], L[j]), inv[i, j])
                for i, j in inv.nonzero_positions(copy=False))


class ReducedIncidenceAlgebra(CombinatorialFreeModule):
    r"""
    The reduced incidence algebra of a poset.

    The reduced incidence algebra `R_P` is a subalgebra of the
    incidence algebra `I_P` where `\alpha(x, y) = \alpha(x', y')` when
    `[x, y]` is isomorphic to `[x', y']` as posets. Thus the delta, Möbius,
    and zeta functions are all elements of `R_P`.
    """
    def __init__(self, I, prefix='R') -> None:
        """
        Initialize ``self``.

        TESTS::

            sage: P = posets.BooleanLattice(3)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: TestSuite(R).run()  # long time
        """
        self._ambient = I
        EC = {}
        P = self._ambient._poset
        if not P.is_finite():
            raise NotImplementedError("only implemented for finite posets")
        for i in self._ambient.basis().keys():
            S = P.subposet(P.interval(*i))
            added = False
            for k in EC:
                if S._hasse_diagram.is_isomorphic(k._hasse_diagram):
                    EC[k].append(i)
                    added = True
                    break
            if not added:
                EC[S] = [i]
        self._equiv_classes = map(sorted, EC.values())
        self._equiv_classes = {cls[0]: cls for cls in self._equiv_classes}
        cat = Algebras(I.base_ring()).FiniteDimensional().WithBasis()
        CombinatorialFreeModule.__init__(self, I.base_ring(),
                                         sorted(self._equiv_classes.keys()),
                                         prefix=prefix, category=cat)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: P.incidence_algebra(QQ).reduced_subalgebra()
            Reduced incidence algebra of Finite lattice containing 16 elements
             over Rational Field
        """
        msg = "Reduced incidence algebra of {} over {}"
        return msg.format(self._ambient._poset, self.base_ring())

    def poset(self):
        """
        Return the defining poset of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.poset()
            Finite lattice containing 16 elements
            sage: R.poset() == P
            True
        """
        return self._ambient._poset

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.some_elements()
            [2*R[(0, 0)] + 2*R[(0, 1)] + 3*R[(0, 3)],
             R[(0, 0)] - R[(0, 1)] + R[(0, 3)] - R[(0, 7)] + R[(0, 15)],
             R[(0, 0)] + R[(0, 1)] + R[(0, 3)] + R[(0, 7)] + R[(0, 15)]]
        """
        return [self.an_element(), self.moebius(), self.zeta()]

    @cached_method
    def one_basis(self):
        """
        Return the index of the element `1` in ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.one_basis()
            (0, 0)
        """
        for A in self.basis().keys():
            if A[0] == A[1]:
                return A

    def delta(self):
        """
        Return the Kronecker delta function in ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.delta()
            R[(0, 0)]
        """
        return self.one()

    @cached_method
    def zeta(self):
        r"""
        Return the `\zeta` function in ``self``.

        The `\zeta` function on a poset `P` is given by

        .. MATH::

            \zeta(x, y) = \begin{cases}
            1 & x \leq y, \\
            0 & x \not\leq y.
            \end{cases}

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.zeta()
            R[(0, 0)] + R[(0, 1)] + R[(0, 3)] + R[(0, 7)] + R[(0, 15)]
        """
        return self.sum(self.basis())

    @cached_method
    def moebius(self):
        """
        Return the Möbius function of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.moebius()
            R[(0, 0)] - R[(0, 1)] + R[(0, 3)] - R[(0, 7)] + R[(0, 15)]
        """
        mu = self._ambient._poset.moebius_function
        R = self.base_ring()
        return self.sum_of_terms((A, R(mu(*A))) for A in self.basis().keys())

    @cached_method
    def _lift_basis(self, x):
        """
        Lift the basis element indexed by ``x`` to the ambient incidence
        algebra of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R._lift_basis((0, 7))
            I[0, 7] + I[0, 11] + I[0, 13] + I[0, 14] + I[1, 15]
             + I[2, 15] + I[4, 15] + I[8, 15]
        """
        return self._ambient.sum_of_monomials(self._equiv_classes[x])

    @lazy_attribute
    def lift(self):
        """
        Return the lift morphism from ``self`` to the ambient space.

        EXAMPLES::

            sage: P = posets.BooleanLattice(2)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R.lift
            Generic morphism:
              From: Reduced incidence algebra of Finite lattice containing 4 elements over Rational Field
              To:   Incidence algebra of Finite lattice containing 4 elements over Rational Field
            sage: R.an_element() - R.one()
            R[(0, 0)] + 2*R[(0, 1)] + 3*R[(0, 3)]
            sage: R.lift(R.an_element() - R.one())
            I[0, 0] + 2*I[0, 1] + 2*I[0, 2] + 3*I[0, 3] + I[1, 1]
             + 2*I[1, 3] + I[2, 2] + 2*I[2, 3] + I[3, 3]
        """
        return self.module_morphism(self._lift_basis, codomain=self._ambient)

    def __getitem__(self, A):
        """
        Return the basis element indexed by ``A``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
            sage: R[0, 0]
            R[(0, 0)]
            sage: R[0, 1]
            R[(0, 1)]
            sage: R[0, 15]
            R[(0, 15)]
            sage: R[0]
            R[(0, 0)]

        This works for any representative of the equivalence class::

            sage: R[3, 3]
            R[(0, 0)]
            sage: R[3, 11]
            R[(0, 1)]

        TESTS:

            sage: R[2, 5]
            Traceback (most recent call last):
            ...
            ValueError: not an interval

            sage: R[-1]
            Traceback (most recent call last):
            ...
            ValueError: not an element of the poset
        """
        if not isinstance(A, (list, tuple)):
            if A not in self._ambient._poset.list():
                raise ValueError("not an element of the poset")
            return self.one()
        else:
            A = tuple(A)
            if len(A) != 2:
                raise ValueError("not an interval")
        for k in self._equiv_classes:
            if A in self._equiv_classes[k]:
                return self.monomial(k)
        raise ValueError("not an interval")

    def _retract(self, x):
        """
        Return the retract of ``x`` from the incidence algebra to ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: R = I.reduced_subalgebra()
            sage: all(R._retract(R.lift(x)) == x for x in R.basis())
            True
            sage: R._retract(I.zeta()) == R.zeta()
            True
            sage: R._retract(I.delta()) == R.delta()
            True
            sage: R._retract(I.moebius()) == R.moebius()
            True
        """
        return self.sum_of_terms((k, x[k]) for k in self.basis().keys())

    class Element(CombinatorialFreeModule.Element):
        """
        An element of a reduced incidence algebra.
        """
        def __call__(self, x, y):
            """
            Return ``self(x, y)``.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
                sage: x = R.an_element()
                sage: x(0, 7)
                0
                sage: x(7, 15)
                2
                sage: x(4, 15)
                0
            """
            return self.parent().lift(self)(x, y)

        def _mul_(self, other):
            """
            Return the product of ``self`` and ``other``.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
                sage: x = R.an_element()
                sage: x * R.zeta()
                2*R[(0, 0)] + 4*R[(0, 1)] + 9*R[(0, 3)] + 17*R[(0, 7)] + 28*R[(0, 15)]
                sage: x * R.moebius()
                2*R[(0, 0)] + R[(0, 3)] - 5*R[(0, 7)] + 12*R[(0, 15)]
                sage: x * R.moebius() * R.zeta() == x
                True
            """
            P = self.parent()
            return P._retract(P.lift(self) * P.lift(other))

        @cached_method
        def to_matrix(self):
            r"""
            Return ``self`` as a matrix.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
                sage: mu = R.moebius()
                sage: mu.to_matrix()
                [ 1 -1 -1  1]
                [ 0  1  0 -1]
                [ 0  0  1 -1]
                [ 0  0  0  1]
            """
            return self.parent().lift(self).to_matrix()

        def is_unit(self):
            """
            Return if ``self`` is a unit.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
                sage: x = R.an_element()
                sage: x.is_unit()
                True
            """
            return self[self.parent().one_basis()].is_unit()

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: R = P.incidence_algebra(QQ).reduced_subalgebra()
                sage: x = R.an_element()
                sage: ~x
                1/2*R[(0, 0)] - 1/2*R[(0, 1)] + 1/4*R[(0, 3)]
                 + 3/2*R[(0, 7)] - 33/4*R[(0, 15)]
            """
            P = self.parent()
            return P._retract(~P.lift(self))

        def lift(self):
            """
            Return the lift of ``self`` to the ambient space.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: R = I.reduced_subalgebra()
                sage: x = R.an_element(); x
                2*R[(0, 0)] + 2*R[(0, 1)] + 3*R[(0, 3)]
                sage: x.lift()
                2*I[0, 0] + 2*I[0, 1] + 2*I[0, 2] + 3*I[0, 3] + 2*I[1, 1]
                 + 2*I[1, 3] + 2*I[2, 2] + 2*I[2, 3] + 2*I[3, 3]
            """
            return self.parent().lift(self)
