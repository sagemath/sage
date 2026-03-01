r"""
Jacobians in the Unique Hess model

This module implements Jacobian arithmetic based on divisor representation by
ideals, much like :mod:`sage.rings.function_field.jacobian_hess`.
This approach to Jacobian arithmetic implementation is attributed to Hess [Hes2002]_.

The Unique Hess variant implemented here is a modification due to Macri [Mac2025]_
of Hess' maximally reduced divisors, but with a linear search reduction algorithm
that performs better in practice when doing a sequence of many Jacobian arithmetic operations.

The Unique Hess model is the only Jacobian model implemented in Sage that uses unique
representatives of divisor classes. To the best of the author's knowledge at time of writing,
this module makes SageMath the only computer algebra system that includes a model of Jacobian
arithmetic that uses unique representatives of divisor classes.

Jacobian
--------

To create a Jacobian in the Unique Hess model, specify ``'unique_hess'`` model and provide a base place of degree 1.

AUTHORS:

- Vincent Macri (2025-10-10): initial version
"""

# ****************************************************************************
#       Copyright (C) 2025 Vincent Macri <vincent.macri@ucaglary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING

import sage.groups.generic as groups_generic
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.rings.function_field import riemann_roch
from sage.structure.richcmp import op_EQ, op_NE, richcmp
from sage.structure.unique_representation import UniqueRepresentation

from .divisor import FunctionFieldDivisor
from .jacobian_base import (
    Jacobian_base,
    JacobianGroup_base,
    JacobianGroup_finite_field_base,
    JacobianPoint_base,
    JacobianPoint_finite_field_base,
)
from .place import FunctionFieldPlace

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Self, TypeAlias

    from sage.rings.integer import Integer

    from .function_field import FunctionField
    from .ideal import FunctionFieldIdeal
    from .ideal import FunctionFieldIdealInfinite as InfiniteIdeal
    FiniteIdeal: TypeAlias = FunctionFieldIdeal  # For readability when we specifically mean a finite ideal


class JacobianPoint(JacobianPoint_base):
    r"""
    A Jacobian point in the "Unique Hess" model.
    Points are represented by pairs of ideals as in the "Hess" model, but we
    use a different reduction algorithm that guarantees that all points have
    a unique representative. Unlike the other Jacobian models in Sage, elements
    of this class are hashable.
    """

    def __init__(self, parent: JacobianGroup, finite_ideal: FiniteIdeal,
                 infinite_ideal: InfiniteIdeal) -> None:
        super().__init__(parent)
        # self._r is determined by the two ideals, but it is computed during
        # reduction and is useful in a few utility methods so we store it.
        self._finite_ideal, self._infinite_ideal, self._r = parent._reduce(finite_ideal, infinite_ideal)

    def __hash__(self) -> int:
        r"""
        Return the hash of ``self``.
        """
        return hash((self._finite_ideal, self._infinite_ideal))

    def _richcmp_(self, other: Self, op) -> bool:
        r"""
        Compare ``self`` with ``other`` with respect to operator ``op``.
        """
        self_data = (self._finite_ideal, self._infinite_ideal)
        other_data = (other._finite_ideal, other._infinite_ideal)
        if op not in (op_EQ, op_NE):
            return richcmp(self_data, other_data, op)
        return (self_data == other_data) == (op is op_EQ)

    def _add_(self, other: Self) -> Self:
        r"""
        Add ``self`` and ``other``.
        """
        G = self.parent()
        finite_ideal = self._finite_ideal * other._finite_ideal
        infinite_ideal = G._infinite_ideal_mult(self._infinite_ideal, other._infinite_ideal)
        return G.element_class(G, finite_ideal, infinite_ideal)

    def _neg_(self) -> Self:
        G = self.parent()
        return G.element_class(G, ~self._finite_ideal, ~self._infinite_ideal)

    def additive_order(self) -> Integer:
        r"""
        Return the order of this point.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='unique_hess', base_div=b).group()
            sage: p = C([-1,2,1]).place()
            sage: pt = G.point(p - b)
            sage: pt.order()
            30
            sage: D = pt.divisor()
            sage: G = C.jacobian(model='unique_hess').group()
            sage: P = G(D)
            sage: P.order()
            30
        """
        bounds = (1, self.parent()._bound_on_order())
        return groups_generic.order_from_bounds(self, bounds)

    def effective_part(self) -> FunctionFieldDivisor:
        r"""
        Return the effective part of this divisor.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='unique_hess', base_div=b, extra_caching=False).group()
            sage: p = C([-1,2,1]).place()
            sage: pt = G.point(p - b)
            sage: pt.effective_part() == p
            True
        """
        return self.divisor() + self.parent().base_divisor() * self._r

    def divisor(self) -> FunctionFieldDivisor:
        r"""
        Return the unique reduced representative of this divisor class
        in the Unique Hess model with respect to the base place.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(29), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: G = C.jacobian(model='unique_hess', base_div=b).group()
            sage: p = C([-1,2,1]).place()
            sage: D = p - b
            sage: pt = G.point(D)
            sage: pt.divisor() == D
            True
            sage: pt == D
            False
        """
        return (~self._finite_ideal).divisor() + (~self._infinite_ideal).divisor()


class JacobianPoint_finite_field(JacobianPoint, JacobianPoint_finite_field_base):
    """
    Points of Jacobians over finite fields
    """


class JacobianGroup(UniqueRepresentation, JacobianGroup_base):
    Element = JacobianPoint

    def __init__(self, parent, function_field: FunctionField, base_div: FunctionFieldDivisor) -> None:
        r"""
        TESTS::

            sage: K = GF(11)
            sage: Kx.<x> = FunctionField(K)
            sage: t = polygen(Kx)
            sage: F.<y> = Kx.extension(t^4 + 9*x*t^3 + (10*x + 7)*t^2 + (7*x^2 + 2*x + 10)*t + 9*x^3 + 3*x^2 + 6*x + 4)
            sage: J = F.jacobian(model='unique_hess')
            sage: G = J.group()
            sage: TestSuite(G).run()
        """
        super().__init__(parent, function_field, base_div)

        A = parent._A
        if not base_div == A:
            raise ValueError('base_div must equal base_place')

        # We use inverted ideals, so _A_inv is actually the uninverted place
        self._A_inv = A.prime_ideal()
        self._A = self._A_inv.inverse()
        self._A_g_minus_1 = self._A ** (self._genus - 1)

        # For faster/more convenient access from the JacobianPoint class
        self._vector_space, self._from_vector_space, self._to_vector_space = self._function_field.free_module(map=True)
        self._A_is_infinite = A.is_infinite_place()
        self._function_field_degree = self._function_field.degree()
        self._maximal_order_finite = self._function_field.maximal_order()
        self._maximal_order_infinite = self._function_field.maximal_order_infinite()

        self._cache_infinite_ideals = parent._cache_infinite_ideals
        self._noncaching_infinite_ideal_mult = lambda J1, J2: J1 * J2
        self._noncaching_inverse_infinite_matrix = lambda J: matrix([self._to_vector_space(b) for b in J.gens_over_base()]).inverse()

        if self._cache_infinite_ideals:
            self._infinite_ideal_mult = self._cached_ideal_mult
            self._inverse_infinite_matrix = self._cached_inverse_infinite_matrix
        else:
            self._infinite_ideal_mult = self._noncaching_infinite_ideal_mult
            self._inverse_infinite_matrix = self._noncaching_inverse_infinite_matrix

        # Ideal multiplication is expensive, so we want to avoid unnecessary multiplication by identity.
        # We define functions to handle this so that we don't need to complicate our reduction logic with branching.
        if self._A_is_infinite:
            self._multiply_pair_by_A = lambda I, J, A: (I, self._infinite_ideal_mult(J, A))
        else:
            self._multiply_pair_by_A = lambda I, J, A: (I * A, J)

    def _reduce(self, I: FiniteIdeal, J: InfiniteIdeal) -> tuple[FiniteIdeal, InfiniteIdeal, int]:
        r"""
        We use Algorithm 4.1 from [Mac2025]_ to find a unique reduced
        representative of the divisor class. This algorithm is optimized
        for typical inputs which will be the vast majority of inputs in
        practice when performing a series of Jacobian computations.

        INPUT:

        - ``I`` -- an ideal of the finite maximal order
        - ``J`` -- an ideal of the infinite maximal order

        OUTPUT:

        A tuple `(I', J', r)` where `(I', J')` is the ideal representation of
        the unique reduced representative of the divisor class. The multiplicity
        of the base place of the output divisor is `-r \leq 0`.
        """
        # We take the input divisor D and find the maximal r such that ℓ(D + rA) = 1
        degree = self._function_field_degree
        to = self._to_vector_space

        def matrices_and_riemann_roch(C, B_inv, tilde_I, tilde_J):
            if self._A_is_infinite:
                B_inv = self._inverse_infinite_matrix(tilde_J)
                return C, B_inv, riemann_roch._short_circuit_riemann_roch_matrices(degree, C, B_inv, tilde_I.gens_over_base())
            C = matrix([to(v) for v in tilde_I.gens_over_base()])
            return C, B_inv, riemann_roch._short_circuit_riemann_roch_matrices(degree, C, B_inv, tilde_I.gens_over_base())

        tilde_I, tilde_J = self._multiply_pair_by_A(I, J, self._A_g_minus_1)
        C = matrix([to(v) for v in tilde_I.gens_over_base()])
        B_inv = self._inverse_infinite_matrix(tilde_J)

        last_basis = riemann_roch._short_circuit_riemann_roch_matrices(degree, C, B_inv, tilde_I.gens_over_base())
        basis = last_basis  # Needed to ensure basis is defined in the g = 1 case

        if last_basis is None:  # ell(D + (g - 1) A) = 0
            r = self._genus
            tilde_I, tilde_J = self._multiply_pair_by_A(tilde_I, tilde_J, self._A)  # D + gA
            C, B_inv, basis = matrices_and_riemann_roch(C, B_inv, tilde_I, tilde_J)
        else:
            r = 0
            last_tilde_I, last_tilde_J = tilde_I, tilde_J
            for n in range(self._genus - 2, -1, -1):
                tilde_I, tilde_J = self._multiply_pair_by_A(tilde_I, tilde_J, self._A_inv)  # D + nA
                C, B_inv, basis = matrices_and_riemann_roch(C, B_inv, tilde_I, tilde_J)
                if basis is None:
                    # ℓ(D + nA) = 0, hence ℓ(D + (n + 1) A) = 1
                    basis = last_basis
                    r = n + 1
                    tilde_I, tilde_J = last_tilde_I, last_tilde_J
                    break
                last_tilde_I, last_tilde_J = tilde_I, tilde_J
                last_basis = basis
            # If we never hit the breakpoint then ℓ(D) = 1, which happens iff D = 0

        f_inv = ~basis

        fI = self._maximal_order_finite.ideal(f_inv)
        fJ = self._maximal_order_infinite.ideal(f_inv)

        return I * fI, self._infinite_ideal_mult(J, fJ), r

    def _use_caching(self, caching: bool):
        r"""
        Change whether we use caching or not.

        Used to internally turn off caching when reducing elements provided by
        the user, as valid input from users exceeds the values that will be
        encountered when doing a series of operations with reduced elements,
        causing potentially unbounded memory use.

        If the Jacobian model was not initialized to use caching, this method does nothing.

        TESTS::

            sage: K = GF(2)
            sage: Kx.<x> = FunctionField(K)
            sage: t = polygen(Kx)
            sage: F.<y> = Kx.extension(t^3 + (x^2 + x + 1)*t^2 + (x^3 + x + 1)*t + x^5 + x^4)
            sage: J = F.jacobian(model='unique_hess', extra_caching=True)
            sage: G = J.group()
            sage: len(G._infinite_ideal_mult.cache)
            0
            sage: len(G._inverse_infinite_matrix.cache)
            0
            sage: zero = J(0)
            sage: len(G._infinite_ideal_mult.cache)
            0
            sage: len(G._inverse_infinite_matrix.cache)
            0
            sage: D1, D2, D3 = G.get_points(3)
            sage: len(G._infinite_ideal_mult.cache)
            0
            sage: len(G._inverse_infinite_matrix.cache)
            0
            sage: _ = D2 + D3
            sage: len(G._infinite_ideal_mult.cache) > 0
            True
            sage: len(G._inverse_infinite_matrix.cache) > 0
            True
        """
        if not self._cache_infinite_ideals:
            pass

        if caching:
            self._infinite_ideal_mult = self._cached_ideal_mult
            self._inverse_infinite_matrix = self._cached_inverse_infinite_matrix
        else:
            self._infinite_ideal_mult = self._noncaching_infinite_ideal_mult
            self._inverse_infinite_matrix = self._noncaching_inverse_infinite_matrix


    @cached_method
    def _cached_inverse_infinite_matrix(self, J):
        r"""
        During the computation of a Riemann-Roch space we need to invert
        the matrix corresponding to the infinite ideal. Matrix inversion
        is somewhat expensive, and there are few distinct infinite ideals
        that we will wish to compute this for, so we cache the results.
        """
        return matrix([self._to_vector_space(b) for b in J.gens_over_base()]).inverse()

    @cached_method(key=lambda self, J1, J2: frozenset((J1, J2)))
    def _cached_ideal_mult(self, J1: InfiniteIdeal, J2: InfiniteIdeal):  # noqa: PLR6301 - this is a method so the cache is deleted when the Jacobian instance is
        r"""
        Cached wrapper around multiplication of infinite ideals.
        Because ideal multiplication is commutative, we use a frozenset as the
        cache key, so we don't cache the same multiplication performed in a
        different order.

        When doing a sequence of computations in the Jacobian of a function field,
        in practice we will only encounter a handful of distinct infinite ideals.
        This means that we can cache the result of infinite ideal arithmetic with
        minimal cost.
        """
        return J1 * J2

    def _element_constructor_(self, x):
        if x == 0:
            return self.zero()

        if isinstance(x, FunctionFieldPlace):
            x = x.divisor()

        if isinstance(x, FunctionFieldDivisor) and x in self._function_field.divisor_group():
            return self.point(x)

        raise ValueError(f'cannot construct a point from {x}')

    def point(self, divisor: FunctionFieldDivisor):
        r"""
        Return the point represented by the divisor of degree zero.
        """
        if divisor.degree() != 0:
            raise ValueError('divisor not of degree zero')

        I, J = riemann_roch._divisor_to_inverted_ideals(divisor)
        self._use_caching(False)
        D = self.element_class(self, I, J)
        self._use_caching(True)  # Does nothing if we weren't already using caching
        return D

    @cached_method
    def zero(self) -> JacobianPoint:
        r"""
        Return the zero element of this group.
        """
        return self.point(self._function_field.divisor_group().zero())

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.
        """
        return f'{super()._repr_()} (Unique Hess model)'


class JacobianGroup_finite_field(JacobianGroup, JacobianGroup_finite_field_base):
    Element = JacobianPoint_finite_field

    def __iter__(self) -> Iterable[JacobianPoint_finite_field]:
        dg = self.function_field().divisor_group()
        g = self._genus
        A = self._base_div
        hits = set()
        for D in dg.effective_divisors(max_degree=g, avoid=[A]):
            pt = D - D.degree() * A
            assert pt.degree() == 0
            point = self.point(pt)
            if point in hits:
                continue
            hits.add(point)
            yield point


class Jacobian(Jacobian_base, UniqueRepresentation):

    def __init__(self, function_field: FunctionField, base_div: FunctionFieldDivisor | FunctionFieldPlace, cache_infinite_ideals: bool = True, **kwds) -> None:
        r"""
        TESTS::

            sage: K = GF(11)
            sage: Kx.<x> = FunctionField(K)
            sage: t = polygen(Kx)
            sage: F.<y> = Kx.extension(t^4 + 9*x*t^3 + (10*x + 7)*t^2 + (7*x^2 + 2*x + 10)*t + 9*x^3 + 3*x^2 + 6*x + 4)
            sage: J = F.jacobian(model='unique_hess')
            sage: TestSuite(J).run()

        Test projective curve::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: J = C.jacobian(model='unique_hess')
            sage: TestSuite(J).run()
        """
        if base_div.degree() != 1:
            raise ValueError('base_div must be degree 1')

        if isinstance(base_div, FunctionFieldPlace):
            super().__init__(function_field, base_div.divisor(), **kwds)
            self._A = base_div
        elif isinstance(base_div, FunctionFieldDivisor):  # Allowed for compatibility with other Jacobian models
            if not base_div.is_prime():
                raise ValueError('base_div must be a prime divisor')
            super().__init__(function_field, base_div, **kwds)
            self._A = base_div.place()
        else:
            raise TypeError('base_div must be a divisor or a place')

        #reveal_type(self._A)

        self._cache_infinite_ideals = cache_infinite_ideals

        if function_field.constant_base_field().is_finite():
            self._group_class = JacobianGroup_finite_field
        else:
            self._group_class = JacobianGroup

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3 - y^2*z, P2)
            sage: C.jacobian(model='unique_hess')
            Jacobian of Projective Plane Curve over Finite Field of size 7
             defined by x^3 - y^2*z - 2*z^3 (Unique Hess model)
        """
        return f'{super()._repr_()} (Unique Hess model)'

    def group(self, k_ext=None):
        r"""
        Raise an error if ``k_ext is not self.function_field().constant_base_field()``,
        otherwise call :meth:`Jacobian_base.group`.
        """
        if k_ext is not None and k_ext is not self._function_field.constant_base_field():
            raise NotImplementedError('Embeddings are not yet supported for the Unique Hess Jacobian model')
        return super().group(k_ext)
