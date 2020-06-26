"""
Askey-Wilson Algebras

AUTHORS:

- Travis Scrimshaw (2018-08): initial version
"""

# ****************************************************************************
#  Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebras import Algebras
from sage.categories.cartesian_product import cartesian_product
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.modules.with_basis.morphism import ModuleMorphismByLinearity
from sage.sets.family import Family
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.sets.non_negative_integers import NonNegativeIntegers


class AskeyWilsonAlgebra(CombinatorialFreeModule):
    r"""
    The (universal) Askey-Wilson algebra.

    Let `R` be a commutative ring. The *universal Askey-Wilson* algebra
    is an associative unital algebra `\Delta_q` over `R[q,q^-1]` given
    by the generators `A, B, C, \alpha, \beta, \gamma` that satisfy the
    following relations:

    .. MATH::

        \begin{aligned}
        (q-q^{-1}) \alpha &= (q^2-q^{-2}) A + qBC - q^{-1}CB, \\
        (q-q^{-1}) \beta  &= (q^2-q^{-2}) B + qCA - q^{-1}AC, \\
        (q-q^{-1}) \gamma &= (q^2-q^{-2}) C + qAB - q^{-1}BA.
        \end{aligned}

    The universal Askey-Wilson contains a
    :meth:`Casimir element <casimir_element>` `\Omega`, and the elements
    `\alpha`, `\beta`, `\gamma`, `\Omega` generate the center of `\Delta_q`,
    which is isomorphic to the polynomial ring
    `(R[q,q^-1])[\alpha,\beta,\gamma,\Omega]` (assuming `q` is not a root
    of unity). Furthermore, the relations imply that `\Delta_q` has a basis
    given by monomials `A^i B^j C^k \alpha^r \beta^s \gamma^t`, where
    `i, j, k, r, s, t \in \ZZ_{\geq 0}`.

    The universal Askey-Wilson algebra also admits a faithful action
    of `PSL_2(\ZZ)` given by the automorphisms `\rho`
    (:meth:`permutation_automorphism`):

    .. MATH::

        A \mapsto B \mapsto C \mapsto A,
        \qquad
        \alpha \mapsto \beta \mapsto \gamma \mapsto \alpha.

    and `\sigma` (:meth:`reflection_automorphism`):

    .. MATH::

        A \mapsto B \mapsto A,
        C \mapsto C + \frac{AB - BA}{q-q^{-1}},
        \qquad
        \alpha \mapsto \beta \mapsto \alpha,
        \gamma \mapsto \gamma.

    Note that `\rho^3 = \sigma^2 = 1` and

    .. MATH::

        \sigma(C) = C - q AB - (1+q^2) C + q \gamma
                  = C - q AB - q^2 C + q \gamma.

    The Askey-Wilson `AW_q(a,b,c)` algebra is a specialization of the
    universal Askey-Wilson algebra by `\alpha = a`, \beta = b`,
    `\gamma = c`, where `a,b,c \in R`. `AW_q(a,b,c)` was first introduced
    by [Zhedanov1991]_ to describe the Askey-Wilson polynomials. The
    Askey-Wilson algebra has a central extension of `\Delta_q`.

    INPUT:

    - ``R`` -- a commutative ring
    - ``q`` -- (optional) the parameter `q`; must be invertible in ``R``

    If ``q`` is not specified, then ``R`` is taken to be the base
    ring of a Laurent polynomial ring with variable `q`. Otherwise
    the element ``q`` must be an element of ``R``.

    .. NOTE::

        No check is performed to ensure ``q`` is not a root of unity,
        which may lead to violations of the results in [Terwilliger2011]_.

    EXAMPLES:

    We create the universal Askey-Wilson algebra and check
    the defining relations::

        sage: AW = algebras.AskeyWilson(QQ)
        sage: AW.inject_variables()
        Defining A, B, C, a, b, g
        sage: q = AW.q()
        sage: (q^2-q^-2)*A + q*B*C - q^-1*C*B == (q-q^-1)*a
        True
        sage: (q^2-q^-2)*B + q*C*A - q^-1*A*C == (q-q^-1)*b
        True
        sage: (q^2-q^-2)*C + q*A*B - q^-1*B*A == (q-q^-1)*g
        True

    Next, we perform some computations::

        sage: C * A
        (q^-2)*A*C + (q^-3-q)*B - (q^-2-1)*b
        sage: B^2 * g^2 * A
        q^4*A*B^2*g^2 - (q^-1-q^7)*B*C*g^2 + (1-q^4)*B*g^3
         + (1-2*q^4+q^8)*A*g^2 - (q-q^3-q^5+q^7)*a*g^2
        sage: (B^3 - A) * (C^2 + q*A*B)
        q^7*A*B^4 + B^3*C^2 - (q^2-q^14)*B^3*C + (q-q^7)*B^3*g - q*A^2*B
         + (3*q^3-4*q^7+q^19)*A*B^2 - A*C^2 - (1-q^6-q^8+q^14)*B^2*a
         - (q^-2-3*q^6+3*q^14-q^22)*B*C
         + (q^-1+q-3*q^3-q^5+2*q^7-q^9+q^13+q^15-q^19)*B*g
         + (2*q^-1-6*q^3+5*q^7-2*q^19+q^23)*A
         - (2-2*q^2-4*q^4+4*q^6+q^8-q^10+q^12-q^14+q^16-q^18-q^20+q^22)*a

    We check the elements `\alpha`, `\beta`, and `\gamma`
    are in the center::

        sage: all(x * gen == gen * x for gen in AW.algebra_generators() for x in [a,b,g])
        True

    We verify that the :meth:`Casimir element <casimir_element>`
    is in the center::

        sage: Omega = AW.casimir_element()
        sage: all(x * Omega == Omega * x for x in [A,B,C])
        True

        sage: x = AW.an_element()
        sage: O2 = Omega^2
        sage: x * O2 == O2 * x
        True

    We prove Lemma 2.1 in [Terwilliger2011]_::

        sage: (q^2-q^-2) * C == (q-q^-1) * g - (q*A*B - q^-1*B*A)
        True
        sage: (q-q^-1) * (q^2-q^-2) * a == (B^2*A - (q^2+q^-2)*B*A*B + A*B^2
        ....:                               + (q^2-q^-2)^2*A + (q-q^-1)^2*B*g)
        True
        sage: (q-q^-1) * (q^2-q^-2) * b == (A^2*B - (q^2+q^-2)*A*B*A + B*A^2
        ....:                               + (q^2-q^-2)^2*B + (q-q^-1)^2*A*g)
        True

    We prove Theorem 2.2 in [Terwilliger2011]_::

        sage: q3 = q^-2 + 1 + q^2
        sage: A^3*B - q3*A^2*B*A + q3*A*B*A^2 - B*A^3 == -(q^2-q^-2)^2 * (A*B - B*A)
        True
        sage: B^3*A - q3*B^2*A*B + q3*B*A*B^2 - A*B^3 == -(q^2-q^-2)^2 * (B*A - A*B)
        True
        sage: (A^2*B^2 - B^2*A^2 + (q^2+q^-2)*(B*A*B*A-A*B*A*B)
        ....:  == -(q^1-q^-1)^2 * (A*B - B*A) * g)
        True

    We construct an Askey-Wilson algebra over `\GF{5}` at `q=2`::

        sage: AW = algebras.AskeyWilson(GF(5), q=2)
        sage: A,B,C,a,b,g = AW.algebra_generators()
        sage: q = AW.q()
        sage: Omega = AW.casimir_element()

        sage: B * A
        4*A*B + 2*g
        sage: C * A
        4*A*C + 2*b
        sage: C * B
        4*B*C + 2*a
        sage: Omega^2
        A^2*B^2*C^2 + A^3*B*C + A*B^3*C + A*B*C^3 + A^4 + 4*A^3*a
         + 2*A^2*B^2 + A^2*B*b + 2*A^2*C^2 + 4*A^2*C*g + 4*A^2*a^2
         + 4*A*B^2*a + 4*A*C^2*a + B^4 + B^3*b + 2*B^2*C^2 + 4*B^2*C*g
         + 4*B^2*b^2 + B*C^2*b + C^4 + 4*C^3*g + 4*C^2*g^2 + 2*a*b*g

        sage: (q^2-q^-2)*A + q*B*C - q^-1*C*B == (q-q^-1)*a
        True
        sage: (q^2-q^-2)*B + q*C*A - q^-1*A*C == (q-q^-1)*b
        True
        sage: (q^2-q^-2)*C + q*A*B - q^-1*B*A == (q-q^-1)*g
        True
        sage: all(x * Omega == Omega * x for x in [A,B,C])
        True

    REFERENCES:

    - [Terwilliger2011]_
    """
    @staticmethod
    def __classcall_private__(cls, R, q=None):
        r"""
        Normalize input to ensure a unique representation.

        TESTS::

            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: AW1 = algebras.AskeyWilson(QQ)
            sage: AW2 = algebras.AskeyWilson(R, q)
            sage: AW1 is AW2
            True

            sage: AW = algebras.AskeyWilson(ZZ, 0)
            Traceback (most recent call last):
            ...
            ValueError: q cannot be 0

            sage: AW = algebras.AskeyWilson(ZZ, 3)
            Traceback (most recent call last):
            ...
            ValueError: q=3 is not invertible in Integer Ring
        """
        if q is None:
            R = LaurentPolynomialRing(R, 'q')
            q = R.gen()
        else:
            q = R(q)
        if q == 0:
            raise ValueError("q cannot be 0")
        if 1/q not in R:
            raise ValueError("q={} is not invertible in {}".format(q, R))
        if R not in Rings().Commutative():
            raise ValueError("{} is not a commutative ring".format(R))
        return super(AskeyWilsonAlgebra, cls).__classcall__(cls, R, q)

    def __init__(self, R, q):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: TestSuite(AW).run()  # long time
        """
        self._q = q
        cat = Algebras(Rings().Commutative()).WithBasis()
        indices = cartesian_product([NonNegativeIntegers()]*6)
        CombinatorialFreeModule.__init__(self, R, indices, prefix='AW',
                                         sorting_key=_basis_key,
                                         sorting_reverse=True,
                                         category=cat)
        self._assign_names('A,B,C,a,b,g')

    def _repr_term(self, t):
        r"""
        Return a string representation of the basis element indexed by ``t``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW._repr_term((0,0,0,0,0,0))
            '1'
            sage: AW._repr_term((5,1,2,3,7,2))
            'A^5*B*C^2*a^3*b^7*g^2'
            sage: AW._repr_term((0,1,0,3,7,2))
            'B*a^3*b^7*g^2'
        """
        def exp(l, e):
            if e == 0:
                return ''
            if e == 1:
                return '*' + l
            return '*' + l + '^{}'.format(e)
        ret = ''.join(exp(l, e) for l, e in zip(['A','B','C','a','b','g'], t))
        if not ret:
            return '1'
        if ret[0] == '*':
            ret = ret[1:]
        return ret

    def _latex_term(self, t):
        r"""
        Return a latex representation of the basis element indexed by ``t``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW._latex_term((0,0,0,0,0,0))
            '1'
            sage: AW._latex_term((5,1,2,3,7,2))
            'A^{5}BC^{2}\\alpha^{3}\\beta^{7}\\gamma^{2}'
            sage: AW._latex_term((0,1,0,3,7,2))
            'B\\alpha^{3}\\beta^{7}\\gamma^{2}'
        """
        if sum(t) == 0:
            return '1'
        def exp(l, e):
            if e == 0:
                return ''
            if e == 1:
                return l
            return l + '^{{{}}}'.format(e)
        var_names = ['A', 'B', 'C', '\\alpha', '\\beta', '\\gamma']
        return ''.join(exp(l, e) for l, e in zip(var_names, t))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.AskeyWilson(QQ)
            Askey-Wilon algebra with q=q over
             Univariate Laurent Polynomial Ring in q over Rational Field
        """
        return "Askey-Wilon algebra with q={} over {}".format(self._q, self.base_ring())

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: G = AW.algebra_generators()
            sage: G['A']
            A
            sage: G['a']
            a
            sage: list(G)
            [A, B, C, a, b, g]
        """
        A = self.variable_names()
        def build_monomial(g):
            exp = [0] * 6
            exp[A.index(g)] = 1
            return self.monomial(self._indices(exp))
        return Family(A, build_monomial)

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.gens()
            (A, B, C, a, b, g)
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element `1` of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.one_basis()
            (0, 0, 0, 0, 0, 0)
        """
        return self._indices([0]*6)

    def q(self):
        r"""
        Return the parameter `q` of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: q = AW.q()
            sage: q
            q
            sage: q.parent()
            Univariate Laurent Polynomial Ring in q over Rational Field
        """
        return self._q

    @cached_method
    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.an_element()
            (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A
        """
        q = self._q
        I = self._indices
        R = self.base_ring()
        elt = {I((1,0,0,0,0,0)): R(1),
               I((1,0,2,0,1,0)): R.an_element(),
               I((0,1,0,2,0,1)): q**2 * R(3),
               I((0,0,0,1,1,3)): q**-3 + R(3) + R(2)*q + q**2}
        return self.element_class(self, elt)

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.some_elements()
            (A, B, C, a, b, g, 1,
             (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A,
             q*A*B*C + q^2*A^2 - q*A*a + (q^-2)*B^2 - (q^-1)*B*b + q^2*C^2 - q*C*g)
        """
        return self.gens() + (self.one(), self.an_element(), self.casimir_element())

    @cached_method
    def casimir_element(self):
        r"""
        Return the Casimir element of ``self``.

        The Casimir element of the Askey-Wilson algebra `\Delta_q` is

        .. MATH::

            \Omega = q ABC + q^2 A^2 + q^{-2} B^2 + q^2 C^2
                     - q A\alpha - q^{-1} B\beta - q C\gamma.

        The center `Z(\Delta_q)` is generated by `\alpha`, `\beta`,
        `\gamma`, and `\Omega`.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.casimir_element()
            q*A*B*C + q^2*A^2 - q*A*a + (q^-2)*B^2 - (q^-1)*B*b + q^2*C^2 - q*C*g

        We check that the Casimir element is in the center::

            sage: Omega = AW.casimir_element()
            sage: all(Omega * gen == gen * Omega for gen in AW.algebra_generators())
            True
        """
        q = self._q
        I = self._indices
        d = {I((1,1,1,0,0,0)): q,      # q ABC
             I((2,0,0,0,0,0)): q**2,   # q^2 A^2
             I((0,2,0,0,0,0)): q**-2,  # q^-2 B^2
             I((0,0,2,0,0,0)): q**2,   # q^2 C^2
             I((1,0,0,1,0,0)): -q,     # -q A\alpha
             I((0,1,0,0,1,0)): -q**-1, # -q^-1 B\beta
             I((0,0,1,0,0,1)): -q}     # -q C\gamma
        return self.element_class(self, d)

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of the basis elements indexed by ``x`` and ``y``.

        INPUT:

        - ``x``, ``y`` -- tuple of length 6

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.product_on_basis((0,0,0,0,0,0), (3,5,2,0,12,3))
            A^3*B^5*C^2*b^12*g^3
            sage: AW.product_on_basis((0,0,0,5,3,5), (3,5,2,0,12,3))
            A^3*B^5*C^2*a^5*b^15*g^8
            sage: AW.product_on_basis((7,0,0,5,3,5), (0,5,2,0,12,3))
            A^7*B^5*C^2*a^5*b^15*g^8
            sage: AW.product_on_basis((7,3,0,5,3,5), (0,2,2,0,12,3))
            A^7*B^5*C^2*a^5*b^15*g^8
            sage: AW.product_on_basis((0,1,0,5,3,5), (2,0,0,0,5,3))
            q^4*A^2*B*a^5*b^8*g^8 - (q^-3-q^5)*A*C*a^5*b^8*g^8
             + (1-q^4)*A*a^5*b^8*g^9 - (q^-4-2+q^4)*B*a^5*b^8*g^8
             + (q^-3-q^-1-q+q^3)*a^5*b^9*g^8
            sage: AW.product_on_basis((0,2,1,0,2,0), (1,1,0,2,1,0))
            q^4*A*B^3*C*a^2*b^3 - (q^5-q^9)*A^2*B^2*a^2*b^3
             + (q^2-q^4)*A*B^2*a^3*b^3 + (q^-3-q)*B^4*a^2*b^3
             - (q^-2-1)*B^3*a^2*b^4 - (q-q^9)*B^2*C^2*a^2*b^3
             + (1-q^4)*B^2*C*a^2*b^3*g + (q^-4+2-5*q^4+2*q^12)*A*B*C*a^2*b^3
             - (q^-1+q-2*q^3-2*q^5+q^7+q^9)*A*B*a^2*b^3*g
             - (q^-3-q^3-2*q^5+q^7+q^9)*B*C*a^3*b^3
             + (q^-2-1-q^2+q^4)*B*a^3*b^3*g
             - (q^-3-2*q+2*q^9-q^13)*A^2*a^2*b^3
             + (2*q^-2-2-3*q^2+3*q^4+q^10-q^12)*A*a^3*b^3
             + (q^-7-2*q^-3+2*q^5-q^9)*B^2*a^2*b^3
             - (q^-6-q^-4-q^-2+1-q^2+q^4+q^6-q^8)*B*a^2*b^4
             - (q^-7-q^-3-2*q+2*q^5+q^9-q^13)*C^2*a^2*b^3
             + (q^-6-3-2*q^2+5*q^4-q^8+q^10-q^12)*C*a^2*b^3*g
             - (q^-1-2*q+2*q^5-q^7)*a^4*b^3
             - (q^-3-q^-1-2*q+2*q^3+q^5-q^7)*a^2*b^3*g^2
        """
        I = self._indices
        # Commute the central parts to the right
        lhs = list(x[:3])
        rhs = list(y)
        for i in range(3,6):
            rhs[i] += x[i]

        # No ABC variables on the RHS to move past
        if sum(rhs[:3]) == 0:
            return self.monomial(I(lhs + rhs[3:]))

        # We recurse using the PBW-type basis property:
        #   that YX = XY + lower order terms (see Theorem 4.1 in Terwilliger).
        q = self._q
        if lhs[2] > 0: # lhs has a C
            if rhs[0] > 0: # rhs has an A to commute with C
                lhs[2] -= 1
                rhs[0] -= 1
                rel = {I((1,0,1,0,0,0)): q**-2,         # q^2 AC
                       I((0,1,0,0,0,0)): q**-3 - q**1,  # q^-1(q^-2-q^2) B
                       I((0,0,0,0,1,0)): 1 - q**-2}     # -q^-1(q^-1-q) b
                rel = self.element_class(self, rel)
                return self.monomial(I(lhs+[0]*3)) * (rel * self.monomial(I(rhs)))
            elif rhs[1] > 0: # rhs has a B to commute with C
                lhs[2] -= 1
                rhs[1] -= 1
                rel = {I((0,1,1,0,0,0)): q**2,          # q^2 BC
                       I((1,0,0,0,0,0)): q**3 - q**-1,  # q(q^2-q^-2) A
                       I((0,0,0,1,0,0)): -q**2 + 1}     # -q(q-q^-1) a
                rel = self.element_class(self, rel)
                return self.monomial(I(lhs+[0]*3)) * (rel * self.monomial(I(rhs)))
            else: # nothing to commute as rhs has no A nor B
                rhs[2] += lhs[2]
                rhs[1] = lhs[1]
                rhs[0] = lhs[0]
                return self.monomial(I(rhs))

        elif lhs[1] > 0: # lhs has a B
            if rhs[0] > 0: # rhs has an A to commute with B
                lhs[1] -= 1
                rhs[0] -= 1
                rel = {I((1,1,0,0,0,0)): q**2,          # q^2 AB
                       I((0,0,1,0,0,0)): q**3 - q**-1,  # q(q^2-q^-2) C
                       I((0,0,0,0,0,1)): -q**2 + 1}     # -q(q-q^-1) g
                rel = self.element_class(self, rel)
                return self.monomial(I(lhs+[0]*3)) * (rel * self.monomial(I(rhs)))
            else: # nothing to commute as rhs has no A
                rhs[1] += lhs[1]
                rhs[0] = lhs[0]
                return self.monomial(I(rhs))

        elif lhs[0] > 0: # lhs has an A
            rhs[0] += lhs[0]
            return self.monomial(I(rhs))

        # otherwise, lhs is just 1
        return self.monomial(I(rhs))

    def permutation_automorphism(self):
        r"""
        Return the permutation automorphism `\rho` of ``self``.

        We define the automorphism `\rho` by

        .. MATH::

            A \mapsto B \mapsto C \mapsto A,
            \qquad
            \alpha \mapsto \beta \mapsto \gamma \mapsto \alpha.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: rho = AW.permutation_automorphism()
            sage: [rho(gen) for gen in AW.algebra_generators()]
            [B, C, A, b, g, a]

            sage: AW.an_element()
            (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A
            sage: rho(AW.an_element())
            (q^-3+3+2*q+q^2)*a^3*b*g + q^5*A^2*B*g + 3*q^2*C*a*b^2
             - (q^-2-q^6)*A*C*g + (q-q^5)*A*g^2 - (q^-3-2*q+q^5)*B*g
             + (q^-2-1-q^2+q^4)*b*g + B

            sage: r3 = rho * rho * rho
            sage: [r3(gen) for gen in AW.algebra_generators()]
            [A, B, C, a, b, g]
            sage: r3(AW.an_element()) == AW.an_element()
            True
        """
        A,B,C,a,b,g = self.gens()
        return AlgebraMorphism(self, [B,C,A,b,g,a], codomain=self)

    rho = permutation_automorphism

    def reflection_automorphism(self):
        r"""
        Return the reflection automorphism `\sigma` of ``self``.

        We define the automorphism `\sigma` by

        .. MATH::

            A \mapsto B \mapsto A,
            \qquad
            C \mapsto C + \frac{AB - BA}{q-q^{-1}}
            = C - qAB - (1+q^2) C + q \gamma,

        .. MATH::

            \alpha \mapsto \beta \mapsto \alpha,
            \gamma \mapsto \gamma.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: sigma = AW.reflection_automorphism()
            sage: [sigma(gen) for gen in AW.algebra_generators()]
            [B, A, -q*A*B - q^2*C + q*g, b, a, g]

            sage: AW.an_element()
            (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A
            sage: sigma(AW.an_element())
            q^9*A^2*B^3*a + (q^10+q^14)*A*B^2*C*a - (q^7+q^9)*A*B^2*a*g
             + (q^-3+3+2*q+q^2)*a*b*g^3 + (q-3*q^9+q^13+q^17)*A^2*B*a
             - (q^2-q^6-q^8+q^14)*A*B*a^2 + 3*q^2*A*b^2*g + (q^5-q^9)*B^3*a
             - (q^6-q^8)*B^2*a*b + q^13*B*C^2*a - 2*q^10*B*C*a*g + q^7*B*a*g^2
             + (q^2-2*q^10+q^18)*A*C*a - (q-q^7-2*q^9+2*q^11-q^15+q^17)*A*a*g
             - (q^3-q^7-q^9+q^13)*C*a^2 + (q^2-q^6-2*q^8+2*q^10)*a^2*g
             + (q-3*q^5+3*q^9-q^13)*B*a - (q^2-q^4-2*q^6+2*q^8+q^10-q^12)*a*b + B

            sage: s2 = sigma * sigma
            sage: [s2(gen) for gen in AW.algebra_generators()]
            [A, B, C, a, b, g]
            sage: s2(AW.an_element()) == AW.an_element()
            True
        """
        A,B,C,a,b,g = self.gens()
        q = self._q
        # Note that sage: (A*B-B*A) / (q-q^-1) == -q*A*B - (1+q^2)*C + q*g
        Cp = C - q*A*B - (1+q**2)*C + q*g
        return AlgebraMorphism(self, [B,A,Cp,b,a,g], codomain=self)

    sigma = reflection_automorphism

    def loop_representation(self):
        r"""
        Return the map `\pi` from ``self`` to `2 \times 2` matrices
        over `R[\lambda,\lambda^{-1}]`, where `F` is the fraction field
        of the base ring of ``self``.

        Let `AW` be the Askey-Wilson algebra over `R`, and let `F` be
        the fraction field of `R`. Let `M` be the space of `2 \times 2`
        matrices over `F[\lambda, \lambda^{-1}]`. Consider the following
        elements of `M`:

        .. MATH::

            \mathcal{A} = \begin{pmatrix}
                \lambda & 1 - \lambda^{-1} \\ 0 & \lambda^{-1}
            \end{pmatrix},
            \qquad
            \mathcal{B} = \begin{pmatrix}
                \lambda^{-1} & 0 \\ \lambda - 1 & \lambda
            \end{pmatrix},
            \qquad
            \mathcal{C} = \begin{pmatrix}
                1 & \lambda - 1 \\ 1 - \lambda^{-1} & \lambda + \lambda^{-1} - 1
            \end{pmatrix}.

        From Lemma 3.11 of [Terwilliger2011]_, we define a
        representation `\pi: AW \to M` by

        .. MATH::

            A \mapsto q \mathcal{A} + q^{-1} \mathcal{A}^{-1},
            \qquad
            B \mapsto q \mathcal{B} + q^{-1} \mathcal{B}^{-1},
            \qquad
            C \mapsto q \mathcal{C} + q^{-1} \mathcal{C}^{-1},

        .. MATH::

            \alpha, \beta, \gamma \mapsto \nu I,

        where `\nu = (q^2 + q^-2)(\lambda + \lambda^{-1})
        + (\lambda + \lambda^{-1})^2`.

        We call this representation the *loop representation* as
        it is a representation using the loop group
        `SL_2(F[\lambda,\lambda^{-1}])`.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: q = AW.q()
            sage: pi = AW.loop_representation()
            sage: A,B,C,a,b,g = [pi(gen) for gen in AW.algebra_generators()]
            sage: A
            [                1/q*lambda^-1 + q*lambda ((-q^2 + 1)/q)*lambda^-1 + ((q^2 - 1)/q)]
            [                                       0                 q*lambda^-1 + 1/q*lambda]
            sage: B
            [             q*lambda^-1 + 1/q*lambda                                     0]
            [((-q^2 + 1)/q) + ((q^2 - 1)/q)*lambda              1/q*lambda^-1 + q*lambda]
            sage: C
            [1/q*lambda^-1 + ((q^2 - 1)/q) + 1/q*lambda      ((q^2 - 1)/q) + ((-q^2 + 1)/q)*lambda]
            [  ((q^2 - 1)/q)*lambda^-1 + ((-q^2 + 1)/q)    q*lambda^-1 + ((-q^2 + 1)/q) + q*lambda]
            sage: a
            [lambda^-2 + ((q^4 + 1)/q^2)*lambda^-1 + 2 + ((q^4 + 1)/q^2)*lambda + lambda^2                                                                             0]
            [                                                                            0 lambda^-2 + ((q^4 + 1)/q^2)*lambda^-1 + 2 + ((q^4 + 1)/q^2)*lambda + lambda^2]
            sage: a == b
            True
            sage: a == g
            True

            sage: AW.an_element()
            (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A
            sage: x = pi(AW.an_element())
            sage: y = (q^-3+3+2*q+q^2)*a*b*g^3 + q*A*C^2*b + 3*q^2*B*a^2*g + A
            sage: x == y
            True

        We check the defining relations of the Askey-Wilson algebra::

            sage: A + (q*B*C - q^-1*C*B) / (q^2 - q^-2) == a / (q + q^-1)
            True
            sage: B + (q*C*A - q^-1*A*C) / (q^2 - q^-2) == b / (q + q^-1)
            True
            sage: C + (q*A*B - q^-1*B*A) / (q^2 - q^-2) == g / (q + q^-1)
            True

        We check Lemma 3.12 in [Terwilliger2011]_::

            sage: M = pi.codomain()
            sage: la = M.base_ring().gen()
            sage: p = M([[0,-1],[1,1]])
            sage: s = M([[0,1],[la,0]])
            sage: rho = AW.rho()
            sage: sigma = AW.sigma()
            sage: all(p*pi(gen)*~p == pi(rho(gen)) for gen in AW.algebra_generators())
            True
            sage: all(s*pi(gen)*~s == pi(sigma(gen)) for gen in AW.algebra_generators())
            True
        """
        from sage.matrix.matrix_space import MatrixSpace
        q = self._q
        base = LaurentPolynomialRing(self.base_ring().fraction_field(), 'lambda')
        la = base.gen()
        inv = ~la
        M = MatrixSpace(base, 2)
        A = M([[la,1-inv],[0,inv]])
        Ai = M([[inv,inv-1],[0,la]])
        B = M([[inv,0],[la-1,la]])
        Bi = M([[la,0],[1-la,inv]])
        C = M([[1,1-la],[inv-1,la+inv-1]])
        Ci = M([[la+inv-1,la-1],[1-inv,1]])
        mu = la + inv
        nu = (self._q**2 + self._q**-2) * mu + mu**2
        nuI = M(nu)
        # After #29374 is fixed, the category can become
        # Algebras(Rings().Commutative()) as it was before #29399.
        category = Rings()
        return AlgebraMorphism(self, [q*A + q**-1*Ai, q*B + q**-1*Bi, q*C + q**-1*Ci,
                                      nuI, nuI, nuI],
                               codomain=M, category=category)

    pi = loop_representation

def _basis_key(t):
    """
    Return a key for the basis element of the Askey-Wilson algebra
    indexed by ``t``.

    EXAMPLES::

        sage: from sage.algebras.askey_wilson import _basis_key
        sage: I = algebras.AskeyWilson(QQ).indices()
        sage: _basis_key(I((0,2,3,1,2,5)))
        (13, (0, 2, 3, 1, 2, 5))
    """
    return (sum(t), t.value)

class AlgebraMorphism(ModuleMorphismByLinearity):
    """
    An algebra morphism of the Askey-Wilson algebra defined by
    the images of the generators.
    """
    def __init__(self, domain, on_generators, position=0, codomain=None,
                 category=None):
        """
        Given a map on the multiplicative basis of a free algebra, this method
        returns the algebra morphism that is the linear extension of its image
        on generators.

        INPUT:

        - ``domain`` -- an Askey-Wilson algebra
        - ``on_generators`` -- a list of length 6 corresponding to
          the images of the generators
        - ``codomain`` -- (optional) the codomain
        - ``position`` -- (default: 0) integer
        - ``category`` -- (optional) category

        OUTPUT:

        - module morphism

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: sigma = AW.sigma()
            sage: TestSuite(sigma).run()
        """
        if category is None:
            category = Algebras(Rings().Commutative()).WithBasis()
        self._on_generators = tuple(on_generators)
        ModuleMorphismByLinearity.__init__(self, domain=domain, codomain=codomain,
                                           position=position, category=category)

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.algebras.askey_wilson import AlgebraMorphism
            sage: AW = algebras.AskeyWilson(QQ)
            sage: rho = AW.rho()
            sage: sigma = AW.sigma()
            sage: id = AlgebraMorphism(AW, AW.gens(), codomain=AW)
            sage: sigma * sigma == id
            True
            sage: id == rho * rho * rho
            True
        """
        return (self.__class__ is other.__class__ and self.parent() == other.parent()
                and self._zero == other._zero
                and self._on_generators == other._on_generators
                and self._position == other._position
                and self._is_module_with_basis_over_same_base_ring
                    == other._is_module_with_basis_over_same_base_ring)

    def _on_basis(self, c):
        r"""
        Computes the image of this morphism on the basis element
        indexed by ``c``.

        INPUT:

        - ``c`` -- a tuple of length 6

        OUTPUT:

        - element of the codomain

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: AW.inject_variables()
            Defining A, B, C, a, b, g
            sage: rho = AW.rho()
            sage: sigma = AW.sigma()
            sage: rho._on_basis((2,1,1,2,0,5))
            q^2*A*B^2*C*a^5*b^2 + (q^-3-q)*B^3*a^5*b^2 - (q^-2-1)*B^2*a^5*b^3
             - (q^-3-q^5)*B*C^2*a^5*b^2 + (q^-2-q^2)*B*C*a^5*b^2*g
             + (q^-2-2*q^2+q^6)*A*C*a^5*b^2 - (q^-1-q-q^3+q^5)*C*a^6*b^2

            sage: sigma._on_basis((2,1,1,2,0,5))
            -q^9*A^2*B^3*b^2*g^5 + (q^2-q^10-q^14)*A*B^2*C*b^2*g^5
             - (q^3-q^7-q^9)*A*B^2*b^2*g^6
             - (2*q-q^5-3*q^9+q^13+q^17)*A^2*B*b^2*g^5
             + (1-q^6-q^8+q^14)*A*B*a*b^2*g^5 + (q^-3-q-q^5+q^9)*B^3*b^2*g^5
             - (q^-2-1-q^6+q^8)*B^2*b^3*g^5 + (q^5-q^13)*B*C^2*b^2*g^5
             - (q^2+q^6-2*q^10)*B*C*b^2*g^6 + (q^3-q^7)*B*b^2*g^7
             + (q^-6-4*q^2+2*q^6+2*q^10-q^18)*A*C*b^2*g^5
             - (q^-3+q^-1-3*q-2*q^3+q^5+2*q^7+2*q^9-2*q^11+q^15-q^17)*A*b^2*g^6
             - (q^-3-2*q-q^3+q^5+q^7+q^9-q^13)*C*a*b^2*g^5
             + (q^-2-1-q^2-q^4+2*q^6+2*q^8-2*q^10)*a*b^2*g^6
             + (q^-7-3*q^-3+2*q+2*q^5-3*q^9+q^13)*B*b^2*g^5
             - (q^-6-q^-4-2*q^-2+2+2*q^6-2*q^8-q^10+q^12)*b^3*g^5

            sage: rho(B*A)
            q^2*B*C - (q^-1-q^3)*A + (1-q^2)*a
            sage: rho(A*B)
            B*C
            sage: rho(A*B*C)
            A*B*C + (q^-3-q)*B^2 - (q^-2-1)*B*b - (q^-3-q)*C^2 + (q^-2-1)*C*g
            sage: rho(B*C*A)
            A*B*C - (q^-3-q)*A^2 + (q^-2-1)*A*a + (q^-3-q)*B^2 - (q^-2-1)*B*b
            sage: rho(C*A*B)
            A*B*C

            sage: rho(C^2*a*b^6*g^2)
            A^2*a^2*b*g^6

            sage: sigma(C^2)
            q^4*A^2*B^2 + (q^3+q^7)*A*B*C - (q^2+q^4)*A*B*g
             - (q^4-q^8)*A^2 + (q^5-q^7)*A*a + (1-q^4)*B^2
             - (q-q^3)*B*b + q^4*C^2 - 2*q^3*C*g + q^2*g^2
            sage: sigma(A*B)
            q^2*A*B - (q^-1-q^3)*C + (1-q^2)*g
            sage: sigma(C + 3*g*A*B)
            3*q^2*A*B*g - q*A*B - (3*q^-1-3*q^3)*C*g
             + (3-3*q^2)*g^2 - q^2*C + q*g
        """
        return self.codomain().prod(self._on_generators[i]**exp
                                    for i,exp in enumerate(c))

    def _composition_(self, right, homset):
        """
        Return the composition of ``self`` and ``right`` in ``homset``.

        EXAMPLES::

            sage: AW = algebras.AskeyWilson(QQ)
            sage: rho = AW.rho()
            sage: sigma = AW.sigma()
            sage: s2 = sigma * sigma
            sage: s2._on_generators
            (A, B, C, a, b, g)
            sage: sr = sigma * rho
            sage: sr._on_generators
            (C, B, -q*B*C - q^2*A + q*a, g, b, a)
            sage: rs = rho * sigma
            sage: rs._on_generators
            (A, -q*A*B - q^2*C + q*g, B, a, g, b)
        """
        if isinstance(right, AlgebraMorphism):
            cat = homset.homset_category()
            return AlgebraMorphism(homset.domain(),
                                   [right(g) for g in self._on_generators],
                                   codomain=homset.codomain(), category=cat)
        return super(self, AlgebraMorphism)._composition_(right, homset)

