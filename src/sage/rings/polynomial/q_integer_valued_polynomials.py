r"""
Quantum-valued polynomial rings

AUTHORS:

- Frédéric Chapoton (2024-03): Initial version
"""
# ***************************************************************************
#  Copyright (C) 2024 Frédéric Chapoton <chapoton-math-unistra-fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from sage.arith.misc import binomial
from sage.categories.algebras import Algebras
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.q_analogues import q_binomial, q_int
from sage.data_structures.blas_dict import linear_combination
from sage.matrix.constructor import matrix
from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


def q_int_x(n):
    """
    Return the interpolating polynomial of `q`-integers.

    INPUT:

    - ``n`` -- a positive integer

    EXAMPLES::

        sage: from sage.rings.polynomial.q_integer_valued_polynomials import q_int_x
        sage: q_int_x(3)
        q^2*x + q + 1
    """
    ring_q = PolynomialRing(ZZ, 'q')
    q = ring_q.gen()
    x = polygen(ring_q, 'x')
    return q_int(n - 1) + q**(n - 1) * x


def q_binomial_x(m, n):
    r"""
    Return a `q`-analogue of ``binomial(m + x, n)``.

    When evaluated at the `q`-integer `[k]_q`, this gives
    the usual `q`-binomial coefficient `[m + k, n]_q`.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_int
        sage: from sage.rings.polynomial.q_integer_valued_polynomials import q_binomial_x, q_int_x
        sage: q_binomial_x(4,2)(0) == q_binomial(4,2)
        True
        sage: q_binomial_x(3,2)(1) == q_binomial(4,2)
        True
        sage: q_binomial_x(3,1) == q_int_x(4)
        True
        sage: q_binomial_x(2,0).parent()
        Univariate Polynomial Ring in x over Fraction Field of
        Univariate Polynomial Ring in q over Integer Ring
    """
    ring = PolynomialRing(PolynomialRing(ZZ, 'q').fraction_field(), 'x')
    if n == 0:
        return ring.one()
    return ring.prod(q_int_x(m + 2 - i) / q_int(i) for i in range(1, n + 1))


class QuantumValuedPolynomialRing(UniqueRepresentation, Parent):
    r"""
    The quantum-valued polynomial ring on some generators over a base ring.

    Quantum-valued polynomial rings are commutative and associative
    algebras, with a basis indexed by integers.

    The basis used here is given by `B[i] = \binom{i+n}{i}` for `i \in \NN`.

    There is a nice formula for the product, see [HaHo2017]_.

    INPUT:

    - ``R`` -- commutative ring

    REFERENCES:

    - [HaHo2017] Nate Harman and Sam Hopkins, *Quantum integer-valued
      polynomials*, J. Alg. Comb. 2017, :doi:`10.1007/s10801-016-0717-3`

    EXAMPLES::

        sage: F = QuantumValuedPolynomialRing(QQ).S(); F
        Quantum-Valued Polynomial Ring over Rational Field
        in the shifted basis

        sage: F.gen()
        S[1]

        sage: S = QuantumValuedPolynomialRing(ZZ); S
        Quantum-Valued Polynomial Ring over Integer Ring
        sage: S.base_ring()
        Univariate Laurent Polynomial Ring in q over Integer Ring

    Quantum-valued polynomial rings commute with their base ring::

        sage: K = QuantumValuedPolynomialRing(QQ).S()
        sage: a = K.gen()
        sage: c = K.monomial(2)

    Quantum-valued polynomial rings are commutative::

        sage: c^3 * a == c * a * c * c
        True

    We can also manipulate elements in the basis and coerce elements from our
    base field::

        sage: F = QuantumValuedPolynomialRing(QQ).S()
        sage: B = F.basis()
        sage: B[2] * B[3]
        (q^-5+q^-4+q^-3)*S[3] - (q^-6+2*q^-5+3*q^-4+3*q^-3+2*q^-2+q^-1)*S[4]
        + (q^-6+q^-5+2*q^-4+2*q^-3+2*q^-2+q^-1+1)*S[5]
        sage: 1 - B[2] * B[2] / 2
        S[0] - (1/2*q^-3)*S[2] + (1/2*q^-4+q^-3+q^-2+1/2*q^-1)*S[3]
        - (1/2*q^-4+1/2*q^-3+q^-2+1/2*q^-1+1/2)*S[4]
    """
    def __init__(self, R) -> None:
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = QuantumValuedPolynomialRing(QQ); F
            Quantum-Valued Polynomial Ring over Rational Field
            sage: TestSuite(F).run()  # not tested

        TESTS::

            sage: QuantumValuedPolynomialRing(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a commutative ring
        """
        if R not in Rings().Commutative():
            msg = "argument R must be a commutative ring"
            raise TypeError(msg)
        laurent = LaurentPolynomialRing(R, 'q')
        self._ground_ring = R
        cat = Algebras(laurent).Commutative().WithBasis()
        Parent.__init__(self, base=laurent, category=cat.WithRealizations())

    _shorthands = ["B", "S"]

    def _repr_(self) -> str:
        r"""
        Return the string representation.

        EXAMPLES::

            sage: F = QuantumValuedPolynomialRing(QQ)
            sage: F  # indirect doctest
            Quantum-Valued Polynomial Ring over Rational Field

            sage: QuantumValuedPolynomialRing(ZZ)
            Quantum-Valued Polynomial Ring over Integer Ring
        """
        base = self.base_ring().base_ring()
        return f"Quantum-Valued Polynomial Ring over {base}"

    def a_realization(self):
        """
        Return a default realization.

        The Shifted realization is chosen.

        EXAMPLES::

            sage: QuantumValuedPolynomialRing(QQ).a_realization()
            Quantum-Valued Polynomial Ring over Rational Field
            in the shifted basis
        """
        return self.Shifted()

    class Bases(Category_realization_of_parent):
        def super_categories(self) -> list:
            r"""
            Return the super-categories of ``self``.

            EXAMPLES::

                sage: A = QuantumValuedPolynomialRing(QQ); A
                Quantum-Valued Polynomial Ring over Rational Field
                sage: C = A.Bases(); C
                Category of bases of Quantum-Valued Polynomial Ring
                over Rational Field
                sage: C.super_categories()
                [Category of realizations of Quantum-Valued Polynomial Ring
                 over Rational Field,
                 Join of Category of algebras with basis
                 over Univariate Laurent Polynomial Ring in q over Rational Field and
                 Category of filtered algebras
                 over Univariate Laurent Polynomial Ring in q over Rational Field and
                 Category of commutative algebras
                 over Univariate Laurent Polynomial Ring in q over Rational Field and
                 Category of realizations of unital magmas]
            """
            A = self.base()
            category = Algebras(A.base_ring()).Commutative().Filtered()
            return [A.Realizations(),
                    category.Realizations().WithBasis()]

        class ParentMethods:
            def ground_ring(self):
                """
                Return the ring of coefficients.

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(QQ).S()
                    sage: A.ground_ring()
                    Rational Field
                """
                return self.realization_of()._ground_ring

            def _repr_(self) -> str:
                r"""
                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(QQ).S()
                    sage: F  # indirect doctest
                    Quantum-Valued Polynomial Ring over Rational Field
                    in the shifted basis
                """
                real = self.realization_of()
                return f"{real} in the {self._realization_name()} basis"

            @cached_method
            def one_basis(self):
                r"""
                Return the number 0, which index the unit of this algebra.

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(QQ).S()
                    sage: A.one_basis()
                    0
                    sage: A.one()
                    S[0]
                """
                return self.basis().keys()(0)

            def degree_on_basis(self, m):
                r"""
                Return the degree of the basis element indexed by ``m``.

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(QQ).S()
                    sage: A.degree_on_basis(4)
                    4
                """
                return ZZ(m)

            def from_polynomial(self, p):
                r"""
                Convert a polynomial into the ring of quantum-valued polynomials.

                This raises a :exc:`ValueError` if this is not possible.

                INPUT:

                - ``p`` -- a polynomial in ``x`` with coefficients in ``QQ(q)``

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).S()
                    sage: S = A.basis()
                    sage: A.from_polynomial((S[1]).polynomial())
                    S[1]
                    sage: A.from_polynomial((S[2]+2*S[3]).polynomial())
                    S[2] + 2*S[3]

                    sage: A = QuantumValuedPolynomialRing(ZZ).B()
                    sage: B = A.basis()
                    sage: A.from_polynomial((B[1]).polynomial())
                    B[1]
                    sage: A.from_polynomial((B[2]+2*B[3]).polynomial())
                    B[2] + 2*B[3]
                """
                B = self.basis()
                poly = self._poly
                laurent_polys = self.base_ring()
                remain = p.change_variable_name('x')
                result = self.zero()
                while remain:
                    N = remain.degree()
                    base_N = poly(N)
                    top_coeff = remain.leading_coefficient() / base_N.lc()
                    denom = top_coeff.denominator()
                    if denom.is_term():
                        numer = top_coeff.numerator()
                        top_coeff_laurent = laurent_polys(numer) / laurent_polys(denom)
                    else:
                        msg = 'not a polynomial with integer values :'
                        msg += f' {top_coeff} is not a Laurent polynomial'
                        raise ValueError(msg)
                    remain += -top_coeff * base_N
                    result += top_coeff_laurent * B[N]
                return result

            def gen(self, i=0):
                r"""
                Return the generator of the algebra.

                The optional argument is ignored.

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: F.gen()
                    S[1]
                """
                return self.basis()[1]

            @cached_method
            def algebra_generators(self):
                r"""
                Return the generators of this algebra.

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).S(); A
                    Quantum-Valued Polynomial Ring over Integer Ring
                    in the shifted basis
                    sage: A.algebra_generators()
                    Family (S[1],)
                """
                return Family([self.basis()[1]])

            gens = algebra_generators

        class ElementMethods:
            def __call__(self, v):
                """
                Return the evaluation at some value ``v``.

                EXAMPLES::

                     sage: F = QuantumValuedPolynomialRing(ZZ).S()
                     sage: B = F.gen()
                     sage: f = B**2+4*B+6
                     sage: f(1/3)
                     (q^2 + 18*q + 99)/9

                     sage: F = QuantumValuedPolynomialRing(ZZ).B()
                     sage: B = F.gen()
                     sage: f = F.monomial(2)+4*B+6
                     sage: f(1/3)
                     (66*q^2 + 66*q - 2)/(9*q^2 + 9*q)
                """
                return self.polynomial()(v)

            def polynomial(self):
                """
                Convert to a polynomial in `x`.

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: S = F.gen()
                    sage: (S+1).polynomial()
                    q*x + 2

                    sage: F = QuantumValuedPolynomialRing(ZZ).B()
                    sage: B = F.gen()
                    sage: (B+1).polynomial()
                    x + 1

                TESTS::

                    sage: F.zero().polynomial().parent()
                    Univariate Polynomial Ring in x over Fraction Field
                    of Univariate Polynomial Ring in q over Integer Ring
                """
                ring = self.parent().ground_ring()
                fractions = PolynomialRing(ring, 'q').fraction_field()
                R = PolynomialRing(fractions, 'x')
                p = self.parent()._poly
                return R.sum(c * p(i) for i, c in self)

            def shift(self, j=1):
                """
                Shift all indices by `j`.

                INPUT:

                - `j` -- integer (default 1)

                In the binomial basis, the shift by 1 corresponds to
                a summation operator from `0` to `x`.

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: B = F.gen()
                    sage: (B+1).shift()
                    S[1] + S[2]
                """
                A = self.parent()
                return A._from_dict({A._indices(i + j): c for i, c in self})

            def sum_of_coefficients(self):
                """
                Return the sum of coefficients.

                In the shifted basis, this is the evaluation at `x=0`.

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: B = F.basis()
                    sage: (B[2]*B[4]).sum_of_coefficients()
                    1
                """
                R = self.parent().base_ring()
                return R.sum(self._monomial_coefficients.values())

    class Shifted(CombinatorialFreeModule, BindableClass):
        r"""
        The quantum-valued polynomial ring in the shifted basis.

        The basis used here is given by `S[i] = \genfrac{[}{]}{0pt}{}{i+x}{i}_q` for `i \in \NN`.

        Assuming `n_1 \leq n_2`, the product of two monomials `S[n_1] \cdot S[n_2]`
        is given by the sum

        .. MATH::

            \sum_{k=0}^{n_1} (-1)^k q^{\binom{k}{2} - n_1 * n_2} \genfrac{[}{]}{0pt}{}{n_1}{k}_q \genfrac{[}{]}{0pt}{}{n_1+n_2-k}{n_1}_q S[n_1 + n_2 - k].
        """
        def __init__(self, A):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(QQ).S(); F
                Quantum-Valued Polynomial Ring over Rational Field
                in the shifted basis
                sage: TestSuite(F).run()  # not tested
            """
            CombinatorialFreeModule.__init__(self, A.base_ring(),
                                             NonNegativeIntegers(),
                                             category=A.Bases(),
                                             prefix="S",
                                             latex_prefix=r"\mathbb{S}")

        def _an_element_(self):
            """
            Return a small element of ``self``.

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(QQ).S()
                sage: F.an_element()
                2*S[0] + 4*S[2]
            """
            NonNeg = self.basis().keys()
            ring = self.base_ring()
            return self.element_class(self, {NonNeg(0): ring(2),
                                             NonNeg(2): ring(4)})

        def _realization_name(self) -> str:
            r"""
            TESTS::

                sage: F = QuantumValuedPolynomialRing(QQ).S()
                sage: F._realization_name()
                'shifted'
            """
            return "shifted"

        def product_on_basis(self, n1, n2):
            r"""
            Return the product of basis elements ``n1`` and ``n2``.

            INPUT:

            - ``n1``, ``n2`` -- integers

            EXAMPLES::

                sage: A = QuantumValuedPolynomialRing(QQ).S()
                sage: A.product_on_basis(0, 1)
                S[1]
            """
            i = ZZ(n1)
            j = ZZ(n2)
            if j < i:
                j, i = i, j
            q = self.base_ring().gen()
            return self._from_dict({i + j - k: (-1)**k
                                    * q_binomial(i, k)
                                    * q_binomial(i + j - k, i)
                                    * q**(binomial(k, 2) - i * j)
                                    for k in range(i + 1)})

        def _from_binomial_basis(self, i):
            """
            Convert from the ``binomial(x,k)`` basis.

            INPUT:

            - ``i`` -- an integer

            EXAMPLES::

                sage: S = QuantumValuedPolynomialRing(ZZ).S()
                sage: B = QuantumValuedPolynomialRing(ZZ).B()
                sage: b = B.basis()
                sage: S(b[3]+1)  # indirect doctest
                -(q^-6-1)*S[0] + (q^-8+q^-7+q^-6)*S[1]
                - (q^-9+q^-8+q^-7)*S[2] + (q^-9)*S[3]
                sage: B(_)
                B[0] + B[3]
            """
            i = ZZ(i)
            R = self.base_ring()
            q = self.base_ring().gen()
            return self._from_dict({k: R((-1)**(i - k) * q_binomial(i, k))
                                    * q**(-i**2 + binomial(i - k, 2))
                                    for k in range(i + 1)})

        def from_h_vector(self, hv):
            """
            Convert from some `h`-vector.

            EXAMPLES::

                sage: A = QuantumValuedPolynomialRing(ZZ).S()
                sage: B = A.basis()
                sage: ex = B[2] + B[3]
                sage: A.from_h_vector(ex.h_vector())
                S[2] + S[3]

                sage: q = A.base_ring().gen()
                sage: ex = B[2] + q*B[3]
                sage: A.from_h_vector(ex.h_vector())
                S[2] + q*S[3]
            """
            B = self.basis()
            ring = self.base_ring()
            q = ring.gen()
            d = len(hv) - 1
            m = matrix(ring, d + 1, d + 1,
                       lambda j, i: (-1)**(d - j) * q_binomial(d - i, d - j, q) *
                       q**(-d * (d - i) + binomial(d - j, 2)))
            v = vector(ring, [hv[i] for i in range(d + 1)])
            return sum(ring(c) * B[i] for i, c in enumerate(m * v))

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            INPUT:

            - ``x`` -- an element of the base ring or something convertible

            EXAMPLES::

                sage: R = QuantumValuedPolynomialRing(QQ).S()
                sage: x = R.gen()
                sage: R(3)
                3*S[0]
                sage: R(x)
                S[1]
            """
            P = x.parent()
            if isinstance(P, QuantumValuedPolynomialRing.Shifted):
                if P is self:
                    return x
                if P is not self.base_ring():
                    return self.element_class(self, x.monomial_coefficients())

            # ok, not a quantum-valued polynomial ring element
            R = self.base_ring()
            # coercion via base ring
            x = R(x)
            if x == 0:
                return self.element_class(self, {})
            return self.from_base_ring_from_one_basis(x)

        def _coerce_map_from_(self, R) -> bool:
            r"""
            Return whether there is a coercion from ``R`` into ``self``.

            INPUT:

            - ``R`` -- a commutative ring

            The things that coerce into ``self`` are

            - Quantum-Valued Polynomial Rings over a base
              with a coercion map into ``self.base_ring()``.

            - Anything with a coercion into ``self.base_ring()``.

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(GF(7)).S(); F
                Quantum-Valued Polynomial Ring over Finite Field of size 7
                in the shifted basis

            Elements of the quantum-valued polynomial ring canonically coerce in::

                sage: x = F.gen()
                sage: F.coerce(x*x) # indirect doctest
                (6*q^-1)*S[1] + (q^-1+1)*S[2]

            Elements of the integers coerce in, since there is a coerce map
            from `\ZZ` to GF(7)::

                sage: F.coerce(1)       # indirect doctest
                S[0]

            There is no coerce map from `\QQ` to `\GF{7}`::

                sage: F.coerce(2/3)  # indirect doctest
                Traceback (most recent call last):
                ...
                TypeError: no canonical coercion from Rational Field to
                Quantum-Valued Polynomial Ring over Finite Field of size 7
                in the shifted basis

            Elements of the base ring coerce in::

                sage: F.coerce(GF(7)(5))
                5*S[0]

            The quantum-valued polynomial ring over `\ZZ` on `x` coerces in, since
            `\ZZ` coerces to `\GF{7}`::

                sage: G = QuantumValuedPolynomialRing(ZZ).S()
                sage: Gx = G.gen()
                sage: z = F.coerce(Gx**2); z
                -(q^-1)*S[1] + (q^-1+1)*S[2]
                sage: z.parent() is F
                True

            However, `\GF{7}` does not coerce to `\ZZ`, so the shuffle
            algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

                sage: G.coerce(x^3+x)
                Traceback (most recent call last):
                ...
                TypeError: no canonical coercion from Quantum-Valued Polynomial
                Ring over Finite Field of size 7 in the shifted basis
                to Quantum-Valued Polynomial Ring over Integer Ring
                in the shifted basis

            TESTS::

                sage: F = QuantumValuedPolynomialRing(ZZ).S()
                sage: G = QuantumValuedPolynomialRing(QQ).S()
                sage: H = QuantumValuedPolynomialRing(ZZ).S()
                sage: F._coerce_map_from_(G)
                False
                sage: G._coerce_map_from_(F)
                True
                sage: F._coerce_map_from_(H)
                True
                sage: F._coerce_map_from_(QQ)
                False
                sage: G._coerce_map_from_(QQ)
                True
                sage: F.has_coerce_map_from(PolynomialRing(ZZ, 'x'))
                False
            """
            # quantum-valued polynomial rings in the same variable
            # over any base that coerces in:
            if isinstance(R, QuantumValuedPolynomialRing.Shifted):
                return self.base_ring().has_coerce_map_from(R.base_ring())
            if isinstance(R, QuantumValuedPolynomialRing.Binomial):
                return R.module_morphism(self._from_binomial_basis,
                                         codomain=self)
            return self.base_ring().has_coerce_map_from(R)

        def _poly(self, i):
            """
            Convert the basis element `S[i]` to a polynomial.

            INPUT:

            - ``i`` -- an integer

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(ZZ).S()
                sage: F._poly(4).factor()
                (1/(q^6 + 3*q^5 + 5*q^4 + 6*q^3 + 5*q^2 + 3*q + 1)) *
                (q*x + 1) * (q^2*x + q + 1) * (q^3*x + q^2 + q + 1) *
                (q^4*x + q^3 + q^2 + q + 1)
            """
            return q_binomial_x(i, i)

        class Element(CombinatorialFreeModule.Element):

            def umbra(self):
                """
                Return the Bernoulli umbra.

                This is the derivative at `-1` of the shift by one.

                This is related to Carlitz's `q`-Bernoulli numbers.

                .. SEEALSO:: :meth:`derivative_at_minus_one`

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: B = F.gen()
                    sage: (B+1).umbra()
                    (q + 2)/(q + 1)
                """
                return self.shift().derivative_at_minus_one()

            def variable_shift(self, k=1):
                r"""
                Return the image by the shift on variables.

                The shift is the substitution operator

                .. MATH::

                    x \mapsto q x + 1.

                INPUT:

                - `k` -- integer (default: 1)

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).S()
                    sage: S = A.basis()
                    sage: S[5].variable_shift()
                    S[0] + q*S[1] + q^2*S[2] + q^3*S[3] + q^4*S[4] + q^5*S[5]

                    sage: S[5].variable_shift(-1)
                    -(q^-5)*S[4] + (q^-5)*S[5]

                TESTS::

                    sage: S[5].variable_shift(0)
                    S[5]
                    sage: S[5].variable_shift().variable_shift(-1)
                    S[5]
                    sage: S[5].variable_shift(2).variable_shift(-2)
                    S[5]
                    sage: S[3].variable_shift(-2)
                    (q^-5)*S[1] - (q^-6+q^-5)*S[2] + (q^-6)*S[3]
                """
                if k == 0:
                    return self

                A = self.parent()
                q = A.base_ring().gen()

                def on_basis(n):
                    return {A._indices(j): q**(k * j)
                            * q_binomial(k + n - 1 - j, n - j)
                            for j in range(n + 1)}

                mc = self._monomial_coefficients
                ret = linear_combination((on_basis(index), coeff)
                                         for index, coeff in mc.items())
                return A.element_class(A, ret)

            def derivative_at_minus_one(self):
                """
                Return the 'derivative' at -1.

                .. SEEALSO:: :meth:`umbra`

                EXAMPLES::

                    sage: F = QuantumValuedPolynomialRing(ZZ).S()
                    sage: B = F.gen()
                    sage: (B+1).derivative_at_minus_one()
                    1
                """
                ring = q_int(1).parent()
                return ring.sum(c / q_int(i) for i, c in self if i > 0)

            def h_vector(self):
                """
                Return the numerator of the generating series of values.

                If ``self`` is an Ehrhart polynomial, this is the h-vector.

                .. SEEALSO:: :meth:`h_polynomial`, :meth:`fraction`

                changement de base vers les (binomial(x+i,d))_{i=0..d}

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).S()
                    sage: ex = A.basis()[4]
                    sage: ex.h_vector()
                    (0, 0, 0, 0, 1)

                    sage: q = polygen(QQ,'q')
                    sage: x = polygen(q.parent(),'x')
                    sage: ex = A.from_polynomial((1+q*x)**3)
                    sage: ex.h_vector()
                    (0, q^3, 2*q + 2*q^2, 1)
                """
                d = max(self.support())
                ring = self.parent().base_ring()
                q = ring.gen()

                def fn(j, i):
                    return ((-1)**(d - j) *
                            q**(binomial(d - j + i + 1, 2) -
                                binomial(i + 1, 2)) *
                            q_binomial(d - i, d - j))
                m = matrix(ring, d + 1, d + 1, fn)
                v = vector(ring, [self.coefficient(i) for i in range(d + 1)])
                return m * v

            def h_polynomial(self):
                """
                Return the `h`-vector as a polynomial.

                .. SEEALSO:: :meth:`h_vector`, :meth:`fraction`

                peut-etre pas dans le bon sens ?

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).S()
                    sage: q = polygen(ZZ,'q')
                    sage: x = polygen(q.parent(),'x')
                    sage: ex = A.from_polynomial((1+q*x)**3)
                    sage: ex.h_polynomial()
                    z^3 + (2*q + 2*q^2)*z^2 + q^3*z
                """
                ring = PolynomialRing(self.parent().base_ring(), 'z')
                return ring(list(self.h_vector()))

            def fraction(self):
                """
                Return the generating series of values as a fraction.

                .. SEEALSO:: :meth:`h_vector`, :meth:`h_polynomial`

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(QQ).S()
                    sage: ex = A.basis()[4]
                    sage: ex.fraction().factor()
                    (-1) * (t - 1)^-1 * (q*t - 1)^-1 * (q^2*t - 1)^-1 * (q^3*t - 1)^-1 * (q^4*t - 1)^-1

                    sage: q = polygen(QQ,'q')
                    sage: x = polygen(q.parent(), 'x')
                    sage: ex = A.from_polynomial((1+q*x)**3)
                    sage: ex.fraction().factor()
                    (t - 1)^-1 * (q*t - 1)^-1 * (q^2*t - 1)^-1 * (q^3*t - 1)^-1 * (q^3*t^2 + 2*q^2*t + 2*q*t + 1)
                    sage: ex.fraction().numerator()
                    q^3*t^2 + 2*q^2*t + 2*q*t + 1
                """
                v = self.h_vector()
                d = len(v)
                R = PolynomialRing(QQ, 'q,t')
                frac_R = R.fraction_field()
                q, t = R.gens()
                denom = R.prod(1 - q**i * t for i in range(d))
                numer = sum(frac_R(v[i]) * t**(d - 1 - i) for i in range(d))
                return numer / denom

    S = Shifted

    # =====     Another basis for the same algebra     =====

    class Binomial(CombinatorialFreeModule, BindableClass):
        r"""
        The quantum-valued polynomial ring in the binomial basis.

        The basis used here is given by `B[i] = \genfrac{[}{]}{0pt}{}{x}{i}_q` for `i \in \NN`.

        Assuming `n_1 \leq n_2`, the product of two monomials `B[n_1] \cdot B[n_2]`
        is given by the sum

        .. MATH::

            \sum_{k=0}^{n_1} q^{(k-n_1)(k-n_2)} \genfrac{[}{]}{0pt}{}{n_1}{k}_q \genfrac{[}{]}{0pt}{}{n_1+n_2-k}{n_1}_q B[n_1 + n_2 - k].
        """
        def __init__(self, A) -> None:
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(QQ).B(); F
                Quantum-Valued Polynomial Ring over Rational Field
                in the binomial basis
                sage: TestSuite(F).run()  # not tested
            """
            CombinatorialFreeModule.__init__(self, A.base_ring(),
                                             NonNegativeIntegers(),
                                             category=A.Bases(),
                                             prefix="B",
                                             latex_prefix=r"\mathbb{B}")

        def _realization_name(self) -> str:
            r"""
            TESTS::

                sage: F = QuantumValuedPolynomialRing(QQ).B()
                sage: F._realization_name()
                'binomial'
            """
            return "binomial"

        def product_on_basis(self, n1, n2):
            r"""
            Return the product of basis elements ``n1`` and ``n2``.

            INPUT:

            - ``n1``, ``n2`` -- integers

            The formula is taken from Theorem 3.4 in Harman-Hopkins.

            EXAMPLES::

                sage: A = QuantumValuedPolynomialRing(QQ).B()
                sage: A.product_on_basis(0, 1)
                B[1]
            """
            i = ZZ(n1)
            j = ZZ(n2)
            if j < i:
                j, i = i, j

            q = self.base_ring().gen()
            return self._from_dict({i + j - k:
                                    q_binomial(i, k)
                                    * q_binomial(i + j - k, i)
                                    * q**((k - i) * (k - j))
                                    for k in range(i + 1)})

        def _from_shifted_basis(self, i):
            """
            Convert from the shifted binomial(x+k,k) basis.

            INPUT:

            - ``i`` -- an integer

            EXAMPLES::

                sage: S = QuantumValuedPolynomialRing(ZZ).S()
                sage: B = QuantumValuedPolynomialRing(ZZ).B()
                sage: s = S.basis()
                sage: B(s[3]+1)  # indirect doctest
                2*B[0] + (q+q^2+q^3)*B[1] + (q^4+q^5+q^6)*B[2] + q^9*B[3]
                sage: S(_)
                S[0] + S[3]
            """
            i = ZZ(i)
            R = self.base_ring()
            q = self.base_ring().gen()
            return self._from_dict({k: R(q_binomial(i, k))
                                    * q**(k**2)
                                    for k in range(i + 1)})

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            EXAMPLES::

                sage: R = QuantumValuedPolynomialRing(QQ).B()
                sage: x = R.gen()
                sage: R(3)
                3*B[0]
                sage: R(x)
                B[1]
            """
            P = x.parent()
            if isinstance(P, QuantumValuedPolynomialRing.Binomial):
                if P is self:
                    return x
                if P is not self.base_ring():
                    return self.element_class(self, x.monomial_coefficients())

            # ok, not a quantum-valued polynomial ring element
            R = self.base_ring()
            # coercion via base ring
            x = R(x)
            if x == 0:
                return self.element_class(self, {})
            return self.from_base_ring_from_one_basis(x)

        def _coerce_map_from_(self, R):
            r"""
            Return whether there is a coercion from ``R`` into ``self``.

            INPUT:

            - ``R`` -- a commutative ring

            The things that coerce into ``self`` are

            - Quantum-Valued Polynomial Rings over a base
              with a coercion map into ``self.base_ring()``.

            - Anything with a coercion into ``self.base_ring()``.

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(GF(7)).B(); F
                Quantum-Valued Polynomial Ring over Finite Field of size 7
                in the binomial basis

            Elements of the integer-valued polynomial ring canonically coerce
            in::

                sage: x = F.gen()
                sage: F.coerce(x*x) # indirect doctest
                B[1] + (q+q^2)*B[2]

            Elements of the integers coerce in, since there is a coerce map
            from `\ZZ` to `\GF(7)`::

                sage: F.coerce(1)       # indirect doctest
                B[0]

            There is no coerce map from `\QQ` to `\GF{7}`::

                sage: F.coerce(2/3)  # indirect doctest
                Traceback (most recent call last):
                ...
                TypeError: no canonical coercion from Rational Field to
                Quantum-Valued Polynomial Ring over Finite Field of size 7
                in the binomial basis

            Elements of the base ring coerce in::

                sage: F.coerce(GF(7)(5))
                5*B[0]

            The integer-valued polynomial ring over `\ZZ` on `x` coerces in,
            since `\ZZ` coerces to `\GF{7}`::

                sage: G = QuantumValuedPolynomialRing(ZZ).B()
                sage: Gx = G.gen()
                sage: z = F.coerce(Gx**2); z
                B[1] + (q+q^2)*B[2]
                sage: z.parent() is F
                True

            However, `\GF{7}` does not coerce to `\ZZ`, so the
            integer-valued polynomial algebra over `\GF{7}` does not
            coerce to the one over `\ZZ`::

                sage: G.coerce(x^3+x)
                Traceback (most recent call last):
                ...
                TypeError: no canonical coercion from
                Quantum-Valued Polynomial Ring over Finite Field of size 7
                in the binomial basis to Quantum-Valued Polynomial Ring
                over Integer Ring in the binomial basis

            TESTS::

                sage: F = QuantumValuedPolynomialRing(ZZ).B()
                sage: G = QuantumValuedPolynomialRing(QQ).B()
                sage: H = QuantumValuedPolynomialRing(ZZ).B()
                sage: F._coerce_map_from_(G)
                False
                sage: G._coerce_map_from_(F)
                True
                sage: F._coerce_map_from_(H)
                True
                sage: F._coerce_map_from_(QQ)
                False
                sage: G._coerce_map_from_(QQ)
                True
                sage: F.has_coerce_map_from(PolynomialRing(ZZ,'x'))
                False
            """
            # quantum-valued polynomial rings over any base
            # that coerces in:
            if isinstance(R, QuantumValuedPolynomialRing.Binomial):
                return self.base_ring().has_coerce_map_from(R.base_ring())
            if isinstance(R, QuantumValuedPolynomialRing.Shifted):
                return R.module_morphism(self._from_shifted_basis,
                                         codomain=self)
            return self.base_ring().has_coerce_map_from(R)

        def _poly(self, i):
            """
            Convert the basis element `B[i]` to a polynomial.

            INPUT:

            - ``i`` -- an integer

            EXAMPLES::

                sage: F = QuantumValuedPolynomialRing(ZZ).B()
                sage: F._poly(4).factor()
                (1/(q^12 + 3*q^11 + 5*q^10 + 6*q^9 + 5*q^8 + 3*q^7 + q^6)) *
                (x - 1) * x * (x - q - 1) * (x - q^2 - q - 1)
            """
            return q_binomial_x(0, i)

        class Element(CombinatorialFreeModule.Element):
            def variable_shift(self, k=1):
                r"""
                Return the image by the shift of variables.

                On polynomials, the action for `k=1` is the shift
                on variables `x \mapsto 1 + qx`.

                This implementation follows formula (5.5) in [HaHo2017]_.

                INPUT:

                - `k` -- nonnegative integer (default: 1)

                EXAMPLES::

                    sage: A = QuantumValuedPolynomialRing(ZZ).B()
                    sage: B = A.basis()
                    sage: B[5].variable_shift()
                    B[4] + q^5*B[5]

                TESTS::

                    sage: B[5].variable_shift(0)
                    B[5]
                """
                if k == 0:
                    return self

                A = self.parent()
                q = A.base_ring().gen()

                def on_basis(n):
                    return {A._indices(j): q**((k + j - n) * j)
                            * q_binomial(k, n - j)
                            for j in range(n + 1)}

                mc = self._monomial_coefficients
                ret = linear_combination((on_basis(index), coeff)
                                         for index, coeff in mc.items())
                return A.element_class(A, ret)

    B = Binomial
