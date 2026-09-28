r"""
Homogeneous and Elementary Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from itertools import combinations_with_replacement, combinations
from sage.combinat.partition import Partitions, _Partitions
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.functions.other import factorial
from sage.misc.misc_c import prod

class SupersymFunctionAlgebra_hom_el(super_sfa.SuperSymAlgebra_multiplicative):
    r"""
    Homogeneous and elementary supersymmetric functions.

    The *homogeneous supersymmetric function* on variables `\mathbf{x}` and
    `\mathbf{y}` is defined as

    .. MATH::

        h_k(\mathbf{x} \mid \mathbf{y}) = \sum_{a+b=k}\sum_{\substack{i_1 \geq \cdots \geq i_b \\ j_1 < \cdots < j_a}}
        y_{j_1} \cdots y_{j_a} x_{i_1} \cdots x_{i_b}.

    The *elementary supersymmetric function* is defined as

    .. MATH::

        e_k(\mathbf{x} \mid \mathbf{y}) = \sum_{a+b=k}\sum_{\substack{i_1 > \cdots > i_b \\ j_1 \leq \cdots \leq j_a}}
        y_{j_1} \cdots y_{j_a} x_{i_1} \cdots x_{i_b}.

    These form multiplicative non-graded bases for the ring of supersymmetric
    functions.

    REFERENCES:

    - [BHS25]_

    INPUT:

    - ``Supersym`` -- the ring of supersymmetric functions
    - ``basis_name`` -- string (default: ``'homogeneous'``); one of the following

      * ``'homogeneous'`` - homogeneous basis of supersymmetric functions
      * ``'elementary'`` - elementary basis of supersymmetric functions

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: h = s.h()
        sage: h
        Supersymmetric functions over Rational Field in the homogeneous basis
        sage: e = s.e()
        sage: e
        Supersymmetric functions over Rational Field in the elementary basis
    """
    def __init__(self, Supersym, basis_name='homogeneous'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: TestSuite(h).run()
        """
        prefix = ''
        if basis_name == 'homogeneous':
            prefix = 'h'
        elif basis_name == 'elementary':
            prefix = 'e'
        else:
            raise ValueError("invalid basis name")
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=True,
                                                   prefix=prefix, basis_name=basis_name)

    def antipode(self, x):
        r"""
        Return the antipode of ``x``.

        INPUT:

        - ``x`` -- element of ``self``

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: e = s.e()
            sage: h = s.h()
            sage: f = e[6,5]

            sage: # Tests
            sage: h.antipode(h[2])
            h[1, 1] - h[2]
            sage: e(h.antipode(h[2]))
            e[2]
            sage: all(e(h[la].antipode()) == e[la] for la in Partitions(6))
            True
        """
        if self.basis_name == 'homogeneous':
            e = self.realization_of().e()
            el = e(x)
        elif self.basis_name == 'elementary':
            h = self.realization_of().h()
            el = h(x)

        return self.sum_of_terms((lam, (-1)**(sum(lam) % 2) * a)
                                 for lam, a in el)

    def coproduct_on_generators(self, n):
        r"""
        Return the coproduct on `h_n`.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: h.coproduct_on_generators(5)
            h[] # h[5] + h[1] # h[4] + h[2] # h[3] + h[3] # h[2] + h[4] # h[1] + h[5] # h[]
        """
        def P(i):
            return _Partitions([i]) if i else _Partitions([])
        T = self.tensor_square()
        return T.sum_of_monomials((P(j), P(n-j)) for j in range(n+1))

    def lift_on_gens(self, n):
        r"""
        Return the homogeneous or elementary basis in terms of powersum
        basis for a given ``n``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: h = s.h()
            sage: h.lift_on_gens(3)
            p[1, 1, 1] + 1/4*p[2, 1] + 1/18*p[3]
        """
        basis_name = self.basis_name
        parts = Partitions(n)
        Supersym = self.realization_of()
        R = self.base_ring()
        res = R.zero()

        def z(part):
            prod = R.one()
            for p in part:
                m = part.to_exp()[p-1]
                prod *= (p ** m) * factorial(p)
            return prod

        ssp = Supersym.p()
        if basis_name == 'homogeneous':
            res = ssp._from_dict({part: ~z(part) for part in parts})
        elif basis_name == 'elementary':
            res = ssp._from_dict({part: (-1) ** (n - len(part)) / z(part) for part in parts})
        return res

    class Element(super_sfa.SuperSymAlgebra_multiplicative.Element):
        def expand(self, n, m, alphabet_x='x', alphabet_y='y'):
            r"""
            Expand the supersymmetric function ``self`` as a supersymmetric
            polynomial in ``n`` variables.

            INPUT:

            - ``n`` -- nonnegative integer
            - ``alphabet_x`` -- (default: ``'x'``) a variable for the expansion `x`
            - ``alphabet_y`` -- (default: ``'y'``) a variable for the expansion `y`

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: Sym = SuperSymmetricFunctions(QQ)
                sage: h = Sym.h()
                sage: h[2,1].expand(1,1)  # corner cases for homogeneous
                x0^3 + 2*x0^2*y0 + x0*y0^2
                sage: h[2,1].expand(0,1)
                0
                sage: h[2,1].expand(0,0)
                0
                sage: h[2,1].expand(1,0)
                x0^3

                sage: # Comparing with Sym
                sage: sym = SymmetricFunctions(QQ)
                sage: h1 = sym.h()
                sage: h[2,1].expand(2,0) == h1[2,1].expand(2)
                True

                sage: # Checking corner cases for elementary
                sage: e = Sym.e()
                sage: e[2,1].expand(1,1)
                x0^2*y0 + 2*x0*y0^2 + y0^3
                sage: e[2,1].expand(1,0)
                0
                sage: e[2,1].expand(0,0)
                0
                sage: e[2,1].expand(0,1)
                y0^3
                sage: # Comparing with sym
                sage: e1 = sym.e()
                sage: e[2,1].expand(2,0) == e1[2,1].expand(2)
                True
            """
            basis_name = self.parent().basis_name
            monomial_coeff = self.monomial_coefficients()
            x_gens = [alphabet_x + str(i) for i in range(n)]
            y_gens = [alphabet_y + str(i) for i in range(m)]
            variables = x_gens + y_gens
            R = PolynomialRing(self.base_ring(), variables)
            R_gens = R.gens_dict()
            x_gens = [R_gens[gen] for gen in x_gens]
            y_gens = [R_gens[gen] for gen in y_gens]

            def el_i(i, X_gens, Y_gens, x_count, y_count):
                req_sum = R.zero()
                for p in range(i+1):
                    for seq1 in combinations(range(y_count), (i - p)):
                        for seq2 in combinations_with_replacement(range(x_count), p):
                                res_prod = prod([Y_gens[i] for i in seq1]) * prod([X_gens[j] for j in seq2])
                                req_sum += res_prod
                return req_sum

            if basis_name == 'homogeneous':

                fin_res = R.zero()
                for part in monomial_coeff:
                    res_prod = R.one()
                    for p in part:
                        res_prod *= el_i(p, x_gens, y_gens, n, m)
                    fin_res += monomial_coeff[part] * res_prod
                return fin_res

            elif basis_name == 'elementary':

                fin_res = R.zero()
                for part in monomial_coeff:
                    res_prod = R.one()
                    for p in part:
                        res_prod *= el_i(p, y_gens, x_gens, m, n)
                        # Same function el_i used for hom and el,
                        # except with gens swapped as per definition
                    fin_res += monomial_coeff[part] * res_prod
                return fin_res

            else:
                raise ValueError("invalid basis name")