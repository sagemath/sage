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

    The homogeneous supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `h_k(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        h_k(\mathbb{x} \mid \mathbb{y}) = \sum_{a+b=k}\sum_{i_1 \geq \ldots \geq i_b \\ j_1 < \ldots < j_a}
        y_{j_1} \ldots y_{j_a} x_{i_1} \ldots x_{i_b}.

    The elementary supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `e_k(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        e_k(\mathbb{x} \mid \mathbb{y}) = \sum_{a+b=k}\sum_{i_1 > \ldots > i_b \\ j_1 \leq \ldots \leq j_a}
        y_{j_1} \ldots y_{j_a} x_{i_1} \ldots x_{i_b}.

    These form a multiplicative non-graded basis for the ring of supersymmetric
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
            raise ValueError("Invalid basis name")
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=True,
                                                   prefix=prefix, basis_name=basis_name)

    def antipode(self, x):
        r"""
        Return the antipode of ``x``.

        INPUT:

        - ``x`` -- element of ``self``

        OUTPUT:

        - the result of the antipode on ``x``

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: e = s.e()
            sage: f = e[6,5]
            sage: e.antipode(f)
            -e[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] +
            9*e[2, 1, 1, 1, 1, 1, 1, 1, 1, 1] -
            29*e[2, 2, 1, 1, 1, 1, 1, 1, 1] + 40*e[2, 2, 2, 1, 1, 1, 1, 1] -
            22*e[2, 2, 2, 2, 1, 1, 1] + 3*e[2, 2, 2, 2, 2, 1] -
            7*e[3, 1, 1, 1, 1, 1, 1, 1, 1] + 39*e[3, 2, 1, 1, 1, 1, 1, 1] -
            64*e[3, 2, 2, 1, 1, 1, 1] + 33*e[3, 2, 2, 2, 1, 1] -
            2*e[3, 2, 2, 2, 2] - 13*e[3, 3, 1, 1, 1, 1, 1] +
            30*e[3, 3, 2, 1, 1, 1] - 15*e[3, 3, 2, 2, 1] - 3*e[3, 3, 3, 1, 1] +
            2*e[3, 3, 3, 2] + 5*e[4, 1, 1, 1, 1, 1, 1, 1] -
            24*e[4, 2, 1, 1, 1, 1, 1] + 29*e[4, 2, 2, 1, 1, 1] -
            8*e[4, 2, 2, 2, 1] + 17*e[4, 3, 1, 1, 1, 1] -
            24*e[4, 3, 2, 1, 1] + 4*e[4, 3, 2, 2] + 2*e[4, 3, 3, 1] -
            6*e[4, 4, 1, 1, 1] + 4*e[4, 4, 2, 1] - 3*e[5, 1, 1, 1, 1, 1, 1] +
            13*e[5, 2, 1, 1, 1, 1] - 12*e[5, 2, 2, 1, 1] + e[5, 2, 2, 2] -
            10*e[5, 3, 1, 1, 1] + 10*e[5, 3, 2, 1] - e[5, 3, 3] +
            7*e[5, 4, 1, 1] - 2*e[5, 4, 2] - 2*e[5, 5, 1] +
            e[6, 1, 1, 1, 1, 1] - 4*e[6, 2, 1, 1, 1] + 3*e[6, 2, 2, 1] +
            3*e[6, 3, 1, 1] - 2*e[6, 3, 2] - 2*e[6, 4, 1] + e[6, 5]
            sage: h = s.h()
            sage: h.antipode(h[2])
            h[1, 1] - h[2]
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
        def expand(self, n, alphabet_x='x', alphabet_y='y'):
            r"""
            Expand the supersymmetric function ``self`` as a supersymmetric
            polynomial in ``n`` variables.

            INPUT:

            - ``n`` -- nonnegative integer
            - ``alphabet_x`` -- (default: ``'x'``) a variable for the expansion `x`
            - ``alphabet_y`` -- (default: ``'y'``) a variable for the expansion `y`

            EXAMPLES::

                sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
                sage: s = SuperSymmetricFunctions(QQ)
                sage: h = s.h()
                sage: h[3,2].expand(3)
                100*x1^2*x2^2*x3^2 + 360*x1^2*x2^2*x3*y1 + 324*x1^2*x2^2*y1^2 +
                180*x1^2*x2*x3*y1*y2 + 324*x1^2*x2*y1^2*y2 +
                81*x1^2*y1^2*y2^2 + 60*x1^2*x2^2*x3 + 108*x1^2*x2^2*y1 +
                90*x1^2*x2*x3*y1 + 162*x1^2*x2*y1^2 + 54*x1^2*x2*y1*y2 +
                81*x1^2*y1^2*y2
            """
            basis_name = self.parent().basis_name
            res = self.base_ring().one()
            monomial_coeff = self.monomial_coefficients()
            x_gens = [alphabet_x+str(i).format(i) for i in range(1, n+1)]
            y_gens = [alphabet_y+str(i).format(i) for i in range(1, n+1)]
            variables = x_gens + y_gens
            R = PolynomialRing(self.base_ring(), variables)
            R_gens = R.gens_dict()
            x_gens = [R_gens[gen] for gen in x_gens]
            y_gens = [R_gens[gen] for gen in y_gens]
            req_sum = R.zero()
            fin_res = R.zero()
            for k in monomial_coeff:
                for ki in k:
                    for p in range(1, ki+1):
                        for seq1 in combinations(range(n), (ki - p)):
                            for seq2 in combinations_with_replacement(range(n), p):
                                if basis_name == 'homogeneous':
                                    res_prod = prod([y_gens[i] for i in range(len(seq1))]) * prod([x_gens[j] for j in range(len(seq2))])
                                    req_sum += res_prod
                                elif basis_name == 'elementary':
                                    res_prod = prod([y_gens[i] for i in range(len(seq2))]) * prod([x_gens[j] for j in range(len(seq1))])
                                    req_sum += res_prod
                    res *= req_sum
                fin_res += monomial_coeff[k] * res
            return fin_res
