r"""
Powersum Basis of Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.partition import _Partitions, Partition
from sage.categories.tensor import tensor
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.misc_c import prod
from sage.misc.lazy_attribute import lazy_attribute

class SupersymFunctionAlgebra_powersum(super_sfa.SuperSymAlgebra_multiplicative):
    r"""
    Powersum supersymmetric functions.

    The powersum supersymmetric function defined on variables `\mathbb{x}` and
    `\mathbb{y}`, `p_i(\mathbb{x} \mid \mathbb{y})`, is given as

    .. MATH::

        p_i(\mathbb{x} \mid \mathbb{y}) = \sum_{k=1}^\infty x_k^i - y_k^i

    These form a non-graded multiplicative basis for the ring of supersymmetric
    functions.

    REFERENCES:

    - [BHS25]_

    INPUT:

    - ``Supersym`` -- ring of supersymmetric functions

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: p = s.p()
        sage: p
        Supersymmetric functions over Rational Field in the powersum basis
    """
    def __init__(self, Supersym):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: TestSuite(p).run()
        """
        super_sfa.SuperSymAlgebra_generic.__init__(self, SuperSym=Supersym, graded=True,
                                                   prefix='p', basis_name='powersum')

    def coproduct_on_generators(self, i):
        r"""
        Return coproduct on generators for power sums `p_i`
        (for integer `i`).

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p.coproduct_on_generators(3)
            p[] # p[3] + p[3] # p[]
        """
        Pi = _Partitions([i])
        P0 = _Partitions([])
        T = self.tensor_square()
        return T.sum_of_monomials([(Pi, P0), (P0, Pi)])

    def antipode_on_basis(self, part):
        r"""
        Return the antipode of ``self[part]``.

        INPUT:

        - ``part`` -- a list or partition

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: from sage.combinat.partition import _Partitions
            sage: part = _Partitions([4,4,2])
            sage: p.antipode_on_basis(part)
            -p[4, 4, 2]
        """
        if not isinstance(part, Partition):
            part = sorted(part, reverse=True)
            part = _Partitions(part)
        if len(part) % 2 == 0:
            return self[part]
        return -self[part]

    def _lift_on_basis(self, part):
        r"""
        Return the value of ``self[part]`` under the plethysm.

        INPUT:

        - ``part`` -- partition

        OUTPUT:

        - value under the plethysm in the tensor product of powersum symmetric functions

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: from sage.combinat.partition import _Partitions
            sage: part = _Partitions([4,4,2])
            sage: p._lift_on_basis(part)
            p[] # p[4, 4, 2] + p[2] # p[4, 4] + 2*p[4] # p[4, 2] +
            2*p[4, 2] # p[4] + p[4, 4] # p[2] + p[4, 4, 2] # p[]
        """
        R = self.base_ring()
        p_sq = SymmetricFunctions(R).p().tensor_square()
        if not part:
            return p_sq.one()
        B = p_sq.basis()
        empty = _Partitions([])
        return p_sq.prod([B[_Partitions([p]), empty] + (-1)**p * B[empty, _Partitions([p])]
                          for p in part])

    @lazy_attribute
    def lift(self):
        r"""
        Return the plethysm of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p.lift
            Generic morphism:
            From: Supersymmetric functions over Rational Field in the powersum basis
            To:   Symmetric Functions over Rational Field in the powersum basis # Symmetric Functions over Rational Field in the powersum basis
        """
        R = self.base_ring()
        p_sym = SymmetricFunctions(R).p()
        T = p_sym.tensor_square()
        return self.module_morphism(self._lift_on_basis, triangular='upper', codomain=T, inverse_on_support=lambda x: x[0])

    @lazy_attribute
    def retract(self):
        r"""
        Return retract of plethysm of ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: p = s.p()
            sage: p.retract
            Generic morphism:
            From: Symmetric Functions over Rational Field in the powersum basis # Symmetric Functions over Rational Field in the powersum basis
            To:   Supersymmetric functions over Rational Field in the powersum basis
        """
        return self.lift.section()

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
                sage: p = s.p()
                sage: p[[4]].expand(3)
                x1^4 + x2^4 + x3^4 - y1^4 - y2^4 - y3^4
            """
            self_parts = self.monomial_coefficients()
            x_gens = [alphabet_x + str(i) for i in range(1, n + 1)]
            y_gens = [alphabet_y + str(i) for i in range(1, n + 1)]
            R = PolynomialRing(self.base_ring(), x_gens + y_gens)
            R_gens = R.gens()
            x_gens1 = R_gens[:n]
            y_gens1 = R_gens[n:]
            return sum(self_parts[part] *
                       prod(sum(x**p - y**p
                                for x, y in zip(x_gens1, y_gens1))
                            for p in part)
                       for part in self_parts)
