r"""
Monomial Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.categories.tensor import tensor

class SupersymFunctionAlgebra_monomial(super_sfa.SuperSymAlgebra_generic):
    r"""
    Monomial supersymmetric functions.

    The *monomial supersymmetric function* defined on variables `\mathbf{x}` and
    `\mathbf{y}`, `m_\lambda(\mathbf{x} \mid \mathbf{y})`, is given in terms of
    monomial symmetric functions and forgotten symmetric functions:

    .. MATH::

       m_\lambda(\mathbf{x} \mid \mathbf{y}) = \sum_{\mu \cup \nu = \lambda} m_\mu(\mathbf{x}) f_\nu(\mathbf{y})

    These form a non-graded multiplicative basis for the ring of supersymmetric
    functions.

    REFERENCES:

    - [Moe2007]_

    INPUT:

    - ``Supersym`` -- ring of supersymmetric functions

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: m = s.m()
        sage: m
        Supersymmetric functions over Rational Field in the monomial basis
    """
    def __init__(self, SuperSym):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: m = s.m()
            sage: TestSuite(m).run()
        """
        super().__init__(SuperSym=SuperSym, graded=False, prefix='m', basis_name='monomial')

    def coproduct_by_coercion(self, elt):
        r"""
        Return the coproduct of ``elt`` by coercing into powersum
        supersymmetric basis and back.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: el = m[1,1,1]
            sage: elt = m.coproduct_by_coercion(el)
            sage: elt
            m[] # m[1, 1, 1] + m[1] # m[1, 1] + m[1, 1] # m[1] + m[1, 1, 1] # m[]
        """
        p = self.realization_of().p()
        return self.tensor_square().sum(coeff * tensor([self(p[x]), self(p[y])])
                                        for (x, y), coeff in p(elt).coproduct())

    def antipode_by_coercion(self, el):
        r"""
        Return the antipode of ``part`` by coercing into powersum
        supersymmetric basis and back.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: elt = m[2,1]
            sage: m.antipode_by_coercion(elt)
            m[2, 1] + 2*m[3]
        """
        p = self.realization_of().p()
        elt = p(el)
        pow_el = p.antipode(elt)
        return self(pow_el)

    def lift_on_basis(self, el):
        r"""
        Return the monomial basis in terms of tensor product of powersum symmetric functions.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: Sym = SuperSymmetricFunctions(QQ)
            sage: m = Sym.m()
            sage: el = m[1,1,1]
            sage: m.lift_on_basis(el)
            -1/6*p[] # p[1, 1, 1] + 1/2*p[] # p[2, 1] - 1/3*p[] # p[3] +
            1/2*p[1] # p[1, 1] - 1/2*p[1] # p[2] - 1/2*p[1, 1] # p[1] +
            1/6*p[1, 1, 1] # p[] + 1/2*p[2] # p[1] - 1/2*p[2, 1] # p[] +
            1/3*p[3] # p[]
        """
        p = self.realization_of().p()
        pow_el = p(el)
        res = p.lift(pow_el)
        return res
