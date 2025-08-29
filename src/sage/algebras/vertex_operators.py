r"""
Vertex Operators

AUTHORS:

- Joseph McDonough (2025-08-04): Initial version

TESTS::

    sage: from sage.algebras.vertex_operators import *
    sage: R = SymmetricFunctions(QQ); p = R.p(); s = R.s()
    sage: P = p.completion()
    sage: F = FermionicFockSpace(s)
    sage: x = 3*F(([3,2,1],0)); x
    3*s[]*|[3, 2, 1], 0>
    sage: J = Current()
    sage: J.act(3, x)
    (-3*s[])*|[1, 1, 1], 0> + (-3*s[])*|[3], 0>
    sage: deg(x)
    6
"""

from sage.algebras.weyl_algebra import DifferentialWeylAlgebra
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.data_structures.stream import Stream_function, Stream_cauchy_compose
from sage.functions.other import factorial
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_function, cached_method
from sage.misc.flatten import flatten
from sage.misc.misc_c import prod
from sage.rings.infinity import PlusInfinity
from sage.rings.integer_ring import ZZ
from sage.rings.lazy_series_ring import LazyLaurentSeriesRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.element import MonoidElement
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class FermionicFockSpace(CombinatorialFreeModule):
    r"""
    A model of the Fermionic Fock space `\mathcal{H}_F`. As a vector space, `\mathcal{H}_F`
    has basis given by pairs `(\lambda, n)` where `\lambda` is a partition, and `n` is an integer.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: F = FermionicFockSpace(QQ)
        sage: F.an_element()
        |[], 0> + 4*|[], 1> + 2*|[1], 0>
        sage: TestSuite(F).run()
    """
    def __init__(self, R):
        r"""
        Initialize ``self``.
        """
        index_set = cartesian_product((Partitions(), ZZ))
        CombinatorialFreeModule.__init__(self, R, index_set, prefix='', bracket='')

    def _repr_term(self, m):
        return '|' + str(m[0]) + ', ' + str(m[1]) + '>'


def BosonicFockSpace(R=SymmetricFunctions(QQ).s(), names=('w',)):
    return LaurentPolynomialRing(R, names=names)


def BFMap():
    r"""
    Returns a map that compute the image of an element of the Fermionic Fock space
    under the Boson-Fermion correspondence. The basis element `|\lambda, n \rangle`
    gets mapped to `w^ns_\lambda`, where `w` is the charge variable and `s_\lambda`
    is the Schur function corresponding to `\lambda`.

    TESTS::

        sage: from sage.algebras.vertex_operators import *
        sage: F = FermionicFockSpace(QQ)
        sage: phi = BFMap()
        sage: phi(F(([2,1],-1)) + F(([],0)))
        s[2, 1]*w^-1 + s[]
    """
    F = FermionicFockSpace(QQ)
    B = BosonicFockSpace()
    w = B.gen()
    s = B.base_ring()
    return F.module_morphism(on_basis=lambda b: w**(b[1])*s(b[0]), codomain=B)


class Current():
    r"""
    The action of the current operators `J_i` on the Fermionic Fock space.
    For `i \neq 0`, we can view `J_i` as moving a particle by `-i` in all possible ways.
    Equivalently, viewed on the bosonic side, the action of `J_i` can be calculated
    using the Murnaghan-Nakayama rule.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: R = SymmetricFunctions(QQ); s = R.s()
        sage: F = FermionicFockSpace(s)
        sage: x = F(([3,2], 0))
        sage: J = Current()
        sage: J.act(-3, x)
        s[]*|[3, 2, 1, 1, 1], 0> + (-s[])*|[3, 2, 2, 1], 0> + (-s[])*|[4, 4], 0>
         + s[]*|[6, 2], 0>
        sage: J.act(1,x)
        s[]*|[2, 2], 0> + s[]*|[3, 1], 0>

    """
    def __init__(self):
        r"""
        Initialize ``self``.
        """
        self._R = SymmetricFunctions(QQ)
        self._s = self._R.s()
        self._p = self._R.p()
        self.fockspace = FermionicFockSpace(self._s)

    def act(self, g, x):
        res = self.fockspace.zero()
        if g > 0:  # \partial p_k
            for (m, c) in x.monomial_coefficients().items():
                par = self._s(m[0])
                for (m2, c2) in (par.skew_by(self._p[g])).monomial_coefficients().items():
                    res += c*c2*self.fockspace((m2, m[1]))
            return res
        elif g < 0:  # -k*p_k
            for (m, c) in x.monomial_coefficients().items():
                par = self._s(m[0])
                for (m2, c2) in (self._s(self._p[-g]*par)).monomial_coefficients().items():
                    res += c*c2*self.fockspace((m2, m[1]))
            return res
        else:  # multiply by charge
            return sum(m[1]*c*self.fockspace(m) for (m, c) in x.monomial_coefficients().items())

    def act_by_partition(self, lam, sign, x):
        for i in lam:
            # TODO: Make this more efficient by acting by the whole partition at once
            x = self.act(i*sign, x)
        return x


@cached_function
def Hamiltonian():
    r"""
    Returns the (exponentiated) Hamiltonian `\exp( \sum_{j \geq 1} p_i J_i)`
    where `J_i` is a current operator and `p_i` is a powersum symmetric function.

    This operator acts on the Fermionic Fock space, and its matrix coefficients are
    the skew Schur functions.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: H = Hamiltonian()
        sage: H.truncate(3)
        p[] + p[1] + (1/2*p[1,1]+1/2*p[2])
    """
    p = SymmetricFunctions(QQ).p()
    P = p.completion()
    return P(lambda n: p[n]/n, valuation=1).exp()


def act_by_H(x):
    r"""
    Computes the action of H on an element ``x`` of the Fermionic Fock Space

    If `x = |\lambda, n\rangle`, then `H\cdot x = \sum_{\mu} s_{\lambda/\mu}|\mu, n\rangle`.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: R = SymmetricFunctions(QQ); s = R.s()
        sage: F = FermionicFockSpace(s)
        sage: act_by_H(F(([2,1],0)))
        (s[2,1])*|[], 0> + (s[1,1]+s[2])*|[1], 0> + s[1]*|[1, 1], 0> + s[1]*|[2], 0> +
         s[]*|[2, 1], 0>
    """
    J = Current()
    H = Hamiltonian().truncate(deg(x) + 1).symmetric_function()
    p = J._p

    res = x.parent().zero()
    for (hm, hc) in H.monomial_coefficients().items():
        res += hc*p(hm)*J.act_by_partition(hm, 1, x)
    return res


def H_mat_coeff(lam, mu):
    r"""
    Compute the coefficient of `|\mu, n\rangle` in `H\cdot |\lambda, n\rangle`.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import *
        sage: lam = [3,2,1]; s = SymmetricFunctions(QQ).s()
        sage: all(H_mat_coeff(lam, mu) == s(lam).skew_by(s(mu)) for mu in Partitions(outer=lam))
        True
    """
    R = SymmetricFunctions(QQ)
    s = R.s()
    F = FermionicFockSpace(s)
    return act_by_H(F((lam, 0))).coefficient((mu, 0))


# degree of element of F
def deg(x):
    return max((m[0].size() for (m, _) in x.monomial_coefficients().items()), default=0)


class HalfVertexOperator():
    def __init__(self, f):
        r"""
        Initialize ``self``.
        """
        self._stream = Stream_cauchy_compose(
            Stream_function(lambda n: ZZ(1)/factorial(n), False, 0),
            Stream_function(f, False, 1),
            False
        )

    def __getitem__(self, i):
        if i in NonNegativeIntegers():
            return self._stream[i]
        raise ValueError("Invalid input")


class AbstractVertexOperator(SageObject):
    """
    Abstract class for Vertex Operators
    """
    def __init__(self, fockspace):
        r"""
        Initialize ``self``.
        """
        self.fockspace = fockspace  # TODO: Check input validity
        super().__init__()
        # super().__init__(VertexOperatorMonoid(fockspace))

    @abstract_method
    def act_on(self, m, f): pass

    @abstract_method
    def full_action(self, f): pass

    @abstract_method
    def matrix_coefficient(self, bra, ket): pass

    def __mul__(self, V):
        return ProductOfVertexOperators([self, V])


class VertexOperator(AbstractVertexOperator):
    r"""
    The action of a Vertex Operator on the Bosonic Fock space. Users should not
    create instances of this class directly, but instead use one of the defined
    subclasses.

    INPUT:

    - ``pos`` -- function taking in nonnegative integers indexing the positive
      half of the vertex operator
    - ``neg`` -- function taking in nonnegative integers indexing the negative
      half of the vertex operator
    - ``cutoff`` -- function taking in symmetric functions which determines how
      far to expand the vertex operator (default: ``lambda x: max(x.degree(), 1)``)
    - ``dcharge`` -- integer (default: ``1``) indicating how this vertex operator
      should change the charge of an element
    - ``fockspace``-- the space that the vertex operators are acting on

    """
    def __init__(self, pos, neg, cutoff=lambda x: max(x.degree(), 1), dcharge=1, fockspace=None):
        r"""
        Initialize ``self``.
        """
        self.pos = HalfVertexOperator(pos)
        self.neg = HalfVertexOperator(neg)
        if fockspace is None:
            fockspace = LaurentPolynomialRing(SymmetricFunctions(QQ).s(), names=('w'))

        # self.spectral = LazyLaurentSeriesRing(self.fockspace, names = ('z',))
        self.cutoff = cutoff
        self.dcharge = dcharge
        super().__init__(fockspace)
    # def act_on_fock_space_element(self, x):
    #     return self.spectral(lambda n: self.act_by_mode(n,x), valuation = -1)

    def act_on(self, i, x):
        r"""
        Action of the ``i``'th Fourier mode of ``self`` on element ``x`` of
        ``self.fockspace``.
        """
        res = self.fockspace.zero()
        # for each component of constant charge, compute the action
        for (charge, fn) in x.monomial_coefficients().items():
            res += self.fockspace.gen()**(charge+self.dcharge)*self._act_on_sym(i, fn, charge)

        return res

    def full_action(self, x, cutoff=None):
        r"""
        The full action of ``self`` on a Fock space element ``x``.

        INPUT:

        - ``x`` -- element of ``self.fockspace``
        """
        F = Family(ZZ, lambda i: self.act_on(i, x))
        return F if cutoff is None else self._approximate(F, cutoff)

    def matrix_coefficient(self, bra, ket, cutoff=None):
        r"""
        Compute the matrix coefficient of ``self`` with vector ``ket``
        and dual vector ``bra``.

        INPUT:

        - ``bra``, ``ket`` -- An ordered pair (`\lambda`, ``c``)
          where `lambda` is an integer partition and ``c`` is an integer.
        """
        w = self.fockspace.gen()
        R = self.fockspace.base_ring()
        f = (w**ket[1])*R(ket[0])
        f_acted_on = self.full_action(f)
        F = Family(ZZ, lambda i: f_acted_on[i].monomial_coefficients().get(
                bra[1], self.fockspace.zero()).monomial_coefficients().get(Partition(bra[0]), self.fockspace.zero()))
        return F if cutoff is None else self._approximate(F, cutoff)

    def _approximate(self, family, cutoff):
        res = {}
        for i in range(-cutoff, cutoff + 1):
            c = family[i]
            if c != 0:
                res[i] = c
        return res
    # act = act_by_mode

    def _act_on_sym(self, i, x, c):
        r"""
        Action of the ``i``'th Fourier mode of ``self`` on a homoegeneous element
        of ``self.fockspace``.
        """
        R = self.fockspace.base_ring()
        op = self._get_operator(i, self.cutoff(x), c)

        # op is a scalar
        if op in R.base_ring():
            return op*x

        res = R.zero()
        for (m, c) in x.monomial_coefficients().items():
            res += c*self._act_on_basis(op, m)
        return res

    @cached_method
    def _act_on_basis(self, op, x):
        """
        Action of a differential operator ``op`` on a basis element of the fock space.
        """
        R = self.fockspace.base_ring()
        p = R.symmetric_function_ring().p()

        res = R.zero()
        for (m, c) in op.monomial_coefficients().items():
            par1 = self._d_to_par(m[0].dict())
            par2 = self._d_to_par(m[1].dict())
            res += c*R(x).skew_by(p(par2)) * p(par1) / prod(par1)
        return res

    @abstract_method
    def _get_operator(self, i, x, c): pass

    def _d_to_par(self, m):
        """
        compute the partition corresponding to a dictionary ``m`` where the value
        of key ``i`` in ``m`` corresponds to the number of parts of size ``i``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.CreationOperator(B)
            sage: V._d_to_par({1:3, 4:1})
            [4, 1, 1, 1]

        """
        res = []
        for i in m:
            res += [i]*int(m[i])
        return Partition(sorted(res, reverse=True))


class ProductOfVertexOperators(AbstractVertexOperator):
    r"""
    Product of vertex operators.

    INPUT:

    - ``vertex_ops`` -- a ``list`` of vertex operators on the same Fock space.

    EXAMPLES::

        sage: B = vertex_operators.BosonicFockSpace()
        sage: Cre = vertex_operators.CreationOperator(B)
        sage: Ann = vertex_operators.AnnihilationOperator(B)
        sage: P1 = Cre*Cre
        sage: P1.act_on([4,3], B.one())
        s[3, 3]*w^2
        sage: P2 = Ann*Cre
        sage: P2.act_on([-3,3], B.one())
        s[]

    """

    def __init__(self, vertex_ops):
        r"""
        Initialize ``self``.
        """
        self.vertex_ops = flatten([x if isinstance(x, VertexOperator) else x.vertex_ops for x in vertex_ops])
        assert all(op.fockspace is vertex_ops[0].fockspace for op in self.vertex_ops)
        self._num_ops = len(vertex_ops)
        super().__init__(self.vertex_ops[0].fockspace)

    def act_on(self, mon, x):
        r"""
        Compute the action of a Fourier mode of ``self`` on ``x``.

        Let `X^i_j` denote the `j`'th Fourier mode of the `i`'th vertex operator of ``self``.
        This method computes `X^1_{mon_1}\cdots X^m_{mon_m}\left|x\right\rangle`.

        INPUT:

        - ``mon`` -- a list of integers of length equal to the number of vertex operators
        - ``x`` -- an element of ``self.fockspace`` to be acted on.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: Cre = vertex_operators.CreationOperator(B)
            sage: P = Cre*Cre*Cre
            sage: P.act_on([3, 2, 1],B.one())
            s[1, 1, 1]*w^3
            sage: P.act_on([3, 1, 2],B.one())
            -s[1, 1, 1]*w^3
            sage: Ann = vertex_operators.AnnihilationOperator(B)
            sage: P = Ann*Ann*Ann
            sage: P.act_on([4, 3, 2],B.one())
            -s[3]*w^-3
            sage: P.act_on([4, 2, 3],B.one())
            s[3]*w^-3
        """
        if len(mon) != len(self.vertex_ops):
            raise ValueError
        for i in range(len(mon) - 1, -1, -1):
            x = self.vertex_ops[i].act_on(mon[i], x)
            if x == 0:
                break
        return self.fockspace(x)

    def full_action(self, f, cutoff=None):
        r"""
        Compute the full action of ``self`` on ``f``.

        INPUT:

        - ``f`` -- Fock space element
        - ``cutoff`` -- (default: None) positive integer

        OUTPUT:

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.CreationOperator(B)
            sage: P = V*V
            sage: f = P.full_action(B.one())
            sage: f[2,1]
            s[1, 1]*w^2
        """
        F = Family(ZZ**self._num_ops, lambda m: self.act_on(m, f))
        return F if cutoff is None else self._approximate(F, cutoff)

    def matrix_coefficient(self, bra, ket, cutoff=None):
        r"""
        Compute the matrix coefficient of ``self`` corresponding to ``bra`` and ``ket``.

        INPUT:

        - ``bra``, ``ket`` -- ordered pair (`\lambda`, c) consisting of an integer
          partition and an integer, indexing a basis element of the Fock space.
        - ``cutoff`` -- (default: None) Nonnegative integer indicating how far to expand the vertex operator.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: Cre = vertex_operators.CreationOperator(B)
            sage: P = Cre*Cre
            sage: P.matrix_coefficient(([],3), ([],1), cutoff=3)  # example of eq 2.21 in [AZ13]
            {(1, 2): -1, (2, 1): 1}
        """
        w = self.fockspace.gen()
        R = self.fockspace.base_ring()
        f = (w**ket[1])*R(ket[0])
        f_acted_on = self.full_action(f)

        F = Family(ZZ**self._num_ops, lambda m: f_acted_on[m].monomial_coefficients().get(
                bra[1], self.fockspace.zero()).monomial_coefficients().get(Partition(bra[0]), self.fockspace.zero()))
        return F if cutoff is None else self._approximate(F, cutoff)

    def _approximate(self, F, cutoff=3):
        r"""
        Approximate the nonzero output of a family ``F`` with integer tuple index set

        Computes the nonzero values of fam for each integer tuple with maximum entry magnitude ``cutoff``.

        INPUT:

        - ``F`` -- Family
        - ``cutoff`` -- (default: 3) Nonnegative integer
        """
        from itertools import product
        res = {}
        for m in product(range(-cutoff, cutoff + 1), repeat=self._num_ops):
            c = F[m]
            if c != 0:
                res[m] = c
        return res

    def vacuum_expectation(self, cutoff=None):
        r"""
        Computes the matrix coefficient `\langle \varnothing | X | \varnothing \rangle`

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: Cre = vertex_operators.CreationOperator(B)
            sage: Ann = vertex_operators.AnnihilationOperator(B)
            sage: P = Ann*Cre
            sage: P.vacuum_expectation(cutoff=4)
            {(-4, 4): 1, (-3, 3): 1, (-2, 2): 1, (-1, 1): 1, (0, 0): 1}
            sage: P = Cre*Ann
            sage: P.vacuum_expectation(cutoff=4)
            {(-4, 4): 1, (-3, 3): 1, (-2, 2): 1, (-1, 1): 1}

        This verifies that, letting `\psi(z), \psi^*(w)` denote the fermionic fields,
        the matrix coefficient of their product is given by `\langle \varnothing | \psi^*(w)\psi^(z)|\varnothing\rangle = \frac{w}{w-z}`
        """
        return self.matrix_coefficient(([], 0), ([], 0), cutoff)

    # TODO: think of a better way to meaningfully represent the product without way too much text.
    def _repr_(self):
        return "Product of " + str(self.vertex_ops)[1:-1]


class CreationOperator(VertexOperator):
    r"""
    The image under the boson-fermion correspond of the fermionic field `\psi^*(z)`.

    Explicitly, the action on the bosonic fock space is given by the series

    .. MATH::

        \exp \left( \sum_{j \geq 1} x_j z^j \right) \exp\left(\sum_{j \geq 1} (-dx_j/j)z^{-j}\right)wz^{Q}

    EXAMPLES::

        sage: B = vertex_operators.BosonicFockSpace(); w = B.gen()
        sage: Cre = vertex_operators.CreationOperator(B)
        sage: Cre.act_on(-1, B.one())
        0
        sage: Cre.act_on(0, B.one())
        s[]*w
        sage: Cre.act_on(1, w)
        s[]*w^2
        sage: Cre.act_on(0, w^-1+ w^-2)
        s[2]*w^-1 + s[1]
        sage: t = Cre.act_on(1, Cre.act_on(-1, w^-2)); t
        s[2, 1]
        sage: Cre.act_on(3, t)
        s[3, 2, 1]*w
        sage: t + Cre.act_on(-1, Cre.act_on(1, w^-2))
        0
    """
    def __init__(self, fockspace):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.CreationOperator(B)
            sage: TestSuite(V).run(skip='_test_pickling')
        """
        self._weyl_algebra = DifferentialWeylAlgebra(QQ, names=('x'), n=PlusInfinity())
        self._x, self._dx = self._weyl_algebra.gens()
        self._spectral = LazyLaurentSeriesRing(fockspace, names=('z',))
        super().__init__(lambda n: self._x[n], lambda n: -self._dx[n]/n, dcharge=1, fockspace=fockspace)

    def _get_operator(self, i, cutoff, c):
        r"""
        Compute the coefficient of `z^{i - c}` in the vertex operator

        .. MATH::

            \exp \left( \sum_{j \geq 1} x_j z^j \right) \exp\left(\sum_{j \geq 1} (-dx_j/j)z^{-j}\right)wz^{Q}

        Thought the coefficient is an infinite sum of weyl algebra elements, we
        only compute the first ``cutoff`` many terms.

        INPUT:

        - ``i`` -- (integer)
        - ``cutoff`` -- (integer) how far to expand the product
        - ``c`` -- (integer)
        """
        op = 0
        for j in range(cutoff+1):
            if i + j - c < 0:
                continue
            op += self.pos[i + j - c]*self.neg[j]
        return op

    # def full_action(self, x):

    #     return self._spectral(lambda i: self.act_on(i, x), valuation=-max(y.degree() for (_, y) in x.monomial_coefficients().items()))

    # def matrix_coefficient(self, bra, ket):
    #     r"""
    #     Compute the matrix coefficient of ``self`` with vector ``ket``
    #     and dual vector ``bra``.

    #     INPUT:

    #     - ``bra``, ``ket`` -- An ordered pair (`\lambda`, ``c``)
    #       where `lambda` is an integer partition and ``c`` is an integer.

    #     EXAMPLES::

    #         sage: from sage.algebras.vertex_operators import *
    #         sage: B = BosonicFockSpace()
    #         sage: V = CreationOperator(B)
    #         sage: V.matrix_coefficient(([1],1),([2,1],0))
    #         s[]*z^-2 + O(s[]*z^3)
    #     """
    #     w = self.fockspace.gen()
    #     R = self.fockspace.base_ring()
    #     f = (w**(ket[1]))*R(ket[0])
    #     zero = self.fockspace.zero()
    #     return self._spectral(lambda i: self.act_on(i, f).monomial_coefficients().get(
    #             bra[1], zero).monomial_coefficients().get(Partition(bra[0]), zero), valuation=-Partition(ket[0]).size()-1)

    def act_by_clifford_gen(self, i, x):
        return self.act_on(i, x)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.CreationOperator(B); V
            Creation Vertex Operator acting on Univariate Laurent
            Polynomial Ring in w over Symmetric Functions over Rational Field
            in the Schur basis

        """
        return f"Creation Vertex Operator acting on {self.fockspace}"


class AnnihilationOperator(VertexOperator):
    r"""
    The image under the boson-fermion correspond of the fermionic field `\psi^*(z)`.

    Explicitly, the action on the bosonic fock space is given by the series

    .. MATH::

        \exp \left( \sum_{j \geq 1} -x_j z^j \right) \exp\left(\sum_{j \geq 1} (dx_j/j)z^{-j}\right)z^{-Q}w^{-1}

    .. WARNING::

        Following the literature, the operator corresponding to the coefficient of
        `z^i` is `\psi_{-i}`, **not** `\psi_i`.

    EXAMPLES::

        sage: B = vertex_operators.BosonicFockSpace(); w = B.gen()
        sage: Ann = vertex_operators.AnnihilationOperator(B)
        sage: Ann.act_on(0, B.one())
        0
        sage: Ann.act_on(1, B.one())
        s[]*w^-1
        sage: Ann.act_on(2, w^-1)
        s[]*w^-2
        sage: Ann.act_on(1, w)
        -s[1]
        sage: Ann.act_on(2, Ann.act_on(0, w^2))
        -s[2, 1]
        sage: Ann.act_on(2, Ann.act_on(0, w^2)) + Ann.act_on(0, Ann.act_on(2, w^2))
        0
    """
    def __init__(self, fockspace):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.AnnihilationOperator(B)
            sage: TestSuite(V).run(skip='_test_pickling')
        """

        self._weyl_algebra = DifferentialWeylAlgebra(QQ, names=('x'), n=PlusInfinity())
        self._x, self._dx = self._weyl_algebra.gens()
        super().__init__(lambda n: -self._x[n], lambda n: self._dx[n]/n, dcharge=-1, fockspace=fockspace)

    def _get_operator(self, i, cutoff, c):
        r"""
        Compute the coefficient of `z^{i + c - 1}` in the vertex operator

        .. MATH::

            \exp \left( \sum_{j \geq 1} -x_j z^j \right) \exp\left(\sum_{j \geq 1} (dx_j/j)z^{-j}\right)z^{-Q}w^{-1}

        Thought the coefficient is an infinite sum of weyl algebra elements, we
        only compute the first ``cutoff`` many terms.

        INPUT:

        - ``i`` -- (integer)
        - ``cutoff`` -- (integer) how far to expand the product
        - ``c`` -- (integer)
        """
        op = 0
        for j in range(cutoff+1):
            if j + i + c - 1 < 0:
                continue

            op += self.pos[j + i + c - 1]*self.neg[j]
        return op

    def act_by_clifford_gen(self, i, x):
        """
        Action of the coefficient of `z^{-i}` on Fock space element ``x``.
        """
        return self.act_on(-i, x)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.AnnihilationOperator(B); V
            Annihilation Vertex Operator acting on Univariate Laurent
            Polynomial Ring in w over Symmetric Functions over Rational Field
            in the Schur basis

        """
        return f"Annihilation Vertex Operator acting on {self.fockspace}"


# class VertexOperatorMonoid(Parent, UniqueRepresentation):
#     def __init__(self, fockspace):
#         from sage.categories.monoids import Monoids
#         self.fockspace = fockspace
#         super().__init__(category=Monoids())

#     def _element_constructor_(self, x):
#         if isinstance(x, AbstractVertexOperator):
#             return self.element_class(self, x)

#     Element = AbstractVertexOperator
