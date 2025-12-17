r"""
Vertex Operators

AUTHORS:

- Joseph McDonough (2025-08-04): Initial version
"""

from sage.algebras.weyl_algebra import DifferentialWeylAlgebra
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.data_structures.stream import Stream_function, Stream_cauchy_compose
from sage.functions.other import factorial
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
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
    """
    def __init__(self, R=QQ):
        r"""
        Initialize ``self``.

        INPUT:

        - ``R`` -- a ring (default: ``QQ``)

        TESTS::

            sage: from sage.algebras.vertex_operators import FermionicFockSpace
            sage: F = FermionicFockSpace(QQ)
            sage: TestSuite(F).run()
        """
        index_set = cartesian_product((Partitions(), ZZ))
        CombinatorialFreeModule.__init__(self, R, index_set, prefix='', bracket='')

    def _repr_term(self, m):
        r"""
        String representation of a basis element ``m`` of ``self``

        EXAMPLES::

            sage: from sage.algebras.vertex_operators import FermionicFockSpace
            sage: F = FermionicFockSpace()
            sage: print(F._repr_term(([3,2,1], -4)))
            |[3, 2, 1], -4>
        """
        return '|' + str(m[0]) + ', ' + str(m[1]) + '>'


def BosonicFockSpace(R=QQ, sym_basis='s', name='w'):
    r"""
    Bosonic Fock space.

    INPUT:

    - ``R`` -- (default: ``QQ``) a ring
    - ``sym_basis`` -- (default: ``'s'``) string or realization of ``SymmetricFunctions(R)``
    - ``name`` -- (default ``'w'``) string indicating the name of the Laurent variable.

    OUTPUT:

    The Laurent polynomial ring over the ring of symmertic functions over ``R`` in
    the specified basis.

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import BosonicFockSpace
        sage: B = BosonicFockSpace(); B
        Univariate Laurent Polynomial Ring in w over Symmetric Functions
        over Rational Field in the Schur basis
        sage: B.an_element()
        s[]*w
        sage: Sym = SymmetricFunctions(QQ); Sym.inject_shorthands(verbose=False)
        sage: BosonicFockSpace(QQ, s)
        Univariate Laurent Polynomial Ring in w over
        Symmetric Functions over Rational Field in the Schur basis
        sage: BosonicFockSpace(FractionField(QQ['q','t']), 'p')
        Univariate Laurent Polynomial Ring in w over Symmetric Functions
        over Fraction Field of Multivariate Polynomial Ring in q, t over
        Rational Field in the powersum basis
    """
    Sym = SymmetricFunctions(R)
    if isinstance(sym_basis, str):
        b = getattr(Sym, sym_basis, None)
        if b is None or b() not in Sym.realizations():
            raise ValueError("Unknown symmetric function basis")
        return LaurentPolynomialRing(b(), names=(name,))
    if sym_basis.base_ring() is not R:
        raise ValueError("")
    return LaurentPolynomialRing(sym_basis, names=(name,))


class HalfVertexOperator(SageObject):
    r"""
    A half vertex operator.

    INPUT:

    - ``f`` -- A function from the positive integers to an infinite weyl algebra

    OUTPUT:

    A stream where the `i`'th value is the coefficient of `z^i` in `\exp \left(\sum_{j \geq 1} f(j) z^j\right)`

    EXAMPLES::

        sage: from sage.algebras.vertex_operators import HalfVertexOperator
        sage: W.<x> = DifferentialWeylAlgebra(QQ, n=oo)
        sage: H = HalfVertexOperator(lambda n: x[n])
        sage: H[3]
        1/6*x[1]^3 + x[1]*x[2] + x[3]
    """
    def __init__(self, f):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.vertex_operators import HalfVertexOperator
            sage: W.<x> = DifferentialWeylAlgebra(QQ, n=oo)
            sage: H = HalfVertexOperator(lambda n: x[n])
            sage: TestSuite(H).run(skip='_test_pickling')
        """
        self._stream = Stream_cauchy_compose(
            Stream_function(lambda n: ZZ(1)/factorial(n), False, 0),
            Stream_function(f, False, 1),
            False
        )

    def __getitem__(self, i):
        r"""
        Get the i'th value of ``self``

        sage: from sage.algebras.vertex_operators import HalfVertexOperator
        sage: W.<x> = DifferentialWeylAlgebra(QQ, n=oo); dx = W.differentials()
        sage: H = HalfVertexOperator(lambda n: x[n] + dx[n])
        sage: H[2]
        1/2*dx[1]^2 + x[1]*dx[1] + dx[2] + 1/2*x[1]^2 + x[2] + 1/2
        sage: H[-1]
        Traceback (most recent call last):
        ...
        ValueError: Invalid input
        """
        if i in NonNegativeIntegers():
            return self._stream[i]
        raise ValueError("Invalid input")


class AbstractVertexOperator(SageObject):
    """
    Abstract class for Vertex Operators
    """
    def __init__(self, fockspace):
        r"""
        Initialize ``self``

        TESTS::

            sage: from sage.algebras.vertex_operators import AbstractVertexOperator
            sage: A = AbstractVertexOperator(vertex_operators.BosonicFockSpace())
            sage: TestSuite(A).run(skip='_test_pickling')
        """
        self.fockspace = fockspace  # TODO: Check input validity
        super().__init__()
        # super().__init__(VertexOperatorMonoid(fockspace))

    @abstract_method(optional=True)
    def act_on(self, m, f):
        r"""
        Act on an element of the Fock space by a mode of ``self``.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: V.act_on(3, vertex_operators.BosonicFockSpace().one())
            s[3]*w
        """

    @abstract_method(optional=True)
    def full_action(self, f):
        r"""
        Act on an element of the Fock space by the full vertex operator.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: V.full_action(vertex_operators.BosonicFockSpace().one())
            Lazy family (<lambda>(i))_{i in Integer Ring}
        """

    @abstract_method(optional=True)
    def matrix_coefficient(self, bra, ket):
        r"""
        Calculate the matrix coefficient of ``self`` for vector ``ket`` and
        dual vector ``bra``.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: mc = V.matrix_coefficient(([],1), ([],0))
            sage: mc[0]
            1
        """

    def __mul__(self, V):
        r"""
        Construct the product of ``self`` with ``V``

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: (V*V).act_on([2,1], vertex_operators.BosonicFockSpace().one())
            s[1, 1]*w^2
        """
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


    EXAMPLES::

        sage: V = vertex_operators.CreationOperator()
        sage: B = V.fockspace; s = B.base_ring()
        sage: V.full_action(B(s[1])*B.gen()^3, cutoff=5)
        {2: -s[]*w^4, 4: s[1, 1]*w^4, 5: s[2, 1]*w^4}
    """
    def __init__(self, pos, neg, cutoff=lambda x: max(x.degree(), 1), dcharge=1, fockspace=None):
        r"""
        Initialize ``self``.

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

        EXAMPLES::

            sage: from sage.algebras.vertex_operators import VertexOperator
            sage: W.<x> = DifferentialWeylAlgebra(QQ, n=oo)
            sage: V = VertexOperator(lambda n: x[n], lambda n: x[n])
            sage: TestSuite(V).run(skip="_test_pickling")
        """
        self.pos = HalfVertexOperator(pos)
        self.neg = HalfVertexOperator(neg)
        if fockspace is None:
            fockspace = BosonicFockSpace()

        self.cutoff = cutoff
        self.dcharge = dcharge
        super().__init__(fockspace)

    def act_on(self, i, f):
        r"""
        Action of the ``i``'th Fourier mode of ``self`` on element ``f`` of
        ``self.fockspace``.

        INPUT:

        - ``i`` -- integer
        - ``f`` -- element of ``self.fockspace``.

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: B = V.fockspace; s = B.base_ring()
            sage: V.act_on(4, B(s[2,1]))
            -s[3, 2, 1]*w^-1
        """
        res = self.fockspace.zero()
        # for each component of constant charge, compute the action
        for (charge, fn) in f.monomial_coefficients().items():
            res += self.fockspace.gen()**(charge+self.dcharge)*self._act_on_sym(i, fn, charge)

        return res

    def full_action(self, f, cutoff=None):
        r"""
        The full action of ``self`` on a Fock space element ``x``.

        If the ``cutoff`` parameter is ``None``, we return a family that
        indexes the action of each mode. Otherwise, returns a dictionary with
        key, value pairs corresponding to the action of modes of magnitude at
        most ``cutoff`` with nonzero value.

        INPUT:

        - ``f`` -- element of ``self.fockspace``
        - ``cutoff`` -- ``None`` or positive integer (default: ``None``)

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: B = V.fockspace
            sage: V.full_action(B.one())
            Lazy family (<lambda>(i))_{i in Integer Ring}
            sage: V.full_action(B.one(), cutoff=5)
            {1: s[]*w^-1, 2: -s[1]*w^-1, 3: s[1, 1]*w^-1,
            4: -s[1, 1, 1]*w^-1, 5: s[1, 1, 1, 1]*w^-1}
        """
        F = Family(ZZ, lambda i: self.act_on(i, f))
        return F if cutoff is None else self._approximate(F, cutoff)

    def matrix_coefficient(self, bra, ket, cutoff=None):
        r"""
        Compute the matrix coefficient of ``self`` with vector ``ket``
        and dual vector ``bra``.

        If the ``cutoff`` parameter is ``None``, we return a family that
        indexes the value of the matrix coefficient for `z^i`. Otherwise,
        returns a dictionary with key, value pairs corresponding to the nonzero
        values with key magnitude at most ``cutoff``.

        INPUT:

        - ``bra``, ``ket`` -- An ordered pair (`\lambda`, ``c``)
          where `lambda` is an integer partition and ``c`` is an integer.
        - ``cutoff`` -- ``None`` or positive integer (default: ``None``)

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: B = V.fockspace
            sage: V.matrix_coefficient(([],-1), ([],0))
            Lazy family (<lambda>(i))_{i in Integer Ring}
            sage: V.matrix_coefficient(([],-1), ([],0), cutoff=5)
            {1: 1}
        """
        w = self.fockspace.gen()
        R = self.fockspace.base_ring()
        f = (w**ket[1])*R(ket[0])
        f_acted_on = self.full_action(f)
        F = Family(ZZ, lambda i: f_acted_on[i].monomial_coefficients().get(
                bra[1], self.fockspace.zero()).monomial_coefficients().get(Partition(bra[0]), self.fockspace.zero()))
        return F if cutoff is None else self._approximate(F, cutoff)

    def _approximate(self, family, cutoff):
        r"""
        Approximates the nonzero values of a ``ZZ`` indexed family

        INPUT:

        - ``family`` -- a family
        - ``cutoff`` -- a nonnegative integer

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: B = V.fockspace
            sage: F = V.full_action(B.gen()^2)
            sage: V._approximate(F, 5)
            {2: s[]*w^3, 3: s[1]*w^3, 4: s[2]*w^3, 5: s[3]*w^3}
        """
        res = {}
        for i in range(-cutoff, cutoff + 1):
            c = family[i]
            if c != 0:
                res[i] = c
        return res

    def _act_on_sym(self, i, f, c):
        r"""
        Action of the ``i``'th Fourier mode of ``self`` on a homoegeneous element
        of ``self.fockspace``.

        INPUT:

        - ``i`` -- an integer indexing the Fourier mode
        - ``f`` -- a ``SymmetricFunction``
        - ``c`` -- an integer indexing the charge of ``f``.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: s = V.fockspace.base_ring()
            sage: V._act_on_sym(1, s[2] + s[1, 1], 0)
            s[1, 1, 1]
            sage: V._act_on_sym(1, s[2] + s[1, 1], -1)
            s[2, 1, 1] + s[2, 2]
            sage: V._act_on_sym(1, s[2] + s[1, 1], 2)
            -s[1]
        """
        R = self.fockspace.base_ring()
        op = self._get_operator(i, self.cutoff(f), c)

        # op is a scalar
        if op in R.base_ring():
            return op*f

        res = R.zero()
        for (m, c) in f.monomial_coefficients().items():
            res += c*self._act_on_basis(op, m)
        return res

    @cached_method
    def _act_on_basis(self, op, x):
        """
        Action of a differential operator ``op`` on a basis element of the fock space.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: s = V.fockspace.base_ring()
            sage: op = V._get_operator(1,5,0)
            sage: V._act_on_basis(op, s[4])
            -s[3, 2]
        """
        R = self.fockspace.base_ring()
        p = R.symmetric_function_ring().p()

        res = R.zero()
        for (m, c) in op.monomial_coefficients().items():
            par1 = self._d_to_par(m[0].dict())
            par2 = self._d_to_par(m[1].dict())
            res += c*R(x).skew_by(p(par2)) * p(par1) / prod(par1)
        return res

    @abstract_method(optional=True)
    def _get_operator(self, i, cutoff, c):
        """
        Compute the Weyl algebra element necessary to find the action of a mode of
        ``self``.

        Since the method for computing this operator varies wildly, we leave it to
        the specific subclass to determine how to implement it.

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: s = V.fockspace.base_ring()
            sage: V._get_operator(1, 1, 1)
            1/2*x[1]^2*dx[1] - x[2]*dx[1] - x[1]
            sage: V._get_operator(1, 2, 1)
            -1/12*x[1]^3*dx[1]^2 + 1/2*x[1]*x[2]*dx[1]^2 - 1/2*x[3]*dx[1]^2 +
            1/2*x[1]^2*dx[1] - x[2]*dx[1] - 1/12*x[1]^3*dx[2] +
            1/2*x[1]*x[2]*dx[2] - 1/2*x[3]*dx[2] - x[1]
        """

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

        TESTS::

            sage: V = vertex_operators.CreationOperator()
            sage: P = V*V
            sage: TestSuite(P).run(skip='_test_pickling')
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

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: P = V*V
            sage: F = P.full_action(P.fockspace.one())
            sage: P._approximate(F, cutoff=2)
            {(0, 1): -s[]*w^2, (0, 2): -s[1]*w^2, (1, 0): s[]*w^2,
            (1, 2): -s[1, 1]*w^2, (2, 0): s[1]*w^2, (2, 1): s[1, 1]*w^2}

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
        r"""
        Return a string representation of ``self``

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: V*V
            Product of Creation Vertex Operator acting on Univariate Laurent
            Polynomial Ring in w over Symmetric Functions over Rational Field
            in the Schur basis, Creation Vertex Operator acting on Univariate
            Laurent Polynomial Ring in w over Symmetric Functions over Rational
            Field in the Schur basis
        """
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
    def __init__(self, fockspace=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = vertex_operators.BosonicFockSpace()
            sage: V = vertex_operators.CreationOperator(B)
            sage: TestSuite(V).run(skip='_test_pickling')
        """
        if fockspace is None:
            fockspace = BosonicFockSpace()
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

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: V._get_operator(2,1,1)
            -1/2*x[1]^2*dx[1] - x[2]*dx[1] + x[1]
        """
        op = 0
        for j in range(cutoff+1):
            if i + j - c < 0:
                continue
            op += self.pos[i + j - c]*self.neg[j]
        return op

    def act_by_clifford_gen(self, i, x):
        r"""
        Act by the Clifford algebra generator `\psi_i` on ``x``.

        EXAMPLES::

            sage: V = vertex_operators.CreationOperator()
            sage: B = V.fockspace
            sage: V.act_by_clifford_gen(6, B.one())
            s[6]*w
        """
        return self.act_on(i, x)

    def _repr_(self):
        """
        Return a string representation of ``self``

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
    def __init__(self, fockspace=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: TestSuite(V).run(skip='_test_pickling')
        """
        if fockspace is None:
            fockspace = BosonicFockSpace()
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

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: V._get_operator(2,1,1)
            -1/6*x[1]^3*dx[1] + x[1]*x[2]*dx[1] - x[3]*dx[1] + 1/2*x[1]^2 - x[2]
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

        EXAMPLES::

            sage: V = vertex_operators.AnnihilationOperator()
            sage: B = V.fockspace
            sage: V.act_by_clifford_gen(-6, B.one())
            -s[1, 1, 1, 1, 1]*w^-1
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
