r"""
Cohomology of coherent sheaves

We provide the cohomology of coherent sheaves over projective schemes by
implementing Maruyama's method. Use the
:meth:`~sage.schemes.projective.coherent_sheaf.CoherentSheaf.cohomology`
method of coherent sheaves on projective schemes to access this.

Maruyama's method is explained and proved in detail in [Kudo2017]_. Here we
summary the main results necessary to understand the implementation.

Let `M` be a graded module finitely generated over the homogeneous coordinate
ring `S` of the projective `r`-space over a field `k`. We aim for computing the
cohomology groups `H^q(\tilde M)` for the coherent sheaf `\tilde M`.

Let `S = k[x_0, x_1, \ldots, x_r]`. Then `M` is a quotient of the
free module `\bigoplus_{i=1}^t S` by a submodule. Let

.. MATH::

    0 \to \bigoplus_{j=1}^{t_{r+1}}S(-m^{(r+1)}_j) \overset{f_{r+1}}{\longrightarrow}
    \cdots \overset{f_1}{\longrightarrow}
    \bigoplus_{j=1}^{t_0}S(-m^{(0)}_j) \overset{f_0}{\longrightarrow} M \to 0

be a minimal free resolution of `M`. Then it induces a complex of (top)
cohomology groups

.. MATH::

    \bigoplus_{j=1}^{t_{i+1}} H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(i+1)}_j))
    \overset{H^r(f_{i+1})}{\longrightarrow}
    \bigoplus_{j=1}^{t_i}H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(i)}_j))
    \overset{H^r(f_{i})}{\longrightarrow}
    \bigoplus_{j=1}^{t_{i-1}} H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(i-1)}_j)),

where `i` runs from `1` to `r`. Now it holds that

.. MATH::

    H^q(\tilde M)\cong \ker H^r(f_{r-q})/\operatorname{im} H^r(f_{r-q+1})

for `1 \le q \le r - 1` and

.. MATH::

    H^r(\tilde M) \cong
    \bigoplus_{j=1}^{t_0}H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(0)}_j))
      / \mathrm{im} H^r(f_1)

and `\dim H^0(\tilde M)` can be computed by the formula

.. MATH::

    \begin{aligned}
    & \dim \bigoplus_{j=1}^{t_{0}} H^0(\mathcal{O}_{\Bold{P}^r}(-m^{(0)}_j))
    - \dim \bigoplus_{j=1}^{t_{r+1}} H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(r+1)}_j))
    + \dim \bigoplus_{j=1}^{t_{r}} H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(r)}_j)) \\
    & \quad - \rank H^0(f_1) - \rank H^r(f_r)
    \end{aligned}

in which the complex of (bottom) cohomology groups

.. MATH::

    \bigoplus_{j=1}^{t_{i+1}} H^0(\mathcal{O}_{\Bold{P}^r}(-m^{(i+1)}_j))
    \overset{H^0(f_{i+1})}{\longrightarrow}
    \bigoplus_{j=1}^{t_i} H^0(\mathcal{O}_{\Bold{P}^r}(-m^{(i)}_j))
    \overset{H^0(f_{i})}{\longrightarrow}
    \bigoplus_{j=1}^{t_{i-1}}H^0(\mathcal{O}_{\Bold{P}^r}(-m^{(i-1)}_j)),

where `i` runs from `1` to `r` is used.

The implemented algorithm works more generally for twisted coherent sheaves
and accepts as input shifted graded module `M(-n)` with shift `n`.

EXAMPLES:

We define the Fermat cubic surface (a curve in `\Bold{P}^2`) and compute
its cohomology groups::

    sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
    sage: X = P2.subscheme([x^4 + y^4 + z^4])
    sage: sh = X.structure_sheaf()
    sage: sh.betti(1)
    3

Internally the cohomology is computed by Maruyama's method::

    sage: c = sh.cohomology(); c
    Sheaf cohomology constructed from S(0) <-- S(-4) <-- 0
    sage: c.top_complex()
    Top Cohomology Complex induced from S(0) <-- S(-4) <-- 0
    sage: c.bottom_complex()
    Bottom Cohomology Complex induced from S(0) <-- S(-4) <-- 0
    sage: c.cohomology_group(1)
    Vector space quotient V/W of dimension 3 over Rational Field where
    V: Vector space of degree 3 and dimension 3 over Rational Field
    Basis matrix:
    [1 0 0]
    [0 1 0]
    [0 0 1]
    W: Vector space of degree 3 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: c.cohomology_group(1).dimension() == sh.betti(1)
    True
    sage: sh.betti(2)
    0
    sage: sh.cohomology(2)
    Vector space quotient V/W of dimension 0 over Rational Field where
    V: Vector space of dimension 0 over Rational Field
    W: Vector space of degree 0 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: c.cohomology_group(2) == sh.cohomology(2)
    True
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.flatten import flatten
from sage.misc.cachefunc import cached_method
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer

class CohomologyGroup(SageObject):
    r"""
    Cohomology group (either top or bottom) of the twisted structure sheaf
    of a projective space.

    INPUT:

    - ``S`` -- direct sum of copies of the coordinate ring of a projective space
    - ``shifts`` -- shifts of the component rings of ``S``
    - ``is_top`` -- boolean

    This represents `\bigoplus_{j=1}^{t_i}
    H^r(\mathcal{O}_{\Bold{P}^r}(-m^{(i)}_j))` for shifts `m^{(i)}_j`.

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: X = P2.subscheme([x^4 + y^4 + z^4])
        sage: c = X.structure_sheaf().cohomology()
        sage: c.bottom_complex().cohomology_group(0)
        Bottom Cohomology Group of dimension 1
        sage: c.top_complex().cohomology_group(0)
        Top Cohomology Group of dimension 0
    """
    def __init__(self, S, shifts, is_top):
        r"""
        Initialize ``self``.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf(twist=1).cohomology()
            sage: B0 = c.bottom_complex().cohomology_group(0)
            sage: TestSuite(B0).run(skip="_test_pickling")
            sage: T1 = c.top_complex().cohomology_group(1)
            sage: TestSuite(T1).run(skip="_test_pickling")
        """
        self._shifts = shifts
        self._is_top = is_top  # if false, then bottom

        n = S.ngens()

        basis = []
        summands_basis = []
        summands_index = []
        rank = 0
        for m in self._shifts:
            # list of integer vectors whose entries are all non-negative integers and sum to -m
            if self._is_top:
                L = [-vector(e) for e in IntegerListsLex(length=n, min_sum=m, max_sum=m, min_part=1)]
            else:
                L = [vector(e) for e in IntegerListsLex(length=n, min_sum=-m, max_sum=-m)]
            basis += L
            summands_basis.append(L)
            summands_index.append(rank)
            rank += len(L)

        self._summands_basis = summands_basis
        self._summands_index = summands_index
        self._basis = tuple(basis)
        self._vector_space = VectorSpace(S.base_ring(), rank)
        self._rank = Integer(rank)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf()._cohomology
            sage: c.bottom_complex().cohomology_group(0)
            Bottom Cohomology Group of dimension 1
            sage: c.top_complex().cohomology_group(1)
            Top Cohomology Group of dimension 3
        """
        start = "Top" if self._is_top else "Bottom"
        return f'{start} Cohomology Group of dimension {self._rank}'


class MaruyamaCohomologyComplex(UniqueRepresentation, SageObject):
    r"""
    The complex of (top/bottom) cohomology groups in Maruyama's method.

    INPUT:

    - ``resolution`` -- a minimal free resolution of a module
    - ``is_top`` -- boolean

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: X = P2.subscheme([x^4 + y^4 + z^4])
        sage: c = X.structure_sheaf(1).cohomology()
        sage: TC = c.top_complex()
        sage: all((TC.differential(i-1) * TC.differential(i)).is_zero() for i in range(4))
        True
        sage: BC = c.bottom_complex()
        sage: all((BC.differential(i-1) * BC.differential(i)).is_zero() for i in range(4))
        True
    """
    def __init__(self, resolution, is_top):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: BC = c.bottom_complex()
            sage: TestSuite(BC).run()
            sage: TC = c.top_complex()
            sage: TestSuite(TC).run()
        """
        self._resolution = resolution
        self._is_top = is_top
        self._base_ring = self._resolution.target().base_ring()
        self._coefficient_field = self._base_ring.base_ring()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.bottom_complex()
            Bottom Cohomology Complex induced from S(0) <-- S(-4) <-- 0
            sage: c.top_complex()
            Top Cohomology Complex induced from S(0) <-- S(-4) <-- 0
        """
        base = "Top" if self._is_top else "Bottom"
        return f"{base} Cohomology Complex induced from {self._resolution}"

    @cached_method
    def cohomology_group(self, i):
        r"""
        Return `i`-th cohomology group `\bigoplus_{j=1}^{t_i}
        H^a(\mathcal{O}_{\Bold{P}^r}(-m^{(i)}_j))` of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.bottom_complex().cohomology_group(1)
            Bottom Cohomology Group of dimension 0
            sage: c.top_complex().cohomology_group(1)
            Top Cohomology Group of dimension 3
        """
        S = self._base_ring
        if i < 0:  # resolution.shifts() raises an error for negative index
            shifts = []
        else:
            shifts = self._resolution.shifts(i)
        return CohomologyGroup(S, shifts, is_top=self._is_top)

    @cached_method
    def differential(self, i):
        r"""
        Return the `i`-th differential map `H^a(f_i)` of ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: BC = c.bottom_complex()
            sage: BC.differential(1)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 0 over Rational Field
            Codomain: Vector space of dimension 1 over Rational Field
            sage: TC = c.top_complex()
            sage: TC.differential(1)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 0 over Rational Field

            sage: BC.differential(-1)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 0 over Rational Field
            Codomain: Vector space of dimension 0 over Rational Field
            sage: BC.differential(0)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 1 over Rational Field
            Codomain: Vector space of dimension 0 over Rational Field
            sage: BC.differential(2)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 0 over Rational Field
            Codomain: Vector space of dimension 0 over Rational Field
            sage: BC.differential(3)
            Vector space morphism represented as left-multiplication by the matrix:
            []
            Domain: Vector space of dimension 0 over Rational Field
            Codomain: Vector space of dimension 0 over Rational Field
        """
        H1 = self.cohomology_group(i)
        H0 = self.cohomology_group(i - 1)
        if i < 0 or i > len(self._resolution):
            return H1._vector_space.hom([], codomain=H0._vector_space, side='right')
        if i == 0:
            rank = H1._vector_space.rank()
            return H1._vector_space.hom([[]]*rank, codomain=H0._vector_space, side='right')
        M = self._resolution.differential(i).matrix()
        K = self._coefficient_field
        zero = K.zero()

        assert M.ncols() == len(H1._summands_basis)
        assert M.nrows() == len(H0._summands_basis)

        A = []
        r = H0._rank
        for i, basis in enumerate(H1._summands_basis):
            for v in basis:
                image = [zero] * r
                for j in range(M.nrows()):
                    f = M[j, i]
                    basis = H0._summands_basis[j]
                    for c, m in zip(f.coefficients(), f.exponents()):
                        u = v + vector(m)
                        assert sum(u) == -H0._shifts[j]
                        if ((self._is_top and any(e < 0 for e in u))
                            or (not self._is_top and any(e >= 0 for e in u))):
                            continue
                        k = H0._summands_index[j] + basis.index(u)
                        image[k] += c
                A.append(vector(K, image))

        return H1._vector_space.hom(A, codomain=H0._vector_space, side='right')


class CoherentSheafCohomology(UniqueRepresentation, SageObject):
    r"""
    The cohomology groups of a coherent sheaf on a projective space.

    ALGORITHM:

    This class implements Maruyama's method to compute the cohomology group
    `H^q(\tilde M(n))` as a vector space over the base field `k` and
    `h^q(\tilde M(n)) = \dim_k H^q(\tilde M(n))`, where `n` denotes the twist.
    See the module :mod:`~sage.schemes.projective.cohomology` documentation
    for more information.

    The cohomology groups are returned as vector space quotients (over `k`).

    INPUT:

    - ``M`` -- a quotient of a free module over `S` by a submodule, where `S`
      is a multivariate polynomial ring
    - ``twist`` -- (default: 0) an integer

    EXAMPLES::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: X = P2.subscheme([x^4 + y^4 + z^4])
        sage: sh = X.structure_sheaf(1)  # twisted sheaf
        sage: sh.cohomology()
        Sheaf cohomology constructed from S(1) <-- S(-3) <-- 0
    """
    def __init__(self, M, twist=0):
        """
        Initialize.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: sh = X.structure_sheaf(1)  # twisted sheaf
            sage: c = sh.cohomology()
            sage: TestSuite(c).run()
        """
        shifts = [-twist] * M.cover().degree()
        self._resolution = M.relations().graded_free_resolution(shifts=shifts)
        self._top_complex = MaruyamaCohomologyComplex(self._resolution, True)
        self._bottom_complex = MaruyamaCohomologyComplex(self._resolution, False)
        base_ring = self._resolution.target().base_ring()
        self._projective_space_dimension = base_ring.ngens() - 1
        self._coefficient_field = base_ring.base_ring()

    def _repr_(self):
        r"""
        Return the string representation.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: X.structure_sheaf().cohomology()
            Sheaf cohomology constructed from S(0) <-- S(-4) <-- 0
        """
        return f"Sheaf cohomology constructed from {self._resolution}"

    def resolution(self):
        r"""
        Return the minimal free resolution used to construct ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.resolution()
            S(0) <-- S(-4) <-- 0
        """
        return self._resolution

    def top_complex(self):
        r"""
        Return the complex of top cohomology groups of the resolution
        used to construct ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.top_complex()
            Top Cohomology Complex induced from S(0) <-- S(-4) <-- 0
        """
        return self._top_complex

    def bottom_complex(self):
        r"""
        Return the complex of bottom cohomology groups of the resolution
        used to construct ``self``.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.bottom_complex()
            Bottom Cohomology Complex induced from S(0) <-- S(-4) <-- 0
        """
        return self._bottom_complex

    def cohomology_group(self, q):
        r"""
        Return the `q`-th cohomology group `H^q(\tilde M)` of ``self``
        as a vector space.

        INPUT:

        - ``q`` -- non-negative integer

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.cohomology_group(1)
            Vector space quotient V/W of dimension 3 over Rational Field where
            V: Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            W: Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: c.cohomology_group(0)
            Vector space of dimension 1 over Rational Field
            sage: c.cohomology_group(-2)
            Vector space of dimension 0 over Rational Field
            sage: c.cohomology_group(6)
            Vector space of dimension 0 over Rational Field
        """
        r = self._projective_space_dimension

        if q == r:
            return self._top_complex.cohomology_group(0)._vector_space.quotient(self._top_complex.differential(1).image())

        if 1 <= q and q < r:
            return self._top_complex.differential(r - q).kernel().quotient(self._top_complex.differential(r - q + 1).image())

        if q == 0:
            # 0th cohomology group of (global sections) is not realized in
            # terms of cohomology groups and differential maps of twisted
            # structure sheaves of the projective space.
            return VectorSpace(self._coefficient_field, self.betti(0))

        return VectorSpace(self._coefficient_field, 0)

    def betti(self, q=None):
        r"""
        Return the `q`-th Betti number, which is the dimension `h^q(\tilde M)`
        of the `q`-th cohomology group `H^q(\tilde M)`.

        INPUT:

        - ``q`` -- (optional) non-negative integer

        OUTPUT:

        If ``q`` is not given, then this returns all Betti numbers
        `h^0, \ldots, h^d`, where `d` is the dimension of the projective
        space. Otherwise returns the `q`-th Betti number `h^q`.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^4 + y^4 + z^4])
            sage: c = X.structure_sheaf().cohomology()
            sage: c.betti(1)
            3
            sage: c.betti()
            [1, 3, 0]
        """
        r = self._projective_space_dimension
        if q is None:
            return [self.betti(k) for k in range(r+1)]

        if q == r:
            return self._top_complex.cohomology_group(0)._rank - self._top_complex.differential(1).rank()

        if 1 <= q and q < r:
            return (self._top_complex.differential(r - q).kernel().dimension()
                    - self._top_complex.differential(r - q + 1).rank())

        if q == 0:
            return (self._bottom_complex.cohomology_group(0)._rank
                    - self._top_complex.cohomology_group(r + 1)._rank
                    + self._top_complex.cohomology_group(r)._rank
                    - self._bottom_complex.differential(1).rank()
                    - self._top_complex.differential(r).rank())

        return Integer(0)
