r"""
Banff cluster algebras

This module provides tools for working with Banff cluster algebras, a
subclass of cluster algebras admitting a finite cover by acyclic cluster charts
via the Banff algorithm.

Every acyclic chart determines finitely many Laurent charts, so a Banff cluster algebra
naturally carries the structure of a finite Laurent intersection ring (FLIR).
This implementation therefore realizes a Banff cluster algebra both as

- a :class:`~sage.algebras.cluster_algebra.ClusterAlgebra`, and
- a :class:`~sage.algebras.flir.FLIR`.

Compared with general cluster algebras, this additional structure provides
effective algorithms for several problems that are difficult or unavailable in
general. In particular, it allows

- explicit membership testing by checking Laurentness in every chart,
- computation of factorizations of elements,
- computation of divisor groups and class groups, and
- other divisor-theoretic functionality inherited from FLIRs.

Constructing the FLIR structure requires running the Banff algorithm and computing
the Laurent cover, which may take a noticeable amount of time and memory.
Consequently, the FLIR structure is initialized lazily and is only constructed when one of the
corresponding methods is first called. It can also be initialized explicitly
during construction by passing ``flir=True``.



.. NOTE::

    For scalar multiplication, coerce scalar constants into
    the Banff cluster algebra explicitly. For example, write
    ``A(1)/(x0*x1)`` rather than ``1/(x0*x1)``.

    This ensures that division is performed inside the Banff cluster
    algebra and avoids Sage coercion issues with Python integers.

.. NOTE::

    For a general cluster algebra, it is not immediate whether an element
    of the ambient Laurent polynomial ring or ambient fraction field
    belongs to the cluster algebra. Membership is in general a subtle
    problem.

    For a Banff cluster algebra, the Banff/FLIR structure provides an
    effective membership test: an ambient expression ``f`` lies in ``A``
    if and only if ``f`` is Laurent in every chart of the associated
    Banff/FLIR cover.

    In this implementation, explicit membership testing is performed by
    ``A._check_membership(f)``.

.. WARNING::

    Coercion via ``A(f)`` should not be used as a membership test. In
    particular, if ``f`` is an element of the ambient Laurent polynomial
    ring or fraction field, the construction ``A(f)`` may succeed without
    performing a membership check.

    To certify that an element lies in the Banff cluster algebra, use
    ``A._check_membership(f)`` explicitly.

EXAMPLES::

    sage: B = Matrix([[0, 1], [-1, 0]])
    sage: A = BanffClusterAlgebra(B)
    sage: A.is_banff()
    True

We can compute generators and a presentation::

    sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
    sage: A = BanffClusterAlgebra(B)
    sage: A.generators()
    [x0^-1*x1 + x0^-1, x0*x1^-1 + x1^-1*x2, x1*x2^-1 + x2^-1]
    sage: Q, R, I = A.presentation()
    sage: R.base_ring() is QQ
    True
    sage: Q
    Quotient of Multivariate Polynomial Ring in X0, X1, X2, T0, T1, T2 over Rational Field by the ideal (X1*T1 - X0 - X2, X0*T0 - X1 - 1, X2*T2 - X1 - 1)


Divisor-theoretic functionality is inherited from the FLIR structure::

    sage: B = Matrix([[0, 1], [-1, 0]])
    sage: A = BanffClusterAlgebra(B)
    sage: f = A(A.gens()[0])
    sage: D = f.divisor()
    sage: D.parent() is A.Div()
    True
    sage: D
    1*PrimeDivisor(chart=('x0_0p', 'x0_1'), p=x0_1 + 1)

Membership in a Banff cluster algebra can be checked explicitly::

    sage: B = Matrix([[0, 1], [-1, 0]])
    sage: A = BanffClusterAlgebra(B)
    sage: x0, x1 = A.gens()
    sage: F = A.ambient().fraction_field()
    sage: f = F(1)/(F(x0)*F(x1))
    sage: A._check_membership(f)
    Traceback (most recent call last):
    ...
    ValueError: Not Laurent in chart ('x0_0p', 'x0_1').
    Substituted: x0_0p/(x0_1^2 + x0_1)

A Banff cluster algebra inherits the divisor-theoretic and factorization
functionality of finite Laurent intersection rings. The FLIR structure is
initialized lazily, when one of these methods is called for the first time.

For example, one can compute the divisor of an element::

    sage: A = ClusterAlgebra(['D', 5], scalars=QQ)
    sage: B = A.b_matrix()
    sage: A = BanffClusterAlgebra(B)
    sage: x0, x1, x2, x3, x4 = A.gens()
    sage: D = x0.divisor()
    sage: D.parent() is A.Div()
    True
    sage: D
    1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2', 'x0_3', 'x0_4'), p=x0_1 + 1)
    sage: E = ((x3*x4)/(x0*x1*x2)).divisor()
    Traceback (most recent call last):
    ...
    ValueError: Not Laurent in chart ('x0_0p', 'x0_1', 'x0_2', 'x0_3', 'x0_4').
    Substituted: x0_0p*x0_3*x0_4/(x0_1^2*x0_2 + x0_1*x0_2)

The divisor group can be accessed directly::

    sage: DivA = A.Div()
    sage: D.parent() is DivA
    True
    sage: DivA.zero()
    0

Divisors are additive. In particular, principal divisors satisfy the expected
multiplicative-to-additive relation::

    sage: (x0*x1).divisor() == x0.divisor() + x1.divisor()
    True
    sage: D1 = x1.divisor()
    sage: D1 + D
    1*PrimeDivisor(chart=('x0_0', 'x0_1p', 'x0_2', 'x0_3', 'x0_4'), p=x0_0*x0_2 + 1) +
     1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2', 'x0_3', 'x0_4'), p=x0_1 + 1)
    sage: D1 - D
    1*PrimeDivisor(chart=('x0_0', 'x0_1p', 'x0_2', 'x0_3', 'x0_4'), p=x0_0*x0_2 + 1) +
     -1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2', 'x0_3', 'x0_4'), p=x0_1 + 1)
    sage: 2*D1
    2*PrimeDivisor(chart=('x0_0', 'x0_1p', 'x0_2', 'x0_3', 'x0_4'), p=x0_0*x0_2 + 1)

The divisor class group is computed from the FLIR cover::

    sage: ClA = A.class_group()
    sage: ClA
    Multiplicative Abelian group isomorphic to Z
    sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
    sage: A = BanffClusterAlgebra(B)
    sage: A.class_group()
    Multiplicative Abelian group isomorphic to Z

One can test whether a divisor is principal and, when it is principal, recover
a generator::

    sage: P1, P2, P3, P4 = A.extra_primes()
    sage: P1, P2, P3, P4   # these are the ht-1 primes containing initial variables
    (PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2'), p=x0_1 + 1),
     PrimeDivisor(chart=('x0_0', 'x0_1', 'x0_2p'), p=x0_1 + 1),
     PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2p'), p=x0_1 + 1),
     PrimeDivisor(chart=('x0_0', 'x0_1p', 'x0_2'), p=x0_0 + x0_2))
    sage: G = A.Div()
    sage: D1 = G({P1: 1})
    sage: D1
    1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2'), p=x0_1 + 1)
    sage: A.is_principal_divisor(D1)
    False
    sage: D2 = G({P1: 1, P3: 1})
    sage: A.is_principal_divisor(D2)
    True
    sage: A.principal_generator(D2)
    x0  

The FLIR structure also provides algorithms for atoms and factorizations.
For an element `f` of the algebra, one can compute its atoms::

    sage: x0, x1, x2 = A.gens()
    sage: f = x1 + 1
    sage: f.atoms()
    [x0, x2, (x1 + 1)/x2, (x1 + 1)/x0]
    sage: f.factor()
    [[((x1 + 1)/x0, 1), (x0, 1)], [((x1 + 1)/x2, 1), (x2, 1)]]

REFERENCES:

- [...] Mara Pompili and Daniel Smertnig, *Factoriality and Class
  Groups of Upper Cluster Algebras and Finite Laurent Intersection Rings:
  A Computational Approach*, 2026. :arxiv:`2601.07520`.

AUTHORS:

- Mara Pompili (2026-06): initial version
- Daniel Smertnig (2026-06): initial version
"""

from copy import copy
from collections.abc import Sequence

from sage.algebras.cluster_algebra import (
    ClusterAlgebra,
    ClusterAlgebraElement,
)
from sage.algebras.flir import FLIR, FLIRChart, FLIRElement
from sage.arith.misc import gcd
from sage.combinat.subset import Subsets
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import Matrix
from sage.misc.classcall_metaclass import typecall
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import CommutativeAlgebraElement

# ============================================================
# Banff algorithm helpers
# ============================================================

def is_seed_acyclic(B, allowed_directions: Sequence[int]) -> bool:
    r"""
    Check whether the principal part of ``B`` on the given mutable indices
    is acyclic.

    INPUT:

    - ``B`` -- an exchange matrix
    - ``allowed_directions`` -- indices of mutable directions to keep

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: B = Matrix([[0, 1], [-1, 0]])
        sage: is_seed_acyclic(B, [0, 1])
        True

        sage: B = Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
        sage: is_seed_acyclic(B, [0, 1, 2])
        False
    """
    P = B[allowed_directions, allowed_directions]
    dg = DiGraph(P.apply_map(lambda x: ZZ.zero() if x <= 0 else ZZ.one()))
    return dg.is_directed_acyclic()



def find_partner_sets(A, allowed_directions: Sequence[int]):
    r"""
    Compute partner sets among the allowed directions.

    Two directions(indices) are partners if their exchange polynomials have a
    non-constant common divisor.

    EXAMPLES::
        sage: B = Matrix([[0, 1], [-1, 0]])
        sage: A = ClusterAlgebra(B, scalars=QQ)
        sage: find_partner_sets(A, [0, 1])
        [(0,), (1,)]
        sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
        sage: A = ClusterAlgebra(B, scalars=QQ)
        sage: find_partner_sets(A, [0, 1, 2])
        [(0, 2), (1,)]
    """
    B = A.b_matrix()
    n = B.nrows()
    R = PolynomialRing(A.base_ring(), n, names='x')
    exch_polys = [R(1)] * n
    gens = R.gens()

    for i in allowed_directions:
        f = (
            prod(gens[h] ** B[h, i] for h in range(n) if B[h, i] > 0)
            + prod(gens[h] ** (-B[h, i]) for h in range(n) if B[h, i] < 0)
        )
        exch_polys[i] = f

    partner_sets = []
    unassigned = set(allowed_directions)

    while unassigned:
        i = unassigned.pop()
        current_set = {i}
        partners = set()
        f_i = exch_polys[i]

        for j in list(unassigned):
            f_j = exch_polys[j]
            common_divisor = gcd(f_i, f_j)
            if not common_divisor.is_unit():
                partners.add(j)

        current_set.update(partners)
        partner_sets.append(tuple(sorted(current_set)))
        unassigned -= partners

    return partner_sets


def find_sink_or_source_covering_pair(B, allowed_directions: Sequence[int]):
    r"""
    Find a sink or source in the principal part.

    Returns a tuple ``(i, j, "sink")`` or ``(i, j, "source")``, where
    ``i`` is a sink/source and ``j`` is a neighboring vertex witnessing it.
    Returns ``None`` if no such vertex exists.

    EXAMPLES::
        sage: B = Matrix([[0, 1], [-1, 0]])
        sage: find_sink_or_source_covering_pair(B, [0, 1])
        (0, 1, 'source')
        sage: B = Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
        sage: find_sink_or_source_covering_pair(B, [0, 1, 2]) is None
        True
    """
    for i in allowed_directions:
        row = [B[i, j] for j in allowed_directions]

        # source: all >=0 and at least one nonzero; pick j with B[i,j] > 0
        if all(x >= 0 for x in row) and any(x != 0 for x in row):
            for j in allowed_directions:
                if B[i, j] > 0:
                    return (i, j, "source")

        # sink: all <=0 and at least one nonzero; pick j with B[i,j] < 0
        if all(x <= 0 for x in row) and any(x != 0 for x in row):
            for j in allowed_directions:
                if B[i, j] < 0:
                    return (i, j, "sink")

    return None


class ClusterAlgebraChart:
    r"""
    Stores a chart algebra A_chart plus morphisms between fraction fields.
    
    - seed: the seed in the original algebra that defines the chart
    - chart: a Sage ClusterAlgebra whose initial seed corresponds to seed
    - to_chart: map from base fraction field -> chart fraction field
    - from_chart: inverse map chart fraction field -> base fraction field
    
    EXAMPLES:
    
    Build a rank-2 acyclic cluster algebra and a chart from one mutation
    step away from it::
    
        sage: B = matrix([[0, 1], [-1, 0]])
        sage: A = ClusterAlgebra(B)
        sage: seed = A.initial_seed()
        sage: seed.mutate(0)
        sage: A_chart = ClusterAlgebra(seed.b_matrix(),
        ....:                          cluster_variable_prefix='x0_')
        sage: chart = ClusterAlgebraChart.from_pair(seed, A_chart, [0, 1])
        sage: chart.chart is A_chart
        True
        sage: chart.seed is seed
        True
        sage: chart.allowed_directions
        [0, 1]
        sage: chart
        ClusterAlgebraChart(
          seed: The seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring obtained from the initial by mutating in direction 0
          chart: A Cluster Algebra with cluster variables x0_0, x0_1 and no coefficients over Integer Ring
          allowed directions: [0, 1]
          )
 
    The two morphisms are inverse to each other on generators, up to the
    identification of the ambient fraction fields::
 
        sage: F = A.ambient().fraction_field()
        sage: Fp = A_chart.ambient().fraction_field()
        sage: chart.to_chart.domain() is F
        True
        sage: chart.from_chart.domain() is Fp
        True
    """
    @classmethod
    def from_pair(cls, seed_chart, A_chart, allowed_directions: Sequence[int]):
        r"""
        Construct a :class:`ClusterAlgebraChart` from a seed of the base
        algebra and an independently-built :class:`ClusterAlgebra` sharing
        the same exchange matrix.
 
        INPUT:
 
        - ``seed_chart`` -- a seed of the *base* cluster algebra, whose
          ``b_matrix()`` agrees with ``A_chart.b_matrix()`` at the initial
          seed
        - ``A_chart`` -- a :class:`ClusterAlgebra` built from that exchange
          matrix (typically with distinct variable names, so Sage does not
          identify it with the base algebra)
        - ``allowed_directions`` -- the mutation directions still active in
          this branch of the recursion
 
        OUTPUT: a :class:`ClusterAlgebraChart` wrapping ``seed_chart`` and
        ``A_chart`` together with the two fraction-field morphisms between
        them
 
        EXAMPLES::
 
            sage: B = matrix([[0, 1], [-1, 0]])
            sage: A = ClusterAlgebra(B)
            sage: seed = A.initial_seed()
            sage: A_chart = ClusterAlgebra(B, cluster_variable_prefix='x0_')
            sage: chart = ClusterAlgebraChart.from_pair(seed, A_chart, [0, 1])
            sage: chart
            ClusterAlgebraChart(
              seed: The initial seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
              chart: A Cluster Algebra with cluster variables x0_0, x0_1 and no coefficients over Integer Ring
              allowed directions: [0, 1]
              )
 
        A chart built at a mutated seed still round-trips its own initial
        cluster variables through ``to_chart`` and back through
        ``from_chart``::
 
            sage: seed2 = A.initial_seed()
            sage: seed2.mutate(0)
            sage: A_chart2 = ClusterAlgebra(seed2.b_matrix(),
            ....:                           cluster_variable_prefix='y0_')
            sage: chart2 = ClusterAlgebraChart.from_pair(seed2, A_chart2, [0, 1])
            sage: v = chart2.from_chart(A_chart2.ambient().fraction_field().gens()[0])
            sage: chart2.to_chart(v) == A_chart2.ambient().fraction_field().gens()[0]
            True
        """
        A_initial = seed_chart.parent()
        n = A_initial.rank()

        L = A_initial.ambient()
        F = L.fraction_field()
        L_prime = A_chart.ambient()
        F_prime = L_prime.fraction_field()

        # Compute cluster variables of seed_chart directly in the ambient ring,
        # bypassing A_initial.retract()/__call__() (and therefore its custom
        # _element_constructor_/_coerce_map_from_ machinery). We only need the
        # raw Laurent expression here, not a checked algebra element - and
        # calling retract() reenters coercion machinery that may not yet be
        # safely usable while charts are still being built.
        def _raw_cluster_variable(seed, i):
            g_vector = tuple(seed.g_vector(i))
            parent = seed.parent()
            Fpoly = parent.F_polynomial(g_vector)
            F_std = Fpoly.subs(parent._yhat)
            g_mon = prod(parent.ambient().gen(k) ** g_vector[k] for k in range(parent.rank()))
            F_trop = parent.ambient()(Fpoly.subs(parent._y))._fraction_pair()[1]
            return parent.ambient()(g_mon * F_std * F_trop) if False else (g_mon * F_std * F_trop)

        mapping = [F(_raw_cluster_variable(seed_chart, i)) for i in range(n)]
        from_chart = F_prime.hom(mapping)

        # Map from base field -> chart field
        reverse_seed = A_chart.initial_seed()
        for i in reversed(seed_chart.path_from_initial_seed()):
            reverse_seed.mutate(i)

        # reverse_seed.parent() is a plain ClusterAlgebra (A_chart), not
        # BanffClusterAlgebra, so its own cluster_variable()/retract() is safe -
        # no override, no risk of reentering Banff/FLIR bootstrap logic.
        reverse_mapping = [F_prime(reverse_seed.cluster_variable(i)) for i in range(n)]
        to_chart = F.hom(reverse_mapping)

        return cls(seed_chart, A_chart, to_chart, from_chart, allowed_directions)

    def __init__(self, seed_chart, A_chart, to_chart, from_chart, allowed_directions: Sequence[int]):
        self.chart = A_chart
        self.to_chart = to_chart
        self.from_chart = from_chart
        self.seed = seed_chart
        self.allowed_directions = list(allowed_directions)

    @property
    def lp_chart(self):
        """The ambient Laurent polynomial ring of the chart."""
        return self.chart.ambient()
    
    def __repr__(self):
        return (
            "ClusterAlgebraChart(\n"
            "  seed: {}\n"
            "  chart: {}\n"
            "  allowed directions: {}\n"
            ")"
        ).format(
            self.seed,
            self.chart,
            self.allowed_directions,
        )
    
    


def _freeze_and_continue(A, allowed_directions: list[int], current_seed, counter: dict, max_steps=None):
    """
    Recursively walk mutations of ``current_seed`` within
    ``allowed_directions`` until every branch reaches an acyclic
    exchange matrix, freezing (removing from ``allowed_directions``) one
    direction of a sink/source-covering pair at each cyclic step.
 
    INPUT:
 
    - ``A`` -- the ambient :class:`ClusterAlgebra` being explored
    - ``allowed_directions`` -- list of directions still mutable in this
      branch
    - ``current_seed`` -- the seed the recursion currently sits at
    - ``counter`` -- a mutable ``dict`` with keys ``'count'`` (recursion
      steps taken, for the ``max_steps`` guard) and ``'chart_num'`` (used
      to generate distinct variable-name prefixes per discovered chart)
    - ``max_steps`` -- optional recursion budget; raises :class:`ValueError`
      if exceeded
 
    OUTPUT: a list of :class:`ClusterAlgebraChart` objects, one per acyclic
    chart discovered along every branch of the recursion

    """
    counter['count'] += 1

    if max_steps is not None and counter['count'] > max_steps:
        raise ValueError("Max steps exceeded in Banff recursion.")

    matrix = current_seed.b_matrix()

    # base case: acyclic principal part
    if is_seed_acyclic(matrix, allowed_directions):
        c = counter['chart_num']
        counter['chart_num'] += 1

        # important: distinct variable names per chart
        A_chart = ClusterAlgebra(  
            matrix,
            scalars=A.base_ring(),
            cluster_variable_prefix=f"x{c}_",
        )

        acyclic_chart = ClusterAlgebraChart.from_pair(
            seed_chart=current_seed,
            A_chart=A_chart,
            allowed_directions=allowed_directions,
        )
        return [acyclic_chart]

    # recursive case
    old_current_seed = A.current_seed()
    A.set_current_seed(current_seed)

    for seed in A.seeds(allowed_directions=allowed_directions, from_current_seed=True):
        B = seed.b_matrix()
        covering_pair = find_sink_or_source_covering_pair(B, allowed_directions)

        if covering_pair:
            i, j, _cover_type = covering_pair

            allowed_i = [d for d in allowed_directions if d != i]
            results_i = _freeze_and_continue(A, allowed_i, seed, counter, max_steps=max_steps)

            allowed_j = [d for d in allowed_directions if d != j]
            results_j = _freeze_and_continue(A, allowed_j, seed, counter, max_steps=max_steps)

            A.set_current_seed(old_current_seed)
            return results_i + results_j

    A.set_current_seed(old_current_seed)
    return []


def _banff_algorithm(A, max_steps=None):
    """
    Compute a set of acyclic charts via Banff recursion.
 
    INPUT:
 
    - ``A`` -- a :class:`ClusterAlgebra`
    - ``max_steps`` -- optional recursion budget passed through to
      :func:`_freeze_and_continue`
 
    OUTPUT: a list of :class:`ClusterAlgebraChart`, covering ``A`` by
    acyclic charts

    """
    B = A.b_matrix()
    allowed_directions = list(range(B.ncols()))
    counter = {'count': 0, 'chart_num': 0}
    return _freeze_and_continue(A, allowed_directions, A.initial_seed(), counter, max_steps=max_steps)


def _FLIR_charts_for_acyclic(acyclic_chart: ClusterAlgebraChart):
    """
    Given an acyclic chart, refine it to FLIR charts using partner sets.
 
    For every nonempty subset ``J`` of every partner set ``S`` found by
    :func:`find_partner_sets`, this mutates the chart's seed along ``J``
    and builds a new :class:`ClusterAlgebraChart` at that mutated seed
    (with variable names decorated by ``"p"`` on the mutated indices so
    Sage keeps the charts distinct).
 
    INPUT:
 
    - ``acyclic_chart`` -- a :class:`ClusterAlgebraChart` whose exchange
      matrix is acyclic on ``acyclic_chart.allowed_directions``
 
    OUTPUT: a list of :class:`ClusterAlgebraChart`, starting with
    ``acyclic_chart`` itself followed by one chart per nonempty subset of
    each partner set.
    """
    A = acyclic_chart.chart
    B = A.b_matrix()
    assert is_seed_acyclic(B, acyclic_chart.allowed_directions)

    charts = [acyclic_chart]

    partner_sets = find_partner_sets(A, acyclic_chart.allowed_directions)

    for S in partner_sets:
        for J in Subsets(S):
            if not J:
                continue

            seed_from_start = copy(acyclic_chart.seed)
            for j in J:
                seed_from_start.mutate(j)

            # rename variables so Sage doesn't identify different charts
            var_names = [
                name + "p" if k in J else name
                for (k, name) in enumerate(A.variable_names())
            ]

            A_prime = ClusterAlgebra(
                seed_from_start.b_matrix(),
                scalars=A.base_ring(),
                cluster_variable_names=var_names,
            )
            chart = ClusterAlgebraChart.from_pair(
                seed_from_start,
                A_prime,
                allowed_directions=acyclic_chart.allowed_directions,
            )
            charts.append(chart)

    return charts



def _system_FLIR_charts_for_acyclic(A, acyclic_charts: list[ClusterAlgebraChart]):
    """
    Apply :func:`_FLIR_charts_for_acyclic` to every chart in
    ``acyclic_charts`` and concatenate the results.
 
    INPUT:
 
    - ``A`` -- the ambient :class:`ClusterAlgebra` (unused by the current
      implementation but kept for API symmetry / future use)
    - ``acyclic_charts`` -- a list of acyclic :class:`ClusterAlgebraChart`
      objects, e.g. as returned by :func:`_banff_algorithm`
 
    OUTPUT: the concatenation of ``_FLIR_charts_for_acyclic(chart)`` over
    every ``chart`` in ``acyclic_charts``
    """
    charts: list[ClusterAlgebraChart] = []
    for chart in acyclic_charts:
        charts += _FLIR_charts_for_acyclic(chart)
    return charts



# ============================================================
# BanffClusterAlgebraElement
# ============================================================

class BanffClusterElement(FLIRElement, ClusterAlgebraElement):
    r"""
    Element of a BanffClusterAlgebra, represented by an element of the base fraction field.

    For expressions with a scalar on the left, such as ``1/f``, users should
    explicitly coerce the scalar into the parent and write ``A(1)/f``.
    In particular, ``1/f`` may fail because of Sage coercion behavior for
    custom parents.

    Scalar constants should be coerced explicitly when dividing by algebra elements::

    EXAMPLES::
        sage: B = Matrix([[0, 1], [-1, 0]])
        sage: A = BanffClusterAlgebra(B)
        sage: x0, x1 = A.gens()
        sage: isinstance(x0, BanffClusterElement)
        True
        sage: isinstance(x0, ClusterAlgebraElement)
        True
        sage: isinstance(x0, FLIRElement)
        True
        sage: x0*x1
        x0*x1
        sage: x0 + x1
        x0 + x1
        sage: (x0 + x1) + (x0 - x1)
        2*x0
        sage: x0 - x1
        x0 - x1

    Adding/subtracting a plain integer coerces it into the algebra first::

        sage: x0 + 1
        x0 + 1
        sage: x0**2
        x0^2
        sage: x1**0
        1

    Due to a current implementation issue, left scalar multiplication by Python
    integers, such as `2*x0`, may cause a segmentation fault. As a workaround,
    the scalar should first be coerced into the parent algebra::

        sage: A(2)*x0
        2*x0
        sage: x1*A(1/2)
        1/2*x1

    This issue should be fixed so that scalar multiplication works directly.


    Division is only defined if the quotient lies in the algebra. For example, ``x0/x1`` is not in the algebra::

        sage: x0/x1
        Traceback (most recent call last):
        ...
        ValueError:  Not Laurent in chart ('x0_0', 'x0_1p').
        Substituted: x0_0*x0_1p/(x0_0 + 1)
    
    """
        
    
    def __init__(self, parent, f, check=False):
        F = parent.ambient().fraction_field()
        L = parent.ambient()

        f = F(f)

        if check:
            parent._check_membership(f)

        CommutativeAlgebraElement.__init__(self, parent)

        self._f = f
        self._value = L(f)
        self.value = self._value
        #FLIRElement.__init__(self, parent, L(f), check=check)

        #ClusterAlgebraElement.__init__(self, parent, L(f))

    def lift(self):
        return self._value

    def _repr_(self):
        return repr(self._f)

    @property
    def f(self):
        r"""
        Return the underlying representative in the ambient fraction
        field.
 
        EXAMPLES::
 
            sage: B = Matrix([[0, 1], [-1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x0, x1 = A.gens()
            sage: (x0 * x1).f
            x0*x1
            sage: (x0 * x1).f.parent() is A.ambient().fraction_field()
            True
            sage: (x0 + x1).f == x0.f + x1.f
            True
            sage: x0.f.parent() is x0.parent()
            False
        """
        return self._f

    # --- arithmetic operations  ---
    def _add_(self, other):
        A = self.parent()
        if not isinstance(other, BanffClusterElement) or other.parent() is not A:
            other = A(other)
        return A(self._f + other._f, check=False)

    def _sub_(self, other):
        A = self.parent()
        if not isinstance(other, BanffClusterElement) or other.parent() is not A:
            other = A(other)
        return A(self._f - other._f, check=False)

    def _mul_(self, other):
        A = self.parent()
        if not isinstance(other, BanffClusterElement) or other.parent() is not A:
            other = A(other)
        return A(self._f * other._f, check=False)
    
    def _lmul_(self, c):
        return self.parent()(c * self._f, check=False)

    def _rmul_(self, c):
        return self.parent()(self._f * c, check=False)

    def _neg_(self):
        return self.parent()(-self._f, check=False)

    def __pow__(self, n):
        return self.parent()(self._f ** int(n), check=False)

    def __truediv__(self, other):
        """
        Division is only defined if the quotient lies in the algebra.
        """
        A = self.parent()
        if not isinstance(other, BanffClusterElement) or other.parent() is not A:
            other = A(other)  # try coercion

        if other._f == 0:
            raise ZeroDivisionError("Division by zero.")

        q = self._f / other._f  # computed in the ambient fraction field
        if not (
            getattr(A, "_flir_initializing", False)
        ):
            A._check_membership(q)
        return A(q, check=False)
    
    
# ============================================================
# BanffClusterAlgebra 
# ============================================================

class BanffClusterAlgebra(ClusterAlgebra, FLIR):
    r"""
    A Banff cluster algebra is a cluster algebra which can be covered by
    acyclic cluster charts using the Banff algorithm.  This class combines
    the cluster algebra structure with the FLIR machinery, so that one can
    compute FLIR charts, divisors, class group data, and factorization-related
    information.

    INPUT:

    - ``data`` -- an exchange matrix or quiver accepted by ``ClusterAlgebra``
    - ``scalars`` -- optional base ring; default is ``QQ``
    - ``term_order`` -- optional monomial order used for FLIR charts
    - ``check_Banff`` -- boolean, default ``True``; whether to certify the
      Banff property during construction
    - ``max_steps`` -- optional recursion budget for the Banff algorithm
    - ``flir`` -- boolean, default ``False``; whether to initialize the FLIR
      structure immediately

    EXAMPLES:

    Construct the rank two Banff cluster algebra::

        sage: B = Matrix([[0, 1], [-1, 0]])
        sage: A = BanffClusterAlgebra(B)
        sage: A
         A Banff Cluster Algebra with initial cluster variables x0, x1 over Rational Field.

    The Banff algorithm is a semi-algorithm: in general, it is not guaranteed to
    terminate. A cluster algebra is called a Banff cluster algebra if the algorithm
    terminates and produces an acyclic cover.

    By default, construction attempts to certify the Banff property::

        sage: B = Matrix([[0, 1, -1, 1], [-1, 0, 1, 1], [1, -1, 0, 1], [-1, -1, -1, 0]])
        sage: A = BanffClusterAlgebra(B)
        sage: A.is_banff()
        True

    The optional `max_steps` parameter bounds the search. It can be used to test
    whether an acyclic cover is found within a prescribed number of steps. Failure
    to find a cover within this bound does not imply that the cluster algebra is
    not Banff; it only means that the semi-algorithm was inconclusive within the
    given budget.

        sage: B = Matrix([[0, 2, -2], [-2, 0, 2], [2, -2, 0]])
        sage: A = BanffClusterAlgebra(B, check_Banff=False) # long time
        sage: A.is_banff(max_steps=5) # long time 
        RuntimeError: Banff semi-algorithm inconclusive within given budget. Try increasing max_steps.

    A Banff Cluster Algebra is a Cluster Algebra and a FLIR::

        sage: B = Matrix([[0, 1, -1, 1], [-1, 0, 1, 1], [1, -1, 0, 1], [-1, -1, -1, 0]])
        sage: A = BanffClusterAlgebra(B)
        sage: isinstance(A, ClusterAlgebra)
        True
        sage: isinstance(A, FLIR)
        True

    Its elements are cluster algebra elements, Banff cluster elements, and
    FLIR elements::

        sage: x0, x1, x2, x3 = A.gens()
        sage: isinstance(x0, BanffClusterElement)
        True
        sage: isinstance(x0, ClusterAlgebraElement)
        True
        sage: isinstance(x0, FLIRElement)
        True

    Basic arithmetic works inside the algebra::

        sage: x0*x1
        x0*x1
        sage: x0 + x1
        x0 + x1
        sage: x0 - x1
        x0 - x1
        sage: (x0 + x1) + (x0 - x1)
        2*x0
        sage: x0 + 1
        x0 + 1

    Scalar multiplication should currently be performed by first coercing the
    scalar into the parent algebra.  Direct left/right multiplication by a Python
    integer, such as ``2*x0``, is a known implementation issue and may cause
    a segmentation fault::

        sage: A(2)*x0
        2*x0

    One can also construct the FLIR data immediately::

        sage: A = BanffClusterAlgebra(B, flir=True)
        sage: A._flir_initialized
        True
    """
    Element = BanffClusterElement  

    @staticmethod
    def __classcall__(cls, data, *args, **kwargs):
        
        kwargs = dict(kwargs)

        # things ClusterAlgebra.__init__ expects
        kwargs.setdefault("scalars", QQ)
        kwargs.setdefault("cluster_variable_prefix", "x")
        kwargs.setdefault("next_free_index", 0)

        # names for mutable cluster variables
        if "cluster_variable_names" not in kwargs:
            if hasattr(data, "ncols"):
                n = data.ncols()
            else:
                n = len(data[0])   # only as fallback if data is list-like
            prefix = kwargs["cluster_variable_prefix"]
            kwargs["cluster_variable_names"] = tuple(f"{prefix}{i}" for i in range(n))

        # coefficients / frozen variables
        if "coefficient_names" not in kwargs:
            if hasattr(data, "nrows") and hasattr(data, "ncols"):
                m = data.nrows() - data.ncols()
            else:
                m = 0
            kwargs["coefficient_names"] = tuple(f"y{i}" for i in range(m))
        if "check_Banff" in kwargs:
            kwargs["check_Banff"] = bool(kwargs["check_Banff"])
        if "flir_charts_recompute" in kwargs:
            kwargs["flir_charts_recompute"] = bool(kwargs["flir_charts_recompute"])

        return typecall(cls, data, **kwargs)

    def __init__(
        self,
        data,
        *,
        scalars=QQ,
        term_order="lex",
        check_Banff=True,
        max_steps=None,
        flir=False,
        **kwargs,
    ):
        kwargs = dict(kwargs)

        kwargs["scalars"] = scalars
        kwargs.setdefault("cluster_variable_prefix", "x")
        kwargs.setdefault("next_free_index", 0)

        if "cluster_variable_names" not in kwargs:
            n = data.ncols()
            prefix = kwargs["cluster_variable_prefix"]
            kwargs["cluster_variable_names"] = tuple(f"{prefix}{i}" for i in range(n))

        if "coefficient_names" not in kwargs:
            m = data.nrows() - data.ncols()
            kwargs["coefficient_names"] = tuple(f"y{i}" for i in range(m))

        super().__init__(data, **kwargs)
        self._populate_coercion_lists_(coerce_list=[self.base_ring()])

        self._banff_ready = False
        self.term_order = term_order
        
        self._cached_flir_charts = None
        self._banff_acyclic_charts_cache = None
        self._banff_certified = None

        self._flir_initialized = False
        self._flir_initializing = False
        self._base_chart = None
        self.charts = None
        self.n = None

        self._extra_primes_cache = None
        self._class_data_cache = None
        self._div_group_cache = None


        self.Element = BanffClusterElement
        self.element_class = BanffClusterElement # type: ignore[assignment]

        if check_Banff:
            self.is_banff(max_steps=max_steps)
        
        self._banff_ready = True

        if flir:
            self._ensure_flir_initialized()

    
    def _ensure_flir_initialized(self, recompute=False):
        if getattr(self, "_flir_initialized", False) and not recompute:
            return

        if getattr(self, "_flir_initializing", False):
            return  # or raise a clearer internal error

        self._flir_initializing = True
        try:
            base_chart = FLIRChart.base(
                self.base_ring(),
                self.variable_names(),
                term_order=self.term_order,
            )
            #NOTE: am i counting base chart double?
            flir_charts = self._build_flir_charts(recompute=recompute)
            all_charts = [base_chart] + list(flir_charts or [])

            self._init_flir_structure(
                base_chart,
                all_charts,
                compute_base_to_charts=True
            )
            self._flir_initialized = True
        finally:
            self._flir_initializing = False

    def is_banff(self, *, max_steps=None, recompute=False) -> bool:
        r"""
        Try to certify that ``self`` is Banff.

        This runs the Banff semi-algorithm.  If a Banff cover is found, the
        corresponding acyclic charts are cached and the method returns ``True``.
        If the algorithm is inconclusive within the given budget, a
        ``RuntimeError`` is raised.

        INPUT:

        - ``max_steps`` -- optional recursion budget
        - ``recompute`` -- boolean, default ``False``; whether to ignore the
        cached result and recompute the Banff cover

        OUTPUT:

        Boolean.  Returns ``True`` if a Banff cover is found.
        """
        if not recompute and getattr(self, "_banff_certified", None) is True:
            return True

        try:
            charts = _banff_algorithm(self, max_steps=max_steps)
        except RuntimeError:
            raise

        if charts:
            self._banff_acyclic_charts_cache = charts
            self._banff_certified = True
            return True
        else:
            raise RuntimeError("Banff semi-algorithm inconclusive within given budget. Try increasing max_steps.")
            
    
    def _element_constructor_(self, x, check=True):
        F = self.ambient().fraction_field()

        if isinstance(x, BanffClusterElement) and x.parent() is self:
            return x
        try:
            fa = F(x)
        except Exception:
            try:
                fa = F(self.ambient()(x))
            except Exception as e:
                raise TypeError(f"Cannot coerce {x!r} into {self}.") from e

        if not getattr(self, "_banff_ready", False):
            return BanffClusterElement(self, fa, check=False)

        if (
            check
            and not getattr(self, "_flir_initializing", False)
        ):
            self._check_membership(fa)

        return BanffClusterElement(self, fa, check=False)



    def _check_membership(self, f):
        r"""
        Check whether ``f`` lies in this Banff cluster algebra.

        The test verifies that ``f`` is Laurent in every FLIR chart associated to
        the Banff cover.  This is the expensive membership test used when coercing
        user-supplied elements into the algebra.

        INPUT:

        - ``f`` -- an element of the ambient fraction field

        EXAMPLES::

            sage: B = Matrix([[0, 1, -1, 1], [-1, 0, 1, 1], [1, -1, 0, 1], [-1, -1, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x0, x1, x2, x3 = A.gens()
            sage: A._check_membership((x0 + x1 + x2)/x3)
            Traceback (most recent call last):
            ...
            ValueError: Not Laurent in chart ('x2_0', 'x2_1', 'x2_2', 'x2_3p').
            Substituted: (x2_0*x2_3p + x2_1*x2_3p + x2_2*x2_3p)/(x2_0*x2_1*x2_2 + 1)


        An element of the ambient fraction field which is not Laurent in all charts
        is rejected::

            sage: F = A.ambient().fraction_field()
            sage: f = F(1)/(F(x0)*F(x1))
            sage: A._check_membership(f)
            Traceback (most recent call last):
            ...
            ValueError: Not Laurent in chart ('x0_0', 'x0_1p', 'x0_2', 'x0_3').
            Substituted: x0_1p/(x0_0*x0_2 + x0_0*x0_3^2)

        """
        if not getattr(self, "_banff_ready", False):
            print("Warning: BanffClusterAlgebra not yet ready; skipping membership check.")
            return
        

        self._ensure_flir_initialized()
        base = self.charts[0] 
        g = base.L(f)
        
        for ch in self.charts:
            expr = ch._substitute_from_base(g)
            try:
                ch.L(expr)
            except Exception:
                raise ValueError(
                    f"Not Laurent in chart {ch.var_names}.\nSubstituted: {expr}"
                )

    def base_extend(self, R):
        """
        Minimal compatibility hook for Sage coercion machinery.
        """
        try:
            if self.base_ring().has_coerce_map_from(R):
                return self
        except TypeError:
            pass

        try:
            parent_R = R.parent()
            if self.base_ring().has_coerce_map_from(parent_R):
                return self
        except Exception:
            pass

        raise NotImplementedError(
            f"Base extension from {self.base_ring()} to {R} is not implemented."
        )
    
    def _coerce_map_from_(self, S): # pyright: ignore[reportIncompatibleMethodOverride]
        if S is self:
            return True
        if S is self.base_ring():
            return True
        if self.base_ring().has_coerce_map_from(S):
            return True
        return ClusterAlgebra._coerce_map_from_(self, S) or FLIR._coerce_map_from_(self, S)
    
    # --- chart construction ---
    
    def _build_flir_charts(self, recompute=False):
        """
        Construct a system of FLIRChart objects:
          - base chart = initial Laurent ring
          - additional charts = obtained from the Banff algorithm and refinements
        """

        if (not recompute) and getattr(self, "_cached_flir_charts", None) is not None:
            return self._cached_flir_charts
        
        K = self.base_ring()

        base_chart = FLIRChart(
            K,
            self.variable_names(),
            term_order=self.term_order,
        )
        
        F_base = base_chart.F

        if recompute or self._banff_acyclic_charts_cache is None:
            self.is_banff(recompute=True)
        acyclic = self._banff_acyclic_charts_cache
        if acyclic is None:
            acyclic = []
        lp_like = _system_FLIR_charts_for_acyclic(self, acyclic)
        charts = []

        for c in lp_like:
            Lc = c.lp_chart
            var_names = tuple(str(v) for v in Lc.gens())

            ch = FLIRChart(
                K,
                var_names,
                term_order=base_chart.term_order,
                base_fraction_field=F_base,
            )

            F_chart = Lc.fraction_field()

            # base -> chart
            base_gens = list(F_base.gens())
            ch.base_to_this = [F_chart(c.to_chart(xi)) for xi in base_gens]

            # chart -> base
            chart_gens = list(F_chart.gens())
            ch.this_to_base = [F_base(c.from_chart(yj)) for yj in chart_gens]

            charts.append(ch)
        

        self._cached_flir_charts = charts
        return charts
    
    
    def _banff_algorithm_with_generators(self, current_seed, allowed_directions):
        """ Find a covering pair, freeze at the two indices, and recursively continue. """

        # If the principal part of the seed is acyclic, we are done. We can explicitly
        # get a generating set by taking the current seed and mutating once in each direction
        B = current_seed.b_matrix()
        
        if is_seed_acyclic(B, allowed_directions):
            F = self.ambient().fraction_field()
            gens = [BanffClusterElement(self, F(x), check=False) for x in current_seed.cluster_variables()]
            for k in allowed_directions:
                current_seed.mutate(k)
                xk_prime = current_seed.cluster_variables()[k]
                gens.append(BanffClusterElement(self, F(xk_prime), check=False))
                gens.append(xk_prime)
                current_seed.mutate(k)
            return gens
            

        # Otherwise, we mutate through all seeds until we find a covering pair
        # By definition, this eventually terminates if A is defined by a Banff quiver;
        # otherwise it may end up in an infinite loop

        A = self
        old_current_seed = A.current_seed()
        A.set_current_seed(current_seed)
       
        seed_count = 0
        
        for seed in A.seeds(allowed_directions=allowed_directions,
                            from_current_seed=True):
            
            if isinstance(seed, KeyboardInterrupt):
                continue
            B = seed.b_matrix()
            covering_pair = find_sink_or_source_covering_pair(B, allowed_directions)

            if covering_pair: # There is a covering pair, we can freeze and recurse
                i, j, cover_type = covering_pair

                # Freeze x_i
                allowed_directions_i = [ d for d in allowed_directions if d != i ]  
                generators_i = self._banff_algorithm_with_generators(seed, allowed_directions_i)  
                # Freeze x_j
                allowed_directions_j = [ d for d in allowed_directions if d != j ]
                generators_j = self._banff_algorithm_with_generators(seed, allowed_directions_j)
                
                # add extra generators ensuring that the presentations patch together correctly
                seed.mutate(i)
                xi_prime = seed.cluster_variables()[i]
                seed.mutate(i)

                cluster = seed.cluster_variables()
                if cover_type == "source":
                    sign = 1
                else:
                    sign = -1
                
                a = A(-1)
                u = A(1)
                for k in range(B.nrows()):
                    x_k = cluster[k]
                    b_ik = B[i,k]
                    if b_ik * sign < 0:
                        u *= x_k**abs(b_ik)
                    elif b_ik * sign > 0:
                        if k != j:
                            a *= x_k**abs(b_ik)
                        else:
                            a *= x_k**(abs(b_ik) - 1)           
                assert cluster[i]*xi_prime + a*cluster[j] == u
                patching_gens = [ cluster[i], cluster[j], xi_prime, a, u ]
                
                A.set_current_seed(old_current_seed)
                return generators_i + generators_j + patching_gens
        A.set_current_seed(old_current_seed)
        return []  # If no sink-source is found, return an empty list

    def generators(self):
        r"""
        Banff Cluster Algebras are finitely generated.  This method returns a list of generators.

        The method recursively applies the Banff algorithm, collects generators
        from the resulting acyclic charts, removes constants and monomials, and
        deduplicates the result deterministically.

        OUTPUT:

        A list of Laurent polynomial generators.

        EXAMPLES::

            sage: B = Matrix([[0, 1], [-1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.generators()
            [x0^-1*x1 + x0^-1, x0*x1^-1 + x1^-1]

        The output is deterministic::

            sage: A.generators() == A.generators()
            True
        """

        A = self
        B=A.b_matrix()
        allowed_directions = list(range(B.ncols()))  # Initially all directions are allowed

        gens = self._banff_algorithm_with_generators(A.initial_seed(), allowed_directions)

        # Throw out constants and monomials (they are products of elemens of the initial seed,
        # which we will add right below)
        L = A.ambient()

        gens = [ L(g) for g in gens if not L(g).is_monomial() and not L(g).is_constant() ]
        
        # Deduplicate deterministically: dedupe by a canonical, content-derived
        # key, then sort by that same key so the result is reproducible
        #regardless of the order _banff_algorithm_with_generators happened to produce them in.
        seen = {}
        for g in gens:
            key = tuple(sorted(g.dict().items()))
            seen[key] = g  # last write wins; fine since duplicates are equal

        gens = [seen[key] for key in sorted(seen)]

        return gens
    
    def presentation(self):
        r"""
        Compute a presentation of this Banff cluster algebra.

        OUTPUT:

        A triple ``(Q, R, I)``, where ``R`` is a polynomial ring, ``I`` is an
        ideal of relations, and ``Q = R.quotient(I)``.

        EXAMPLES::

            sage: B = Matrix([[0, 1, -1, 1], [-1, 0, 1, 1], [1, -1, 0, 1], [-1, -1, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: Q, R, I = A.presentation()
            sage: R.base_ring() is QQ
            True
            sage: Q
             Quotient of Multivariate Polynomial Ring in 
             ...
             X1^2*T0*T3 + X1^2*T5*T6 + T0*T3*T5*T6 + T5^2*T6^2 - X1*T3*T4 - X1*T4*T5, T4^4*T5 + X1*T0*T4^2 + T4^2*T5^2 + T0*T1*T5*T6 - T4^3 + X1*T0*T5 - T4*T5, X2*T2*T4*T5*T6 - T2*T3*T4*T5 + T1^2*T3*T6 + T2*T4^2 - X2*T5*T6 - T1*T6^2 + T3*T5 - T4, T2^2*T4^2*T5*T6 - T1^3*T3*T6 + T2^2*T5^2*T6 - T1*T2*T4^2 - T2*T4*T5*T6 + T1^2*T6^2 - T1*T2*T5 - T1*T3*T5 + T1*T4 + T5*T6, X2*T2^2*T5^2*T6^2 + T1^3*T3^2*T6 - T2^2*T3*T5^2*T6 + 2*T2^2*T4*T5*T6 + T2*T3*T4*T5*T6 - 2*T1^2*T3*T6^2 + T1*T2*T3*T5 + T1*T3^2*T5 + X2*T5*T6^2 + T1*T6^3 - T1*T2*T4 - 2*T2*T5*T6 - 2*T3*T5*T6 + T1, T2^3*T4*T5^2*T6^2 - T1^4*T3^2*T6 + T1*T2^2*T3*T5^2*T6 - 2*T1*T2^2*T4*T5*T6 + 2*T1^3*T3*T6^2 - 2*T2^2*T5^2*T6^2 - T1^2*T2*T3*T5 - T1^2*T3^2*T5 - T1^2*T6^3 + T1^2*T2*T4 + 3*T1*T2*T5*T6 + 2*T1*T3*T5*T6 - T5*T6^2 - T1^2, T2^4*T5^3*T6^3 - T1^5*T3^3*T6 + T1^2*T2^2*T3^2*T5^2*T6 + 3*T1^4*T3^2*T6^2 - 3*T1*T2^3*T5^2*T6^2 - 3*T1*T2^2*T3*T5^2*T6^2 - T1^3*T2*T3^2*T5 - T1^3*T3^3*T5 - 3*T1^3*T3*T6^3 + 2*T2^2*T5^2*T6^3 + 3*T1^2*T2^2*T5*T6 + 4*T1^2*T2*T3*T5*T6 + 3*T1^2*T3^2*T5*T6 + T1^2*T6^4 - 3*T1*T2*T5*T6^2 - 3*T1*T3*T5*T6^2 - T1^3*T2 - T1^3*T3 + T5*T6^3 + T1^2*T6)
        """

        generators = self.generators()
        L = self.ambient()
        scalars = L.base_ring()
        n = L.ngens()

        yvars = [ f"Y{i}" for i in range(n) ]
        xvars = [ f"X{i}" for i in range(n) ]

        R = PolynomialRing(scalars, xvars, order='degrevlex')
        I = R.ideal(0)

        numgens = len(generators)

        for j in range(numgens):
            f = generators[j]

            curgens  = list(R.gens())
            currels = I.gens()

            curgens = yvars + [f"T{j}"] + curgens

            S = PolynomialRing(scalars, curgens)

            Xsub = dict([ (L(f"x{i}"), S(f"X{i}")) for i in range(n) ])
            Ysub = dict([ (L(f"x{i}"), S(f"Y{i}")) for i in range(n) ])

            
            xyrels = [ S(f"X{i}*Y{i} - 1") for i in range (n) ]
            currels = [ S(r) for r in currels ]
            currels += xyrels
                
            # construct the denominator exponents of f
            dexp = [0]*(n+1)
            for e in f.exponents():
                for i in range(n):
                    if e[i] < 0:
                        dexp[i] = max(-e[i], dexp[i])

            d = prod([ L.gen(i) ** dexp[i] for i in range(n) ])
            g = f*d
            G=g.substitute(Xsub)
            D=d.substitute(Ysub)

            currels.append(G*D - S(f"T{j}"))
    
            Iprime = S.ideal(currels)
            
            J = Iprime.elimination_ideal([ S(f"Y{i}") for i in range(n) ])  

            R = PolynomialRing(scalars, list(R.gens()) + [f"T{j}"])
            I = R.ideal(J.gens())

        return R.quo(I), R, I
        
    # ---------------------------
    # methods from FLIR
    # --------------------------

    def extra_primes(self, recompute=False):
        r"""
        Return the list P1,...,Pr of height-one prime ideals that contain x1*...*xn.
        
        Notice that the height-one spectrum of self is given by P1,...,Pr 
        and the height-one sprectrum of self.ambient().

        EXAMPLES::
            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.extra_primes()
            [PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2'), p=x0_1 + 1),
             PrimeDivisor(chart=('x0_0', 'x0_1', 'x0_2p'), p=x0_1 + 1),
             PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2p'), p=x0_1 + 1),
             PrimeDivisor(chart=('x0_0', 'x0_1p', 'x0_2'), p=x0_0 + x0_2)]
        """
        self._ensure_flir_initialized(recompute=recompute)
        return FLIR.extra_primes(self, recompute=recompute)

    def divisor_group(self):
        r""" Return the divisor group of this Banff cluster algebra. 
        The divisor group is the free abelian group generated by the height-one prime divisors of ``self``. 
        It is represented by an instance of :class:`FLIRDivisorGroup`.
        The group is constructed during FLIR initialization and cached, so repeated calls return the same object. 
        
        OUTPUT: A :class:`FLIRDivisorGroup`.
        """
        self._ensure_flir_initialized()
        return FLIR.divisor_group(self)

    def Div(self):
        r""" Return the divisor group `\operatorname{Div}(A)`. 
        
        This is a shorthand for :meth:`divisor_group`. 
        
        OUTPUT: A :class:`FLIRDivisorGroup`, whose elements are instances of :class:`FLIRDivisor`.

        EXAMPLES::
            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.Div()
            Divisor group Div(A) of  A Banff Cluster Algebra with initial cluster variables x0, x1, x2 over Rational Field. with coefficients in Integer Ring

        """

        self._ensure_flir_initialized()
        return FLIR.Div(self)

    def divisor(self, data=None):
        r""" Construct a divisor (= an element of self.Div()) of this Banff cluster algebra. 
        
        INPUT: - ``data`` -- optional input describing a divisor. 
        The accepted forms are:
          - ``None`` or ``0``: the zero divisor; 
          - a :class:`FLIRDivisor` belonging to ``self.Div()``; 
          - a :class:`FLIRPrimeDivisor`, interpreted with coefficient one; 
          - a dictionary mapping prime divisors to integer coefficients; 
          - an iterable of ``(prime, coefficient)`` pairs;
          - a :class:`BanffClusterElement`, or an object coercible to the ambient fraction field, in which case its principal divisor is returned. 
          
        OUTPUT: A :class:`FLIRDivisor` in ``self.Div()``. 
          
        EXAMPLES: 
        
        Construct the zero divisor::
            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)

            sage: A.divisor() 
            0 
            sage: A.divisor(0) 
            0 
            
        Compute a principal divisor from an element:: 
            sage: x1, x2, x3 = A.gens()
            sage: D = A.divisor(x1) 
            sage: D
            1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2'), p=x0_1 + 1) +
             1*PrimeDivisor(chart=('x0_0p', 'x0_1', 'x0_2p'), p=x0_1 + 1)
            sage: D.parent() is A.Div() 
            True 
            
        Existing divisors in the same divisor group are returned unchanged:: 
            sage: A.divisor(D) is D 
            True 
        
        """
        self._ensure_flir_initialized()
        return FLIR.divisor(self, data)

    def class_data(self, recompute=False):
        r""" Return the data used to compute the divisor class group. 
        
        This includes the distinguished prime divisors, their valuation matrix, its Smith normal form, 
        and the resulting presentation of the divisor class group. 
        
        By default, previously computed data is reused. 
        Set ``recompute`` to ``True`` to discard the cached data and compute it again. 
        
        INPUT: - ``recompute`` -- boolean (default: ``False``); whether to recompute the class-group data 
        
        OUTPUT: A :class:`ClassGroupData` instance.
        """
        self._ensure_flir_initialized(recompute=recompute)
        return FLIR.class_data(self, recompute=recompute)

    def class_group(self):
        r""" Return the divisor class group of this Banff cluster algebra. 
        The class group is the group of Weil divisors modulo principal divisors.
        It is computed from the valuation matrix of the distinguished height-one primes. 
        
        OUTPUT: A Sage abelian group representing ``Cl(self)``
        
        EXAMPLES::
            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.class_group()
            Multiplicative Abelian group isomorphic to Z

            sage: B = Matrix([[0, 1], [-1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.class_group()
            Trivial Abelian group

            sage: B = Matrix([[0, 1, -1,1,0], [-1, 0, 1,0,-1], [1, -1, 0,-1,1],[-1,0,1,0,-1],[0,1,-1,1,0]])
            sage: A = BanffClusterAlgebra(B)
            sage: A.class_group()
            Multiplicative Abelian group isomorphic to Z x Z x Z x Z x Z x Z x Z x Z x Z x Z x Z
        """
        self._ensure_flir_initialized()
        return FLIR.class_group(self)

    def divisor_class(self, D):
        r"""Compute the class of a divisor in the divisor class group.

            Accept:
            - BanffClusterElement (the divisor of the element is used)
            - FLIRDivisor
            - FLIRPrimeDivisor (interpreted as 1*P)
            - dict {prime: exponent}

            Return: class [D] in Cl(A) = Z^r / im(M).

            EXAMPLES::

                sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
                sage: A = BanffClusterAlgebra(B)
                sage: x1, x2, x3 = A.gens()
                sage: z = A(x1*x3)
                sage: D = A.divisor(z)
                sage: A.divisor_class(D)
                (0)
                sage: P,Q,_,_ = A.extra_primes()
                sage: P = A.divisor(P)
                sage: Q = A.divisor(Q)
                sage: A.divisor_class(P)
                (1)
                sage: A.divisor_class((P+Q))
                (2)
        
            """
        self._ensure_flir_initialized()
        return FLIR.divisor_class(self, D)

    def is_principal_divisor(self, D):
        """
        Check if the divisor D is principal, i.e. if [D] = 0 in Cl(A).

        EXAMPLES:

        The divisor of an element of the FLIR is always principal::

            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x1, x2, x3 = A.gens()
            sage: z = A(x1*x3)
            sage: D = A.divisor(z)
            sage: A.is_principal_divisor(D)
            True
            sage: P,_,_,_ = A.extra_primes()
            sage: P = A.divisor(P)
            sage: A.is_principal_divisor(A.divisor(P))
            False
        """
        self._ensure_flir_initialized()
        return FLIR.is_principal_divisor(self, D)

    def principal_generator(self, D):
        """
        If a divisor D is principal, attempt to produce ``f \\in K(x)`` such that ``div(f) = D``,
        where K=self.fraction_field().

        Returns ``None`` if ``D`` is not principal.

        EXAMPLES:

        Recovering a generator for the divisor of a known element::

            sage: B = Matrix([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x0, x1, x2 = A.gens()
            sage: z = A(x0*x2)
            sage: D = A.divisor(z)
            sage: g = A.principal_generator(D)
            sage: g
            x0*x2
            sage: A.divisor(g) == D
            True
            sage: C = A.charts
            sage: c = C[0]
            sage: P = FLIRPrimeDivisor(c, x1 + 1)
            sage: Q = A.extra_primes()[1]
            sage: D1 = A.divisor_group()({P: 1, Q: 1})
            sage: D1.is_principal()
            True
            sage: D1.gen()
            (x1 + 1)/x0

        If the divisor is not principal, ``None`` is returned::

            sage: P = list(D.support())[0]
            sage: A.principal_generator(A.divisor({P: 1})) is None
            True

        
        """
        self._ensure_flir_initialized()
        return FLIR.principal_generator(self, D)

    def find_atoms(self, a):
        r"""
        Compute the atoms dividing a BanffClusterElement

        EXAMPLES::

            sage: B = Matrix([[0, 1, -1,1,0], [-1, 0, 1,0,-1], [1, -1, 0,-1,1],[-1,0,1,0,-1],[0,1,-1,1,0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x1, x2, x3, x4, x5 = A.gens()
            sage: f = x2*x4+x3
            sage: A.find_atoms(f)
            [x0, x4, (x1*x3 + x2)/x4, (x1*x3 + x2)/x0]
        """
        self._ensure_flir_initialized()
        return FLIR._find_atoms(self, a)

    def factorizations(self, a):
        r"""
        Compute all the factorizations in atoms of a BanffClusterElement

        EXAMPLES::

            sage: B = Matrix([[0, 1, -1,1,0], [-1, 0, 1,0,-1], [1, -1, 0,-1,1],[-1,0,1,0,-1],[0,1,-1,1,0]])
            sage: A = BanffClusterAlgebra(B)
            sage: x1, x2, x3, x4, x5 = A.gens()
            sage: f = x2*x4+x3
            sage: A.factorizations(f)
            ([x0, x4, (x1*x3 + x2)/x4, (x1*x3 + x2)/x0],
             [[(x4, 1), ((x1*x3 + x2)/x4, 1)], [(x0, 1), ((x1*x3 + x2)/x0, 1)]],
             [1, 1])
        """
        self._ensure_flir_initialized()
        return FLIR.factorizations(self, a)

    def __repr__(self):
        var_names = self.initial_cluster_variable_names()
        var_names_str = (" " if len(var_names) == 1 else "s ") + ", ".join(var_names)
        return (f" A Banff Cluster Algebra with initial cluster variable"
                f"{var_names_str}"
                f" over {self.base_ring()}.")
    


