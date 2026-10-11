import pytest
from random import randint

from itertools import chain
from sage.geometry.cone import Cone, random_cone
from sage.geometry.toric_lattice import ToricLattice
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.symbolic.ring import SR


def _gen_ops_test_cones():
    r"""
    Generate two cones that are suitable for testing (cross)
    positive, Lyapunov-like, and Z-operators.

    This extends :func:`sage.geometry.cone.random_cone` by choosing
    conservative parameters, and by rejecting any cone that looks too
    complicated for these methods to be computed in a reasonable
    amount of time. We generate two cones simultaneously because,
    while the cross-positive, Lyapunov-like, and Z-operators are
    defined on a single cone, the regular (non-cross) positive
    operators involve *two*.  We need to ensure that the various
    combinations are all computationally feasible; "success" is thus a
    property of the pair.
    """
    J2 = random_cone(max_ambient_dim=5, max_rays=10)
    J3 = random_cone(max_ambient_dim=5, max_rays=10)

    worst_dim = J2.dim() * max(J2.dual().dim(), J3.dual().dim())
    worst_npairs = J2.nrays() * max(J2.dual().nrays(), J2.dual().nrays())

    # How many pairs is too many pairs, in a given dimension? Note
    # that the maximum ambient dimension for both cones is 5, so there
    # are only a few cases to worry about.  "Too many" means that with
    # Cone(..., check=True), constructing the cone of operators itself
    # will be slow. The numbers below were obtained by hand.
    #
    # (product of dimensions => max number of ray pairs)
    limits = {
        9: 250,
        12: 100,
        15: 45,
        16: 36,
        20: 32,
        25: 30,
    }
    if not J2.is_proper() or not J3.is_proper():
        # When both J2 and J3 are proper, check=False gets used.
        if worst_dim in limits:
            # Ignore really low-dimensional cases.
            if worst_npairs >= limits[worst_dim]:
                return _gen_ops_test_cones()


    # We _also_ need to check that the dual we're going to take a dual
    # of is not too complicated: the complexity of dual() is bad in
    # general. It's also not easy to predict... so we are pretty much
    # guessing here.
    J2_J2_pi_dual = J2._positive_operators_dual(J2)
    J2_J3_pi_dual = J2._positive_operators_dual(J3)

    worst_nrays = max(J2_J2_pi_dual.nrays(), J2_J3_pi_dual.nrays())
    if worst_dim in limits and (worst_nrays >= limits[worst_dim]):
        # Again: do not expect this comparison to mean anything.
        # It's a bound; it's better than nothing.
        return _gen_ops_test_cones()

    # If we've made it this far, we monkey-patch cache a few more
    # intermediate results, starting with the posops-duals that
    # we already have handy. Note that there are two of them! The
    # tests had better not involve any others.
    def J2_pi_dual_method(s,o):
        if o is J2:
            return J2_J2_pi_dual
        return J2_J3_pi_dual  # "o" must be J3

    J2._positive_operators_dual = J2_pi_dual_method.__get__(J2)

    # Do the same for the dual of the cross-positive operators, but note
    # in this case that there is no "other" argument; the cross-positive
    # operators are on a single cone.
    J2_cp_dual = J2._cross_positive_operators_dual()
    J2_cp_dual_method = lambda s: J2_cp_dual
    J2._cross_positive_operators_dual = J2_cp_dual_method.__get__(J2)

    return (J2,J3)

@pytest.fixture
def K():
    r"""
    The baseline random cone fixture, meant to be used in tests
    that do not explode computationally.
    """
    return random_cone()

@pytest.fixture
def P():
    r"""
    Like :func:`K`, but guaranteed to be proper.
    """
    return random_cone(solid=True, strictly_convex=True)

@pytest.fixture
def J():
    r"""
    Another :func:`K`, for tests that require two cones of a
    similar nature.

    """
    return random_cone()

@pytest.fixture
def Q():
    r"""
    Like :func:`J`, but guaranteed to be proper.
    """
    return random_cone(solid=True, strictly_convex=True)

@pytest.fixture
def K2_and_K3():
    return _gen_ops_test_cones()

@pytest.fixture
def K2(K2_and_K3):
    r"""
    Cone (1 of 2) used to test positive, cross-positive,
    Lyapunov-like, and Z-operators.
    """
    return K2_and_K3[0]

@pytest.fixture
def K3(K2_and_K3):
    r"""
    Cone (2 of 2) used to test positive, cross-positive, Lyapunov-like,
    and Z-operators. Used mainly for positive operators tests, where
    two different cones can be involved.
    """
    return K2_and_K3[1]

@pytest.fixture
def K2_cp_gens(K2):
    return K2.cross_positive_operators_gens()

@pytest.fixture
def K2_cp_cone(K2, K2_cp_gens):
    L = ToricLattice(K2.lattice_dim()**2)
    return Cone((g.list() for g in K2_cp_gens),
                lattice=L,
                check=False)

@pytest.fixture
def K2_posops_gens(K2):
    return K2.positive_operators_gens()

@pytest.fixture
def K2_K3_posops_gens(K2, K3):
    return K2.positive_operators_gens(K3)

@pytest.fixture
def K2_posops_cone(K2, K2_posops_gens):
    L = ToricLattice(K2.lattice_dim()**2)
    return Cone((g.list() for g in K2_posops_gens),
                lattice=L,
                check=False)

@pytest.fixture
def K2_K3_posops_cone(K2, K3, K2_K3_posops_gens):
    L = ToricLattice(K2.lattice_dim()*K3.lattice_dim())
    return Cone((g.list() for g in K2_K3_posops_gens),
                lattice=L,
                check=False)

@pytest.fixture
def K2_ll_basis(K2):
    return K2.lyapunov_like_basis()

def test_dual_lattice_dual_is_original(K):
    r"""
    The dual lattice of a dual lattice is the original.
    """
    assert K.dual_lattice().dual() is K.lattice()

def test_codimension_bounds(K):
    r"""
    Codimension is an integer between zero and the ambient
    (lattice or vector space) dimension.
    """
    c = K.codim()
    assert c.is_integral() and 0 <= c <= K.lattice_dim()

def test_codim_is_zero_iff_cone_is_solid(K):
    r"""
    A cone is solid if and only if it has codimension zero.
    """
    assert K.is_solid() == K.codim().is_zero()

def test_codim_is_dual_lineality(K):
    r"""
    The codimension of a cone is equal to the lineality of its
    dual.
    """
    assert K.codim() == K.dual().lineality()

def test_cone_is_solid_iff_it_spans_its_lattice(K):
    r"""
    The span of a cone is the entire ambient space if and only if
    that cone is solid.
    """
    V1 = K.lattice().vector_space()
    V2 = K.span().vector_space()
    assert K.is_solid() == (V1 == V2)

def test_dual_of_dual_is_original(K):
    r"""
    The dual cone of a dual cone is the original cone.
    """
    assert K.dual().dual() is K

def test_cones_are_equivalent_to_themselves(K):
    r"""
    Sanity check: a cone is always equivalent to itself.
    """
    assert K.is_equivalent(K)

def test_cones_are_faces_of_themselves(K):
    r"""
    Every cone is (trivially) a face of itself.
    """
    assert K.is_face_of(K)

def pointed_cones_are_isomorphic_to_themselves(K):
    r"""
    Sanity check: a cone is always isomorphic to itself.

    Isomorphism is implemented only for pointed cones, so
    we use the strict quotient of the arbitrary (random)
    cone ``K``.
    """
    K_p = K.strict_quotient()
    assert K_p.is_isomorphic(K_p)

def test_lineality_bounds(K):
    r"""
    Lineality is an integer between zero and the ambient
    dimension.
    """
    l = K.lineality()
    assert l.is_integral() and 0 <= l <= K.lattice_dim()

def test_lineality_zero_iff_pointed(K):
    r"""
    The lineality of a pointed cone is zero, essentially by
    definition.
    """
    assert K.is_pointed() == K.lineality().is_zero()

def test_lineality_space_is_dual_perp(K):
    r"""
    The linear subspace of any closed convex cone can be identified
    with the orthogonal complement of the span of its dual.
    """
    expected = K.dual().span().vector_space().complement()
    assert K.linear_subspace() == expected

def test_strict_quotient_is_pointed(K):
    r"""
    The strict quotient of any cone should be pointed (that's the
    point of the method).
    """
    assert K.strict_quotient().is_strictly_convex()

def test_strict_quotient_of_pointed_is_itself(K):
    r"""
    The strict quotient of a cone is itself if and only if that
    cone is pointed (strictly convex).
    """
    assert K.is_strictly_convex() == (K.strict_quotient() is K)

def test_strict_quotient_dim_bound(K):
    r"""
    The orthogonal complement of a cone's linear subspace has the
    same dimension as its dual, so the strict quotient cannot have a
    larger dimension than the dual.
    """
    assert K.strict_quotient().dim() <= K.dual().dim()

def test_strict_quotient_is_idempotent(K):
    r"""
    As corollary of the strict quotient being pointed and the
    strict quotient of a pointed cone being itself, the strict
    quotient is idempotent.
    """
    C = K.strict_quotient()
    assert C.strict_quotient() == C

def test_solid_restriction_is_solid(K):
    r"""
    The solid restriction of a cone is itself if and only if that
    cone is solid.
    """
    assert K.is_solid() == (K.solid_restriction() is K)

def test_solid_restriction_has_same_properties_as_original(K):
    r"""
    The solid restriction of a cone has the same dimension,
    number of rays, etc. as the the original cone. (The restriction
    changes only its representation and ambient space.)
    """
    sr = K.solid_restriction()
    assert all((sr.dim() == K.dim(),
                sr.n_rays() == K.n_rays(),
                sr.lineality() == K.lineality(),
                len(sr.facets()) == len(K.facets())))

def test_solid_iff_dual_pointed(K):
    r"""
    A closed convex cone is solid if and only if its dual is
    strictly convex.
    """
    assert K.is_solid() == K.dual().is_strictly_convex()

def test_random_element_membership(K):
    r"""
    Random elements of a cone belong to the cone, its ambient
    vector space, and its lattice. The same is true for conic
    combinations of random elements.
    """
    L = K.lattice()
    V = L.vector_space()
    F = L.base_field()

    # Note that after the change of ring (from ZZ to QQ), we lose
    # membership in the lattice.
    e1 = sum((K.random_element()
              for _ in range(randint(0,10))),
             start=L.zero())
    e2 = sum((K.random_element(ring=F)
              for _ in range(randint(0,10))),
             start=L.zero())

    assert all(x in K and x in V for x in (e1,e2))
    assert e1 in L

def test_pointed_cones_contain_no_random_line(K):
    r"""
    A pointed cone contains no lines, and thus no negative
    multiples of any of its elements (besides the origin).

    The input cone ``K`` is arbitrary (random), so we take
    its strict quotient to obtain a pointed cone.
    """
    K_p = K.strict_quotient()
    x = K_p.random_element()
    assert x.is_zero() or not K_p.contains(-x)

@pytest.mark.long
def test_reducibility_preserved_by_isomorphism(K):
    r"""
    Reducibility is preserved under linear isomorphisms. We take
    the strict quotient of ``K`` to ensure that we have a pointed cone
    for which reducibility is well-defined.
    """
    C = K.strict_quotient()

    n = C.ambient_dim()
    F = C.lattice().base_field()
    q = F._random_nonzero_element()
    A = q*matrix.random(F, n, algorithm='unimodular')
    AC = Cone([ r*A for r in C.rays() ], lattice=C.lattice())
    assert C.is_reducible() == AC.is_reducible()

@pytest.mark.long
def test_reducibility_criteria(K,J,P,Q,K2,K3):
    r"""
    In [GT2014]_ it is shown that a (nontrivial) proper polyhedral
    cone is irreducible if and only if its Lyapunov rank is one.
    A related test combines Theorem 4.7 of [HFP1976]_ with the
    Z-operator algorithm in [Or2018b]_.

    Loop through several fixtures to increase the chances that one of
    them is both proper AND nontrivial. If not, just make a new one.
    """
    winners = ( C for C in (P,Q,K,J,K2,K3)
                if not C.is_trivial() and C.is_proper() )
    try:
        C = next(winners)
    except StopIteration:
        C = random_cone(strictly_convex=True,
                        solid=True,
                        max_ambient_dim=5,
                        min_rays=1,
                        max_rays=7)

    assert C.is_reducible() == (C.lyapunov_rank() != 1)
    d = C._cross_positive_operators_dual().dim()
    assert C.is_reducible() == (d < C.dim()**2 - 1)

@pytest.mark.longlong
def test_is_cross_positive(K2, K2_cp_gens):
    r"""
    The cross-positive property is possessed by every
    cross-positive operator.
    """
    F = K2.lattice().base_field()
    n = K2.lattice_dim()

    # Check over SR, too.
    assert all(A.is_cross_positive_on(K2) for A in K2_cp_gens)
    assert all(A.change_ring(SR).is_cross_positive_on(K2) for A in K2_cp_gens)

    # The identity matrix is always cross-positive.
    assert matrix.identity(F, n).is_cross_positive_on(K2)

    # The zero matrix is always cross-positive.
    assert matrix.zero(F, n).is_cross_positive_on(K2)

@pytest.mark.longlong
def test_is_Z_operator(K2, K2_cp_gens):
    r"""
    The Z property is possessed by every Z-operator.
    """
    F = K2.lattice().base_field()
    n = K2.lattice_dim()

    # Obviously cheating, but the Z_operators_gens() method is simple
    # enough to justify faking it here; otherwise we'd have to
    # recompute the cross-positive gens that we already have cached.
    Z_gens = [ -cp for cp in K2_cp_gens ]

    # Check over SR, too.
    assert all(A.is_Z_operator_on(K2) for A in Z_gens)
    assert all(A.change_ring(SR).is_Z_operator_on(K2) for A in Z_gens)

    # The identity matrix is always a Z-operator.
    assert matrix.identity(F, n).is_Z_operator_on(K2)

    # The zero matrix is always a Z-operator.
    assert matrix.zero(F, n).is_Z_operator_on(K2)

@pytest.mark.longlong
def test_is_positive(K2, K3, K2_K3_posops_gens):
    r"""
    Every positive operator is positive.
    """
    F = K2.lattice().base_field()
    m = K2.lattice_dim()
    n = K3.lattice_dim()

    # Check over SR, too.
    assert all(A.is_positive_operator_on(K2,K3) for A in K2_K3_posops_gens)
    assert all(A.change_ring(SR).is_positive_operator_on(K2,K3)
               for A in K2_K3_posops_gens)

    # The identity matrix is a positive-operator on every cone.
    assert matrix.identity(F, m).is_positive_operator_on(K2)

    # The zero matrix is always a positive operator.
    assert matrix.zero(F, n, m).is_positive_operator_on(K2, K3)

@pytest.mark.longlong
def test_is_lyapunov_like(K2, K2_ll_basis):
    r"""
    Every Lyapunov-like operator is Lyapunov-like.
    """
    F = K2.lattice().base_field()
    n = K2.lattice_dim()

    # Check over SR, too.
    assert all(A.is_lyapunov_like_on(K2) for A in K2_ll_basis)
    assert all(A.change_ring(SR).is_lyapunov_like_on(K2)
               for A in K2_ll_basis)

    # The identity matrix is always Lyapunov-like.
    assert matrix.identity(F, n).is_lyapunov_like_on(K2)

    # The zero matrix is always Lyapunov-like.
    assert matrix.zero(F, n).is_lyapunov_like_on(K2)

@pytest.mark.longlong
def test_lyapunov_like_is_plus_minus_Z(K2):
    r"""
    A matrix is Lyapunov-like on a cone if and only if both the
    matrix and its negation are cross-positive (or Z) on the cone.
    """
    R = K2.lattice().base_field()
    A = matrix.random(R, K2.lattice_dim())
    actual = A.is_lyapunov_like_on(K2)
    expected = all( B.is_cross_positive_on(K2) for  B in (A, -A) )
    assert actual == expected

@pytest.mark.longlong
def test_cross_positive_linspace(K2, K2_cp_cone, K2_ll_basis):
    r"""
    The lineality space of the cone of cross-positive operators is
    the space of Lyapunov-like operators [Or2018b]_.
    """
    L = ToricLattice(K2.lattice_dim()**2)
    V = L.vector_space()
    long_basis = (V(l.list()) for l in K2_ll_basis)
    lls = V.span(long_basis)
    assert K2_cp_cone.linear_subspace() == lls

@pytest.mark.longlong
def test_cross_positive_permutation(K2, K2_cp_gens):
    r"""
    The cross-positive operators of a permuted cone can be
    obtained by conjugation.
    """
    L = ToricLattice(K2.lattice_dim()**2)
    p = SymmetricGroup(K2.lattice_dim()).random_element().matrix()
    pK2 = Cone((p*k for k in K2), K2.lattice(), check=False)
    actual = Cone((g.list() for g in pK2.cross_positive_operators_gens()),
                  lattice=L,
                  check=False)
    expected = Cone(((p*g*p.inverse()).list() for g in K2_cp_gens),
                    lattice=L,
                    check=False)
    assert actual.is_equivalent(expected)

@pytest.mark.longlong
def test_cross_positive_adjoint(K2, K2_cp_cone):
    r"""
    An operator is cross-positive on a cone if and only if its
    adjoint is cross-positive on the dual of that cone [Or2018b]_.
    """
    n = K2.lattice_dim()
    L = ToricLattice(n**2)
    F = K2.lattice().base_field()
    M = MatrixSpace(F, n)
    W = L.vector_space()

    dual_gens = K2.dual().cross_positive_operators_gens()
    cp_star = Cone((g.list() for g in dual_gens),
                   lattice=L,
                   check=False)

    A = M(K2_cp_cone.random_element(ring=F).list())
    assert cp_star.contains(W(A.transpose().list()))

    A = M(cp_star.random_element(ring=F).list())
    assert K2_cp_cone.contains(W(A.transpose().list()))

@pytest.mark.longlong
def test_cross_positive_dim(K2, K2_cp_cone, K2_posops_cone):
    r"""
    The lineality spaces of the duals of the positive and cross-
    positive operator cones are equal. From this it follows that
    the dimensions of the cross-positive operator cone and positive
    operator cone are equal [Or2018b]_.
    """
    pi_star = K2._positive_operators_dual(K2)
    cp_star = K2._cross_positive_operators_dual()

    assert pi_star.linear_subspace() == cp_star.linear_subspace()
    assert K2_posops_cone.dim() == K2_cp_cone.dim()

@pytest.mark.longlong
def test_posops_dual_linspace(K2):
    r"""
    The lineality space of the dual of the positive operators
    can be computed from the lineality spaces of the cone and
    its dual [Or2018b]_.
    """
    pi_dual = K2._positive_operators_dual(K2)
    V = pi_dual.ambient_vector_space()
    actual = pi_dual.linear_subspace()
    U1 = ( V((s.tensor_product(x)).list())
           for x in K2.lines()
           for s in K2.dual() )
    U2 = ( V((s.tensor_product(x)).list())
           for x in K2
           for s in K2.dual().lines() )
    expected = V.span(chain(U1,U2))
    assert actual == expected

@pytest.mark.longlong
def test_posops_dual_lineality(K2):
    r"""
    The lineality of the dual of the positive operators is known
    from the lineality space of the original [Or2018b]_.
    """
    n = K2.lattice_dim()
    m = K2.dim()
    l = K2.lineality()
    actual = K2._positive_operators_dual(K2).lineality()
    expected = l*(m - l) + m*(n - m)
    assert actual == expected

@pytest.mark.longlong
def test_posops_dimension(K2, K2_posops_cone):
    r"""
    The dimension of the positive operators on a cone depends on the
    dimension and lineality of that cone [Or2018b]_.
    """
    n = K2.lattice_dim()
    m = K2.dim()
    l = K2.lineality()
    actual = K2_posops_cone.dim()
    expected = n**2 - l*(m - l) - (n - m)*m
    assert actual == expected

@pytest.mark.longlong
def test_posops_lineality_from_gens(K2, K2_posops_cone):
    r"""
    The lineality of the positive operators follows from the
    description of its generators [Or2018b]_.
    """
    n = K2.lattice_dim()
    actual = K2_posops_cone.lineality()
    expected = n**2 - K2.dim()*K2.dual().dim()
    assert actual == expected

@pytest.mark.longlong
def test_posops_proper(K2, K2_posops_cone):
    r"""
    A cone is proper if and only if its positive operators form a
    proper cone [Or2018b]_.
    """
    assert K2.is_proper() == K2_posops_cone.is_proper()

@pytest.mark.longlong
def test_posops_on_element(K2, K3, K2_K3_posops_cone):
    r"""
    A random positive operator should send a random element of one
    cone into the other cone.
    """
    F = K2.lattice().base_field()
    P = matrix(K3.lattice_dim(),
               K2.lattice_dim(),
               K2_K3_posops_cone.random_element(F).list())
    assert K3.contains(P*K2.random_element(ring=F))

@pytest.mark.longlong
def test_posops_permutation(K2, K2_posops_gens):
    r"""
    The positive operators on a permuted cone can be obtained by
    conjugation.
    """
    L = ToricLattice(K2.lattice_dim()**2)
    p = SymmetricGroup(K2.lattice_dim()).random_element().matrix()
    pK = Cone((p*k for k in K2), K2.lattice(), check=False)
    actual = Cone((g.list() for g in pK.positive_operators_gens()),
                  lattice=L,
                  check=False)
    expected = Cone(((p*g*p.inverse()).list() for g in K2_posops_gens),
                    lattice=L,
                    check=False)
    assert actual.is_equivalent(expected)

@pytest.mark.longlong
def test_posops_adjoint(K2, K3, K2_K3_posops_cone):
    r"""
    An operator is positive from one cone to another if and only if
    its adjoint is positive from the dual of the second cone to the
    dual of the first.
    """
    F = K2.lattice().base_field()
    n = K2.lattice_dim()
    m = K3.lattice_dim()
    L = ToricLattice(n*m)
    W = L.vector_space()
    M_fwd = MatrixSpace(F, m, n)
    M_back = MatrixSpace(F, n, m)

    pi_fwd = K2_K3_posops_cone
    pi_back_gens = K3.dual().positive_operators_gens(K2.dual())
    pi_back = Cone((g.list() for g in pi_back_gens),
                   lattice=L,
                   check=False)

    A = M_fwd(pi_fwd.random_element(ring=F).list())
    assert pi_back.contains(W(A.transpose().list()))
    A = M_back(pi_back.random_element(ring=F).list())
    assert pi_fwd.contains(W(A.transpose().list()))

@pytest.mark.longlong
def test_posops_lyapunov_rank(K2, K3, K2_K3_posops_cone):
    r"""
    The Lyapunov rank of the positive operators is the product of
    the Lyapunov ranks of the associated cones if both are proper
    [Or2018a]_.
    """
    if not K2_K3_posops_cone.is_proper():
        return

    beta1 = K2.lyapunov_rank()
    beta2 = K3.lyapunov_rank()
    assert K2_K3_posops_cone.lyapunov_rank() == beta1*beta2

@pytest.mark.longlong
def test_posops_LL(K2, K3, K2_K3_posops_cone):
    r"""
    Lyapunov-like operators on a proper polyhedral positive operator
    cone can be computed from the Lyapunov-like operators on the cones
    with respect to which the operators are positive [Or2018a]_.
    """
    if not K2_K3_posops_cone.is_proper():
        return

    F = K2.lattice().base_field()
    m = K2.lattice_dim()
    n = K3.lattice_dim()
    M1 = MatrixSpace(F, m)
    M2 = MatrixSpace(F, n)
    W = VectorSpace(F, (m**2)*(n**2))

    tps = (M2(s.list()).tensor_product(M1(x.list()))
           for x in K2.dual().lyapunov_like_basis()
           for s in K3.lyapunov_like_basis())
    expected = W.span( W(x.list()) for x in tps )
    LL_pi = K2_K3_posops_cone.lyapunov_like_basis()
    actual = W.span( W(x.list()) for x in LL_pi )
    assert actual == expected

@pytest.mark.longlong
def test_lyapunov_like_adjoint(K2, K2_ll_basis):
    r"""
    The Lyapunov-like transformations on a cone and its dual are
    transposes of one another. However, there's no reason to expect
    that one basis will consist of transposes of the other.
    """
    F = K2.lattice().base_field()
    V = VectorSpace(F, K2.lattice_dim()**2)
    LL2 = (A.transpose() for A in K2.dual().lyapunov_like_basis())
    LL1_vecs = (V(m.list()) for m in K2_ll_basis)
    LL2_vecs = (V(m.list()) for m in LL2)
    assert V.span(LL1_vecs) == V.span(LL2_vecs)

@pytest.mark.longlong
def test_LL_closed_under_lie_bracket(K2, K2_ll_basis):
    r"""
    The space of all Lyapunov-like transformations is a Lie algebra
    and should therefore be closed under the lie bracket.
    """
    W = VectorSpace(K2.lattice().base_field(), K2.lattice_dim()**2)
    LL_W = W.span( W(m.list()) for m in K2_ll_basis )
    brackets = (W((A1*A2 - A2*A1).list()) for A1 in K2_ll_basis
                                          for A2 in K2_ll_basis)
    assert all(b in LL_W for b in brackets)

@pytest.mark.long
def test_lyapunov_rank_of_direct_sum(P, Q):
    r"""
    Lyapunov rank should be additive on a product of proper cones
    [RNPA2011]_.
    """
    PQ = P.cartesian_product(Q)
    assert PQ.lyapunov_rank() == P.lyapunov_rank() + Q.lyapunov_rank()

@pytest.mark.long
def test_lyapunov_rank_automorphism_invariance(K):
    r"""
    Lyapunov rank should be invariant under a linear isomorphism
    [Or2017]_.
    """
    F = K.lattice().base_field()
    n = K.lattice_dim()
    A = matrix.random(F, n, algorithm='unimodular')
    AK = Cone((A*r for r in K), lattice=K.lattice())
    assert K.lyapunov_rank() == AK.lyapunov_rank()

@pytest.mark.long
def test_lyapunov_rank_dual_invariance(K):
    r"""
    Lyapunov rank is invariant under the duality operator
    [RNPA2011]_.
    """
    assert K.lyapunov_rank() == K.dual().lyapunov_rank()

@pytest.mark.long
def test_lyapunov_rank_valid_range(K, P):
    r"""
    The Lyapunov rank of a proper polyhedral cone in a
    non-trivial `n`-dimensional space can be any number between `1`
    and `n` inclusive, excluding `n-1` [GT2014]_. No polyhedral closed
    convex cone in `n` dimensions (proper or otherwise) has Lyapunov
    rank `n-1` [Or2017]_.
    """
    b = P.lyapunov_rank()
    n = P.lattice_dim()
    if not P.is_trivial():
        # The paper overlooks the trivial case!
        assert 1 <= b <= n
    assert b != n-1

    assert K.lyapunov_rank() != K.lattice_dim()-1

@pytest.mark.long
def test_lyapunov_rank_from_proper_subcone(K):
    r"""
    The calculation of the Lyapunov rank of an improper cone can
    be reduced to that of a proper cone [Or2017]_.
    """
    K_SP = K.solid_restriction().strict_quotient()
    l = K.lineality()
    c = K.codim()
    actual = K.lyapunov_rank()
    expected = K_SP.lyapunov_rank() + K.dim()*(l + c) + c**2
    assert actual == expected

@pytest.mark.long
def test_lyapunov_rank_agrees_with_naive_definition(K):
    r"""
    The Lyapunov rank of a cone is the length of a
    :meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.lyapunov_like_basis`
    for it.
    """
    assert K.lyapunov_rank() == len(K.lyapunov_like_basis())

@pytest.mark.long
def test_lyapunov_rank_of_perfect_cone(K):
    r"""
    A "perfect" cone has Lyapunov rank `n` or more in `n`
    dimensions. We can make any cone perfect by adding a slack
    variable.
    """
    L = ToricLattice(K.lattice_dim() + 1)
    K_perf = Cone([r.list() + [0] for r in K], lattice=L)
    assert K_perf.lyapunov_rank() >= K_perf.lattice_dim()

@pytest.mark.long
def test_discrete_complementarity_set_of_dual(K):
    r"""
    A discrete complementarity set for the dual can be obtained by
    switching components in a discrete complementarity set of the
    original cone.
    """
    dcs_dual = K.dual().discrete_complementarity_set()
    expected = tuple((x,s) for (s,x) in dcs_dual)
    actual = K.discrete_complementarity_set()
    assert sorted(actual) == sorted(expected)

@pytest.mark.long
def test_discrete_complementarity_set_is_complementary(K):
    r"""
    The pairs in a discrete complementarity set are in fact
    complementary.
    """
    dcs = K.discrete_complementarity_set()
    assert sum((s*x).abs() for (x,s) in dcs) == 0


def test_random_cone_params_are_respected(K):
    r"""
    The dim/ray bounds and pointed/solid constraints should
    actually be met by the result.
    """
    solid = bool(randint(0,1))
    pointed = bool(randint(0,1))
    min_ambient_dim = randint(5,10)
    min_rays = randint(0,8)

    # Start with random maxima, but fudge them a bit so that they
    # don't contradict the minima.
    max_rays = randint(min_rays,2*min_rays)
    max_rays = max(max_rays, min_ambient_dim+3)
    max_ambient_dim = randint(min_ambient_dim, 2*min_ambient_dim)
    max_ambient_dim = max(max_ambient_dim, min_ambient_dim+4)

    C = random_cone(min_rays=min_rays,
                    max_rays=max_rays,
                    solid=solid,
                    strictly_convex=pointed,
                    min_ambient_dim=min_ambient_dim,
                    max_ambient_dim=max_ambient_dim)

    assert all((
        min_rays <= C.n_rays() <= max_rays,
        min_ambient_dim <= C.lattice_dim() <= max_ambient_dim,
        solid == C.is_solid(),
        pointed == C.is_strictly_convex()
    ))

    # Again with a fixed lattice (in which we know a satisfactory cone
    # exists, because we just found one) rather than dimension bounds.
    L = C.lattice()
    C = random_cone(min_rays=min_rays,
                    max_rays=max_rays,
                    solid=solid,
                    strictly_convex=pointed,
                    lattice=L)

    assert all((
        min_rays <= C.n_rays() <= max_rays,
        C.lattice() is L,
        solid == C.is_solid(),
        pointed == C.is_strictly_convex()
    ))
