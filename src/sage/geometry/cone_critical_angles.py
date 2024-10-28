r"""
Find maximal angles between polyhedral convex cones

.. WARNING::

    This module is considered internal and its contents are subject to
    change at any time without (deprecation) warning. The stable
    interface is
    :meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.max_angle`.

Finding the maximal (or equivalently, the minimal) angle between two
polyhedral convex cones is a hard nonconvex optimization problem. The
problem for a single cone was introduced in [IS2005]_, and was later
extended in [SS2016]_ to two cones as a generalization of the
principal angle between two vector subspaces.

Seeger and Sossa proposed an algorithm in [SS2016]_ to find maximal
angles, and [Or2020]_ elaborates on that algorithm. It is this latest
improvement that is implemented (more or less) by this module. The
fact that perturbations of pointed cones may not change the answer too
much [Or2024]_ is taken into consideration to avoid pathological
cases.

This module is internal to SageMath; the interface presented to users
consists of a public method,
:meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.max_angle` for
polyhedral convex cones. Even though all of the functions in this
module are internal, some are more internal than others. There are a
few functions that are used only in doctests, and not by any code that
an end-user would run. Breaking somewhat with tradition, only those
methods have been prefixed with an underscore.
"""

from sage.functions.trig import arccos, cos
from sage.matrix.constructor import matrix
from sage.misc.misc import powerset
from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import AA
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.symbolic.constants import pi


def _normalize_gevp_solution(gevp_solution):
    r"""
    Normalize the results of :func:`solve_gevp_nonzero` and
    :func:`_solve_gevp_naive`.

    Those two functions return solutions (pairs of vectors) to an
    eigenvalue problem, but eigenvectors are only unique up to a
    scalar multiple. This function normalizes those results so that
    every eigenvector has a leading entry of positive one. This allows
    us to identity equivalent solutions to those problems.

    INPUT:

    A quartet ``gevp_solution`` whose components are, in order:

    - ``eigenvalue`` -- ignored

    - ``xi`` -- first component of the `( \xi, \eta )` eigenvector

    - ``eta`` -- second component of the `( \xi, \eta )` eigenvector

    - ``multiplicity`` -- ignored

    OUTPUT:

    If `c` is the first nonzero component of the concatenated `( \xi,
    \eta )` vector, then a quartet whose components are, in order:

    - ``eigenvalue`` -- the unmodified ``eigenvalue`` argument

    - ``xi*(1/c)`` -- the `\xi` component normalized so that the first
      nonzero component of the concatenated vector
      `( \xi, \eta )` is positive one

    - ``eta*(1/c)`` -- the `\eta` component normalized so that the
      first nonzero component of the concatenated vector
      `( \xi, \eta )` is positive one

    - ``multiplicity`` -- the unmodified ``multiplicity`` argument

    If there is no such `c` (that is, if both `\xi` and `\eta` are
    zero), then the entire input quartet is returned unmodified.

    EXAMPLES::

        sage: from sage.geometry.cone_critical_angles import (
        ....:   _normalize_gevp_solution)
        sage: s1 = (-1, vector(QQ,[0,-2]), vector(QQ,[1]), 1)
        sage: _normalize_gevp_solution(s1)
        (-1, (0, 1), (-1/2), 1)
        sage: s2 = (1, vector(QQ,[0,0]), vector(QQ,[0,0,-1]), 2)
        sage: _normalize_gevp_solution(s2)
        (1, (0, 0), (0, 0, 1), 2)
    """
    eigenvalue, xi, eta, multiplicity = gevp_solution
    from itertools import chain

    # We'll use this default of zero as the scaling factor if we don't
    # find a better one; that is, if xi = eta = 0 -- in which case the
    # additional multiplication by zero is a no-op.
    scale = 0
    for c in chain(xi, eta):
        if c != 0:
            scale = ~c
            break

    xi *= scale
    xi.set_immutable()
    eta *= scale
    eta.set_immutable()
    return (eigenvalue, xi, eta, multiplicity)


def _random_admissible_cone(ambient_dim):
    r"""
    Generate a random cone in a lattice of dimension
    ``ambient_dim`` that isn't trivial.

    This is a convenience method used to simplify some test cases. The
    number of rays that the cone possesses is limited to two more than
    ``ambient_dim``; so, for example, you will not get more than five
    rays in a three-dimensional space. This limits the amount of time
    spent in any one test case. In contrast with the definition in
    [Or2020]_ we consider the full ambient space to be admissible (it
    doesn't hurt anything, and [Or2024]_ was forced to allow it).

    INPUT:

    - ``ambient_dim`` -- positive integer representing the dimension
      of the ambient lattice in which the returned cone lives

    OUTPUT:

    A "random" nontrivial closed convex cone in a lattice of dimension
    ``ambient_dim``.

    A :exc:`ValueError` is raised if ``ambient_dim`` is not
    positive.

    EXAMPLES:

    The result has all of the desired properties::

        sage: from sage.geometry.cone_critical_angles import (
        ....:   _random_admissible_cone )
        sage: K = _random_admissible_cone(5)
        sage: K.lattice_dim()
        5
        sage: K.is_trivial()
        False

    Unless the ``ambient_dim`` argument is nonsense::

        sage: from sage.geometry.cone_critical_angles import (
        ....:   _random_admissible_cone )
        sage: K = _random_admissible_cone(0)
        Traceback (most recent call last):
        ...
        ValueError: there are no nontrivial cones in dimension 0
    """
    if ambient_dim < 1 or ambient_dim not in ZZ:
        # The random_cone() method already crashes if we ask the
        # impossible of it, but having this here emits a more sensible
        # error message.
        raise ValueError("there are no nontrivial cones in dimension %d"
                         % ambient_dim)

    args = { 'min_ambient_dim': ambient_dim,
             'max_ambient_dim': ambient_dim,
             'min_rays': 1,
             'max_rays': ambient_dim+2 }

    from sage.geometry.cone import random_cone
    return random_cone(**args)

    return K


def gevp_licis(G):
    r"""
    Return all nonempty subsets of indices for the columns of
    ``G`` that correspond to linearly independent sets (of columns of
    ``G``).

    Mnemonic: linearly independent column-index subsets (LICIS).

    The returned lists are all sorted in the same (the natural) order;
    and are returned as lists so that they may be used to index into
    the rows/columns of matrices.

    INPUT:

    - ``G`` -- the matrix whose linearly independent column index sets
      we want

    OUTPUT:

    A generator that returns sorted lists of natural numbers. Each
    generated list ``I`` is a set of indices corresponding to columns
    of ``G`` that, when considered as a set, is linearly independent.

    EXAMPLES:

    The linearly independent subsets of the matrix corresponding to a
    line (with two generators pointing in opposite directions) are the
    one-element subsets, since the only two-element subset isn't
    linearly independent::

        sage: from sage.geometry.cone_critical_angles import gevp_licis
        sage: K = Cone([(1,0),(-1,0)])
        sage: G = matrix.column(K.rays())
        sage: list(gevp_licis(G))
        [[0], [1]]

    The matrix for the trivial cone has no linearly independent
    subsets, since we require them to be nonempty::

        sage: from sage.geometry.cone_critical_angles import gevp_licis
        sage: trivial_cone = cones.trivial(0)
        sage: trivial_cone.is_trivial()
        True
        sage: list(gevp_licis(matrix.column(trivial_cone.rays())))
        []

    All rays in the nonnegative orthant of `R^{n}` are
    linearly independent, so we should get back `2^{n} - 1` subsets
    after accounting for the absence of the empty set::

        sage: from sage.geometry.cone_critical_angles import gevp_licis
        sage: K = cones.nonnegative_orthant(3)
        sage: G = matrix.column(K.rays())
        sage: len(list(gevp_licis(G))) == 2^(K.nrays()) - 1
        True

    TESTS:

    All sets corresponding to the returned indices should be linearly
    independent::

        sage: from sage.geometry.cone_critical_angles import gevp_licis
        sage: K = random_cone(max_rays=8)
        sage: G = matrix.column(K.rays())
        sage: all( len(s) == K.rays(s).dimension() for s in gevp_licis(G) )
        True
    """
    from sage.matroids.linear_matroid import LinearMatroid

    # There's a fast implementation of this for matroids, but we need
    # to drop the empty set from its output and convert the rest to
    # lists that are all sorted in the same order.
    return map(sorted, filter(bool, LinearMatroid(G).independent_sets()))


def _solve_gevp_naive(GG, HH, M, I, J):
    r"""
    Solve the generalized eigenvalue problem in Theorem 3
    [Or2020]_ in a very naive way, by (slowly) inverting the matrices
    and finding the eigenvalues in the product space.

    This is used only for testing, to ensure that the smart way of
    solving the generalized eigenvalue problem via
    :func:`solve_gevp_zero` and :func:`solve_gevp_nonzero` returns the
    same answers as the dumb way.

    This returns a generator, like :func:`solve_gevp_zero` and
    :func:`solve_gevp_nonzero`.

    INPUT:

    See the arguments for :func:`solve_gevp_nonzero`.

    ALGORITHM:

    We construct the two matrices `A` and `B` in Theorem 3 [Or2020]_
    in block form, and then use the naive "inverse" method on `B` to
    move it to the left and obtain `M = B^{-1}A`. We then compute the
    right eigenvectors of the whole big matrix `M`.

    EXAMPLES:

    A simple usage example, that also appears as Example 3 in
    [Or2020]_::

        sage: from sage.geometry.cone_critical_angles import _solve_gevp_naive
        sage: K = cones.nonnegative_orthant(2)
        sage: G = matrix.column(K.rays())
        sage: GG = G.transpose() * G
        sage: I = [0]
        sage: J = [1]
        sage: list(_solve_gevp_naive(GG,GG,GG,I,J))
        [(0, (1), (0), 2), (0, (0), (1), 2)]

    Check Example 4 [Or2020]_ symbolically to ensure that we get
    eigenspaces of dimension `n=2` corresponding to the eigenvalues
    `\cos\theta = -1` and `\cos\theta = 1`::

        sage: from sage.geometry.cone_critical_angles import _solve_gevp_naive
        sage: g11,g12,g21,g22 = SR.var('g11,g12,g21,g22', domain='real')
        sage: h11,h12,h21,h22 = SR.var('h11,h12,h21,h22', domain='real')
        sage: gs = [[g11,g12], [g21,g22]]
        sage: hs = [[h11,h12], [h21,h22]]
        sage: G = matrix.column(gs)
        sage: H = matrix.column(hs)
        sage: GG = G.transpose() * G
        sage: HH = H.transpose() * H
        sage: M = G.transpose() * H
        sage: I = [0, 1]
        sage: J = [0, 1]
        sage: all( v in [-1,1] and m == 2
        ....:      for (v,_,_,m) in _solve_gevp_naive(GG,HH,M,I,J) )
        True
    """
    A = matrix.block([
        [ZZ.zero(), M[I,J]],
        [M.transpose()[J,I], ZZ.zero()]
    ])
    B = matrix.block([
        [GG[I,I], ZZ.zero()],
        [ZZ.zero(), HH[J,J]]
    ])
    M = B.inverse() * A

    # We'll format the result to match the solve_gevp_nonzero() return value.
    for (evalue, evectors, multiplicity) in M.eigenvectors_right():
        for z in evectors:
            xi = z[0:len(I)]
            xi.set_immutable()
            eta = z[len(I):]
            eta.set_immutable()
            yield (evalue, xi, eta, multiplicity)


def solve_gevp_zero(M, I, J):
    r"""
    Solve the generalized eigenvalue problem in Theorem 3
    [Or2020]_ for a zero eigenvalue using Propositions 3 and 4
    [Or2020]_.

    INPUT:

    - ``M`` -- the matrix whose `(i,j)`-th entry is the inner product
      of `g_{i}` and `h_{j}` as in Proposition 6 [Or2020]_

    - ``I`` -- a linearly independent column-index set for the matrix
      `G` that appears in Theorem 3 [Or2020]_

    - ``J`` -- a linearly independent column-index set for the matrix
      `H` that appears in Theorem 3 [Or2020]_

    OUTPUT:

    A generator of ``(eigenvalue, xi, eta, multiplicity)`` quartets
    where

    - ``eigenvalue`` is zero (the eigenvalue of the system)

    - ``xi`` is the first (length ``len(I)``) component of an
      eigenvector associated with ``eigenvalue``

    - ``eta`` is the second (length ``len(J)``) component of an
      eigenvector associated with ``eigenvalue``

    - ``multiplicity`` is the dimension of the eigenspace associated
      with ``eigenvalue``

    ALGORITHM:

    Proposition 4 in [Or2020]_ is used.

    EXAMPLES:

    This particular configuration results in the zero matrix in the
    eigenvalue problem, so the only solutions correspond to the
    eigenvalue zero::

        sage: from sage.geometry.cone_critical_angles import solve_gevp_zero
        sage: K = cones.nonnegative_orthant(2)
        sage: G = matrix.column(K.rays())
        sage: GG = G.transpose() * G
        sage: I = [0]
        sage: J = [1]
        sage: list(solve_gevp_zero(GG, I, J))
        [(0, (1), (0), 2), (0, (0), (1), 2)]
    """
    # A Cartesian product would be more appropriate here, but Sage
    # isn't smart enough to figure out a basis for the product. So,
    # we use the direct sum and then chop it up.
    M_IJ = M[I,J]
    xi_space = M_IJ.left_kernel()
    eta_space = M_IJ.right_kernel()

    fake_cartprod = xi_space.direct_sum(eta_space)
    multiplicity = fake_cartprod.dimension()

    for z in fake_cartprod.basis():
        z1 = z[0:len(I)]
        z1.set_immutable()
        z2 = z[len(I):]
        z2.set_immutable()

        # The base ring of M will either be RDF or AA, which is enough
        # to contain any eigenvalues that will arise... meaning that
        # if we use the corresponding "zero" here, it will match the
        # field of the eigenvalues returned by the nonzero function.
        yield (M.base_ring().zero(), z1, z2, multiplicity)


def solve_gevp_nonzero(GG, HH, M, I, J):
    r"""
    Solve the generalized eigenvalue problem in Theorem 3
    [Or2020]_ for a nonzero eigenvalue using Propositions 3 and 5
    [Or2020]_.

    INPUT:

    - ``GG`` -- the matrix whose `(i,j)`-th entry is the inner product
      of `g_{i}` and `g_{j}`, which are in turn the `i`-th and `j`-th
      columns of the matrix `G` in Theorem 3 [Or2020]_

    - ``HH`` -- the matrix whose `(i,j)`-th entry is the inner product
      of `h_{i}` and `h_{j}`, which are in turn the `i`-th and `j`-th
      columns of the matrix `H` in Theorem 3 [Or2020]_

    - ``M`` -- the matrix whose `(i,j)`-th entry is the inner product
      of `g_{i}` and `h_{j}` as in Proposition 6 in [Or2020]_

    - ``I`` -- a linearly independent column-index set for the matrix
      `G` that appears in Theorem 3 [Or2020]_

    - ``J`` -- a linearly independent column-index set for the matrix
      `H` that appears in Theorem 3 [Or2020]_

    OUTPUT:

    A generator of ``(eigenvalue, xi, eta, multiplicity)`` quartets
    where

    - ``eigenvalue`` is a real eigenvalue of the system

    - ``xi`` is the first (length ``len(I)``) component of an
      eigenvector associated with ``eigenvalue``

    - ``eta`` is the second (length ``len(J)``) component of an
      eigenvector associated with ``eigenvalue``

    - ``multiplicity`` is the dimension of the eigenspace associated
      with ``eigenvalue``

    Note that we do not return a basis for each eigenspace along with
    its eigenvalue. For the application we have in mind, an eigenspace
    of dimension greater than one (so, ``multiplicity > 1``) is an
    error. As such, our return value is optimized for convenience in
    the non-error case, where there is only one eigenvector (spanning
    a one-dimensional eigenspace) associated with each eigenvalue.

    ALGORITHM:

    According to Proposition 5 [Or2020]_, the solutions corresponding
    to nonzero eigenvalues can be found by solving a smaller
    eigenvalue problem in only the variable `\xi`. So, we do that, and
    then solve for `\eta` in terms of `\xi` as described in the
    proposition.

    EXAMPLES:

    When the zero solutions are included, this function returns the
    same solutions as the naive method on the Schur cone in three
    dimensions::

        sage: from itertools import chain
        sage: from sage.geometry.cone_critical_angles import (
        ....:   _normalize_gevp_solution,
        ....:   _solve_gevp_naive,
        ....:   gevp_licis,
        ....:   solve_gevp_nonzero,
        ....:   solve_gevp_zero)
        sage: K = cones.schur(3)
        sage: gs = [g.change_ring(AA).normalized() for g in K]
        sage: G = matrix.column(gs)
        sage: GG = G.transpose() * G
        sage: G_index_sets = list(gevp_licis(G))
        sage: all(
        ....:   set(
        ....:     _normalize_gevp_solution(s)
        ....:     for s in
        ....:     chain(
        ....:       solve_gevp_zero(GG, I, J),
        ....:       solve_gevp_nonzero(GG, GG, GG, I, J)
        ....:     )
        ....:   )
        ....:   ==
        ....:   set(
        ....:     _normalize_gevp_solution(s)
        ....:     for s in
        ....:     _solve_gevp_naive(GG,GG,GG,I,J)
        ....:   )
        ....:   for I in G_index_sets
        ....:   for J in G_index_sets
        ....: )
        True

    TESTS:

    This function should return the same solutions (with zero included,
    of course) as the naive implementation even for random cones::

        sage: # long time
        sage: from itertools import chain
        sage: from sage.geometry.cone_critical_angles import (
        ....:   _normalize_gevp_solution,
        ....:   _random_admissible_cone,
        ....:   _solve_gevp_naive,
        ....:   gevp_licis,
        ....:   solve_gevp_nonzero,
        ....:   solve_gevp_zero)
        sage: n = ZZ.random_element(1,3)
        sage: P = _random_admissible_cone(ambient_dim=n)
        sage: Q = _random_admissible_cone(ambient_dim=n)
        sage: gs = [g.change_ring(AA).normalized() for g in P]
        sage: G = matrix.column(gs)
        sage: GG = G.transpose() * G
        sage: hs = [h.change_ring(AA).normalized() for h in Q]
        sage: H = matrix.column(hs)
        sage: HH = H.transpose() * H
        sage: M = G.transpose() * H
        sage: G_index_sets = list(gevp_licis(G))
        sage: H_index_sets = list(gevp_licis(H))
        sage: all(
        ....:   set(
        ....:     _normalize_gevp_solution(s)
        ....:     for s in
        ....:     chain(
        ....:       solve_gevp_zero(M, I, J),
        ....:       solve_gevp_nonzero(GG, HH, M, I, J)
        ....:     )
        ....:   )
        ....:   ==
        ....:   set(
        ....:     _normalize_gevp_solution(s)
        ....:     for s in
        ....:     _solve_gevp_naive(GG, HH, M, I, J)
        ....:   )
        ....:   for I in G_index_sets
        ....:   for J in H_index_sets
        ....: )
        True

    According to Proposition 7 [Or2020]_, the only eigenvalues that
    arise when either ``G`` or ``H`` is invertible are `-1`, `0`, and
    `1`::

        sage: # long time
        sage: from sage.geometry.cone_critical_angles import (
        ....:   _random_admissible_cone,
        ....:   gevp_licis,
        ....:   solve_gevp_nonzero)
        sage: n = ZZ.random_element(1,3)
        sage: P = _random_admissible_cone(ambient_dim=n)
        sage: Q = _random_admissible_cone(ambient_dim=n)
        sage: gs = [g.change_ring(AA).normalized() for g in P]
        sage: hs = [h.change_ring(AA).normalized() for h in Q]
        sage: G = matrix.column(gs)
        sage: GG = G.transpose() * G
        sage: H = matrix.column(hs)
        sage: HH = H.transpose() * H
        sage: M = G.transpose() * H
        sage: from itertools import product
        sage: all(
        ....:  (v in [-1,0,1]
        ....:   for (v,_,_,_) in solve_gevp_nonzero(GG, HH, M, I, J))
        ....:   for (I,J) in product(gevp_licis(G),gevp_licis(H))
        ....:   if len(I) == n or len(J) == n )
        True
    """
    if len(J) < len(I):
        # We can always opt to solve the smaller problem. Reading the
        # first three assignments below, you should be able to
        # convince yourself that switching GG <-> HH, I <-> J, and
        # transposing M does in fact switch from the "xi problem" to
        # the "eta problem."
        yield from ((l, xi, eta, m)
                    for (l, eta, xi, m)
                    in solve_gevp_nonzero(HH, GG, M.transpose(), J, I))
    else:
        M_IJ = M[I,J]
        G_I_pinv_H_J = GG[I,I].inverse_positive_definite() * M_IJ
        H_J_pinv_G_I = HH[J,J].inverse_positive_definite() * M_IJ.transpose()
        L = (G_I_pinv_H_J * H_J_pinv_G_I)

        for (sigma, xis, m) in L.eigenvectors_right():
            if sigma > 0:
                # Avoid recomputing these for each xi in xis
                sigma_sqrt = sigma.sqrt()
                inv_sqrt = ~sigma_sqrt
                pm_sqrt_inv_pairs = [
                    (-sigma_sqrt, -inv_sqrt),
                    (sigma_sqrt, inv_sqrt)
                ]

                for xi in xis:
                    for l, li in pm_sqrt_inv_pairs:
                        eta = li * H_J_pinv_G_I*xi
                        eta.set_immutable()
                        yield (l, xi, eta, m)


def compute_gevp_M(gs, hs):
    r"""
    Compute the matrix `M` whose `(i,j)`-th entry is the inner
    product of ``gs[i]`` and ``hs[j]``.

    This is the "generalized Gram matrix" appearing in Proposition 6
    in [Or2020]_. For efficiency, we also return the minimal pair,
    whose inner product is minimal among the entries of `M`. This
    allows our consumer to bail out immediately (knowing the optimal
    pair!) if it turns out that the maximal angle is acute; i.e. if
    the smallest entry of `M` is nonnegative.

    INPUT:

    - ``gs`` -- a linearly independent list of unit-norm generators
      for the cone `P`

    - ``hs`` -- a linearly independent list of unit-norm generators
      for the cone `Q`

    OUTPUT: a tuple containing four elements, in order:

    - The matrix `M` described in Proposition 6

    - The minimal entry in the matrix `M`

    - A vector in ``gs`` that achieves that minimal inner product
      along with the next element of the tuple

    - A vector in ``hs`` that achieves the minimal inner product
      along with the previous element in the tuple

    EXAMPLES::

        sage: from sage.geometry.cone_critical_angles import compute_gevp_M
        sage: P = Cone([ (1,2,0), (3,4,0) ])
        sage: Q = Cone([ (-1,4,1), (5,-2,-1),  (-1,-1,5) ])
        sage: gs = [g.change_ring(QQ) for g in P]
        sage: hs = [h.change_ring(QQ) for h in Q]
        sage: M = compute_gevp_M(gs, hs)[0]
        sage: all( M[i][j] == gs[i].inner_product(hs[j])
        ....:       for i in range(P.nrays())
        ....:       for j in range(Q.nrays()) )
        True

    TESTS:

    The products `(G_{I})^{T}H_{J}` correspond to
    submatrices of the "generalized Gram matrix" `M` in Proposition
    6. Note that SageMath does (row,column) indexing but [Or2020]_
    does (column,row) indexing::

        sage: from sage.geometry.cone_critical_angles import (
        ....:   _random_admissible_cone,
        ....:   compute_gevp_M,
        ....:   gevp_licis)
        sage: n = ZZ.random_element(1,4)
        sage: n = ZZ.random_element(1,8) # long time
        sage: P = _random_admissible_cone(ambient_dim=n)
        sage: Q = _random_admissible_cone(ambient_dim=n)
        sage: gs = [g.change_ring(QQ) for g in P]
        sage: hs = [h.change_ring(QQ) for h in Q]
        sage: M = compute_gevp_M(gs,hs)[0]
        sage: f = lambda i,j: gs[i].inner_product(hs[j])
        sage: expected_M = matrix(QQ, P.nrays(), Q.nrays(), f)
        sage: M == expected_M
        True
        sage: G = matrix.column(gs)
        sage: H = matrix.column(hs)
        sage: def _test_indexing(I,J):
        ....:      G_I = G.matrix_from_columns(I)
        ....:      H_J = H.matrix_from_columns(J)
        ....:      return (G_I.transpose()*H_J == M[I,J]
        ....:              and
        ....:              H_J.transpose()*G_I == M.transpose()[J,I])
        sage: G_index_sets = list(gevp_licis(G))
        sage: H_index_sets = list(gevp_licis(H))
        sage: all( _test_indexing(I,J) for I in G_index_sets
        ....:                          for J in H_index_sets )
        True
    """
    min_u = gs[0]
    min_v = hs[0]
    min_ip = min_u.inner_product(min_v)

    M = []
    for g in gs:
        M_i = []
        for h in hs:
            val = g.inner_product(h)
            M_i.append(val)
            if (val < min_ip):
                min_ip = val
                min_u = g
                min_v = h
        M.append(M_i)

    return (matrix(M), min_ip, min_u, min_v)


def check_gevp_feasibility(cos_theta, xi, eta, G_I, G_I_c_T,
                           H_J, H_J_c_T, epsilon):
    r"""
    Determine if a solution to the generalized eigenvalue problem
    in Theorem 3 [Or2020]_ is feasible.

    Implementation detail: we take four matrices that we are capable
    of computing as parameters instead, because we will be called in a
    nested loop "for all `I`... and for all `J`..." The data
    corresponding to `I` should be computed only once, which means
    that we can't do it here -- it needs to be done outside of the `J`
    loop. For symmetry (and to avoid relying on too many
    cross-function implementation details), we also insist that the
    `J` data be passed in.

    INPUT:

    - ``cos_theta`` -- an eigenvalue corresponding to
      `( \xi, \eta )`

    - ``xi`` -- first component of the `( \xi, \eta )` eigenvector

    - ``eta`` -- second component of the `( \xi, \eta )` eigenvector

    - ``G_I`` -- the submatrix of `G` with columns indexed by `I`

    - ``G_I_c_T`` -- a matrix whose rows are the non-`I` columns of `G`

    - ``H_J`` -- the submatrix of `H` with columns indexed by `J`

    - ``H_J_c_T`` -- a matrix whose rows are the non-`J` columns of `H`

    - ``epsilon`` -- the tolerance to use when making comparisons

    OUTPUT:

    A triple containing (in order),

    - a boolean,
    - a vector in the cone `P` (of the same length as ``xi``), and
    - a vector in the cone `Q` (of the same length as ``eta``).

    If `( \xi, \eta )` is feasible, we return ``(True, u, v)`` where `u`
    and `v` are the vectors in `P` and `Q` respectively that form the
    the angle `\theta`.

    If `( \xi, \eta )` is **not** feasible, then we return ``(False, 0, 0)``
    where ``0`` should be interpreted to mean the zero vector in the
    appropriate space.

    EXAMPLES:

    If `\xi` has any components less than "zero," it isn't feasible::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [-1,1])
        sage: eta = vector(QQ, [1,1,1])
        sage: check_gevp_feasibility(0,xi,eta,None,None,None,None,0)
        (False, (0, 0), (0, 0, 0))

    If `\eta` has any components less than "zero," it isn't feasible::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [2])
        sage: eta = vector(QQ, [1,-4,4,5])
        sage: check_gevp_feasibility(0,xi,eta,None,None,None,None,0)
        (False, (0), (0, 0, 0, 0))

    If `\xi` and `\eta` are equal and if `G_{I}` and `H_{J}` are not,
    then the copy of `\eta` that's been scaled by the norm of `G_{I}\xi`
    generally won't satisfy its norm-equality constraint::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [1,1])
        sage: eta = xi
        sage: G_I = matrix.identity(QQ,2)
        sage: H_J = 2*G_I
        sage: check_gevp_feasibility(0,xi,eta,G_I,None,H_J,None,0)
        (False, (0, 0), (0, 0))

    When `\cos\theta` is zero, the inequality (42) in Theorem 7.3
    [SS2016]_ is just an inner product with `v` which we can make
    positive by ensuring that all of the entries of `H_{J}` are
    positive. So, if any of the rows of ``G_I_c_T`` contain a negative
    entry, (42) will fail::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [1/2,1/2,1/2,1/2])
        sage: eta = xi
        sage: G_I = matrix.identity(QQ,4)
        sage: G_I_c_T = matrix(QQ, [[0,-1,0,0]])
        sage: H_J = G_I
        sage: check_gevp_feasibility(0,xi,eta,G_I,G_I_c_T,H_J,None,0)
        (False, (0, 0, 0, 0), (0, 0, 0, 0))

    Likewise we can make (43) fail in exactly the same way::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [1/2,1/2,1/2,1/2])
        sage: eta = xi
        sage: G_I = matrix.identity(QQ,4)
        sage: G_I_c_T = matrix(QQ, [[0,1,0,0]])
        sage: H_J = G_I
        sage: H_J_c_T = matrix(QQ, [[0,-1,0,0]])
        sage: check_gevp_feasibility(0,xi,eta,G_I,G_I_c_T,H_J,H_J_c_T,0)
        (False, (0, 0, 0, 0), (0, 0, 0, 0))

    Finally, if we ensure that everything works, we get back a feasible
    result along with the vectors (scaled `\xi` and `\eta`) that worked::

        sage: from sage.geometry.cone_critical_angles import(
        ....:   check_gevp_feasibility)
        sage: xi = vector(QQ, [1/2,1/2,1/2,1/2])
        sage: eta = xi
        sage: G_I = matrix.identity(QQ,4)
        sage: G_I_c_T = matrix(QQ, [[0,1,0,0]])
        sage: H_J = G_I
        sage: H_J_c_T = matrix(QQ, [[0,1,0,0]])
        sage: check_gevp_feasibility(0,xi,eta,G_I,G_I_c_T,H_J,H_J_c_T,0)
        (True, (1/2, 1/2, 1/2, 1/2), (1/2, 1/2, 1/2, 1/2))
    """
    infeasible_result = (False, 0*xi, 0*eta)
    if min(xi) <= -epsilon or min(eta) <= -epsilon:
        # xi or eta isn't in the interior of the nonnegative orthant,
        # so skip this (non-)solution.
        return infeasible_result

    # Rescale xi to satisfy (44), and rescale eta by the same amount,
    # because (xi,eta) needs to remain in the same one-dimensional
    # eigenspace.
    scale = ~((G_I*xi).norm())
    xi_hat = xi * scale
    eta_hat = eta * scale

    # Now check that (45) is satisfied.
    if ((H_J*eta_hat).norm() - 1).abs() > epsilon:
        return infeasible_result

    # And check that (42,43) are satisfied.
    v = H_J * eta_hat
    rhs = v - cos_theta*G_I*xi_hat

    if any(x < -epsilon for x in G_I_c_T * rhs):
        return infeasible_result

    u = G_I * xi_hat
    rhs = u - cos_theta*H_J*eta_hat
    if any(x < -epsilon for x in H_J_c_T * rhs):
        return infeasible_result

    return (True, u, v)


def max_angle(P, Q, exact, epsilon):
    r"""
    Find the maximal angle between the cones `P` and `Q`.

    This implements
    :meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.max_angle`,
    which should be fully documented.

    EXAMPLES:

    For the sake of the user interface, the argument validation for
    this function is performed in the associated cone method; we can
    therefore crash it by feeding it invalid input like an
    inadmissible cone::

        sage: from sage.geometry.cone_critical_angles import max_angle
        sage: K = cones.trivial(3)
        sage: max_angle(K,K,True,0)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range
    """
    # The lattice dimensions of P and Q are guaranteed to be equal
    # because the cone method checks it before calling us.
    n = P.lattice_dim()

    ring = RDF
    if exact:
        ring = AA
        # For some reason we can go RR -> QQ -> AA, but not
        # straight from RR to AA.
        epsilon = QQ(epsilon)
    epsilon = ring(epsilon)

    # First check if P is contained in the dual of Q. Keep track of
    # the minimum inner product (and associated vectors) while doing
    # so; then if P is contained in dual(Q), we just return the pair
    # with the smallest inner product.
    gs = [g.change_ring(ring).normalized() for g in P]
    Q_is_P = (P == Q) # This is used again later
    if Q_is_P:
        hs = gs
    else:
        hs = [h.change_ring(ring).normalized() for h in Q]

    (M, min_ip, min_u, min_v) = compute_gevp_M(gs,hs)

    if min_ip >= 0: # The maximal angle is acute!
        return (arccos(min_ip), min_u, min_v)

    # Also check to see if the maximal angle is pi. In particular this
    # is true when either P or Q is the entire ambient space.
    P_and_negative_Q = P.intersection(-Q)
    if not (P_and_negative_Q.is_trivial()):
        u = P_and_negative_Q.ray(0).change_ring(ring).normalized()
        v = -u
        return (pi, u, v)

    # When P == Q, GG and HH are both just M. We rule out the
    # cardinality ``n`` in the index sets because it will eventually
    # result in a GEVP whose only solutions are lambda in {-1,0,1};
    # none of which we want! We rule out 0 and 1 with the acute check,
    # and -1 with the pi (P_and_negative_Q) check.
    #
    # It's VERY IMPORTANT that we construct lists from the index set
    # generators, because we're going to use them in a nested loop!
    G = matrix.column(gs)
    G_index_sets = [s for s in gevp_licis(G) if not len(s) == n]

    if Q_is_P:
        GG = M
        H = G
        HH = M
        H_index_sets = G_index_sets
    else:
        GG = G.transpose() * G
        H = matrix.column(hs)
        HH = H.transpose() * H
        H_index_sets = [s for s in gevp_licis(H) if not len(s) == n]

    # Keep track of the (cos-theta, xi, eta, multiplicity) tuples with
    # multiplicity > 1. These are only a problem if they could
    # potentially be maximal. Therefore, we want to inspect them AFTER
    # we've checked all of the multiplicity=1 angles and found the
    # largest. This allows us to ignore most of the problematic
    # eigenspaces, independent of the order in which we run through I,J.
    big_eigenspaces = []

    for I in G_index_sets:
        G_I = G.matrix_from_columns(I)
        I_complement = [i for i in range(P.nrays()) if i not in I]
        G_I_c_T = G.matrix_from_columns(I_complement).transpose()

        for J in H_index_sets:
            J_complement = [j for j in range(Q.nrays()) if j not in J]
            H_J = H.matrix_from_columns(J)
            H_J_c_T = H.matrix_from_columns(J_complement).transpose()

            for (cos_theta,xi,eta,mult) in solve_gevp_nonzero(GG, HH, M, I, J):

                if cos_theta >= min_ip:
                    # This potential critical angle is smaller than or
                    # equal to one that we've already found. Why
                    # bother?
                    continue

                if cos_theta == -1:
                    # We already ruled this case out with the
                    # "P_and_negative_Q" trick.
                    continue

                (is_feasible, u, v) = check_gevp_feasibility(cos_theta,
                                                             xi,
                                                             eta,
                                                             G_I,
                                                             G_I_c_T,
                                                             H_J,
                                                             H_J_c_T,
                                                             epsilon)

                if is_feasible:
                    min_ip = cos_theta
                    min_u = u
                    min_v = v
                elif mult > 1:
                    # Save this for later. The eigenvalue cos_theta
                    # might be so big that we can ignore it.
                    big_eigenspaces.append((cos_theta, xi, eta, mult))
                    continue

    for (cos_theta, xi, eta, mult) in big_eigenspaces:
        if cos_theta < min_ip:
            # The existence of a big eigenspace is only a problem if
            # cos_theta could actually be minimal.

            if exact and P.is_strictly_convex():
                if Q_is_P or Q.is_strictly_convex():
                    # The maximal angle composed with the conic hull
                    # is continuous at (gs,hs), so we can retry with
                    # inexact arithmetic. That will "perturb"
                    # everything, hopefully eliminating any larger
                    # eigenspaces and without changing the answer too
                    # much.
                    return max_angle(P, Q, False, epsilon)

            # Either we don't know that the maximal angle of the conic
            # hull is continuous at (gs,hs), or we're already using
            # inexact arithmetic. There's nothing left to try. (Note
            # that the case where either P or Q is the ambient space
            # was handled much earlier, since in that case the maximal
            # angle is obviously pi.)
            raise ValueError('eigenspace of dimension %d > 1 '
                             'corresponding to eigenvalue %s'
                              % (mult, cos_theta))

    return (arccos(min_ip), min_u, min_v)
