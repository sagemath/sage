# sage.doctest: needs sage.geometry.polyhedron
"""
Diamond cutting implementation

AUTHORS:

- Jan Poeschko (2012-07-02): initial version
"""
# ****************************************************************************
#       Copyright (C) 2012 Jan Poeschko <jan@poeschko.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from math import sqrt, floor, ceil


def plane_inequality(v) -> list:
    """
    Return the inequality for points on the same side as the origin
    with respect to the plane through ``v`` normal to ``v``.

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import plane_inequality
        sage: ieq = plane_inequality([1, -1]); ieq
        [2, -1, 1]
        sage: ieq[0] + vector(ieq[1:]) * vector([1, -1])
        0
    """
    v = vector(v)
    c = -v * v
    if c < 0:
        c, v = -c, -v
    return [c] + list(v)


def jacobi(M):
    r"""
    Compute the upper-triangular part of the Cholesky/Jacobi
    decomposition of the symmetric matrix ``M``.

    Let `M` be a symmetric `n \times n`-matrix over a field `F`.
    Let `m_{i,j}` denote the `(i,j)`-th entry of `M` for any
    `1 \leq i \leq n` and `1 \leq j \leq n`. Then, the
    upper-triangular part computed by this method is the
    upper-triangular `n \times n`-matrix `Q` whose
    `(i,j)`-th entry `q_{i,j}` satisfies

    .. MATH::

        q_{i,j} =
        \begin{cases}
            \frac{1}{q_{i,i}} \left( m_{i,j} - \sum_{r<i} q_{r,r} q_{r,i} q_{r,j} \right) & i < j, \\
            m_{i,j} - \sum_{r<i} q_{r,r} q_{r,i}^2 & i = j, \\
            0 & i > j,
        \end{cases}

    for all `1 \leq i \leq n` and `1 \leq j \leq n`. (These
    equalities determine the entries of `Q` uniquely by
    recursion.) This matrix `Q` is defined for every invertible
    `n \times n`-matrix `M`. Its definition is taken from (2.3)
    of [FP1985]_.

    .. NOTE::

        This should be a method of matrices.

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import jacobi
        sage: jacobi(identity_matrix(3) * 4)
        [4 0 0]
        [0 4 0]
        [0 0 4]

        sage: def testall(M):
        ....:      Q = jacobi(M)
        ....:      for j in range(3):
        ....:          for i in range(j):
        ....:              if Q[i,j] * Q[i,i] != M[i,j] - sum(Q[r,i] * Q[r,j] * Q[r,r] for r in range(i)):
        ....:                  return False
        ....:      for i in range(3):
        ....:          if Q[i,i] != M[i,i] - sum(Q[r,i] ** 2 * Q[r,r] for r in range(i)):
        ....:              return False
        ....:          for j in range(i):
        ....:              if Q[i,j] != 0:
        ....:                  return False
        ....:      return True

        sage: M = Matrix(QQ, [[8,1,5], [1,6,0], [5,0,3]])
        sage: Q = jacobi(M); Q
        [    8   1/8   5/8]
        [    0  47/8 -5/47]
        [    0     0 -9/47]
        sage: testall(M)
        True

        sage: M = Matrix(QQ, [[3,6,-1,7],[6,9,8,5],[-1,8,2,4],[7,5,4,0]])
        sage: testall(M)
        True
    """
    if not M.is_square():
        raise ValueError("the matrix must be square")
    dim = M.nrows()
    q = [list(row) for row in M]
    for i in range(dim - 1):
        for j in range(i + 1, dim):
            q[j][i] = q[i][j]
            q[i][j] = q[i][j] / q[i][i]
        for k in range(i + 1, dim):
            for l in range(k, dim):
                q[k][l] -= q[k][i] * q[i][l]
    for i in range(1, dim):
        for j in range(i):
            q[i][j] = 0
    return matrix(q)


def diamond_cut(V, GM, C, verbose=False) -> Polyhedron:
    r"""
    Perform diamond cutting on polyhedron ``V`` with basis matrix ``GM``
    and squared radius ``C``.

    INPUT:

    - ``V`` -- polyhedron to cut from

    - ``GM`` -- half of the basis matrix of the lattice

    - ``C`` -- square of the radius to use in cutting algorithm

    - ``verbose`` -- boolean (default: ``False``); whether to print
      debug information

    OUTPUT: a :class:`Polyhedron` instance

    ALGORITHM:

    Use the algorithm in (2.8) of [FP1985]_ to iterate through the nonzero
    vectors ``hv`` of length at most `\sqrt{C}` in the lattice spanned by
    ``GM``. (Actually, the algorithm only constructs one vector from each pair
    ``{hv, -hv}``.) For each such vector ``hv``, intersect ``V`` with the
    half-spaces defined by ``plane_inequality(hv)`` and
    ``plane_inequality(-hv)``.

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import diamond_cut
        sage: V = Polyhedron([[0], [2]])
        sage: GM = matrix([2])
        sage: V = diamond_cut(V, GM, 4)
        sage: V.vertices()
        (A vertex at (2), A vertex at (0))

    TESTS:

    Verify that code works when no cuts are performed::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: v = vector(ZZ, [1,1,-1])
        sage: L = IntegerLattice([v])
        sage: C = L.voronoi_cell(radius=0.1)
    """
    if verbose:
        print("Cut\n{}\nwith squared radius {}".format(GM, C))

    dim = GM.dimensions()
    if dim[0] != dim[1]:
        raise ValueError("the matrix must be square")
    dim = dim[0]
    T = [0] * dim
    U = [0] * dim
    x = [0] * dim
    L = [0] * dim

    # calculate the Gram matrix
    q = matrix([[sum(GM[i][k] * GM[j][k] for k in range(dim))
                 for j in range(dim)] for i in range(dim)])
    if verbose:
        print("q:\n{}".format(q.n()))
    # apply Cholesky/Jacobi decomposition
    q = jacobi(q)
    if verbose:
        print("q:\n{}".format(q.n()))

    i = dim - 1
    T[i] = C
    U[i] = 0

    new_dimension = True
    cut_count = 0
    inequalities = []
    while True:
        if verbose:
            print(f"Dimension: {i}")
        if new_dimension:
            Z = sqrt(T[i] / q[i][i])
            if verbose:
                print("Z: {}".format(Z))
            L[i] = int(floor(Z - U[i]))
            if verbose:
                print("L: {}".format(L))
            x[i] = int(ceil(-Z - U[i]) - 1)
            new_dimension = False

        x[i] += 1
        if verbose:
            print(f"x: {x}")
        if x[i] > L[i]:
            i += 1
        elif i > 0:
            T[i - 1] = T[i] - q[i][i] * (x[i] + U[i]) ** 2
            i -= 1
            U[i] = 0
            for j in range(i + 1, dim):
                U[i] += q[i][j] * x[j]
            new_dimension = True
        else:
            if all(elmt == 0 for elmt in x):
                break
            hv = [0] * dim
            for k in range(dim):
                for j in range(dim):
                    hv[k] += x[j] * GM[j][k]
            hv = vector(hv)

            for hv in [hv, -hv]:
                cut_count += 1
                if verbose:
                    print("\n%d) Cut using normal vector %s" % (cut_count, hv))
                inequalities.append(plane_inequality(hv))

    if verbose:
        print("Final cut")
    if inequalities:
        cut = Polyhedron(ieqs=inequalities)
        V = V.intersection(cut)

    if verbose:
        print("End")

    return V


def calculate_voronoi_cell(basis, radius=None, verbose=False) -> Polyhedron:
    """
    Calculate the Voronoi cell of the lattice defined by basis.

    INPUT:

    - ``basis`` -- embedded basis matrix of the lattice

    - ``radius`` -- square of radius of basis vectors to consider

    - ``verbose`` -- whether to print debug information

    OUTPUT: a :class:`Polyhedron` instance

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import calculate_voronoi_cell
        sage: V = calculate_voronoi_cell(matrix([[1, 0], [0, 1]]))
        sage: V.volume()
        1

    TESTS:

    Verify that :issue:`39507` is fixed::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: v = vector(ZZ, [1,1,1,-1])
        sage: L = IntegerLattice([v])
        sage: print(v in L)
        True
        sage: print(L.closest_vector(v))
        (1, 1, 1, -1)
        sage: C = L.voronoi_cell()
        sage: C.Hrepresentation()
        (An inequality (-1, -1, -1, 1) x + 2 >= 0,
         An inequality (1, 1, 1, -1) x + 2 >= 0)
        sage: v = vector(ZZ, [1,1,-1])
        sage: L = IntegerLattice([v])
        sage: C = L.voronoi_cell()
        sage: C.Hrepresentation()
        (An inequality (-2, -2, 2) x + 3 >= 0,
         An inequality (2, 2, -2) x + 3 >= 0)
        sage: C.Vrepresentation()
        (A line in the direction (0, 1, 1),
         A line in the direction (1, 0, 1),
         A vertex at (0, 0, -3/2),
         A vertex at (0, 0, 3/2))

    Verify that :issue:`37086` is fixed::

        sage: from sage.modules.free_module_integer import IntegerLattice
        sage: l  = [7, 0, -1, -2, -1, -2, 7, -2, 0, 0, -2,
        ....:       0, 7, -2, 0, -1, -2, -1, 7, 0 , -1, -1, 0, -2, 7]
        sage: M = matrix(5, 5, l)
        sage: C = IntegerLattice(M).voronoi_cell()
        sage: C
        A 5-dimensional polyhedron in QQ^5 defined as the
        convex hull of 720 vertices
    """
    dim = basis.dimensions()
    # LLL-reduce for efficiency.
    basis = basis.LLL()
    if radius is None:
        # Convert the basis matrix to use RDF numbers for efficiency when we
        # calculate the triangular matrix of the QR decomposition.
        from sage.rings.real_double import RDF
        transposed_RDF_matrix = (basis.transpose()).change_ring(RDF)
        R = transposed_RDF_matrix.QR()[1]
        # The length of the vector formed by the diagonal entries of R is an
        # upper bound for twice the covering radius, so it is an upper bound
        # on the length of the lattice vectors that need to be considered for
        # diamond cutting. However, the value of the `radius` keyword is
        # actually a squared length, so there is no square root in the
        # following formula.
        radius = sum(R[i, i]**2 for i in range(dim[0]))
        # We then divide by 4 as we will divide the basis by 2 later on.
        radius = ceil(radius / 4)
    artificial_length = None
    if dim[0] < dim[1]:
        F = basis.base_ring().fraction_field()
        # Introduce "artificial" basis points (representing infinity).
        additional_vectors = (F**dim[1]).subspace(basis).complement().basis()
        additional_vectors = matrix(additional_vectors)
        # LLL-reduce for efficiency.
        additional_vectors = additional_vectors.LLL()

        from sage.rings.real_double import RDF
        # Convert the basis matrix to use RDF numbers for efficiency when we
        # perform the QR decomposition.
        transposed_RDF_matrix = additional_vectors.transpose().change_ring(RDF)
        R = transposed_RDF_matrix.QR()[1]
        # Since R is triangular, its smallest diagonal entry provides a
        # lower bound on the length of the shortest nonzero vector in the
        # lattice spanned by the artificial points. We square it because
        # value of `radius` is a squared length.
        shortest_vector_lower_bound = min(R[i, i]**2
                                          for i in range(dim[1] - dim[0]))
        # We will multiply our artificial points by the following scalar in
        # order to make sure the squared length of the shortest
        # nonzero vector is greater than radius, even after the vectors
        # are divided by 2.
        artificial_length = ceil(2.001 * sqrt(radius / shortest_vector_lower_bound))
        additional_vectors *= artificial_length
        basis = basis.stack(additional_vectors)
        basis = matrix([v for v in basis if v])
        dim = basis.dimensions()
    if dim[0] != dim[1]:
        raise ValueError("invalid matrix")
    basis = basis / 2

    ieqs = []
    for v in basis:
        ieqs.append(plane_inequality(v))
        ieqs.append(plane_inequality(-v))
    Q = Polyhedron(ieqs=ieqs)

    V = diamond_cut(Q, basis, radius, verbose=verbose)

    if artificial_length is not None:
        # Remove inequalities introduced by artificial basis points.
        H = V.Hrepresentation()
        H = [v for v in H if all(not V._is_zero(v.A() * w / 2 - v.b()) and
                                 not V._is_zero(v.A() * (-w) / 2 - v.b())
                                 for w in additional_vectors)]
        V = Polyhedron(ieqs=H)

    return V
