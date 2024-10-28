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
from sage.matrix.constructor import matrix, identity_matrix
from sage.modules.free_module_element import vector

from math import sqrt, floor, ceil


def plane_inequality(v):
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
            a_{i,j} - \sum_{r<i} q_{r,r} q_{r,i}^2 & i = j, \\
            0 & i > j,
        \end{cases}

    for all `1 \leq i \leq n` and `1 \leq j \leq n`. (These
    equalities determine the entries of `Q` uniquely by
    recursion.) This matrix `Q` is defined for all `M` in a
    certain Zariski-dense open subset of the set of all
    `n \times n`-matrices.

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


def diamond_cut(V, GM, C, verbose=False):
    r"""
    Perform diamond cutting on polyhedron ``V`` with basis matrix ``GM``
    and radius ``C``.

    INPUT:

    - ``V`` -- polyhedron to cut from

    - ``GM`` -- half of the basis matrix of the lattice

    - ``C`` -- radius to use in cutting algorithm

    - ``verbose`` -- boolean (default: ``False``); whether to print debug information

    OUTPUT: a :class:`Polyhedron` instance

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import diamond_cut
        sage: V = Polyhedron([[0], [2]])
        sage: GM = matrix([2])
        sage: V = diamond_cut(V, GM, 4)
        sage: V.vertices()
        (A vertex at (2), A vertex at (0))
    """
    if verbose:
        print("Cut\n{}\nwith radius {}".format(GM, C))

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
            print("Dimension: {}".format(i))
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
            print("x: {}".format(x))
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
    cut = Polyhedron(ieqs=inequalities)
    V = V.intersection(cut)

    if verbose:
        print("End")

    return V


def calculate_voronoi_cell(basis, radius=None, verbose=False):
    """
    Calculate the Voronoi cell of the lattice defined by basis.

    INPUT:

    - ``basis`` -- embedded basis matrix of the lattice

    - ``radius`` -- radius of basis vectors to consider

    - ``verbose`` -- whether to print debug information

    OUTPUT: a :class:`Polyhedron` instance

    EXAMPLES::

        sage: from sage.modules.diamond_cutting import calculate_voronoi_cell
        sage: V = calculate_voronoi_cell(matrix([[1, 0], [0, 1]]))
        sage: V.volume()
        1
    """
    dim = basis.dimensions()
    artificial_length = None
    if dim[0] < dim[1]:
        # introduce "artificial" basis points (representing infinity)
        def approx_norm(v):
            r, r1 = (v.inner_product(v)).sqrtrem()
            return r + (r1 > 0)
        artificial_length = max(approx_norm(v) for v in basis) * 2
        additional_vectors = identity_matrix(dim[1]) * artificial_length
        basis = basis.stack(additional_vectors)
        # LLL-reduce to get quadratic matrix
        basis = basis.LLL()
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

    # twice the length of longest vertex in Q is a safe choice
    if radius is None:
        radius = 2 * max(v.inner_product(v) for v in basis)

    V = diamond_cut(Q, basis, radius, verbose=verbose)

    if artificial_length is not None:
        # remove inequalities introduced by artificial basis points
        H = V.Hrepresentation()
        H = [v for v in H if all(not V._is_zero(v.A() * w / 2 - v.b()) and
                                 not V._is_zero(v.A() * (-w) / 2 - v.b())
                                 for w in additional_vectors)]
        V = Polyhedron(ieqs=H)

    return V
