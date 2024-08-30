# sage.doctest: needs sage.libs.singular
r"""
Helper Functions For Freeness Of Hyperplane Arrangements

This contains the algorithms to check for freeness of a hyperplane
arrangement. See
:meth:`sage.geometry.hyperplane_arrangement.HyperplaneArrangementElement.is_free`
for details.

.. NOTE::

    This could be extended to a freeness check for more general modules
    over a polynomial ring.
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
import sage.libs.singular.function_factory as fun_fact

def less_generators(X):
    """
    Reduce the generator matrix of the module defined by ``X``.

    This is Algorithm 6.4 in [BC2012]_ and relies on the row syzygies of
    the matrix ``X``.

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.check_freeness import less_generators
        sage: R.<x,y,z> = QQ[]
        sage: m = matrix([[1, 0, 0], [0, z, -1], [0, 0, 0], [0, y, 1]])
        sage: less_generators(m)
        [ 1  0  0]
        [ 0  z -1]
        [ 0  y  1]
    """
    R = X.base_ring()
    syz = fun_fact.ff.syz
    zero = R.zero()
    while True:
        Z = matrix(R, syz(X.transpose()))
        # get_column_independent_unit_positions
        J = range(Z.ncols())
        K = []
        for r in Z.rows():
            for j in J:
                elt = r[j]
                if elt.is_unit():
                    K.append(j)
                    J = [l for l in J if r[l] == zero]
                    break
        if not K: # K is empty
            return X
        Kd = set(range(X.nrows())).difference(K)
        X = X.matrix_from_rows(sorted(Kd))

def construct_free_chain(A):
    """
    Construct the free chain for the hyperplanes ``A``.

    ALGORITHM:

    We follow Algorithm 6.5 in [BC2012]_.

    INPUT:

    - ``A`` -- a hyperplane arrangement

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.check_freeness import construct_free_chain
        sage: H.<x,y,z> = HyperplaneArrangements(QQ)
        sage: A = H(z, y+z, x+y+z)
        sage: construct_free_chain(A)
        [
        [1 0 0]  [ 1  0  0]  [    0     1     0]
        [0 1 0]  [ 0  z -1]  [y + z     0    -1]
        [0 0 z], [ 0  y  1], [    x     0     1]
        ]
    """
    AL = list(A)
    if not AL: # Empty arrangement
        return []

    S = A.parent().ambient_space().symmetric_space()
    # Compute the morphisms \phi_{H_j} : S^{1xR} \to S / <b_j>
    B = [H.to_symmetric_space() for H in AL]
    phi = [matrix(S, [[beta.derivative(x)] for x in S.gens()]) for beta in B]

    # Setup variables
    syz = fun_fact.ff.syz
    G = S.gens()
    r = len(G)
    indices = list(range(len(B)))

    # Helper function
    def next_step(indices, prev, T):
        for pos,i in enumerate(indices):
            U = prev * T
            mu = U * phi[i]
            mu = mu.stack(matrix.diagonal([B[i]]).dense_matrix())
            row_syzygy = matrix(S, syz(mu.transpose())).matrix_from_columns(range(r))
            Y = less_generators(row_syzygy)
            if not Y.is_square():
                continue

            if len(indices) == 1:
                return [prev, Y]

            I = list(indices)
            I.pop(pos)
            ret = next_step(I, Y, U)
            if ret is not None:
                return [prev] + ret
        return None

    T = matrix.identity(S, r)
    for i in indices:
        mu = phi[i].stack(matrix.diagonal([B[i]]).dense_matrix())
        row_syzygy = matrix(S, syz(mu.transpose())).matrix_from_columns(range(r))
        Y = less_generators(row_syzygy)
        if not Y.is_square():
            continue

        if len(indices) == 1:
            return [Y]

        I = list(indices)
        I.pop(i)
        ret = next_step(I, Y, T)
        if ret is not None:
            return ret
    return None
