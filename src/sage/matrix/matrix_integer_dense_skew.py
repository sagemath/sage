"""
Skew Smith normal form over ZZ
"""

from sage.matrix.constructor import Matrix
from sage.rings.integer_ring import ZZ
from sage.arith.misc import xgcd


def skew_form(M, transformation=False):
    r"""
    Compute the skew Smith normal form of this square integer matrix as defined
    in Theorem IV.1 in [New1972]_.

    Given an integer skew-symmetric matrix `M`, compute a unimodular matrix `S`
    such that `S^T M S` is block diagonal with `2x2` blocks `(0, d_i; -d_i, 0)`
    where `d_i` are positive integers with `d_i | d_{i+1}`. If ``transformation``
    is True, also return ``S``.

    INPUT:
    
    - ``M`` -- integer n x n skew-symmetric matrix

    - ``transformation`` -- boolean (default: ``False``); if ``True``,
      return transformation matrix ``S``

    OUTPUT:

    - ``M'`` -- integer n x n matrix in skew Smith normal form
    - ``S``  -- integer n x n unimodular matrix with `S^T M S = M'`
      (if ``transformation`` is True)

    EXAMPLES::

        sage: M = matrix(ZZ, [
        ....:     [   0,    8,    0,   -4,   -2,   12 ],
        ....:     [  -8,    0,   12,    0,    0,    0 ],
        ....:     [   0,  -12,    0,    6,    0,  -18 ],
        ....:     [   4,    0,   -6,    0,    0,    0 ],
        ....:     [   2,    0,    0,    0,    0,    0 ],
        ....:     [ -12,    0,   18,    0,    0,    0 ]
        ....: ])
        sage: M_prime, S = M.skew_form(transformation=True)
        sage: M_prime
        [ 0  2  0  0  0  0]
        [-2  0  0  0  0  0]
        [ 0  0  0  6  0  0]
        [ 0  0 -6  0  0  0]
        [ 0  0  0  0  0  0]
        [ 0  0  0  0  0  0]
        sage: bool(S.transpose() * M * S == M_prime)
        True
    """
    if M.ncols() != M.nrows():
        raise TypeError("Input must be a square matrix.")
    if M != -M.transpose():
        raise ValueError("Input matrix must be skew-symmetric.")

    A = M.change_ring(ZZ)
    n = A.nrows()
    S = Matrix.identity(ZZ, n) if transformation else None

    k = 0
    while k < n - 1:
        # We work on the submatrix A[k:, k:]
        m = n - k

        # Find a non-zero entry in the strict upper triangle of the submatrix
        piv_row, piv_col = None, None
        min_abs_val = float("inf")
        for r in range(m):
            for c in range(r + 1, m):
                val = A[k + r, k + c]
                if val != 0 and abs(val) < min_abs_val:
                    min_abs_val = abs(val)
                    piv_row, piv_col = r, c

        # If the submatrix is zero, we are done
        if piv_row is None:
            break

        # Move the pivot to the (k, k+1) position using swaps
        if piv_row != 0:
            E = Matrix.identity(ZZ, n)
            E.swap_columns(k, k + piv_row)
            A = E.transpose() * A * E
            if transformation:
                S = S * E

        if piv_col != 1:
            E = Matrix.identity(ZZ, n)
            E.swap_columns(k + 1, k + piv_col)
            A = E.transpose() * A * E
            if transformation:
                S = S * E

        # GCD descent loop
        while True:
            # Step 1: Clear the first row of the submatrix (row k)
            for j in range(2, m):
                g, h = A[k, k + 1], A[k, k + j]
                if h == 0:
                    continue
                d, x, y = xgcd(g, h)
                g_inv, h_inv = g // d, h // d

                # Create the elementary matrix for the GCD step
                E = Matrix.identity(ZZ, n)
                i1, i2 = k + 1, k + j
                E[i1, i1], E[i2, i1] = x, y
                E[i1, i2], E[i2, i2] = -h_inv, g_inv

                # Apply the congruence transformation
                A = E.transpose() * A * E
                if transformation:
                    S = S * E

            # Step 2: Check for divisibility
            h = A[k, k + 1]
            if h == 0:
                break

            bad_entry_row = None
            for r in range(1, m):
                for c in range(r, m):
                    if A[k + r, k + c] % h != 0:
                        bad_entry_row = r
                        break
                if bad_entry_row is not None:
                    break

            if bad_entry_row is None:
                break

            # Perturb to improve divisibility
            r_idx = k + bad_entry_row
            E = Matrix.identity(ZZ, n)
            E[r_idx, k] = 1

            # Apply the congruence transformation
            A = E.transpose() * A * E
            if transformation:
                S = S * E

        h = A[k, k + 1]
        if h == 0:
            continue

        # Step 3: Extract the 2x2 block by zeroing out the rest of rows/cols k, k+1
        for j in range(2, m):
            val = A[k + 1, k + j]
            if val == 0:
                continue
            q = val // h

            E = Matrix.identity(ZZ, n)
            E[k, k + j] = q

            # Apply the congruence transformation
            A = E.transpose() * A * E
            if transformation:
                S = S * E

        k += 2

    # Final normalisation to make diagonal entries positive
    for i in range(n // 2):
        if A[2 * i, 2 * i + 1] < 0:
            idx = 2 * i + 1
            # Create elementary matrix to rescale a column by -1
            E = Matrix.identity(ZZ, n)
            E[idx, idx] = -1

            # Apply the congruence transformation
            A = E.transpose() * A * E
            if transformation:
                S = S * E

    if transformation:
        return A, S
    return A