# sage_setup: distribution = sagemath-modules
r"""
Special methods for matrices over discrete valuation rings/fields.
"""

# ****************************************************************************
#       Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#                          Raphaël Pagès <raphael.pages@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#          https://www.gnu.org/licenses/
# ****************************************************************************


from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.element cimport RingElement
from sage.rings.infinity import Infinity


# We assume that H is square
cpdef hessenbergize_cdvf(Matrix_generic_dense H):
    r"""
    Replace `H` with a Hessenberg form of it.

    .. NOTE::

        This function assumes that H is a matrix over
        a complete discrete valuation field.

        The pivot on each column is always chosen
        with maximal relative precision, which ensures
        the numerical stability of the algorithm.

    TESTS::

        sage: # needs sage.rings.padics
        sage: K = Qp(5, print_mode="digits", prec=5)
        sage: H = matrix(K, 3, 3, range(9))
        sage: H
        [        0  ...00001  ...00002]
        [ ...00003  ...00004 ...000010]
        [ ...00011  ...00012  ...00013]
        sage: H.hessenbergize()
        sage: H
        [        0  ...00010  ...00002]
        [ ...00003  ...00024 ...000010]
        [ ...00000  ...44440  ...44443]

    ::

        sage: M = random_matrix(K, 6, 6)                                                # needs sage.rings.padics
        sage: M.charpoly()[0] == M.determinant()                                        # needs sage.rings.padics
        True

    We check that :issue:`31753` is resolved::

        sage: R.<t> = GF(5)[[]]
        sage: M = matrix(3, 3, [ 1, t + O(t^3), t^2,
        ....:                    1 + t + O(t^3), 2 + t^2, 3 + 2*t + O(t^3),
        ....:                    t - t^2, 2*t, 1 + t ])
        sage: M.charpoly()
        x^3 + (1 + 4*t + 4*t^2 + O(t^3))*x^2 + (t + 2*t^2 + O(t^3))*x + 3 + 2*t^2 + O(t^3)
    """
    cdef Py_ssize_t n, i, j, k
    cdef RingElement entry, pivot, inv, scalar

    n = H.nrows()
    for j in range(n-1):
        k = j + 1
        #maxi = H.get_unsafe(k, j).precision_relative()
        maxi = - H.get_unsafe(k, j).valuation()
        i = j + 2
        while maxi is not Infinity and i < n:
            entry = H.get_unsafe(i, j)
            if entry:
                #m = entry.precision_relative()
                m = -entry.valuation()
                if m > maxi:
                    maxi = m
                    k = i
            i += 1

        if k != j + 1:
            H.swap_rows_c(j+1, k)
            H.swap_columns_c(j+1, k)
        pivot = H.get_unsafe(j+1, j)
        # print(pivot)
        if pivot:
            inv = (~pivot).lift_to_precision()
            for i in range(j+2, n):
                scalar = inv * H.get_unsafe(i, j)
                scalar = scalar.lift_to_precision()
                # print("  ", scalar)
                H.add_multiple_of_row_c(i, j+1, -scalar, 0)
                H.add_multiple_of_column_c(j+1, i, scalar, 0)
