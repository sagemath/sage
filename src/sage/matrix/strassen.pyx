"""
Generic Asymptotically Fast Strassen Algorithms

This implements asymptotically fast echelon form and matrix
multiplication algorithms.
"""

# ***************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.matrix.matrix_window cimport MatrixWindow

from cysignals.signals cimport sig_on, sig_off


def strassen_window_multiply(C, A,B, cutoff):
    """
    Multiply the submatrices specified by A and B, places result in
    C. Assumes that A and B have compatible dimensions to be
    multiplied, and that C is the correct size to receive the product,
    and that they are all defined over the same ring.

    Uses Strassen multiplication at high levels and then uses
    MatrixWindow methods at low levels.

    EXAMPLES: The following matrix
    dimensions are chosen especially to exercise the eight possible
    parity combinations that could occur while subdividing the matrix
    in the Strassen recursion. The base case in both cases will be a
    (4x5) matrix times a (5x6) matrix.

    ::

        sage: A = MatrixSpace(Integers(2^65), 64, 83).random_element()
        sage: B = MatrixSpace(Integers(2^65), 83, 101).random_element()
        sage: A._multiply_classical(B) == A._multiply_strassen(B, 3) #indirect doctest
        True

    AUTHORS:

    - David Harvey
    - Simon King (2011-07): Improve memory efficiency; :issue:`11610`
    """
    strassen_window_multiply_c(C, A, B, cutoff)


cdef strassen_window_multiply_c(MatrixWindow C, MatrixWindow A,
                                MatrixWindow B, Py_ssize_t cutoff):
    # todo -- I'm not sure how to interpret "cutoff". Should it be...
    # (a) the minimum side length of the matrices (currently implemented below)
    # (b) the maximum side length of the matrices
    # (c) the total number of entries being multiplied
    # (d) something else entirely?

    cdef Py_ssize_t A_nrows, A_ncols, B_ncols
    A_nrows = A._nrows
    A_ncols = A._ncols        # this should also be the number of rows of B
    B_ncols = B._ncols

    if (A_nrows <= cutoff) or (A_ncols <= cutoff) or (B_ncols <= cutoff):
        # note: this code is only reached if the TOP level is already beneath
        # the cutoff. In a typical large multiplication, the base case is
        # handled directly (see below).
        C.set_to_prod(A, B)
        return

    # Construct windows for the four quadrants of each matrix.
    # Note that if the side lengths are odd we're ignoring the
    # final row/column for the moment.

    cdef Py_ssize_t A_sub_nrows, A_sub_ncols, B_sub_ncols
    A_sub_nrows = A_nrows >> 1
    A_sub_ncols = A_ncols >> 1     # this is also like B_sub_nrows
    B_sub_ncols = B_ncols >> 1

    cdef bint have_cutoff = (A_sub_nrows <= cutoff) or (A_sub_ncols <= cutoff) or (B_sub_ncols <= cutoff)

    cdef MatrixWindow A00, A01, A10, A11, B00, B01, B10, B11
    A00 = A.matrix_window(0,           0,           A_sub_nrows, A_sub_ncols)
    A01 = A.matrix_window(0,           A_sub_ncols, A_sub_nrows, A_sub_ncols)
    A10 = A.matrix_window(A_sub_nrows, 0,           A_sub_nrows, A_sub_ncols)
    A11 = A.matrix_window(A_sub_nrows, A_sub_ncols, A_sub_nrows, A_sub_ncols)
    B00 = B.matrix_window(0,           0,           A_sub_ncols, B_sub_ncols)
    B01 = B.matrix_window(0,           B_sub_ncols, A_sub_ncols, B_sub_ncols)
    B10 = B.matrix_window(A_sub_ncols, 0,           A_sub_ncols, B_sub_ncols)
    B11 = B.matrix_window(A_sub_ncols, B_sub_ncols, A_sub_ncols, B_sub_ncols)

    # Allocate temp space.
    # S.K: We already have allocated C, so, we should use it for temporary results.
    # We use the schedule from Douglas-Heroux-Slishman-Smith (see also Boyer-Pernet-Zhou,
    # "Memory efficient scheduling of Strassen-Winograd's matrix multiplication algorithm",
    # Table 1).

    cdef MatrixWindow S0, S1, S2, S3, T0, T1, T2, T3, P0, P1, P2, P3, P4, P5, P6, U0, U1, U2, U3, U4, U5, U6
    cdef MatrixWindow X, Y
    X = A.new_empty_window(A_sub_nrows, max(A_sub_ncols, B_sub_ncols))
    Y = B.new_empty_window(A_sub_ncols, B_sub_ncols)

    # 1 S2 = A00-A10 in X
    S2 = X.matrix_window(0, 0, A_sub_nrows, A_sub_ncols)
    S2.set_to_diff(A00, A10)

    # 2 T2 = B11-B01 in Y
    T2 = Y
    T2.set_to_diff(B11, B01)

    # 3 P6 = S2*T2 in C10
    P6 = C.matrix_window(A_sub_nrows, 0, A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P6.set_to_prod(S2, T2)
    else:
        strassen_window_multiply_c(P6, S2, T2, cutoff)

    # 4 S0 = A10+A11 in X
    S0 = X.matrix_window(0, 0, A_sub_nrows, A_sub_ncols)
    S0.set_to_sum(A10, A11)

    # 5 T0 = B01-B00 in Y
    T0 = Y
    T0.set_to_diff(B01, B00)

    # 6 P4 = S0*T0 in C11
    P4 = C.matrix_window(A_sub_nrows, B_sub_ncols, A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P4.set_to_prod(S0, T0)
    else:
        strassen_window_multiply_c(P4, S0, T0, cutoff)

    # 7 S1 = S0-A00 in X
    S1 = X.matrix_window(0, 0, A_sub_nrows, A_sub_ncols)
    S1.set_to_diff(S0, A00)

    # 8 T1 = B11-T0 in Y
    T1 = Y
    T1.set_to_diff(B11,T0)

    # 9 P5 = S1*T1 in C01
    P5 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P5.set_to_prod(S1, T1)
    else:
        strassen_window_multiply_c(P5, S1, T1, cutoff)

    #10 S3 = A01-S1 in X
    S3 = X.matrix_window(0, 0, A_sub_nrows, A_sub_ncols)
    S3.set_to_diff(A01,S1)

    #11 P2 = S3*B11 in C00
    P2 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P2.set_to_prod(S3, B11)
    else:
        strassen_window_multiply_c(P2, S3, B11, cutoff)

    #12 P0 = A00*B00 in X
    P0 = X.matrix_window(0, 0, A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P0.set_to_prod(A00, B00)
    else:
        strassen_window_multiply_c(P0, A00, B00, cutoff)

    #13 U1 = P0+P5 in C01
    U1 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U1.set_to_sum(P0, P5)

    #14 U2 = U1+P6 in C10
    U2 = C.matrix_window(A_sub_nrows, 0,           A_sub_nrows, B_sub_ncols)
    U2.set_to_sum(U1, P6)

    #15 U3 = U1+P4 in C01
    U3 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U3.set_to_sum(U1, P4)

    #16 U6 = U2+P4 in C11 (final)
    U6 = C.matrix_window(A_sub_nrows, B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U6.set_to_sum(U2, P4)

    #17 U4 = U3+P2 in C01 (final)
    U4 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U4.set_to_sum(U3, P2)

    #18 T3 = T1-B10 in Y
    T3 = Y
    T3.set_to_diff(T1, B10)

    #19 P3 = A11*T3 in C00
    P3 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P3.set_to_prod(A11, T3)
    else:
        strassen_window_multiply_c(P3, A11, T3, cutoff)

    #20 U5 = U2-P3 in C10 (final)
    U5 = C.matrix_window(A_sub_nrows, 0,           A_sub_nrows, B_sub_ncols)
    U5.set_to_diff(U2, P3)

    #21 P1 = A01*B10 in C00
    P1 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    if have_cutoff:
        P1.set_to_prod(A01, B10)
    else:
        strassen_window_multiply_c(P1, A01, B10, cutoff)

    #22 U0 = P0+P1 in C00 (final)
    U0 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    U0.set_to_sum(P0, P1)

    # Now deal with the leftover row and/or column (if they exist).

    cdef MatrixWindow B_last_col, C_last_col, B_bulk, A_last_row, C_last_row, B_last_row, A_last_col, C_bulk

    if B_ncols & 1:
        B_last_col = B.matrix_window(0, B_ncols-1, A_ncols, 1)
        C_last_col = C.matrix_window(0, B_ncols-1, A_nrows, 1)
        C_last_col.set_to_prod(A, B_last_col)

    if A_nrows & 1:
        A_last_row = A.matrix_window(A_nrows-1, 0, 1, A_ncols)
        if B_ncols & 1:
            B_bulk = B.matrix_window(0, 0, A_ncols, B_ncols-1)
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols-1)
        else:
            B_bulk = B
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols)
        C_last_row.set_to_prod(A_last_row, B_bulk)

    if A_ncols & 1:
        A_last_col = A.matrix_window(0, A_ncols-1, A_sub_nrows << 1, 1)
        B_last_row = B.matrix_window(A_ncols-1, 0, 1, B_sub_ncols << 1)
        C_bulk = C.matrix_window(0, 0, A_sub_nrows << 1, B_sub_ncols << 1)
        C_bulk.add_prod(A_last_col, B_last_row)

cdef subtract_strassen_product(MatrixWindow result, MatrixWindow A, MatrixWindow B, Py_ssize_t cutoff):
    cdef MatrixWindow to_sub
    if (cutoff == -1 or result.ncols() <= cutoff or result.nrows() <= cutoff):
        result.subtract_prod(A, B)
    else:
        to_sub = A.new_empty_window(result.nrows(), result.ncols())
        strassen_window_multiply_c(to_sub, A, B, cutoff)
        result.subtract(to_sub)


def strassen_echelon(MatrixWindow A, cutoff):
    """
    Compute echelon form, in place. Internal function, call with
    M.echelonize(algorithm='strassen')

    Based on work of Robert Bradshaw
    and David Harvey at MSRI workshop in 2006.

    INPUT:

    - ``A`` -- matrix window

    - ``cutoff`` -- size at which algorithm reverts to
      naive Gaussian elimination and multiplication must be at least 1

    OUTPUT: the list of pivot columns

    EXAMPLES::

        sage: A = matrix(QQ, 7, [5, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 3, 1, 0, -1, 0, 0, -1, 0, 1, 2, -1, 1, 0, -1, 0, 1, 3, -1, 1, 0, 0, -2, 0, 2, 0, 1, 0, 0, -1, 0, 1, 0, 1])
        sage: B = A.__copy__(); B._echelon_strassen(1); B
        [ 1  0  0  0  0  0  0]
        [ 0  1  0 -1  0  1  0]
        [ 0  0  1  0  0  0  0]
        [ 0  0  0  0  1  0  0]
        [ 0  0  0  0  0  0  1]
        [ 0  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0]
        sage: C = A.__copy__(); C._echelon_strassen(2); C == B
        True
        sage: C = A.__copy__(); C._echelon_strassen(4); C == B
        True

    ::

        sage: n = 32; A = matrix(Integers(389),n,range(n^2))
        sage: B = A.__copy__(); B._echelon_in_place_classical()
        sage: C = A.__copy__(); C._echelon_strassen(2)
        sage: B == C
        True

    TESTS::

        sage: A = matrix(Integers(7), 4, 4, [1,2,0,3,0,0,1,0,0,1,0,0,0,0,0,1])
        sage: B = A.__copy__(); B._echelon_in_place_classical()
        sage: C = A.__copy__(); C._echelon_strassen(2)
        sage: B == C
        True

    ::

        sage: A = matrix(Integers(7), 4, 4, [1,0,5,0,2,0,3,6,5,1,2,6,4,6,1,1])
        sage: B = A.__copy__(); B._echelon_in_place_classical()
        sage: C = A.__copy__(); C._echelon_strassen(2)   #indirect doctest
        sage: B == C
        True

    AUTHORS:

    - Robert Bradshaw
    """
    if cutoff < 1:
        raise ValueError("cutoff must be at least 1")
    sig_on()
    strassen_echelon_c(A, cutoff, A._matrix._strassen_default_cutoff(A._matrix))
    sig_off()


cdef strassen_echelon_c(MatrixWindow A, Py_ssize_t cutoff, Py_ssize_t mul_cutoff):
    # The following notation will be used in the comments below, which should be understood to give
    # the general idea of what's going on, as if there were no inconvenient non-pivot columns.
    # The original matrix is given by [ A B ]
    #                                 [ C D ]
    # For compactness, let A' denote the inverse of A
    # top_left, top_right, bottom_left, and bottom_right loosely correspond to A, B, C, and D respectively,
    # however, the "cut" between the top and bottom rows need not be the same.

    cdef Py_ssize_t nrows, ncols
    nrows = A.nrows()
    ncols = A.ncols()

    if nrows <= cutoff or ncols <= cutoff:
        return list(A.echelon_in_place())

    cdef Py_ssize_t top_h, bottom_cut, bottom_h, bottom_start, top_cut
    cdef Py_ssize_t prev_pivot_count
    cdef Py_ssize_t split
    split = nrows // 2

    cdef MatrixWindow top, bottom, top_left, top_right, bottom_left, bottom_right, clear

    top = A.matrix_window(0, 0, split, ncols)
    bottom = A.matrix_window(split, 0, nrows-split, ncols)

    top_pivots = strassen_echelon_c(top, cutoff, mul_cutoff)
    # effectively "multiplied" top row by A^{-1}
    #     [  I  A'B ]
    #     [  C   D  ]

    top_pivot_intervals = int_range(top_pivots)
    top_h = len(top_pivots)

    if top_h == 0:
        #                                 [ 0 0 ]
        # the whole top is a zero matrix, [ C D ]. Run echelon on the bottom
        bottom_pivots = strassen_echelon_c(bottom, cutoff, mul_cutoff)
        #             [  0   0  ]
        # we now have [  I  C'D ], proceed to sorting

    else:
        bottom_cut = max(top_pivots) + 1
        bottom_left = bottom.matrix_window(0, 0, nrows-split, bottom_cut)

        if top_h == ncols:
            bottom.set_to_zero()
            # [ I ]
            # [ 0 ]
            # proceed to sorting

        else:
            if bottom_cut == top_h:
                clear = bottom_left
            else:
                clear = bottom_left.to_matrix().matrix_from_columns(top_pivots).matrix_window() # TODO: read only, can I do this faster? Also below
            # Subtract off C time top from the bottom_right
            if bottom_cut < ncols:
                bottom_right = bottom.matrix_window(0, bottom_cut, nrows-split, ncols-bottom_cut)
                subtract_strassen_product(bottom_right, clear, top.matrix_window(0, bottom_cut, top_h, ncols-bottom_cut), mul_cutoff)
            # [  I      A'B   ]
            # [  *   D - CA'B ]

            # Now subtract off C times the top from the bottom_left (pivots -> 0)
            if bottom_cut == top_h:
                bottom_left.set_to_zero()
                bottom_start = bottom_cut

            else:
                for cols in top_pivot_intervals:
                    bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]).set_to_zero()
                non_pivots = int_range(0, bottom_cut) - top_pivot_intervals
                for cols in non_pivots:
                    if cols[0] == 0: continue
                    prev_pivot_count = len(top_pivot_intervals - int_range(cols[0]+cols[1], bottom_cut - cols[0]+cols[1]))
                    subtract_strassen_product(bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]),
                                              clear.matrix_window(0, 0, nrows-split, prev_pivot_count),
                                              top.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                              mul_cutoff)
                bottom_start = non_pivots._intervals[0][0]
            # [  I      A'B   ]
            # [  0   D - CA'B ]

            # Now recursively do echelon form on the bottom
            bottom_pivots_rel = strassen_echelon_c(bottom.matrix_window(0, bottom_start, nrows-split, ncols-bottom_start), cutoff, mul_cutoff)
            # [  I  A'B ]
            # [  0  I F ]
            bottom_pivots = []
            for pivot in bottom_pivots_rel:
                bottom_pivots.append(pivot + bottom_start)
            bottom_h = len(bottom_pivots)

            if bottom_h == 0:
                pass
                # [  I  A'B ]
                # [  0   0  ]
                # proceed to sorting

            else:
                #     [  I  A'B ]  =  [  I  E  G  ]
                # let [  0  I F ]  =  [  0  I  F  ]
                top_cut = max(max(bottom_pivots) + 1, bottom_cut)

                # Note: left with respect to leftmost nonzero column of bottom
                top_left = top.matrix_window(0, bottom_start, top_h, top_cut - bottom_start)

                if bottom_h + top_h < ncols:

                    if top_cut - bottom_start == bottom_h:
                        clear = top_left
                    else:
                        clear = top_left.to_matrix().matrix_from_columns(bottom_pivots_rel).matrix_window()

                # subtract off E times bottom from top right
                if top_cut < ncols:

                    top_right = top.matrix_window(0, top_cut, top_h, ncols - top_cut)
                    subtract_strassen_product(top_right, clear, bottom.matrix_window(0, top_cut, bottom_h, ncols - top_cut), mul_cutoff)

                # [  I  *  G - EF ]
                # [  0  I     F   ]

                # Now subtract of E times bottom from top left
                if top_cut - bottom_start == bottom_h:
                    top_left.set_to_zero()

                else:
                    bottom_pivot_intervals = int_range(bottom_pivots)
                    non_pivots = int_range(bottom_start, top_cut - bottom_start) - bottom_pivot_intervals - top_pivot_intervals
                    for cols in non_pivots:
                        if cols[0] == 0: continue
                        prev_pivot_count = len(bottom_pivot_intervals - int_range(cols[0]+cols[1], top_cut - cols[0]+cols[1]))
                        subtract_strassen_product(top.matrix_window(0, cols[0], top_h, cols[1]),
                                                  clear.matrix_window(0, 0, top_h, prev_pivot_count),
                                                  bottom.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                                  mul_cutoff)
                    for cols in bottom_pivot_intervals:
                        top.matrix_window(0, cols[0], top_h, cols[1]).set_to_zero()

                # [  I  0  G - EF ]
                # [  0  I     F   ]
                # proceed to sorting

    # subrows already sorted...maybe I could do this more efficiently in cases with few pivot columns (e.g. merge sort)

    pivots = top_pivots
    pivots.extend(bottom_pivots)
    pivots.sort()

    cdef Py_ssize_t i, cur_row
    for cur_row from 0 <= cur_row < len(pivots):
        pivot = pivots[cur_row]
        for i from cur_row <= i < nrows:
            if not A.element_is_zero(i, pivot):
                break
        if i > cur_row and i < nrows:
            A.swap_rows(i, cur_row)

    return pivots


################################
# lots of room for optimization....
# eventually, should I just pass these around
# rather than lists of ints for pivots?
# would need new from_cols
class int_range:
    r"""
    Represent a list of integers as a list of integer intervals.

    .. NOTE::

        Repetitions are not considered.

    Useful class for dealing with pivots in the Strassen echelon, could
    have much more general application

    INPUT:

    It can be one of the following:

    - ``indices`` -- integer; start of the unique interval
    - ``range`` -- integer; length of the unique interval

    OR

    - ``indices`` -- list of integers, the integers to wrap into intervals

    OR

    - ``indices`` -- ``None`` (default), shortcut for an empty list

    OUTPUT:

    An instance of ``int_range``, i.e. a list of pairs ``(start, length)``.

    EXAMPLES:

    From a pair of integers::

        sage: from sage.matrix.strassen import int_range
        sage: int_range(2, 4)
        [(2, 4)]

    Default::

        sage: int_range()
        []

    From a list of integers::

        sage: int_range([1,2,3,4])
        [(1, 4)]
        sage: int_range([1,2,3,4,6,7,8])
        [(1, 4), (6, 3)]
        sage: int_range([1,2,3,4,100,101,102])
        [(1, 4), (100, 3)]
        sage: int_range([1,1000,2,101,3,4,100,102])
        [(1, 4), (100, 3), (1000, 1)]

    Repetitions are not considered::

        sage: int_range([1,2,3])
        [(1, 3)]
        sage: int_range([1,1,1,1,2,2,2,3])
        [(1, 3)]

    AUTHORS:

    - Robert Bradshaw
    """
    def __init__(self, indices=None, range=None):
        r"""
        See ``sage.matrix.strassen.int_range`` for full documentation.

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: int_range(2, 4)
            [(2, 4)]
        """
        if indices is None:
            self._intervals = []
            return
        if range is not None:
            self._intervals = [(int(indices), int(range))]
        else:
            self._intervals = []
            if not indices:
                return
            indices.sort()
            start = None
            last = None
            for ix in indices:
                if last is None:
                    start = ix
                elif ix-last > 1:
                    self._intervals.append((start, last-start+1))
                    start = ix
                last = ix
            self._intervals.append((start, last-start+1))

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: int_range([4,5,6,20,21,22,23])
            [(4, 3), (20, 4)]
            sage: int_range([])
            []
        """
        return str(self._intervals)

    def intervals(self):
        r"""
        Return the list of intervals.

        OUTPUT: list of pairs of integers

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([4,5,6,20,21,22,23])
            sage: I.intervals()
            [(4, 3), (20, 4)]
            sage: type(I.intervals())
            <... 'list'>
        """
        return self._intervals

    def to_list(self):
        r"""
        Return the (sorted) list of integers represented by this object.

        OUTPUT: list of integers

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([6,20,21,4,5,22,23])
            sage: I.to_list()
            [4, 5, 6, 20, 21, 22, 23]

        ::

            sage: I = int_range(34, 9)
            sage: I.to_list()
            [34, 35, 36, 37, 38, 39, 40, 41, 42]

        Repetitions are not considered::

            sage: I = int_range([1,1,1,1,2,2,2,3])
            sage: I.to_list()
            [1, 2, 3]
        """
        return [i for iv0, iv1 in self._intervals
                for i in range(iv0, iv0 + iv1)]

    def __iter__(self):
        r"""
        Return an iterator over the intervals.

        OUTPUT: iterator

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([6,20,21,4,5,22,23])
            sage: it = iter(I)
            sage: next(it)
            (4, 3)
            sage: next(it)
            (20, 4)
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        return iter(self._intervals)

    def __len__(self):
        r"""
        Return the number of integers represented by this object.

        OUTPUT: Python integer

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([6,20,21,4,5,22,23])
            sage: len(I)
            7

        ::

            sage: I = int_range([1,1,1,1,2,2,2,3])
            sage: len(I)
            3
        """
        return sum(iv1 for _, iv1 in self._intervals)

    def __add__(self, right):
        r"""
        Return the union of ``self`` and ``right``.

        INPUT:

        - ``right`` -- an instance of ``int_range``

        OUTPUT: an instance of ``int_range``

        .. NOTE::

            Yes, this two could be a lot faster...
            Basically, this class is for abstracting away what I was trying
            to do by hand in several places

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([1,1,1,1,2,2,2,3])
            sage: J = int_range([6,20,21,4,5,22,23])
            sage: I + J
            [(1, 6), (20, 4)]
        """
        full_list = self.to_list()
        full_list.extend(right.to_list())
        return int_range(full_list)

    def __sub__(self, right):
        r"""
        Return the set difference of ``self`` and ``right``.

        INPUT:

        - ``right`` -- an instance of ``int_range``

        OUTPUT: an instance of ``int_range``

        .. NOTE::

            Yes, this two could be a lot faster...
            Basically, this class is for abstracting away what I was trying
            to do by hand in several places

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([1,2,3,4,5])
            sage: J = int_range([6,20,21,4,5,22,23])
            sage: J - I
            [(6, 1), (20, 4)]
        """
        right_list = set(right.to_list())
        diff = [i for i in self.to_list() if i not in right_list]
        return int_range(diff)

    def __mul__(self, right):
        r"""
        Return the intersection of ``self`` and ``right``.

        INPUT:

        - ``right`` -- an instance of ``int_range``

        OUTPUT: an instance of ``int_range``

        EXAMPLES::

            sage: from sage.matrix.strassen import int_range
            sage: I = int_range([1,2,3,4,5])
            sage: J = int_range([6,20,21,4,5,22,23])
            sage: J * I
            [(4, 2)]
        """
        self_list = set(self.to_list())
        intersection = [i for i in right.to_list() if i in self_list]
        return int_range(intersection)


# Useful test code:
def test(n, m, R, c=2):
    r"""
    Test code for the Strassen algorithm.

    INPUT:

    - ``n`` -- integer
    - ``m`` -- integer
    - ``R`` -- ring
    - ``c`` -- integer (default: 2)

    EXAMPLES::

        sage: from sage.matrix.strassen import test
        sage: for n in range(5):
        ....:     print("{} {}".format(n, test(2*n,n,Frac(QQ['x']),2)))
        0 True
        1 True
        2 True
        3 True
        4 True
    """
    from sage.matrix.constructor import matrix
    A = matrix(R, n, m, range(n * m))
    B = A.__copy__()
    B._echelon_in_place_classical()
    C = A.__copy__()
    C._echelon_strassen(c)
    return B == C


# This stuff gets tested extensively elsewhere, and the functions
# below aren't callable now without using Pyrex.


#     todo: doc cutoff parameter as soon as I work out what it really means

#     EXAMPLES:

#         The following matrix dimensions are chosen especially to exercise the
#         eight possible parity combinations that could occur while subdividing
#         the matrix in the strassen recursion. The base case in both cases will
#         be a (4x5) matrix times a (5x6) matrix.

#         TODO -- the doctests below are currently not
#         tested/enabled/working -- enable them when linear algebra
#         restructing gets going.

#         sage: dim1 = 64; dim2 = 83; dim3 = 101
#         sage: R = MatrixSpace(QQ, dim1, dim2)
#         sage: S = MatrixSpace(QQ, dim2, dim3)
#         sage: T = MatrixSpace(QQ, dim1, dim3)


#         sage: A = R.random_element(range(-30, 30))
#         sage: B = S.random_element(range(-30, 30))
#         sage: C = T(0)
#         sage: D = T(0)

#         sage: A_window = A.matrix_window(0, 0, dim1, dim2)
#         sage: B_window = B.matrix_window(0, 0, dim2, dim3)
#         sage: C_window = C.matrix_window(0, 0, dim1, dim3)
#         sage: D_window = D.matrix_window(0, 0, dim1, dim3)

#         sage: from sage.matrix.strassen import strassen_window_multiply
#         sage: strassen_window_multiply(C_window, A_window, B_window, 2)   # use strassen method
#         sage: D_window.set_to_prod(A_window, B_window)             # use naive method
#         sage: C_window == D_window
#         True

#         sage: dim1 = 79; dim2 = 83; dim3 = 101
#         sage: R = MatrixSpace(QQ, dim1, dim2)
#         sage: S = MatrixSpace(QQ, dim2, dim3)
#         sage: T = MatrixSpace(QQ, dim1, dim3)

#         sage: A = R.random_element(range(30))
#         sage: B = S.random_element(range(30))
#         sage: C = T(0)
#         sage: D = T(0)

#         sage: A_window = A.matrix_window(0, 0, dim1, dim2)
#         sage: B_window = B.matrix_window(0, 0, dim2, dim3)
#         sage: C_window = C.matrix_window(0, 0, dim1, dim3)

#         sage: strassen_window_multiply(C, A, B, 2)   # use strassen method
#         sage: D.set_to_prod(A, B)                    # use naive method

#         sage: C == D
#         True
