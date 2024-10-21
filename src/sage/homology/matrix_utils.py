# sage_setup: distribution = sagemath-modules
"""
Utility Functions for Matrices

The actual computation of homology groups ends up being linear algebra
with the differentials thought of as matrices. This module contains
some utility functions for this purpose.
"""
########################################################################
#       Copyright (C) 2013 John H. Palmieri <palmieri@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
########################################################################

# TODO: this module is a clear candidate for cythonizing. Need to
# evaluate speed benefits.

from sage.matrix.constructor import matrix


def dhsw_snf(mat, verbose=False):
    """
    Preprocess a matrix using the "Elimination algorithm" described by
    Dumas et al. [DHSW2003]_, and then call ``elementary_divisors`` on the
    resulting (smaller) matrix.

    .. NOTE::

        'snf' stands for 'Smith Normal Form'.

    INPUT:

    - ``mat`` -- integer matrix, either sparse or dense

    (They use the transpose of the matrix considered here, so they use
    rows instead of columns.)

    ALGORITHM:

    Go through ``mat`` one column at a time.  For each
    column, add multiples of previous columns to it until either

    - it's zero, in which case it should be deleted.
    - its first nonzero entry is 1 or -1, in which case it should be kept.
    - its first nonzero entry is something else, in which case it is
      deferred until the second pass.

    Then do a second pass on the deferred columns.

    At this point, the columns with 1 or -1 in the first entry
    contribute to the rank of the matrix, and these can be counted and
    then deleted (after using the 1 or -1 entry to clear out its row).
    Suppose that there were `N` of these.

    The resulting matrix should be much smaller; we then feed it
    to Sage's ``elementary_divisors`` function, and prepend `N` 1s to
    account for the rows deleted in the previous step.

    EXAMPLES::

        sage: from sage.homology.matrix_utils import dhsw_snf
        sage: mat = matrix(ZZ, 3, 4, range(12))
        sage: dhsw_snf(mat)
        [1, 4, 0]
        sage: mat = random_matrix(ZZ, 20, 20, x=-1, y=2)
        sage: mat.elementary_divisors() == dhsw_snf(mat)
        True
    """
    ring = mat.base_ring()
    rows = mat.nrows()
    cols = mat.ncols()
    new_data = {}
    new_mat = matrix(ring, rows, cols, new_data)
    add_to_rank = 0
    zero_cols = 0
    if verbose:
        print("old matrix: %s by %s" % (rows, cols))
    # leading_positions: dictionary of lists indexed by row: if first
    # nonzero entry in column c is in row r, then leading_positions[r]
    # should contain c
    leading_positions = {}
    # pass 1:
    if verbose:
        print("starting pass 1")
    for j in range(cols):
        # new_col is a matrix with one column: sparse matrices seem to
        # be less buggy than sparse vectors (#5184, #5185), and
        # perhaps also faster.
        new_col = mat.matrix_from_columns([j])
        if new_col.is_zero():
            zero_cols += 1
        else:
            check_leading = True
            while check_leading:
                i = new_col.nonzero_positions_in_column(0)[0]
                entry = new_col[i, 0]
                check_leading = False
                if i in leading_positions:
                    for c in leading_positions[i]:
                        earlier = new_mat[i, c]
                        # right now we don't check to see if entry divides
                        # earlier, because we don't want to modify the
                        # earlier columns of the matrix.  Deal with this
                        # in pass 2.
                        if entry and earlier.divides(entry):
                            quo = entry.divide_knowing_divisible_by(earlier)
                            new_col = new_col - quo * new_mat.matrix_from_columns([c])
                            entry = 0
                            if not new_col.is_zero():
                                check_leading = True
            if not new_col.is_zero():
                new_mat.set_column(j-zero_cols, new_col.column(0))
                i = new_col.nonzero_positions_in_column(0)[0]
                if i in leading_positions:
                    leading_positions[i].append(j-zero_cols)
                else:
                    leading_positions[i] = [j-zero_cols]
            else:
                zero_cols += 1
    # pass 2:
    # first eliminate the zero columns at the end
    cols = cols - zero_cols
    zero_cols = 0
    new_mat = new_mat.matrix_from_columns(range(cols))
    if verbose:
        print("starting pass 2")
    keep_columns = list(range(cols))
    check_leading = True
    while check_leading:
        check_leading = False
        new_leading = leading_positions.copy()
        for i in leading_positions:
            if len(leading_positions[i]) > 1:
                j = leading_positions[i][0]
                jth = new_mat[i, j]
                for n in leading_positions[i][1:]:
                    nth = new_mat[i, n]
                    if jth.divides(nth):
                        quo = nth.divide_knowing_divisible_by(jth)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    elif nth.divides(jth):
                        quo = jth.divide_knowing_divisible_by(nth)
                        jth = nth
                        new_mat.swap_columns(n, j)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    else:
                        g, r, s = jth.xgcd(nth)
                        unit, A, B = r.xgcd(-s)  # unit ought to be 1 here
                        jth_col = new_mat.column(j)
                        nth_col = new_mat.column(n)
                        new_mat.set_column(j, r*jth_col + s*nth_col)
                        new_mat.set_column(n, B*jth_col + A*nth_col)
                        nth = B*jth + A*nth
                        jth = g
                        # at this point, jth should divide nth
                        quo = nth.divide_knowing_divisible_by(jth)
                        new_mat.add_multiple_of_column(n, j, -quo)
                    new_leading[i].remove(n)
                    if new_mat.column(n).is_zero():
                        keep_columns.remove(n)
                        zero_cols += 1
                    else:
                        new_r = new_mat.column(n).nonzero_positions()[0]
                        if new_r in new_leading:
                            new_leading[new_r].append(n)
                        else:
                            new_leading[new_r] = [n]
                        check_leading = True
        leading_positions = new_leading
    # pass 3: get rid of columns which start with 1 or -1
    if verbose:
        print("starting pass 3")
    max_leading = 1
    for i in leading_positions:
        j = leading_positions[i][0]
        entry = new_mat[i, j]
        if entry.abs() == 1:
            add_to_rank += 1
            keep_columns.remove(j)
            for c in new_mat.nonzero_positions_in_row(i):
                if c in keep_columns:
                    new_mat.add_multiple_of_column(c, j, -entry * new_mat[i, c])
        else:
            max_leading = max(max_leading, new_mat[i, j].abs())
    # form the new matrix
    if max_leading != 1:
        new_mat = new_mat.matrix_from_columns(keep_columns)
        if verbose:
            print("new matrix: %s by %s" % (new_mat.nrows(), new_mat.ncols()))
        if new_mat.is_sparse():
            ed = [1]*add_to_rank + new_mat.dense_matrix().elementary_divisors()
        else:
            ed = [1]*add_to_rank + new_mat.elementary_divisors()
    else:
        if verbose:
            print("new matrix: all pivots are 1 or -1")
        ed = [1]*add_to_rank

    if len(ed) < rows:
        return ed + [0]*(rows - len(ed))
    return ed[:rows]
