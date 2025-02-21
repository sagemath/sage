r"""
Covering Arrays (CA)

A Covering Array, denoted `CA(N;t,k,v)`, is an `N` by `k` array with
entries from a set of `v` elements with the property that in every
selection of `t` columns, every sequence of `t` elements appears in at
least one row.

An Orthogonal Array, denoted `OA(N;t,k,v)` is a covering array with the
property that every sequence of `t`-elements appears in exactly one row.
(See :mod:`sage.combinat.designs.orthogonal_arrays`).

This module collects methods relating to covering arrays, some of which
are inherited from orthogonal array methods. This module defines the
following functions:

.. csv-table::
    :class: contentstable
    :widths: 50, 50
    :delim: |

    :meth:`~sage.combinat.designs.designs_pyx.is_covering_array` | Check that an input list of lists is a `CA(N;t,k,v)`.
    :meth:`~sage.combinat.designs.covering_array.CA_relabel` | Return a relabelled version of the `CA`.
    :meth:`~sage.combinat.designs.covering_array.CA_standard_label` | Return a version of the `CA` relabelled to symbols `(0,\dots,n-1)`.
    :meth:`~sage.combinat.designs.covering_array.truncate_columns` | Return an array with `k` columns from a larger one.
    :meth:`~sage.combinat.designs.covering_array.Kleitman_Spencer_Katona` | Return a `CA(N; 2, k, 2)` using N as input.
    :meth:`~sage.combinat.designs.covering_array.column_Kleitman_Spencer_Katona` | Return a `CA(N; 2, k, 2)` using k as input.
    :meth:`~sage.combinat.designs.covering_array.database_check` | Check if CA can be made from the database of combinatorial designs.
    :meth:`~sage.combinat.designs.covering_array.covering_array` | Return a `CA` with given parameters.

REFERENCES:

- [Colb2004]_

- [SMC2006]_

- [WC2007]_

AUTHORS:

- Aaron Dwyer and brett stevens (2022): initial version
"""

# **********************************************************************
#       Copyright (C) 2022 Aaron Dwyer <aarondwyer@cmail.carleton.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************

from .orthogonal_arrays import OA_relabel, OA_standard_label
CA_relabel = OA_relabel
CA_standard_label = OA_standard_label


def truncate_columns(array, k):
    r"""
    Return a covering array with `k` columns, obtained by removing excess
    columns from a larger covering array.

    INPUT:

    - ``array`` -- the array to be truncated.

    - ``k`` -- the number of columns desired. Must be less than the
      number of columns in ``array``.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: from sage.combinat.designs.covering_array import truncate_columns
        sage: from sage.combinat.designs.database import ca_11_2_5_3
        sage: C = ca_11_2_5_3()
        sage: D = truncate_columns(C,7)
        Traceback (most recent call last):
        ...
        ValueError: array only has 5 columns
        sage: E = truncate_columns(C,4)
        sage: is_covering_array(E,parameters=True)
        (True, (11, 2, 4, 3))

    """
    oldk = len(array[0])

    if oldk == k:
        return array

    if oldk < k:
        raise ValueError("array only has {} columns".format(oldk))

    return [row[:k] for row in array]


def Kleitman_Spencer_Katona(N):
    r"""
    Return a `CA(N; 2, k, 2)` where `k = \binom {N-1}{\lceil N/2 \rceil}`.

    INPUT:

    - ``N`` -- the number of rows in the array, must be an integer greater
      than 3 since any smaller would not produce enough columns for a
      strength 2 array.

    This construction is referenced in [Colb2004]_ from [KS1973]_ and [Kat1973]_

    **Construction**

    Take all distinct binary `N`-tuples of weight `N/2` that have a 0
    in the first position and place them as columns in an array.

    EXAMPLES::

        sage: from sage.combinat.designs.covering_array import Kleitman_Spencer_Katona
        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: C = Kleitman_Spencer_Katona(2)
        Traceback (most recent call last):
        ...
        ValueError: N must be greater than 3
        sage: C = Kleitman_Spencer_Katona(5)
        sage: is_covering_array(C,parameters=True)
        (True, (5, 2, 4, 2))

    """
    from itertools import combinations
    from sage.arith.misc import integer_ceil
    if N < 4:
        raise ValueError("N must be greater than 3")

    col_list = []
    for p in combinations(range(N-1), integer_ceil(N/2)):
        S = [0]*N
        for i in p:
            S[i] = 1
        col_list.append(S)
    return [[col_list[j][i] for j in range(len(col_list))] for i in range(N)]


def column_Kleitman_Spencer_Katona(k):
    r"""
    Return a covering array with `k` columns using the Kleitman-Spencer-Katona
    method.

    See :func:`~sage.combinat.designs.covering_array.Kleitman_Spencer_Katona`

    INPUT:

    - ``k`` -- the number of columns in the array, must be an integer
      greater than 3 since any smaller is a trivial array for strength 2.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: from sage.combinat.designs.covering_array import column_Kleitman_Spencer_Katona
        sage: C = column_Kleitman_Spencer_Katona(20)
        sage: is_covering_array(C,parameters=True)
        (True, (8, 2, 20, 2))
        sage: column_Kleitman_Spencer_Katona(25000)
        Traceback (most recent call last):
        ...
        NotImplementedError: not implemented for k > 24310

    """
    kdict = [(3, 4), (4, 5), (10, 6), (15, 7), (35, 8), (56, 9),
             (126, 10), (210, 11), (462, 12), (792, 13), (1716, 14),
             (3003, 15), (6435, 16), (11440, 17), (24310, 18)]

    if k > kdict[-1][0]:
        raise NotImplementedError("not implemented for k > {}".format(kdict[-1][0]))

    for (ki, N) in kdict:
        if k <= ki:
            return truncate_columns(Kleitman_Spencer_Katona(N), k)


def database_check(number_columns, strength, levels):
    r"""
    Check if the database can be used to build a CA with the given parameters.
    If so return the CA, if not return False.

    INPUT:

    - ``strength`` (integer) -- the parameter `t` of the covering array,
      such that in any selection of `t` columns of the array, every
      `t`-tuple appears at least once.

    - ``levels`` (integer) -- the parameter `v` which is the number of
      unique symbols that appear in the covering array.

    - ``number_columns`` (integer) -- the number of columns desired for
      the covering array.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: from sage.combinat.designs.covering_array import database_check
        sage: C = database_check(6, 2, 3)
        sage: is_covering_array(C, parameters=True)
        (True, (12, 2, 6, 3))
        sage: database_check(6, 3, 3)
        False

    """
    import sage.combinat.designs.database as DB

    if (strength, levels) in DB.CA_constructions:
        for i in DB.CA_constructions[(strength, levels)]:
            if number_columns <= i[1]:
                CA = "ca_{}_{}_{}_{}".format(i[0], strength, i[1], levels)
                f = getattr(DB, CA)
                return truncate_columns(f(), number_columns)
        return False
    else:
        return False


def covering_array(strength, number_columns, levels):
    r"""
    Build a `CA(N; t, k, v)` using direct constructions, where `N` is the
    smallest size known.

    INPUT:

    - ``strength`` (integer) -- the parameter `t` of the covering array,
      such that in any selection of `t` columns of the array, every
      `t`-tuple appears at least once.

    - ``levels`` (integer) -- the parameter `v` which is the number of
      unique symbols that appear in the covering array.

    - ``number_columns`` (integer) -- the number of columns desired for
      the covering array.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: from sage.combinat.designs.covering_array import covering_array
        sage: C1 = covering_array(2, 7, 3)
        sage: is_covering_array(C1,parameters=True)
        (True, (12, 2, 7, 3))
        sage: C2 = covering_array(2, 11, 2)
        sage: is_covering_array(C2,parameters=True)
        (True, (7, 2, 11, 2))
        sage: C3 = covering_array(2, 8, 7)
        sage: is_covering_array(C3,parameters=True)
        (True, (49, 2, 8, 7))
        sage: C4 = covering_array(2, 50, 7)
        No direct construction known and/or implemented for a CA(N; 2, 50, 7)

    """
    from sage.combinat.designs.orthogonal_arrays import orthogonal_array

    if levels == 2 and strength == 2:
        return column_Kleitman_Spencer_Katona(number_columns)

    in_database = database_check(number_columns, strength, levels)
    if in_database:
        return in_database

    if orthogonal_array(number_columns, levels, strength, existence=True) is True:
        return orthogonal_array(number_columns, levels, strength)

    else:
        print("No direct construction known and/or implemented for a CA(N; {}, {}, {})".format(
            strength, number_columns, levels))
        return
