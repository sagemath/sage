r"""
Hash Families


Hash Families are a way to look at hash functions as an array of
integers. Both Perfect Hash Families and Covering Perfect Hash Families
are hash families with special properties that can be used to build
covering arrays.


A perfect hash family denoted `PHF(N;t,k,v)` is an array with
`N` rows and `k` columns with entries from a set of `v` elements with
the property that in every selection of `t` columns, there exists at
least one row where every symbol is unique.


Let `v` be a prime power, and `F_v` be the finite field order
`v`. Let `R_{t,v}` be the `v^t \times t` matrix consisting of length `t`
row vectors whose entries are from `F_v`. Let `V_{t,v}` be the
set of non-zero length `t` column vectors with entries from `F_v`
with first non-zero coordinate equal to 1. Finally let `U_{t,v}` be the
vectors of `V_{t,v}` whose first coordinate is 1.


A Covering Perfect Hash Family denoted `CPHF(N;t,k,v)` is an
array, C, with `N` rows and `k` columns with entries from `V_{t,v}`
with the property that for every selection of `t` column indices,
`\{ \gamma_1,...,\gamma_t\}` there is at least one row index, `\rho`,
for which the matrix `[C_{\rho\gamma_1},... C_{\rho\gamma_t}]` is non
singular. The array is a Sherwood Covering Perfect Hash Family (SCPHF)
if it uses entries from `U_{t,v}` instead of `V_{t,v}`.


This module collects methods relating to PHFs and CPHFs


.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |


    :meth:`~sage.combinat.designs.hash_family.is_PHF` | Check that an input list of lists is a `PHF(N;t,k,v)`.
    :meth:`~sage.combinat.designs.hash_family.strength2_PHF` | Create a `PHF(N;2,k,2)`.
    :meth:`~sage.combinat.designs.hash_family.is_CPHF` | Check that an input list of lists is a `CPHF(N;t,k,v)`.
    
   
REFERENCES:


- [WC2007]_


- [SMC2006]_


- [Colb2004]_


AUTHORS:


- Aaron Dwyer and brett stevens (2024): initial version
"""


# ********************************************************************** #
#       Copyright (C) 2024 Aaron Dwyer <aarondwyer@cmail.carleton.ca>    #
#                                                                        #
# This program is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 2 of the License, or      # 
# (at your option) any later version.                                    #
#                  https://www.gnu.org/licenses/                         #
# ********************************************************************** #


def is_PHF(array, strength):
    r"""
    Check if the input is a perfect hash family with given strength.


    INPUT:


    - ``array`` -- the perfect hash family to be tested, is a list of
      lists of positive integers including 0.


    - ``strength`` -- integer; the parameter `t` of the perfect hash
      family such that in any selection of `t` columns of the array,
      there is at least one row where all entries are unique


    EXAMPLES::


        sage: from sage.combinat.designs.hash_family import is_PHF
        sage: P = [[0, 0, 0, 0, 1, 1, 1, 1],
        ....:      [0, 0, 1, 1, 0, 0, 1, 1],  
        ....:      [0, 1, 0, 1, 0, 1, 0, 1]]
        sage: is_PHF(P,  3)
        Traceback (most recent call last):
        ...
        ValueError: "PHF must have v >= t, this array only has 2 symbols"


        sage: P = [[1, 1, 1, 0],
        ....:      [1, 1, 0, 0],
        ....:      [0, 0, 0]]
        sage: is_PHF(P,  2)
        Traceback (most recent call last):
        ...
        ValueError: Not all rows are the same length, row 2 is not the same length as row 0


        sage: P = [[0, 0, 0, 1, 1, 1, 2, 2, 2],
        ...:       [0, 1, 2, 0, 1, 2, 0, 1, 2]]
        sage: is_PHF(P,  2)
        True
    """
    from itertools import combinations


    symbol_list = list({x for l in array for x in l})


    if len(symbol_list) < strength:
        raise ValueError("PHF must have v >= t, this array only has {} symbols".format(len(symbol_list)))


    k = len(array[0])
    for row in array:
        if len(row) != k:
            raise ValueError("Not all rows are the same length, row {} is not the same length as row 0.format(row)")
           
    for comb in combinations(range(len(array[0])), strength):
        found = False
        for row in array:
            # Check if the row has all entries distinct, if we cannot
            # find one then this is not a PHF
            rowlist = [row[ti] for ti in comb]
            if len(set(rowlist)) == len(rowlist):
                found = True
                break
        if found == False:
            return False
    return True


def strength2_PHF(num_rows, num_symbols):
    r"""
    Return a `PHF(N;2,v^N,v)` by placing all `v^N` tuples of length
    `N` as columns.


    INPUT:


    - ``num_rows`` -- positive integer; the parameter `N` of the PHF,
      which is the number of rows


    - ``num_symbols`` -- positive integer; the parameter `v` of the PHF,
      which is the number of symbols


    EXAMPLES::
        sage: from sage.combinat.designs.hash_family import strength2_PHF
        sage: from sage.combinat.designs.hash_family import is_PHF
        sage: P = strength_2_PHF(3,3)
        sage: is_PHF(P,2)
        True
    """
    from itertools import product


    col_list = list(product(range(num_symbols), repeat = num_rows))
    return [[col_list[j][i] for j in range(len(col_list))] for i in range(num_rows)]


def is_CPHF(array, strength, field, compact = False):
    r"""
    Check if the input is a covering perfect hash family with given
    strength and entries are vectors with entries in the given field.


    INPUT:


    - ``array`` -- the covering perfect hash family to be tested, is a
      list of lists of lists which represents an array where each entry
      is a vector. These vectors must have entries that lie in `field`.


    - ``strength`` -- integer; the parameter `t` of the covering perfect
      hash family such that in any selection of `t` columns of the array,
      there is at least one row where the submatrix formed by the `t`
      vectors in that row is non singular.


    - ``field`` -- a sage.rings.finite_rings.finite_field object.


    EXAMPLES::
        sage: from sage.combinat.designs.designs_pyx import is_CPHF
        sage: C = [[[0,0,0],[1,1,1],[0,0,0]],
        ....:      [[1,0,1],[1,0,1],[1,1]]]
        ValueError("All vectors must be the same length")
    """
    from itertools import combinations
    from sage.matrix.constructor import Matrix


    # Ensure all vectors same length and all entries lie in `field`
    dim = len(array[0][0])
    for row in array:
        for vector in row:
            if len(vector) != dim:
                raise ValueError("All vectors must be the same length")
            for i in vector:
                if i.parent() != field:
                    print(i.parent())
                    raise ValueError("All entries must be in 'field'")


    for comb in combinations(range(len(array[0])), strength):
        # Iterate over every possible selections of t columns, and check
        # that in some row of the array, the t columns form a
        # non-singular matrix over the given field
        for row in array:
            found = False
            M = Matrix(field, [row[ti] for ti in comb])
            if M.rank() == strength:
                found = True
                break
        if not found:
            return False
    return True