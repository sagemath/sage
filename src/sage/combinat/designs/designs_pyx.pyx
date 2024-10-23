r"""
Cython functions for combinatorial designs

This module implements the design methods that need to be somewhat efficient.

Functions
---------
"""

from sage.data_structures.bitset_base cimport *

from libc.string cimport memset

from cysignals.memory cimport sig_malloc, sig_calloc, sig_realloc, sig_free

from sage.misc.unknown import Unknown


def is_covering_array(array, strength=None, levels=None, verbose=False, parameters=False):
    r"""
    Check if the input is a covering array with given strength.

    See :mod:`sage.combinat.designs.covering_array` for a definition.

    INPUT:

    - ``array`` -- the Covering Array to be tested

    - ``strength`` -- integer; the parameter `t` of the covering array,
      such that in any selection of `t` columns of the array, every `t`
      -tuple appears at least once. If set to None then all t > 0 are
      tested to and the maximal strength is used.

    - ``levels`` -- the number of symbols that appear in ``array``
      If set to None, then each unique entry in ``array`` is counted

    - ``verbose`` -- boolean; whether to display some information about
      the covering array

    - ``parameters`` -- boolean; whether to return the parameters of
      the Covering Array. If set to ``True``, the function returns a
      pair ``(boolean_answer,(N,t,k,v))``.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_covering_array
        sage: C = [[1, 1, 1, 0],
        ....:      [1, 1, 0, 0],
        ....:      [0, 0, 0]]
        sage: is_covering_array(C)
        Traceback (most recent call last):
        ...
        ValueError: Not all rows are the same length, row 2 is not the same length as row 0

        sage: C = [[0, 1, 1],
        ....:      [1, 1, 0],
        ....:      [1, 0, 1],
        ....:      [0, 0, 0,]]
        sage: is_covering_array(C,strength=4)
        Traceback (most recent call last):
        ...
        ValueError: Strength must be equal or less than number of columns

        sage: C = [[0, 1, 1],
        ....:      [1, 1, 1],
        ....:      [1, 0, 1]]
        sage: is_covering_array(C,verbose=True)
        A 3 by 3 Covering Array with strength 0 with entries from a symbol set of size 2
        True

        sage: C = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ....:      [0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        ....:      [0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
        ....:      [1, 0, 1, 1, 0, 1, 1, 0, 0, 1],
        ....:      [1, 1, 0, 1, 1, 0, 1, 0, 1, 0],
        ....:      [1, 1, 1, 0, 1, 1, 0, 1, 2, 0]]
        sage: is_covering_array(C,levels=2)
        Traceback (most recent call last):
        ...
        ValueError: Array should contain integer symbols from 0 to 1

        sage: C = [[1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2],
        ....:      [1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2],
        ....:      [1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0],
        ....:      [0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1],
        ....:      [2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2],
        ....:      [1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2],
        ....:      [1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1],
        ....:      [2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2],
        ....:      [1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0],
        ....:      [0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2],
        ....:      [1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0],
        ....:      [0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0],
        ....:      [0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1],
        ....:      [2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1],
        ....:      [2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1],
        ....:      [2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0],
        ....:      [0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2],
        ....:      [1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1],
        ....:      [2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1],
        ....:      [2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2],
        ....:      [1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1],
        ....:      [2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0],
        ....:      [0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1],
        ....:      [2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0],
        ....:      [0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0],
        ....:      [0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2],
        ....:      [1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2],
        ....:      [1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0],
        ....:      [0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0],
        ....:      [0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1],
        ....:      [2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0],
        ....:      [0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1],
        ....:      [2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2],
        ....:      [1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1],
        ....:      [2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1],
        ....:      [2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2],
        ....:      [1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0],
        ....:      [0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1],
        ....:      [2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1],
        ....:      [2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1],
        ....:      [2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0],
        ....:      [0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0],
        ....:      [0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2],
        ....:      [1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0],
        ....:      [0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2],
        ....:      [1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1],
        ....:      [2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2],
        ....:      [1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2],
        ....:      [1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1],
        ....:      [2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0],
        ....:      [0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2],
        ....:      [1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2],
        ....:      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        sage: is_covering_array(C,parameters=True)
        (True, (53, 3, 13, 3))

        sage: C = [[1, 0, 1, 1, 2, 0, 2, 2],
        ....:      [2, 1, 0, 1, 1, 2, 0, 2],
        ....:      [2, 2, 1, 0, 1, 1, 2, 0],
        ....:      [0, 2, 2, 1, 0, 1, 1, 2],
        ....:      [2, 0, 2, 2, 1, 0, 1, 1],
        ....:      [1, 2, 0, 2, 2, 1, 0, 1],
        ....:      [1, 1, 2, 0, 2, 2, 1, 0],
        ....:      [0, 1, 1, 2, 0, 2, 2, 1]]
        sage: is_covering_array(C,strength=2,parameters=True)
        (False, (8, 0, 8, 3))
    """
    from itertools import product, combinations

    if levels is None:
        symbol_list = list({x for l in array for x in l})
        levels = len(symbol_list)
    else:
        symbol_list = range(levels)

    number_rows = len(array)
    number_columns = len(array[0])

    for row in array:
        if len(row) != number_columns:
            raise ValueError("Not all rows are the same length, row {} is not the same length as row 0".format(array.index(row)))
        else:
            for entry in row:
                if int(entry) != entry or entry < -1 or entry >= levels:
                    raise ValueError("Array should contain integer symbols from 0 to {}".format(levels-1))

    result = True

    # If strength t is inputted, check that for every selection of t
    # columns, each v^t t-tuple is found in some row.
    if strength:
        if strength > number_columns:
            raise ValueError("Strength must be equal or less than number of columns")
        wstrength = strength
        for comb in combinations(range(number_columns), wstrength):
            existing_ttuples = set(tuple([row[ti] for ti in comb]) for row in array)
            if len(existing_ttuples) != levels ** wstrength:
                wstrength = 0
                result = False
                break

    # If no strength t is inputted, starting at t=1 check all t until
    # one of the v^t t-tuples does not appear.
    else:
        wstrength = 1
        finished = False
        do_iterate = True
        while finished is False:
            for comb in combinations(range(number_columns), wstrength):
                tuple_dictionary = {item: 0 for item in product(symbol_list, repeat=wstrength)}
                for row in array:
                    tuple_dictionary[tuple([row[ti] for ti in comb])] += 1
                if 0 in tuple_dictionary.values():
                    wstrength -= 1
                    finished = True
                    break
                elif do_iterate and any(value < levels for value in tuple_dictionary.values()):
                    do_iterate = False
                    finished = True
            if finished is False and wstrength < number_columns and do_iterate:
                wstrength += 1

    if verbose:
            print('A {} by {} Covering Array with strength {} with entries from a symbol set of size {}'.format(number_rows, number_columns, wstrength, levels))

    if parameters:
        return (result, (number_rows, wstrength, number_columns, levels))
    else:
        return result


def is_orthogonal_array(OA, int k, int n, int t=2, verbose=False, terminology='OA'):
    r"""
    Check that the integer matrix `OA` is an `OA(k,n,t)`.

    See :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    for a definition.

    INPUT:

    - ``OA`` -- the Orthogonal Array to be tested

    - ``k``, ``n``, ``t`` -- integers; only implemented for `t=2`

    - ``verbose`` -- boolean; whether to display some information when ``OA``
      is not an orthogonal array `OA(k,n)`

    - ``terminology`` -- string; how to phrase the information when ``verbose =
      True``. Possible values are `"OA"`, `"MOLS"`

    EXAMPLES::

        sage: # needs sage.schemes
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: OA = designs.orthogonal_arrays.build(8,9)
        sage: is_orthogonal_array(OA,8,9)
        True
        sage: is_orthogonal_array(OA,8,10)
        False
        sage: OA[4][3] = 1
        sage: is_orthogonal_array(OA,8,9)
        False
        sage: is_orthogonal_array(OA,8,9,verbose=True)
        Columns 0 and 3 are not orthogonal
        False
        sage: is_orthogonal_array(OA,8,9, verbose=True, terminology='MOLS')
        Squares 0 and 3 are not orthogonal
        False

    TESTS::

        sage: # needs sage.schemes
        sage: is_orthogonal_array(OA,8,9, t=3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for t=2
        sage: is_orthogonal_array([[3]*8],8,9, verbose=True)
        The number of rows is 1 instead of 9^2=81
        False
        sage: is_orthogonal_array([[3]*8],8,9, verbose=True, terminology='MOLS')
        All squares do not have dimension n^2=9^2
        False
        sage: is_orthogonal_array([[3]*7],8,9, verbose=True)
        Some row does not have length 8
        False
        sage: is_orthogonal_array([[3]*7],8,9, verbose=True, terminology='MOLS')
        The number of squares is not 6
        False

    Up to relabelling, there is a unique `OA(3,2)`. So their number is just the
    cardinality of the relabeling group which is `S_2^3 \times S_3` and has
    cardinality `48`::

        sage: from itertools import product
        sage: n = 0
        sage: for a in product(product((0,1), repeat=3), repeat=4):                     # needs sage.schemes
        ....:     if is_orthogonal_array(a,3,2):
        ....:          n += 1
        sage: n
        48
    """
    cdef int n2 = n*n
    cdef int x

    if t != 2:
        raise NotImplementedError("only implemented for t=2")

    for R in OA:
        if len(R) != k:
            if verbose:
                print({"OA": "Some row does not have length "+str(k),
                       "MOLS": "The number of squares is not "+str(k-2)}[terminology])
            return False

    if len(OA) != n2:
        if verbose:
            print({"OA": "The number of rows is {} instead of {}^2={}".format(len(OA),n,n2),
                   "MOLS": "All squares do not have dimension n^2={}^2".format(n)}[terminology])
        return False

    if n == 0:
        return True

    cdef int i,j,l

    # A copy of OA
    cdef unsigned short * OAc = <unsigned short *> sig_malloc(k*n2*sizeof(unsigned short))

    cdef unsigned short * C1
    cdef unsigned short * C2

    # failed malloc ?
    if OAc is NULL:
        raise MemoryError

    # Filling OAc
    for i,R in enumerate(OA):
        for j,x in enumerate(R):
            if x < 0 or x >= n:
                if verbose:
                    print({"OA": "{} is not in the interval [0..{}]".format(x,n-1),
                           "MOLS": "Entry {} was expected to be in the interval [0..{}]".format(x,n-1)}[terminology])
                sig_free(OAc)
                return False
            OAc[j*n2+i] = x

    # A bitset to keep track of pairs of values
    cdef bitset_t seen
    bitset_init(seen, n2)

    for i in range(k): # For any column C1
        C1 = OAc+i*n2
        for j in range(i+1,k): # For any column C2 > C1
            C2 = OAc+j*n2
            bitset_set_first_n(seen, 0) # No pair has been seen yet
            for l in range(n2):
                bitset_add(seen,n*C1[l]+C2[l])

            if bitset_len(seen) != n2: # Have we seen all pairs ?
                sig_free(OAc)
                bitset_free(seen)
                if verbose:
                    print({"OA": "Columns {} and {} are not orthogonal".format(i,j),
                           "MOLS": "Squares {} and {} are not orthogonal".format(i,j)}[terminology])
                return False

    sig_free(OAc)
    bitset_free(seen)
    return True


def is_group_divisible_design(groups,blocks,v,G=None,K=None,lambd=1,verbose=False):
    r"""
    Check that input is a Group Divisible Design on `\{0, \ldots, v-1\}`.

    For more information on Group Divisible Designs, see
    :class:`~sage.combinat.designs.group_divisible_designs.GroupDivisibleDesign`.

    INPUT:

    - ``groups`` -- a partition of `X`. If set to ``None`` the groups are
      guessed automatically, and the function returns ``(True, guessed_groups)``
      instead of ``True``

    - ``blocks`` -- collection of blocks

    - ``v`` -- integers; size of the ground set assumed to be `X=\{0,...,v-1\}`

    - ``G`` -- list of integers of which the sizes of the groups must be
      elements. Set to ``None`` (automatic guess) by default

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements. Set to ``None`` (automatic guess) by default

    - ``lambd`` -- value of `\lambda`. Set to `1` by default

    - ``verbose`` -- boolean; whether to display some information when the
      design is not a GDD

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_group_divisible_design
        sage: TD = designs.transversal_design(4,10)                                     # needs sage.modules
        sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
        sage: is_group_divisible_design(groups,TD,40,lambd=1)                           # needs sage.modules
        True

    TESTS::

        sage: TD = designs.transversal_design(4,10)                                     # needs sage.modules
        sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
        sage: is_group_divisible_design(groups, TD, 40, lambd=2, verbose=True)          # needs sage.modules
        the pair (0,10) has been seen 1 times but lambda=2
        False
        sage: is_group_divisible_design([[1,2],[3,4]],[[1,2]],40,lambd=1,verbose=True)
        groups is not a partition of [0,...,39]
        False
        sage: is_group_divisible_design([list(range(40))],[[1,2]],40,lambd=1,verbose=True)
        the pair (1,2) belongs to a group but appears in some block
        False
        sage: is_group_divisible_design([list(range(40))],[[2,2]],40,lambd=1,verbose=True)
        The following block has repeated elements: [2, 2]
        False
        sage: is_group_divisible_design([list(range(40))],[["e",2]],40,lambd=1,verbose=True)
        e does not belong to [0,...,39]
        False
        sage: is_group_divisible_design([list(range(40))],[list(range(40))],40,G=[5],lambd=1,verbose=True)
        a group has size 40 while G=[5]
        False
        sage: is_group_divisible_design([list(range(40))],[["e",2]],40,K=[1],lambd=1,verbose=True)
        a block has size 2 while K=[1]
        False

        sage: # needs sage.schemes
        sage: p = designs.projective_plane(3)
        sage: is_group_divisible_design(None, p.blocks(), 13)
        (True, [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])
        sage: is_group_divisible_design(None, p.blocks()*2, 13, verbose=True)
        the pair (0,1) has been seen 2 times but lambda=1
        False
        sage: is_group_divisible_design(None, p.blocks()*2, 13, lambd=2)
        (True, [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])
    """
    cdef int n = v
    cdef int i,ii,j,jj,s
    cdef int l = lambd
    cdef bint guess_groups = groups is None

    if v < 0 or lambd < 0:
        if verbose:
            print("v={} and lambda={} must be nonnegative integers".format(v,l))
        return False

    # Block sizes are element of K
    if K is not None:
        K = set(K)
        for b in blocks:
            if not len(b) in K:
                if verbose:
                    print("a block has size {} while K={}".format(len(b),list(K)))
                return False

    # Check that "groups" consists of disjoint sets whose union has length n
    if (groups is not None and
        (sum(len(g) for g in groups) != n or
         len(set().union(*groups)) != n)):
        if verbose:
            print("groups is not a partition of [0,...,{}]".format(n-1))
        return False

    # Checks that the blocks are indeed sets and do not repeat elements
    for b in blocks:
        if len(b) != len(set(b)):
            if verbose:
                print("The following block has repeated elements: {}".format(b))
            return False

    # Check that the groups/blocks belong to [0,...,n-1]
    from itertools import chain
    for b in chain(groups if groups is not None else [],blocks):
        for x in b:
            try:
                i = x
            except TypeError:
                i = -1
            if i < 0 or i >= n:
                if verbose:
                    print("{} does not belong to [0,...,{}]".format(x, n-1))
                return False

    cdef unsigned short * matrix = <unsigned short *> sig_calloc(n*n,sizeof(unsigned short))
    if matrix is NULL:
        raise MemoryError

    # Counts the number of occurrences of each pair of points
    for b in blocks:
        s = len(b)
        for i in range(s):
            ii = b[i]
            for j in range(i+1,s):
                jj = b[j]
                matrix[ii*n+jj] += 1
                matrix[jj*n+ii] += 1

    # Guess the groups (if necessary)
    if groups is None:
        from sage.sets.disjoint_set import DisjointSet_of_integers
        groups = DisjointSet_of_integers(n)
        for i in range(n):
            for j in range(i + 1, n):
                if matrix[i * n + j] == 0:
                    groups.union(i, j)
        groups = list(groups.root_to_elements_dict().values())

    # Group sizes are element of G
    if G is not None:
        G = set(G)
        for g in groups:
            if not len(g) in G:
                if verbose:
                    print("a group has size {} while G={}".format(len(g),list(G)))
                sig_free(matrix)
                return False

    # Checks that two points of the same group were never covered
    for g in groups:
        s = len(g)
        for i in range(s):
            ii = g[i]
            for j in range(i+1,s):
                jj = g[j]
                if matrix[ii*n+jj] != 0:
                    if verbose:
                        print("the pair ({},{}) belongs to a group but appears in some block".format(ii, jj))
                    sig_free(matrix)
                    return False

                # We fill the entries with what is expected by the next loop
                matrix[ii*n+jj] = l
                matrix[jj*n+ii] = l

    # Checking that what should be equal to lambda IS equal to lambda
    for i in range(n):
        for j in range(i+1,n):
            if matrix[i*n+j] != l:
                if verbose:
                    print("the pair ({},{}) has been seen {} times but lambda={}".format(i,j,matrix[i*n+j],l))
                sig_free(matrix)
                return False

    sig_free(matrix)

    return True if not guess_groups else (True, groups)


def is_pairwise_balanced_design(blocks,v,K=None,lambd=1,verbose=False):
    r"""
    Check that input is a Pairwise Balanced Design (PBD) on `\{0, \ldots, v-1\}`.

    For more information on Pairwise Balanced Designs (PBD), see
    :class:`~sage.combinat.designs.bibd.PairwiseBalancedDesign`.

    INPUT:

    - ``blocks`` -- collection of blocks

    - ``v`` -- integers; size of the ground set assumed to be `X=\{0,...,v-1\}`

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements; set to ``None`` (automatic guess) by default

    - ``lambd`` -- value of `\lambda` (default: `1`)

    - ``verbose`` -- boolean; whether to display some information when the
      design is not a PBD

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_pairwise_balanced_design
        sage: sts = designs.steiner_triple_system(9)
        sage: is_pairwise_balanced_design(sts,9,[3],1)
        True
        sage: TD = designs.transversal_design(4,10).blocks()                            # needs sage.modules
        sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
        sage: is_pairwise_balanced_design(TD + groups, 40, [4,10], 1, verbose=True)     # needs sage.modules
        True

    TESTS::

        sage: from sage.combinat.designs.designs_pyx import is_pairwise_balanced_design
        sage: is_pairwise_balanced_design(TD + groups, 40, [4,10], 2, verbose=True)     # needs sage.modules
        the pair (0,1) has been seen 1 times but lambda=2
        False
        sage: is_pairwise_balanced_design(TD + groups, 40, [10], 1, verbose=True)       # needs sage.modules
        a block has size 4 while K=[10]
        False
        sage: is_pairwise_balanced_design([[2,2]], 40, [2], 1, verbose=True)
        The following block has repeated elements: [2, 2]
        False
        sage: is_pairwise_balanced_design([["e",2]], 40, [2], 1, verbose=True)
        e does not belong to [0,...,39]
        False
    """
    return is_group_divisible_design([[i] for i in range(v)],
                                     blocks,
                                     v,
                                     K=K,
                                     lambd=lambd,
                                     verbose=verbose)


def is_projective_plane(blocks, verbose=False):
    r"""
    Test whether the blocks form a projective plane on `\{0,...,v-1\}`.

    A *projective plane* is an incidence structure that has the following properties:

    1. Given any two distinct points, there is exactly one line incident with both of them.
    2. Given any two distinct lines, there is exactly one point incident with both of them.
    3. There are four points such that no line is incident with more than two of them.

    For more informations, see :wikipedia:`Projective_plane`.

    :meth:`~IncidenceStructure.is_t_design` can also check if an incidence structure is a projective plane
    with the parameters `v=k^2+k+1`, `t=2` and `l=1`.

    INPUT:

    - ``blocks`` -- collection of blocks

    - ``verbose`` -- whether to print additional information

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_projective_plane
        sage: p = designs.projective_plane(4)                                           # needs sage.schemes
        sage: b = p.blocks()                                                            # needs sage.schemes
        sage: is_projective_plane(b, verbose=True)                                      # needs sage.schemes
        True

        sage: # needs sage.schemes
        sage: p = designs.projective_plane(2)
        sage: b = p.blocks()
        sage: is_projective_plane(b)
        True
        sage: b[0][2] = 5
        sage: is_projective_plane(b, verbose=True)
        the pair (0,5) has been seen 2 times but lambda=1
        False

        sage: is_projective_plane([[0,1,2],[1,2,4]], verbose=True)
        the pair (0,3) has been seen 0 times but lambda=1
        False

        sage: is_projective_plane([[1]], verbose=True)
        First block has less than 3 points.
        False

        sage: # needs sage.schemes
        sage: p = designs.projective_plane(2)
        sage: b = p.blocks()
        sage: b[2].append(4)
        sage: is_projective_plane(b, verbose=True)
        a block has size 4 while K=[3]
        False
    """
    if not blocks:
        if verbose:
            print('There is no block.')
        return False
    k = len(blocks[0])-1
    if k < 2:
        if verbose:
            print('First block has less than 3 points.')
        return False
    v = k**2 + k + 1
    return is_group_divisible_design([[i] for i in range(v)],
                                     blocks,
                                     v,
                                     K=[k+1],
                                     lambd=1,
                                     verbose=verbose)


def is_difference_matrix(M,G,k,lmbda=1,verbose=False):
    r"""
    Test if `M` is a `(G,k,\lambda)`-difference matrix.

    A matrix `M` is a `(G,k,\lambda)`-difference matrix if its entries are
    element of `G`, and if for any two rows `R,R'` of `M` and `x\in G` there
    are exactly `\lambda` values `i` such that `R_i-R'_i=x`.

    INPUT:

    - ``M`` -- a matrix with entries from ``G``

    - ``G`` -- a group

    - ``k`` -- integer

    - ``lmbda`` -- integer (default: `1`)

    - ``verbose`` -- boolean; whether to print some information when the answer
      is ``False``

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_difference_matrix
        sage: q = 3**3
        sage: F = GF(q,'x')                                                             # needs sage.rings.finite_rings
        sage: M = [[x*y for y in F] for x in F]                                         # needs sage.rings.finite_rings
        sage: is_difference_matrix(M,F,q,verbose=1)                                     # needs sage.rings.finite_rings
        True

        sage: B = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ....:      [0, 1, 2, 3, 4, 2, 3, 4, 0, 1],
        ....:      [0, 2, 4, 1, 3, 3, 0, 2, 4, 1]]
        sage: G = GF(5)
        sage: B = [[G(b) for b in R] for R in B]
        sage: is_difference_matrix(list(zip(*B)),G,3,2)
        True

    Bad input::

        sage: # needs sage.rings.finite_rings
        sage: for R in M: R.append(None)
        sage: is_difference_matrix(M,F,q,verbose=1)
        The matrix has 28 columns but k=27
        False
        sage: for R in M: _=R.pop(-1)
        sage: M.append([None]*3**3)
        sage: is_difference_matrix(M,F,q,verbose=1)
        The matrix has 28 rows instead of lambda(|G|-1+2u)+mu=1(27-1+2.0)+1=27
        False
        sage: _= M.pop(-1)
        sage: for R in M: R[-1] = 0
        sage: is_difference_matrix(M,F,q,verbose=1)
        Columns 0 and 26 generate 0 exactly 27 times instead of the expected mu(=1)
        False
        sage: for R in M: R[-1] = 1
        sage: M[-1][-1] = 0
        sage: is_difference_matrix(M,F,q,verbose=1)
        Columns 0 and 26 do not generate all elements of G exactly lambda(=1) times.
        The element x appeared 0 times as a difference.
        False
    """
    return is_quasi_difference_matrix(M,G,k,lmbda=lmbda,mu=lmbda,u=0,verbose=verbose)


def is_quasi_difference_matrix(M,G,int k,int lmbda,int mu,int u,verbose=False):
    r"""
    Test if the matrix is a `(G,k;\lambda,\mu;u)`-quasi-difference matrix.

    Let `G` be an abelian group of order `n`. A
    `(n,k;\lambda,\mu;u)`-quasi-difference matrix (QDM) is a matrix `Q_{ij}`
    with `\lambda(n-1+2u)+\mu` rows and `k` columns, with each entry either
    equal to ``None`` (i.e. the 'missing entries') or to an element of `G`. Each
    column contains exactly `\lambda u` empty entries, and each row contains at
    most one ``None``. Furthermore, for each `1\leq i<j\leq k`, the multiset

    .. MATH::

        \{q_{li}-q_{lj}:1\leq l\leq \lambda (n-1+2u)+\mu, \text{ with } q_{li}\text{ and }q_{lj}\text{ not empty}\}

    contains `\lambda` times every nonzero element of `G` and contains `\mu`
    times `0`.

    INPUT:

    - ``M`` -- a matrix with entries from ``G`` (or equal to ``None`` for
      missing entries)

    - ``G`` -- a group

    - ``k``, ``lmbda``, ``mu``, ``u`` -- integers

    - ``verbose`` -- boolean; whether to print some information when the answer
      is ``False``

    EXAMPLES:

    Differences matrices::

        sage: from sage.combinat.designs.designs_pyx import is_quasi_difference_matrix
        sage: q = 3**3
        sage: F = GF(q,'x')                                                             # needs sage.rings.finite_rings
        sage: M = [[x*y for y in F] for x in F]                                         # needs sage.rings.finite_rings
        sage: is_quasi_difference_matrix(M,F,q,1,1,0,verbose=1)                         # needs sage.rings.finite_rings
        True

        sage: B = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ....:      [0, 1, 2, 3, 4, 2, 3, 4, 0, 1],
        ....:      [0, 2, 4, 1, 3, 3, 0, 2, 4, 1]]
        sage: G = GF(5)
        sage: B = [[G(b) for b in R] for R in B]
        sage: is_quasi_difference_matrix(list(zip(*B)),G,3,2,2,0)
        True

    A quasi-difference matrix from the database::

        sage: from sage.combinat.designs.database import QDM
        sage: G,M = QDM[38,1][37,1,1,1][1]()
        sage: is_quasi_difference_matrix(M,G,k=6,lmbda=1,mu=1,u=1)
        True

    Bad input::

        sage: is_quasi_difference_matrix(M,G,k=6,lmbda=1,mu=1,u=3,verbose=1)
        The matrix has 39 rows instead of lambda(|G|-1+2u)+mu=1(37-1+2.3)+1=43
        False
        sage: is_quasi_difference_matrix(M,G,k=6,lmbda=1,mu=2,u=1,verbose=1)
        The matrix has 39 rows instead of lambda(|G|-1+2u)+mu=1(37-1+2.1)+2=40
        False
        sage: M[3][1] = None
        sage: is_quasi_difference_matrix(M,G,k=6,lmbda=1,mu=1,u=1,verbose=1)
        Row 3 contains more than one empty entry
        False
        sage: M[3][1] = 1
        sage: M[6][1] = None
        sage: is_quasi_difference_matrix(M,G,k=6,lmbda=1,mu=1,u=1,verbose=1)
        Column 1 contains 2 empty entries instead of the expected lambda.u=1.1=1
        False
    """
    from sage.combinat.designs.difference_family import group_law

    assert k>=2
    assert lmbda >=1
    assert mu>=0
    assert u>=0

    cdef int n = G.cardinality()
    cdef int M_nrows = len(M)
    cdef int i,j,ii
    cdef bint bit

    # Height of the matrix
    if lmbda*(n-1+2*u)+mu != M_nrows:
        if verbose:
            print("The matrix has {} rows instead of lambda(|G|-1+2u)+mu={}({}-1+2.{})+{}={}".format(M_nrows,lmbda,n,u,mu,lmbda*(n-1+2*u)+mu))
        return False

    # Width of the matrix
    for R in M:
        if len(R)!=k:
            if verbose:
                print("The matrix has {} columns but k={}".format(len(R),k))
            return False

    # When |G|=0
    if M_nrows == 0:
        return True

    # Map group element with integers
    cdef list int_to_group = list(G)
    cdef dict group_to_int = {v:i for i,v in enumerate(int_to_group)}

    # Allocations
    cdef int ** x_minus_y     = <int **> sig_malloc((n+1)*sizeof(int *))
    cdef int * x_minus_y_data = <int *>  sig_malloc((n+1)*(n+1)*sizeof(int))
    cdef int * M_c            = <int *>  sig_malloc(k*M_nrows*sizeof(int))
    cdef int * G_seen         = <int *>  sig_malloc((n+1)*sizeof(int))
    if (x_minus_y == NULL or x_minus_y_data == NULL or M_c == NULL or G_seen == NULL):
        sig_free(x_minus_y)
        sig_free(x_minus_y_data)
        sig_free(G_seen)
        sig_free(M_c)
        raise MemoryError

    # The "x-y" table. If g_i, g_j \in G, then x_minus_y[i][j] is equal to
    # group_to_int[g_i-g_j].
    #
    # In order to handle empty values represented by n, we have
    # x_minus_y[?][n]=x_minus_y[n][?]=n
    zero, op, inv = group_law(G)
    x_minus_y[0] = x_minus_y_data
    for i in range(1,n+1):
        x_minus_y[i] = x_minus_y[i-1] + n+1

    # Elements of G
    for j,Gj in enumerate(int_to_group):
        minus_Gj = inv(Gj)
        assert op(Gj, minus_Gj) == zero
        for i,Gi in enumerate(int_to_group):
            x_minus_y[i][j] = group_to_int[op(Gi,minus_Gj)]

    # Empty values
    for i in range(n+1):
        x_minus_y[n][i]=n
        x_minus_y[i][n]=n

    # A copy of the matrix
    for i,R in enumerate(M):
        for j,x in enumerate(R):
            M_c[i*k+j] = group_to_int[G(x)] if x is not None else n

    # Each row contains at most one empty entry
    if u:
        for i in range(M_nrows):
            bit = False
            for j in range(k):
                if M_c[i*k+j] == n:
                    if bit:
                        if verbose:
                            print("Row {} contains more than one empty entry".format(i))
                        sig_free(x_minus_y_data)
                        sig_free(x_minus_y)
                        sig_free(G_seen)
                        sig_free(M_c)
                        return False
                    bit = True

    # Each column contains lmbda*u empty entries
        for j in range(k):
            ii = 0
            for i in range(M_nrows):
                if M_c[i*k+j] == n:
                    ii += 1
            if ii!=lmbda*u:
                if verbose:
                    print("Column {} contains {} empty entries instead of the expected "
                          "lambda.u={}.{}={}".format(j, ii, lmbda, u, lmbda*u))
                sig_free(x_minus_y_data)
                sig_free(x_minus_y)
                sig_free(G_seen)
                sig_free(M_c)
                return False

    # We are now ready to test every pair of columns
    for i in range(k):
        for j in range(i+1,k):
            memset(G_seen, 0, (n+1)*sizeof(int))
            for ii in range(M_nrows):
                G_seen[x_minus_y[M_c[ii*k+i]][M_c[ii*k+j]]] += 1

            if G_seen[0] != mu: # Bad number of 0
                if verbose:
                    print("Columns {} and {} generate 0 exactly {} times "
                          "instead of the expected mu(={})".format(i,j,G_seen[0],mu))
                sig_free(x_minus_y_data)
                sig_free(x_minus_y)
                sig_free(G_seen)
                sig_free(M_c)
                return False

            for ii in range(1,n): # bad number of g_ii\in G
                if G_seen[ii] != lmbda:
                    if verbose:
                        print("Columns {} and {} do not generate all elements of G "
                         "exactly lambda(={}) times. The element {} appeared {} "
                         "times as a difference.".format(i,j,lmbda,int_to_group[ii],G_seen[ii]))
                    sig_free(x_minus_y_data)
                    sig_free(x_minus_y)
                    sig_free(G_seen)
                    sig_free(M_c)
                    return False

    sig_free(x_minus_y_data)
    sig_free(x_minus_y)
    sig_free(G_seen)
    sig_free(M_c)
    return True


# Cached information for OA constructions (see .pxd file for more info)

_OA_cache = <cache_entry *> sig_malloc(2*sizeof(cache_entry))
if (_OA_cache == NULL):
    sig_free(_OA_cache)
    raise MemoryError
_OA_cache[0].max_true = -1
_OA_cache[1].max_true = -1
_OA_cache_size = 2

cpdef _OA_cache_set(int k, int n, truth_value):
    r"""
    Set a value in the OA cache of existence results.

    INPUT:

    - ``k``, ``n`` -- integers

    - ``truth_value`` -- one of ``True``, ``False``, ``Unknown``
    """
    global _OA_cache, _OA_cache_size
    cdef int i
    if _OA_cache_size <= n:
        new_cache_size = n+100
        _OA_cache = <cache_entry *> sig_realloc(_OA_cache,new_cache_size*sizeof(cache_entry))
        if _OA_cache == NULL:
            sig_free(_OA_cache)
            raise MemoryError

        for i in range(_OA_cache_size,new_cache_size):
            _OA_cache[i].max_true = 0
            _OA_cache[i].min_unknown = -1
            _OA_cache[i].max_unknown = 0
            _OA_cache[i].min_false = -1

        _OA_cache_size = new_cache_size

    if truth_value is True:
        _OA_cache[n].max_true    = k if k>_OA_cache[n].max_true    else _OA_cache[n].max_true
    elif truth_value is Unknown:
        _OA_cache[n].min_unknown = k if k<_OA_cache[n].min_unknown else _OA_cache[n].min_unknown
        _OA_cache[n].max_unknown = k if k>_OA_cache[n].max_unknown else _OA_cache[n].max_unknown
    else:
        _OA_cache[n].min_false   = k if k<_OA_cache[n].min_false   else _OA_cache[n].min_false

cpdef _OA_cache_get(int k,int n):
    r"""
    Get a value from the OA cache of existence results.

    INPUT:

    - ``k``, ``n`` -- integers
    """
    if n>=_OA_cache_size:
        return None
    if k <= _OA_cache[n].max_true:
        return True
    elif (k >= _OA_cache[n].min_unknown and k <= _OA_cache[n].max_unknown):
        return Unknown
    elif k >= _OA_cache[n].min_false:
        return False

    return None

cpdef _OA_cache_construction_available(int k,int n):
    r"""
    Test if a construction is implemented using the cache's information.

    INPUT:

    - ``k``, ``n`` -- integers
    """
    if n>=_OA_cache_size:
        return Unknown
    if k <= _OA_cache[n].max_true:
        return True
    if k >= _OA_cache[n].min_unknown:
        return False
    else:
        return Unknown
