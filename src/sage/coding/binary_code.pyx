# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Optimized low-level binary code representation

Some computations with linear binary codes. Fix a basis for `\GF{2}^n`.
A linear binary code is a linear subspace of `\GF{2}^n`, together with
this choice of basis. A permutation `g \in S_n` of the fixed basis
gives rise to a permutation of the vectors, or words, in `\GF{2}^n`,
sending `(w_i)` to `(w_{g(i)})`. The permutation automorphism group of
the code `C` is the set of permutations of the basis that bijectively
map `C` to itself. Note that if `g` is such a permutation, then

.. MATH::

    g(a_i) + g(b_i) = (a_{g(i)} + b_{g(i)}) = g((a_i) + (b_i)).

Over other fields, it is also required that the map be linear, which
as per above boils down to scalar multiplication. However, over
`\GF{2},` the only scalars are 0 and 1, so the linearity condition has
trivial effect.

AUTHOR:

- Robert L Miller (Oct-Nov 2007)

* compiled code data structure
* union-find based orbit partition
* optimized partition stack class
* NICE-based partition refinement algorithm
* canonical generation function
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.string cimport memcpy
from cpython.mem cimport *
from cpython.object cimport PyObject_RichCompare
from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from sage.structure.element cimport Matrix
from sage.misc.timing import cputime
from sage.rings.integer cimport Integer
from copy import copy
from sage.data_structures.bitset_base cimport *

WORD_SIZE = sizeof(codeword) << 3

cdef enum:
    chunk_size = 8

cdef inline int min(int a, int b) noexcept:
    if a > b:
        return b
    else:
        return a

# NOTE - Since most of the functions are used from within the module, cdef'd
# functions come without an underscore, and the def'd equivalents, which are
# essentially only for doctesting and debugging, have underscores.

cdef int *hamming_weights() noexcept:
    cdef int *ham_wts
    cdef int i
    ham_wts = <int *> sig_malloc( 65536 * sizeof(int) )
    if ham_wts is NULL:
        sig_free(ham_wts)
        raise MemoryError("Memory.")
    ham_wts[0] = 0
    ham_wts[1] = 1
    ham_wts[2] = 1
    ham_wts[3] = 2
    for i from 4 <= i < 16:
        ham_wts[i] = ham_wts[i & 3] + ham_wts[(i>>2) & 3]
    for i from 16 <= i < 256:
        ham_wts[i] = ham_wts[i & 15] + ham_wts[(i>>4) & 15]
    for i from 256 <= i < 65536:
        ham_wts[i] = ham_wts[i & 255] + ham_wts[(i>>8) & 255]
    return ham_wts


def weight_dist(M):
    """
    Compute the weight distribution of the row space of `M`.

    EXAMPLES::

        sage: from sage.coding.binary_code import weight_dist
        sage: M = Matrix(GF(2), [
        ....:  [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
        ....:  [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],
        ....:  [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
        ....:  [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
        ....:  [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]])
        sage: weight_dist(M)
        [1, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0, 0, 0, 0, 0, 0, 1]
        sage: M = Matrix(GF(2), [
        ....:  [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
        ....:  [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0],
        ....:  [0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1],
        ....:  [0,0,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1]])
        sage: weight_dist(M)
        [1, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 4, 0, 0, 0, 0, 0]
        sage: M = Matrix(GF(2), [
        ....:  [1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0],
        ....:  [0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0],
        ....:  [0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0],
        ....:  [0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0],
        ....:  [0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0],
        ....:  [0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0],
        ....:  [0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0],
        ....:  [0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1]])
        sage: weight_dist(M)
        [1, 0, 0, 0, 0, 0, 68, 0, 85, 0, 68, 0, 34, 0, 0, 0, 0, 0]
    """
    cdef bitset_t word
    cdef int i,j,k, dim=M.nrows(), deg=M.ncols()
    cdef list L
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *LL = <int *> mem.malloc((deg+1) * sizeof(int))
    cdef bitset_s *basis = <bitset_s *> mem.malloc(dim * sizeof(bitset_s))
    for i from 0 <= i < dim:
        bitset_init(&basis[i], deg)
        bitset_zero(&basis[i])
        for j in M.row(i).nonzero_positions():
            bitset_set(&basis[i], j)
    for i from 0 <= i < deg+1: LL[i] = 0
    bitset_init(word, deg)
    bitset_zero(word)
    i = 0
    j = 0
    while True:
        LL[bitset_hamming_weight(word)] += 1
        i ^= 1
        k = 0
        if not i:
            while not j & (1 << k): k += 1
            k += 1
        if k == dim: break
        else:
            j ^= (1 << k)
            bitset_xor(word, word, &basis[k])
    bitset_free(word)
    L = [int(LL[i]) for i from 0 <= i < deg+1]
    for i from 0 <= i < dim:
        bitset_free(&basis[i])
    return L


def test_word_perms(t_limit=5.0):
    r"""
    Test the :class:`WordPermutation` structs for at least ``t_limit`` seconds.

    These are structures written in pure C for speed, and are tested from this
    function, which performs the following tests:

    1.  Tests :func:`create_word_perm`, which creates a :class:`WordPermutation`
        from a Python list `L` representing a permutation `i \mapsto
        L[i]`. Takes a random word and permutes it by a random list permutation,
        and tests that the result agrees with doing it the slow way.

    1b. Tests :func:`create_array_word_perm`, which creates a
        :class:`WordPermutation` from a C array. Does the same as above.

    2.  Tests :func:`create_comp_word_perm`, which creates a
        :class:`WordPermutation` as a composition of two
        :class:`WordPermutation` objects. Takes a random word and two random
        permutations, and tests that the result of permuting by the composition
        is correct.

    3.  Tests :func:`create_inv_word_perm` and :func:`create_id_word_perm`,
        which create a :class:`WordPermutation` as the inverse and identity
        permutations, resp.  Takes a random word and a random permutation, and
        tests that the result permuting by the permutation and its inverse in
        either order, and permuting by the identity both return the original
        word.

    .. NOTE::

        The functions :func:`permute_word_by_wp` and :func:`dealloc_word_perm`
        are implicitly involved in each of the above tests.

    TESTS::

        sage: from sage.coding.binary_code import test_word_perms
        sage: test_word_perms()  # long time (5s on sage.math, 2011)
    """
    cdef WordPermutation *g
    cdef WordPermutation *h
    cdef WordPermutation *i
    cdef codeword cw1, cw2, cw3
    cdef int n = sizeof(codeword) << 3
    cdef int j
    cdef int *arr = <int*> sig_malloc(n * sizeof(int))
    if arr is NULL:
        raise MemoryError("Error allocating memory.")
    from sage.misc.prandom import randint
    from sage.combinat.permutation import Permutations
    S = Permutations(list(range(n)))
    t = cputime()
    while cputime(t) < t_limit:
        word = [randint(0, 1) for _ in range(n)]
        cw1 = 0
        for j from 0 <= j < n:
            cw1 += (<codeword>word[j]) << (<codeword>j)
        # 1. test create_word_perm
        gg = S.random_element()
        g = create_word_perm(gg)
        word2 = [0]*n
        for j from 0 <= j < n:
            word2[gg[j]] = word[j]
        cw2 = permute_word_by_wp(g, cw1)
        cw3 = 0
        for j from 0 <= j < n:
            cw3 += (<codeword>word2[j]) << (<codeword>j)
        if cw3 != cw2:
            print("ERROR1")
        dealloc_word_perm(g)
        # 1b. test create_array_word_perm
        gg = S.random_element()
        for j from 0 <= j < n:
            arr[j] = gg[j]
        g = create_array_word_perm(arr, 0, n)
        word2 = [0]*n
        for j from 0 <= j < n:
            word2[gg[j]] = word[j]
        cw2 = permute_word_by_wp(g, cw1)
        cw3 = 0
        for j from 0 <= j < n:
            cw3 += (<codeword>word2[j]) << (<codeword>j)
        if cw3 != cw2:
            print("ERROR1b")
        dealloc_word_perm(g)
        # 2. test create_comp_word_perm
        gg = S.random_element()
        hh = S.random_element()
        g = create_word_perm(gg)
        h = create_word_perm(hh)
        i = create_comp_word_perm(g, h)
        word2 = [0]*n
        for j from 0 <= j < n:
            word2[gg[hh[j]]] = word[j]
        cw2 = permute_word_by_wp(i, cw1)
        cw3 = 0
        for j from 0 <= j < n:
            cw3 += (<codeword>word2[j]) << (<codeword>j)
        if cw3 != cw2:
            print("ERROR2")
        dealloc_word_perm(g)
        dealloc_word_perm(h)
        dealloc_word_perm(i)
        # 3. test create_inv_word_perm and create_id_word_perm
        gg = S.random_element()
        g = create_word_perm(gg)
        h = create_inv_word_perm(g)
        i = create_id_word_perm(n)
        cw2 = permute_word_by_wp(g, cw1)
        cw2 = permute_word_by_wp(h, cw2)
        if cw1 != cw2:
            print("ERROR3a")
        cw2 = permute_word_by_wp(h, cw1)
        cw2 = permute_word_by_wp(g, cw2)
        if cw1 != cw2:
            print("ERROR3b")
        cw2 = permute_word_by_wp(i, cw1)
        if cw1 != cw2:
            print("ERROR3c")
        dealloc_word_perm(g)
        dealloc_word_perm(h)
        dealloc_word_perm(i)
    sig_free(arr)


cdef WordPermutation *create_word_perm(object list_perm) noexcept:
    r"""
    Create a word permutation from a Python list permutation L, i.e. such that
    `i \mapsto L[i]`.
    """
    cdef int i, j, parity, comb, words_per_chunk, num_chunks = 1
    cdef codeword *images_i
    cdef codeword image
    cdef WordPermutation *word_perm = <WordPermutation *> sig_malloc( sizeof(WordPermutation) )
    if word_perm is NULL:
        raise RuntimeError("Error allocating memory.")
    word_perm.degree = len(list_perm)
    list_perm = copy(list_perm)
    while num_chunks*chunk_size < word_perm.degree:
        num_chunks += 1
    word_perm.images = <codeword **> sig_malloc(num_chunks * sizeof(codeword *))
    if word_perm.images is NULL:
        sig_free(word_perm)
        raise RuntimeError("Error allocating memory.")
    word_perm.chunk_num = num_chunks
    words_per_chunk = 1 << chunk_size
    word_perm.gate = ( (<codeword>1) << chunk_size ) - 1
    list_perm += list(range(len(list_perm), chunk_size*num_chunks))
    word_perm.chunk_words = words_per_chunk
    for i from 0 <= i < num_chunks:
        images_i = <codeword *> sig_malloc(words_per_chunk * sizeof(codeword))
        if images_i is NULL:
            for j from 0 <= j < i:
                sig_free(word_perm.images[j])
            sig_free(word_perm.images)
            sig_free(word_perm)
            raise RuntimeError("Error allocating memory.")
        word_perm.images[i] = images_i
        for j from 0 <= j < chunk_size:
            images_i[1 << j] = (<codeword>1) << list_perm[chunk_size*i + j]
        image = <codeword> 0
        parity = 0
        comb = 0
        while True:
            images_i[comb] = image
            parity ^= 1
            j = 0
            if not parity:
                while not comb & (1 << j): j += 1
                j += 1
            if j == chunk_size: break
            else:
                comb ^= (1 << j)
                image ^= images_i[1 << j]
    return word_perm

cdef WordPermutation *create_array_word_perm(int *array, int start, int degree) noexcept:
    """
    Create a word permutation of a given degree from a C array, starting at start.
    """
    cdef int i, j, cslim, parity, comb, words_per_chunk, num_chunks = 1
    cdef codeword *images_i
    cdef codeword image
    cdef WordPermutation *word_perm = <WordPermutation *> sig_malloc( sizeof(WordPermutation) )
    if word_perm is NULL:
        raise RuntimeError("Error allocating memory.")
    word_perm.degree = degree
    while num_chunks*chunk_size < word_perm.degree:
        num_chunks += 1
    word_perm.images = <codeword **> sig_malloc(num_chunks * sizeof(codeword *))
    if word_perm.images is NULL:
        sig_free(word_perm)
        raise RuntimeError("Error allocating memory.")
    word_perm.chunk_num = num_chunks
    words_per_chunk = 1 << chunk_size
    word_perm.gate = ( (<codeword>1) << chunk_size ) - 1
    word_perm.chunk_words = words_per_chunk
    for i from 0 <= i < num_chunks:
        images_i = <codeword *> sig_malloc(words_per_chunk * sizeof(codeword))
        if images_i is NULL:
            for j from 0 <= j < i:
                sig_free(word_perm.images[j])
            sig_free(word_perm.images)
            sig_free(word_perm)
            raise RuntimeError("Error allocating memory.")
        word_perm.images[i] = images_i
        cslim = min(chunk_size, degree - i*chunk_size)
        for j from 0 <= j < cslim:
            images_i[1 << j] = (<codeword>1) << array[start + chunk_size*i + j]
        image = <codeword> 0
        parity = 0
        comb = 0
        while True:
            images_i[comb] = image
            parity ^= 1
            j = 0
            if not parity:
                while not comb & (1 << j): j += 1
                j += 1
            if j == chunk_size: break
            else:
                comb ^= (1 << j)
                image ^= images_i[1 << j]
    return word_perm

cdef WordPermutation *create_id_word_perm(int degree) noexcept:
    """
    Create the identity word permutation of degree degree.
    """
    cdef int i, j, parity, comb, words_per_chunk, num_chunks = 1
    cdef codeword *images_i
    cdef codeword image
    cdef WordPermutation *word_perm = <WordPermutation *> sig_malloc( sizeof(WordPermutation) )
    if word_perm is NULL:
        raise RuntimeError("Error allocating memory.")
    word_perm.degree = degree
    while num_chunks*chunk_size < degree:
        num_chunks += 1
    word_perm.images = <codeword **> sig_malloc(num_chunks * sizeof(codeword *))
    if word_perm.images is NULL:
        sig_free(word_perm)
        raise RuntimeError("Error allocating memory.")
    word_perm.chunk_num = num_chunks
    words_per_chunk = 1 << chunk_size
    word_perm.gate = ( (<codeword>1) << chunk_size ) - 1
    word_perm.chunk_words = words_per_chunk
    for i from 0 <= i < num_chunks:
        images_i = <codeword *> sig_malloc(words_per_chunk * sizeof(codeword))
        if images_i is NULL:
            for j from 0 <= j < i:
                sig_free(word_perm.images[j])
            sig_free(word_perm.images)
            sig_free(word_perm)
            raise RuntimeError("Error allocating memory.")
        word_perm.images[i] = images_i
        for j from 0 <= j < chunk_size:
            images_i[1 << j] = (<codeword>1) << (chunk_size*i + j)
        image = <codeword> 0
        parity = 0
        comb = 0
        while True:
            images_i[comb] = image
            parity ^= 1
            j = 0
            if not parity:
                while not comb & (1 << j): j += 1
                j += 1
            if j == chunk_size: break
            else:
                comb ^= (1 << j)
                image ^= images_i[1 << j]
    return word_perm

cdef WordPermutation *create_comp_word_perm(WordPermutation *g, WordPermutation *h) noexcept:
    r"""
    Create the composition of word permutations `g \circ h`.
    """
    cdef int i, j, parity, comb, words_per_chunk, num_chunks = 1
    cdef codeword *images_i
    cdef codeword image
    cdef WordPermutation *word_perm = <WordPermutation *> sig_malloc( sizeof(WordPermutation) )
    if word_perm is NULL:
        raise RuntimeError("Error allocating memory.")
    word_perm.degree = g.degree
    while num_chunks*chunk_size < word_perm.degree:
        num_chunks += 1
    word_perm.images = <codeword **> sig_malloc(num_chunks * sizeof(codeword *))
    if word_perm.images is NULL:
        sig_free(word_perm)
        raise RuntimeError("Error allocating memory.")
    word_perm.chunk_num = num_chunks
    words_per_chunk = 1 << chunk_size
    word_perm.gate = ( (<codeword>1) << chunk_size ) - 1
    word_perm.chunk_words = words_per_chunk
    for i from 0 <= i < num_chunks:
        images_i = <codeword *> sig_malloc(words_per_chunk * sizeof(codeword))
        if images_i is NULL:
            for j from 0 <= j < i:
                sig_free(word_perm.images[j])
            sig_free(word_perm.images)
            sig_free(word_perm)
            raise RuntimeError("Error allocating memory.")
        word_perm.images[i] = images_i
        for j from 0 <= j < chunk_size:
            image = (<codeword>1) << (chunk_size*i + j)
            image = permute_word_by_wp(h, image)
            image = permute_word_by_wp(g, image)
            images_i[1 << j] = image
        image = <codeword> 0
        parity = 0
        comb = 0
        while True:
            images_i[comb] = image
            parity ^= 1
            j = 0
            if not parity:
                while not comb & (1 << j): j += 1
                j += 1
            if j == chunk_size: break
            else:
                comb ^= (1 << j)
                image ^= images_i[1 << j]
    return word_perm

cdef WordPermutation *create_inv_word_perm(WordPermutation *g) noexcept:
    r"""
    Create the inverse `g^{-1}` of the word permutation of `g`.
    """
    cdef int i, j
    cdef int *array = <int *> sig_malloc( g.degree * sizeof(int) )
    cdef codeword temp
    cdef WordPermutation *w
    for i from 0 <= i < g.degree:
        j = 0
        temp = permute_word_by_wp(g, (<codeword>1) << i)
        while not ((<codeword>1) << j) & temp:
            j += 1
        array[j] = i
    w = create_array_word_perm(array, 0, g.degree)
    sig_free(array)
    return w

cdef int dealloc_word_perm(WordPermutation *wp) noexcept:
    """
    Free the memory used by a word permutation.
    """
    cdef int i
    for i from 0 <= i < wp.chunk_num:
        sig_free(wp.images[i])
    sig_free(wp.images)
    sig_free(wp)

cdef codeword permute_word_by_wp(WordPermutation *wp, codeword word) noexcept:
    """
    Return the codeword obtained by applying the permutation wp to word.
    """
    cdef int num_chunks = wp.chunk_num
    cdef int i
    cdef codeword gate = wp.gate
    cdef codeword image = 0
    cdef codeword **images = wp.images
    for i from 0 <= i < num_chunks:
        image += images[i][(word >> i*chunk_size) & gate]
    return image


def test_expand_to_ortho_basis(B=None):
    r"""
    This function is written in pure C for speed, and is tested from this
    function.

    INPUT:

    - ``B`` -- a :class:`BinaryCode` in standard form

    OUTPUT:

    An array of codewords which represent the expansion of a basis for `B` to a
    basis for `(B^\prime)^\perp`, where `B^\prime = B` if the all-ones vector 1
    is in `B`, otherwise `B^\prime = \text{span}(B,1)` (note that this guarantees
    that all the vectors in the span of the output have even weight).

    TESTS::

        sage: from sage.coding.binary_code import test_expand_to_ortho_basis, BinaryCode
        sage: M = Matrix(GF(2), [[1,1,1,1,1,1,0,0,0,0], [0,0,1,1,1,1,1,1,1,1]])
        sage: B = BinaryCode(M)
        sage: B.put_in_std_form()
        0
        sage: test_expand_to_ortho_basis(B=B)
        INPUT CODE:
        Binary [10,2] linear code, generator matrix
        [1010001111]
        [0101111111]
        Expanding to the basis of an orthogonal complement...
        Basis:
        0010000010
        0001000010
        0000100001
        0000010001
        0000001001
    """
    cdef codeword *output
    cdef int k=0, i
    cdef BinaryCode C
    if not isinstance(B, BinaryCode):
        raise TypeError()
    C = B
    print("INPUT CODE:")
    print(C)
    print("Expanding to the basis of an orthogonal complement...")
    output = expand_to_ortho_basis(C, C.ncols)
    print("Basis:")
    while output[k]:
        k += 1
    for i from 0 <= i < k:
        print(''.join(reversed(Integer(output[i]).binary().zfill(C.ncols))))
    sig_free(output)


cdef codeword *expand_to_ortho_basis(BinaryCode B, int n) noexcept:
    r"""
    INPUT:

    - ``B`` -- a BinaryCode in standard form
    - ``n`` -- the degree

    OUTPUT:

    An array of codewords which represent the expansion of a basis for `B` to a
    basis for `(B^\prime)^\perp`, where `B^\prime = B` if the all-ones vector 1
    is in `B`, otherwise `B^\prime = \text{span}(B,1)` (note that this guarantees
    that all the vectors in the span of the output have even weight).
    """
    # assumes B is already in standard form
    cdef codeword *basis
    cdef codeword word = 0, temp, new, pivots = 0, combo, parity
    cdef codeword n_gate = (~<codeword>0) >> ( (sizeof(codeword)<<3) - n)
    cdef int i, j, m, k = B.nrows, dead, d
    cdef WordPermutation *wp
    basis = <codeword *> sig_malloc( (n+1) * sizeof(codeword) )
    if basis is NULL:
        raise MemoryError()
    for i from 0 <= i < k:
        basis[i] = B.basis[i]
        word ^= basis[i]
    # If 11...1 is already a word of the code,
    # then being orthogonal to the code guarantees
    # being even weight. Otherwise, add this in.
    word = (~word) & n_gate
    if word:
        basis[k] = word
        temp = (<codeword>1) << k
        i = k
        while not word & temp:
            temp = temp << 1
            i += 1
        for j from 0 <= j < k:
            if temp & basis[j]:
                basis[j] ^= word
        temp += (<codeword>1 << k) - 1
        i = k
        word = <codeword>1 << k
        k += 1
    else:  # NOTE THIS WILL NEVER HAPPEN AS CURRENTLY SET UP!
        temp = (<codeword>1 << k) - 1
        i = k
        word = <codeword>1 << k
    # Now:
    # k is the length of the basis so far
    j = k
    while i < n:
        # we are now looking at the ith free variable,
        # word has a 1 in the ith place, and
        # j is the current row we are putting in basis
        new = 0
        for m from 0 <= m < k:
            if basis[m] & word:
                new ^= basis[m]
        basis[j] = (new & temp) + word
        j += ((word ^ temp) >> i) & 1
        i += 1
        word = word << 1
    temp = (<codeword>1 << B.nrows) - 1
    for i from k <= i < n:
        basis[i-k] = basis[i] ^ B.words[basis[i] & temp]
    k = n-k
    i = 0
    word = (<codeword>1 << B.nrows)
    while i < k and (word & n_gate):
        m = i
        while m < k and not basis[m] & word:
            m += 1
        if m < k:
            pivots += word
            if m != i:
                new = basis[i]
                basis[i] = basis[m]
                basis[m] = new
            for j from 0 <= j < i:
                if basis[j] & word:
                    basis[j] ^= basis[i]
            for j from i < j < k:
                if basis[j] & word:
                    basis[j] ^= basis[i]
            i += 1
        word = word << 1
    for j from i <= j < n:
        basis[j] = 0
    # now basis is length i
    perm = list(range(B.nrows))
    perm_c = []
    for j from B.nrows <= j < B.ncols:
        if (<codeword>1 << j) & pivots:
            perm.append(j)
        else:
            perm_c.append(j)
    perm.extend(perm_c)
    perm.extend(list(range(B.ncols, n)))
    perm_c = [0]*n
    for j from 0 <= j < n:
        perm_c[perm[j]] = j
    wp = create_word_perm(perm_c)
    for j from 0 <= j < i:
        basis[j] = permute_word_by_wp(wp, basis[j])
    for j from 0 <= j < B.nrows:
        B.basis[j] = permute_word_by_wp(wp, B.basis[j])
    dealloc_word_perm(wp)
    word = 0
    parity = 0
    combo = 0
    while True:
        B.words[combo] = word
        parity ^= 1
        j = 0
        if not parity:
            while not combo & (1 << j): j += 1
            j += 1
        if j == B.nrows: break
        else:
            combo ^= (1 << j)
            word ^= B.basis[j]
    return basis

cdef class BinaryCode:
    """
    Minimal, but optimized, binary code object.

    EXAMPLES::

        sage: import sage.coding.binary_code
        sage: from sage.coding.binary_code import *
        sage: M = Matrix(GF(2), [[1,1,1,1]])
        sage: B = BinaryCode(M)     # create from matrix
        sage: C = BinaryCode(B, 60) # create using glue
        sage: D = BinaryCode(C, 240)
        sage: E = BinaryCode(D, 85)
        sage: B
        Binary [4,1] linear code, generator matrix
        [1111]
        sage: C
        Binary [6,2] linear code, generator matrix
        [111100]
        [001111]
        sage: D
        Binary [8,3] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        sage: E
        Binary [8,4] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        [10101010]

        sage: M = Matrix(GF(2), [[1]*32])
        sage: B = BinaryCode(M)
        sage: B
        Binary [32,1] linear code, generator matrix
        [11111111111111111111111111111111]
    """
    def __cinit__(self, arg1, arg2=None):
        """
        Initialize.

        TESTS::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: TestSuite(B).run()
        """
        cdef int nrows, i, j, size
        cdef int nwords, other_nwords, parity, combination
        cdef codeword word, glue_word
        cdef BinaryCode other
        cdef codeword *self_words
        cdef codeword *self_basis
        cdef codeword *other_basis

        self.radix = sizeof(int) << 3

        if isinstance(arg1, Matrix):
            self.ncols = arg1.ncols()
            self.nrows = arg1.nrows()
            nrows = self.nrows
            self.nwords = 1 << nrows
            nwords = self.nwords
        elif isinstance(arg1, BinaryCode):
            other = <BinaryCode> arg1
            self.nrows = other.nrows + 1
            glue_word = <codeword> arg2
            size = 0
            while 0 < ((<codeword> 1) << size) <= glue_word:
                size += 1
            if other.ncols > size:
                self.ncols = other.ncols
            else:
                self.ncols = size
            other_nwords = other.nwords
            self.nwords = 2 * other_nwords
            nrows = self.nrows
            nwords = self.nwords
        else:
            raise NotImplementedError

        if self.nrows >= self.radix or self.ncols > self.radix:
            raise NotImplementedError("Columns and rows are stored as ints. This code is too big.")

        self.mem = MemoryAllocator()
        self.words = <codeword *> self.mem.malloc(nwords * sizeof(int))
        self.basis = <codeword *> self.mem.malloc(nrows * sizeof(int))
        self_words = self.words
        self_basis = self.basis

        if isinstance(arg1, Matrix):
            rows = arg1.rows()
            for i from 0 <= i < nrows:
                word = <codeword> 0
                for j in rows[i].nonzero_positions():
                    word += (1<<j)
                self_basis[i] = word

            word = <codeword> 0
            parity = 0
            combination = 0
            while True:
                self_words[combination] = word
                parity ^= 1
                j = 0
                if not parity:
                    while not combination & (1 << j): j += 1
                    j += 1
                if j == nrows: break
                else:
                    combination ^= (1 << j)
                    word ^= self_basis[j]

        else:  # isinstance(arg1, BinaryCode)
            other_basis = other.basis
            for i from 0 <= i < nrows-1:
                self_basis[i] = other_basis[i]
            i = nrows - 1
            self_basis[i] = glue_word

            memcpy(self_words, other.words, other_nwords*(self.radix>>3))

            for combination from 0 <= combination < other_nwords:
                self_words[combination+other_nwords] = self_words[combination] ^ glue_word

    def __reduce__(self):
        """
        Method for pickling and unpickling BinaryCodes.

        TESTS::

            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: loads(dumps(B)) == B
            True
        """
        return BinaryCode, (self.matrix(),)

    def __richcmp__(self, other, int op):
        """
        Comparison of BinaryCodes.

        TESTS::

            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: C = BinaryCode(B.matrix())
            sage: B == C
            True
        """
        if type(self) is not type(other):
            return NotImplemented
        return PyObject_RichCompare(self.matrix(), other.matrix(), op)

    def matrix(self):
        """
        Return the generator matrix of the :class:`BinaryCode`, i.e. the code is
        the rowspace of ``B.matrix()``.

        EXAMPLES::

            sage: M = Matrix(GF(2), [[1,1,1,1,0,0], [0,0,1,1,1,1]])
            sage: from sage.coding.binary_code import *
            sage: B = BinaryCode(M)
            sage: B.matrix()
            [1 1 1 1 0 0]
            [0 0 1 1 1 1]
        """
        cdef int i, j
        from sage.matrix.constructor import matrix
        from sage.rings.finite_rings.finite_field_constructor import GF
        rows = []
        for i from 0 <= i < self.nrows:
            row = [0]*self.ncols
            for j from 0 <= j < self.ncols:
                if self.basis[i] & ((<codeword>1) << j):
                    row[j] = 1
            rows.append(row)
        return matrix(GF(2), self.nrows, self.ncols, rows)

    def print_data(self):
        """
        Print all data for ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: C = BinaryCode(B, 60)
            sage: D = BinaryCode(C, 240)
            sage: E = BinaryCode(D, 85)
            sage: B.print_data() # random - actually "print(P.print_data())"
            ncols: 4
            nrows: 1
            nwords: 2
            radix: 32
            basis:
            1111
            words:
            0000
            1111
            sage: C.print_data() # random - actually "print(P.print_data())"
            ncols: 6
            nrows: 2
            nwords: 4
            radix: 32
            basis:
            111100
            001111
            words:
            000000
            111100
            001111
            110011
            sage: D.print_data() # random - actually "print(P.print_data())"
            ncols: 8
            nrows: 3
            nwords: 8
            radix: 32
            basis:
            11110000
            00111100
            00001111
            words:
            00000000
            11110000
            00111100
            11001100
            00001111
            11111111
            00110011
            11000011
            sage: E.print_data() # random - actually "print(P.print_data())"
            ncols: 8
            nrows: 4
            nwords: 16
            radix: 32
            basis:
            11110000
            00111100
            00001111
            10101010
            words:
            00000000
            11110000
            00111100
            11001100
            00001111
            11111111
            00110011
            11000011
            10101010
            01011010
            10010110
            01100110
            10100101
            01010101
            10011001
            01101001
        """
        from sage.graphs.generic_graph_pyx import int_to_binary_string
        cdef int ui
        cdef int i
        s = ''
        s += "ncols:" + str(self.ncols)
        s += "\nnrows:" + str(self.nrows)
        s += "\nnwords:" + str(self.nwords)
        s += "\nradix:" + str(self.radix)
        s += "\nbasis:\n"
        for i from 0 <= i < self.nrows:
            b = list(int_to_binary_string(self.basis[i]).zfill(self.ncols))
            b.reverse()
            b.append('\n')
            s += ''.join(b)
        s += "\nwords:\n"
        for ui from 0 <= ui < self.nwords:
            b = list(int_to_binary_string(self.words[ui]).zfill(self.ncols))
            b.reverse()
            b.append('\n')
            s += ''.join(b)

    def __repr__(self):
        """
        String representation of ``self``.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]
        """
        cdef int i, j
        s = 'Binary [%d,%d] linear code, generator matrix\n'%(self.ncols, self.nrows)
        for i from 0 <= i < self.nrows:
            s += '[' + self._word((<codeword> 1)<<i) + ']\n'
        return s

    def _word(self, coords):
        """
        Considering ``coords`` as an integer in binary, think of the 0s and 1s
        as coefficients of the basis given by ``self.matrix()``. This function
        returns a string representation of that word.

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B._word(0)
            '0000'
            sage: B._word(1)
            '1111'

        Note that behavior under input which does not represent a word in
        the code is unspecified (gives nonsense).
        """
        s = ''
        for j from 0 <= j < self.ncols:
            s += '%d'%self.is_one(coords,j)
        return s

    def _is_one(self, word, col):
        """
        Return the col-th letter of word, i.e. 0 or 1. Words are expressed
        as integers, which represent linear combinations of the rows of the
        generator matrix of the code.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]
            sage: B._is_one(7, 4)
            0
            sage: B._is_one(8, 4)
            1
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 13, 12, 15, 14])
            1
        """
        return self.is_one(word, col) != 0

    cdef int is_one(self, int word, int column) noexcept:
        return (self.words[word] & (<codeword> 1 << column)) >> column

    def _is_automorphism(self, col_gamma, word_gamma):
        """
        Check whether a given permutation is an automorphism of the code.

        INPUT:

        - ``col_gamma`` -- permutation sending i |--> col_gamma[i] acting
          on the columns
        - ``word_gamma`` -- permutation sending i |--> word_gamma[i] acting
          on the words

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 13, 12, 15, 14])
            1
        """
        cdef int i
        cdef int *_col_gamma
        cdef int *_word_gamma
        _word_gamma = <int *> sig_malloc(self.nwords * sizeof(int))
        _col_gamma = <int *> sig_malloc(self.ncols * sizeof(int))
        if _col_gamma is NULL or _word_gamma is NULL:
            if _word_gamma is not NULL: sig_free(_word_gamma)
            if _col_gamma is not NULL: sig_free(_col_gamma)
            raise MemoryError("Memory.")
        for i from 0 <= i < self.nwords:
            _word_gamma[i] = word_gamma[i]
        for i from 0 <= i < self.ncols:
            _col_gamma[i] = col_gamma[i]
        result = self.is_automorphism(_col_gamma, _word_gamma)
        sig_free(_col_gamma)
        sig_free(_word_gamma)
        return result

    cdef int is_automorphism(self, int *col_gamma, int *word_gamma) noexcept:
        cdef int i, j, self_nwords = self.nwords, self_ncols = self.ncols
        i = 1
        while i < self_nwords:
            for j from 0 <= j < self_ncols:
                if self.is_one(i, j) != self.is_one(word_gamma[i], col_gamma[j]):
                    return 0
            i = i << 1
        return 1

    def apply_permutation(self, labeling):
        """
        Apply a column permutation to the code.

        INPUT:

        - ``labeling`` -- list permutation of the columns

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: B = BinaryCode(codes.GolayCode(GF(2)).generator_matrix())
            sage: B
            Binary [24,12] linear code, generator matrix
            [100000000000101011100011]
            [010000000000111110010010]
            [001000000000110100101011]
            [000100000000110001110110]
            [000010000000110011011001]
            [000001000000011001101101]
            [000000100000001100110111]
            [000000010000101101111000]
            [000000001000010110111100]
            [000000000100001011011110]
            [000000000010101110001101]
            [000000000001010111000111]
            sage: B.apply_permutation(list(range(11,-1,-1)) + list(range(12, 24)))
            sage: B
            Binary [24,12] linear code, generator matrix
            [000000000001101011100011]
            [000000000010111110010010]
            [000000000100110100101011]
            [000000001000110001110110]
            [000000010000110011011001]
            [000000100000011001101101]
            [000001000000001100110111]
            [000010000000101101111000]
            [000100000000010110111100]
            [001000000000001011011110]
            [010000000000101110001101]
            [100000000000010111000111]
        """
        # Tests for this function implicitly test _apply_permutation_to_basis
        # and _update_words_from_basis. These functions should not be used
        # individually by the user, so they remain cdef'd.
        self._apply_permutation_to_basis(labeling)
        self._update_words_from_basis()

    cdef void _apply_permutation_to_basis(self, object labeling) noexcept:
        cdef WordPermutation *wp
        cdef int i
        wp = create_word_perm(labeling)
        for i from 0 <= i < self.nrows:
            self.basis[i] = permute_word_by_wp(wp, self.basis[i])
        dealloc_word_perm(wp)

    cdef void _update_words_from_basis(self) noexcept:
        cdef codeword word
        cdef int j, parity, combination
        word = 0
        parity = 0
        combination = 0
        while True:
            self.words[combination] = word
            parity ^= 1
            j = 0
            if not parity:
                while not combination & (1 << j): j += 1
                j += 1
            if j == self.nrows: break
            else:
                combination ^= (1 << j)
                word ^= self.basis[j]

    cpdef int put_in_std_form(self) noexcept:
        """
        Put the code in binary form, which is defined by an identity matrix on
        the left, augmented by a matrix of data.

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0], [0,0,1,1,1,1]])
            sage: B = BinaryCode(M); B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: B.put_in_std_form(); B
            0
            Binary [6,2] linear code, generator matrix
            [101011]
            [010111]
        """
        cdef codeword swap, current = 1, pivots = 0
        cdef int i, j, k, row = 0
        cdef object perm
        while row < self.nrows:
            i = row
            while i < self.nrows and not self.basis[i] & current:
                i += 1
            if i < self.nrows:
                pivots += current
                if i != row:
                    swap = self.basis[row]
                    self.basis[row] = self.basis[i]
                    self.basis[i] = swap
                for j from 0 <= j < row:
                    if self.basis[j] & current:
                        self.basis[j] ^= self.basis[row]
                for j from row < j < self.nrows:
                    if self.basis[j] & current:
                        self.basis[j] ^= self.basis[row]
                row += 1
            current = current << 1
        perm = [0]*self.ncols
        j = 0
        k = self.nrows
        for i from 0 <= i < self.ncols:
            if ((<codeword>1) << i) & pivots:
                perm[i] = j
                j += 1
            else:
                perm[i] = k
                k += 1
        self._apply_permutation_to_basis(perm)
        self._update_words_from_basis()


cdef class OrbitPartition:
    """
    Structure which keeps track of which vertices are equivalent
    under the part of the automorphism group that has already been
    seen, during search. Essentially a disjoint-set data structure*,
    which also keeps track of the minimum element and size of each
    cell of the partition, and the size of the partition.

    See :wikipedia:`Disjoint-set_data_structure`
    """
    def __cinit__(self, int nrows, int ncols):
        """
        Initialize.

        TESTS::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: TestSuite(O).run(skip='_test_pickling')
        """
        cdef int col
        cdef int nwords, word
        nwords = (1 << nrows)
        self.nwords = nwords
        self.ncols = ncols
        self.mem = MemoryAllocator()
        self.wd_parent = <int *> self.mem.malloc(nwords * sizeof(int))
        self.wd_rank = <int *> self.mem.malloc(nwords * sizeof(int))
        self.wd_min_cell_rep = <int *> self.mem.malloc(nwords * sizeof(int))
        self.wd_size = <int *> self.mem.malloc(nwords * sizeof(int))
        self.col_parent = <int *> self.mem.malloc(ncols * sizeof(int))
        self.col_rank = <int *> self.mem.malloc(ncols * sizeof(int))
        self.col_min_cell_rep = <int *> self.mem.malloc(ncols * sizeof(int))
        self.col_size = <int *> self.mem.malloc(ncols * sizeof(int))
        for word from 0 <= word < nwords:
            self.wd_parent[word] = word
            self.wd_rank[word] = 0
            self.wd_min_cell_rep[word] = word
            self.wd_size[word] = 1
        for col from 0 <= col < ncols:
            self.col_parent[col] = col
            self.col_rank[col] = 0
            self.col_min_cell_rep[col] = col
            self.col_size[col] = 1

    def __repr__(self):
        """
        Return a string representation of the orbit partition.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
        """
        cdef int i
        cdef int j
        s = 'OrbitPartition on %d words and %d columns. Data:\n'%(self.nwords, self.ncols)
#        s += 'Parents::\n'
        s += 'Words:\n'
        for i from 0 <= i < self.nwords:
            s += '%d,'%self.wd_parent[i]
        s = s[:-1] + '\nColumns:\n'
        for j from 0 <= j < self.ncols:
            s += '%d,'%self.col_parent[j]
#        s = s[:-1] + '\n'
#        s += 'Min Cell Reps::\n'
#        s += 'Words:\n'
#        for i from 0 <= i < self.nwords:
#            s += '%d,'%self.wd_min_cell_rep[i]
#        s = s[:-1] + '\nColumns:\n'
#        for j from 0 <= j < self.ncols:
#            s += '%d,'%self.col_min_cell_rep[j]
        return s[:-1]

    def _wd_find(self, word):
        """
        Return the root of word.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_find(12)
            12
        """
        return self.wd_find(word)

    cdef int wd_find(self, int word) noexcept:
        if self.wd_parent[word] == word:
            return word
        else:
            self.wd_parent[word] = self.wd_find(self.wd_parent[word])
        return self.wd_parent[word]

    def _wd_union(self, x, y):
        """
        Join the cells containing x and y.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_union(1,10)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,1,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_find(10)
            1
        """
        self.wd_union(x, y)

    cdef void wd_union(self, int x, int y) noexcept:
        cdef int x_root, y_root
        x_root = self.wd_find(x)
        y_root = self.wd_find(y)
        if self.wd_rank[x_root] > self.wd_rank[y_root]:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[x_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[x_root] += self.wd_size[y_root]
        elif self.wd_rank[x_root] < self.wd_rank[y_root]:
            self.wd_parent[x_root] = y_root
            self.wd_min_cell_rep[y_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[y_root] += self.wd_size[x_root]
        elif x_root != y_root:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[x_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[x_root] += self.wd_size[y_root]
            self.wd_rank[x_root] += 1

    def _col_find(self, col):
        """
        Return the root of col.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._col_find(6)
            6
        """
        return self.col_find(col)

    cdef int col_find(self, int col) noexcept:
        if self.col_parent[col] == col:
            return col
        else:
            self.col_parent[col] = self.col_find(self.col_parent[col])
            return self.col_parent[col]

    def _col_union(self, x, y):
        """
        Join the cells containing x and y.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._col_union(1,4)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,1,5,6,7
            sage: O._col_find(4)
            1
        """
        self.col_union(x, y)

    cdef void col_union(self, int x, int y) noexcept:
        cdef int x_root, y_root
        x_root = self.col_find(x)
        y_root = self.col_find(y)
        if self.col_rank[x_root] > self.col_rank[y_root]:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[x_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[x_root] += self.col_size[y_root]
        elif self.col_rank[x_root] < self.col_rank[y_root]:
            self.col_parent[x_root] = y_root
            self.col_min_cell_rep[y_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[y_root] += self.col_size[x_root]
        elif x_root != y_root:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[x_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[x_root] += self.col_size[y_root]
            self.col_rank[x_root] += 1

    def _merge_perm(self, col_gamma, wd_gamma):
        """
        Merges the cells of ``self`` under the given permutation. If gamma[a] = b,
        then after merge_perm, a and b will be in the same cell. Returns 0 if
        nothing was done, otherwise returns 1.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._merge_perm([1,0,3,2,4,5,6,7], [0,1,2,3,4,5,6,7,9,8,11,10,13,12,15,14])
            1
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,8,10,10,12,12,14,14
            Columns:
            0,0,2,2,4,5,6,7
        """
        cdef int i
        cdef int *_col_gamma
        cdef int *_wd_gamma
        _wd_gamma = <int *> sig_malloc(self.nwords * sizeof(int))
        _col_gamma = <int *> sig_malloc(self.ncols * sizeof(int))
        if _col_gamma is NULL or _wd_gamma is NULL:
            if _wd_gamma is not NULL: sig_free(_wd_gamma)
            if _col_gamma is not NULL: sig_free(_col_gamma)
            raise MemoryError("Memory.")
        for i from 0 <= i < self.nwords:
            _wd_gamma[i] = wd_gamma[i]
        for i from 0 <= i < self.ncols:
            _col_gamma[i] = col_gamma[i]
        result = self.merge_perm(_col_gamma, _wd_gamma)
        sig_free(_col_gamma)
        sig_free(_wd_gamma)
        return result

    cdef int merge_perm(self, int *col_gamma, int *wd_gamma) noexcept:
        cdef int i, gamma_i_root
        cdef int j, gamma_j_root, return_value = 0
        cdef int *self_wd_parent = self.wd_parent
        cdef int *self_col_parent = self.col_parent
        for i from 0 <= i < self.nwords:
            gamma_i_root = self.wd_find(wd_gamma[i])
            if gamma_i_root != i:
                return_value = 1
                self.wd_union(i, gamma_i_root)
        for j from 0 <= j < self.ncols:
            gamma_j_root = self.col_find(col_gamma[j])
            if gamma_j_root != j:
                return_value = 1
                self.col_union(j, gamma_j_root)
        return return_value

cdef class PartitionStack:
    """
    Partition stack structure for traversing the search tree during automorphism
    group computation.
    """
    def __cinit__(self, arg1, arg2=None):
        """
        Initialize.

        TESTS::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: TestSuite(P).run(skip='_test_pickling')
        """
        cdef int k, nwords, ncols, sizeof_int
        cdef PartitionStack other = None
        cdef int *wd_ents
        cdef int *wd_lvls
        cdef int *col_ents
        cdef int *col_lvls
        cdef int *col_degs
        cdef int *col_counts
        cdef int *col_output
        cdef int *wd_degs
        cdef int *wd_counts
        cdef int *wd_output
        sizeof_int = sizeof(int)

        try:
            self.nrows = <int> arg1
            self.nwords = 1 << self.nrows
            self.ncols = <int> arg2
        except Exception:
            other = arg1
            self.nrows = other.nrows
            self.nwords = other.nwords
            self.ncols = other.ncols

        self.radix = sizeof_int << 3
        self.flag = (1 << (self.radix-1))

        # data
        self.mem = MemoryAllocator()
        self.wd_ents = <int *> self.mem.malloc(self.nwords * sizeof_int)
        self.wd_lvls = <int *> self.mem.malloc(self.nwords * sizeof_int)
        self.col_ents = <int *> self.mem.malloc(self.ncols  * sizeof_int)
        self.col_lvls = <int *> self.mem.malloc(self.ncols  * sizeof_int)

        # scratch space
        self.col_degs = <int *> self.mem.malloc( self.ncols  * sizeof_int )
        self.col_counts = <int *> self.mem.malloc( self.nwords * sizeof_int )
        self.col_output = <int *> self.mem.malloc( self.ncols  * sizeof_int )
        self.wd_degs = <int *> self.mem.malloc( self.nwords * sizeof_int )
        self.wd_counts = <int *> self.mem.malloc( (self.ncols+1)  * sizeof_int )
        self.wd_output = <int *> self.mem.malloc( self.nwords * sizeof_int )

        nwords = self.nwords
        ncols = self.ncols

        if other:
            memcpy(self.wd_ents,  other.wd_ents, self.nwords * sizeof_int)
            memcpy(self.wd_lvls,  other.wd_lvls, self.nwords * sizeof_int)
            memcpy(self.col_ents, other.col_ents, self.ncols * sizeof_int)
            memcpy(self.col_lvls, other.col_lvls, self.ncols * sizeof_int)
        else:
            wd_ents = self.wd_ents
            wd_lvls = self.wd_lvls
            col_ents = self.col_ents
            col_lvls = self.col_lvls
            for k from 0 <= k < nwords-1:
                wd_ents[k] = k
                wd_lvls[k] = 2*ncols
            for k from 0 <= k < ncols-1:
                col_ents[k] = k
                col_lvls[k] = 2*ncols
            wd_ents[nwords-1] = nwords-1
            wd_lvls[nwords-1] = -1
            col_ents[ncols-1] = ncols-1
            col_lvls[ncols-1] = -1

        col_degs = self.col_degs
        col_counts = self.col_counts
        col_output = self.col_output
        wd_degs = self.wd_degs
        wd_counts = self.wd_counts
        wd_output = self.wd_output
        for k from 0 <= k < ncols:
            col_degs[k]=0
            col_output[k]=0
            wd_counts[k]=0
        wd_counts[ncols]=0
        for k from 0 <= k < nwords:
            col_counts[k]=0
            wd_degs[k]=0
            wd_output[k]=0

    def __dealloc__(self):
        if self.basis_locations:
            sig_free(self.basis_locations)

    def print_data(self):
        """
        Print all data for ``self``.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: print(P.print_data())
            nwords:4
            nrows:2
            ncols:6
            radix:32
            wd_ents:
            0
            1
            2
            3
            wd_lvls:
            12
            12
            12
            -1
            col_ents:
            0
            1
            2
            3
            4
            5
            col_lvls:
            12
            12
            12
            12
            12
            -1
            col_degs:
            0
            0
            0
            0
            0
            0
            col_counts:
            0
            0
            0
            0
            col_output:
            0
            0
            0
            0
            0
            0
            wd_degs:
            0
            0
            0
            0
            wd_counts:
            0
            0
            0
            0
            0
            0
            0
            wd_output:
            0
            0
            0
            0
        """
        cdef int i, j
        s = ''
        s += "nwords:" + str(self.nwords) + '\n'
        s += "nrows:" + str(self.nrows) + '\n'
        s += "ncols:" + str(self.ncols) + '\n'
        s += "radix:" + str(self.radix) + '\n'
        s += "wd_ents:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_ents[i]) + '\n'
        s += "wd_lvls:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_lvls[i]) + '\n'
        s += "col_ents:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_ents[i]) + '\n'
        s += "col_lvls:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_lvls[i]) + '\n'
        s += "col_degs:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_degs[i]) + '\n'
        s += "col_counts:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.col_counts[i]) + '\n'
        s += "col_output:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_output[i]) + '\n'
        s += "wd_degs:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_degs[i]) + '\n'
        s += "wd_counts:" + '\n'
        for i from 0 <= i < self.ncols + 1:
            s += str(self.wd_counts[i]) + '\n'
        s += "wd_output:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_output[i]) + '\n'
        if self.basis_locations:
            s += "basis_locations:" + '\n'
            j = 1
            while (1 << j) < self.nwords:
                j += 1
            for i from 0 <= i < j:
                s += str(self.basis_locations[i]) + '\n'
        return s

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
        """
        cdef int i, j, k
        s = ''
        last = ''
        current = ''
        for k from 0 <= k < 2*self.ncols:
            current = self._repr_at_k(k)
            if current == last: break
            s += current
            last = current
        return s

    def _repr_at_k(self, k):
        """
        Give a string representing the partition at level k.

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6); P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            sage: P._repr_at_k(0)
            '({0,1,2,3})  ({0,1,2,3,4,5})\n'
        """
        s = '({'
        for j from 0 <= j < self.nwords:
            s += str(self.wd_ents[j])
            if self.wd_lvls[j] <= k:
                s += '},{'
            else:
                s += ','
        s = s[:-2] + ')  '
        s += '({'
        for j from 0 <= j < self.ncols:
            s += str(self.col_ents[j])
            if self.col_lvls[j] <= k:
                s += '},{'
            else:
                s += ','
        s = s[:-2] + ')\n'
        return s

    def _is_discrete(self, k):
        """
        Return whether the partition at level k is discrete.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sort_wds(0, [0,2,3,1], 5)
            0
            sage: P
            ({0,3,1,2})  ({0,1,2,3,4,5})
            ({0,3,1,2})  ({0},{1,2,3,4,5})
            ({0,3,1,2})  ({0},{1},{2,3,4,5})
            ({0,3,1,2})  ({0},{1},{2},{3,4,5})
            ({0,3,1,2})  ({0},{1},{2},{3},{4,5})
            ({0},{3},{1},{2})  ({0},{1},{2},{3},{4},{5})
            sage: P._is_discrete(4)
            0
            sage: P._is_discrete(5)
            1
        """
        return self.is_discrete(k)

    cdef int is_discrete(self, int k) noexcept:
        cdef int i, self_ncols = self.ncols, self_nwords = self.nwords
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_wd_lvls = self.wd_lvls
        for i from 0 <= i < self_ncols:
            if self_col_lvls[i] > k:
                return 0
        for i from 0 <= i < self_nwords:
            if self_wd_lvls[i] > k:
                return 0
        return 1

    def _num_cells(self, k):
        """
        Return the number of cells in the partition at level k.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._num_cells(3)
            5
        """
        return self.num_cells(k)

    cdef int num_cells(self, int k) noexcept:
        cdef int i, j = 0
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_lvls = self.col_lvls
        for i from 0 <= i < self.nwords:
            if self_wd_lvls[i] <= k:
                j += 1
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] <= k:
                j += 1
        return j

    def _sat_225(self, k):
        """
        Return whether the partition at level k satisfies the hypotheses of
        Lemma 2.25 in Brendan McKay's Practical Graph Isomorphism paper (see
        sage/graphs/graph_isom.pyx.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sat_225(3)
            0
            sage: P._sat_225(4)
            1
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
        """
        return self.sat_225(k)

    cdef int sat_225(self, int k) noexcept:
        cdef int i, n = self.nwords + self.ncols, in_cell = 0
        cdef int nontrivial_cells = 0, total_cells = self.num_cells(k)
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_lvls = self.col_lvls
        if n <= total_cells + 4:
            return 1
        for i from 0 <= i < self.nwords:
            if self_wd_lvls[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        in_cell = 0
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        if n == total_cells + nontrivial_cells:
            return 1
        if n == total_cells + nontrivial_cells + 1:
            return 1
        return 0

#    def _new_min_cell_reps(self, k): #TODO
#        """
#        Returns an integer whose bits represent which columns are minimal cell
#        representatives.
#
#        EXAMPLES::
#
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._min_cell_reps(2)
#            sage: Integer(a).binary()
#            '111'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.min_cell_reps(k)
#
#    cdef int min_cell_reps(self, int k):
#        cdef int i
#        cdef int reps = 1
#        cdef int *self_col_lvls = self.col_lvls
#        for i from 0 < i < self.ncols:
#            if self_col_lvls[i-1] <= k:
#                reps += (1 << i)
#        return reps
#
    cdef void new_min_cell_reps(self, int k, unsigned int *Omega, int start) noexcept:
        cdef int i, j
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents
        cdef int *self_wd_ents = self.wd_ents
        cdef int reps = (1 << self_col_ents[0]), length, word
        cdef int radix = self.radix, nwords = self.nwords, ncols = self.ncols
        length = 1 + nwords/radix
        if nwords%radix:
            length += 1
        for i from 0 <= i < length:
            Omega[start+i] = 0
        for i from 0 < i < ncols:
            Omega[start] += ((self_col_lvls[i-1] <= k) << self_col_ents[i])
        Omega[start+1] = (1 << self_wd_ents[0])
        for i from 0 < i < nwords:
            if self_wd_lvls[i-1] <= k:
                word = self_wd_lvls[i-1]
                Omega[start+1+word/radix] += (1 << word%radix)

#    def _fixed_cols(self, mcrs, k): #TODO
#        """
#        Returns an integer whose bits represent which columns are fixed. For
#        efficiency, mcrs is the output of min_cell_reps.
#
#        EXAMPLES::
#
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._fixed_cols(7, 2)
#            sage: Integer(a).binary()
#            '11'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.fixed_cols(mcrs, k)
#
#    cdef int fixed_cols(self, int mcrs, int k):
#        cdef int i
#        cdef int fixed = 0
#        cdef int *self_col_lvls = self.col_lvls
#        for i from 0 <= i < self.ncols:
#            if self_col_lvls[i] <= k:
#                fixed += (1 << i)
#        return fixed & mcrs
#
    cdef void fixed_vertices(self, int k, unsigned int *Phi, unsigned int *Omega, int start) noexcept:
        cdef int i, j, length, ell, fixed = 0
        cdef int radix = self.radix, nwords = self.nwords, ncols = self.ncols
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents
        cdef int *self_wd_ents = self.wd_ents
        for i from 0 <= i < ncols:
            fixed += ((self_col_lvls[i] <= k) << self_col_ents[i])
        Phi[start] = fixed & Omega[start]
        # zero out the rest of Phi
        length = 1 + nwords/self.radix
        if nwords%self.radix:
            length += 1
        for i from 0 < i < length:
            Phi[start+i] = 0
        for i from 0 <= i < nwords:
            ell = self_wd_ents[i]
            Phi[start+1+ell/radix] = ((self_wd_lvls[i] <= k) << ell%radix)
        for i from 0 < i < length:
            Phi[i] &= Omega[i]

#    def _first_smallest_nontrivial(self, k): #TODO
#        """
#        Returns an integer representing the first, smallest nontrivial cell of columns.
#
#        EXAMPLES::
#
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._first_smallest_nontrivial(2)
#            sage: Integer(a).binary().zfill(32)
#            '00000000000000000000000000111100'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.first_smallest_nontrivial(k)
#
#    cdef int first_smallest_nontrivial(self, int k):
#        cdef int cell
#        cdef int i = 0, j = 0, location = 0, ncols = self.ncols
#        cdef int *self_col_lvls = self.col_lvls
#        while True:
#            if self_col_lvls[i] <= k:
#                if i != j and ncols > i - j + 1:
#                    ncols = i - j + 1
#                    location = j
#                j = i + 1
#            if self_col_lvls[i] == -1: break
#            i += 1
#        # location now points to the beginning of the first, smallest,
#        # nontrivial cell
#        j = location
#        self.v = self.col_ents[j]
#        while True:
#            if self_col_lvls[j] <= k: break
#            j += 1
#        # j now points to the last element of the cell
#        i = self.radix - j - 1                 # the cell is represented in binary, reading from the right:
#        cell = (~0 << location) ^ (~0 << j+1)  # <-------            self.radix               ----->
#        return cell                            # [0]*(radix-j-1) + [1]*(j-location+1) + [0]*location
#
    cdef int new_first_smallest_nontrivial(self, int k, unsigned int *W, int start) noexcept:
        cdef int ell
        cdef int i = 0, j = 0, location = 0, min = self.ncols, nwords = self.nwords
        cdef int min_is_col = 1, radix = self.radix
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents
        cdef int *self_wd_ents = self.wd_ents
        while True:
            if self_col_lvls[i] <= k:
                if i != j and min > i - j + 1:
                    min = i - j + 1
                    location = j
                j = i + 1
            if self_col_lvls[i] == -1: break
            i += 1
#        i = 0; j = 0
#        while True:
#            if self_wd_lvls[i] <= k:
#                if i != j and min > i - j + 1:
#                    min = i - j + 1
#                    min_is_col = 0
#                    location = j
#                j = i + 1
#            if self_wd_lvls[i] == -1: break
#            i += 1
        # location now points to the beginning of the first, smallest,
        # nontrivial cell
        j = location
        # zero out this level of W:
        ell = 1 + nwords/radix
        if nwords%radix:
            ell += 1
        for i from 0 <= i < ell:
            W[start+i] = 0
        if min_is_col:
            while True:
                if self_col_lvls[j] <= k: break
                j += 1
            # j now points to the last element of the cell
            i = location
            while i <= j:
                W[start] ^= (1 << self_col_ents[i])
                i += 1
            return self_col_ents[location]
        else:
            while True:
                if self_wd_lvls[j] <= k: break
                j += 1
            # j now points to the last element of the cell
            i = location
            while i <= j:
                ell = self_wd_ents[i]
                W[start+1+ell/radix] ^= (1 << ell%radix)
                i += 1
            return self_wd_ents[location] ^ self.flag

    def _dangerous_dont_use_set_ents_lvls(self, col_ents, col_lvls, wd_ents, wd_lvls):
        """
        For debugging only.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            sage: P._dangerous_dont_use_set_ents_lvls([99]*6, [0,3,2,3,5,-1], [4,3,5,6], [3,2,1,-1])
            sage: P
            ({4,3,5,6})  ({99},{99,99,99,99,99})
            ({4,3,5},{6})  ({99},{99,99,99,99,99})
            ({4,3},{5},{6})  ({99},{99,99},{99,99,99})
            ({4},{3},{5},{6})  ({99},{99},{99},{99},{99,99})
        """
        cdef int i
        for i from 0 <= i < len(col_ents):
            self.col_ents[i] = col_ents[i]
        for i from 0 <= i < len(col_lvls):
            self.col_lvls[i] = col_lvls[i]
        for i from 0 <= i < len(wd_ents):
            self.wd_ents[i] = wd_ents[i]
        for i from 0 <= i < len(wd_lvls):
            self.wd_lvls[i] = wd_lvls[i]

    def _col_percolate(self, start, end):
        """
        Do one round of bubble sort on ents.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(list(range(5,-1,-1)), [1,2,2,3,3,-1], list(range(3,-1,-1)), [1,1,2,-1])
            sage: P
            ({3,2,1,0})  ({5,4,3,2,1,0})
            ({3},{2},{1,0})  ({5},{4,3,2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})  ({0,5,4,3,2,1})
            ({0},{3},{2,1})  ({0},{5,4,3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3},{2},{1})
        """
        self.col_percolate(start, end)

    cdef void col_percolate(self, int start, int end) noexcept:
        cdef int i, temp
        cdef int *self_col_ents = self.col_ents
        for i from end >= i > start:
            if self_col_ents[i] < self_col_ents[i-1]:
                temp = self.col_ents[i]
                self_col_ents[i] = self_col_ents[i-1]
                self_col_ents[i-1] = temp

    def _wd_percolate(self, start, end):
        """
        Do one round of bubble sort on ents.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(list(range(5,-1,-1)), [1,2,2,3,3,-1], list(range(3,-1,-1)), [1,1,2,-1])
            sage: P
            ({3,2,1,0})  ({5,4,3,2,1,0})
            ({3},{2},{1,0})  ({5},{4,3,2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})  ({0,5,4,3,2,1})
            ({0},{3},{2,1})  ({0},{5,4,3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3},{2},{1})
        """
        self.wd_percolate(start, end)

    cdef void wd_percolate(self, int start, int end) noexcept:
        cdef int i, temp
        cdef int *self_wd_ents = self.wd_ents
        for i from end >= i > start:
            if self_wd_ents[i] < self_wd_ents[i-1]:
                temp = self.wd_ents[i]
                self_wd_ents[i] = self_wd_ents[i-1]
                self_wd_ents[i-1] = temp

#    def _split_column(self, int v, int k): #TODO
#        """
#        Split column v out, placing it before the rest of the cell it was in.
#        Returns the location of the split column.
#
#        EXAMPLES::
#
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#            sage: P = PartitionStack(2, 6)
#            sage: P._split_column(0,1)
#            0
#            sage: P._split_column(2,2)
#            1
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,2,1,3,4,5})
#            ({0},{2,1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#
#        """
#        return self.split_column(v, k)
#
#    cdef int split_column(self, int v, int k):
#        cdef int i = 0, j
#        cdef int *self_col_ents = self.col_ents
#        cdef int *self_col_lvls = self.col_lvls
#        while self_col_ents[i] != v: i += 1
#        j = i
#        while self_col_lvls[i] > k: i += 1
#        if j == 0 or self_col_lvls[j-1] <= k:
#            self.col_percolate(j+1, i)
#        else:
#            while j != 0 and self_col_lvls[j-1] > k:
#                self_col_ents[j] = self_col_ents[j-1]
#                j -= 1
#            self_col_ents[j] = v
#        self_col_lvls[j] = k
#        return j
#

    def _split_vertex(self, v, k):
        """
        Split vertex v out, placing it before the rest of the cell it was in.
        Returns the location of the split vertex.

        .. NOTE::

            There is a convention regarding whether a vertex is a word or a
            column. See the 'flag' attribute of the PartitionStack object:
            If vertex&flag is not zero, it is a word.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
        """
        return self.split_vertex(v, k)

    cdef int split_vertex(self, int v, int k) noexcept:
        cdef int i = 0, j, flag = self.flag
        cdef int *ents
        cdef int *lvls
        if v & flag:
            ents = self.wd_ents
            lvls = self.wd_lvls
            v = v ^ flag
            while ents[i] != v: i += 1
            v = v ^ flag
        else:
            ents = self.col_ents
            lvls = self.col_lvls
            while ents[i] != v: i += 1
        j = i
        while lvls[i] > k: i += 1
        if j == 0 or lvls[j-1] <= k:
            if v & self.flag:
                self.wd_percolate(j+1, i)
            else:
                self.col_percolate(j+1, i)
        else:
            while j != 0 and lvls[j-1] > k:
                ents[j] = ents[j-1]
                j -= 1
            if v & flag:
                ents[j] = v ^ flag
            else:
                ents[j] = v
        lvls[j] = k
        if v & flag:
            return j ^ flag
        else:
            return j

    def _col_degree(self, C, col, wd_ptr, k):
        """
        Return the number of words in the cell specified by wd_ptr that have a
        1 in the col-th column.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._col_degree(B, 2, 0, 2)
            2
        """
        return self.col_degree(C, col, wd_ptr, k)

    cdef int col_degree(self, BinaryCode CG, int col, int wd_ptr, int k) noexcept:
        cdef int i = 0
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_wd_ents = self.wd_ents
        while True:
            if CG.is_one(self_wd_ents[wd_ptr], col):
                i += 1
            if self_wd_lvls[wd_ptr] > k:
                wd_ptr += 1
            else:
                break
        return i

    def _wd_degree(self, C, wd, col_ptr, k):
        """
        Return the number of columns in the cell specified by col_ptr that are
        1 in wd.

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._wd_degree(B, 1, 1, 1)
            3
        """
        cdef int *ham_wts = hamming_weights()
        result = self.wd_degree(C, wd, col_ptr, k, ham_wts)
        sig_free(ham_wts)
        return result

    cdef int wd_degree(self, BinaryCode CG, int wd, int col_ptr, int k, int *ham_wts) noexcept:

        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_col_ents = self.col_ents
        cdef int mask = (1 << self_col_ents[col_ptr])
        while self_col_lvls[col_ptr] > k:
            col_ptr += 1
            mask += (1 << self_col_ents[col_ptr])
        mask &= CG.words[wd]
        return ham_wts[mask & 65535] + ham_wts[(mask >> 16) & 65535]

    def _sort_cols(self, start, degrees, k):
        """
        Essentially a counting sort, but on only one cell of the partition.

        INPUT:

        - ``start`` -- location of the beginning of the cell
        - ``k`` -- at what level of refinement the partition of interest lies
        - ``degrees`` -- the counts to sort by

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sort_cols(1, [0,2,2,1,1], 1)
            2
            sage: P
            ({0,1,2,3})  ({0,1,4,5,2,3})
            ({0,1,2,3})  ({0},{1},{4,5},{2,3})
        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.col_degs[i] = degrees[i]
        return self.sort_cols(start, k)

    cdef int sort_cols(self, int start, int k) noexcept:
        cdef int i, j, max, max_location, self_ncols = self.ncols
        cdef int self_nwords = self.nwords, ii
        cdef int *self_col_counts = self.col_counts
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_col_degs = self.col_degs
        cdef int *self_col_ents = self.col_ents
        cdef int *self_col_output = self.col_output
        for ii from 0 <= ii < self_nwords:
            self_col_counts[ii] = 0
        i = 0
        while self_col_lvls[i+start] > k:
            self_col_counts[self_col_degs[i]] += 1
            i += 1
        self_col_counts[self_col_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_col_counts[0]
        max_location = 0
        for ii from 0 < ii < self_nwords:
            if self_col_counts[ii] > max:
                max = self_col_counts[ii]
                max_location = ii
            self_col_counts[ii] += self_col_counts[ii-1]

        for j from i >= j >= 0:
            self_col_counts[self_col_degs[j]] -= 1
            self_col_output[self_col_counts[self_col_degs[j]]] = self_col_ents[start+j]

        max_location = self_col_counts[max_location] + start

        for j from 0 <= j <= i:
            self_col_ents[start+j] = self_col_output[j]

        ii = 1
        while ii < self_nwords and self_col_counts[ii] <= i:
            if self_col_counts[ii] > 0:
                self_col_lvls[start + self_col_counts[ii] - 1] = k
            self.col_percolate(start + self_col_counts[ii-1], start + self_col_counts[ii] - 1)
            ii += 1

        return max_location

    def _sort_wds(self, start, degrees, k):
        """
        Essentially a counting sort, but on only one cell of the partition.

        INPUT:

        - ``start`` -- location of the beginning of the cell
        - ``k`` -- at what level of refinement the partition of interest lies
        - ``degrees`` -- the counts to sort by

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(3, 6)
            sage: P._sort_wds(0, [0,0,3,3,3,3,2,2], 1)
            4
            sage: P
            ({0,1,6,7,2,3,4,5})  ({0,1,2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})  ({0,1,2,3,4,5})
        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.wd_degs[i] = degrees[i]
        return self.sort_wds(start, k)

    cdef int sort_wds(self, int start, int k) noexcept:
        cdef int i, j, max, max_location, self_nwords = self.nwords
        cdef int ii, self_ncols = self.ncols
        cdef int *self_wd_counts = self.wd_counts
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_wd_degs = self.wd_degs
        cdef int *self_wd_ents = self.wd_ents
        cdef int *self_wd_output = self.wd_output

        for ii from 0 <= ii < self_ncols+1:
            self_wd_counts[ii] = 0
        i = 0
        while self_wd_lvls[i+start] > k:
            self_wd_counts[self_wd_degs[i]] += 1
            i += 1
        self_wd_counts[self_wd_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_wd_counts[0]
        max_location = 0
        for ii from 0 < ii < self_ncols+1:
            if self_wd_counts[ii] > max:
                max = self_wd_counts[ii]
                max_location = ii
            self_wd_counts[ii] += self_wd_counts[ii-1]

        for j from i >= j >= 0:
            if j > i: break # cython bug with ints...
            self_wd_counts[self_wd_degs[j]] -= 1
            self_wd_output[self_wd_counts[self_wd_degs[j]]] = self_wd_ents[start+j]

        max_location = self_wd_counts[max_location] + start

        for j from 0 <= j <= i:
            self_wd_ents[start+j] = self_wd_output[j]

        ii = 1
        while ii < self_ncols+1 and self_wd_counts[ii] <= i:
            if self_wd_counts[ii] > 0:
                self_wd_lvls[start + self_wd_counts[ii] - 1] = k
            self.wd_percolate(start + self_wd_counts[ii-1], start + self_wd_counts[ii] - 1)
            ii += 1

        return max_location

    def _refine(self, k, alpha, CG):
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(1, [[0,0],[1,0]], B)
            181
            sage: P._split_vertex(0, 2)
            0
            sage: P._refine(2, [[0,0]], B)
            290
            sage: P._split_vertex(1, 3)
            1
            sage: P._refine(3, [[0,1]], B)
            463
            sage: P._split_vertex(2, 4)
            2
            sage: P._refine(4, [[0,2]], B)
            1500
            sage: P._split_vertex(3, 5)
            3
            sage: P._refine(5, [[0,3]], B)
            641
            sage: P._split_vertex(4, 6)
            4
            sage: P._refine(6, [[0,4]], B)
            1224
            sage: P._is_discrete(5)
            0
            sage: P._is_discrete(6)
            1
            sage: P
            ({0,4,6,2,13,9,11,15,10,14,12,8,7,3,1,5})  ({0,1,2,3,4,7,6,5})
            ({0},{4,6,2,13,9,11,15,10,14,12,8,7,3,1},{5})  ({0,1,2,3,4,7,6,5})
            ({0},{4,6,2,13,9,11,15},{10,14,12,8,7,3,1},{5})  ({0},{1,2,3,4,7,6,5})
            ({0},{4,6,2},{13,9,11,15},{10,14,12,8},{7,3,1},{5})  ({0},{1},{2,3,4,7,6,5})
            ({0},{4},{6,2},{13,9},{11,15},{10,14},{12,8},{7,3},{1},{5})  ({0},{1},{2},{3,4,7,6,5})
            ({0},{4},{6,2},{13,9},{11,15},{10,14},{12,8},{7,3},{1},{5})  ({0},{1},{2},{3},{4,7,6,5})
            ({0},{4},{6},{2},{13},{9},{11},{15},{10},{14},{12},{8},{7},{3},{1},{5})  ({0},{1},{2},{3},{4},{7},{6},{5})
        """
        cdef int i, alpha_length = len(alpha)
        cdef int *_alpha = <int *> sig_malloc( (self.nwords + self.ncols) * sizeof(int) )
        cdef int *ham_wts = hamming_weights()
        if _alpha is NULL:
            raise MemoryError("Memory.")
        for i from 0 <= i < alpha_length:
            if alpha[i][0]:
                _alpha[i] = alpha[i][1] ^ self.flag
            else:
                _alpha[i] = alpha[i][1]
        result = self.refine(k, _alpha, alpha_length, CG, ham_wts)
        sig_free(_alpha)
        sig_free(ham_wts)
        return result

    cdef int refine(self, int k, int *alpha, int alpha_length, BinaryCode CG, int *ham_wts) noexcept:
        cdef int q, r, s, t, flag = self.flag, self_ncols = self.ncols
        cdef int t_w, self_nwords = self.nwords, invariant = 0, i, j, m = 0
        cdef int *self_wd_degs = self.wd_degs
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_wd_ents = self.wd_ents
        cdef int *self_col_degs = self.col_degs
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_col_ents = self.col_ents
        while not self.is_discrete(k) and m < alpha_length:
            invariant += 1
            j = 0
            if alpha[m] & flag:
                while j < self_ncols:
                    i = j
                    s = 0
                    invariant += 8
                    while True:
                        self_col_degs[i-j] = self.col_degree(CG, self_col_ents[i], alpha[m]^flag, k)
                        if s == 0 and self_col_degs[i-j] != self_col_degs[0]: s = 1
                        i += 1
                        if self_col_lvls[i-1] <= k: break
                    if s:
                        invariant += 8
                        t = self.sort_cols(j, k)
                        invariant += t
                        q = m
                        while q < alpha_length:
                            if alpha[q] == j:
                                alpha[q] = t
                                break
                            q += 1
                        r = j
                        while True:
                            if r == j or self.col_lvls[r-1] == k:
                                if r != t:
                                    alpha[alpha_length] = r
                                    alpha_length += 1
                            r += 1
                            if r >= i: break
                        invariant += self.col_degree(CG, self_col_ents[i-1], alpha[m]^flag, k)
                        invariant += (i-j)
                    j = i
            else:
                while j < self.nwords:
                    i = j
                    s = 0
                    invariant += 64
                    while True:
                        self_wd_degs[i-j] = self.wd_degree(CG, self_wd_ents[i], alpha[m], k, ham_wts)
                        if s == 0 and self_wd_degs[i-j] != self_wd_degs[0]: s = 1
                        i += 1
                        if self_wd_lvls[i-1] <= k: break
                    if s:
                        invariant += 64
                        t_w = self.sort_wds(j, k)
                        invariant += t_w
                        q = m
                        j ^= flag
                        while q < alpha_length:
                            if alpha[q] == j:
                                alpha[q] = t_w ^ flag
                                break
                            q += 1
                        j ^= flag
                        r = j
                        while True:
                            if r == j or self.wd_lvls[r-1] == k:
                                if r != t_w:
                                    alpha[alpha_length] = r^flag
                                    alpha_length += 1
                            r += 1
                            if r >= i: break
                        invariant += self.wd_degree(CG, self_wd_ents[i-1], alpha[m], k, ham_wts)
                        invariant += (i-j)
                    j = i
            m += 1
        return invariant

    def _clear(self, k):
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._clear(2)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
        """
        self.clear(k)

    cdef void clear(self, int k) noexcept:
        cdef int i, j = 0, nwords = self.nwords, ncols = self.ncols
        cdef int *wd_lvls = self.wd_lvls
        cdef int *col_lvls = self.col_lvls
        for i from 0 <= i < nwords:
            if wd_lvls[i] >= k:
                wd_lvls[i] += 1
            if wd_lvls[i] < k:
                self.wd_percolate(j, i)
                j = i + 1
        j = 0
        for i from 0 <= i < ncols:
            if col_lvls[i] >= k:
                col_lvls[i] += 1
            if col_lvls[i] < k:
                self.col_percolate(j, i)
                j = i + 1

    cpdef int cmp(self, PartitionStack other, BinaryCode CG) noexcept:
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0], [0,0,1,1,1,1,0,0],
            ....:                    [0,0,0,0,1,1,1,1], [1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(0, [[0,0],[1,0]], B)
            181
            sage: P._split_vertex(0, 1)
            0
            sage: P._refine(1, [[0,0]], B)
            290
            sage: P._split_vertex(1, 2)
            1
            sage: P._refine(2, [[0,1]], B)
            463
            sage: P._split_vertex(2, 3)
            2
            sage: P._refine(3, [[0,2]], B)
            1500
            sage: P._split_vertex(4, 4)
            4
            sage: P._refine(4, [[0,4]], B)
            1224
            sage: P._is_discrete(4)
            1
            sage: Q = PartitionStack(P)
            sage: Q._clear(4)
            sage: Q._split_vertex(5, 4)
            4
            sage: Q._refine(4, [[0,4]], B)
            1224
            sage: Q._is_discrete(4)
            1
            sage: Q.cmp(P, B)
            0
        """
        cdef int *self_wd_ents = self.wd_ents
        cdef codeword *CG_words = CG.words
        cdef int i, j, l, m, span = 1, ncols = self.ncols, nwords = self.nwords
        for i from 0 < i < nwords:
            for j from 0 <= j < ncols:
                l = CG.is_one(self.wd_ents[i], self.col_ents[j])
                m = CG.is_one(other.wd_ents[i], other.col_ents[j])
                if l != m:
                    return l - m
        return 0

    def print_basis(self):
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 8)
            sage: P._dangerous_dont_use_set_ents_lvls(
            ....:     list(range(8)),
            ....:     list(range(7)) + [-1],
            ....:     [4,7,12,11,1,9,3,0,2,5,6,8,10,13,14,15],
            ....:     [0]*16)
            sage: P
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1,2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6},{7})
            sage: P._find_basis()
            sage: P.print_basis()
            basis_locations:
            4
            8
            0
            11
        """
        cdef int i, j
        if self.basis_locations:
            print("basis_locations:")
            j = 1
            while (1 << j) < self.nwords:
                j += 1
            for i from 0 <= i < j:
                print(self.basis_locations[i])

    def _find_basis(self):
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 8)
            sage: P._dangerous_dont_use_set_ents_lvls(list(range(8)), list(range(7))+[-1], [4,7,12,11,1,9,3,0,2,5,6,8,10,13,14,15], [0]*16)
            sage: P
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1,2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6},{7})
            sage: P._find_basis()
            sage: P.print_basis()
            basis_locations:
            4
            8
            0
            11
        """
        cdef int i
        cdef int *ham_wts = hamming_weights()
        self.find_basis(ham_wts)
        sig_free(ham_wts)

    cdef int find_basis(self, int *ham_wts) noexcept:
        cdef int i = 0, j, k, nwords = self.nwords, weight, basis_elts = 0, nrows = self.nrows
        cdef int *self_wd_ents = self.wd_ents
        if self.basis_locations is NULL:
            self.basis_locations = <int *> sig_malloc( 2 * nrows * sizeof(int) )
        if self.basis_locations is NULL:
            raise MemoryError("Memory.")
        while i < nwords:
            j = self_wd_ents[i]
            weight = ham_wts[j & 65535] + ham_wts[(j>>16) & 65535]
            if weight == 1:
                basis_elts += 1
                k = 0
                while not (1<<k) & j:
                    k += 1
                self.basis_locations[k] = i
                if basis_elts == nrows: break
            i += 1
        for i from 0 <= i < nrows:
            self.basis_locations[nrows + i] = self_wd_ents[1 << i]

    def _get_permutation(self, other):
        """
        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(0, [[0,0],[1,0]], B)
            181
            sage: P._split_vertex(0, 1)
            0
            sage: P._refine(1, [[0,0]], B)
            290
            sage: P._split_vertex(1, 2)
            1
            sage: P._refine(2, [[0,1]], B)
            463
            sage: P._split_vertex(2, 3)
            2
            sage: P._refine(3, [[0,2]], B)
            1500
            sage: P._split_vertex(4, 4)
            4
            sage: P._refine(4, [[0,4]], B)
            1224
            sage: P._is_discrete(4)
            1
            sage: Q = PartitionStack(P)
            sage: Q._clear(4)
            sage: Q._split_vertex(5, 4)
            4
            sage: Q._refine(4, [[0,4]], B)
            1224
            sage: Q._is_discrete(4)
            1
            sage: P._get_permutation(Q)
            ([0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 8, 9, 10, 11], [0, 1, 2, 3, 5, 4, 7, 6])
        """
        cdef int i
        cdef MemoryAllocator loc_mem = MemoryAllocator()
        cdef int *word_g = <int *> loc_mem.malloc(self.nwords * sizeof(int))
        cdef int *col_g = <int *> loc_mem.malloc(self.ncols * sizeof(int))
        self.get_permutation(other, word_g, col_g)
        word_l = [word_g[i] for i from 0 <= i < self.nwords]
        col_l = [col_g[i] for i from 0 <= i < self.ncols]
        return word_l, col_l

    cdef void get_permutation(self, PartitionStack other, int *word_gamma, int *col_gamma) noexcept:
        cdef int i
        cdef int *self_wd_ents = self.wd_ents
        cdef int *other_wd_ents = other.wd_ents
        cdef int *self_col_ents = self.col_ents
        cdef int *other_col_ents = other.col_ents
        # word_gamma[i] := image of the ith row as linear comb of rows
        for i from 0 <= i < self.nwords:
            word_gamma[other_wd_ents[i]] = self_wd_ents[i]
        for i from 0 <= i < self.ncols:
            col_gamma[other_col_ents[i]] = self_col_ents[i]

cdef class BinaryCodeClassifier:

    def __cinit__(self):
        """
        Initialize.

        TESTS::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: BC = BinaryCodeClassifier()
            sage: TestSuite(BC).run(skip='_test_pickling')
        """
        self.radix = sizeof(codeword) << 3
        self.ham_wts = hamming_weights()
        self.L = 100 # memory limit for Phi and Omega- multiply by 8KB
        self.aut_gens_size = self.radix * 100

        self.w_gamma_size = 1 << (self.radix/2)
        self.alpha_size = self.w_gamma_size + self.radix
        self.Phi_size = self.w_gamma_size/self.radix + 1

        self.mem = MemoryAllocator()
        self.w_gamma = <int *> self.mem.malloc(self.w_gamma_size * sizeof(int))
        self.alpha = <int *> self.mem.malloc(self.alpha_size * sizeof(int))
        self.Phi = <unsigned int *> self.mem.malloc(self.Phi_size * (self.L+1) * sizeof(unsigned int))
        self.Omega = <unsigned int *> self.mem.malloc(self.Phi_size * self.L * sizeof(unsigned int))
        self.W = <unsigned int *> self.mem.malloc(self.Phi_size * self.radix * 2 * sizeof(unsigned int))

        self.base = <int *> self.mem.malloc(self.radix * sizeof(int))
        self.aut_gp_gens = <int *> self.mem.malloc(self.aut_gens_size  * sizeof(int))
        self.c_gamma = <int *> self.mem.malloc(self.radix * sizeof(int))
        self.labeling = <int *> self.mem.malloc(self.radix * 3 * sizeof(int))
        self.Lambda1 = <int *> self.mem.malloc(self.radix * 2 * sizeof(int))
        self.Lambda2 = <int *> self.mem.malloc(self.radix * 2 * sizeof(int))
        self.Lambda3 = <int *> self.mem.malloc(self.radix * 2 * sizeof(int))
        self.v = <int *> self.mem.malloc(self.radix * 2 * sizeof(int))
        self.e = <int *> self.mem.malloc(self.radix * 2 * sizeof(int))

    cdef void record_automorphism(self, int *gamma, int ncols) noexcept:
        cdef int i, j
        if self.aut_gp_index + ncols > self.aut_gens_size:
            self.aut_gens_size *= 2
            self.aut_gp_gens = <int *> self.mem.realloc(self.aut_gp_gens, self.aut_gens_size * sizeof(int))
        j = self.aut_gp_index
        for i from 0 <= i < ncols:
            self.aut_gp_gens[i+j] = gamma[i]
        self.aut_gp_index += ncols

    def _aut_gp_and_can_label(self, CC, verbosity=0):
        """
        Compute the automorphism group and canonical label of the code ``CC``.

        INPUT:

        - ``CC`` -- a BinaryCode object
        - ``verbosity`` -- nonnegative integer

        OUTPUT:

        a tuple, (gens, labeling, size, base)
        - gens; list of permutations (in list form) representing generators
          of the permutation automorphism group of the code CC
        - labeling; a permutation representing the canonical labeling of the
          code. mostly for internal use; entries describe the relabeling
          on the columns.
        - size; the order of the automorphism group
        - base; a set of cols whose action determines the action on all cols

        EXAMPLES::

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: BC = BinaryCodeClassifier()

            sage: M = Matrix(GF(2),[
            ....:  [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
            ....:  [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
            ....:  [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size, base = BC._aut_gp_and_can_label(B)

            sage: # needs sage.groups
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            322560
            sage: size
            322560

            sage: M = Matrix(GF(2),[
            ....:  [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0],
            ....:  [0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1],
            ....:  [0,0,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size, base = BC._aut_gp_and_can_label(B)

            sage: # needs sage.groups
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            2304
            sage: size
            2304

            sage: M = Matrix(GF(2),[
            ....:  [1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0],
            ....:  [0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0],
            ....:  [0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0],
            ....:  [0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0],
            ....:  [0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0],
            ....:  [0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0],
            ....:  [0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0],
            ....:  [0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size, base = BC._aut_gp_and_can_label(B)

            sage: # needs sage.groups
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            136
            sage: size
            136

            sage: M = Matrix(GF(2),[
            ....:  [0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1],
            ....:  [1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,0],
            ....:  [0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0],
            ....:  [1,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1],
            ....:  [1,1,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,0],
            ....:  [1,0,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0],
            ....:  [0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,0,0,0],
            ....:  [0,0,1,0,0,1,0,1,1,1,0,1,1,0,1,0,0,1,0,0,0,1],
            ....:  [0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1],
            ....:  [1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,1],
            ....:  [0,0,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0,1,1,1,0]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size, base = BC._aut_gp_and_can_label(B)

            sage: # needs sage.groups
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            887040
            sage: size
            887040

            sage: B = BinaryCode(Matrix(GF(2),[[1,0,1],[0,1,1]]))
            sage: BC._aut_gp_and_can_label(B)
            ([[0, 2, 1], [1, 0, 2]], [0, 1, 2], 6, [0, 1])

            sage: B = BinaryCode(Matrix(GF(2),[[1,1,1,1]]))
            sage: BC._aut_gp_and_can_label(B)
            ([[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]], [0, 1, 2, 3], 24, [0, 1, 2])

            sage: B = BinaryCode(Matrix(GF(2),[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]))
            sage: gens, labeling, size, base = BC._aut_gp_and_can_label(B)
            sage: size
            87178291200

            sage: M = Matrix(GF(2),[
            ....:  [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],
            ....:  [0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1],
            ....:  [0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1],
            ....:  [0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1],
            ....:  [0,1,0,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0]])
            sage: B = BinaryCode(M)
            sage: BC._aut_gp_and_can_label(B)[2]
            2160

            sage: M = Matrix(GF(2),[
            ....:  [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
            ....:  [1,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0],
            ....:  [1,1,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,0,0],
            ....:  [1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: BC._aut_gp_and_can_label(B)[2]
            294912

            sage: M = Matrix(GF(2), [
            ....:  [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
            ....:  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0],
            ....:  [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0],
            ....:  [0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0],
            ....:  [0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1],
            ....:  [0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,1],
            ....:  [0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1]])
            sage: B = BinaryCode(M)
            sage: BC = BinaryCodeClassifier()
            sage: BC._aut_gp_and_can_label(B)[2]
            442368
        """
        cdef int i, j
        cdef BinaryCode C = CC
        self.aut_gp_and_can_label(C, verbosity)
        i = 0
        py_aut_gp_gens = []
        while i < self.aut_gp_index:
            gen = [self.aut_gp_gens[i+j] for j from 0 <= j < C.ncols]
            py_aut_gp_gens.append(gen)
            i += C.ncols
        py_labeling = [self.labeling[i] for i from 0 <= i < C.ncols]
        base = []
        for i from 0 <= i < self.radix:
            if self.base[i] == -1:
                break
            base.append(self.base[i])
        aut_gp_size = self.aut_gp_size
        return py_aut_gp_gens, py_labeling, aut_gp_size, base

    cdef void aut_gp_and_can_label(self, BinaryCode C, int verbosity) noexcept:

        # declare variables:
        cdef int i, j, ii, jj, iii, jjj, iiii # local variables

        cdef PartitionStack nu, zeta, rho # nu is the current position in the tree,
                                          # zeta the first terminal position,
                                          # and rho the best-so-far guess at canonical labeling position
        cdef int k = 0  # the number of partitions in nu
        cdef int k_rho  # the number of partitions in rho
        cdef int *v = self.v    # list of vertices determining nu
        cdef int h = -1 # longest common ancestor of zeta and nu: zeta[h] == nu[h], zeta[h+1] != nu[h+1]
                        # -1 indicates that zeta is not yet defined
        cdef int hb     # longest common ancestor of rho and nu:
                        # rho[hb] == nu[hb], rho[hb+1] != nu[hb+1]
        cdef int hh = 1 # the height of the oldest ancestor of nu satisfying Lemma 2.25 in [1]:
                        # if nu does not satisfy it at k, then hh = k
        cdef int ht # smallest such that all descendants of zeta[ht] are equivalent under
                    # the portion of the automorphism group so far discovered
        cdef int *alpha # for storing pointers to cells of nu[k]
        cdef int tvc    # tvc keeps track of which vertex is the first where nu and zeta differ-
                        # zeta was defined by splitting one vertex, and nu was defined by splitting tvc

        cdef OrbitPartition Theta # keeps track of which vertices have been discovered to be equivalent
        cdef unsigned int *Phi      # Phi stores the fixed point sets of each automorphism
        cdef unsigned int *Omega    # Omega stores the minimal elements of each cell of the orbit partition
        cdef int l = -1    # current index for storing values in Phi and Omega- we start at -1 so that when
                           # we increment first, the first place we write to is 0.
        cdef unsigned int *W    # for each k, W[k] is a list (as int mask) of the vertices to be searched down from
                       # the current partition, at k. Phi and Omega are ultimately used to make the size of
                       # W as small as possible
        cdef int *e  # 0 or 1, whether or not we have used Omega and Phi to narrow down W[k] yet: see states 12 and 17

        cdef int index = 0 # Define $\Gamma^{(-1)} := \text{Aut}(C)$, and
                           # $\Gamma^{(i)} := \Gamma^{(-1)}_{v_0,...,v_i}$.
                           # Then index = $|\Gamma^{(k-1)}|/|\Gamma^{(k)}|$ at (POINT A)
                           # and size = $|\Gamma^{(k-1)}|$ at (POINT A) and (POINT B).

        cdef int *Lambda = self.Lambda1             # for tracking indicator values- zf and zb are
        cdef int *zf__Lambda_zeta = self.Lambda2    # indicator vectors remembering Lambda[k] for
        cdef int *zb__Lambda_rho = self.Lambda3     # zeta and rho, respectively
        cdef int qzb              # keeps track of Lambda[k] {>,<,=} zb[k]
        cdef int hzf__h_zeta      # the max height for which Lambda and zf agree
        cdef int hzb__h_rho = -1  # the max height for which Lambda and zb agree

        cdef int *word_gamma
        cdef int *col_gamma = self.c_gamma # used for storing permutations
        cdef int nwords = C.nwords, ncols = C.ncols, nrows = C.nrows
        cdef int *ham_wts = self.ham_wts
        cdef int state  # keeps track of position in algorithm - see sage/graphs/graph_isom.pyx, search for "STATE DIAGRAM"

        self.aut_gp_index = 0
        self.aut_gp_size = Integer(1)

        if self.w_gamma_size < nwords:
            while self.w_gamma_size < nwords:
                self.w_gamma_size *= 2
            self.alpha_size = self.w_gamma_size + self.radix
            self.Phi_size = self.w_gamma_size/self.radix + 1
            self.w_gamma = <int *> self.mem.realloc(self.w_gamma, self.w_gamma_size * sizeof(int))
            self.alpha = <int *> self.mem.realloc(self.alpha, self.alpha_size * sizeof(int))
            self.Phi = <unsigned int *> self.mem.realloc(self.Phi, self.Phi_size * self.L * sizeof(int))
            self.Omega = <unsigned int *> self.mem.realloc(self.Omega, self.Phi_size * self.L * sizeof(int))
            self.W = <unsigned int *> self.mem.realloc(self.W, self.Phi_size * self.radix * 2 * sizeof(int))

        for i from 0 <= i < self.Phi_size * self.L:
            self.Omega[i] = 0
        word_gamma = self.w_gamma
        alpha = self.alpha # think of alpha as of length exactly nwords + ncols
        Phi = self.Phi
        Omega = self.Omega
        W = self.W
        e = self.e
        nu = PartitionStack(nrows, ncols)
        Theta = OrbitPartition(nrows, ncols)

        # trivial case
        if ncols == 0 or nrows == 0:
            raise NotImplementedError("Must supply a nontrivial code.")

        state = 1
        while state != -1:

            if state == 1: # Entry point: once only
                alpha[0] = 0
                alpha[1] = nu.flag
                nu.refine(k, alpha, 2, C, ham_wts)
                if nu.sat_225(k): hh = k
                if nu.is_discrete(k):
                    state = 18
                    continue

                # store the first smallest nontrivial cell in W[k], and set v[k]
                # equal to its minimum element
                v[k] = nu.new_first_smallest_nontrivial(k, W, self.Phi_size * k)

                Lambda[k] = 0
                e[k] = 0
                state = 2

            elif state == 2: # Move down the search tree one level by refining nu:
                             # split out a vertex, and refine nu against it
                k += 1
                nu.clear(k)

                alpha[0] = nu.split_vertex(v[k-1], k)
                Lambda[k] = nu.refine(k, alpha, 1, C, ham_wts) # store the invariant to Lambda[k]
                # only if this is the first time moving down the search tree:
                if h == -1:
                    state = 5
                    continue

                # update hzf__h_zeta
                if hzf__h_zeta == k-1 and Lambda[k] == zf__Lambda_zeta[k]: hzf__h_zeta = k
                # update qzb
                if qzb == 0:
                    if zb__Lambda_rho[k] == -1 or Lambda[k] < zb__Lambda_rho[k]:
                        qzb = -1
                    elif Lambda[k] > zb__Lambda_rho[k]:
                        qzb = 1
                    else:
                        qzb = 0
                # update hzb
                if hzb__h_rho == k-1 and qzb == 0: hzb__h_rho = k
                # if Lambda[k] > zb[k], then zb[k] := Lambda[k]
                # (zb keeps track of the indicator invariants corresponding to
                # rho, the closest canonical leaf so far seen- if Lambda is
                # bigger, then rho must be about to change
                if qzb > 0: zb__Lambda_rho[k] = Lambda[k]
                state = 3

            elif state == 3:  # attempt to rule out automorphisms while moving down the tree
                # if k > hzf, then we know that nu currently does not look like zeta, the first
                # terminal node encountered, thus there is no automorphism to discover. If qzb < 0,
                # i.e. Lambda[k] < zb[k], then the indicator is not maximal, and we can't reach a
                # canonical leaf. If neither of these is the case, then proceed to state 4.
                if hzf__h_zeta <= k or qzb >= 0:
                    state = 4
                else:
                    state = 6

            elif state == 4:  # at this point we have -not- ruled out the presence of automorphisms
                if nu.is_discrete(k):
                    state = 7
                    continue  # we have a terminal node, so process it

                # otherwise, prepare to split out another column:
                # store the first smallest nontrivial cell in W[k], and set v[k]
                # equal to its minimum element
                v[k] = nu.new_first_smallest_nontrivial(k, W, self.Phi_size * k)
                if not nu.sat_225(k): hh = k + 1
                e[k] = 0 # see state 12 and 17
                state = 2 # continue down the tree

            elif state == 5:
                # same as state 3, but in the case where we haven't yet defined zeta
                # i.e. this is our first time down the tree. Once we get to the bottom,
                # we will have zeta = nu = rho, so we do:
                zf__Lambda_zeta[k] = Lambda[k]
                zb__Lambda_rho[k] = Lambda[k]
                state = 4

            elif state == 6:
                # at this stage, there is no reason to continue downward, so backtrack
                j = k

                # return to the longest ancestor nu[i] of nu that could have a
                # descendant equivalent to zeta or could improve on rho.
                # All terminal nodes descending from nu[hh] are known to be
                # equivalent, so i < hh. Also, if i > hzb, none of the
                # descendants of nu[i] can improve rho, since the indicator is
                # off (Lambda(nu) < Lambda(rho)). If i >= ht, then no descendant
                # of nu[i] is equivalent to zeta (see [1, p67]).
                if ht-1 > hzb__h_rho:
                    if ht-1 < hh-1:
                        k = ht-1
                    else:
                        k = hh-1
                else:
                    if hzb__h_rho < hh-1:
                        k = hzb__h_rho
                    else:
                        k = hh-1
                # TODO: is the following line necessary?
                if k == -1:
                    k = 0

                if hb > k:  # update hb since we are backtracking
                    hb = k
                # if j == hh, then all nodes lower than our current position are equivalent, so bail out
                if j == hh:
                    state = 13
                    continue

                # recall hh: the height of the oldest ancestor of zeta for which Lemma 2.25 is
                # satisfied, which implies that all terminal nodes descended from there are equivalent.
                # If we are looking at such a node, then the partition at nu[hh] can be used for later
                # pruning, so we store its fixed set and a set of representatives of its cells.
                if l < self.L-1: l += 1
                nu.new_min_cell_reps(hh, Omega, self.Phi_size*l)
                nu.fixed_vertices(hh, Phi, Omega, self.Phi_size*l)

                state = 12

            elif state == 7: # we have just arrived at a terminal node of the search tree T(G, Pi)
                # if this is the first terminal node, go directly to 18, to
                # process zeta
                if h == -1:
                    state = 18
                    continue

                # hzf is the extremal height of ancestors of both nu and zeta, so if k < hzf, nu is not
                # equivalent to zeta, i.e. there is no automorphism to discover.
                if k < hzf__h_zeta:
                    state = 8
                    continue

                nu.get_permutation(zeta, word_gamma, col_gamma)

                # if C^gamma == C, the permutation is an automorphism, goto 10
                if C.is_automorphism(col_gamma, word_gamma):
                    state = 10
                else:
                    state = 8

            elif state == 8: # we have just ruled out the presence of automorphism and have not yet
                             # considered whether nu improves on rho
                # if qzb < 0, then rho already has larger indicator tuple
                if qzb < 0:
                    state = 6
                    continue

                # if Lambda[k] > zb[k] or nu is shorter than rho, then we have an improvement for rho
                if (qzb > 0) or (k < k_rho):
                    state = 9
                    continue

                # now Lambda[k] == zb[k] and k == k_rho, so we appeal to an enumeration:
                j = nu.cmp(rho, C)
                # if C(nu) > C(rho), we have a new label, goto 9
                if j > 0:
                    state = 9
                    continue

                # if C(nu) < C(rho), no new label, goto 6
                if j < 0:
                    state = 6
                    continue

                # if C(nu) == C(rho), get the automorphism and goto 10
                rho.get_permutation(nu, word_gamma, col_gamma)

                state = 10

            elif state == 9: # nu is a better guess at the canonical label than rho
                rho = PartitionStack(nu)
                k_rho = k
                qzb = 0
                hb = k
                hzb__h_rho = k
                # set zb[k+1] = Infinity
                zb__Lambda_rho[k+1] = -1
                state = 6

            elif state == 10: # we have an automorphism to process
                # increment l
                if l < self.L-1: l += 1
                # store information about the automorphism to Omega and Phi
                ii = self.Phi_size*l
                jj = 1 + nwords/self.radix
#                Omega[ii] = ~(~0 << ncols)
                for i from 0 <= i < jj:
                    Omega[ii+i] = ~0
                    Phi[ii+i] = 0
                if nwords%self.radix:
                    jj += 1
#                Omega[ii+jj-1] = ~((1 << nwords%self.radix) - 1)
                # Omega stores the minimum cell representatives
                i = 0
                while i < ncols:
                    j = col_gamma[i]         # i is a minimum
                    while j != i:            # cell rep,
                        Omega[ii] ^= (1<<j)  # so cancel
                        j = col_gamma[j]     # cellmates
                    i += 1
                    while i < ncols and not Omega[ii]&(1<<i): # find minimal element
                        i += 1                                # of next cell
                i = 0
                jj = self.radix
                while i < nwords:
                    j = word_gamma[i]
                    while j != i:
                        Omega[ii+1+j/jj] ^= (1<<(j%jj))
                        j = word_gamma[j]
                    i += 1
                    while i < nwords and not Omega[ii+1+i/jj]&(1<<(i%jj)):
                        i += 1
                # Phi stores the columns fixed by the automorphism
                for i from 0 <= i < ncols:
                    if col_gamma[i] == i:
                        Phi[ii] ^= (1 << i)
                for i from 0 <= i < nwords:
                    if word_gamma[i] == i:
                        Phi[ii+1+i/jj] ^= (1<<(i%jj))

                # Now incorporate the automorphism into Theta
                j = Theta.merge_perm(col_gamma, word_gamma)

                # j stores whether anything happened or not- if not, then the automorphism we have
                # discovered is already in the subgroup spanned by the generators we have output
                if not j:
                    state = 11
                    continue

                # otherwise, we have a new generator, so record it:
                self.record_automorphism(col_gamma, ncols)
                # The variable tvc was set to be the minimum element of W[k] the last time the
                # algorithm came up to meet zeta. At this point, we were considering the new
                # possibilities for descending away from zeta at this level.
                # if this is still a minimum cell representative of Theta, even in light
                # of this new automorphism, then the current branch off of zeta hasn't been
                # found equivalent to one already searched yet, so there may still be a
                # better canonical label downward.
                if tvc & nu.flag:
                    i = tvc^nu.flag
                    if Theta.wd_min_cell_rep[Theta.wd_find(i)] == i:
                        state = 11
                        continue
                elif Theta.col_min_cell_rep[Theta.col_find(tvc)] == tvc:
                    state = 11
                    continue

                # Otherwise, proceed to where zeta meets nu:
                k = h
                state = 13

            elif state == 11: # We have just found a new automorphism, and deduced that there may
                # be a better canonical label below the current branch off of zeta. So go to where
                # nu meets rho
                k = hb
                state = 12

            elif state == 12: # Coming here from either state 6 or 11, the algorithm has discovered
                              # some new information. 11 came from 10, where a new line in Omega and
                              # Phi was just recorded, and 6 stored information about implicit auto-
                              # morphisms in Omega and Phi
                if e[k] == 1:
                    # this means that the algorithm has come upward to this position (in state 17)
                    # before, so we have already intersected W[k] with the bulk of Omega and Phi, but
                    # we should still catch up with the latest ones
                    ii = self.Phi_size*l
                    jj = self.Phi_size*k
                    j = 1 + nwords/self.radix
                    if nwords%self.radix:
                        j += 1
                    W[jj] &= Omega[ii]
                    for i from 0 < i < j:
                        W[jj+i] &= Omega[ii+i]
                state = 13

            elif state == 13: # hub state
                if k == -1:
                    state = -1
                    continue  # exit point

                if k > h:
                    state = 17
                    continue # we are still on the same principal branch from zeta

                if k == h:
                    state = 14
                    continue # update the stabilizer index and check for new splits,
                                                # since we have returned to a partition of zeta
                # otherwise k < h, hence we have just backtracked up zeta, and are one level closer to done
                h = k
                tvc = 0
                jj = self.Phi_size*k
                if W[jj]:
                    while not (1 << tvc) & W[jj]:
                        tvc += 1
                else:
                    ii = 0
                    while not W[jj+1+ii]:
                        ii += 1
                    while not W[jj+1+ii] & (1 << tvc):
                        tvc += 1
                    tvc = (ii*self.radix + tvc) ^ nu.flag
                # now tvc points to the minimal cell representative of W[k]
                state = 14

            elif state == 14: # see if there are any more splits to make from this level of zeta (see state 17)
                if v[k]&nu.flag == tvc&nu.flag:
                    if tvc&nu.flag:
                        if Theta.wd_find(v[k]^nu.flag) == Theta.wd_find(tvc^nu.flag):
                            index += 1
                    else:
                        if Theta.col_find(v[k]) == Theta.col_find(tvc):
                            index += 1

                            # keep tabs on how many elements are in the same cell of Theta as tvc
                # find the next split
                jj = self.Phi_size*k
                if v[k]&nu.flag:
                    ii = self.radix
                    i = (v[k]^nu.flag) + 1
                    while i < nwords and not (1 << i%ii) & W[jj+1+i/ii]:
                        i += 1
                    if i < nwords:
                        v[k] = i^nu.flag
                    else:
                        # there is no new split at this level
                        state = 16
                        continue
                    # new split column better be a minimal representative in Theta, or wasted effort
                    if Theta.wd_min_cell_rep[Theta.wd_find(i)] == i:
                        state = 15
                    else:
                        state = 14
                else:
                    i = v[k] + 1
                    while i < ncols and not (1 << i) & W[jj]:
                        i += 1
                    if i < ncols:
                        v[k] = i
                    else:
                        # there is no new split at this level
                        state = 16
                        continue
                    # new split column better be a minimal representative in Theta, or wasted effort
                    if Theta.col_min_cell_rep[Theta.col_find(v[k])] == v[k]:
                        state = 15
                    else:
                        state = 14

            elif state == 15: # split out the column v[k]
                # hh is smallest such that nu[hh] satisfies Lemma 2.25. If it is larger than k+1,
                # it must be modified, since we are changing that part
                if k + 1 < hh:
                    hh = k + 1
                # hzf is maximal such that indicators line up for nu and zeta
                if k < hzf__h_zeta:
                    hzf__h_zeta = k
                # hzb is longest such that nu and rho have the same indicators
                if hzb__h_rho >= k:
                    hzb__h_rho = k
                    qzb = 0
                state = 2

            elif state == 16: # backtrack up zeta, updating information about stabilizer vector
                jj = self.Phi_size*k
                if W[jj]:
                    i = W[jj]
                    j = ham_wts[i & 65535] + ham_wts[(i >> 16) & 65535]
                else:
                    i = 0
                    j = 0
                    ii = self.radix
                    while i*ii < nwords:
                        iii = W[jj+1+i]
                        j += ham_wts[iii & 65535] + ham_wts[(iii >> 16) & 65535]
                        i += 1
                if j == index and ht == k + 1: ht = k
                self.aut_gp_size *= index
                # (POINT A)
                index = 0
                k -= 1
                if hb > k: # update hb since we are backtracking
                    hb = k
                state = 13

            elif state == 17: # see if there are any more splits to make from this level of nu (and not zeta)

                jjj = self.Phi_size*k
                if e[k] == 0: # now is the time to narrow down W[k] by Omega and Phi
                    # intersect W[k] with each Omega[i] such that v[0]...v[k-1] is in Phi[i]
                    jj = self.Phi_size*self.L
                    iii = nwords/self.radix
                    if nwords%self.radix:
                        iii += 1
                    for ii from 0 <= ii < iii:
                        Phi[jj+ii] = 0
                    for ii from 0 <= ii < k:
                        if v[ii]&nu.flag:
                            i = v[ii]^nu.flag
                            Phi[jj+1+i/self.radix] ^= (1 << i%self.radix)
                        else:
                            Phi[jj] ^= (1 << v[ii])
                    for i from 0 <= i <= l:
                        ii = self.Phi_size*i
                        iiii = 1
                        for j from 0 <= j < iii:
                            if Phi[ii + j] & Phi[jj + j] != Phi[jj + j]:
                                iiii = 0
                                break
                        if iiii:
                            for j from 0 <= j < iii:
                                W[jjj + j] &= Omega[ii + j]
                e[k] = 1

                # see if there is a vertex to split out
                if nu.flag&v[k]:
                    i = (v[k]^nu.flag)
                    while i < nwords:
                        i += 1
                        if (1 << i%self.radix) & W[jjj+1+i/self.radix]: break
                    if i < nwords:
                        v[k] = i^nu.flag
                        state = 15
                        continue
                else:
                    i = v[k]
                    while i < ncols:
                        i += 1
                        if (1 << i) & W[jjj]: break
                    if i < ncols:
                        v[k] = i
                        state = 15
                        continue

                k -= 1
                state = 13

            elif state == 18: # the first time nu becomes a discrete partition: set up zeta, our "identity" leaf
                # initialize counters for zeta:
                h = k # zeta[h] == nu[h]
                ht = k # nodes descended from zeta[ht] are all equivalent
                hzf__h_zeta = k # max such that indicators for zeta and nu agree
                zeta = PartitionStack(nu)
                for i from 0 <= i < k:
                    self.base[i] = v[i]
                self.base_size = k
                if k != self.radix:
                    self.base[k] = -1
                # (POINT B)
                k -= 1
                rho = PartitionStack(nu)
                # initialize counters for rho:
                k_rho = k+1 # number of partitions in rho
                hzb__h_rho = k # max such that indicators for rho and nu agree - BDM had k+1
                hb = k # rho[hb] == nu[hb] - BDM had k+1
                qzb = 0 # Lambda[k] == zb[k], so...
                state = 13

        # end big while loop
        rho.find_basis(ham_wts)
        for i from 0 <= i < ncols:
            self.labeling[rho.col_ents[i]] = i
        for i from 0 <= i < 2*nrows:
            self.labeling[i+ncols] = rho.basis_locations[i]

    def put_in_canonical_form(self, BinaryCode B):
        """
        Puts the code into canonical form.

        Canonical form is obtained by performing row reduction, permuting the
        pivots to the front so that the generator matrix is of the form: the
        identity matrix augmented to the right by arbitrary data.

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: BC = BinaryCodeClassifier()
            sage: B = BinaryCode(codes.GolayCode(GF(2)).generator_matrix())
            sage: B.apply_permutation(list(range(24,-1,-1)))
            sage: B
            Binary [24,12] linear code, generator matrix
            [011000111010100000000000]
            [001001001111100000000001]
            [011010100101100000000010]
            [001101110001100000000100]
            [010011011001100000001000]
            [010110110011000000010000]
            [011101100110000000100000]
            [000011110110100001000000]
            [000111101101000010000000]
            [001111011010000100000000]
            [010110001110101000000000]
            [011100011101010000000000]
            sage: BC.put_in_canonical_form(B)
            sage: B
            Binary [24,12] linear code, generator matrix
            [100000000000001100111001]
            [010000000000001010001111]
            [001000000000001111010010]
            [000100000000010110101010]
            [000010000000010110010101]
            [000001000000010001101101]
            [000000100000011000110110]
            [000000010000011111001001]
            [000000001000010101110011]
            [000000000100010011011110]
            [000000000010001011110101]
            [000000000001001101101110]
        """
        aut_gp_gens, labeling, size, base = self._aut_gp_and_can_label(B)
        B._apply_permutation_to_basis(labeling)
        B.put_in_std_form()

    def generate_children(self, BinaryCode B, int n, int d=2):
        """
        Use canonical augmentation to generate children of the code `B`.

        INPUT:

        - ``B`` -- a :class:`BinaryCode`

        - ``n`` -- limit on the degree of the code

        - ``d`` -- test whether new vector has weight divisible by `d`. If
          `d=4`, this ensures that all doubly-even canonically augmented
          children are generated.

        EXAMPLES::

            sage: from sage.coding.binary_code import *
            sage: BC = BinaryCodeClassifier()
            sage: B = BinaryCode(Matrix(GF(2), [[1,1,1,1]]))
            sage: BC.generate_children(B, 6, 4)                                         # needs sage.groups
            [
            [1 1 1 1 0 0]
            [0 1 0 1 1 1]
            ]

        .. NOTE::

            The function ``codes.databases.self_orthogonal_binary_codes`` makes heavy
            use of this function.

        MORE EXAMPLES::

            sage: # needs sage.groups
            sage: soc_iter = codes.databases.self_orthogonal_binary_codes(12, 6, 4)
            sage: L = list(soc_iter)
            sage: for n in range(13):
            ....:   s = 'n=%2d : '%n
            ....:   for k in range(1,7):
            ....:       s += '%3d '%len([C for C in L
            ....:                        if C.length() == n and C.dimension() == k])
            ....:   print(s)
            n= 0 :   0   0   0   0   0   0
            n= 1 :   0   0   0   0   0   0
            n= 2 :   0   0   0   0   0   0
            n= 3 :   0   0   0   0   0   0
            n= 4 :   1   0   0   0   0   0
            n= 5 :   0   0   0   0   0   0
            n= 6 :   0   1   0   0   0   0
            n= 7 :   0   0   1   0   0   0
            n= 8 :   1   1   1   1   0   0
            n= 9 :   0   0   0   0   0   0
            n=10 :   0   1   1   1   0   0
            n=11 :   0   0   1   1   0   0
            n=12 :   1   2   3   4   2   0
        """
        cdef BinaryCode m
        cdef codeword *ortho_basis
        cdef codeword *B_can_lab
        cdef codeword current, swap
        cdef codeword word, temp, gate, nonzero_gate, orbit, bwd, k_gate
        cdef codeword *temp_basis
        cdef codeword *orbit_checks
        cdef codeword orb_chx_size, orb_chx_shift, radix_gate
        cdef WordPermutation *gwp
        cdef WordPermutation *hwp
        cdef WordPermutation *can_lab
        cdef WordPermutation *can_lab_inv
        cdef WordPermutation **parent_generators
        cdef BinaryCode B_aug
        cdef int i, ii, j, jj, ij, k = 0, parity, combo, num_gens
        cdef int base_size, row
        cdef int *multimod2_index
        cdef int *ham_wts = self.ham_wts
        cdef int *num_inner_gens
        cdef int *num_outer_gens
        cdef int *v
        cdef int log_2_radix
        cdef bint bingo, bingo2, bingo3

        B.put_in_std_form()
        ortho_basis = expand_to_ortho_basis(B, n) # modifies B!

        aut_gp_gens, labeling, size, base = self._aut_gp_and_can_label(B)
        B_can_lab = <codeword *> sig_malloc(B.nrows * sizeof(codeword))
        can_lab = create_word_perm(labeling[:B.ncols])
        if B_can_lab is NULL or can_lab is NULL:
            sig_free(ortho_basis)
            if B_can_lab is not NULL:
                sig_free(B_can_lab)
            if can_lab is not NULL:
                sig_free(can_lab)
            raise MemoryError()
        for i from 0 <= i < B.nrows:
            B_can_lab[i] = permute_word_by_wp(can_lab, B.basis[i])
        dealloc_word_perm(can_lab)
        row = 0
        current = 1
        while row < B.nrows:
            i = row
            while i < B.nrows and not B_can_lab[i] & current:
                i += 1
            if i < B.nrows:
                if i != row:
                    swap = B_can_lab[row]
                    B_can_lab[row] = B_can_lab[i]
                    B_can_lab[i] = swap
                for j from 0 <= j < row:
                    if B_can_lab[j] & current:
                        B_can_lab[j] ^= B_can_lab[row]
                for j from row < j < B.nrows:
                    if B_can_lab[j] & current:
                        B_can_lab[j] ^= B_can_lab[row]
                row += 1
            current = current << 1
        num_gens = len(aut_gp_gens)
        base_size = len(base)

        parent_generators = <WordPermutation **> sig_malloc( len(aut_gp_gens) * sizeof(WordPermutation*) )
        temp_basis = <codeword *> sig_malloc( self.radix * sizeof(codeword) )

        output = []

        for i from 0 <= i < len(aut_gp_gens):
            parent_generators[i] = create_word_perm(aut_gp_gens[i] + list(range(B.ncols, n)))

        word = 0
        while ortho_basis[k] & (((<codeword>1) << B.ncols) - 1):
            k += 1
        j = k
        while ortho_basis[j]:
            word ^= ortho_basis[j]
            j += 1

        log_2_radix = 0
        while ((<codeword>1) << log_2_radix) < <codeword>self.radix:
            log_2_radix += 1
        # now we assume (<codeword>1 << log_2_radix) == self.radix
        if k < log_2_radix:
            orb_chx_size = 0
        else:
            orb_chx_size = k - log_2_radix
        orbit_checks = <codeword *> sig_malloc( ((<codeword>1) << orb_chx_size) * sizeof(codeword) )
        if orbit_checks is NULL:
            raise MemoryError()
        for temp from 0 <= temp < ((<codeword>1) << orb_chx_size):
            orbit_checks[temp] = 0

        combo = 0
        parity = 0
        gate = (<codeword>1 << B.nrows) - 1
        k_gate = (<codeword>1 << k) - 1
        nonzero_gate = ( (<codeword>1 << (n-B.ncols)) - 1 ) << B.ncols
        radix_gate = (((<codeword>1) << log_2_radix) - 1)

        while True:
            if nonzero_gate & word == nonzero_gate and \
              (ham_wts[word & 65535] + ham_wts[(word >> 16) & 65535])%d == 0:
                temp = (word >> B.nrows) & ((<codeword>1 << k) - 1)
                if not orbit_checks[temp >> log_2_radix] & ((<codeword>1) << (temp & radix_gate)):
                    B_aug = BinaryCode(B, word)
                    aug_aut_gp_gens, aug_labeling, aug_size, aug_base = self._aut_gp_and_can_label(B_aug)

                    # check if (B, B_aug) ~ (m(B_aug), B_aug)

                    can_lab = create_word_perm(aug_labeling[:n])

                    can_lab_inv = create_inv_word_perm(can_lab)
                    for j from 0 <= j < B_aug.nrows:
                        temp_basis[j] = permute_word_by_wp(can_lab, B_aug.basis[j])

                    # row reduce to get canonical label
                    i = 0
                    j = 0
                    while j < B_aug.nrows:
                        ii = j
                        while ii < B_aug.nrows and not temp_basis[ii] & (<codeword>1 << i):
                            ii += 1
                        if ii != B_aug.nrows:
                            if ii != j:
                                swap = temp_basis[ii]
                                temp_basis[ii] = temp_basis[j]
                                temp_basis[j] = swap
                            for jj from 0 <= jj < j:
                                if temp_basis[jj] & (<codeword>1 << i):
                                    temp_basis[jj] ^= temp_basis[j]
                            for jj from j < jj < B_aug.nrows:
                                if temp_basis[jj] & (<codeword>1 << i):
                                    temp_basis[jj] ^= temp_basis[j]
                            j += 1
                        i += 1
                    # done row reduction

                    for j from 0 <= j < B.nrows:
                        temp_basis[j] = permute_word_by_wp(can_lab_inv, temp_basis[j])
                    from sage.matrix.constructor import matrix
                    from sage.rings.integer_ring import ZZ
                    from sage.groups.perm_gps.permgroup import PermutationGroup
                    from sage.groups.perm_gps.constructor import PermutationGroupElement
                    from sage.libs.gap.libgap import libgap
                    rs = []
                    for i from 0 <= i < B.nrows:
                        r = []
                        for j from 0 <= j < n:
                            r.append((((<codeword>1)<<j)&temp_basis[i])>>j)
                        rs.append(r)
                    m = BinaryCode(matrix(ZZ, rs))

                    m_aut_gp_gens, m_labeling, m_size, m_base = self._aut_gp_and_can_label(m)
                    if True:  # size*factorial(n-B.ncols) == m_size:

                        if len(m_aut_gp_gens) == 0:
                            aut_m = PermutationGroup([()])
                        else:
                            aut_m = PermutationGroup([PermutationGroupElement([a+1 for a in g]) for g in m_aut_gp_gens])

                        if len(aug_aut_gp_gens) == 0:
                            aut_B_aug = libgap(PermutationGroup([()]))
                        else:
                            aut_B_aug = libgap(PermutationGroup([PermutationGroupElement([a+1 for a in g]) for g in aug_aut_gp_gens]))
                        H = libgap(aut_m).Intersection2(aut_B_aug)
                        rt_transversal = [[int(a) - 1 for a in g.ListPerm(n)] for g in aut_B_aug.RightTransversal(H) if not g.IsOne()]
                        rt_transversal.append(list(range(n)))

                        bingo2 = 0
                        for coset_rep in rt_transversal:
                            hwp = create_word_perm(coset_rep)
                            # dealloc_word_perm(gwp)
                            bingo2 = 1
                            for j from 0 <= j < B.nrows:
                                temp = permute_word_by_wp(hwp, temp_basis[j])
                                if temp != B.words[temp & gate]:
                                    bingo2 = 0
                                    dealloc_word_perm(hwp)
                                    break
                            if bingo2:
                                dealloc_word_perm(hwp)
                                break
                        if bingo2:
                            from sage.matrix.constructor import matrix
                            from sage.rings.finite_rings.finite_field_constructor import GF
                            M = matrix(GF(2), B_aug.nrows, B_aug.ncols)
                            for i from 0 <= i < B_aug.ncols:
                                for j from 0 <= j < B_aug.nrows:
                                    M[j,i] = B_aug.is_one(1 << j, i)
                            output.append(M)
                    dealloc_word_perm(can_lab)
                    dealloc_word_perm(can_lab_inv)
                #...

                    orbits = [word]
                    j = 0
                    while j < len(orbits):
                        for i from 0 <= i < len(aut_gp_gens):
                            temp = <codeword> orbits[j]
                            temp = permute_word_by_wp(parent_generators[i], temp)
                            temp ^= B.words[temp & gate]
                            if temp not in orbits:
                                orbits.append(temp)
                        j += 1
                    for temp in orbits:
                        temp = (temp >> B.nrows) & k_gate
                        orbit_checks[temp >> log_2_radix] |= ((<codeword>1) << (temp & radix_gate))

            parity ^= 1
            i = 0
            if not parity:
                while not combo & (1 << i): i += 1
                i += 1
            if i == k: break
            else:
                combo ^= (1 << i)
                word ^= ortho_basis[i]

        for i from 0 <= i < len(aut_gp_gens):
            dealloc_word_perm(parent_generators[i])
        sig_free(B_can_lab)
        sig_free(parent_generators)
        sig_free(orbit_checks)
        sig_free(ortho_basis)
        sig_free(temp_basis)
        return output
