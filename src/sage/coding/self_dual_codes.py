# sage.doctest: needs sage.groups sage.modules
r"""
Enumerating binary self-dual codes

This module implements functions useful for studying binary self-dual codes.
The main function is ``self_dual_binary_codes``, which is a case-by-case list
of entries, each represented by a Python dictionary.

Format of each entry: a Python dictionary with keys ``'order autgp'``, ``'spectrum'``,
``'code'``, ``'Comment'``, ``'Type'``, where

- ``'code'`` -- a sd code `C` of length `n`, dim `n/2`, over `\GF{2}`

- ``'order autgp'`` -- order of the permutation automorphism group of `C`

- ``'Type'`` -- the type of `C` (which can be ``'I'`` or ``'II'``, in the binary case)

- ``'spectrum'`` -- the spectrum `[A_0,A_1,...,A_n]`

- ``'Comment'`` -- possibly an empty string

Python dictionaries were used since they seemed to be both
human-readable and allow others to update the database easiest.

- The following double ``for`` loop can be time-consuming but should
  be run once in a while for testing purposes. It should only print
  ``True`` and have no trace-back errors::

    sage: for n in [4,6,8,10,12,14,16,18,20,22]:                # not tested
    ....:    C = self_dual_binary_codes(n); m = len(C.keys())
    ....:    for i in range(m):
    ....:        C0 = C["%s"%n]["%s"%i]["code"]
    ....:        print([n,i,C["%s"%n]["%s"%i]["spectrum"] == C0.spectrum()])
    ....:        print(C0 == C0.dual_code())
    ....:        G = C0.automorphism_group_binary_code()
    ....:        print(C["%s" % n]["%s" % i]["order autgp"] == G.order())

- To check if the "Riemann hypothesis" holds, run the following
  code::

    sage: R = PolynomialRing(CC,"T")                            # not tested
    sage: T = R.gen()                                           # not tested
    sage: for n in [4,6,8,10,12,14,16,18,20,22]:                # not tested
    ....:      C = self_dual_binary_codes(n); m = len(C["%s"%n].keys())
    ....:      for i in range(m):
    ....:          C0 = C["%s"%n]["%s"%i]["code"]
    ....:          if C0.minimum_distance()>2:
    ....:              f = R(C0.sd_zeta_polynomial())
    ....:              print([n,i,[z[0].abs() for z in f.roots()]])

You should get lists of numbers equal to 0.707106781186548.

Here's a rather naive construction of self-dual codes in the binary
case:

For even `m`, let `A_m` denote the `m\times m` matrix over `\GF{2}`
given by adding the all 1s matrix to the identity matrix (in
``MatrixSpace(GF(2),m,m)`` of course). If `M_1, ..., M_r` are square
matrices, let `diag(M_1,M_2,...,M_r)` denote the "block diagonal"
matrix with the matrices `M_i` on the diagonal and 0s elsewhere. Let
`C(m_1,...,m_r,s)` denote the linear code with generator matrix
having block form `G = (I, A)`, where
`A = diag(A_{m_1},A_{m_2},...,A_{m_r},I_s)`, for some
(even) `m_i`'s and `s`, where
`m_1+m_2+...+m_r+s=n/2`. Note: Such codes
`C(m_1,...,m_r,s)` are SD.

SD codes not of this form will be called (for the purpose of
documenting the code below) "exceptional". Except when `n` is
"small", most sd codes are exceptional (based on a counting
argument and table 9.1 in the Huffman+Pless [HP2003]_, page 347).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AUTHORS:

- David Joyner (2007-08-11)

REFERENCES:

- [HP2003] \W. C. Huffman, V. Pless, Fundamentals of
  Error-Correcting Codes, Cambridge Univ. Press, 2003.

- [P] \V. Pless, *A classification of self-orthogonal codes over GF(2)*,
  Discrete Math 3 (1972) 209-246.
"""

from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.matrix.constructor import block_diagonal_matrix
from sage.rings.integer_ring import ZZ
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.cachefunc import cached_function

from sage.coding.linear_code import LinearCode

_F = GF(2)


def _MS(n):
    r"""
    For internal use; returns the floor(n/2) x n matrix space over GF(2).

    EXAMPLES::

        sage: import sage.coding.self_dual_codes as self_dual_codes
        sage: self_dual_codes._MS(2)
        Full MatrixSpace of 1 by 2 dense matrices over Finite Field of size 2
        sage: self_dual_codes._MS(3)
        Full MatrixSpace of 1 by 3 dense matrices over Finite Field of size 2
        sage: self_dual_codes._MS(8)
        Full MatrixSpace of 4 by 8 dense matrices over Finite Field of size 2
    """
    n2 = ZZ(n)/2
    return MatrixSpace(_F, n2, n)


def _matA(n):
    r"""
    For internal use; returns a list of square matrices over GF(2) `(a_{ij})`
    of sizes 0 x 0, 1 x 1, ..., n x n which are of the form
    `(a_{ij} = 1) + (a_{ij} = \delta_{ij})`.

    EXAMPLES::

        sage: import sage.coding.self_dual_codes as self_dual_codes
        sage: self_dual_codes._matA(4)
        [
                        [0 1 1]
                 [0 1]  [1 0 1]
        [], [0], [1 0], [1 1 0]
        ]
    """
    A = []
    n2 = n.quo_rem(2)[0]
    for j in range(n2+2):
        MS0 = MatrixSpace(_F, j, j)
        I = MS0.identity_matrix()
        O = MS0(j*j*[1])
        A.append(I+O)
    return A


def _matId(n):
    r"""
    For internal use; returns a list of identity matrices over GF(2)
    of sizes (floor(n/2)-j) x (floor(n/2)-j) for j = 0 ... (floor(n/2)-1).

    EXAMPLES::

        sage: import sage.coding.self_dual_codes as self_dual_codes
        sage: self_dual_codes._matId(6)
        [
        [1 0 0]
        [0 1 0]  [1 0]
        [0 0 1], [0 1], [1]
        ]
    """
    Id = []
    n2 = n.quo_rem(2)[0]
    for j in range(n2):
        MSn = MatrixSpace(_F, n2-j, n2-j)
        Id.append(MSn.identity_matrix())
    return Id


def _MS2(n):
    r"""
    For internal use; returns the floor(n/2) x floor(n/2) matrix space over GF(2).

    EXAMPLES::

        sage: import sage.coding.self_dual_codes as self_dual_codes
        sage: self_dual_codes._MS2(8)
        Full MatrixSpace of 4 by 4 dense matrices over Finite Field of size 2
    """
    n2 = n.quo_rem(2)[0]
    return MatrixSpace(_F, n2, n2)


def _I2(n):
    r"""
    Internal function.

    EXAMPLES::

        sage: from sage.coding.self_dual_codes import _I2
        sage: _I2(3)
        [1]
        sage: _I2(5)
        [1 0]
        [0 1]
        sage: _I2(7)
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    return _MS2(n).identity_matrix()


@cached_function
def _And7():
    """
    Auxiliary matrix And7.

    EXAMPLES::

        sage: from sage.coding.self_dual_codes import _And7
        sage: _And7()
        [1 1 1 0 0 1 1]
        [1 1 1 0 1 0 1]
        [1 1 1 0 1 1 0]
        [0 0 0 0 1 1 1]
        [0 1 1 1 0 0 0]
        [1 0 1 1 0 0 0]
        [1 1 0 1 0 0 0]
    """
    return matrix(_F, [[1, 1, 1, 0, 0, 1, 1],
                       [1, 1, 1, 0, 1, 0, 1],
                       [1, 1, 1, 0, 1, 1, 0],
                       [0, 0, 0, 0, 1, 1, 1],
                       [0, 1, 1, 1, 0, 0, 0],
                       [1, 0, 1, 1, 0, 0, 0],
                       [1, 1, 0, 1, 0, 0, 0]])


@cached_function
def _H8():
    """
    Auxiliary matrix H8.

    EXAMPLES::

        sage: from sage.coding.self_dual_codes import _H8
        sage: _H8()
        [ 1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1]
    """
    return matrix(ZZ, [[1, 1, 1, 1, 1, 1, 1, 1],
                       [1, -1, 1, -1, 1, -1, 1, -1],
                       [1, 1, -1, -1, 1, 1, -1, -1],
                       [1, -1, -1, 1, 1, -1, -1, 1],
                       [1, 1, 1, 1, -1, -1, -1, -1],
                       [1, -1, 1, -1, -1, 1, -1, 1],
                       [1, 1, -1, -1, -1, -1, 1, 1],
                       [1, -1, -1, 1, -1, 1, 1, -1]]) # from Guava's Hadamard matrices database

# Remark: The above matrix constructions aid in computing some "small" self-dual codes.

############## main functions ##############


def self_dual_binary_codes(n):
    r"""
    Return the dictionary of inequivalent binary self dual codes of length `n`.

    For `n=4` even, returns the sd codes of a given length, up to (perm)
    equivalence, the (perm) aut gp, and the type.

    The number of inequivalent "diagonal" sd binary codes in the database of
    length n is ("diagonal" is defined by the conjecture above) is the
    same as the restricted partition number of `n`, where only integers
    from the set 1, 4, 6, 8, ... are allowed. This is the coefficient of
    `x^n` in the series expansion
    `(1-x)^{-1}\prod_{2^\infty (1-x^{2j})^{-1}}`. Typing the
    command ``f = (1-x)(-1)\*prod([(1-x(2\*j))(-1) for j in range(2,18)])``
    into Sage, we obtain for the coeffs of `x^4`,
    `x^6`, ... [1, 1, 2, 2, 3, 3, 5, 5, 7, 7, 11, 11, 15, 15,
    22, 22, 30, 30, 42, 42, 56, 56, 77, 77, 101, 101, 135, 135, 176,
    176, 231] These numbers grow too slowly to account for all the sd
    codes (see Huffman+Pless' Table 9.1, referenced above). In fact, in
    Table 9.10 of [HP2003]_, the number `B_n` of inequivalent sd binary codes
    of length `n` is given::

        n   2 4 6 8 10 12 14 16 18 20 22 24  26  28  30
        B_n 1 1 1 2  2  3  4  7  9 16 25 55 103 261 731

    According to http://oeis.org/classic/A003179,
    the next 2 entries are: 3295, 24147.

    EXAMPLES::

        sage: C = codes.databases.self_dual_binary_codes(10)
        sage: C["10"]["0"]["code"] == C["10"]["0"]["code"].dual_code()
        True
        sage: C["10"]["1"]["code"] == C["10"]["1"]["code"].dual_code()
        True
        sage: len(C["10"].keys())  # number of inequiv sd codes of length 10
        2
        sage: C = codes.databases.self_dual_binary_codes(12)
        sage: C["12"]["0"]["code"] == C["12"]["0"]["code"].dual_code()
        True
        sage: C["12"]["1"]["code"] == C["12"]["1"]["code"].dual_code()
        True
        sage: C["12"]["2"]["code"] == C["12"]["2"]["code"].dual_code()
        True
    """
    self_dual_codes = {}

    if n == 4:
        # this code is Type I
        # [4,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup([ "(2,4)",  "(1,2)(3,4)" ])
        spectrum = [1, 0, 2, 0, 1]
        self_dual_codes_4_0 = {"order autgp":8,"code":LinearCode(genmat),"spectrum":spectrum,
                        "Type":"I","Comment":"Unique."}
        self_dual_codes["4"] = {"0":self_dual_codes_4_0}
        return self_dual_codes

    if n == 6:
        # this is Type I
        # [6,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(3,6)", "(2,3)(5,6)", "(1,2)(4,5)"] )
        spectrum = [1, 0, 3, 0, 3, 0, 1]
        self_dual_codes_6_0 = {"order autgp":48,"code":LinearCode(genmat),"spectrum":spectrum,
                "Type":"I","Comment":"Unique"}
        self_dual_codes["6"] = {"0":self_dual_codes_6_0}
        return self_dual_codes

    if n == 8:
        # the first code is Type I, the second is Type II
        # the second code is equiv to the extended Hamming [8,4,4] code.
        # [8,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(4,8)", "(3,4)(7,8)", "(2,3)(6,7)", "(1,2)(5,6)"] )
        spectrum = [1, 0, 4, 0, 6, 0, 4, 0, 1]
        self_dual_codes_8_0 = {"order autgp":384,"code":LinearCode(genmat),"spectrum":spectrum,
               "Type":"I","Comment":"Unique Type I of this length."}
        # [8,1]:
        genmat = _I2(n).augment(_matA(n)[4])
        # G = PermutationGroup( ["(4,5)(6,7)", "(4,6)(5,7)", "(3,4)(7,8)",\
        #                    "(2,3)(6,7)", "(1,2)(5,6)"] )
        spectrum = [1, 0, 0, 0, 14, 0, 0, 0, 1]
        self_dual_codes_8_1 = {"order autgp":1344,"code":LinearCode(genmat),"spectrum":spectrum,
                "Type":"II","Comment":"Unique Type II of this length."}
        self_dual_codes["8"] = {"0":self_dual_codes_8_0,"1":self_dual_codes_8_1}
        return self_dual_codes

    if n == 10:
        # Both of these are Type I; one has a unique lowest weight codeword
        # [10,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(5,10)", "(4,5)(9,10)", "(3,4)(8,9)",\
        #                       "(2,3)(7,8)", "(1,2)(6,7)"] )
        spectrum = [1, 0, 5, 0, 10, 0, 10, 0, 5, 0, 1]
        self_dual_codes_10_0 = {"order autgp":3840,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"No Type II of this length."}
        # [10,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        # G = PermutationGroup( ["(5,10)", "(4,6)(7,8)", "(4,7)(6,8)", "(3,4)(8,9)",\
        #                       "(2,3)(7,8)", "(1,2)(6,7)"] )
        spectrum = [1, 0, 1, 0, 14, 0, 14, 0, 1, 0, 1]
        self_dual_codes_10_1 = {"order autgp":2688,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Unique lowest weight nonzero codeword."}
        self_dual_codes["10"] = {"0":self_dual_codes_10_0,"1":self_dual_codes_10_1}
        return self_dual_codes

    if n == 12:
        # all of these are Type I
        # [12,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(6,12)", "(5,6)(11,12)", "(4,5)(10,11)", "(3,4)(9,10)",\
        #                       "(2,3)(8,9)", "(1,2)(7,8)"] )
        spectrum = [1, 0, 6, 0, 15, 0, 20, 0, 15, 0, 6, 0, 1]
        self_dual_codes_12_0 = {"order autgp":48080,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"No Type II of this length."}
        # [12,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        # G = PermutationGroup( ["(2,3)(4,7)", "(2,4)(3,7)", "(2,4,9)(3,7,8)", "(2,4,8,10)(3,9)",\
        #       "(1,2,4,7,8,10)(3,9)", "(2,4,8,10)(3,9)(6,12)", "(2,4,8,10)(3,9)(5,6,11,12)"] )
        spectrum = [1, 0, 2, 0, 15, 0, 28, 0, 15, 0, 2, 0, 1]
        self_dual_codes_12_1 = {"order autgp":10752,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Smallest automorphism group of these."}
        # [12,2]:
        genmat = _I2(n).augment(_matA(n)[6])
        # G = PermutationGroup( ["(5,6)(11,12)", "(5,11)(6,12)", "(4,5)(10,11)", "(3,4)(9,10)",\
        #                     "(2,3)(8,9)", "(1,2)(7,8)"] )
        spectrum = [1, 0, 0, 0, 15, 0, 32, 0, 15, 0, 0, 0, 1]
        self_dual_codes_12_2 = {"order autgp":23040,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Largest minimum distance of these."}
        self_dual_codes["12"] = {"0":self_dual_codes_12_0,"1":self_dual_codes_12_1,"2":self_dual_codes_12_2}
        return self_dual_codes

    if n == 14:
        # all of these are Type I; one has a unique lowest weight codeword
        # (there are 4 total inequiv sd codes of n = 14, by Table 9.10 [HP2003])
        # [14,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(7,14)", "(6,7)(13,14)", "(5,6)(12,13)", "(4,5)(11,12)",\
        #            "(3,4)(10,11)", "(2,3)(9,10)", "(1,2)(8,9)"] )
        spectrum = [1, 0, 7, 0, 21, 0, 35, 0, 35, 0, 21, 0, 7, 0, 1]
        self_dual_codes_14_0 = {"order autgp":645120,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"No Type II of this length. Huge aut gp."}
        # [14,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        # G = PermutationGroup( ["(7,14)", "(6,7)(13,14)", "(5,6)(12,13)", "(4,8)(9,10)",\
        #              "(4,9)(8,10)", "(3,4)(10,11)", "(2,3)(9,10)", "(1,2)(8,9)"] )
        spectrum = [1, 0, 3, 0, 17, 0, 43, 0, 43, 0, 17, 0, 3, 0, 1]
        self_dual_codes_14_1 = {"order autgp":64512,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Automorphism group has order 64512."}
        # [14,2]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matId(n)[6]]))
        # G = PermutationGroup( ["(7,14)", "(5,6)(12,13)", "(5,12)(6,13)", "(4,5)(11,12)",\
        #                        "(3,4)(10,11)", "(2,3)(9,10)", "(1,2)(8,9)"] )
        spectrum = [1, 0, 1, 0, 15, 0, 47, 0, 47, 0, 15, 0, 1, 0, 1]
        self_dual_codes_14_2 = {"order autgp":46080,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Unique codeword of weight 2."}
        # [14,3]:
        genmat = _I2(n).augment(_And7())
        # G = PermutationGroup( ["(7,11)(12,13)", "(7,12)(11,13)", "(6,9)(10,14)",\
        #      "(6,10)(9,14)", "(5,6)(8,9)", "(4,5)(9,10), (2,3)(11,12)", "(2,7)(3,13)",\
        #      "(1,2)(12,13)", "(1,4)(2,5)(3,8)(6,7)(9,13)(10,12)(11,14)"])
        spectrum = [1, 0, 0, 0, 14, 0, 49, 0, 49, 0, 14, 0, 0, 0, 1]
        self_dual_codes_14_3 = {"order autgp":56448,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Largest minimum distance of these."}
        self_dual_codes["14"] = {"0":self_dual_codes_14_0,"1":self_dual_codes_14_1,"2":self_dual_codes_14_2,
                  "3":self_dual_codes_14_3}
        return self_dual_codes

    if n == 16:
        # 4 of these are Type I, 2 are Type II. The 2 Type II codes
        # are formally equivalent but with different automorphism groups
        # [16,0]:
        genmat = _I2(n).augment(_I2(n))
        #  G = PermutationGroup( [ "(8,16)", "(7,8)(15,16)", "(6,7)(14,15)", "(5,6)(13,14)",
        #                       "(4,5)(12,13)", "(3,4)(11,12)", "(2,3)(10,11)", "(1,2)(9,10)"] )
        spectrum = [1, 0, 8, 0, 28, 0, 56, 0, 70, 0, 56, 0, 28, 0, 8, 0, 1]
        self_dual_codes_16_0 = {"order autgp":10321920,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Huge aut gp."}
        # [16,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        #  G = PermutationGroup( [ "(8,16)", "(7,8)(15,16)", "(6,7)(14,15)", "(5,6)(13,14)",\
        #        "(4,9)(10,11)", "(4,10)(9,11)", "(3,4)(11,12)", "(2,3)(10,11)", "(1,2)(9,10)"] )
        spectrum = [1, 0, 4, 0, 20, 0, 60, 0, 86, 0, 60, 0, 20, 0, 4, 0, 1]
        self_dual_codes_16_1 = {"order autgp":516096,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [16,2]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matA(n)[4]]))
        #  G = PermutationGroup( [ "(8,13)(14,15)", "(8,14)(13,15)", "(7,8)(15,16)", "(6,7)(14,15)",\
        #     "(5,6)(13,14)", "(4,9)(10,11)", "(4,10)(9,11)", "(3,4)(11,12)", "(2,3)(10,11)",\
        #     "(1,2)(9,10)","(1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)"] )
        spectrum = [1, 0, 0, 0, 28, 0, 0, 0, 198, 0, 0, 0, 28, 0, 0, 0, 1]
        self_dual_codes_16_2 = {"order autgp":3612672,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"II","Comment":"Same spectrum as the other Type II code."}
        # [16,3]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matId(n)[6]]))
        # G = PermutationGroup( [ "(8,16)", "(7,8)(15,16)", "(5,6)(13,14)", "(5,13)(6,14)",\
        #             "(4,5)(12,13)", "(3,4)(11,12)", "(2,3)(10,11)", "(1,2)(9,10)"] )
        spectrum = [1, 0, 2, 0, 16, 0, 62, 0, 94, 0, 62, 0, 16, 0, 2, 0, 1]
        self_dual_codes_16_3 = {"order autgp":184320,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [16,4]:
        genmat = _I2(n).augment(_matA(n)[8])
        # an equivalent form: See also [20,8] using A[10]
        # [(1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1),
        #  (0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1),
        #  (0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0),
        #  (0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0),
        #  (0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0),
        #  (0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0),
        #  (0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0),
        #  (0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1)]
        #  G = PermutationGroup( [ "(7,8)(15,16)", "(7,15)(8,16)", "(6,7)(14,15)",\
        #      "(5,6)(13,14)","(4,5)(12,13)","(3,4)(11,12)", "(2,3)(10,11)", "(1,2)(9,10)"] )
        spectrum = [1, 0, 0, 0, 28, 0, 0, 0, 198, 0, 0, 0, 28, 0, 0, 0, 1]
        self_dual_codes_16_4 = {"order autgp":5160960,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"II","Comment":"Same spectrum as the other Type II code. Large aut gp."}
        # [16,5]:
        genmat = _I2(n).augment(block_diagonal_matrix([_And7(),_matId(n)[7]]))
        #  G = PermutationGroup( [ "(8,16)", "(7,12)(13,14)", "(7,13)(12,14)",\
        #      "(6,10)(11,15)", "(6,11)(10,15)", "(5,6)(9,10)", "(4,5)(10,11)",\
        #      "(2,3)(12,13)", "(2,7)(3,14)", "(1,2)(13,14)",\
        #      "(1,4)(2,5)(3,9)(6,7)(10,14)(11,13)(12,15)" ] )
        spectrum = [1, 0, 1, 0, 14, 0, 63, 0, 98, 0, 63, 0, 14, 0, 1, 0, 1]
        self_dual_codes_16_5 = {"order autgp":112896,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"'Exceptional' construction."}
        # [16,6]:
        J8 = MatrixSpace(ZZ,8,8)(64*[1])
        genmat = _I2(n).augment(_I2(n)+_MS2(n)((_H8()+J8)/2))
        #  G = PermutationGroup( [ "(7,9)(10,16)", "(7,10)(9,16)", "(6,7)(10,11)",\
        #       "(4,6)(11,13)", "(3,5)(12,14)", "(3,12)(5,14)", "(2,3)(14,15)",\
        #       "(1,2)(8,15)", "(1,4)(2,6)(3,7)(5,16)(8,13)(9,12)(10,14)(11,15)" ] )
        spectrum = [1, 0, 0, 0, 12, 0, 64, 0, 102, 0, 64, 0, 12, 0, 0, 0, 1]
        self_dual_codes_16_6 = {"order autgp":73728,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"'Exceptional' construction. Min dist 4."}
        self_dual_codes["16"] = {"0":self_dual_codes_16_0,"1":self_dual_codes_16_1,"2":self_dual_codes_16_2,
                  "3":self_dual_codes_16_3,"4":self_dual_codes_16_4,"5":self_dual_codes_16_5,"6":self_dual_codes_16_6}
        return self_dual_codes

    if n == 18:
        # all of these are Type I, all are "extensions" of the n=16 codes
        # [18,3] and [18,4] each has a unique lowest weight codeword. Also, they
        # are formally equivalent but with different automorphism groups
        # [18,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( [ "(9,18)", "(8,9)(17,18)", "(7,8)(16,17)", "(6,7)(15,16)",\
        #     "(5,6)(14,15)", "(4,5)(13,14)", "(3,4)(12,13)", "(2,3)(11,12)", "(1,2)(10,11)" ] )
        spectrum = [1, 0, 9, 0, 36, 0, 84, 0, 126, 0, 126, 0, 84, 0, 36, 0, 9, 0, 1]
        self_dual_codes_18_0 = {"order autgp":185794560,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Huge aut gp. S_9x(ZZ/2ZZ)^9?"}
        # [18,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        #   G = PermutationGroup( [ "(9,18)", "(8,9)(17,18)", "(7,8)(16,17)", "(6,7)(15,16)",\
        #       "(5,6)(14,15)", "(4,10)(11,12)", "(4,11)(10,12)", "(3,4)(12,13)",\
        #       "(2,3)(11,12)", "(1,2)(10,11)" ] )
        spectrum = [1, 0, 5, 0, 24, 0, 80, 0, 146, 0, 146, 0, 80, 0, 24, 0, 5, 0, 1]
        self_dual_codes_18_1 = {"order autgp":5160960,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Large aut gp."}
        # [18,2]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matId(n)[6]]))
        #  G = PermutationGroup( [ "(9,18)", "(8,9)(17,18)", "(7,8)(16,17)", "(5,6)(14,15)",\
        #       "(5,14)(6,15)","(4,5)(13,14)", "(3,4)(12,13)", "(2,3)(11,12)", "(1,2)(10,11)"] )
        spectrum = [1, 0, 3, 0, 18, 0, 78, 0, 156, 0, 156, 0, 78, 0, 18, 0, 3, 0, 1]
        self_dual_codes_18_2 = {"order autgp":1105920,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": ""}
        # [18,3]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matA(n)[4],_matId(n)[8]]))
        #   G = PermutationGroup( [ "(9,18)", "(8,14)(15,16)", "(8,15)(14,16)", "(7,8)(16,17)",\
        #      "(6,7)(15,16)","(5,6)(14,15)", "(4,10)(11,12)", "(4,11)(10,12)",\
        #      "(3,4)(12,13)", "(2,3)(11,12)","(1,2)(10,11)",\
        #      "(1,5)(2,6)(3,7)(4,8)(10,14)(11,15)(12,16)(13,17)" ] )
        spectrum = [1, 0, 1, 0, 28, 0, 28, 0, 198, 0, 198, 0, 28, 0, 28, 0, 1, 0, 1]
        self_dual_codes_18_3 = {"order autgp":7225344,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Large aut gp. Unique codeword of smallest nonzero wt.\
                 Same spectrum as '[18,4]' sd code."}
        # [18,4]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[8],_matId(n)[8]]))
        # G = PermutationGroup( [ "(9,18)", "(7,8)(16,17)", "(7,16)(8,17)", "(6,7)(15,16)", \
        #     "(5,6)(14,15)", "(4,5)(13,14)", "(3,4)(12,13)", "(2,3)(11,12)", "(1,2)(10,11)" ] )
        spectrum = [1, 0, 1, 0, 28, 0, 28, 0, 198, 0, 198, 0, 28, 0, 28, 0, 1, 0, 1]
        self_dual_codes_18_4 = {"order autgp":10321920,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Huge aut gp. Unique codeword of smallest nonzero wt.\
                 Same spectrum as '[18,3]' sd code."}
        # [18,5]:
        C = self_dual_binary_codes(n-2)["%s" % (n-2)]["5"]["code"]
        A0 = C.redundancy_matrix()
        genmat = _I2(n).augment(block_diagonal_matrix([A0,_matId(n)[8]]))
        # G = PermutationGroup( [ "(5,10)(6,11)", "(5,11)(6,10)", "(5,11,12)(6,7,10)",\
        #     "(5,11,10,7,12,6,13)", "(2,15)(3,16)(5,11,10,7,12,6,13)",\
        #     "(2,16)(3,15)(5,11,10,7,12,6,13)", "(2,16,14)(3,15,4)(5,11,10,7,12,6,13)",\
        #     "(1,2,16,15,4,3,14)(5,11,10,7,12,6,13)", "(1,5,14,6,16,11,15,7,3,10,4,12,2,13)",\
        #     "(2,16,14)(3,15,4)(5,11,10,7,12,6,13)(9,18)",\
        #     "(2,16,14)(3,15,4)(5,11,10,7,12,6,13)(8,9,17,18)" ] )
        spectrum = [1, 0, 2, 0, 15, 0, 77, 0, 161, 0, 161, 0, 77, 0, 15, 0, 2, 0, 1]
        self_dual_codes_18_5 = {"order autgp":451584,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional' construction."}
        # [18,6]:
        C = self_dual_binary_codes(n-2)["%s" % (n-2)]["6"]["code"]
        A0 = C.redundancy_matrix()
        genmat = _I2(n).augment(block_diagonal_matrix([A0,_matId(n)[8]]))
        G = PermutationGroup( [ "(9,18)", "(7,10)(11,17)", "(7,11)(10,17)", "(6,7)(11,12)",
              "(4,6)(12,14)", "(3,5)(13,15)", "(3,13)(5,15)", "(2,3)(15,16)", "(1,2)(8,16)",
              "(1,4)(2,6)(3,7)(5,17)(8,14)(10,13)(11,15)(12,16)" ] )
        spectrum = [1, 0, 1, 0, 12, 0, 76, 0, 166, 0, 166, 0, 76, 0, 12, 0, 1, 0, 1]
        self_dual_codes_18_6 = {"order autgp":147456,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional'. Unique codeword of smallest nonzero wt."}
        # [18,7] (equiv to H18 in [P])
        genmat = _MS(n)([[1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0],
                         [0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1],
                         [0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1],
                         [0,0,0,1,0,0,0,0,0,1,1,1,1,0,0,0,0,1],
                         [0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0],
                         [0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0],
                         [0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0],
                         [0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1],
                         [0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1]])
        # G = PermutationGroup( [ "(9,10)(16,18)", "(9,16)(10,18)", "(8,9)(14,16)",\
        #          "(7,11)(12,17)", "(7,12)(11,17)", "(5,6)(11,12)", "(5,7)(6,17)",\
        #          "(4,13)(5,8)(6,14)(7,9)(10,12)(11,18)(16,17)", "(3,4)(13,15)",\
        #          "(1,2)(5,8)(6,14)(7,9)(10,12)(11,18)(16,17)", "(1,3)(2,15)",\
        #          "(1,5)(2,6)(3,7)(4,11)(10,18)(12,13)(15,17)" ] )
        spectrum = [1, 0, 0, 0, 9, 0, 75, 0, 171, 0, 171, 0, 75, 0, 9, 0, 0, 0, 1]
        self_dual_codes_18_7 = {"order autgp":82944,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional' construction. Min dist 4."}
        # [18, 8] (equiv to I18 in [P])
        I18 = _MS(n)([[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
                      [1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0],
                      [0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1],
                      [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]])
        genmat = _MS(n)([[1,0,0,0,0,0,0,0,0, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                         [0,1,0,0,0,0,0,0,0, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                         [0,0,1,0,0,0,0,0,0, 0, 1, 1, 0, 0, 0, 1, 1, 1],
                         [0,0,0,1,0,0,0,0,0, 0, 1, 0, 0, 1, 0, 1, 1, 1],
                         [0,0,0,0,1,0,0,0,0, 0, 1, 0, 1, 0, 0, 1, 1, 1],
                         [0,0,0,0,0,1,0,0,0, 1, 1, 0, 0, 0, 0, 1, 1, 1],
                         [0,0,0,0,0,0,1,0,0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
                         [0,0,0,0,0,0,0,1,0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
                         [0,0,0,0,0,0,0,0,1, 0, 0, 0, 0, 0, 1, 1, 1, 0]])
        G = PermutationGroup( [ "(9,15)(16,17)", "(9,16)(15,17)", "(8,9)(17,18)",
                       "(7,8)(16,17)", "(5,6)(10,13)", "(5,10)(6,13)", "(4,5)(13,14)",
                      "(3,4)(12,14)", "(1,2)(6,10)", "(1,3)(2,12)" ] )
        spectrum = [1, 0, 0, 0, 17, 0, 51, 0, 187, 0, 187, 0, 51, 0, 17, 0, 0, 0, 1]
        self_dual_codes_18_8 = {"order autgp":322560,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional' construction. Min dist 4."}
        self_dual_codes["18"] = {"0":self_dual_codes_18_0,"1":self_dual_codes_18_1,"2":self_dual_codes_18_2,
                  "3":self_dual_codes_18_3,"4":self_dual_codes_18_4,"5":self_dual_codes_18_5,
                  "6":self_dual_codes_18_6,"7":self_dual_codes_18_7,"8":self_dual_codes_18_8}
        return self_dual_codes

    if n == 20:
    # all of these of these are Type I; 2 of these codes
    # are formally equivalent but with different automorphism groups;
    # one of these has a unique codeword of lowest weight
        A10 = MatrixSpace(_F,10,10)([[1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                                     [1, 1, 1, 0, 1, 0, 1, 0, 1, 1],
                                     [1, 0, 0, 1, 0, 1, 0, 1, 0, 1],
                                     [0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
                                     [0, 0, 1, 1, 0, 1, 0, 1, 0, 1],
                                     [0, 0, 0, 1, 0, 1, 1, 1, 0, 1],
                                     [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                                     [0, 0, 0, 1, 0, 0, 0, 0, 1, 1],
                                     [0, 0, 0, 0, 0, 1, 0, 0, 1, 1],
                                     [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]])
        # [20,0]:
        genmat = _I2(n).augment(_I2(n))
        # G = PermutationGroup( ["(10,20)", "(9,10)(19,20)", "(8,9)(18,19)", "(7,8)(17,18)", "(6,7)(16,17)",\
        #            "(5,6)(15,16)", "(4,5)(14,15)", "(3,4)(13,14)", "(2,3)(12,13)", "(1,2)(11,12)"] )
        spectrum = [1, 0, 10, 0, 45, 0, 120, 0, 210, 0, 252, 0, 210, 0, 120, 0, 45, 0, 10, 0, 1]
        self_dual_codes_20_0 = {"order autgp":3715891200,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Huge aut gp"}
        # [20,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        # G = PermutationGroup( [ "(10,20)", "(9,10)(19,20)", "(8,9)(18,19)", "(7,8)(17,18)", "(6,7)(16,17)",\
        #         "(5,6)(15,16)", "(4,11)(12,13)", "(4,12)(11,13)", "(3,4)(13,14)",\
        #         "(2,3)(12,13)", "(1,2)(11,12)"] )
        spectrum = [1, 0, 6, 0, 29, 0, 104, 0, 226, 0, 292, 0, 226, 0, 104, 0, 29, 0, 6, 0, 1]
        self_dual_codes_20_1 = {"order autgp":61931520,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [20,2]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matId(n)[6]]))
        #  G = PermutationGroup( [ "(10,20)", "(9,10)(19,20)", "(8,9)(18,19)", "(7,8)(17,18)",\
        #          "(5,6)(15,16)", "(5,15)(6,16)", "(4,5)(14,15)", "(3,4)(13,14)",\
        #          "(2,3)(12,13)", "(1,2)(11,12)"] )
        spectrum = [1, 0, 4, 0, 21, 0, 96, 0, 234, 0, 312, 0, 234, 0, 96, 0, 21, 0, 4, 0, 1]
        self_dual_codes_20_2 = {"order autgp":8847360,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [20,3]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matA(n)[4]]))
        # G = PermutationGroup( [ "(5,6)(15,16)", "(5,15)(6,16)", "(4,5)(14,15)", "(3,4)(13,14)",\
        #             "(2,3)(12,13)", "(1,2)(11,12)", "(8,17)(9,10)", "(8,10)(9,17)", "(8,10,20)(9,19,17)",\
        #             "(8,19,20,9,17,10,18)", "(7,8,19,20,9,18)(10,17)"] )
        spectrum = [1, 0, 0, 0, 29, 0, 32, 0, 226, 0, 448, 0, 226, 0, 32, 0, 29, 0, 0, 0, 1]
        self_dual_codes_20_3 = {"order autgp":30965760,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Min dist 4."}
        # [20,4]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matA(n)[4],_matId(n)[8]]))
        #  G = PermutationGroup( [ "(5,15)(6,16)", "(5,16)(6,15)", "(5,16,7)(6,17,15)", "(5,15,8)(6,17,7)",\
        #              "(5,17,18)(6,15,8), (3,14)(4,13)(5,17,18)(6,15,8)", "(3,13)(4,14)(5,17,18)(6,15,8)",\
        #              "(2,3,14)(4,13,11)(5,17,18)(6,15,8)"," (2,3,12)(4,11,14)(5,17,18)(6,15,8)",\
        #              "(1,2,3,11,14,4,12)(5,17,18)(6,15,8)", "(1,5,13,17,14,8,2,7,3,16,12,6,11,18)(4,15)",\
        #               "(2,3,12)(4,11,14)(5,17,18)(6,15,8)(10,20)",\
        #               "(2,3,12)(4,11,14)(5,17,18)(6,15,8)(9,10,19,20)"] )
        spectrum = [1, 0, 2, 0, 29, 0, 56, 0, 226, 0, 396, 0, 226, 0, 56, 0, 29, 0, 2, 0, 1]
        self_dual_codes_20_4 = {"order autgp":28901376,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [20,5]:
        genmat = _I2(n).augment(block_diagonal_matrix([_And7(),_matId(n)[7]]))
        # G = PermutationGroup( [ "(10,20)", "(9,10)(19,20)", "(8,9)(18,19)",\
        #        "(7,11)(12,14)", "(7,12)(11,14)", "(6,7)(12,13)", "(5,6)(11,12)",\
        #       "(4,15)(16,17)", "(4,16)(15,17)", "(2,3)(16,17)", "(2,4)(3,15)",\
        #        "(1,2)(15,16)", "(1,5)(2,6)(3,13)(4,7)(11,16)(12,15)(14,17)" ] ) # order 2709504
        spectrum = [1, 0, 3, 0, 17, 0, 92, 0, 238, 0, 322, 0, 238, 0, 92, 0, 17, 0, 3, 0, 1]
        self_dual_codes_20_5 = {"order autgp":2709504,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional' construction."}
        # [20,6]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[8],_matId(n)[8]]))
        # G = PermutationGroup( [ "(7,8)(17,18)", "(7,17)(8,18)", "(6,7)(16,17)", "(5,6)(15,16)",\
        #        "(4,5)(14,15)", "(3,4)(13,14)", "(2,3)(12,13)", "(1,2)(11,12)",\
        #        "(10,20)", "(9,10,19,20)"] )
        spectrum = [1, 0, 2, 0, 29, 0, 56, 0, 226, 0, 396, 0, 226, 0, 56, 0, 29, 0, 2, 0, 1]
        self_dual_codes_20_6 = {"order autgp":41287680,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [20,7]:
        A0 = self_dual_binary_codes(n-4)["16"]["6"]["code"].redundancy_matrix()
        genmat = _I2(n).augment(block_diagonal_matrix([A0,_matId(n)[8]]))
        # G = PermutationGroup( [ "(10,20)", "(9,10)(19,20)", "(7,11)(12,18)",\
        #    "(7,12)(11,18)", "(6,7)(12,13)", "(4,6)(13,15)", "(3,5)(14,16)",\
        #    "(3,14)(5,16)", "(2,3)(16,17)", "(1,2)(8,17)",\
        #    "(1,4)(2,6)(3,7)(5,18)(8,15)(11,14)(12,16)(13,17)" ] )
        spectrum = [1,0,2,0,13,0,88,0,242,0,332,0,242,0,88,0,13,0,2,0,1]
        self_dual_codes_20_7 = {"order autgp":589824,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"'Exceptional' construction."}
        # [20,8]: (genmat, J20, and genmat2 are all equiv)
        genmat = _I2(n).augment(_matA(n)[10])
        J20 = _MS(n)([[1,1,1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0,0,1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0,0,0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0,0,0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0,0,0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0,0,0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                     [0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
                     [0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
                     [0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
                     [1,0,1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]])
        genmat2 = _MS(n)([[1,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [0,1,0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1],
                         [0,0,1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0,0,0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0,0,0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0,0,0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0,0,0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0,0,0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0,0,0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0,0,0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1]])
        #  G = PermutationGroup( [ "(9,10)(19,20)", "(9,19)(10,20)", "(8,9)(18,19)", "(7,8)(17,18)",\
        #        "(6,7)(16,17)", "(5,6)(15,16)", "(4,5)(14,15)", "(3,4)(13,14)",\
        #        "(2,3)(12,13)", "(1,2)(11,12)"] )
        spectrum = [1, 0, 0, 0, 45, 0, 0, 0, 210, 0, 512, 0, 210, 0, 0, 0, 45, 0, 0, 0, 1]
        self_dual_codes_20_8 = {"order autgp":1857945600,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Huge aut gp. Min dist 4."}
        # [20,9]: (genmat, K20 are equiv)
        genmat = _I2(n).augment(A10)
        K20 = _MS(n)([[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
                      [1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0]])
        #genmat = K20 # not in standard form
        #  G = PermutationGroup( [ "(4,13)(5,15)", "(4,15)(5,13)", "(3,4,13)(5,11,15)",
        #   "(3,4,6,11,15,17)(5,13)", "(3,5,17,4,12)(6,15,7,11,13)",
        #   "(1,2)(3,5,17,4,7,11,13,6,15,12)", "(1,3,5,17,4,12)(2,11,13,6,15,7)",
        #   "(3,5,17,4,12)(6,15,7,11,13)(10,18)(19,20)", "(3,5,17,4,12)(6,15,7,11,13)(10,19)(18,20)",
        #   "(3,5,17,4,12)(6,15,7,11,13)(9,10)(16,18)",
        #   "(3,5,17,4,12)(6,15,7,11,13)(8,9)(14,16)" ] )
        spectrum = [1, 0, 0, 0, 21, 0, 48, 0, 234, 0, 416, 0, 234, 0, 48, 0, 21, 0, 0, 0, 1]
        self_dual_codes_20_9 = {"order autgp":4423680,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Min dist 4."}
        # [20,10]
        L20 = _MS(n)([[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
                      [0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,1,0,0,0,0],
                      [0,1,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0]])
        genmat = L20 # not in standard form
        # G = PermutationGroup( [ "(17,18)(19,20)", "(17,19)(18,20)", "(15,16)(19,20)",
        #        "(15,17)(16,18)", "(10,11)(12,13)", "(10,12)(11,13)", "(9,10)(13,14)",
        #        "(8,9)(12,13)", "(3,4)(5,6)", "(3,5)(4,6)", "(2,3)(6,7)", "(1,2)(5,6)",
        #        "(1,8)(2,9)(3,10)(4,11)(5,12)(6,13)(7,14)(19,20)" ] ) # order 1354752
        spectrum = [1, 0, 0, 0, 17, 0, 56, 0, 238, 0, 400, 0, 238, 0, 56, 0, 17, 0, 0, 0, 1]
        self_dual_codes_20_10 = {"order autgp":1354752,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Min dist 4."}
        # [20,11]
        S20 = _MS(n)([[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
                      [1,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0],
                      [1,1,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,0,0],
                      [1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0]] )
        genmat = S20 # not in standard form
        # G = PermutationGroup( [ "(17,18)(19,20)", "(17,19)(18,20)", "(13,14)(15,16)",
        #    "(13,15)(14,16)", "(11,12)(15,16)", "(11,13)(12,14)", "(9,10)(15,16)",
        #    "(9,11)(10,12)", "(5,6)(7,8)", "(5,7)(6,8)", "(3,4)(7,8)", "(3,5)(4,6)",
        #    "(1,2)(7,8)", "(1,3)(2,4)", "(1,9)(2,10)(3,11)(4,12)(5,13)(6,14)(7,15)(8,16)" ] )
        # G.order() = 294912
        spectrum = [1, 0, 0, 0, 13, 0, 64, 0, 242, 0, 384, 0, 242, 0, 64, 0, 13, 0, 0, 0, 1]
        self_dual_codes_20_11 = {"order autgp":294912,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Min dist 4."}
        # [20,12]
        R20 = _MS(n)([[0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
                      [0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,1,0],
                      [1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0],
                      [1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1],
                      [1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1]])
        genmat = R20 # not in standard form
        #  G = PermutationGroup( [ "(17,18)(19,20)", "(17,19)(18,20)", "(15,16)(19,20)",
        #    "(15,17)(16,18)", "(11,12)(13,14)", "(11,13)(12,14)", "(9,10)(13,14)",
        #    "(9,11)(10,12)", "(5,6)(7,8)", "(5,7)(6,8)", "(3,4)(7,8)", "(3,5)(4,6)",
        #    "(3,9,15)(4,10,16)(5,11,17)(6,12,18)(7,14,19)(8,13,20)",
        #    "(1,2)(7,8)(9,15)(10,16)(11,17)(12,18)(13,19)(14,20)" ] ) # order 82944
        spectrum = [1, 0, 0, 0, 9, 0, 72, 0, 246, 0, 368, 0, 246, 0, 72, 0, 9, 0, 0, 0, 1]
        self_dual_codes_20_12 = {"order autgp":82944,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Min dist 4."}
        # [20,13]
        M20 = _MS(n)([[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
                      [0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
                      [1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0],
                      [0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
                      [0,0,1,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,1,1,0,1,1,0,1,0,0,1,0,0,0,0]])
        genmat = M20 # not in standard form
        #  G = PermutationGroup( [ "(17,18)(19,20)", "(17,19)(18,20)", "(13,14)(15,16)",
        #            "(13,15)(14,16)", "(9,10)(11,12)", "(9,11)(10,12)", "(5,6)(7,8)",
        #            "(5,7)(6,8)", "(5,9)(6,11)(7,12)(8,10)(13,17)(14,19)(15,18)(16,20)",
        #            "(5,13)(6,15)(7,14)(8,16)(9,17)(10,20)(11,18)(12,19)",
        #            "(3,4)(6,7)(11,12)(13,17)(14,18)(15,19)(16,20)",
        #            "(2,3)(7,8)(9,13)(10,14)(11,15)(12,16)(19,20)",
        #            "(1,2)(6,7)(11,12)(13,17)(14,18)(15,19)(16,20)",
        #            "(1,5)(2,6)(3,7)(4,8)(9,17)(10,18)(11,19)(12,20)" ] )
        spectrum = [1, 0, 0, 0, 5, 0, 80, 0, 250, 0, 352, 0, 250, 0, 80, 0, 5, 0, 0, 0, 1]
        self_dual_codes_20_13 = {"order autgp":122880,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "Min dist 4."}
        # [20,14]:  # aut gp of this computed using a program by Robert Miller
        A0 = self_dual_binary_codes(n-2)["18"]["8"]["code"].redundancy_matrix()
        genmat = _I2(n).augment(block_diagonal_matrix([A0,_matId(n)[9]]))
        # [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        #  [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
        #  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0],
        #  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0],
        #  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0],
        #  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0],
        #  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
        #  G = PermutationGroup( [ "(8,19)(16,17)", "(8,16)(17,19)", "(9,18)(16,17)", "(8,9)(18,19)",
        #                 "(7,8)(17,18)", "(4,15)(5,14)", "(4,5)(14,15)", "(4,15)(6,11)", "(5,6)(11,14)",
        #                 "(3,13)(4,15)", "(3,15)(4,13)", "(1,2)(4,15)", "(1,4)(2,15)(3,5)(13,14)", "(10,20)" ] )
        spectrum = [1, 0, 1, 0, 17, 0, 68, 0, 238, 0, 374, 0, 238, 0, 68, 0, 17, 0, 1, 0, 1]
        self_dual_codes_20_14 = {"order autgp":645120,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment": "'Exceptional' construction."}
        # [20,15]:
        A0 = self_dual_binary_codes(n-2)["18"]["7"]["code"].redundancy_matrix()
        genmat = _I2(n).augment(block_diagonal_matrix([A0,_matId(n)[9]]))
        # [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
        #  [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0],
        #  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0],
        #  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0],
        #  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0],
        #  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0],
        #  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        #  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
        #  G = PermutationGroup( [ "(10,20)", "(9,11)(17,19)", "(9,17)(11,19)", "(8,9)(15,17)",
        #     "(7,12)(13,18)", "(7,13)(12,18)", "(5,6)(12,13)", "(5,7)(6,18)",
        #     "(4,14)(5,8)(6,15)(7,9)(11,13)(12,19)(17,18)", "(3,4)(14,16)",
        #     "(1,2)(5,8)(6,15)(7,9)(11,13)(12,19)(17,18)", "(1,3)(2,16)",
        #     "(1,5)(2,6)(3,7)(4,12)(11,19)(13,14)(16,18)" ] ) # order 165888
        spectrum = [1, 0, 1, 0, 9, 0, 84, 0, 246, 0, 342, 0, 246, 0, 84, 0, 9, 0, 1, 0, 1]
        self_dual_codes_20_15 = {"order autgp":165888,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"'Exceptional' construction. Unique lowest wt codeword."}
        self_dual_codes["20"] = {"0":self_dual_codes_20_0,"1":self_dual_codes_20_1,"2":self_dual_codes_20_2,
                  "3":self_dual_codes_20_3,"4":self_dual_codes_20_4,"5":self_dual_codes_20_5,
                  "6":self_dual_codes_20_6,"7":self_dual_codes_20_7,"8":self_dual_codes_20_8,
                  "9":self_dual_codes_20_9,"10":self_dual_codes_20_10,"11":self_dual_codes_20_11,
                  "12":self_dual_codes_20_12,"13":self_dual_codes_20_13,"14":self_dual_codes_20_14,
                  "15":self_dual_codes_20_15}
        return self_dual_codes

    if n == 22:
        # all of these of these are Type I; 2 of these codes
        # are formally equivalent but with different automorphism groups
        #    *** Incomplete ***   (7 out of 25)
        # [22,0]:
        genmat = _I2(n).augment(_I2(n))
        #    G = PermutationGroup( [ "(11,22)", "(10,11)(21,22)", "(9,10)(20,21)",\
        #        "(8,9)(19,20)", "(7,8)(18,19)", "(6,7)(17,18)", "(5,6)(16,17)",\
        #        "(4,5)(15,16)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] ) # S_11x(ZZ/2ZZ)^11??
        spectrum = [1, 0, 11, 0, 55, 0, 165, 0, 330, 0, 462, 0, 462, 0, 330, 0, 165, 0, 55, 0, 11, 0, 1]
        self_dual_codes_22_0 = {"order autgp":81749606400,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Huge aut gp."}
        # [22,1]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matId(n)[4]]))
        #   G = PermutationGroup( [ "(11,22)", "(10,11)(21,22)", "(9,10)(20,21)",\
        #         "(8,9)(19,20)", "(7,8)(18,19)", "(6,7)(17,18)", "(5,6)(16,17)",\
        #         "(4,12)(13,14)", "(4,13)(12,14)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] )
        spectrum = [1, 0, 7, 0, 35, 0, 133, 0, 330, 0, 518, 0, 518, 0, 330, 0, 133, 0, 35, 0, 7, 0, 1]
        self_dual_codes_22_1 = {"order autgp":867041280,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [22,2]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matId(n)[6]]))
        #   G = PermutationGroup( [ "(11,22)", "(10,11)(21,22)", "(9,10)(20,21)",\
        #         "(8,9)(19,20)", "(7,8)(18,19)", "(5,6)(16,17)", "(5,16)(6,17)",\
        #         "(4,5)(15,16)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] )
        spectrum = [1, 0, 5, 0, 25, 0, 117, 0, 330, 0, 546, 0, 546, 0, 330, 0, 117, 0, 25, 0, 5, 0, 1]
        self_dual_codes_22_2 = {"order autgp":88473600,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":""}
        # [22,3]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[8],_matId(n)[8]]))
        #   G = PermutationGroup( [ "(11,22)", "(10,11)(21,22)", "(9,10)(20,21)",\
        #          "(7,8)(18,19)", "(7,18)(8,19)", "(6,7)(17,18)", "(5,6)(16,17)",\
        #          "(4,5)(15,16)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] )
        spectrum = [1, 0, 3, 0, 31, 0, 85, 0, 282, 0, 622, 0, 622, 0, 282, 0, 85, 0, 31, 0, 3, 0, 1]
        self_dual_codes_22_3 = {"order autgp":247726080,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Same spectrum as the '[20,5]' code."}
        # [22,4]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[10],_matId(n)[10]]))
        #   G = PermutationGroup( [ "(11,22)", "(9,10)(20,21)", "(9,20)(10,21)",\
        #        "(8,9)(19,20)", "(7,8)(18,19)", "(6,7)(17,18)", "(5,6)(16,17)",\
        #        "(4,5)(15,16)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] )
        spectrum = [1, 0, 1, 0, 45, 0, 45, 0, 210, 0, 722, 0, 722, 0, 210, 0, 45, 0, 45, 0, 1, 0, 1]
        self_dual_codes_22_4 = {"order autgp":3715891200,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Unique lowest weight codeword."}
        # [22,5]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[4],_matA(n)[4],_matId(n)[8]]))
        #   G = PermutationGroup( [ "(11,22)", "(10,11)(21,22)", "(9,10)(20,21)",\
        #         "(8,16)(17,18)", "(8,17)(16,18)", "(7,8)(18,19)", "(6,7)(17,18)",\
        #         "(5,6)(16,17)", "(4,12)(13,14)", "(4,13)(12,14)", "(3,4)(14,15)",\
        #         "(2,3)(13,14)", "(1,2)(12,13)", "(1,5)(2,6)(3,7)(4,8)(12,16)(13,17)(14,18)(15,19)" ] )
        spectrum = [1, 0, 3, 0, 31, 0, 85, 0, 282, 0, 622, 0, 622, 0, 282, 0, 85, 0, 31, 0, 3, 0, 1]
        self_dual_codes_22_5 = {"order autgp":173408256,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Same spectrum as the '[20,3]' code."}
        # [22,6]:
        genmat = _I2(n).augment(block_diagonal_matrix([_matA(n)[6],_matA(n)[4],_matId(n)[10]]))
        #   G = PermutationGroup( [ "(11,22)", "(10,18)(19,20)", "(10,19)(18,20)",\
        #         "(9,10)(20,21)", "(8,9)(19,20)", "(7,8)(18,19)", "(5,6)(16,17)",\
        #         "(5,16)(6,17)", "(4,5)(15,16)", "(3,4)(14,15)", "(2,3)(13,14)", "(1,2)(12,13)" ] )
        spectrum = [1, 0, 1, 0, 29, 0, 61, 0, 258, 0, 674, 0, 674, 0, 258, 0, 61, 0, 29, 0, 1, 0, 1]
        self_dual_codes_22_6 = {"order autgp":61931520,"code":LinearCode(genmat),"spectrum":spectrum,
                 "Type":"I","Comment":"Unique lowest weight codeword."}
        self_dual_codes["22"] = {"0":self_dual_codes_22_0,"1":self_dual_codes_22_1,"2":self_dual_codes_22_2,
                          "3":self_dual_codes_22_3,"4":self_dual_codes_22_4,"5":self_dual_codes_22_5,
                          "6":self_dual_codes_22_6}
        return self_dual_codes
