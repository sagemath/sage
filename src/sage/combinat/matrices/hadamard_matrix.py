r"""
Hadamard matrices

A Hadamard matrix is an `n\times n` matrix `H` whose entries are either `+1` or `-1`
and whose rows are mutually orthogonal. For example, the matrix `H_2`
defined by

.. MATH::

    \left(\begin{array}{rr}
    1 & 1 \\
    1 & -1
    \end{array}\right)

is a Hadamard matrix. An `n\times n` matrix `H` whose entries are either `+1` or
`-1` is a Hadamard matrix if and only if:

(a) `|det(H)|=n^{n/2}` or

(b)  `H*H^t = n\cdot I_n`, where `I_n` is the identity matrix.

In general, the tensor product of an `m\times m` Hadamard matrix and an
`n\times n` Hadamard matrix is an `(mn)\times (mn)` matrix. In
particular, if there is an `n\times n` Hadamard matrix then there is a
`(2n)\times (2n)` Hadamard matrix (since one may tensor with `H_2`).
This particular case is sometimes called the Sylvester construction.

The Hadamard conjecture (possibly due to Paley) states that a Hadamard
matrix of order `n` exists if and only if `n= 1, 2` or `n` is a multiple
of `4`.

The module below implements the Paley constructions (see for example
[Hora]_) and the Sylvester construction. It also allows you to pull a
Hadamard matrix from the database at [HadaSloa]_.

AUTHORS:

- David Joyner (2009-05-17): initial version

REFERENCES:

.. [HadaSloa] N.J.A. Sloane's Library of Hadamard Matrices, at
   http://neilsloane.com/hadamard/
.. [HadaWiki] Hadamard matrices on Wikipedia, :wikipedia:`Hadamard_matrix`
.. [Hora] K. J. Horadam, Hadamard Matrices and Their Applications,
   Princeton University Press, 2006.
"""
from sage.rings.arith import kronecker_symbol
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix, block_matrix, block_diagonal_matrix, diagonal_matrix
from urllib import urlopen
from sage.misc.functional import is_even
from sage.rings.arith import is_prime, is_square, is_prime_power, divisors
from math import sqrt
from sage.matrix.constructor import identity_matrix as I
from sage.matrix.constructor import ones_matrix     as J

def H1(i, j, p):
    """
    Returns the i,j-th entry of the Paley matrix, type I case.
    The Paley type I case corresponds to the case `p \cong 3 \mod{4}`
    for a prime `p`.

    .. TODO::

        This construction holds more generally for prime powers `q`
        congruent to `3 \mod{4}`. We should implement these but we
        first need to implement Quadratic character for `GF(q)`.

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.H1(1,2,3)
        1
    """
    if i == 0 or j == 0:
        return 1
    # what follows will not be executed for (i, j) = (0, 0).
    if i == j:
        return -1
    return -kronecker_symbol(i - j, p)

def H2(i, j, p):
    """
    Returns the i,j-th entry of the Paley matrix, type II case.

    The Paley type II case corresponds to the case `p \cong 1 \mod{4}`
    for a prime `p` (see [Hora]_).

    .. TODO::

        This construction holds more generally for prime powers `q`
        congruent to `1 \mod{4}`. We should implement these but we
        first need to implement Quadratic character for `GF(q)`.

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.H2(1,2,5)
        1
    """
    if i == 0 and j == 0:
        return 0
    if i == 0 or j == 0:
        return 1
    if i == j:
        return 0
    return kronecker_symbol(i - j, p)

def normalise_hadamard(H):
    """
    Return the normalised Hadamard matrix corresponding to ``H``.

    The normalised Hadamard matrix corresponding to a Hadamard matrix `H` is a
    matrix whose every entry in the first row and column is +1.

    EXAMPLES::

        sage: H = sage.combinat.matrices.hadamard_matrix.normalise_hadamard(hadamard_matrix(4))
        sage: H == hadamard_matrix(4)
        True
    """
    Hc1 = H.column(0)
    Hr1 = H.row(0)
    for i in range(H.ncols()):
        if Hc1[i] < 0:
            H.rescale_row(i, -1)
    for i in range(H.nrows()):
        if Hr1[i] < 0:
            H.rescale_col(i, -1)
    return H

def hadamard_matrix_paleyI(n):
    """
    Implements the Paley type I construction.

    The Paley type I case corresponds to the case `p \cong 3 \mod{4}` for a
    prime `p` (see [Hora]_).

    EXAMPLES:

    We note that this method returns a normalised Hadamard matrix ::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyI(4)
        [ 1  1  1  1]
        [ 1 -1 -1  1]
        [ 1  1 -1 -1]
        [ 1 -1  1 -1]
    """
    p = n - 1
    if not(is_prime(p) and (p % 4 == 3)):
        raise ValueError("The order %s is not covered by the Paley type I construction." % n)
    H = matrix(ZZ, [[H1(i, j, p) for i in range(n)] for j in range(n)])
    # normalising H so that first row and column have only +1 entries.
    return normalise_hadamard(H)

def hadamard_matrix_paleyII(n):
    """
    Implements the Paley type II construction.

    The Paley type II case corresponds to the case `p \cong 1 \mod{4}` for a
    prime `p` (see [Hora]_).

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12).det()
        2985984
        sage: 12^6
        2985984

    We note that the method returns a normalised Hadamard matrix ::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12)
        [ 1  1  1  1  1  1| 1  1  1  1  1  1]
        [ 1  1  1 -1 -1  1|-1 -1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1|-1  1 -1  1 -1 -1]
        [ 1 -1  1  1  1 -1|-1 -1  1 -1  1 -1]
        [ 1 -1 -1  1  1  1|-1 -1 -1  1 -1  1]
        [ 1  1 -1 -1  1  1|-1  1 -1 -1  1 -1]
        [-----------------+-----------------]
        [ 1 -1 -1 -1 -1 -1|-1  1  1  1  1  1]
        [ 1 -1  1 -1 -1  1| 1 -1 -1  1  1 -1]
        [ 1  1 -1  1 -1 -1| 1 -1 -1 -1  1  1]
        [ 1 -1  1 -1  1 -1| 1  1 -1 -1 -1  1]
        [ 1 -1 -1  1 -1  1| 1  1  1 -1 -1 -1]
        [ 1  1 -1 -1  1 -1| 1 -1  1  1 -1 -1]
    """
    N = Integer(n/2)
    p = N - 1
    if not(is_prime(p) and (p % 4 == 1)):
        raise ValueError("The order %s is not covered by the Paley type II construction." % n)
    S = matrix(ZZ, [[H2(i, j, p) for i in range(N)] for j in range(N)])
    H = block_matrix([[S + 1, S - 1], [1 - S, S + 1]])
    # normalising H so that first row and column have only +1 entries.
    return normalise_hadamard(H)

def hadamard_matrix(n):
    """
    Tries to construct a Hadamard matrix using a combination of Paley
    and Sylvester constructions.

    EXAMPLES::

        sage: hadamard_matrix(12).det()
        2985984
        sage: 12^6
        2985984
        sage: hadamard_matrix(1)
        [1]
        sage: hadamard_matrix(2)
        [ 1  1]
        [ 1 -1]
        sage: hadamard_matrix(8)
        [ 1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1]
        sage: hadamard_matrix(8).det() == 8^4
        True

    We note that the method `hadamard_matrix()` returns a normalised Hadamard matrix
    (the entries in the first row and column are all +1) ::

        sage: hadamard_matrix(12)
        [ 1  1  1  1  1  1| 1  1  1  1  1  1]
        [ 1  1  1 -1 -1  1|-1 -1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1|-1  1 -1  1 -1 -1]
        [ 1 -1  1  1  1 -1|-1 -1  1 -1  1 -1]
        [ 1 -1 -1  1  1  1|-1 -1 -1  1 -1  1]
        [ 1  1 -1 -1  1  1|-1  1 -1 -1  1 -1]
        [-----------------+-----------------]
        [ 1 -1 -1 -1 -1 -1|-1  1  1  1  1  1]
        [ 1 -1  1 -1 -1  1| 1 -1 -1  1  1 -1]
        [ 1  1 -1  1 -1 -1| 1 -1 -1 -1  1  1]
        [ 1 -1  1 -1  1 -1| 1  1 -1 -1 -1  1]
        [ 1 -1 -1  1 -1  1| 1  1  1 -1 -1 -1]
        [ 1  1 -1 -1  1 -1| 1 -1  1  1 -1 -1]
    """
    if not(n % 4 == 0) and (n > 2):
        raise ValueError("The Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        return matrix([[1, 1], [1, -1]])
    if is_even(n):
        N = Integer(n / 2)
    elif n == 1:
        return matrix([1])
    if is_prime(N - 1) and (N - 1) % 4 == 1:
        return hadamard_matrix_paleyII(n)
    elif n == 4 or n % 8 == 0:
        had = hadamard_matrix(Integer(n / 2))
        chad1 = matrix([list(r) + list(r) for r in had.rows()])
        mhad = (-1) * had
        R = len(had.rows())
        chad2 = matrix([list(had.rows()[i]) + list(mhad.rows()[i])
                       for i in range(R)])
        return chad1.stack(chad2)
    elif is_prime(N - 1) and (N - 1) % 4 == 3:
        return hadamard_matrix_paleyI(n)
    else:
        raise ValueError("The Hadamard matrix of order %s is not yet implemented." % n)

def hadamard_matrix_www(url_file, comments=False):
    """
    Pulls file from Sloane's database and returns the corresponding Hadamard
    matrix as a Sage matrix.

    You must input a filename of the form "had.n.xxx.txt" as described
    on the webpage http://neilsloane.com/hadamard/, where
    "xxx" could be empty or a number of some characters.

    If comments=True then the "Automorphism..." line of the had.n.xxx.txt
    file is printed if it exists. Otherwise nothing is done.

    EXAMPLES::

        sage: hadamard_matrix_www("had.4.txt")      # optional - internet
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1  1 -1 -1]
        [ 1 -1 -1  1]
        sage: hadamard_matrix_www("had.16.2.txt",comments=True)   # optional - internet
        Automorphism group has order = 49152 = 2^14 * 3
        [ 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1]
        [ 1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1]
        [ 1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1]
        [ 1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1]
        [ 1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1 -1]
        [ 1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1]
        [ 1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1]
        [ 1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1]
    """
    n = eval(url_file.split(".")[1])
    rws = []
    url = "http://neilsloane.com/hadamard/" + url_file
    f = urlopen(url)
    s = f.readlines()
    for i in range(n):
        r = []
        for j in range(n):
            if s[i][j] == "+":
                r.append(1)
            else:
                r.append(-1)
        rws.append(r)
    f.close()
    if comments:
        lastline = s[-1]
        if lastline[0] == "A":
            print lastline
    return matrix(rws)

_rshcd_cache = {}

def regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e,existence=False):
    r"""
    Return a Regular Symmetric Hadamard Matrix with Constant Diagonal.

    A Hadamard matrix is said to be *regular* if its rows all sum to the same
    value.

    When `\epsilon\in\{-1,+1\}`, we say that `M` is a `(n,\epsilon)-RSHCD` if
    `M` is a regular symmetric Hadamard matrix with constant diagonal
    `\delta\in\{-1,+1\}` and row values all equal to `\delta \epsilon
    \sqrt(n)`. For more information, see [HX10]_ or 10.5.1 in
    [BH12]_.

    INPUT:

    - ``n`` (integer) -- side of the matrix

    - ``e`` -- one of `-1` or `+1`, equal to the value of `\epsilon`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,1)
        [ 1  1  1 -1]
        [ 1  1 -1  1]
        [ 1 -1  1  1]
        [-1  1  1  1]
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,-1)
        [ 1 -1 -1 -1]
        [-1  1 -1 -1]
        [-1 -1  1 -1]
        [-1 -1 -1  1]

    Other hardcoded values::

        sage: for n,e in [(36,1),(36,-1),(100,1),(100,-1),(196, 1)]:
        ....:     print regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e)
        36 x 36 dense matrix over Integer Ring
        36 x 36 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        196 x 196 dense matrix over Integer Ring

    From two close prime powers::

        sage: print regular_symmetric_hadamard_matrix_with_constant_diagonal(64,-1)
        64 x 64 dense matrix over Integer Ring

    Recursive construction::

        sage: print regular_symmetric_hadamard_matrix_with_constant_diagonal(144,-1)
        144 x 144 dense matrix over Integer Ring

    REFERENCE:

    .. [BH12] A. Brouwer and W. Haemers,
      Spectra of graphs,
      Springer, 2012.

    .. [HX10] W. Haemers and Q. Xiang,
      Strongly regular graphs with parameters `(4m^4,2m^4+m^2,m^4+m^2,m^4+m^2)` exist for all `m>1`,
      European Journal of Combinatorics,
      Volume 31, Issue 6, August 2010, Pages 1553-1559,
      http://dx.doi.org/10.1016/j.ejc.2009.07.009.
    """
    if existence and (n,e) in _rshcd_cache:
        return _rshcd_cache[n,e]

    from sage.graphs.strongly_regular_db import strongly_regular_graph

    def true():
        _rshcd_cache[n,e] = True
        return True

    M = None
    if abs(e) != 1:
        raise ValueError
    if n<0:
        if existence:
            return False
        raise ValueError
    elif n == 4:
        if existence:
            return true()
        if e == 1:
            M = J(4)-2*matrix(4,[[int(i+j == 3) for i in range(4)] for j in range(4)])
        else:
            M = -J(4)+2*I(4)
    elif n ==  36:
        if existence:
            return true()
        if e == 1:
            M = strongly_regular_graph(36, 15, 6, 6).adjacency_matrix()
            M = J(36) - 2*M
        else:
            M = strongly_regular_graph(36,14,4,6).adjacency_matrix()
            M =  -J(36) + 2*M + 2*I(36)
    elif n == 100:
        if existence:
            return true()
        if e == -1:
            M = strongly_regular_graph(100,44,18,20).adjacency_matrix()
            M = 2*M - J(100) + 2*I(100)
        else:
            M = strongly_regular_graph(100,45,20,20).adjacency_matrix()
            M = J(100) - 2*M
    elif n == 196 and e == 1:
        if existence:
            return true()
        M = strongly_regular_graph(196,91,42,42).adjacency_matrix()
        M = J(196) - 2*M
    elif (  e  == 1                 and
          n%16 == 0                 and
          is_square(n)              and
          is_prime_power(sqrt(n)-1) and
          is_prime_power(sqrt(n)+1)):
        if existence:
            return true()
        M = -rshcd_from_close_prime_powers(int(sqrt(n)))

    # Recursive construction: the kronecker product of two RSHCD is a RSHCD
    else:
        from itertools import product
        for n1,e1 in product(divisors(n)[1:-1],[-1,1]):
            e2 = e1*e
            n2 = n//n1
            if (regular_symmetric_hadamard_matrix_with_constant_diagonal(n1,e1,existence=True) and
                regular_symmetric_hadamard_matrix_with_constant_diagonal(n2,e2,existence=True)):
                if existence:
                    return true()
                M1 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n1,e1)
                M2 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n2,e2)
                M  = M1.tensor_product(M2)
                break

    if M is None:
        from sage.misc.unknown import Unknown
        _rshcd_cache[n,e] = Unknown
        if existence:
            return Unknown
        raise ValueError("I do not know how to build a {}-RSHCD".format((n,e)))

    assert M*M.transpose() == n*I(n)
    assert set(map(sum,M)) == {e*sqrt(n)}

    return M

def _helper_payley_matrix(n):
    r"""
    Return the marix constructed in Lemma 1.19 page 291 of [SWW72]_.

    This function return a `n^2` matrix `M` whose rows/columns are indexed by
    the element of a finite field on `n` elements `x_1,...,x_n`. The value
    `M_{i,j}` is equal to `\chi(x_i-x_j)`. Note that `n` must be an odd prime power.

    The elements `x_1,...,x_n` are ordered in such a way that the matrix is
    symmetric with respect to its second diagonal. The matrix is symmetric if
    n==4k+1, and skew-symmetric if n=4k-1.

    INPUT:

    - ``n`` -- a prime power

    .. SEEALSO::

        :func:`rshcd_from_close_prime_powers`

    EXAMPLE::

        sage: from sage.combinat.matrices.hadamard_matrix import _helper_payley_matrix
        sage: _helper_payley_matrix(5)
        [ 0  1 -1 -1  1]
        [ 1  0  1 -1 -1]
        [-1  1  0  1 -1]
        [-1 -1  1  0  1]
        [ 1 -1 -1  1  0]
    """
    from sage.rings.finite_rings.constructor import FiniteField as GF
    K = GF(n,conway=True,prefix='x')

    # Order the elements of K in K_list
    # so that K_list[i] = -K_list[n-i-1]
    K_pairs = set(frozenset([x,-x]) for x in K)
    K_pairs.discard(frozenset([0]))
    K_list = [None]*n
    for i,(x,y) in enumerate(K_pairs):
        K_list[i]   = x
        K_list[-i-1] = y
    K_list[n//2] = K(0)

    M = matrix(n,[[2*((x-y).is_square())-1
                   for x in K_list]
                  for y in K_list])
    M = M-I(n)
    assert (M*J(n)).is_zero()
    assert (M*M.transpose()) == n*I(n)-J(n)
    return M

def rshcd_from_close_prime_powers(n):
    r"""
    Return a `(n^2,1)`-RSHCD when `n-1` and `n+1` are odd prime powers and `n=0\pmod{4}`.

    The construction implemented here appears in Theorem 4.3 from [GS14]_.

    Note that the authors of [SWW72]_ claim in Corollary 5.12 (page 342) to have
    proved the same result without the `n=0\pmod{4}` restriction with a *very*
    similar construction. So far, however, I (Nathann Cohen) have not been able
    to make it work.

    INPUT:

    - ``n`` -- an integer congruent to `0\pmod{4}`

    .. SEEALSO::

        :func:`regular_symmetric_hadamard_matrix_with_constant_diagonal`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import rshcd_from_close_prime_powers
        sage: rshcd_from_close_prime_powers(4)
        [-1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1]
        [-1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1]
        [ 1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1]
        [-1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1]
        [ 1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1]
        [-1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1]
        [-1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1]
        [ 1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1]
        [-1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1  1]
        [ 1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1]
        [-1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1]
        [-1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1]
        [-1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1]
        [-1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1]

    REFERENCE:

    .. [GS14] J.M. Goethals, and J. J. Seidel,
      Strongly regular graphs derived from combinatorial designs,
      Canadian Journal of Mathematics 22(1970), 597-614,
      http://dx.doi.org/10.4153/CJM-1970-067-9

    .. [SWW72] A Street, W. Wallis, J. Wallis,
      Combinatorics: Room squares, sum-free sets, Hadamard matrices.
      Lecture notes in Mathematics 292 (1972).
    """
    if n%4:
        raise ValueError("n(={}) must be congruent to 0 mod 4")

    a,b = sorted([n-1,n+1],key=lambda x:-x%4)
    Sa  = _helper_payley_matrix(a)
    Sb  = _helper_payley_matrix(b)
    U   = matrix(a,[[int(i+j == a-1) for i in range(a)] for j in range(a)])

    K = (U*Sa).tensor_product(Sb) + U.tensor_product(J(b)-I(b)) - J(a).tensor_product(I(b))

    F = lambda x:diagonal_matrix([-(-1)**i for i in range(x)])
    G = block_diagonal_matrix([J(1),I(a).tensor_product(F(b))])
    e = matrix(a*b,[1]*(a*b))
    H = block_matrix(2,[-J(1),e.transpose(),e,K])

    HH = G*H*G
    assert len(set(map(sum,HH))) == 1
    assert HH**2 == n**2*I(n**2)
    return HH
