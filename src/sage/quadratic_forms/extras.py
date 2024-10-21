# sage_setup: distribution = sagemath-modules
"""
Extra functions for quadratic forms
"""

from sage.matrix.constructor import matrix
from sage.structure.element import Matrix
from sage.arith.misc import legendre_symbol
from sage.rings.integer_ring import ZZ


def is_triangular_number(n, return_value=False):
    """
    Return whether ``n`` is a triangular number.

    A *triangular number* is a number of the form `k(k+1)/2` for some
    nonnegative integer `n`. See :wikipedia:`Triangular_number`. The sequence
    of triangular number is references as A000217 in the Online encyclopedia of
    integer sequences (OEIS).

    If you want to get the value of `k` for which `n=k(k+1)/2` set the
    argument ``return_value`` to ``True`` (see the examples below).

    INPUT:

    - ``n`` -- integer

    - ``return_value`` -- boolean (default: ``False``); if set to
      ``True`` the function returns a pair made of a boolean and the value `v`
      such that `v(v+1)/2 = n`

    EXAMPLES::

        sage: is_triangular_number(3)
        True
        sage: is_triangular_number(3, return_value=True)
        (True, 2)

        sage: is_triangular_number(2)
        False
        sage: is_triangular_number(2, return_value=True)
        (False, None)

        sage: is_triangular_number(25*(25+1)/2)
        True

        sage: is_triangular_number(10^6 * (10^6 +1)/2, return_value=True)
        (True, 1000000)

    TESTS::

        sage: F1 = [n for n in range(1,100*(100+1)/2)
        ....:       if is_triangular_number(n)]
        sage: F2 = [n*(n+1)/2 for n in range(1,100)]
        sage: F1 == F2
        True

        sage: for n in range(1000):
        ....:     res,v = is_triangular_number(n,return_value=True)
        ....:     assert res == is_triangular_number(n)
        ....:     if res: assert v*(v+1)/2 == n
    """
    n = ZZ(n)

    if return_value:
        if n < 0:
            return (False, None)
        if n == 0:
            return (True, ZZ.zero())
        s, r = (8 * n + 1).sqrtrem()
        if r:
            return (False, None)
        return (True, (s - 1) // 2)

    return (8 * n + 1).is_square()


def extend_to_primitive(A_input):
    r"""
    Given a matrix (resp. list of vectors), extend it to a square
    matrix (resp. list of vectors), such that its determinant is the
    gcd of its minors (i.e. extend the basis of a lattice to a
    "maximal" one in `\ZZ^n`).

    Author(s): Gonzalo Tornaria and Jonathan Hanke.

    INPUT:

    - ``A_input`` -- a matrix or a list of length n vectors (in the same space)

    OUTPUT: a square matrix or a list of n vectors (resp.)

    EXAMPLES::

        sage: A = Matrix(ZZ, 3, 2, range(6))
        sage: extend_to_primitive(A)
        [ 0  1 -1]
        [ 2  3  0]
        [ 4  5  0]

        sage: extend_to_primitive([vector([1,2,3])])
        [(1, 2, 3), (0, 1, 1), (-1, 0, 0)]
    """
    # Deal with a list of vectors
    if not isinstance(A_input, Matrix):
        A = matrix(A_input)      # Make a matrix A with the given rows.
        vec_output_flag = True
    else:
        A = A_input
        vec_output_flag = False

    # Arrange for A  to have more columns than rows.
    if A.is_square():
        return A
    if A.nrows() > A.ncols():
        return extend_to_primitive(A.transpose()).transpose()

    # Setup
    k = A.nrows()
    n = A.ncols()
    R = A.base_ring()

    # Smith normal form transformation, assuming more columns than rows
    V = A.smith_form()[2]

    # Extend the matrix in new coordinates, then switch back.
    B = A * V
    B_new = matrix(R, n - k, n)
    for i in range(n - k):
        B_new[i, n - i - 1] = 1
    C = B.stack(B_new)
    D = C * V**(-1)

    # Normalize for a positive determinant
    if D.det() < 0:
        D.rescale_row(n - 1, -1)

    # Return the current information
    if vec_output_flag:
        return D.rows()
    else:
        return D


def least_quadratic_nonresidue(p):
    r"""
    Return the smallest positive integer quadratic non-residue in `\ZZ/p\ZZ` for primes `p>2`.

    EXAMPLES::

        sage: least_quadratic_nonresidue(5)
        2
        sage: [least_quadratic_nonresidue(p) for p in prime_range(3, 100)]              # needs sage.libs.pari
        [2, 2, 3, 2, 2, 3, 2, 5, 2, 3, 2, 3, 2, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5]

    TESTS:

    This raises an error if input is a positive composite integer. ::

        sage: least_quadratic_nonresidue(20)
        Traceback (most recent call last):
        ...
        ValueError: p must be a prime number > 2

    This raises an error if input is 2. This is because every integer is a
    quadratic residue modulo 2. ::

        sage: least_quadratic_nonresidue(2)
        Traceback (most recent call last):
        ...
        ValueError: there are no quadratic non-residues in Z/2Z
    """
    p1 = abs(p)

    # Deal with the prime p = 2 and |p| <= 1.
    if p1 == 2:
        raise ValueError("there are no quadratic non-residues in Z/2Z")
    if p1 < 2:
        raise ValueError("p must be a prime number > 2")

    # Find the smallest non-residue mod p
    # For 7/8 of primes the answer is 2, 3 or 5:
    if p % 8 in (3, 5):
        return ZZ(2)
    if p % 12 in (5, 7):
        return ZZ(3)
    if p % 5 in (2, 3):
        return ZZ(5)
    # default case (first needed for p=71):
    if not p.is_prime():
        raise ValueError("p must be a prime number > 2")
    from sage.arith.srange import xsrange
    for r in xsrange(7, p):
        if legendre_symbol(r, p) == -1:
            return ZZ(r)
