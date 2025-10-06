"""
Counting Congruence Solutions

This file provides more user-friendly Python front-ends to the Cython code in
:mod:`sage.quadratic_forms.count_local`.
"""
################################################################
# Methods for counting/computing the number of representations #
# of a number by a quadratic form in Z/(p^k)Z of various types #
################################################################

from sage.quadratic_forms.count_local_2 import CountAllLocalTypesNaive


def count_congruence_solutions_as_vector(self, p, k, m, zvec, nzvec):
    r"""
    Return the number of integer solution vectors `x` satisfying the
    congruence `Q(x) = m` (mod `p^k`) of each solution type (i.e. All,
    Good, Zero, Bad, BadI, BadII) which satisfy the additional
    congruence conditions of having certain coefficients = 0 (mod `p`)
    and certain collections of coefficients not congruent to the zero
    vector (mod `p`).

    A solution vector `x` satisfies the additional congruence conditions
    specified by ``zvec`` and ``nzvec`` (and therefore is counted) iff both of
    the following conditions hold:

    1. `x_i = 0` (mod `p`) for all `i` in ``zvec``

    2. `x_i \neq 0` (mod `p`) for all `i` in ``nzvec``


    REFERENCES:

    See Hanke's (????) paper "Local Densities and explicit bounds...", p??? for
    the definitions of the solution types and congruence conditions.

    INPUT:

    - ``p`` -- prime number > 0
    - ``k`` -- integer > 0
    - ``m`` -- integer (depending only on mod `p^k`)
    - ``zvec``, ``nzvec`` -- lists of integers in ``range(self.dim())``, or ``None``

    OUTPUT:

    a list of six integers `\geq 0` representing the solution types:
    [All, Good, Zero, Bad, BadI, BadII]

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions_as_vector(3, 1, 1, [], [])
        [0, 0, 0, 0, 0, 0]
        sage: Q.count_congruence_solutions_as_vector(3, 1, 1, None, [])
        [0, 0, 0, 0, 0, 0]
        sage: Q.count_congruence_solutions_as_vector(3, 1, 1, [], None)
        [6, 6, 0, 0, 0, 0]
        sage: Q.count_congruence_solutions_as_vector(3, 1, 1, None, None)
        [6, 6, 0, 0, 0, 0]
        sage: Q.count_congruence_solutions_as_vector(3, 1, 2, None, None)
        [6, 6, 0, 0, 0, 0]
        sage: Q.count_congruence_solutions_as_vector(3, 1, 0, None, None)
        [15, 12, 1, 2, 0, 2]
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)


# //////////////////////////////////////////
# // Front-ends for our counting routines //
# //////////////////////////////////////////

def count_congruence_solutions(self, p, k, m, zvec, nzvec):
    r"""
    Count all solutions of `Q(x) = m` (mod `p^k`) satisfying the
    additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers in ``range(self.dim())``, or ``None``

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions(3, 1, 0, None, None)
        15
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[0]


def count_congruence_solutions__good_type(self, p, k, m, zvec, nzvec):
    r"""
    Count the good-type solutions of `Q(x) = m` (mod `p^k`) satisfying the
    additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions__good_type(3, 1, 0, None, None)
        12
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[1]


def count_congruence_solutions__zero_type(self, p, k, m, zvec, nzvec):
    """
    Count the zero-type solutions of `Q(x) = m` (mod `p^k`) satisfying the
    additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions__zero_type(3, 1, 0, None, None)
        1
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[2]


def count_congruence_solutions__bad_type(self, p, k, m, zvec, nzvec):
    """
    Count the bad-type solutions of `Q(x) = m` (mod `p^k`) satisfying the
    additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions__bad_type(3, 1, 0, None, None)
        2
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[3]


def count_congruence_solutions__bad_type_I(self, p, k, m, zvec, nzvec):
    """
    Count the bad-typeI solutions of `Q(x) = m` (mod `p^k`) satisfying
    the additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions__bad_type_I(3, 1, 0, None, None)
        0
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[4]


def count_congruence_solutions__bad_type_II(self, p, k, m, zvec, nzvec):
    """
    Count the bad-typeII solutions of `Q(x) = m` (mod `p^k`) satisfying
    the additional congruence conditions described in
    :meth:`QuadraticForm.count_congruence_solutions_as_vector`.

    INPUT:

    - ``p`` -- prime number > 0

    - ``k`` -- integer > 0

    - ``m`` -- integer (depending only on mod `p^k`)

    - ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.count_congruence_solutions__bad_type_II(3, 1, 0, None, None)
        2
    """
    return CountAllLocalTypesNaive(self, p, k, m, zvec, nzvec)[5]
