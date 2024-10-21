# sage_setup: distribution = sagemath-modules
"""
Random quadratic forms

This file contains a set of routines to create a random quadratic form.
"""
from sage.categories.rings import Rings
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.quadratic_forms.ternary_qf import TernaryQF
from sage.rings.integer_ring import ZZ


################################################
# Routines to create a random quadratic form ##
################################################

def random_quadraticform(R, n, rand_arg_list=None):
    r"""
    Create a random quadratic form in `n` variables defined over the ring `R`.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    INPUT:

    - ``R`` -- a ring
    - ``n`` -- integer `\ge 0`
    - ``rand_arg_list`` -- list of at most 3 arguments which can be taken by
      ``R.random_element()``

    OUTPUT: a quadratic form over the ring `R`

    EXAMPLES::

        sage: random_quadraticform(ZZ, 3, [1,5])    # random
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 2 3 ]
        [ * 1 4 ]
        [ * * 3 ]

        sage: random_quadraticform(ZZ, 3, [-5,5])    # random
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 2 -5 ]
        [ * 2 -2 ]
        [ * * -5 ]

        sage: random_quadraticform(ZZ, 3, [-50,50])    # random
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 8 -23 ]
        [ * 0 0 ]
        [ * * 6 ]

    TESTS::

        sage: random_quadraticform(ZZ, 3, [1,2,3,4])
        Traceback (most recent call last):
        ...
        TypeError: the list of randomness arguments can have at most 3 elements
    """
    if rand_arg_list is None:
        rand_arg_list = []
    if len(rand_arg_list) > 3:
        raise TypeError("the list of randomness arguments can have "
                        "at most 3 elements")
    if R not in Rings():
        raise TypeError("the first argument must be a ring")
    # Create a list of upper-triangular entries for the quadratic form
    n2 = (n * (n + 1)) // 2
    if not rand_arg_list:
        rand_list = [R.random_element() for _ in range(n2)]
    else:
        rand_list = [R.random_element(*rand_arg_list) for _ in range(n2)]
    return QuadraticForm(R, n, rand_list)


def random_quadraticform_with_conditions(R, n, condition_list=[],
                                         rand_arg_list=None):
    """
    Create a random quadratic form in `n` variables defined over the ring `R`
    satisfying a list of boolean (i.e. True/False) conditions.

    The conditions `c` appearing in the list must be boolean functions which
    can be called either as ``Q.c()`` or ``c(Q)``, where ``Q`` is the random
    quadratic form.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    EXAMPLES::

        sage: check = QuadraticForm.is_positive_definite
        sage: Q = random_quadraticform_with_conditions(ZZ, 3, [check], [-5, 5])
        sage: Q    # random
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 -2 -5 ]
        [ * 2 2 ]
        [ * * 3 ]
    """
    if rand_arg_list is None:
        rand_arg_list = []

    Q = random_quadraticform(R, n, rand_arg_list)
    done_flag = True

    # Check that all conditions are satisfied
    while done_flag:
        done_flag = False
        for c in condition_list:

            # Check if condition c is satisfied
            try:
                bool_ans = Q.c()
            except Exception:
                bool_ans = c(Q)

            # Create a new quadratic form if a condition fails
            if not bool_ans:
                Q = random_quadraticform(R, n, rand_arg_list)
                done_flag = True
                break

    # Return the quadratic form
    return Q


def random_ternaryqf(rand_arg_list=None):
    """
    Create a random ternary quadratic form.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    INPUT:

    - ``rand_arg_list`` -- list of at most 3 arguments which can be taken by
      ``R.random_element()``

    OUTPUT: a ternary quadratic form

    EXAMPLES::

        sage: random_ternaryqf()  # random
        Ternary quadratic form with integer coefficients:
        [1 1 4]
        [-1 1 -1]
        sage: random_ternaryqf([-1, 2])  # random
        Ternary quadratic form with integer coefficients:
        [1 0 1]
        [-1 -1 -1]
        sage: random_ternaryqf([-10, 10, "uniform"])  # random
        Ternary quadratic form with integer coefficients:
        [7 -8 2]
        [0 3 -6]
    """
    if rand_arg_list is None:
        rand_arg_list = []
    if not rand_arg_list:
        rand_list = [ZZ.random_element() for _ in range(6)]
    else:
        rand_list = [ZZ.random_element(*rand_arg_list) for _ in range(6)]
    return TernaryQF(rand_list)


def random_ternaryqf_with_conditions(condition_list=[], rand_arg_list=None):
    """
    Create a random ternary quadratic form satisfying a list of boolean
    (i.e. True/False) conditions.

    The conditions `c` appearing in the list must be boolean functions which
    can be called either as ``Q.c()`` or ``c(Q)``, where ``Q`` is the random
    ternary quadratic form.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    EXAMPLES::

        sage: check = TernaryQF.is_positive_definite
        sage: Q = random_ternaryqf_with_conditions([check], [-5, 5])
        sage: Q    # random
        Ternary quadratic form with integer coefficients:
        [3 4 2]
        [2 -2 -1]
    """
    if rand_arg_list is None:
        rand_arg_list = []
    Q = random_ternaryqf(rand_arg_list)
    done_flag = True

    # Check that all conditions are satisfied
    while done_flag:
        done_flag = False
        for c in condition_list:
            # Check if condition c is satisfied
            try:
                bool_ans = Q.c()
            except Exception:
                bool_ans = c(Q)

            # Create a new quadratic form if a condition fails
            if not bool_ans:
                Q = random_ternaryqf(rand_arg_list)
                done_flag = True
                break
    return Q
