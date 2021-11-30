"""
Creating A Random Quadratic Form
"""
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.quadratic_forms.ternary_qf import TernaryQF
from sage.rings.ring import is_Ring
from sage.rings.integer_ring import ZZ

################################################
## Routines to create a random quadratic form ##
################################################

def random_quadraticform(R, n, rand_arg_list=[]):
    r"""
    Create a random quadratic form in `n` variables defined over the ring `R`.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    INPUT:

    - `R` -- a ring.
    - `n` -- an integer `\ge 0`
    - ``rand_arg_list`` -- a list of at most 3 arguments which can be taken by
      ``R.random_element()``.

    OUTPUT:

    A quadratic form over the ring `R`.

    EXAMPLES::

        sage: random_quadraticform(ZZ, 3, [1,5])    ## RANDOM
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 2 3 ]
        [ * 1 4 ]
        [ * * 3 ]

    ::

        sage: random_quadraticform(ZZ, 3, [-5,5])    ## RANDOM
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 2 -5 ]
        [ * 2 -2 ]
        [ * * -5 ]

    ::

        sage: random_quadraticform(ZZ, 3, [-50,50])    ## RANDOM
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 8 -23 ]
        [ * 0 0 ]
        [ * * 6 ]
    """
    ## Sanity Checks: We have a ring and there are at most 3 parameters for randomness!
    if len(rand_arg_list) > 3:
        raise TypeError("Oops!  The list of randomness arguments can have at most 3 elements.")
    if not is_Ring(R):
        raise TypeError("Oops!  The first argument must be a ring.")

    ## Create a list of upper-triangular entries for the quadratic form
    L = len(rand_arg_list)
    nn = int(n*(n+1)/2)
    if L == 0:
        rand_list = [R.random_element()  for _ in range(nn)]
    elif L == 1:
        rand_list = [R.random_element(rand_arg_list[0])  for _ in range(nn)]
    elif L == 2:
        rand_list = [R.random_element(rand_arg_list[0], rand_arg_list[1])  for _ in range(nn)]
    elif L == 3:
        rand_list = [R.random_element(rand_arg_list[0], rand_arg_list[1], rand_arg_list[2])  for _ in range(nn)]

    ## Return  the Quadratic Form
    return QuadraticForm(R, n, rand_list)


def random_quadraticform_with_conditions(R, n, condition_list=[], rand_arg_list=[]):
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

        sage: Q = random_quadraticform_with_conditions(ZZ, 3, [QuadraticForm.is_positive_definite], [-5, 5])
        sage: Q    ## RANDOM
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 -2 -5 ]
        [ * 2 2 ]
        [ * * 3 ]

    """
    Q = random_quadraticform(R, n, rand_arg_list)
    Done_Flag = True

    ## Check that all conditions are satisfied
    while Done_Flag:
        Done_Flag = False
        for c in condition_list:

            ## Check if condition c is satisfied
            try:
                bool_ans = Q.c()
            except Exception:
                bool_ans = c(Q)

            ## Create a new quadratic form if a condition fails
            if not bool_ans:
                Q = random_quadraticform(R, n, rand_arg_list)
                Done_Flag = True
                break

    ## Return the quadratic form
    return Q

def random_ternaryqf(rand_arg_list = []):
    """
    Create a random ternary quadratic form.

    The last (and optional) argument ``rand_arg_list`` is a list of at most 3
    elements which is passed (as at most 3 separate variables) into the method
    ``R.random_element()``.

    INPUT:

    - ``rand_arg_list`` -- a list of at most 3 arguments which can be taken by
      ``R.random_element()``.

    OUTPUT:

    A ternary quadratic form.

    EXAMPLES::

        sage: random_ternaryqf()  ##RANDOM
        Ternary quadratic form with integer coefficients:
        [1 1 4]
        [-1 1 -1]
        sage: random_ternaryqf([-1, 2])  ##RANDOM
        Ternary quadratic form with integer coefficients:
        [1 0 1]
        [-1 -1 -1]
        sage: random_ternaryqf([-10, 10, "uniform"])  ##RANDOM
        Ternary quadratic form with integer coefficients:
        [7 -8 2]
        [0 3 -6]
    """


    R = ZZ
    n = 6
    L = len(rand_arg_list)
    if L == 0:
        rand_list = [ R.random_element() for _ in range(n)]
    elif L == 1:
        rand_list = [ R.random_element(rand_arg_list[0]) for _ in range(6)]
    elif L == 2:
        rand_list = [ R.random_element(rand_arg_list[0], rand_arg_list[1]) for _ in range(6)]
    elif L == 3:
        rand_list = [ R.random_element(rand_arg_list[0], rand_arg_list[1], rand_arg_list[2]) for _ in range(6)]

    return TernaryQF(rand_list)



def random_ternaryqf_with_conditions(condition_list=[], rand_arg_list=[]):
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

        sage: Q = random_ternaryqf_with_conditions([TernaryQF.is_positive_definite], [-5, 5])
        sage: Q    ## RANDOM
        Ternary quadratic form with integer coefficients:
        [3 4 2]
        [2 -2 -1]
    """

    Q = random_ternaryqf(rand_arg_list)
    Done_Flag = True

    ## Check that all conditions are satisfied
    while Done_Flag:
        Done_Flag = False
        for c in condition_list:

            ## Check if condition c is satisfied
            try:
                bool_ans = Q.c()
            except Exception:
                bool_ans = c(Q)

            ## Create a new quadratic form if a condition fails
            if not bool_ans:
                Q = random_ternaryqf(rand_arg_list)
                Done_Flag = True
                break
    return Q
