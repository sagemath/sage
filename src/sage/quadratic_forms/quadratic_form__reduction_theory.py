"""
Reduction Theory
"""
from copy import deepcopy
from sage.matrix.constructor import matrix
from sage.functions.all import floor
from sage.misc.mrange import mrange
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ


def reduced_binary_form1(self):
    r"""
    Reduce the form `ax^2 + bxy+cy^2` to satisfy the reduced condition `|b| \le
    a \le c`, with `b \ge 0` if `a = c`. This reduction occurs within the
    proper class, so all transformations are taken to have determinant 1.

    EXAMPLES::

        sage: QuadraticForm(ZZ,2,[5,5,2]).reduced_binary_form1()
        (
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 2 -1 ]
        [ * 2 ]                                                            ,
        <BLANKLINE>
        [ 0 -1]
        [ 1  1]
        )
    """
    if self.dim() != 2:
        raise TypeError("This must be a binary form for now...")

    R = self.base_ring()
    interior_reduced_flag = False
    Q = deepcopy(self)
    M = matrix(R, 2, 2, [1,0,0,1])

    while not interior_reduced_flag:
        interior_reduced_flag = True

        ## Arrange for a <= c
        if Q[0,0] > Q[1,1]:
            M_new = matrix(R,2,2,[0, -1, 1, 0])
            Q = Q(M_new)
            M = M * M_new
            interior_reduced_flag = False

        ## Arrange for |b| <= a
        if abs(Q[0,1]) > Q[0,0]:
            r = R(floor(round(Q[0,1]/(2*Q[0,0]))))
            M_new = matrix(R,2,2,[1, -r, 0, 1])
            Q = Q(M_new)
            M = M * M_new
            interior_reduced_flag = False

    return Q, M




def reduced_ternary_form__Dickson(self):
    """
    Find the unique reduced ternary form according to the conditions
    of Dickson's "Studies in the Theory of Numbers", pp164-171.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1, 1, 1])
        sage: Q.reduced_ternary_form__Dickson()
        Traceback (most recent call last):
        ...
        NotImplementedError: TO DO

    """
    raise NotImplementedError("TO DO")



def reduced_binary_form(self):
    """
    Find a form which is reduced in the sense that no further binary
    form reductions can be done to reduce the original form.

    EXAMPLES::

        sage: QuadraticForm(ZZ,2,[5,5,2]).reduced_binary_form()
        (
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 2 -1 ]
        [ * 2 ]                                                            ,
        <BLANKLINE>
        [ 0 -1]
        [ 1  1]
        )
    """
    R = self.base_ring()
    n = self.dim()
    interior_reduced_flag = False
    Q = deepcopy(self)
    M = matrix(R, n, n)
    for i in range(n):
        M[i,i] = 1


    while not interior_reduced_flag:
        interior_reduced_flag = True

        ## Arrange for (weakly) increasing diagonal entries
        for i in range(n):
            for j in range(i+1,n):
                if Q[i,i] > Q[j,j]:
                    M_new = matrix(R,n,n)
                    for k in range(n):
                        M_new[k,k] = 1
                    M_new[i,j] = -1
                    M_new[j,i] = 1
                    M_new[i,i] = 0
                    M_new[j,j] = 1

                    Q = Q(M_new)
                    M = M * M_new
                    interior_reduced_flag = False

                ## Arrange for |b| <= a
                if abs(Q[i,j]) > Q[i,i]:
                    r = R(floor(round(Q[i,j]/(2*Q[i,i]))))

                    M_new = matrix(R,n,n)
                    for k in range(n):
                        M_new[k,k] = 1
                    M_new[i,j] = -r

                    Q = Q(M_new)
                    M = M * M_new
                    interior_reduced_flag = False

    return Q, M


def minkowski_reduction(self):
    """
    Find a Minkowski-reduced form equivalent to the given one.
    This means that

    .. MATH::

            Q(v_k) <= Q(s_1 * v_1 + ... + s_n * v_n)

    for all `s_i` where GCD`(s_k, ... s_n) = 1`.

    Note: When Q has dim <= 4 we can take all `s_i` in {1, 0, -1}.

    References:
        Schulze-Pillot's paper on "An algorithm for computing genera
            of ternary and quaternary quadratic forms", p138.
        Donaldson's 1979 paper "Minkowski Reduction of Integral
            Matrices", p203.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ,4,[30, 17, 11, 12, 29, 25, 62, 64, 25, 110])
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 30 17 11 12 ]
        [ * 29 25 62 ]
        [ * * 64 25 ]
        [ * * * 110 ]
        sage: Q.minkowski_reduction()
        (
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 30 17 11 -5 ]
        [ * 29 25 4 ]
        [ * * 64 0 ]
        [ * * * 77 ]                                                       ,
        <BLANKLINE>
        [ 1  0  0  0]
        [ 0  1  0 -1]
        [ 0  0  1  0]
        [ 0  0  0  1]
        )

    ::

        sage: Q=QuadraticForm(ZZ,4,[1, -2, 0, 0, 2, 0, 0, 2, 0, 2])
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 -2 0 0 ]
        [ * 2 0 0 ]
        [ * * 2 0 ]
        [ * * * 2 ]
        sage: Q.minkowski_reduction()
        (
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 0 0 ]
        [ * 1 0 0 ]
        [ * * 2 0 ]
        [ * * * 2 ]                                                        ,
        <BLANKLINE>
        [1 1 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        )

    ::

        sage: Q=QuadraticForm(ZZ,5,[2,2,0,0,0,2,2,0,0,2,2,0,2,2,2])
        sage: Q.Gram_matrix()
        [2 1 0 0 0]
        [1 2 1 0 0]
        [0 1 2 1 0]
        [0 0 1 2 1]
        [0 0 0 1 2]
        sage: Q.minkowski_reduction()
        Traceback (most recent call last):
        ...
        NotImplementedError: This algorithm is only for dimensions less than 5
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm
    from sage.quadratic_forms.quadratic_form import matrix
    if not self.is_positive_definite():
        raise TypeError("Minkowski reduction only works for positive definite forms")
    if self.dim() > 4:
        raise NotImplementedError("This algorithm is only for dimensions less than 5")

    R = self.base_ring()
    n = self.dim()
    Q = deepcopy(self)
    M = matrix(R, n, n)
    for i in range(n):
        M[i, i] = 1

    ## Begin the reduction
    done_flag = False
    while not done_flag:

        ## Loop through possible shorted vectors until
        done_flag = True
        for j in range(n-1, -1, -1):
            for a_first in mrange([3  for i in range(j)]):
                y = [x-1 for x in a_first] + [1] + [0 for k in range(n-1-j)]
                e_j = [0  for k in range(n)]
                e_j[j] = 1

                ## Reduce if a shorter vector is found
                if Q(y) < Q(e_j):

                    ## Create the transformation matrix
                    M_new = matrix(R, n, n)
                    for k in range(n):
                        M_new[k,k] = 1
                    for k in range(n):
                        M_new[k,j] = y[k]

                    ## Perform the reduction and restart the loop
                    Q = QuadraticForm(M_new.transpose()*Q.matrix()*M_new)
                    M = M * M_new
                    done_flag = False

                if not done_flag:
                    break

            if not done_flag:
                break

    ## Return the results
    return Q, M




def minkowski_reduction_for_4vars__SP(self):
    """
    Find a Minkowski-reduced form equivalent to the given one.
    This means that

        Q(`v_k`) <= Q(`s_1 * v_1 + ... + s_n * v_n`)

    for all `s_i` where GCD(`s_k, ... s_n`) = 1.

    Note: When Q has dim <= 4 we can take all `s_i` in {1, 0, -1}.

    References:
        Schulze-Pillot's paper on "An algorithm for computing genera
            of ternary and quaternary quadratic forms", p138.
        Donaldson's 1979 paper "Minkowski Reduction of Integral
            Matrices", p203.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ,4,[30,17,11,12,29,25,62,64,25,110])
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 30 17 11 12 ]
        [ * 29 25 62 ]
        [ * * 64 25 ]
        [ * * * 110 ]
        sage: Q.minkowski_reduction_for_4vars__SP()
        (
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 29 -17 25 4 ]
        [ * 30 -11 5 ]
        [ * * 64 0 ]
        [ * * * 77 ]                                                       ,
        <BLANKLINE>
        [ 0  1  0  0]
        [ 1  0  0 -1]
        [ 0  0  1  0]
        [ 0  0  0  1]
        )
    """
    R = self.base_ring()
    n = self.dim()
    Q = deepcopy(self)
    M = matrix(R, n, n)
    for i in range(n):
        M[i, i] = 1

    ## Only allow 4-variable forms
    if n != 4:
        raise TypeError("Oops!  The given quadratic form has " + str(n) +  \
                " != 4 variables. =|")


    ## Step 1: Begin the reduction
    done_flag = False
    while not done_flag:

        ## Loop through possible shorter vectors
        done_flag = True
        for j in range(n-1, -1, -1):
            for a_first in mrange([2  for i in range(j)]):
                y = [x-1 for x in a_first] + [1] + [0 for k in range(n-1-j)]
                e_j = [0  for k in range(n)]
                e_j[j] = 1

                ## Reduce if a shorter vector is found
                if Q(y) < Q(e_j):

                    ## Further n=4 computations
                    B_y_vec = Q.matrix() * vector(ZZ, y)
                        ## SP's B = our self.matrix()/2
                        ## SP's A = coeff matrix of his B
                        ## Here we compute the double of both and compare.
                    B_sum = sum([abs(B_y_vec[i])  for i in range(4)  if i != j])
                    A_sum = sum([abs(Q[i,j])  for i in range(4)  if i != j])
                    B_max = max([abs(B_y_vec[i])  for i in range(4)  if i != j])
                    A_max = max([abs(Q[i,j])  for i in range(4)  if i != j])

                    if (B_sum < A_sum) or ((B_sum == A_sum) and (B_max < A_max)):

                        ## Create the transformation matrix
                        M_new = matrix(R, n, n)
                        for k in range(n):
                            M_new[k,k] = 1
                        for k in range(n):
                            M_new[k,j] = y[k]

                        ## Perform the reduction and restart the loop
                        Q = Q(M_new)
                        M = M * M_new
                        done_flag = False

                if not done_flag:
                    break

            if not done_flag:
                break

    ## Step 2: Order A by certain criteria
    for i in range(4):
        for j in range(i+1,4):

            ## Condition (a)
            if (Q[i,i] > Q[j,j]):
                Q.swap_variables(i,j,in_place=True)
                M_new = matrix(R,n,n)
                M_new[i,j] = -1
                M_new[j,i] = 1
                for r in range(4):
                    if (r == i) or (r == j):
                        M_new[r,r] = 0
                    else:
                        M_new[r,r] = 1
                M = M * M_new

            elif (Q[i,i] == Q[j,j]):
                i_sum = sum([abs(Q[i,k])  for k in range(4)  if k != i])
                j_sum = sum([abs(Q[j,k])  for k in range(4)  if k != j])

                ## Condition (b)
                if (i_sum > j_sum):
                    Q.swap_variables(i,j,in_place=True)
                    M_new = matrix(R,n,n)
                    M_new[i,j] = -1
                    M_new[j,i] = 1
                    for r in range(4):
                        if (r == i) or (r == j):
                            M_new[r,r] = 0
                        else:
                            M_new[r,r] = 1
                    M = M * M_new

                elif (i_sum == j_sum):
                    for k in [2,1,0]:   ## TO DO: These steps are a little redundant...
                        Q1 = Q.matrix()

                        c_flag = True
                        for l in range(k+1,4):
                            c_flag = c_flag and (abs(Q1[i,l]) == abs(Q1[j,l]))

                        ## Condition (c)
                        if c_flag and (abs(Q1[i,k]) > abs(Q1[j,k])):
                            Q.swap_variables(i,j,in_place=True)
                            M_new = matrix(R,n,n)
                            M_new[i,j] = -1
                            M_new[j,i] = 1
                            for r in range(4):
                                if (r == i) or (r == j):
                                    M_new[r,r] = 0
                                else:
                                    M_new[r,r] = 1
                            M = M * M_new


    ## Step 3: Order the signs
    for i in range(4):
        if Q[i,3] < 0:
            Q.multiply_variable(-1, i, in_place=True)
            M_new = matrix(R,n,n)
            for r in range(4):
                if r == i:
                    M_new[r,r] = -1
                else:
                    M_new[r,r] = 1
            M = M * M_new

    for i in range(4):
        j = 3
        while (Q[i,j] == 0):
            j += -1
        if (Q[i,j] < 0):
            Q.multiply_variable(-1, i, in_place=True)
            M_new = matrix(R,n,n)
            for r in range(4):
                if r == i:
                    M_new[r,r] = -1
                else:
                    M_new[r,r] = 1
            M = M * M_new

    if Q[1,2] < 0:
        ## Test a row 1 sign change
        if (Q[1,3] <= 0 and \
            ((Q[1,3] < 0) or (Q[1,3] == 0 and Q[1,2] < 0)  \
                or (Q[1,3] == 0 and Q[1,2] == 0 and Q[1,1] < 0))):
            Q.multiply_variable(-1, i, in_place=True)
            M_new = matrix(R,n,n)
            for r in range(4):
                if r == i:
                    M_new[r,r] = -1
                else:
                    M_new[r,r] = 1
            M = M * M_new

        elif (Q[2,3] <= 0 and \
            ((Q[2,3] < 0) or (Q[2,3] == 0 and Q[2,2] < 0)  \
                or (Q[2,3] == 0 and Q[2,2] == 0 and Q[2,1] < 0))):
            Q.multiply_variable(-1, i, in_place=True)
            M_new = matrix(R,n,n)
            for r in range(4):
                if r == i:
                    M_new[r,r] = -1
                else:
                    M_new[r,r] = 1
            M = M * M_new


    ## Return the results
    return Q, M


