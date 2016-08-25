r"""
Optimised Cython code for counting congruence solutions
"""
from __future__ import print_function

from sage.arith.all import valuation, kronecker_symbol, is_prime
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

from sage.rings.integer_ring import ZZ

from sage.rings.finite_rings.integer_mod cimport IntegerMod_gmp
from sage.sets.set import Set




def extract_sublist_indices(Biglist, Smalllist):
    """
    Returns the indices of Biglist which index the entries of
    Smalllist appearing in Biglist.  (Note that Smalllist may not be a
    sublist of Biglist.)

    NOTE 1: This is an internal routine which deals with re-indexing
    lists, and is not exported to the QuadraticForm namespace!

    NOTE 2: This should really by applied only when BigList has no
    repeated entries.

    TO DO: *** Please revisit this routine, and eliminate it! ***

    INPUT:
        Biglist, Smalllist -- two lists of a common type, where Biglist has no
        repeated entries.

    OUTPUT:
        a list of integers >= 0

    EXAMPLES::

        sage: from sage.quadratic_forms.count_local_2 import extract_sublist_indices

        sage: biglist = [1,3,5,7,8,2,4]
        sage: sublist = [5,3,2]
        sage: sublist == [biglist[i]  for i in extract_sublist_indices(biglist, sublist)]  ## Ok whenever Smalllist is a sublist of Biglist
        True

        sage: extract_sublist_indices([1,2,3,6,9,11], [1,3,2,9])
        [0, 2, 1, 4]

        sage: extract_sublist_indices([1,2,3,6,9,11], [1,3,10,2,9,0])
        [0, 2, 1, 4]

        sage: extract_sublist_indices([1,3,5,3,8], [1,5])
        Traceback (most recent call last):
        ...
        TypeError: Biglist must not have repeated entries!
    """
    ## Check that Biglist has no repeated entries
    Big_set = Set(Biglist)
    if len(Set(Biglist)) != len(Biglist):
        raise TypeError("Biglist must not have repeated entries!")

    ## Extract the indices of Biglist needed to make Sublist
    index_list = []
    for x in Smalllist:
        try:
            index_list.append(Biglist.index(x))
        except ValueError:                                   ## This happens when an entry of Smalllist is not contained in Biglist
            None

    ## Return the list if indices
    return index_list




def count_modp__by_gauss_sum(n, p, m, Qdet):
    """
    Returns the number of solutions of Q(x) = m over the finite field
    Z/pZ, where p is a prime number > 2 and Q is a non-degenerate
    quadratic form of dimension n >= 1 and has Gram determinant Qdet.

    REFERENCE:
        These are defined in Table 1 on p363 of Hanke's "Local
        Densities..." paper.

    INPUT:

    - n -- an integer >= 1
    - p -- a prime number > 2
    - m -- an integer
    - Qdet -- a integer which is non-zero mod p

    OUTPUT:
        an integer >= 0

    EXAMPLES::

        sage: from sage.quadratic_forms.count_local_2 import count_modp__by_gauss_sum

        sage: count_modp__by_gauss_sum(3, 3, 0, 1)    ## for Q = x^2 + y^2 + z^2  => Gram Det = 1 (mod 3)
        9
        sage: count_modp__by_gauss_sum(3, 3, 1, 1)    ## for Q = x^2 + y^2 + z^2  => Gram Det = 1 (mod 3)
        6
        sage: count_modp__by_gauss_sum(3, 3, 2, 1)    ## for Q = x^2 + y^2 + z^2  => Gram Det = 1 (mod 3)
        12

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: [Q.count_congruence_solutions(3, 1, m, None, None) == count_modp__by_gauss_sum(3, 3, m, 1)  for m in range(3)]
        [True, True, True]


        sage: count_modp__by_gauss_sum(3, 3, 0, 2)    ## for Q = x^2 + y^2 + 2*z^2  => Gram Det = 2 (mod 3)
        9
        sage: count_modp__by_gauss_sum(3, 3, 1, 2)    ## for Q = x^2 + y^2 + 2*z^2  => Gram Det = 2 (mod 3)
        12
        sage: count_modp__by_gauss_sum(3, 3, 2, 2)    ## for Q = x^2 + y^2 + 2*z^2  => Gram Det = 2 (mod 3)
        6

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,2])
        sage: [Q.count_congruence_solutions(3, 1, m, None, None) == count_modp__by_gauss_sum(3, 3, m, 2)  for m in range(3)]
        [True, True, True]


    """
    ## Check that Qdet is non-degenerate
    if Qdet % p == 0:
        raise RuntimeError("Qdet must be non-zero.")

    ## Check that p is prime > 0
    if not is_prime(p) or p == 2:
        raise RuntimeError("p must be a prime number > 2.")

    ## Check that n >= 1
    if n < 1:
        raise RuntimeError("the dimension n must be >= 1.")

    ## Compute the Gauss sum
    neg1 = -1
    if (m % p == 0):
        if (n % 2 != 0):
            count = (p**(n-1))
        else:
            count = (p**(n-1)) + (p-1) * (p**((n-2)/2)) * kronecker_symbol(((neg1**(n/2)) * Qdet) % p, p)
    else:
        if (n % 2 != 0):
            count = (p**(n-1)) + (p**((n-1)/2)) * kronecker_symbol(((neg1**((n-1)/2)) * Qdet * m) % p, p)
        else:
            count = (p**(n-1)) - (p**((n-2)/2)) * kronecker_symbol(((neg1**(n/2)) * Qdet) % p, p)

    ## Return the result
    return count






cdef CountAllLocalTypesNaive_cdef(Q, p, k, m, zvec, nzvec):
    """
    This Cython routine is documented in its Python wrapper method
    QuadraticForm.count_congruence_solutions_by_type().
    """
    cdef long n, i
    cdef long a, b    ## Used to quickly evaluate Q(v)
    cdef long ptr     ## Used to increment the vector
    cdef long solntype    ## Used to store the kind of solution we find


    ## Some shortcuts and definitions
    n = Q.dim()
    R = p ** k
    Q1 = Q.base_change_to(IntegerModRing(R))


    ## Cython Variables
    cdef IntegerMod_gmp zero, one
    zero = IntegerMod_gmp(IntegerModRing(R), 0)
    one = IntegerMod_gmp(IntegerModRing(R), 1)



    ## Initialize the counting vector
    count_vector = [0  for i in range(6)]

    ## Initialize v = (0, ... , 0)
    v = [Mod(0, R)  for i in range(n)]


    ## Some declarations to speed up the loop
    R_n = R ** n
    m1 = Mod(m, R)

    ## Count the local solutions
    for i from 0 <= i < R_n:

        ## Perform a carry (when value = R-1) until we can increment freely
        ptr = len(v)
        while ((ptr > 0) and (v[ptr-1] == R-1)):
            v[ptr-1] += 1
            ptr += -1

        ## Only increment if we're not already at the zero vector =)
        if (ptr > 0):
            v[ptr-1] += 1


        ## Evaluate Q(v) quickly
        tmp_val = Mod(0, R)
        for a from 0 <= a < n:
            for b from a <= b < n:
                tmp_val += Q1[a,b] * v[a] * v[b]

        ## Sort the solution by it's type
        #if (Q1(v) == m1):
        if (tmp_val == m1):
            solntype = local_solution_type_cdef(Q1, p, v, zvec, nzvec)
            if (solntype != 0):
                count_vector[solntype] += 1


    ## Generate the Bad-type and Total counts
    count_vector[3] = count_vector[4] + count_vector[5]
    count_vector[0] = count_vector[1] + count_vector[2] + count_vector[3]

    ## Return the solution counts
    return count_vector




def CountAllLocalTypesNaive(Q, p, k, m, zvec, nzvec):
    """
    This is an internal routine, which is called by
    :meth:`sage.quadratic_forms.quadratic_form.QuadraticForm.count_congruence_solutions_by_type
    QuadraticForm.count_congruence_solutions_by_type`. See the documentation of
    that method for more details.

    INPUT:

    - `Q` -- quadratic form over `\ZZ`
    - `p` -- prime number > 0
    - `k` -- an integer > 0
    - `m` -- an integer (depending only on mod `p^k`)
    - ``zvec``, ``nzvec`` -- a list of integers in ``range(Q.dim())``, or ``None``

    OUTPUT:
        a list of six integers `\ge 0` representing the solution types: ``[All,
        Good, Zero, Bad, BadI, BadII]``


    EXAMPLES::

        sage: from sage.quadratic_forms.count_local_2 import CountAllLocalTypesNaive
        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: CountAllLocalTypesNaive(Q, 3, 1, 1, None, None)
        [6, 6, 0, 0, 0, 0]
        sage: CountAllLocalTypesNaive(Q, 3, 1, 2, None, None)
        [6, 6, 0, 0, 0, 0]
        sage: CountAllLocalTypesNaive(Q, 3, 1, 0, None, None)
        [15, 12, 1, 2, 0, 2]

    """
    return CountAllLocalTypesNaive_cdef(Q, p, k, m, zvec, nzvec)







cdef local_solution_type_cdef(Q, p, w, zvec, nzvec):
    """
    Internal routine to check if a given solution vector w (of Q(w) =
    m mod p^k) is of a certain local type and satisfies certain
    congruence conditions mod p.

    NOTE: No internal checking is done to test if p is a prime >=2, or
    that Q has the same size as w.

    """
    cdef long i
    cdef long n

    n = Q.dim()


    ## Check if the solution satisfies the zvec "zero" congruence conditions
    ## (either zvec is empty or its components index the zero vector mod p)
    if (zvec is None) or (len(zvec) == 0):
        zero_flag = True
    else:
        zero_flag = False
        i = 0
        while ( (i < len(zvec)) and ((w[zvec[i]] % p) == 0) ):  ## Increment so long as our entry is zero (mod p)
            i += 1
        if (i == len(zvec)):      ## If we make it through all entries then the solution is zero (mod p)
            zero_flag = True


    ## DIAGNOSTIC
    #print("IsLocalSolutionType: Finished the Zero congruence condition test \n")

    if (zero_flag is False):
        return <long> 0

    ## DIAGNOSTIC
    #print("IsLocalSolutionType: Passed the Zero congruence condition test \n")

    ## Check if the solution satisfies the nzvec "nonzero" congruence conditions
    ## (nzvec is non-empty and its components index a non-zero vector mod p)
    if (nzvec is None):
        nonzero_flag = True
    elif (len(nzvec) == 0):
        nonzero_flag = False           ## Trivially no solutions in this case!
    else:
        nonzero_flag = False
        i = 0
        while ((nonzero_flag == False) and (i < len(nzvec))):
            if ((w[nzvec[i]] % p) != 0):
                nonzero_flag = True           ## The non-zero condition is satisfied when we find one non-zero entry
            i += 1

    if (nonzero_flag == False):
        return <long> 0


    ## Check if the solution has the appropriate (local) type:
    ## -------------------------------------------------------

    ## 1: Check Good-type
    for i from 0 <= i < n:
        if (((w[i] % p) != 0)  and ((Q[i,i] % p) != 0)):
            return <long> 1
    if (p == 2):
        for i from 0 <= i < (n - 1):
            if (((Q[i,i+1] % p) != 0) and (((w[i] % p) != 0) or ((w[i+1] % p) != 0))):
                return <long> 1


    ## 2: Check Zero-type
    Zero_flag = True
    for i from 0 <= i < n:
        if ((w[i] % p) != 0):
            Zero_flag = False
    if (Zero_flag == True):
        return <long> 2


    ## Check if wS1 is zero or not
    wS1_nonzero_flag = False
    for i from 0 <= i < n:

        ## Compute the valuation of each index, allowing for off-diagonal terms
        if (Q[i,i] == 0):
            if (i == 0):
                val = valuation(Q[i,i+1], p)    ## Look at the term to the right
            elif (i == n - 1):
                val = valuation(Q[i-1,i], p)    ## Look at the term above
            else:
                val = valuation(Q[i,i+1] + Q[i-1,i], p)    ## Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(Q[i,i], p)

        ## Test each index
        if ((val == 1) and ((w[i] % p) != 0)):
            wS1_nonzero_flag = True

    ## 4: Check Bad-type I
    if (wS1_nonzero_flag is True):
        return <long> 4

    ## 5: Check Bad-type II
    if (wS1_nonzero_flag is False):
        return <long> 5

    ## Error if we get here! =o
    print("   Solution vector is " + str(w))
    print("   and Q is \n" + str(Q) + "\n")
    raise RuntimeError("Error in IsLocalSolutionType: Should not execute this line... =( \n")
