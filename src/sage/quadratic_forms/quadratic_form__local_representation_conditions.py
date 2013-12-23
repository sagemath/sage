"""
Local Representation Conditions
"""
##########################################################################
## Class for keeping track of the local conditions for representability ##
## of numbers by a quadratic form over ZZ (and eventually QQ also).     ##
##########################################################################

from copy import deepcopy

from sage.rings.integer_ring import ZZ
from sage.rings.arith import prime_divisors, valuation, is_square
from sage.quadratic_forms.extras import least_quadratic_nonresidue
from sage.rings.infinity import infinity
from sage.misc.functional import numerator, denominator
from sage.rings.rational_field import QQ



class QuadraticFormLocalRepresentationConditions():
    """
    Creates a class for dealing with the local conditions of a
    quadratic form, and checking local representability of numbers.

    EXAMPLES::

        sage: Q4 = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q4.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [].  For these and the reals, we have:
             Reals:   [0, +Infinity]
        sage: Q4.is_locally_represented_number(1)
        True
        sage: Q4.is_locally_universal_at_all_primes()
        True
        sage: Q4.is_locally_universal_at_all_places()
        False
        sage: L = [m  for m in range(-5, 100)  if Q4.is_locally_represented_number(m)]
        sage: L == range(100)
        True

    ::

        sage: Q3 = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q3.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [2].  For these and the reals, we have:
             Reals:   [0, +Infinity]
             p = 2:   [0, 0, 0, +Infinity, 0, 0, 0, 0]
        sage: E = [m  for m in range(100)  if not Q3.is_locally_represented_number(m)]
        sage: E1 = [m  for m in range(1,100)  if m / 2**(2*floor(valuation(m,2)/2)) % 8 == 7]
        sage: E == E1
        True
        sage: E
        [7, 15, 23, 28, 31, 39, 47, 55, 60, 63, 71, 79, 87, 92, 95]

    ::

        sage: Q2 = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q2.local_representation_conditions()
        This 2-dimensional form represents the p-adic integers of even
        valuation for all primes p except [2].
        For these and the reals, we have:
             Reals:   [0, +Infinity]
             p = 2:   [0, +Infinity, 0, +Infinity, 0, +Infinity, 0, +Infinity]
        sage: Q2.is_locally_universal_at_all_places()
        False
        sage: Q2.is_locally_universal_at_all_primes()
        False
        sage: L = [m  for m in range(-5, 25)  if Q2.is_locally_represented_number(m)]
        sage: L1 = [0] + [m  for m in range(1,25)  \
              if len([p  for p in prime_factors(squarefree_part(ZZ(m)))  if (p % 4) == 3]) % 2 == 0]
        sage: L == L1
        True
        sage: L
        [0, 1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20, 21]

    ::

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1])
        sage: Q1.local_representation_conditions()
        This 1-dimensional form only represents square multiples of 1.
        sage: L = [m  for m in range(100)  if Q1.is_locally_represented_number(m)]
        sage: L
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    ::

        sage: Q0 = DiagonalQuadraticForm(ZZ, [])
        sage: Q0.local_representation_conditions()
        This 0-dimensional form only represents zero.
        sage: L = [m  for m in range(100)  if Q0.is_locally_represented_number(m)]
        sage: L
        [0]

    """

    def __init__(self, Q):
        """
        Takes a QuadraticForm and computes its local conditions (if
        they don't already exist).  The recompute_flag overrides the
        previously computed conditions if they exist, and stores the
        new conditions.

        INPUT:
            Q -- Quadratic form over ZZ

        OUTPUT:
            a  QuadraticFormLocalRepresentationConditions object

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: QuadraticFormLocalRepresentationConditions(Q)
            This form represents the p-adic integers Z_p for all primes p except
            [].  For these and the reals, we have:
                 Reals:   [0, +Infinity]

        """

        ## Check that the form Q is integer-valued (we can relax this later)
        if Q.base_ring() != ZZ:
            raise TypeError, "We require that the quadratic form be defined over ZZ (integer-values) for now."


        ## Basic structure initialization
        self.local_repn_array = []    ## List of all local conditions
        self.dim = Q.dim()       ## We allow this to be any non-negative integer.
        self.exceptional_primes = [infinity]

        ## Deal with the special cases of 0 and 1-dimensional forms
        if self.dim == 0:
            self.coeff = None
            return
        elif self.dim == 1:
            self.coeff = Q[0,0]
            return
        else:
            self.coeff = None


        ## Compute the local conditions at the real numbers (i.e. "p = infinity")
        ## ----------------------------------------------------------------------
        M = Q.matrix()
        E = M.eigenspaces_left()
        M_eigenvalues = [E[i][0]  for i in range(len(E))]

        pos_flag = infinity
        neg_flag = infinity

        for e in M_eigenvalues:
            if e > 0:
                pos_flag = 0
            elif e < 0:
                neg_flag = 0

        real_vec = [infinity, pos_flag, neg_flag, None, None, None, None, None, None]
        self.local_repn_array.append(real_vec)


        ## Compute the local conditions for representability:
        ## --------------------------------------------------
        N = Q.level()
        level_primes = prime_divisors(N)
        prime_repn_modulus_list = [p**(valuation(4*N, p) + 2)  for p in level_primes]

        ## Make a table of local normal forms for each p | N
        local_normal_forms = [Q.local_normal_form(p)  for p in level_primes]

        ## Check local representability conditions for each prime
        for i in range(len(level_primes)):
            p = level_primes[i]
            tmp_local_repn_vec = [p, None, None, None, None, None, None, None, None]
            sqclass = self.squareclass_vector(p)

            ## Check the representability in each Z_p squareclass
            for j in range(len(sqclass)):
                m = sqclass[j]
                k = 0
                repn_flag = False

                while ((repn_flag == False) and (m < 4 * N * p * p)):
                    if (local_normal_forms[i].local_density(p, m) > 0):
                        tmp_local_repn_vec[j+1] = k
                        repn_flag = True
                    k = k + 1
                    m = m * p * p

                ## If we're not represented, write "infinity" to signify
                ## that this squareclass is fully obstructed
                if (repn_flag == False):
                    tmp_local_repn_vec[j+1] = infinity

            ## Test if the conditions at p give exactly Z_p when dim >=3, or
            ## if we represent the elements of even valuation >= 2 when dim = 2.
            omit_flag = True
            if self.dim >= 2:
                ## Check that all entries are zero or 'None'
                for x in tmp_local_repn_vec[1:]:
                    if not ((x == 0) or (x == None)):
                        omit_flag = False

            ## Add the results for this prime if there is a congruence obstruction
            if omit_flag == False:
                self.local_repn_array.append(tmp_local_repn_vec)
                self.exceptional_primes.append(p)


    def __repr__(self):
        """
        Print the local conditions.

        INPUT:
            none

        OUTPUT:
            string

        TO DO:  Improve the output for the real numbers, and special output for locally universality.
        Also give names to the squareclasses, so it's clear what the output means! =)

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.__repr__()
            'This 2-dimensional form represents the p-adic integers of even\nvaluation for all primes p except [2].\nFor these and the reals, we have:\n     Reals:   [0, +Infinity]\n     p = 2:   [0, +Infinity, 0, +Infinity, 0, +Infinity, 0, +Infinity]\n'

        """
        if self.dim == 0:
            out_str = "This 0-dimensional form only represents zero."
        elif self.dim == 1:
            out_str = "This 1-dimensional form only represents square multiples of " + str(self.coeff) + "."
        elif self.dim == 2:
            out_str = "This 2-dimensional form represents the p-adic integers of even\n"
            out_str += "valuation for all primes p except " + str(self.exceptional_primes[1:]) + ".\n"
            out_str += "For these and the reals, we have:\n"
        else:
            out_str = "This form represents the p-adic integers Z_p for all primes p except \n"
            out_str += str(self.exceptional_primes[1:]) + ".  For these and the reals, we have:\n"

        for v in self.local_repn_array:
            if v[0] == infinity:
                out_str += "     " + "Reals:   " + str(v[1:3]) + "\n"
            elif v[0] == 2:
                out_str += "     " + "p = 2:   " + str(v[1:]) + "\n"
            else:
                out_str += "     " + "p = " + str(v[0]) + ":   " + str(v[1:5]) + "\n"

        return out_str



    def __eq__(self, right):
        """
        Determines if two sets of local conditions are equal.

        INPUT:
            right -- a QuadraticFormLocalRepresentationConditions object

        OUTPUT:
            boolean

        EXAMPLES::

             sage: Q1 = DiagonalQuadraticForm(ZZ, [1,1])
             sage: Q2 = DiagonalQuadraticForm(ZZ, [1,1,1])
             sage: Q3 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
             sage: Q4 = DiagonalQuadraticForm(ZZ, [1,1,1,1])

             sage: Q1.local_representation_conditions() == Q2.local_representation_conditions()
             False
             sage: Q1.local_representation_conditions() == Q3.local_representation_conditions()
             False
             sage: Q1.local_representation_conditions() == Q4.local_representation_conditions()
             False
             sage: Q2.local_representation_conditions() == Q3.local_representation_conditions()
             False
             sage: Q3.local_representation_conditions() == Q4.local_representation_conditions()
             True
        """
        if not isinstance(right, QuadraticFormLocalRepresentationConditions):
            return False

        ## Check the dimensions agree when they affect the kind of representation conditions.
        if ((self.dim <= 2) or (right.dim <= 2)) and self.dim != right.dim:
            return False

        ## Check equality by dimension
        if self.dim == 0:
            return True
        elif self.dim == 1:
            return self.coeff == right.coeff     ## Compare coefficients in dimension 1 (since ZZ has only one unit square)
        else:
            return (self.exceptional_primes == right.exceptional_primes) \
                and (self.local_repn_array == right.local_repn_array)


    def squareclass_vector(self, p):
        """
        Gives a vector of integers which are normalized
        representatives for the `p`-adic rational squareclasses
        (or the real squareclasses) at the prime `p`.

        INPUT:
            `p` -- a positive prime number or "infinity".

        OUTPUT:
            a list of integers

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.squareclass_vector(5)
            [1, 2, 5, 10]

        """
        if p == infinity:
            return [1, -1]
        elif p == 2:
            return [1, 3, 5, 7, 2, 6, 10, 14]
        else:
            r = least_quadratic_nonresidue(p)
            return [1, r, p, p*r]



    def local_conditions_vector_for_prime(self, p):
        """
        Returns a local representation vector for the (possibly infinite) prime `p`.

        INPUT:
            `p` -- a positive prime number.  (Is 'infinity' allowed here?)

        OUTPUT:
            a list of integers

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.local_conditions_vector_for_prime(2)
            [2, 0, 0, 0, +Infinity, 0, 0, 0, 0]
            sage: C.local_conditions_vector_for_prime(3)
            [3, 0, 0, 0, 0, None, None, None, None]

        """
        ## Check if p is non-generic
        if p in self.exceptional_primes:
            return deepcopy(self.local_repn_array[self.exceptional_primes.index(p)])

        ## Otherwise, generate a vector at this (finite) prime
        if self.dim >= 3:
            if p == 2:
                return [2, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                return [p, 0, 0, 0, 0, None, None, None, None]

        elif self.dim == 2:
            if p == 2:
                return [2, 0, 0, 0, 0, infinity, infinity, infinity, infinity]
            else:
                return [p, 0, 0, infinity, infinity, None, None, None, None]

        elif self.dim == 1:
            v = [p, None, None, None, None, None, None, None, None]
            sqclass = self.squareclass_vector(p)

            for i in range(len(sq_class)):
                if QQ(self.coeff / sqclass[i]).is_padic_square(p):    ## Note:This should happen only once!
                    nu = valuation(self.coeff / sqclass[i], p) / 2
                else:
                    v[i+1] = infinity

        elif self.dim == 0:
            if p == 2:
                return [2, infinity, infinity, infinity, infinity, infinity, infinity, infinity, infinity]
            else:
                return [p, infinity, infinity, infinity, infinity, None, None, None, None]

        raise RuntimeError, "Error... The dimension stored should be a non-negative integer!"



    def is_universal_at_prime(self, p):
        """
        Determines if the (integer-valued/rational) quadratic form represents all of `Z_p`.

        INPUT:
            `p` -- a positive prime number or "infinity".

        OUTPUT:
            boolean

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_universal_at_prime(2)
            False
            sage: C.is_universal_at_prime(3)
            True
            sage: C.is_universal_at_prime(infinity)
            False

        """
        ## Check if the prime behaves generically for n >= 3.
        if (self.dim >= 3) and not (p in self.exceptional_primes):
            return True

        ## Check if the prime behaves generically for n <= 2.
        if (self.dim <= 2) and not (p in self.exceptional_primes):
            return False

        ## Check if the prime is "infinity" (for the reals)
        if p == infinity:
            v = self.local_repn_array[0]
            if p != v[0]:
                raise RuntimeError, "Error... The first vector should be for the real numbers!"
            return (v[1:3] == [0,0])     ## True iff the form is indefinite

        ## Check non-generic "finite" primes
        v = self.local_conditions_vector_for_prime(p)
        Zp_univ_flag = True
        for nu in v[1:]:
            if (nu != None) and ((nu != 0) or (nu == infinity)):
                Zp_univ_flag = False
        return Zp_univ_flag


    def is_universal_at_all_finite_primes(self):
        """
        Determines if the quadratic form represents `Z_p` for all finite/non-archimedean primes.

        INPUT:
            none

        OUTPUT:
            boolean

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_universal_at_all_finite_primes()
            False

        ::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_universal_at_all_finite_primes()
            True

        """
        ## Check if dim <= 2.
        if self.dim <= 2:
            return False

        ## Check that all non-generic finite primes are universal
        univ_flag = True
        for p in self.exceptional_primes[1:]:      ## Omit p = "infinity" here
            univ_flag = univ_flag and self.is_universal_at_prime(p)
        return univ_flag


    def is_universal_at_all_places(self):
        """
        Determines if the quadratic form represents `Z_p` for all
        finite/non-archimedean primes, and represents all real numbers.

        INPUT:
            none

        OUTPUT:
            boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions

        ::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_universal_at_all_places()
            False

        ::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_universal_at_all_places()
            False

        ::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,-1])
            sage: C = QuadraticFormLocalRepresentationConditions(Q)     # long time (8.5 s)
            sage: C.is_universal_at_all_places()                        # long time
            True

        """
        ## Check if dim <= 2.
        if self.dim <= 2:
            return False

        ## Check that all non-generic finite primes are universal
        for p in self.exceptional_primes:
            if not self.is_universal_at_prime(p):
                return False
        return True



    def is_locally_represented_at_place(self, m, p):
        """
        Determines if the rational number m is locally represented by the
        quadratic form at the (possibly infinite) prime `p`.

        INPUT:

            `m` -- an integer

            `p` -- a positive prime number or "infinity".

        OUTPUT:
            boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_locally_represented_at_place(7, 2)
            False
            sage: C.is_locally_represented_at_place(1, 3)
            True
            sage: C.is_locally_represented_at_place(-1, infinity)
            False
            sage: C.is_locally_represented_at_place(1, infinity)
            True
            sage: C.is_locally_represented_at_place(0, infinity)
            True

        """
        ## Sanity Check
        if not m in QQ:
            raise TypeError, "Oops!  m = " + str(m) +  " is not a rational number!"

       ## Representing zero
        if m == 0:
            return True

        ## 0-dim'l forms
        if self.dim == 0:    ## Here m != 0
            return False

        ## 1-dim'l forms
        if self.dim == 1:
            m1 = QQ(m) / self.coeff
            if p == infinity:
                return (m1 > 0)
            else:
                return (valuation(m1, p) >= 0) and m1.is_padic_square(p)

        ## >= 2-dim'l forms
        local_vec = self.local_conditions_vector_for_prime(p)

        ## Check the real place
        if p == infinity:
            if m > 0:
                return local_vec[1] == 0
            elif m < 0:
                return local_vec[2] == 0
            else:   ## m == 0
                return True

        ## Check at a finite place
        sqclass = self.squareclass_vector(p)
        for s in sqclass:
            #print "m =", m, "   s =", s, "   m/s =", (QQ(m)/s)
            if (QQ(m)/s).is_padic_square(p):
                nu = valuation(m/s, p)
                return local_vec[sqclass.index(s) + 1] <= (nu / 2)



    def is_locally_represented(self, m):
        """
        Determines if the rational number `m` is locally represented by
        the quadratic form (allowing vectors with coefficients in `Z_p` at all
        places).

        INPUT:
            `m` -- an integer

        OUTPUT:
            boolean

        EXAMPLES::

            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.is_locally_represented(7)
            False
            sage: C.is_locally_represented(28)
            False
            sage: C.is_locally_represented(11)
            True
            sage: C.is_locally_represented(QQ(1)/QQ(2))
            False

        """
        ## Representing zero
        if m == 0:
            return True

        ## 0-dim'l forms
        if self.dim == 0:    ## Here m != 0
            return False

        ## 1-dim'l forms
        if self.dim == 1:
            m1 = m / self.coeff
            return (m1 in ZZ) and is_square(m1)

        ## Check the generic primes (when n = 2 or n >= 3)
        m_primes = prime_divisors(numerator(m) * denominator(m))
        for p in m_primes:
            if not p in self.exceptional_primes:
               val = valuation(m, p)
               if (val < 0):
                   return False

        ## Check the non-generic primes (when n = 2 or n >= 3)
        for p in self.exceptional_primes:
            if not self.is_locally_represented_at_place(m, p):
                return False

        ## If we got here, we're locally represented!
        return True



## --------------------  End of QuadraticFormLocalRepresentationConditions Class   ----------------------



def local_representation_conditions(self, recompute_flag=False, silent_flag=False):
    """
    WARNING: THIS ONLY WORKS CORRECTLY FOR FORMS IN >=3 VARIABLES,
        WHICH ARE LOCALLY UNIVERSAL AT ALMOST ALL PRIMES!

    This class finds the local conditions for a number to be integrally
    represented by an integer-valued quadratic form.  These conditions
    are stored in "self.__local_representability_conditions" and
    consist of a list of 9 element vectors, with one for each prime
    with a local obstruction (though only the first 5 are meaningful
    unless `p=2` ).  The first element is always the prime `p` where the
    local obstruction occurs, and the next 8 (or 4) entries represent
    square-classes in the `p`-adic integers `Z_p`, and are labeled by the
    `Q_p` square-classes `t*(Q_p)^2` with `t` given as follows:

        `p > 2`  ==>  [ *  1  u  p  u p  *  *  *  * ]

        `p = 2`  ==>  [ *  1  3  5  7  2  6  10  14 ]

    The integer appearing in each place tells us how `p`-divisible a
    number needs to be in that square-class in order to be locally
    represented by Q.  A negative number indicates that the entire `Q_p`
    square-class is not represented, while a positive number `x` indicates
    that `t*p^{(2*x)} (Z_p)^2` is locally represented but `t*p^{(2*(x-1))}`
    `(Z_p)^2` is not.

    As an example, the vector

        [2  3  0  0  0  0  2  0  infinity]

    tells us that all positive integers are locally represented at p=2
    except those of the forms:

        `2^6 * u * r^2` with  `u = 1 (mod 8)`

        `2^5 * u * r^2` with `u = 3 (mod 8)`

        `2 * u * r^2` with `u = 7 (mod 8)`

    At the real numbers, the vector which looks like

        [infinity, 0, infinity, None, None, None, None, None, None]

    means that Q is negative definite (i.e. the 0 tells us all
    positive reals are represented).  The real vector always appears,
    and is listed before the other ones.

    INPUT:
        none

    OUTPUT:
        A list of 9-element vectors describing the representation
        obstructions at primes dividing the level.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [])
        sage: Q.local_representation_conditions()
        This 0-dimensional form only represents zero.

        sage: Q = DiagonalQuadraticForm(ZZ, [5])
        sage: Q.local_representation_conditions()
        This 1-dimensional form only represents square multiples of 5.

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q1.local_representation_conditions()
        This 2-dimensional form represents the p-adic integers of even
        valuation for all primes p except [2].
        For these and the reals, we have:
         Reals:   [0, +Infinity]
         p = 2:   [0, +Infinity, 0, +Infinity, 0, +Infinity, 0, +Infinity]


        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q1.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [2].  For these and the reals, we have:
         Reals:   [0, +Infinity]
         p = 2:   [0, 0, 0, +Infinity, 0, 0, 0, 0]

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q1.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [].  For these and the reals, we have:
         Reals:   [0, +Infinity]

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,3,3,3])
        sage: Q1.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [3].  For these and the reals, we have:
         Reals:   [0, +Infinity]
         p = 3:   [0, 1, 0, 0]

        sage: Q2 = DiagonalQuadraticForm(ZZ, [2,3,3,3])
        sage: Q2.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [3].  For these and the reals, we have:
         Reals:   [0, +Infinity]
         p = 3:   [1, 0, 0, 0]

        sage: Q3 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q3.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [].  For these and the reals, we have:
         Reals:   [0, +Infinity]

    """
    ## Recompute the local conditions if they don't exist or the recompute_flag is set.
    if (not hasattr(self, "__local_representability_conditions")) or (recompute_flag == True):
        self.__local_representability_conditions = QuadraticFormLocalRepresentationConditions(self)

    ## Return the local conditions if the silent_flag is not set.
    if not silent_flag:
        return self.__local_representability_conditions



def is_locally_universal_at_prime(self, p):
    """
    Determines if the (integer-valued/rational) quadratic form represents all of `Z_p`.

    INPUT:
        `p` -- a positive prime number or "infinity".

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.is_locally_universal_at_prime(2)
        True
        sage: Q.is_locally_universal_at_prime(3)
        True
        sage: Q.is_locally_universal_at_prime(5)
        True
        sage: Q.is_locally_universal_at_prime(infinity)
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.is_locally_universal_at_prime(2)
        False
        sage: Q.is_locally_universal_at_prime(3)
        True
        sage: Q.is_locally_universal_at_prime(5)
        True
        sage: Q.is_locally_universal_at_prime(infinity)
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,-1])
        sage: Q.is_locally_universal_at_prime(infinity)
        True

    """
    self.local_representation_conditions(silent_flag=True)
    return self.__local_representability_conditions.is_universal_at_prime(p)



def is_locally_universal_at_all_primes(self):
    """
    Determines if the quadratic form represents `Z_p` for all finite/non-archimedean primes.

    INPUT:
        none

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.is_locally_universal_at_all_primes()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.is_locally_universal_at_all_primes()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.is_locally_universal_at_all_primes()
        False

    """
    self.local_representation_conditions(silent_flag=True)
    return self.__local_representability_conditions.is_universal_at_all_finite_primes()



def is_locally_universal_at_all_places(self):
    """
    Determines if the quadratic form represents `Z_p` for all
    finite/non-archimedean primes, and represents all real numbers.

    INPUT:
        none

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.is_locally_universal_at_all_places()
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.is_locally_universal_at_all_places()
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,-1])
        sage: Q.is_locally_universal_at_all_places()        # long time (8.5 s)
        True

    """
    self.local_representation_conditions(silent_flag=True)
    return self.__local_representability_conditions.is_universal_at_all_places()



def is_locally_represented_number_at_place(self, m, p):
    """
    Determines if the rational number m is locally represented by the
    quadratic form at the (possibly infinite) prime `p`.

    INPUT:

        `m` -- an integer

        `p` -- a prime number > 0 or 'infinity'

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.is_locally_represented_number_at_place(7, infinity)
        True
        sage: Q.is_locally_represented_number_at_place(7, 2)
        False
        sage: Q.is_locally_represented_number_at_place(7, 3)
        True
        sage: Q.is_locally_represented_number_at_place(7, 5)
        True
        sage: Q.is_locally_represented_number_at_place(-1, infinity)
        False
        sage: Q.is_locally_represented_number_at_place(-1, 2)
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,-1])
        sage: Q.is_locally_represented_number_at_place(7, infinity)     # long time (8.5 s)
        True
        sage: Q.is_locally_represented_number_at_place(7, 2)            # long time
        True
        sage: Q.is_locally_represented_number_at_place(7, 3)            # long time
        True
        sage: Q.is_locally_represented_number_at_place(7, 5)            # long time
        True

    """
    self.local_representation_conditions(silent_flag=True)
    return self.__local_representability_conditions.is_locally_represented_at_place(m, p)



def is_locally_represented_number(self, m):
    """
    Determines if the rational number m is locally represented by the quadratic form.

    INPUT:
        `m` -- an integer

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.is_locally_represented_number(2)
        True
        sage: Q.is_locally_represented_number(7)
        False
        sage: Q.is_locally_represented_number(-1)
        False
        sage: Q.is_locally_represented_number(28)
        False
        sage: Q.is_locally_represented_number(0)
        True

    """
    self.local_representation_conditions(silent_flag=True)
    return self.__local_representability_conditions.is_locally_represented(m)
