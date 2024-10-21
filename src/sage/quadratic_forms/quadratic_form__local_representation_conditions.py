# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.pari sage.rings.number_field
"""
Local Representation Conditions
"""
########################################################################
# Class for keeping track of the local conditions for representability #
# of numbers by a quadratic form over ZZ (and eventually QQ also).     #
########################################################################
from __future__ import annotations
from copy import deepcopy

from sage.arith.misc import is_square, prime_divisors, valuation
from sage.misc.functional import numerator, denominator
from sage.quadratic_forms.extras import least_quadratic_nonresidue
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class QuadraticFormLocalRepresentationConditions:
    """
    A class for dealing with the local conditions of a
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
        sage: L == list(range(100))
        True

    ::

        sage: Q3 = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q3.local_representation_conditions()
        This form represents the p-adic integers Z_p for all primes p except
        [2].  For these and the reals, we have:
             Reals:   [0, +Infinity]
             p = 2:   [0, 0, 0, +Infinity, 0, 0, 0, 0]
        sage: E = [m  for m in range(100)  if not Q3.is_locally_represented_number(m)]
        sage: E1 = [m  for m in range(1,100)  if m / 2**(2 * (valuation(m,2) // 2)) % 8 == 7]
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
        sage: L1 = [0] + [m for m in range(1, 25)
        ....:             if len([p for p in prime_factors(squarefree_part(ZZ(m)))
        ....:                       if (p % 4) == 3]) % 2 == 0]
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
        r"""
        Take a :class:`QuadraticForm` and computes its local conditions (if
        they do not already exist).  The ``recompute_flag`` overrides the
        previously computed conditions if they exist, and stores the
        new conditions.

        INPUT:

        - ``Q`` -- Quadratic form over `\ZZ`

        OUTPUT: a :class:`QuadraticFormLocalRepresentationConditions` object

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: QuadraticFormLocalRepresentationConditions(Q)
            This form represents the p-adic integers Z_p for all primes p except
            [].  For these and the reals, we have:
                 Reals:   [0, +Infinity]
        """
        # Check that the form Q is integer-valued (we can relax this later)
        if Q.base_ring() != ZZ:
            raise TypeError("We require that the quadratic form be defined over ZZ (integer-values) for now.")

        # Basic structure initialization
        self.local_repn_array = []    # List of all local conditions
        self.dim = Q.dim()       # We allow this to be any nonnegative integer.
        self.exceptional_primes = [infinity]

        # Deal with the special cases of 0 and 1-dimensional forms
        if self.dim == 0:
            self.coeff = None
            return
        if self.dim == 1:
            self.coeff = Q[0, 0]
            return

        self.coeff = None

        # Compute the local conditions at the real numbers (i.e. "p = infinity")
        # ----------------------------------------------------------------------
        M = Q.matrix()
        E = M.eigenspaces_left()
        M_eigenvalues = [E[i][0] for i in range(len(E))]

        pos_flag = infinity
        neg_flag = infinity

        for e in M_eigenvalues:
            if e > 0:
                pos_flag = 0
            elif e < 0:
                neg_flag = 0

        real_vec = [infinity, pos_flag, neg_flag, None, None, None, None, None, None]
        self.local_repn_array.append(real_vec)

        # Compute the local conditions for representability:
        # --------------------------------------------------
        N = Q.level()
        level_primes = prime_divisors(N)

        # Make a table of local normal forms for each p | N
        local_normal_forms = [Q.local_normal_form(p) for p in level_primes]

        # Check local representability conditions for each prime
        for i, p in enumerate(level_primes):
            tmp_local_repn_vec = [p, None, None, None, None, None, None, None, None]
            sqclass = self.squareclass_vector(p)

            # Check the representability in each Z_p squareclass
            for j, m in enumerate(sqclass):
                k = 0
                repn_flag = False

                while ((not repn_flag) and (m < 4 * N * p * p)):
                    if local_normal_forms[i].local_density(p, m) > 0:
                        tmp_local_repn_vec[j + 1] = k
                        repn_flag = True
                    k = k + 1
                    m = m * p * p

                # If we're not represented, write "infinity" to signify
                # that this squareclass is fully obstructed
                if not repn_flag:
                    tmp_local_repn_vec[j + 1] = infinity

            # Test if the conditions at p give exactly Z_p when dim >=3, or
            # if we represent the elements of even valuation >= 2 when dim = 2.
            omit_flag = True
            if self.dim >= 2:
                # Check that all entries are zero or 'None'
                for x in tmp_local_repn_vec[1:]:
                    if not ((x == 0) or (x is None)):
                        omit_flag = False

            # Add the results for this prime if there is a congruence obstruction
            if not omit_flag:
                self.local_repn_array.append(tmp_local_repn_vec)
                self.exceptional_primes.append(p)

    def __repr__(self) -> str:
        r"""
        Print the local conditions.

        OUTPUT: string

        .. TODO::

            Improve the output for the real numbers, and special output for locally universality.
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

    def __eq__(self, right) -> bool:
        """
        Determine if two sets of local conditions are equal.

        INPUT:

        - ``right`` -- a QuadraticFormLocalRepresentationConditions object

        OUTPUT: boolean

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

        # Check the dimensions agree when they affect the kind of representation conditions.
        if ((self.dim <= 2) or (right.dim <= 2)) and self.dim != right.dim:
            return False

        # Check equality by dimension
        if self.dim == 0:
            return True
        if self.dim == 1:
            return self.coeff == right.coeff     # Compare coefficients in dimension 1 (since ZZ has only one unit square)
        return ((self.exceptional_primes == right.exceptional_primes)
                and (self.local_repn_array == right.local_repn_array))

    def squareclass_vector(self, p) -> list:
        """
        Return a list of integers which are normalized
        representatives for the `p`-adic rational squareclasses
        (or the real squareclasses) at the prime `p`.

        INPUT:

        - ``p`` -- a positive prime number or "infinity"

        OUTPUT: list of integers

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.squareclass_vector(5)
            [1, 2, 5, 10]
        """
        if p == infinity:
            return [1, -1]
        if p == 2:
            return [1, 3, 5, 7, 2, 6, 10, 14]
        r = least_quadratic_nonresidue(p)
        return [1, r, p, p * r]

    def local_conditions_vector_for_prime(self, p) -> list:
        """
        Return a local representation vector for the (possibly infinite) prime `p`.

        INPUT:

        - ``p`` -- a positive prime number.  (Is 'infinity' allowed here?)

        OUTPUT: list of integers

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: from sage.quadratic_forms.quadratic_form__local_representation_conditions import QuadraticFormLocalRepresentationConditions
            sage: C = QuadraticFormLocalRepresentationConditions(Q)
            sage: C.local_conditions_vector_for_prime(2)
            [2, 0, 0, 0, +Infinity, 0, 0, 0, 0]
            sage: C.local_conditions_vector_for_prime(3)
            [3, 0, 0, 0, 0, None, None, None, None]
        """
        # Check if p is non-generic
        if p in self.exceptional_primes:
            return deepcopy(self.local_repn_array[self.exceptional_primes.index(p)])

        # Otherwise, generate a vector at this (finite) prime
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

            for i, sqi in enumerate(sqclass):
                if QQ(self.coeff / sqi).is_padic_square(p):    # Note:This should happen only once!
                    nu = valuation(self.coeff / sqi, p) / 2    # UNUSED VARIABLE !
                else:
                    v[i + 1] = infinity

        elif self.dim == 0:
            if p == 2:
                return [2, infinity, infinity, infinity, infinity, infinity, infinity, infinity, infinity]
            return [p, infinity, infinity, infinity, infinity, None, None, None, None]

        raise RuntimeError("the stored dimension should be a nonnegative integer")

    def is_universal_at_prime(self, p) -> bool:
        r"""
        Determine if the (integer-valued/rational) quadratic form represents all of `\ZZ_p`.

        INPUT:

        - ``p`` -- a positive prime number or "infinity"

        OUTPUT: boolean

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
        # Check if the prime behaves generically for n >= 3.
        if self.dim >= 3 and p not in self.exceptional_primes:
            return True

        # Check if the prime behaves generically for n <= 2.
        if self.dim <= 2 and p not in self.exceptional_primes:
            return False

        # Check if the prime is "infinity" (for the reals)
        if p == infinity:
            v = self.local_repn_array[0]
            if p != v[0]:
                raise RuntimeError("Error... The first vector should be for the real numbers!")
            return (v[1:3] == [0, 0])     # True iff the form is indefinite

        # Check non-generic "finite" primes
        v = self.local_conditions_vector_for_prime(p)
        Zp_univ_flag = True
        for nu in v[1:]:
            if (nu is not None) and ((nu != 0) or (nu == infinity)):
                Zp_univ_flag = False
        return Zp_univ_flag

    def is_universal_at_all_finite_primes(self) -> bool:
        r"""
        Determine if the quadratic form represents `\ZZ_p` for all finite/non-archimedean primes.

        OUTPUT: boolean

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
        # Check if dim <= 2.
        if self.dim <= 2:
            return False

        # Check that all non-generic finite primes are universal
        # Omit p = "infinity" here
        return all(self.is_universal_at_prime(p)
                   for p in self.exceptional_primes[1:])

    def is_universal_at_all_places(self) -> bool:
        r"""
        Determine if the quadratic form represents `\ZZ_p` for all
        finite/non-archimedean primes, and represents all real numbers.

        OUTPUT: boolean

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
        # Check if dim <= 2.
        if self.dim <= 2:
            return False

        # Check that all non-generic finite primes are universal
        return all(self.is_universal_at_prime(p)
                   for p in self.exceptional_primes)

    def is_locally_represented_at_place(self, m, p) -> bool:
        """
        Determine if the rational number `m` is locally represented by the
        quadratic form at the (possibly infinite) prime `p`.

        INPUT:

        - ``m`` -- integer

        - ``p`` -- a positive prime number or "infinity"

        OUTPUT: boolean

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
        # Sanity Check
        if m not in QQ:
            raise TypeError(f"m = {m} is not a rational number")

        # Representing zero
        if m == 0:
            return True

        # 0-dim'l forms
        if self.dim == 0:   # Here m != 0
            return False

        # 1-dim'l forms
        if self.dim == 1:
            m1 = QQ(m) / self.coeff
            if p == infinity:
                return m1 > 0
            return (valuation(m1, p) >= 0) and m1.is_padic_square(p)

        # >= 2-dim'l forms
        local_vec = self.local_conditions_vector_for_prime(p)

        # Check the real place
        if p == infinity:
            if m > 0:
                return local_vec[1] == 0
            if m < 0:
                return local_vec[2] == 0
            # m == 0
            return True

        # Check at a finite place
        sqclass = self.squareclass_vector(p)
        for s in sqclass:
            if (QQ(m) / s).is_padic_square(p):
                nu = valuation(m // s, p)
                return local_vec[sqclass.index(s) + 1] <= (nu / 2)

    def is_locally_represented(self, m) -> bool:
        r"""
        Determine if the rational number `m` is locally represented by
        the quadratic form (allowing vectors with coefficients in `\ZZ_p` at all
        places).

        INPUT:

        - ``m`` -- integer

        OUTPUT: boolean

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
        # Representing zero
        if m == 0:
            return True

        # 0-dim'l forms
        if self.dim == 0:    # Here m != 0
            return False

        # 1-dim'l forms
        if self.dim == 1:
            m1 = m / self.coeff
            return (m1 in ZZ) and is_square(m1)

        # Check the generic primes (when n = 2 or n >= 3)
        m_primes = prime_divisors(numerator(m) * denominator(m))
        for p in m_primes:
            if p not in self.exceptional_primes:
                val = valuation(m, p)
                if val < 0:
                    return False

        # Check the non-generic primes (when n = 2 or n >= 3)
        for p in self.exceptional_primes:
            if not self.is_locally_represented_at_place(m, p):
                return False

        # If we got here, we're locally represented!
        return True

# ---  End of QuadraticFormLocalRepresentationConditions Class ---


def local_representation_conditions(self, recompute_flag=False, silent_flag=False):
    r"""
    .. WARNING::

        This only works correctly for forms in >=3 variables,
        which are locally universal at almost all primes!

    This class finds the local conditions for a number to be integrally
    represented by an integer-valued quadratic form.  These conditions
    are stored in ``self.__local_representability_conditions`` and
    consist of a list of 9 element vectors, with one for each prime
    with a local obstruction (though only the first 5 are meaningful
    unless `p=2`).  The first element is always the prime `p` where the
    local obstruction occurs, and the next 8 (or 4) entries represent
    square-classes in the `p`-adic integers `\ZZ_p`, and are labeled by the
    `\QQ_p` square-classes `t\cdot (\QQ_p)^2` with `t` given as follows:

    - for `p > 2`, ``[ *  1  u  p  u p  *  *  *  * ]``,

    - for `p = 2`, ``[ *  1  3  5  7  2  6  10  14 ]``.

    The integer appearing in each place tells us how `p`-divisible a
    number needs to be in that square-class in order to be locally
    represented by `Q`.  A negative number indicates that the entire `\QQ_p`
    square-class is not represented, while a positive number `x` indicates
    that `t\cdot p^{(2\cdot x)} (\ZZ_p)^2` is locally represented but `t\cdot p^{(2\cdot (x-1))}`
    `(\ZZ_p)^2` is not.

    As an example, the vector ``[2  3  0  0  0  0  2  0  infinity]``
    tells us that all positive integers are locally represented at `p=2`
    except those of the forms:

    - `2^6\cdot u\cdot r^2` with `u = 1` (mod 8)

    - `2^5\cdot u\cdot r^2` with `u = 3` (mod 8)

    - `2\cdot u\cdot r^2` with `u = 7` (mod 8)

    At the real numbers, the vector which looks like ``[infinity, 0, infinity, None, None, None, None, None, None]``
    means that `Q` is negative definite (i.e., the 0 tells us all
    positive reals are represented).  The real vector always appears,
    and is listed before the other ones.

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
    # Recompute the local conditions if they do not exist or the recompute_flag is set.
    if not hasattr(self, "__local_representability_conditions") or recompute_flag:
        self.__local_representability_conditions = QuadraticFormLocalRepresentationConditions(self)

    # Return the local conditions if the silent_flag is not set.
    if not silent_flag:
        return self.__local_representability_conditions


def is_locally_universal_at_prime(self, p) -> bool:
    r"""
    Determine if the (integer-valued/rational) quadratic form represents all of `\ZZ_p`.

    INPUT:

    - ``p`` -- a positive prime number or "infinity"

    OUTPUT: boolean

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


def is_locally_universal_at_all_primes(self) -> bool:
    r"""
    Determine if the quadratic form represents `\ZZ_p` for all finite/non-archimedean primes.

    OUTPUT: boolean

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


def is_locally_universal_at_all_places(self) -> bool:
    r"""
    Determine if the quadratic form represents `\ZZ_p` for all
    finite/non-archimedean primes, and represents all real numbers.

    OUTPUT: boolean

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


def is_locally_represented_number_at_place(self, m, p) -> bool:
    """
    Determine if the rational number `m` is locally represented by the
    quadratic form at the (possibly infinite) prime `p`.

    INPUT:

    - ``m`` -- integer

    - ``p`` -- a prime number > 0 or 'infinity'

    OUTPUT: boolean

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


def is_locally_represented_number(self, m) -> bool:
    """
    Determine if the rational number `m` is locally represented
    by the quadratic form.

    INPUT:

    - ``m`` -- integer

    OUTPUT: boolean

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
