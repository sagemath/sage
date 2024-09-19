"""
Local Density Congruence
"""
##########################################################################
#  Methods which compute the local densities for representing a number
#  by a quadratic form at a prime (possibly subject to additional
#  congruence conditions).
##########################################################################
from copy import deepcopy

from sage.sets.set import Set
from sage.rings.rational_field import QQ
from sage.arith.misc import valuation
from sage.misc.verbose import verbose

from sage.quadratic_forms.count_local_2 import count_modp__by_gauss_sum, CountAllLocalGoodTypesNormalForm


def count_modp_solutions__by_Gauss_sum(self, p, m):
    r"""
    Return the number of solutions of `Q(x) = m` (mod `p`) of a
    non-degenerate quadratic form over the finite field `\ZZ/p\ZZ`,
    where `p` is a prime number > 2.

    .. NOTE::

        We adopt the useful convention that a zero-dimensional
        quadratic form has exactly one solution always (i.e. the empty
        vector).

    These are defined in Table 1 on p363 of Hanke's "Local Densities..." paper.

    INPUT:

    - ``p`` -- a prime number > 2

    - ``m`` -- integer

    OUTPUT: integer `\geq 0`

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: [Q.count_modp_solutions__by_Gauss_sum(3, m)  for m in range(3)]
        [9, 6, 12]

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,2])
        sage: [Q.count_modp_solutions__by_Gauss_sum(3, m)  for m in range(3)]
        [9, 12, 6]
    """
    if self.dim() == 0:
        return 1
    return count_modp__by_gauss_sum(self.dim(), p, m, self.Gram_det())


def local_good_density_congruence_odd(self, p, m, Zvec, NZvec):
    """
    Find the Good-type local density of `Q` representing `m` at `p`.
    (Assuming that `p > 2` and `Q` is given in local diagonal form.)

    The additional congruence condition arguments ``Zvec`` and ``NZvec`` can
    be either a list of indices or None.  ``Zvec=[]`` is equivalent to
    ``Zvec=None``, which both impose no additional conditions, but
    ``NZvec=[]`` returns no solutions always while ``NZvec=None`` imposes no
    additional condition.

    .. TODO::

        Add type checking for ``Zvec``, ``NZvec``, and that `Q` is in local
        normal form.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_good_density_congruence_odd(3, 1, None, None)
        2/3

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_good_density_congruence_odd(3, 1, None, None)
        8/9
    """
    n = self.dim()

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    #  Assuming Q is diagonal, find the indices of the p-unit (diagonal) entries
    UnitVec = Set(i for i in range(n) if self[i, i] % p)
    NonUnitVec = Set(range(n)) - UnitVec

    #  Take cases on the existence of additional nonzero congruence conditions (mod p)
    UnitVec_minus_Zvec = list(UnitVec - Set(Zvec))
    NonUnitVec_minus_Zvec = list(NonUnitVec - Set(Zvec))
    Q_Unit_minus_Zvec = self.extract_variables(UnitVec_minus_Zvec)

    if NZvec is None:
        if m % p:
            total = Q_Unit_minus_Zvec.count_modp_solutions__by_Gauss_sum(p, m) * p**len(NonUnitVec_minus_Zvec)
        else:
            total = (Q_Unit_minus_Zvec.count_modp_solutions__by_Gauss_sum(p, m) - 1) * p**len(NonUnitVec_minus_Zvec)

    else:
        UnitVec_minus_ZNZvec = list(UnitVec - (Set(Zvec) + Set(NZvec)))
        NonUnitVec_minus_ZNZvec = list(NonUnitVec - (Set(Zvec) + Set(NZvec)))
        Q_Unit_minus_ZNZvec = self.extract_variables(UnitVec_minus_ZNZvec)

        if m % p:
            total = Q_Unit_minus_Zvec.count_modp_solutions__by_Gauss_sum(p, m) * p**len(NonUnitVec_minus_Zvec) \
                - Q_Unit_minus_ZNZvec.count_modp_solutions__by_Gauss_sum(p, m) * p**len(NonUnitVec_minus_ZNZvec)
        else:
            total = (Q_Unit_minus_Zvec.count_modp_solutions__by_Gauss_sum(p, m) - 1) * p**len(NonUnitVec_minus_Zvec) \
                - (Q_Unit_minus_ZNZvec.count_modp_solutions__by_Gauss_sum(p, m) - 1) * p**len(NonUnitVec_minus_ZNZvec)

    #  Return the Good-type representation density
    good_density = QQ(total) / p**(n - 1)
    return good_density


def local_good_density_congruence_even(self, m, Zvec, NZvec):
    r"""
    Find the Good-type local density of `Q` representing `m` at `p=2`.
    (Assuming `Q` is given in local diagonal form.)

    The additional congruence condition arguments ``Zvec`` and ``NZvec`` can
    be either a list of indices or None.  ``Zvec=[]`` is equivalent to
    ``Zvec=None`` which both impose no additional conditions, but
    ``NZvec=[]`` returns no solutions always while ``NZvec=None`` imposes no
    additional condition.

    .. WARNING::

        Here the indices passed in ``Zvec`` and ``NZvec`` represent
        indices of the solution vector `x` of `Q(x) = m` (mod `p^k`), and *not*
        the Jordan components of `Q`.  They therefore are required (and
        assumed) to include either all or none of the indices of a given
        Jordan component of `Q`.  This is only important when `p=2` since
        otherwise all Jordan blocks are `1 \times 1`, and so there the indices and
        Jordan blocks coincide.

    .. TODO::

        Add type checking for ``Zvec`` and ``NZvec``, and that `Q` is in local
        normal form.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and 2-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_good_density_congruence_even(1, None, None)
        1

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_good_density_congruence_even(1, None, None)
        1
        sage: Q.local_good_density_congruence_even(2, None, None)
        3/2
        sage: Q.local_good_density_congruence_even(3, None, None)
        1
        sage: Q.local_good_density_congruence_even(4, None, None)
        1/2

    ::

        sage: Q = QuadraticForm(ZZ, 4, range(10))
        sage: Q[0,0] = 5
        sage: Q[1,1] = 10
        sage: Q[2,2] = 15
        sage: Q[3,3] = 20
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 5 1 2 3 ]
        [ * 10 5 6 ]
        [ * * 15 8 ]
        [ * * * 20 ]
        sage: Q.theta_series(20)                                                        # needs sage.libs.pari
        1 + 2*q^5 + 2*q^10 + 2*q^14 + 2*q^15 + 2*q^16 + 2*q^18 + O(q^20)
        sage: Q.local_normal_form(2)                                                    # needs sage.libs.pari sage.rings.padics
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 0 0 ]
        [ * 0 0 0 ]
        [ * * 0 1 ]
        [ * * * 0 ]
        sage: Q.local_good_density_congruence_even(1, None, None)
        3/4
        sage: Q.local_good_density_congruence_even(2, None, None)
        1
        sage: Q.local_good_density_congruence_even(5, None, None)
        3/4
    """
    n = self.dim()

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    #  Find the indices of x for which the associated Jordan blocks are nonzero mod 8    TODO: Move this to special Jordan block code separately!
    #  -------------------------------------------------------------------------------
    Not8vec = []
    for i in range(n):

        #  DIAGNOSTIC
        verbose(" i = " + str(i))
        verbose(" n = " + str(n))
        verbose(" Not8vec = " + str(Not8vec))

        nz_flag = False

        #  Check if the diagonal entry isn't divisible 8
        if self[i, i] % 8:
            nz_flag = True

        #  Check appropriate off-diagonal entries aren't divisible by 8
        else:

            #  Special check for first off-diagonal entry
            if i == 0 and self[i, i + 1] % 8:
                nz_flag = True

            #  Special check for last off-diagonal entry
            elif i == n - 1 and self[i - 1, i] % 8:
                nz_flag = True

            #  Check for the middle off-diagonal entries
            else:
                if (i > 0) and (i < n - 1) and (self[i, i + 1] % 8 or
                                                self[i - 1, i] % 8):
                    nz_flag = True

        #  Remember the (vector) index if it's not part of a Jordan block of norm divisible by 8
        if nz_flag:
            Not8vec += [i]

    #  Compute the number of Good-type solutions mod 8:
    #  ------------------------------------------------

    #  Setup the indexing sets for additional zero congruence solutions
    Q_Not8 = self.extract_variables(Not8vec)
    Not8 = Set(Not8vec)
    Is8 = Set(range(n)) - Not8

    Z = Set(Zvec)
    Z_Not8 = Not8.intersection(Z)
    Z_Is8 = Is8.intersection(Z)
    Is8_minus_Z = Is8 - Z_Is8

    # DIAGNOSTIC
    verbose("Z = " + str(Z))
    verbose("Z_Not8 = " + str(Z_Not8))
    verbose("Z_Is8 = " + str(Z_Is8))
    verbose("Is8_minus_Z = " + str(Is8_minus_Z))

    # Take cases on the existence of additional nonzero congruence conditions (mod 2)
    if NZvec is None:
        total = (4 ** len(Z_Is8)) * (8 ** len(Is8_minus_Z)) \
            * CountAllLocalGoodTypesNormalForm(Q_Not8,2, 3, m, list(Z_Not8), None)
    else:
        ZNZ = Z + Set(NZvec)
        ZNZ_Not8 = Not8.intersection(ZNZ)
        ZNZ_Is8 = Is8.intersection(ZNZ)
        Is8_minus_ZNZ = Is8 - ZNZ_Is8

        #  DIAGNOSTIC
        verbose("ZNZ = " + str(ZNZ))
        verbose("ZNZ_Not8 = " + str(ZNZ_Not8))
        verbose("ZNZ_Is8 = " + str(ZNZ_Is8))
        verbose("Is8_minus_ZNZ = " + str(Is8_minus_ZNZ))

        total = (4 ** len(Z_Is8)) * (8 ** len(Is8_minus_Z)) \
            * CountAllLocalGoodTypesNormalForm(Q_Not8, 2, 3, m, list(Z_Not8), None) \
            - (4 ** len(ZNZ_Is8)) * (8 ** len(Is8_minus_ZNZ)) \
            * CountAllLocalGoodTypesNormalForm(Q_Not8, 2, 3, m, list(ZNZ_Not8), None)

    # DIAGNOSTIC
    verbose("total = " + str(total))

    # Return the associated Good-type representation density
    return QQ(total) / 8**(n - 1)


def local_good_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the Good-type local density of `Q` representing `m` at `p`.
    (Front end routine for parity specific routines for `p`.)

    .. TODO::

        Add documentation about the additional congruence
        conditions ``Zvec`` and ``NZvec``.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_good_density_congruence(2, 1, None, None)
        1
        sage: Q.local_good_density_congruence(3, 1, None, None)
        2/3

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_good_density_congruence(2, 1, None, None)
        1
        sage: Q.local_good_density_congruence(3, 1, None, None)
        8/9
    """
    #  DIAGNOSTIC
    verbose(" In local_good_density_congruence with ")
    verbose(" Q is: \n" + str(self))
    verbose(" p = " + str(p))
    verbose(" m = " + str(m))
    verbose(" Zvec = " + str(Zvec))
    verbose(" NZvec = " + str(NZvec))

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    n = self.dim()

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    # There was here a commented-out check that Q is in local normal form
    # (it often may not be since the reduction procedure
    # often mixes up the order of the valuations...)
    # This commented-out code was removed in issue #32960

    # Decide which routine to use to compute the Good-type density
    if p > 2:
        return self.local_good_density_congruence_odd(p, m, Zvec, NZvec)

    if p == 2:
        return self.local_good_density_congruence_even(m, Zvec, NZvec)

    raise RuntimeError("\n Error in Local_Good_Density: The 'prime' p = " + str(p) + " is < 2. \n")


def local_zero_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the Zero-type local density of `Q` representing `m` at `p`,
    allowing certain congruence conditions mod `p`.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_zero_density_congruence(2, 2, None, None)
        0
        sage: Q.local_zero_density_congruence(2, 4, None, None)
        1/2
        sage: Q.local_zero_density_congruence(3, 6, None, None)
        0
        sage: Q.local_zero_density_congruence(3, 9, None, None)
        2/9

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_zero_density_congruence(2, 2, None, None)
        0
        sage: Q.local_zero_density_congruence(2, 4, None, None)
        1/4
        sage: Q.local_zero_density_congruence(3, 6, None, None)
        0
        sage: Q.local_zero_density_congruence(3, 9, None, None)
        8/81
    """
    #  DIAGNOSTIC
    verbose(" In local_zero_density_congruence with ")
    verbose(" Q is: \n" + str(self))
    verbose(" p = " + str(p))
    verbose(" m = " + str(m))
    verbose(" Zvec = " + str(Zvec))
    verbose(" NZvec = " + str(NZvec))

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    n = self.dim()

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    p2 = p * p

    #  Check some conditions for no zero-type solutions to exist
    if m % p2 or NZvec is not None:
        return 0

    #  Use the reduction procedure to return the result
    return self.local_density_congruence(p, m / p2, None, None) / p**(self.dim() - 2)


def local_badI_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the Bad-type I local density of `Q` representing `m` at `p`.
    (Assuming that `p > 2` and `Q` is given in local diagonal form.)

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_badI_density_congruence(2, 1, None, None)
        0
        sage: Q.local_badI_density_congruence(2, 2, None, None)
        1
        sage: Q.local_badI_density_congruence(2, 4, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 1, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 6, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 9, None, None)
        0

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_badI_density_congruence(2, 1, None, None)
        0
        sage: Q.local_badI_density_congruence(2, 2, None, None)
        0
        sage: Q.local_badI_density_congruence(2, 4, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 2, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 6, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 9, None, None)
        0

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,3,9])
        sage: Q.local_badI_density_congruence(3, 1, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 3, None, None)
        4/3
        sage: Q.local_badI_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_badI_density_congruence(3, 9, None, None)
        0
        sage: Q.local_badI_density_congruence(3, 18, None, None)
        0
    """
    #  DIAGNOSTIC
    verbose(" In local_badI_density_congruence with ")
    verbose(" Q is: \n" + str(self))
    verbose(" p = " + str(p))
    verbose(" m = " + str(m))
    verbose(" Zvec = " + str(Zvec))
    verbose(" NZvec = " + str(NZvec))

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    n = self.dim()

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    #  Define the indexing set S_0, and determine if S_1 is empty:
    #  -----------------------------------------------------------
    S0 = []
    S1_empty_flag = True
    # This is used to check if we should be computing BI solutions at all!
    # (We should really to this earlier, but S1 must be nonzero to proceed.)

    #  Find the valuation of each variable (which will be the same over 2x2 blocks),
    #  remembering those of valuation 0 and if an entry of valuation 1 exists.
    for i in range(n):

        #  Compute the valuation of each index, allowing for off-diagonal terms
        if self[i, i] == 0:
            if i == 0:
                val = valuation(self[i, i + 1], p)    # Look at the term to the right
            else:
                if i == n - 1:
                    val = valuation(self[i - 1, i], p)    # Look at the term above
                else:
                    val = valuation(self[i, i + 1] + self[i - 1, i], p)    # Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(self[i, i], p)

        if val == 0:
            S0 += [i]
        elif val == 1:
            S1_empty_flag = False    # Need to have a non-empty S1 set to proceed with Bad-type I reduction...

    #  Check that S1 is non-empty and p|m to proceed, otherwise return no solutions.
    if S1_empty_flag or m % p:
        return 0

    #  Check some conditions for no bad-type I solutions to exist
    if (NZvec is not None) and (len(Set(S0).intersection(Set(NZvec))) != 0):
        return 0

    #  Check that the form is primitive...     WHY DO WE NEED TO DO THIS?!?
    if not S0:
        print(" Using Q = " + str(self))
        print(" and p = " + str(p))
        raise RuntimeError("the form is not primitive")

    #  DIAGNOSTIC
    verbose(" m = " + str(m) + "   p = " + str(p))
    verbose(" S0 = " + str(S0))
    verbose(" len(S0) = " + str(len(S0)))

    #  Make the form Qnew for the reduction procedure:
    #  -----------------------------------------------
    Qnew = deepcopy(self)    # TODO:  DO THIS WITHOUT A copy()
    for i in range(n):
        if i in S0:
            Qnew[i, i] = p * Qnew[i, i]
            if ((p == 2) and (i < n - 1)):
                Qnew[i, i + 1] = p * Qnew[i, i + 1]
        else:
            Qnew[i, i] = Qnew[i, i] / p
            if ((p == 2) and (i < n - 1)):
                Qnew[i, i + 1] = Qnew[i, i + 1] / p

    #  DIAGNOSTIC
    verbose("\n\n Check of Bad-type I reduction: \n")
    verbose(" Q is " + str(self))
    verbose(" Qnew is " + str(Qnew))
    verbose(" p = " + str(p))
    verbose(" m / p = " + str(m / p))
    verbose(" NZvec " + str(NZvec))

    #  Do the reduction
    Zvec_geq_1 = list(Set([i for i in Zvec if i not in S0]))
    if NZvec is None:
        NZvec_geq_1 = NZvec
    else:
        NZvec_geq_1 = list(Set([i for i in NZvec if i not in S0]))

    return QQ(p**(1 - len(S0))) * Qnew.local_good_density_congruence(p, m / p, Zvec_geq_1, NZvec_geq_1)


def local_badII_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the Bad-type II local density of `Q` representing `m` at `p`.
    (Assuming that `p > 2` and `Q` is given in local diagonal form.)

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_badII_density_congruence(2, 1, None, None)
        0
        sage: Q.local_badII_density_congruence(2, 2, None, None)
        0
        sage: Q.local_badII_density_congruence(2, 4, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 1, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 6, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 9, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 27, None, None)
        0

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,3,9,9])
        sage: Q.local_badII_density_congruence(3, 1, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 3, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 6, None, None)
        0
        sage: Q.local_badII_density_congruence(3, 9, None, None)
        4/27
        sage: Q.local_badII_density_congruence(3, 18, None, None)
        4/9
    """
    #  DIAGNOSTIC
    verbose(" In local_badII_density_congruence with ")
    verbose(" Q is: \n" + str(self))
    verbose(" p = " + str(p))
    verbose(" m = " + str(m))
    verbose(" Zvec = " + str(Zvec))
    verbose(" NZvec = " + str(NZvec))

    #  Put the Zvec congruence condition in a standard form
    if Zvec is None:
        Zvec = []

    n = self.dim()

    #  Sanity Check on Zvec and NZvec:
    #  -------------------------------
    Sn = Set(range(n))
    if (Zvec is not None) and (len(Set(Zvec) + Sn) > n):
        raise RuntimeError("Zvec must be a subset of {0, ..., n-1}.")
    if (NZvec is not None) and (len(Set(NZvec) + Sn) > n):
        raise RuntimeError("NZvec must be a subset of {0, ..., n-1}.")

    #  Define the indexing sets S_i:
    #  -----------------------------
    S0 = []
    S1 = []
    S2plus = []

    for i in range(n):

        #  Compute the valuation of each index, allowing for off-diagonal terms
        if self[i, i] == 0:
            if i == 0:
                val = valuation(self[i, i + 1], p)    # Look at the term to the right
            elif i == n - 1:
                val = valuation(self[i - 1, i], p)    # Look at the term above
            else:
                val = valuation(self[i, i + 1] + self[i - 1, i], p)    # Finds the valuation of the off-diagonal term since only one isn't zero
        else:
            val = valuation(self[i, i], p)

        #  Sort the indices into disjoint sets by their valuation
        if (val == 0):
            S0 += [i]
        elif (val == 1):
            S1 += [i]
        elif (val >= 2):
            S2plus += [i]

    #  Check that S2 is non-empty and p^2 divides m to proceed, otherwise return no solutions.
    p2 = p * p
    if not S2plus or m % p2:
        return 0

    #  Check some conditions for no bad-type II solutions to exist
    if (NZvec is not None) and (len(Set(S2plus).intersection(Set(NZvec))) == 0):
        return 0

    #  Check that the form is primitive...      WHY IS THIS NECESSARY?
    if not S0:
        print(" Using Q = " + str(self))
        print(" and p = " + str(p))
        raise RuntimeError("the form is not primitive")

    #  DIAGNOSTIC
    verbose("\n Entering BII routine ")
    verbose(" S0 is " + str(S0))
    verbose(" S1 is " + str(S1))
    verbose(" S2plus is " + str(S2plus))
    verbose(" m = " + str(m) + "   p = " + str(p))

    #  Make the form Qnew for the reduction procedure:
    #  -----------------------------------------------
    Qnew = deepcopy(self)    # TODO:  DO THIS WITHOUT A copy()
    for i in range(n):
        if i in S2plus:
            Qnew[i, i] = Qnew[i, i] / p2
            if (p == 2) and (i < n - 1):
                Qnew[i, i + 1] = Qnew[i, i + 1] / p2

    #  DIAGNOSTIC
    verbose("\n\n Check of Bad-type II reduction: \n")
    verbose(" Q is " + str(self))
    verbose(" Qnew is " + str(Qnew))

    #  Perform the reduction formula
    Zvec_geq_2 = list(Set([i for i in Zvec if i in S2plus]))
    if NZvec is None:
        NZvec_geq_2 = NZvec
    else:
        NZvec_geq_2 = list(Set([i for i in NZvec if i in S2plus]))

    diff = Qnew.local_density_congruence(p, m / p2, Zvec_geq_2, NZvec_geq_2)
    diff -= Qnew.local_density_congruence(p, m / p2, S2plus, NZvec_geq_2)
    return QQ(p**(len(S2plus) + 2 - n)) * diff


def local_bad_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the Bad-type local density of `Q` representing
    `m` at `p`, allowing certain congruence conditions mod `p`.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_bad_density_congruence(2, 1, None, None)
        0
        sage: Q.local_bad_density_congruence(2, 2, None, None)
        1
        sage: Q.local_bad_density_congruence(2, 4, None, None)
        0
        sage: Q.local_bad_density_congruence(3, 1, None, None)
        0
        sage: Q.local_bad_density_congruence(3, 6, None, None)
        0
        sage: Q.local_bad_density_congruence(3, 9, None, None)
        0
        sage: Q.local_bad_density_congruence(3, 27, None, None)
        0

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,3,9,9])
        sage: Q.local_bad_density_congruence(3, 1, None, None)
        0
        sage: Q.local_bad_density_congruence(3, 3, None, None)
        4/3
        sage: Q.local_bad_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_bad_density_congruence(3, 9, None, None)
        4/27
        sage: Q.local_bad_density_congruence(3, 18, None, None)
        4/9
        sage: Q.local_bad_density_congruence(3, 27, None, None)
        8/27
    """
    return self.local_badI_density_congruence(p, m, Zvec, NZvec) + self.local_badII_density_congruence(p, m, Zvec, NZvec)

########################################################
#  local_density and local_density_congruence routines #
########################################################


def local_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the local density of `Q` representing `m` at `p`,
    allowing certain congruence conditions mod `p`.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_density_congruence(p=2, m=1, Zvec=None, NZvec=None)
        1
        sage: Q.local_density_congruence(p=3, m=1, Zvec=None, NZvec=None)
        8/9
        sage: Q.local_density_congruence(p=5, m=1, Zvec=None, NZvec=None)
        24/25
        sage: Q.local_density_congruence(p=7, m=1, Zvec=None, NZvec=None)
        48/49
        sage: Q.local_density_congruence(p=11, m=1, Zvec=None, NZvec=None)
        120/121

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_density_congruence(2, 1, None, None)
        1
        sage: Q.local_density_congruence(2, 2, None, None)
        1
        sage: Q.local_density_congruence(2, 4, None, None)
        3/2
        sage: Q.local_density_congruence(3, 1, None, None)
        2/3
        sage: Q.local_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_density_congruence(3, 9, None, None)
        14/9
        sage: Q.local_density_congruence(3, 27, None, None)
        2

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,3,9,9])
        sage: Q.local_density_congruence(3, 1, None, None)
        2
        sage: Q.local_density_congruence(3, 3, None, None)
        4/3
        sage: Q.local_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_density_congruence(3, 9, None, None)
        2/9
        sage: Q.local_density_congruence(3, 18, None, None)
        4/9
    """
    return self.local_good_density_congruence(p, m, Zvec, NZvec) \
        + self.local_zero_density_congruence(p, m, Zvec, NZvec) \
        + self.local_bad_density_congruence(p, m, Zvec, NZvec)


def local_primitive_density_congruence(self, p, m, Zvec=None, NZvec=None):
    """
    Find the primitive local density of `Q` representing
    `m` at `p`, allowing certain congruence conditions mod `p`.

    .. NOTE::

        The following routine is not used internally, but is included for consistency.

    INPUT:

    - ``self`` -- quadratic form `Q`, assumed to be block diagonal and `p`-integral

    - ``p`` -- a prime number

    - ``m`` -- integer

    - ``Zvec``, ``NZvec`` -- non-repeating lists of integers in ``range(self.dim())`` or ``None``

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.local_primitive_density_congruence(p=2, m=1, Zvec=None, NZvec=None)
        1
        sage: Q.local_primitive_density_congruence(p=3, m=1, Zvec=None, NZvec=None)
        8/9
        sage: Q.local_primitive_density_congruence(p=5, m=1, Zvec=None, NZvec=None)
        24/25
        sage: Q.local_primitive_density_congruence(p=7, m=1, Zvec=None, NZvec=None)
        48/49
        sage: Q.local_primitive_density_congruence(p=11, m=1, Zvec=None, NZvec=None)
        120/121

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
        sage: Q.local_primitive_density_congruence(2, 1, None, None)
        1
        sage: Q.local_primitive_density_congruence(2, 2, None, None)
        1
        sage: Q.local_primitive_density_congruence(2, 4, None, None)
        1
        sage: Q.local_primitive_density_congruence(3, 1, None, None)
        2/3
        sage: Q.local_primitive_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_primitive_density_congruence(3, 9, None, None)
        4/3
        sage: Q.local_primitive_density_congruence(3, 27, None, None)
        4/3

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,3,9,9])
        sage: Q.local_primitive_density_congruence(3, 1, None, None)
        2
        sage: Q.local_primitive_density_congruence(3, 3, None, None)
        4/3
        sage: Q.local_primitive_density_congruence(3, 6, None, None)
        4/3
        sage: Q.local_primitive_density_congruence(3, 9, None, None)
        4/27
        sage: Q.local_primitive_density_congruence(3, 18, None, None)
        4/9
        sage: Q.local_primitive_density_congruence(3, 27, None, None)
        8/27
        sage: Q.local_primitive_density_congruence(3, 81, None, None)
        8/27
        sage: Q.local_primitive_density_congruence(3, 243, None, None)
        8/27
    """
    return self.local_good_density_congruence(p, m, Zvec, NZvec) \
        + self.local_bad_density_congruence(p, m, Zvec, NZvec)
