"""
Binary Recurrence Sequences

This class implements several methods relating to general linear binary
recurrence sequences, including a sieve to find perfect powers in integral
linear binary recurrence sequences.

EXAMPLES::

    sage: R = BinaryRecurrenceSequence(1,1)        #the Fibonacci sequence
    sage: R(137)        #the 137th term of the Fibonacci sequence
    19134702400093278081449423917
    sage: R(137) == fibonacci(137)
    True
    sage: [R(i) % 4 for i in range(12)]
    [0, 1, 1, 2, 3, 1, 0, 1, 1, 2, 3, 1]
    sage: R.period(4)        #the period of the fibonacci sequence modulo 4
    6
    sage: R.pthpowers(2, 10**10)        # long time (7 seconds) -- in fact these are all squares, c.f. [BMS06]
    [0, 1, 2, 12]

    sage: S = BinaryRecurrenceSequence(8,1) #a Lucas sequence
    sage: S.period(73)
    148
    sage: S(5) % 73 == S(5 +148) %73
    True
    sage: S.pthpowers(3,10**10)    # long time (3 seconds) -- provably finds the indices of all 3rd powers less than 10^10
    [0, 1, 2]

    sage: T = BinaryRecurrenceSequence(2,0,1,2)
    sage: [T(i) for i in range(10)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    sage: T.is_degenerate()
    True
    sage: T.is_geometric()
    True
    sage: T.pthpowers(7,10**30)
    Traceback (most recent call last):
    ...
    ValueError: The degenerate binary recurrence sequence is geometric or quasigeometric and has many pth powers.


AUTHORS:

- Isabel Vogt (2013): initial version

See [SV2013]_, [BMS2006]_, and [SS1983]_.
"""

# ****************************************************************************
#       Copyright (C) 2013 Isabel Vogt <ivogt161@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.arith.all import lcm, next_prime, is_prime, next_prime_power, legendre_symbol
from sage.functions.log import log
from sage.misc.functional import sqrt


class BinaryRecurrenceSequence(SageObject):

    """
    Create a linear binary recurrence sequence defined by initial conditions
    `u_0` and `u_1` and recurrence relation `u_{n+2} = b*u_{n+1}+c*u_n`.

    INPUT:

    - ``b`` -- an integer (partially determining the recurrence relation)

    - ``c`` -- an integer (partially determining the recurrence relation)

    - ``u0`` -- an integer (the 0th term of the binary recurrence sequence)

    - ``u1`` -- an integer (the 1st term of the binary recurrence sequence)


    OUTPUT:

    - An integral linear binary recurrence sequence defined by ``u0``, ``u1``, and `u_{n+2} = b*u_{n+1}+c*u_n`

    .. SEEALSO::

        :func:`fibonacci`, :func:`lucas_number1`, :func:`lucas_number2`

    EXAMPLES::

        sage: R = BinaryRecurrenceSequence(3,3,2,1)
        sage: R
        Binary recurrence sequence defined by: u_n = 3 * u_{n-1} + 3 * u_{n-2};
        With initial conditions: u_0 = 2, and u_1 = 1

    """

    def __init__(self, b, c, u0=0, u1=1):

        """
        See :class:`BinaryRecurrenceSequence` for full documentation.

        EXAMPLES::

            sage: R = BinaryRecurrenceSequence(3,3,2,1)
            sage: R
            Binary recurrence sequence defined by: u_n = 3 * u_{n-1} + 3 * u_{n-2};
            With initial conditions: u_0 = 2, and u_1 = 1

            sage: R = BinaryRecurrenceSequence(1,1)
            sage: loads(R.dumps()) == R
            True

        """
        self.b = b
        self.c = c
        self.u0 = u0
        self.u1 = u1
        self._period_dict = {}    #dictionary to cache the period of a sequence for future lookup
        self._PGoodness = {} #dictionary to cache primes that are "good" by some prime power
        self._ell = 1 #variable that keeps track of the last prime power to be used as a goodness

    def __repr__(self):
        """
        Give string representation of the class.

        EXAMPLES::

            sage: R = BinaryRecurrenceSequence(3,3,2,1)
            sage: R
            Binary recurrence sequence defined by: u_n = 3 * u_{n-1} + 3 * u_{n-2};
            With initial conditions: u_0 = 2, and u_1 = 1

        """
        return 'Binary recurrence sequence defined by: u_n = ' + str(self.b) + ' * u_{n-1} + ' + str(self.c) + ' * u_{n-2};\nWith initial conditions: u_0 = ' + str(self.u0) + ', and u_1 = ' + str(self.u1)

    def __eq__(self, other):
        """
        Compare two binary recurrence sequences.

        EXAMPLES::

            sage: R = BinaryRecurrenceSequence(3,3,2,1)
            sage: S = BinaryRecurrenceSequence(3,3,2,1)
            sage: R == S
            True

            sage: T = BinaryRecurrenceSequence(3,3,2,2)
            sage: R == T
            False
        """
        return (self.u0 == other.u0) and (self.u1 == other.u1) and (self.b == other.b) and (self.c == other.c)

    def __call__(self, n, modulus=0):
        """
        Give the nth term of a binary recurrence sequence, possibly mod some modulus.

        INPUT:

        - ``n`` -- an integer (the index of the term in the binary recurrence sequence)

        - ``modulus`` -- a natural number (optional --  default value is 0)

        OUTPUT:

        - An integer (the nth term of the binary recurrence sequence modulo ``modulus``)

        EXAMPLES::

            sage: R = BinaryRecurrenceSequence(3,3,2,1)
            sage: R(2)
            9
            sage: R(101)
            16158686318788579168659644539538474790082623100896663971001
            sage: R(101,12)
            9
            sage: R(101)%12
            9
        """
        R = Integers(modulus)
        F = matrix(R, [[0, 1], [self.c, self.b]])
        # F*[u_{n}, u_{n+1}]^T = [u_{n+1}, u_{n+2}]^T (T indicates transpose).
        v = vector(R, [self.u0, self.u1])
        return list(F**n * v)[0]

    def is_degenerate(self):
        """
        Decide whether the binary recurrence sequence is degenerate.

        Let `\\alpha` and `\\beta` denote the roots of the characteristic polynomial
        `p(x) = x^2-bx -c`.  Let `a = u_1-u_0\\beta/(\\beta - \\alpha)` and
        `b = u_1-u_0\\alpha/(\\beta - \\alpha)`.  The sequence is, thus, given by
        `u_n = a \\alpha^n - b\\beta^n`.  Then we say that the sequence is nondegenerate
        if and only if `a*b*\\alpha*\\beta \\neq 0` and `\\alpha/\\beta` is not a
        root of unity.

        More concretely, there are 4 classes of degeneracy, that can all be formulated
        in terms of the matrix `F = [[0,1], [c, b]]`.

        - `F` is singular --  this corresponds to ``c`` = 0, and thus `\\alpha*\\beta = 0`. This sequence is geometric after term ``u0`` and so we call it ``quasigeometric``.

        - `v = [[u_0], [u_1]]` is an eigenvector of `F` -- this corresponds to a ``geometric`` sequence with `a*b = 0`.

        - `F` is nondiagonalizable -- this corresponds to `\\alpha = \\beta`.  This sequence will be the point-wise product of an arithmetic and geometric sequence.

        - `F^k` is scaler, for some `k>1` -- this corresponds to `\\alpha/\\beta` a `k` th root of unity. This sequence is a union of several geometric sequences, and so we again call it ``quasigeometric``.

        EXAMPLES::

            sage: S = BinaryRecurrenceSequence(0,1)
            sage: S.is_degenerate()
            True
            sage: S.is_geometric()
            False
            sage: S.is_quasigeometric()
            True

            sage: R = BinaryRecurrenceSequence(3,-2)
            sage: R.is_degenerate()
            False

            sage: T = BinaryRecurrenceSequence(2,-1)
            sage: T.is_degenerate()
            True
            sage: T.is_arithmetic()
            True
        """
        D = self.b**2 + 4 * self.c
        if D != 0:
            if D.is_square():
                A = sqrt(D)
            else:
                K = QuadraticField(D, 'x')
                A = K.gen()

            aa = (self.u1 - self.u0*(self.b + A)/2)/(A)        #called `a` in Docstring
            bb = (self.u1 - self.u0*(self.b - A)/2)/(A)        #called `b` in Docstring

            #(b+A)/2 is called alpha in Docstring, (b-A)/2 is called beta in Docstring

            if self.b != A:
                if ((self.b+A)/(self.b-A))**6 == 1:
                    return True
            else:
                return True

            return aa*bb*(self.b + A)*(self.b - A) == 0

        return True

    def is_geometric(self):
        """
        Decide whether the binary recurrence sequence is geometric - ie a geometric sequence.

        This is a subcase of a degenerate binary recurrence sequence, for which `ab=0`, i.e.
        `u_{n}/u_{n-1}=r` for some value of `r`.

        See :meth:`is_degenerate` for a description of
        degeneracy and definitions of `a` and `b`.

        EXAMPLES::

            sage: S = BinaryRecurrenceSequence(2,0,1,2)
            sage: [S(i) for i in range(10)]
            [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
            sage: S.is_geometric()
            True
        """
        #If [u_0, u_1]^T is an eigenvector for the incrementation matrix F = [[0,1],[c,b]], then the sequence
        #is geometric, ie we can write u_n = a*r^n for some a and r.

        #We decide if u0, u1, u2 = b*u1+c*u0 are in geometric progression by whether u1^2 = (b*u1+c*u0)*u0

        return (self.u1)**2 == (self.b*self.u1 + self.c*self.u0)*self.u0

    def is_quasigeometric(self):
        """
        Decide whether the binary recurrence sequence is degenerate and similar to a geometric sequence,
        i.e. the union of multiple geometric sequences, or geometric after term ``u0``.

        If `\\alpha/\\beta` is a `k` th root of unity, where `k>1`, then necessarily `k = 2, 3, 4, 6`.
        Then `F = [[0,1],[c,b]` is diagonalizable, and `F^k = [[\\alpha^k, 0], [0,\\beta^k]]` is scaler
        matrix.  Thus for all values of `j` mod `k`, the `j` mod `k` terms of `u_n` form a geometric
        series.

        If `\\alpha` or `\\beta` is zero, this implies that `c=0`.  This is the case when `F` is
        singular.  In this case, `u_1, u_2, u_3, ...` is geometric.

        EXAMPLES::

            sage: S = BinaryRecurrenceSequence(0,1)
            sage: [S(i) for i in range(10)]
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
            sage: S.is_quasigeometric()
            True

            sage: R = BinaryRecurrenceSequence(3,0)
            sage: [R(i) for i in range(10)]
            [0, 1, 3, 9, 27, 81, 243, 729, 2187, 6561]
            sage: R.is_quasigeometric()
            True
        """
        # First test if F is singular... i.e. beta = 0
        if self.c == 0:
            return True

        # Otherwise test if alpha/beta is a root of unity that is not 1
        D = self.b**2 + 4 * self.c
        if D != 0:    # thus alpha/beta != 1
            if D.is_square():
                A = sqrt(D)
            else:
                K = QuadraticField(D, 'x')
                A = K.gen()
            if ((self.b+A)/(self.b-A))**6 == 1:
                return True

        return False

    def is_arithmetic(self):
        """
        Decide whether the sequence is degenerate and an arithmetic sequence.

        The sequence is arithmetic if and only if `u_1 - u_0 = u_2 - u_1 = u_3 - u_2`.

        This corresponds to the matrix `F = [[0,1],[c,b]]` being nondiagonalizable
        and `\\alpha/\\beta = 1`.

        EXAMPLES::

            sage: S = BinaryRecurrenceSequence(2,-1)
            sage: [S(i) for i in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: S.is_arithmetic()
            True
        """
        return (self(1) - self(0) == self(2) - self(1) == self(3) - self(2))

    def period(self, m):
        """
        Return the period of the binary recurrence sequence modulo
        an integer ``m``.

        If `n_1` is congruent to `n_2` modulo ``period(m)``, then `u_{n_1}` is
        is congruent to `u_{n_2}` modulo ``m``.

        INPUT:

        - ``m`` -- an integer (modulo which the period of the recurrence relation is calculated).

        OUTPUT:

        - The integer (the period of the sequence modulo m)

        EXAMPLES:

        If `p = \\pm 1 \\mod 5`, then the period of the Fibonacci sequence
        mod `p` is `p-1` (c.f. Lemma 3.3 of [BMS2006]_).

        ::

            sage: R = BinaryRecurrenceSequence(1,1)
            sage: R.period(31)
            30

            sage: [R(i) % 4 for i in range(12)]
            [0, 1, 1, 2, 3, 1, 0, 1, 1, 2, 3, 1]
            sage: R.period(4)
            6

        This function works for degenerate sequences as well.

        ::

            sage: S = BinaryRecurrenceSequence(2,0,1,2)
            sage: S.is_degenerate()
            True
            sage: S.is_geometric()
            True
            sage: [S(i) % 17 for i in range(16)]
            [1, 2, 4, 8, 16, 15, 13, 9, 1, 2, 4, 8, 16, 15, 13, 9]
            sage: S.period(17)
            8

        .. NOTE:: The answer is cached.
        """

        #If we have already computed the period mod m, then we return the stored value.

        if m in self._period_dict:
            return self._period_dict[m]

        else:
            R = Integers(m)
            A = matrix(R, [[0, 1], [self.c, self.b]])
            w = vector(R, [self.u0, self.u1])
            Fac = list(m.factor())
            Periods = {}

            #To compute the period mod m, we compute the least integer n such that A^n*w == w.  This necessarily
            #divides the order of A as a matrix in GL_2(Z/mZ).

            #We compute the period modulo all distinct prime powers dividing m, and combine via the lcm.
            #To compute the period mod p^e, we first compute the order mod p.  Then the period mod p^e
            #must divide p^{4e-4}*period(p), as the subgroup of matrices mod p^e, which reduce to
            #the identity mod p is of order (p^{e-1})^4.  So we compute the period mod p^e by successively
            #multiplying the period mod p by powers of p.

            for i in Fac:
                p = i[0]
                e = i[1]
                #first compute the period mod p
                if p in self._period_dict:
                    perp = self._period_dict[p]
                else:
                    F = A.change_ring(GF(p))
                    v = w.change_ring(GF(p))
                    FF = F**(p-1)
                    p1fac = list((p-1).factor())

                    #The order of any matrix in GL_2(F_p) either divides p(p-1) or (p-1)(p+1).
                    #The order divides p-1 if it is diagonalizable.  In any case, det(F^(p-1))=1,
                    #so if tr(F^(p-1)) = 2, then it must be triangular of the form [[1,a],[0,1]].
                    #The order of the subgroup of matrices of this form is p, so the order must divide
                    #p(p-1) -- in fact it must be a multiple of p.  If this is not the case, then the
                    #order divides (p-1)(p+1).  As the period divides the order of the matrix in GL_2(F_p),
                    #these conditions hold for the period as well.

                    #check if the order divides (p-1)
                    if FF*v == v:
                        M = p-1
                        Mfac = p1fac

                    #check if the trace is 2, then the order is a multiple of p dividing p*(p-1)
                    elif FF.trace() == 2:
                        M = p-1
                        Mfac = p1fac
                        F = F**p        #replace F by F^p as now we only need to determine the factor dividing (p-1)

                    #otherwise it will divide (p+1)(p-1)
                    else:
                        M = (p+1)*(p-1)
                        p2fac = list((p+1).factor())        #factor the (p+1) and (p-1) terms separately and then combine for speed
                        Mfac_dic = {}
                        for i0, i1 in list(p1fac + p2fac):
                            if i0 not in Mfac_dic:
                                Mfac_dic[i0] = i1
                            else:
                                Mfac_dic[i0] += i1
                        Mfac = list(Mfac_dic.items())

                    #Now use a fast order algorithm to compute the period.  We know that the period divides
                    #M = i_1*i_2*...*i_l where the i_j denote not necessarily distinct prime factors.  As
                    #F^M*v == v, for each i_j, if F^(M/i_j)*v == v, then the period divides (M/i_j).  After
                    #all factors have been iterated over, the result is the period mod p.

                    Mfac = list(Mfac)
                    C = []

                    #expand the list of prime factors so every factor is with multiplicity 1

                    for i0, i1 in Mfac:
                        for j in range(i1):
                            C.append(i0)

                    Mfac = C
                    n = M
                    for ii in Mfac:
                        b = n // ii
                        if F**b * v == v:
                            n = b
                    perp = n

                #Now compute the period mod p^e by stepping up by multiples of p
                F = A.change_ring(Integers(p**e))
                v = w.change_ring(Integers(p**e))
                FF = F**perp
                if FF*v == v:
                    perpe = perp
                else:
                    tries = 0
                    while True:
                        tries += 1
                        FF = FF**p
                        if FF*v == v:
                            perpe = perp*p**tries
                            break
                Periods[p] = perpe

            #take the lcm of the periods mod all distinct primes dividing m
            period = 1
            for p in Periods:
                period = lcm(Periods[p], period)

            self._period_dict[m] = period        #cache the period mod m
            return period

    def pthpowers(self, p, Bound):
        """
        Find the indices of proveably all pth powers in the recurrence sequence bounded by Bound.

        Let `u_n` be a binary recurrence sequence.  A ``p`` th power in `u_n` is a solution
        to `u_n = y^p` for some integer `y`.  There are only finitely many ``p`` th powers in
        any recurrence sequence [SS1983]_.

        INPUT:

        - ``p`` - a rational prime integer (the fixed p in `u_n = y^p`)

        - ``Bound`` - a natural number (the maximum index `n` in `u_n = y^p` that is checked).

        OUTPUT:

        - A list of the indices of all ``p`` th powers less bounded by ``Bound``.  If the sequence is degenerate and there are many ``p`` th powers, raises ``ValueError``.

        EXAMPLES::

            sage: R = BinaryRecurrenceSequence(1,1)        #the Fibonacci sequence
            sage: R.pthpowers(2, 10**10)        # long time (7 seconds) -- in fact these are all squares, c.f. [BMS2006]_
            [0, 1, 2, 12]

            sage: S = BinaryRecurrenceSequence(8,1) #a Lucas sequence
            sage: S.pthpowers(3,10**10)    # long time (3 seconds) -- provably finds the indices of all 3rd powers less than 10^10
            [0, 1, 2]

            sage: Q = BinaryRecurrenceSequence(3,3,2,1)
            sage: Q.pthpowers(11,10**10)          # long time (7.5 seconds)
            [1]

        If the sequence is degenerate, and there are no ``p`` th powers, returns `[]`.  Otherwise, if
        there are many ``p`` th powers, raises ``ValueError``.

        ::

            sage: T = BinaryRecurrenceSequence(2,0,1,2)
            sage: T.is_degenerate()
            True
            sage: T.is_geometric()
            True
            sage: T.pthpowers(7,10**30)
            Traceback (most recent call last):
            ...
            ValueError: The degenerate binary recurrence sequence is geometric or quasigeometric and has many pth powers.

            sage: L = BinaryRecurrenceSequence(4,0,2,2)
            sage: [L(i).factor() for i in range(10)]
            [2, 2, 2^3, 2^5, 2^7, 2^9, 2^11, 2^13, 2^15, 2^17]
            sage: L.is_quasigeometric()
            True
            sage: L.pthpowers(2,10**30)
            []

        .. NOTE::

            This function is primarily optimized in the range where
            ``Bound`` is much larger than ``p``.
        """
        #Thanks to Jesse Silliman for helpful conversations!

        #Reset the dictionary of good primes, as this depends on p
        self._PGoodness = {}
        #Starting lower bound on good primes
        self._ell = 1

        #If the sequence is geometric, then the `n`th term is `a*r^n`.  Thus the
        #property of being a ``p`` th power is periodic mod ``p``.  So there are either
        #no ``p`` th powers if there are none in the first ``p`` terms, or many if there
        #is at least one in the first ``p`` terms.

        if self.is_geometric() or self.is_quasigeometric():
            no_powers = True
            for i in range(1, 6*p+1):
                if _is_p_power(self(i), p):
                    no_powers = False
                    break
            if no_powers:
                if _is_p_power(self.u0, p):
                    return [0]
                return []
            else:
                raise ValueError("The degenerate binary recurrence sequence is geometric or quasigeometric and has many pth powers.")

        #If the sequence is degenerate without being geometric or quasigeometric, there
        #may be many ``p`` th powers or no ``p`` th powers.

        elif (self.b**2+4*self.c) == 0:

            #This is the case if the matrix F is not diagonalizable, ie b^2 +4c = 0, and alpha/beta = 1.

            alpha = self.b/2

            #In this case, u_n = u_0*alpha^n + (u_1 - u_0*alpha)*n*alpha^(n-1) = alpha^(n-1)*(u_0 +n*(u_1 - u_0*alpha)),
            #that is, it is a geometric term (alpha^(n-1)) times an arithmetic term (u_0 + n*(u_1-u_0*alpha)).

            #Look at classes n = k mod p, for k = 1,...,p.

            for k in range(1, p + 1):

                #The linear equation alpha^(k-1)*u_0 + (k+pm)*(alpha^(k-1)*u1 - u0*alpha^k)
                #must thus be a pth power.  This is a linear equation in m, namely, A + B*m, where

                A = (alpha**(k-1)*self.u0 + k*(alpha**(k-1)*self.u1 - self.u0*alpha**k))
                B = p*(alpha**(k-1)*self.u1 - self.u0*alpha**k)

                #This linear equation represents a pth power iff A is a pth power mod B.

                if _is_p_power_mod(A, p, B):
                    raise ValueError("The degenerate binary recurrence sequence has many pth powers.")
            return []

        #We find ``p`` th powers using an elementary sieve.  Term `u_n` is a ``p`` th
        #power if and only if it is a ``p`` th power modulo every prime `\\ell`.  This condition
        #gives nontrivial information if ``p`` divides the order of the multiplicative group of
        #`\\Bold(F)_{\\ell}`, i.e. if `\\ell` is ` 1 \mod{p}`, as then only `1/p` terms are ``p`` th
        #powers modulo `\\ell``.

        #Thus, given such an `\\ell`, we get a set of necessary congruences for the index modulo the
        #the period of the sequence mod `\\ell`.  Then we intersect these congruences for many primes
        #to get a tight list modulo a growing modulus.  In order to keep this step manageable, we
        #only use primes `\\ell` that have particularly smooth periods.

        #Some congruences in the list will remain as the modulus grows.  If a congruence remains through
        #7 rounds of increasing the modulus, then we check if this corresponds to a perfect power (if
        #it does, we add it to our list of indices corresponding to ``p`` th powers).  The rest of the congruences
        #are transient and grow with the modulus.  Once the smallest of these is greater than the bound,
        #the list of known indices corresponding to ``p`` th powers is complete.

        else:

            if Bound < 3 * p:

                powers = []
                ell = p + 1

                while not is_prime(ell):
                    ell += p

                F = GF(ell)
                a0 = F(self.u0)
                a1 = F(self.u1) #a0 and a1 are variables for terms in sequence
                bf, cf = F(self.b), F(self.c)

                for n in range(Bound): # n is the index of the a0

                    #Check whether a0 is a perfect power mod ell
                    if _is_p_power_mod(a0, p, ell):
                        #if a0 is a perfect power mod ell, check if nth term is ppower
                        if _is_p_power(self(n), p):
                            powers.append(n)

                    a0, a1 = a1, bf*a1 + cf*a0        #step up the variables

            else:

                powers = []        #documents the indices of the sequence that provably correspond to pth powers
                cong = [0]        #list of necessary congruences on the index for it to correspond to pth powers
                Possible_count = {}    #keeps track of the number of rounds a congruence lasts in cong

                #These parameters are involved in how we choose primes to increase the modulus
                qqold = 1        #we believe that we know complete information coming from primes good by qqold
                M1 = 1            #we have congruences modulo M1, this may not be the tightest list
                M2 = p            #we want to move to have congruences mod M2
                qq = 1            #the largest prime power divisor of M1 is qq

                #This loop ups the modulus.
                while True:

                    #Try to get good data mod M2

                    #patience of how long we should search for a "good prime"
                    patience = 0.01 * _estimated_time(lcm(M2, p*next_prime_power(qq)), M1, len(cong), p)
                    tries = 0

                    #This loop uses primes to get a small set of congruences mod M2.
                    while True:

                        #only proceed if took less than patience time to find the next good prime
                        ell = _next_good_prime(p, self, qq, patience, qqold)
                        if ell:

                            #gather congruence data for the sequence mod ell, which will be mod period(ell) = modu
                            cong1, modu = _find_cong1(p, self, ell)

                            CongNew = []        #makes a new list from cong that is now mod M = lcm(M1, modu) instead of M1
                            M = lcm(M1, modu)
                            for k in range(M // M1):
                                for i in cong:
                                    CongNew.append(k * M1 + i)
                            cong = set(CongNew)

                            M1 = M

                            killed_something = False        #keeps track of when cong1 can rule out a congruence in cong

                            #CRT by hand to gain speed
                            for i in list(cong):
                                if not (i % modu in cong1):        #congruence in cong is inconsistent with any in cong1
                                    cong.remove(i)            #remove that congruence
                                    killed_something = True

                            if M1 == M2:
                                if not killed_something:
                                    tries += 1
                                    if tries == 2:            #try twice to rule out congruences
                                        cong = list(cong)
                                        qqold = qq
                                        qq = next_prime_power(qq)
                                        M2 = lcm(M2, p * qq)
                                        break

                        else:
                            qq = next_prime_power(qq)
                            M2 = lcm(M2, p * qq)
                            cong = list(cong)
                            break

                    #Document how long each element of cong has been there
                    for i in cong:
                        if i in Possible_count:
                            Possible_count[i] += 1
                        else:
                            Possible_count[i] = 1

                    #Check how long each element has persisted, if it is for at least 7 cycles,
                    #then we check to see if it is actually a perfect power
                    for i in Possible_count:
                        if Possible_count[i] == 7:
                            n = Integer(i)
                            if n < Bound:
                                if _is_p_power(self(n), p):
                                    powers.append(n)

                    #check for a contradiction
                    if len(cong) > len(powers):
                        if cong[len(powers)] > Bound:
                            break
                    elif M1 > Bound:
                        break

            return powers


def _prime_powers(N):
    """
    Find the prime powers dividing ``N``.

    In other words, if `N = q_1^(e_1)q_2^(e_2)...q_n^(e_n)`, it returns
    `[q_1^(e_1),q_2^(e_2),...,q_n^(e_n)]`.

    INPUT:

    - ``N`` -- an integer

    OUTPUT:

    - A list of the prime powers dividing N.

    EXAMPLES::

        sage: sage.combinat.binary_recurrence_sequences._prime_powers(124656)
        [3, 16, 49, 53]

        sage: sage.combinat.binary_recurrence_sequences._prime_powers(65537)
        [65537]
    """
    return sorted(i**j for i, j in N.factor())


def _largest_ppower_divisor(N):
    """
    Find the largest prime power divisor of N.

    INPUT:

    - ``N`` -- an integer

    OUTPUT:

    The largest prime power dividing ``N``.

    EXAMPLES::

        sage: sage.combinat.binary_recurrence_sequences._largest_ppower_divisor(124656)
        53
        sage: sage.combinat.binary_recurrence_sequences._largest_ppower_divisor(65537)
        65537
    """
    return _prime_powers(N)[-1]


def _goodness(n, R, p):
    """
    Return the goodness of ``n`` for the sequence ``R`` and the prime ``p`` -- that is the largest
    non-``p`` prime power dividing ``period(n)``.

    INPUT:

    - ``n`` --  an integer

    - ``R`` -- an object in the class ``BinaryRecurrenceSequence``

    - ``p`` -- a rational prime

    OUTPUT:

    - An integer which is the "goodness" of ``n``, i.e. the largest non-``p`` prime power dividing ``period(n)``.

    EXAMPLES::

        sage: R = BinaryRecurrenceSequence(11,2)
        sage: sage.combinat.binary_recurrence_sequences._goodness(89,R,7)
        11

        sage: R = BinaryRecurrenceSequence(1,1)
        sage: sage.combinat.binary_recurrence_sequences._goodness(13,R,7)
        4
        sage: R.period(13)        #the period of R mod 13 is divisible by 7
        28
    """
    # The period of R mod ell
    K = R.period(n)
    return _largest_ppower_divisor(K // K.gcd(p))


def _next_good_prime(p, R, qq, patience, qqold):

    """
    Find the next prime `\\ell` which is good by ``qq`` but not by ``qqold``, 1 mod ``p``, and for which
    ``b^2+4*c`` is a square mod `\\ell`, for the sequence ``R`` if it is possible in runtime patience.

    INPUT:

    - ``p`` -- a prime

    - ``R`` -- an object in the class ``BinaryRecurrenceSequence``

    - ``qq`` -- a perfect power

    - ``patience`` -- a real number

    - ``qqold`` --  a perfect power less than or equal to ``qq``

    OUTPUT:

    - A prime `\\ell` such that `\\ell` is 1 mod ``p``, ``b^2+4*c`` is a square mod `\\ell` and the period of `\\ell` has ``goodness`` by ``qq`` but not ``qqold``, if patience has not be surpased.  Otherwise ``False``.


    EXAMPLES::

        sage: R = BinaryRecurrenceSequence(1,1)
        sage: sage.combinat.binary_recurrence_sequences._next_good_prime(7,R,1,100,1)        #ran out of patience to search for good primes
        False
        sage: sage.combinat.binary_recurrence_sequences._next_good_prime(7,R,2,100,1)
        29
        sage: sage.combinat.binary_recurrence_sequences._next_good_prime(7,R,2,100,2)        #ran out of patience, as qqold == qq, so no primes work
        False

    """

    #We are looking for pth powers in R.
    #Our primes must be good by qq, but not qqold.
    #We only allow patience number of iterations to find a good prime.

    #The variable _ell for R keeps track of the last "good" prime returned
    #that was not found from the dictionary _PGoodness

    #First, we check to see if we have already computed the goodness of a prime that fits
    #our requirement of being good by qq but not by qqold.  This is stored in the _PGoodness
    #dictionary.

    #Then if we have, we return the smallest such prime and delete it from the list.  If not, we
    #search through patience number of primes R._ell to find one good by qq but not qqold.  If it is
    #not good by either qqold or qq, then we add this prime to R._PGoodness under its goodness.

    #Possible_Primes keeps track of possible primes satisfying our goodness requirements we might return
    Possible_Primes = []

    #check to see if anything in R._PGoodness fits our goodness requirements
    for j in R._PGoodness:
        if (qqold < j <= qq) and len(R._PGoodness[j]):
            Possible_Primes.append(R._PGoodness[j][0])

    #If we found good primes, we take the smallest
    if Possible_Primes:
        q = min(Possible_Primes)
        n = _goodness(q, R, p)
        del R._PGoodness[n][0]    #if we are going to use it, then we delete it from R._PGoodness
        return q

    #If nothing is already stored in R._PGoodness, we start (from where we left off at R._ell) checking
    #for good primes.  We only tolerate patience number of tries before giving up.
    else:
        i = 0
        while i < patience:
            i += 1
            R._ell = next_prime(R._ell)

            #we require that R._ell is 1 mod p, so that p divides the order of the multiplicative
            #group mod R._ell, so that not all elements of GF(R._ell) are pth powers.
            if R._ell % p == 1:

                #requiring that b^2 + 4c is a square in GF(R._ell) ensures that the period mod R._ell
                #divides R._ell - 1
                if legendre_symbol(R.b**2+4*R.c, R._ell) == 1:

                    N = _goodness(R._ell, R, p)

                    #proceed only if R._ell satisfies the goodness requirements
                    if qqold < N <= qq:
                        return R._ell

                    #if we do not use the prime, we store it in R._PGoodness
                    else:
                        if N in R._PGoodness:
                            R._PGoodness[N].append(R._ell)
                        else:
                            R._PGoodness[N] = [R._ell]

        return False


def _is_p_power_mod(a, p, N):
    """
    Determine if ``a`` is a ``p`` th power modulo ``N``.

    By the CRT, this is equivalent to the condition that ``a`` be a ``p`` th power mod all
    distinct prime powers dividing ``N``.  For each of these, we use the strong statement of
    Hensel's lemma to lift ``p`` th powers mod `q` or `q^2` or `q^3` to ``p`` th powers mod `q^e`.

    INPUT:

    - ``a`` -- an integer

    - ``p`` -- a rational prime number

    - ``N`` -- a positive integer

    OUTPUT:

    - True if ``a`` is a ``p`` th power modulo ``N``; False otherwise.

    EXAMPLES::

        sage: sage.combinat.binary_recurrence_sequences._is_p_power_mod(2**3,7,29)
        False
        sage: sage.combinat.binary_recurrence_sequences._is_p_power_mod(2**3,3,29)
        True

    """

    #By the chinese remainder theorem, we can answer this question by examining whether
    #a is a pth power mod q^e, for all distinct prime powers q^e dividing N.

    for q, e in N.factor():

        #If a = q^v*x, with

        v = a.valuation(q)

        #then if v>=e, a is congruent to 0 mod q^e and is thus a pth power trivially.

        if v >= e:
            continue

        # otherwise, it can only be a pth power if v is a multiple of p.
        if v % p:
            return False

        #in this cse it is a pth power if x is a pth power mod q^(e-v), so let x = aa,
        #and (e-v) = ee:

        aa = a/q**v
        ee = e - v

        #The above steps are equivalent to the statement that we may assume a and qq are
        #relatively prime, if we replace a with aa and e with ee.  Now we must determine when
        #aa is a pth power mod q^ee for (aa,q)=1.

        #If q != p, then by Hensel's lemma, we may lift a pth power mod q, to a pth power
        #mod q^2, etc.

        if q != p:

            #aa is necessarily a pth power mod q if p does not divide the order of the multiplicative
            #group mod q, ie if q is not 1 mod p.

            if q % p == 1:

                #otherwise aa if a pth power mod q iff aa^(q-1)/p == 1

                if GF(q)(aa)**((q-1)/p) != 1:
                    return False

        #If q = p and ee = 1, then everything is a pth power p by Fermat's little theorem.

        elif ee > 1:

            #We use the strong statement of Hensel's lemma, which implies that if p is odd
            #and aa is a pth power mod p^2, then aa is a pth power mod any higher power of p

            if p % 2:

                #ZZ/(p^2)ZZ^\times is abstractly isomorphic to ZZ/(p)ZZ cross ZZ/(p-1)ZZ. then
                #aa is a pth power mod p^2 if (aa)^(p*(p-1)/p) == 1, ie if aa^(p-1) == 1.

                if Integers(p**2)(aa)**(p-1) != 1:
                    return False

            #Otherwise, p=2.  By the strong statement of Hensel's lemma, if aa is a pth power
            #mod p^3, then it is a pth power mod higher powers of p.  So we need only check if it
            #is a pth power mod p^2 and p^3.

            elif ee == 2:

                #all odd squares a 1 mod 4

                if aa % 4 != 1:
                    return False

            #all odd squares are 1 mod 8

            elif aa % 8 != 1:
                return False

    return True


def _estimated_time(M2, M1, length, p):

    """
    Find the estimated time to extend congruences mod M1 to consistent congruences mod M2.

    INPUT:

    - ``M2`` -- an integer (the new modulus)

    - ``M1`` -- an integer (the old modulus)

    - ``length`` -- a list (the current length of the list of congruences mod ``M1``)

    - ``p`` --  a prime

    OUTPUT:

    - The estimated run time of the "CRT" step to combine consistent congruences.

    EXAMPLES::

        sage: sage.combinat.binary_recurrence_sequences._estimated_time(2**4*3**2*5*7*11*13*17, 2**4*3**2*5*7*11*13, 20, 7)
        106.211159309421

    """

    #The heuristic run time of the CRT step to go from modulus M1 to M2

    #length is the current length of cong

    Q = p * log(M2) #Size of our primes.
    NPrimes = log(M2/M1) / log(Q) #The number of primes

    return (length * (Q/p)**NPrimes).n()


#Find the list of necessary congruences for the index n of binary recurrence
#sequence R using the fact that the reduction mod ell must be a pth power
def _find_cong1(p, R, ell):
    """
    Find the list of permissible indices `n` for which `u_n = y^p` mod ``ell``.

    INPUT:

    - ``p`` -- a prime number

    - ``R`` -- an object in class :class:`BinaryRecurrenceSequence`

    - ``ell`` -- a prime number

    OUTPUT:

    - A list of permissible values of `n` modulo ``period(ell)`` and the integer ``period(ell)``.

    EXAMPLES::

        sage: R = BinaryRecurrenceSequence(1,1)
        sage: sage.combinat.binary_recurrence_sequences._find_cong1(7, R, 29)
        ([0, 1, 2, 12, 13], 14)
    """
    F = GF(ell)
    u0 = F(R.u0)
    u1 = F(R.u1)
    bf, cf = F(R.b), F(R.c)
    a0 = u0
    a1 = u1        #a0 and a1 are variables for terms in sequence

    #The set of pth powers mod ell
    PPowers = set(i**p for i in F)

    #The period of R mod ell
    modu = R.period(ell)

    #cong1 keeps track of congruences mod modu for the sequence mod ell
    cong1 = []

    for n in range(modu): # n is the index of the a0

        #Check whether a0 is a perfect power mod ell
        if a0 in PPowers:
            #if a0 is a perfect power mod ell, add the index
            #to the list of necessary congruences
            cong1.append(n)

        a0, a1 = a1, bf*a1 + cf*a0        #step up the variables

    cong1.sort()

    return cong1, modu


def _is_p_power(a, p):
    """
    Determine whether ``a`` is a perfect ``p`` th power.

    INPUT:

    - ``a`` -- an integer

    - ``p`` -- a prime number

    OUTPUT:

    - True if ``a`` is a ``p`` th power; else False.

    EXAMPLES::

        sage: sage.combinat.binary_recurrence_sequences._is_p_power(2**7,7)
        True
        sage: sage.combinat.binary_recurrence_sequences._is_p_power(2**7*3**2,7)
        False
    """
    return int(a**(1/p))**p == a
    # slower tentative ?
    # _, test = Integer(a).nth_root(p, truncate_mode=True)
    # return test
