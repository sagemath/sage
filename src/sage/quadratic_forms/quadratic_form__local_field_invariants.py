"""
Local Field Invariants

This contains routines to compute local (p-adic) invariants of
quadratic forms over the rationals.
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein and Jonathan Hanke
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

###########################################################################
## TO DO: Add routines for hasse invariants at all places, anisotropic
## places, is_semi_definite, and support for number fields.
###########################################################################


from copy import deepcopy
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.arith.all import prime_divisors, hilbert_symbol
from sage.functions.all import sgn
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method


def rational_diagonal_form(self, return_matrix=False):
    """
    Return a diagonal form equivalent to the given quadratic from
    over the fraction field of its defining ring.

    INPUT:

    - ``return_matrix`` -- (boolean, default: False) also return the
      transformation matrix

    OUTPUT: either ``D`` (if ``return_matrix`` is false) or ``(D,T)``
    (if ``return_matrix`` is true) where

    - ``D`` -- the diagonalized form of this quadratic form.

    - ``T`` -- transformation matrix. This is such that
      ``T.transpose() * self.matrix() * T`` gives ``D.matrix()``.

    Both ``D`` and ``T`` are defined over the fraction field of the
    base ring of the given form.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [0,1,-1])
        sage: Q
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 0 1 ]
        [ * -1 ]
        sage: Q.rational_diagonal_form()
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1/4 0 ]
        [ * -1 ]

    If we start with a diagonal form, we get back the same form defined
    over the fraction field::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.rational_diagonal_form()
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]

    In the following example, we check the consistency of the
    transformation matrix::

        sage: Q = QuadraticForm(ZZ, 4, range(10))
        sage: D, T = Q.rational_diagonal_form(return_matrix=True)
        sage: D
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ -1/16 0 0 0 ]
        [ * 4 0 0 ]
        [ * * 13 0 ]
        [ * * * 563/52 ]
        sage: T
        [     1      0     11 149/26]
        [  -1/8      1     -2 -10/13]
        [     0      0      1 -29/26]
        [     0      0      0      1]
        sage: T.transpose() * Q.matrix() * T
        [  -1/8      0      0      0]
        [     0      8      0      0]
        [     0      0     26      0]
        [     0      0      0 563/26]
        sage: D.matrix()
        [  -1/8      0      0      0]
        [     0      8      0      0]
        [     0      0     26      0]
        [     0      0      0 563/26]

    ::

        sage: Q1 = QuadraticForm(ZZ, 4, [1, 1, 0, 0, 1, 0, 0, 1, 0, 18])
        sage: Q1
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 1 0 0 ]
        [ * 1 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]
        sage: Q1.rational_diagonal_form(return_matrix=True)
        (
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3/4 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]                                                         ,
        <BLANKLINE>
        [   1 -1/2    0    0]
        [   0    1    0    0]
        [   0    0    1    0]
        [   0    0    0    1]
        )

    PARI returns a singular transformation matrix for this case::

        sage: Q = QuadraticForm(QQ, 2, [1/2, 1, 1/2])
        sage: Q.rational_diagonal_form()
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1/2 0 ]
        [ * 0 ]

    This example cannot be computed by PARI::

        sage: Q = QuadraticForm(RIF, 4, range(10))
        sage: Q.__pari__()
        Traceback (most recent call last):
        ...
        TypeError
        sage: Q.rational_diagonal_form()
        Quadratic form in 4 variables over Real Interval Field with 53 bits of precision with coefficients:
        [ 5 0.?e-14 0.?e-13 0.?e-13 ]
        [ * -0.05000000000000? 0.?e-12 0.?e-12 ]
        [ * * 13.00000000000? 0.?e-10 ]
        [ * * * 10.8269230769? ]

    TESTS:

    Changing the output quadratic form does not affect the caching::

        sage: Q, T = Q1.rational_diagonal_form(return_matrix=True)
        sage: Q[0,0] = 13
        sage: Q1.rational_diagonal_form()
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3/4 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]

    The transformation matrix is immutable::

        sage: T[0,0] = 13
        Traceback (most recent call last):
        ...
        ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
    """
    Q, T = self._rational_diagonal_form_and_transformation()
    T.set_immutable()

    # Quadratic forms do not support immutability, so we need to make
    # a copy to be safe.
    Q = deepcopy(Q)

    if return_matrix:
        return Q, T
    else:
        return Q


@cached_method
def _rational_diagonal_form_and_transformation(self):
    """
    Return a diagonal form equivalent to the given quadratic from and
    the corresponding transformation matrix. This is over the fraction
    field of the base ring of the given quadratic form.

    OUTPUT: a tuple ``(D,T)`` where

    - ``D`` -- the diagonalized form of this quadratic form.

    - ``T`` -- transformation matrix. This is such that
      ``T.transpose() * self.matrix() * T`` gives ``D.matrix()``.

    Both ``D`` and ``T`` are defined over the fraction field of the
    base ring of the given form.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, [1, 1, 0, 0, 1, 0, 0, 1, 0, 18])
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 1 0 0 ]
        [ * 1 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]
        sage: Q._rational_diagonal_form_and_transformation()
        (
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3/4 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]                                                         ,
        <BLANKLINE>
        [   1 -1/2    0    0]
        [   0    1    0    0]
        [   0    0    1    0]
        [   0    0    0    1]
        )
    """
    n = self.dim()
    K = self.base_ring().fraction_field()
    Q = self.base_change_to(K)
    MS = MatrixSpace(K, n, n)

    try:
        # Try PARI if the type is supported
        pariself = self.__pari__()
        # Check that conversion back works
        MS(pariself.sage())
    except Exception:
        pass
    else:
        R = pariself.qfgaussred()
        # Diagonal matrix
        D = MS()
        for i in range(n):
            D[i,i] = R[i,i]
        Q = Q.parent()(D)
        # Transformation matrix (inverted)
        T = MS(R.sage())
        for i in range(n):
            T[i,i] = K.one()
        try:
            return Q, ~T
        except ZeroDivisionError:
            # Singular case is not fully supported by PARI
            pass

    # General case if conversion to/from PARI failed
    T = MS(1)

    ## Clear the entries one row at a time.
    for i in range(n):

        ## Deal with rows where the diagonal entry is zero.
        if Q[i,i] == 0:

            ## Look for a non-zero entry and use it to make the diagonal non-zero (if it exists)
            for j in range(i+1, n):
                if Q[i,j] != 0:
                    temp = MS(1)
                    if Q[i,j] + Q[j,j] == 0:
                        temp[j, i] = -1
                    else:
                        temp[j, i] = 1

                    ## Apply the transformation
                    Q = Q(temp)
                    T = T * temp
                    break

        ## Create a matrix which deals with off-diagonal entries (all at once for each row)
        temp = MS(1)
        for j in range(i+1, n):
            if Q[i,j] != 0:
                temp[i,j] = -Q[i,j] / (Q[i,i] * 2)    ## This should only occur when Q[i,i] != 0, which the above step guarantees.

        Q = Q(temp)
        T = T * temp

    return Q, T


def signature_vector(self):
    """
    Returns the triple `(p, n, z)` of integers where

    - `p` = number of positive eigenvalues
    - `n` = number of negative eigenvalues
    - `z` = number of zero eigenvalues

    for the symmetric matrix associated to Q.

    INPUT:

    (none)

    OUTPUT:

    a triple of integers >= 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,0,0,-4])
        sage: Q.signature_vector()
        (1, 1, 2)

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,-3,-4])
        sage: Q.signature_vector()
        (2, 2, 0)

    ::

        sage: Q = QuadraticForm(ZZ, 4, range(10)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 2 3 ]
        [ * 4 5 6 ]
        [ * * 7 8 ]
        [ * * * 9 ]
        sage: Q.signature_vector()
        (3, 1, 0)

    """
    diag = self.rational_diagonal_form()
    p = 0
    n = 0
    z = 0
    for i in range(diag.dim()):
        if diag[i,i] > 0:
            p += 1
        elif diag[i,i] < 0:
            n += 1
        else:
            z += 1

    return (p, n, z)



def signature(self):
    """
    Returns the signature of the quadratic form, defined as:

       number of positive eigenvalues - number of negative eigenvalues

    of the matrix of the quadratic form.

    INPUT:

        None

    OUTPUT:

        an integer

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,0,0,-4,3,11,3])
        sage: Q.signature()
        3

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,2,-3,-4])
        sage: Q.signature()
        0

    ::

        sage: Q = QuadraticForm(ZZ, 4, range(10)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 2 3 ]
        [ * 4 5 6 ]
        [ * * 7 8 ]
        [ * * * 9 ]
        sage: Q.signature()
        2

    """
    (p, n, z) = self.signature_vector()
    return p - n


def hasse_invariant(self, p):
    r"""
    Computes the Hasse invariant at a prime `p` or at infinity, as given on p55 of
    Cassels's book.  If Q is diagonal with coefficients `a_i`, then the
    (Cassels) Hasse invariant is given by

    .. MATH::

        c_p = \prod_{i < j} (a_i, a_j)_p

    where `(a,b)_p` is the Hilbert symbol at `p`.  The underlying
    quadratic form must be non-degenerate over `Q_p` for this to make
    sense.

    .. WARNING::

        This is different from the O'Meara Hasse invariant, which
        allows `i <= j` in the product.  That is given by the method
        hasse_invariant__OMeara(p).

    .. NOTE::

        We should really rename this hasse_invariant__Cassels(), and
        set hasse_invariant() as a front-end to it.

    INPUT:

    - `p` -- a prime number > 0 or `-1` for the infinite place

    OUTPUT:

    1 or -1

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
        sage: Q.rational_diagonal_form()
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1 0 ]
        [ * 2 ]
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1])
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1,5])
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: K.<a>=NumberField(x^2-23)
        sage: Q=DiagonalQuadraticForm(K,[-a,a+2])
        sage: [Q.hasse_invariant(p) for p in K.primes_above(19)]
        [-1, 1]
    """
    ## TO DO: Need to deal with the case n=1 separately somewhere!

    Diag = self.rational_diagonal_form()
    R = Diag.base_ring()

    ## DIAGNOSTIC
    #print "\n Q = " + str(self)
    #print "\n Q diagonalized at p = " + str(p) + " gives " + str(Diag)

    hasse_temp = 1
    n = Diag.dim()

    if R == QQ:
        for j in range(n-1):
            for k in range(j+1, n):
                hasse_temp = hasse_temp * hilbert_symbol(Diag[j,j], Diag[k,k], p)

    else:
        for j in range(n-1):
            for k in range(j+1, n):
                hasse_temp = hasse_temp * R.hilbert_symbol(Diag[j,j], Diag[k,k], p)

    return hasse_temp


def hasse_invariant__OMeara(self, p):
    r"""
    Compute the O'Meara Hasse invariant at a prime `p`.

    This is defined on
    p167 of O'Meara's book. If Q is diagonal with coefficients `a_i`,
    then the (Cassels) Hasse invariant is given by

    .. MATH::

        c_p = \prod_{i <= j} (a_i, a_j)_p

    where `(a,b)_p` is the Hilbert symbol at `p`.

    .. WARNING::

        This is different from the (Cassels) Hasse invariant, which
        only allows `i < j` in the product.  That is given by the method
        hasse_invariant(p).

    INPUT:

    - `p` -- a prime number > 0 or `-1` for the infinite place

    OUTPUT:

    1 or -1

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
        sage: Q.rational_diagonal_form()
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1 0 ]
        [ * 2 ]
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1])
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ,[1,-1,-1])
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: K.<a>=NumberField(x^2-23)
        sage: Q = DiagonalQuadraticForm(K,[-a,a+2])
        sage: [Q.hasse_invariant__OMeara(p) for p in K.primes_above(19)]
        [1, 1]
    """
    ## TO DO: Need to deal with the case n=1 separately somewhere!

    Diag = self.rational_diagonal_form()
    R = Diag.base_ring()

    ## DIAGNOSTIC
    #print "\n Q = " + str(self)
    #print "\n Q diagonalized at p = " + str(p) + " gives " + str(Diag)

    hasse_temp = 1
    n = Diag.dim()
    if R == QQ:
        for j in range(n):
            for k in range(j, n):
                hasse_temp = hasse_temp * hilbert_symbol(Diag[j,j], Diag[k,k], p)

    else:
        for j in range(n):
            for k in range(j, n):
                hasse_temp = hasse_temp * R.hilbert_symbol(Diag[j,j], Diag[k,k], p)

    return hasse_temp


def is_hyperbolic(self, p):
    r"""
    Check if the quadratic form is a sum of hyperbolic planes over
    the `p`-adic numbers `\QQ_p` or over the real numbers `\RR`.

    REFERENCES:

    This criteria follows from Cassels's "Rational Quadratic Forms":

    - local invariants for hyperbolic plane (Lemma 2.4, p58)
    - direct sum formulas (Lemma 2.3, p58)

    INPUT:

    - `p` -- a prime number > 0 or `-1` for the infinite place

    OUTPUT:

    boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.is_hyperbolic(-1)
        False
        sage: Q.is_hyperbolic(2)
        False
        sage: Q.is_hyperbolic(3)
        False
        sage: Q.is_hyperbolic(5)     ## Here -1 is a square, so it's true.
        True
        sage: Q.is_hyperbolic(7)
        False
        sage: Q.is_hyperbolic(13)    ## Here -1 is a square, so it's true.
        True
    """
    ## False for odd-dim'l forms
    if self.dim() % 2:
        return False

    ## True for the zero form
    if not self.dim():
        return True

    ## Compare local invariants
    ## Note: since the dimension is even, the extra powers of 2 in
    ##        self.det() := Det(2*Q) don't affect the answer!
    m = ZZ(self.dim() // 2)
    if p == -1:
        return self.signature() == 0

    if p == 2:
        return (QQ(self.det() * (-1) ** m).is_padic_square(p) and
                self.hasse_invariant(p) ==
                (-1) ** m.binomial(2))  # here -1 is hilbert_symbol(-1,-1,2)

    return (QQ(self.det() * (-1) ** m).is_padic_square(p) and
            self.hasse_invariant(p) == 1)


def is_anisotropic(self, p):
    r"""
    Check if the quadratic form is anisotropic over the p-adic numbers `\QQ_p` or `\RR`.

    INPUT:

    - `p` -- a prime number > 0 or `-1` for the infinite place

    OUTPUT:

    boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.is_anisotropic(2)
        True
        sage: Q.is_anisotropic(3)
        True
        sage: Q.is_anisotropic(5)
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1])
        sage: Q.is_anisotropic(2)
        False
        sage: Q.is_anisotropic(3)
        False
        sage: Q.is_anisotropic(5)
        False

    ::

        sage: [DiagonalQuadraticForm(ZZ, [1, -least_quadratic_nonresidue(p)]).is_anisotropic(p)  for p in prime_range(3, 30)]
        [True, True, True, True, True, True, True, True, True]

    ::

        sage: [DiagonalQuadraticForm(ZZ, [1, -least_quadratic_nonresidue(p), p, -p*least_quadratic_nonresidue(p)]).is_anisotropic(p)  for p in prime_range(3, 30)]
        [True, True, True, True, True, True, True, True, True]
    """
    ## TO DO: Should check that p is prime
    if p == -1:
        return self.is_definite()

    n = self.dim()
    D = self.det()

    if n >= 5:
        return False

    if n == 4:
        return (QQ(D).is_padic_square(p) and
                (self.hasse_invariant(p) == - hilbert_symbol(-1, -1, p)))

    if n == 3:
        return self.hasse_invariant(p) != hilbert_symbol(-1, -D, p)

    if n == 2:
        return not QQ(-D).is_padic_square(p)

    if n == 1:
        return self[0, 0] != 0

    raise NotImplementedError("Oops!  We haven't established a convention for 0-dim'l quadratic forms... =(")


def is_isotropic(self, p):
    """
    Checks if Q is isotropic over the p-adic numbers `Q_p` or `RR`.

    INPUT:

    - `p` -- a prime number > 0 or `-1` for the infinite place

    OUTPUT:

    boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.is_isotropic(2)
        False
        sage: Q.is_isotropic(3)
        False
        sage: Q.is_isotropic(5)
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1])
        sage: Q.is_isotropic(2)
        True
        sage: Q.is_isotropic(3)
        True
        sage: Q.is_isotropic(5)
        True

    ::

        sage: [DiagonalQuadraticForm(ZZ, [1, -least_quadratic_nonresidue(p)]).is_isotropic(p)  for p in prime_range(3, 30)]
        [False, False, False, False, False, False, False, False, False]

    ::

        sage: [DiagonalQuadraticForm(ZZ, [1, -least_quadratic_nonresidue(p), p, -p*least_quadratic_nonresidue(p)]).is_isotropic(p)  for p in prime_range(3, 30)]
        [False, False, False, False, False, False, False, False, False]

    """
    return not self.is_anisotropic(p)


def anisotropic_primes(self):
    """
    Return a list with all of the anisotropic primes of the quadratic form.

    The infinite place is denoted by `-1`.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.anisotropic_primes()
        [2, -1]

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.anisotropic_primes()
        [2, -1]

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.anisotropic_primes()
        [-1]
    """
    # Look at all prime divisors of 2 * Det(Q) to find the
    # anisotropic primes...
    possible_primes = prime_divisors(2 * self.det()) + [-1]
    return [p for p in possible_primes if self.is_anisotropic(p)]


def compute_definiteness(self):
    """
    Computes whether the given quadratic form is positive-definite,
    negative-definite, indefinite, degenerate, or the zero form.

    This caches one of the following strings in self.__definiteness_string:
    "pos_def", "neg_def", "indef", "zero", "degenerate".  It is called
    from all routines like:

        is_positive_definite(), is_negative_definite(), is_indefinite(), etc.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:

        QuadraticForm

    OUTPUT:

        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.compute_definiteness()
        sage: Q.is_positive_definite()
        True
        sage: Q.is_negative_definite()
        False
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [])
        sage: Q.compute_definiteness()
        sage: Q.is_positive_definite()
        True
        sage: Q.is_negative_definite()
        True
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,0,-1])
        sage: Q.compute_definiteness()
        sage: Q.is_positive_definite()
        False
        sage: Q.is_negative_definite()
        False
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        False

    """
    ## Sanity Check
    if not ((self.base_ring() == ZZ) or (self.base_ring() == QQ) or (self.base_ring() == RR)):
        raise NotImplementedError("Oops!  We can only check definiteness over ZZ, QQ, and RR for now.")

    ## Some useful variables
    n = self.dim()

    ## Deal with the zero-diml form
    if n == 0:
        self.__definiteness_string = "zero"
        return

    sig_pos, sig_neg, sig_zer = self.signature_vector()

    ## Determine and cache the definiteness string
    if sig_zer > 0:
        self.__definiteness_string = "degenerate"
        return
    elif sig_neg == n:
        self.__definiteness_string = "neg_def"
        return
    elif sig_pos == n:
        self.__definiteness_string = "pos_def"
        return
    else:
        self.__definiteness_string = "indefinite"
        return


def compute_definiteness_string_by_determinants(self):
    """
    Compute the (positive) definiteness of a quadratic form by looking
    at the signs of all of its upper-left subdeterminants.  See also
    self.compute_definiteness() for more documentation.

    INPUT:

        None

    OUTPUT:

        string describing the definiteness

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.compute_definiteness_string_by_determinants()
        'pos_def'

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [])
        sage: Q.compute_definiteness_string_by_determinants()
        'zero'

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,0,-1])
        sage: Q.compute_definiteness_string_by_determinants()
        'degenerate'

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-1])
        sage: Q.compute_definiteness_string_by_determinants()
        'indefinite'

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-1])
        sage: Q.compute_definiteness_string_by_determinants()
        'neg_def'

    """
    ## Sanity Check
    if not ((self.base_ring() == ZZ) or (self.base_ring() == QQ) or (self.base_ring() == RR)):
        raise NotImplementedError("Oops!  We can only check definiteness over ZZ, QQ, and RR for now.")

    ## Some useful variables
    n = self.dim()
    M = self.matrix()


    ## Deal with the zero-diml form
    if n == 0:
        return "zero"



    ## Deal with degenerate forms
    if self.det() == 0:
        return "degenerate"


    ## Check the sign of the ratios of consecutive determinants of the upper triangular r x r submatrices
    first_coeff = self[0,0]
    for r in range(1,n+1):
        I = list(range(r))
        new_det = M.matrix_from_rows_and_columns(I, I).det()

        ## Check for a (non-degenerate) zero -- so it's indefinite
        if new_det == 0:
            return "indefinite"


        ## Check for a change of signs in the upper r x r submatrix -- so it's indefinite
        if sgn(first_coeff)**r != sgn(new_det):
            return "indefinite"


    ## Here all ratios of determinants have the correct sign, so the matrix is (pos or neg) definite.
    if first_coeff > 0:
        return "pos_def"
    else:
        return "neg_def"





def is_positive_definite(self):
    """
    Determines if the given quadratic form is positive-definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:

        None

    OUTPUT:

        boolean -- True or False

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5])
        sage: Q.is_positive_definite()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_positive_definite()
        False

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except AttributeError:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "pos_def") or (def_str == "zero")




def is_negative_definite(self):
    """
    Determines if the given quadratic form is negative-definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:

        None

    OUTPUT:

        boolean -- True or False

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_negative_definite()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_negative_definite()
        False

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except AttributeError:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "neg_def") or (def_str == "zero")



def is_indefinite(self):
    """
    Determines if the given quadratic form is indefinite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is not considered indefinite.

    INPUT:

        None

    OUTPUT:

        boolean -- True or False

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_indefinite()
        False

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_indefinite()
        True

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except AttributeError:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return def_str == "indefinite"


def is_definite(self):
    """
    Determines if the given quadratic form is (positive or negative) definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered indefinite.

    INPUT:

        None

    OUTPUT:

        boolean -- True or False

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_definite()
        True

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_definite()
        False

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except AttributeError:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "pos_def") or (def_str == "neg_def") or (def_str == "zero")


