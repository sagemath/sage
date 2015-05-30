"""
Local Field Invariants
"""
#*****************************************************************************
#       Copyright (C) 2007 William Stein and Jonathan Hanke
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


###########################################################################
## TO DO: Add routines for hasse invariants at all places, anisotropic
## places, is_semi_definite, and support for number fields.
###########################################################################


import copy

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.rings.arith import prime_divisors, hilbert_symbol
from sage.functions.all import sgn
from sage.rings.fraction_field import FractionField
from sage.matrix.matrix_space import MatrixSpace


## Routines to compute local (p-adic) invariants of a quadratic form Q:
## (Note: Here Q is the matrix so that Q(x) = x^t * Q * x.)
## --------------------------------------------------------------------

def rational_diagonal_form(self, return_matrix=False):
    """
    Returns a diagonal form equivalent to Q over the fraction field of
    its defining ring.  If the return_matrix is True, then we return
    the transformation matrix performing the diagonalization as the
    second argument.

    INPUT:
        none

    OUTPUT:
        Q -- the diagonalized form of this quadratic form
        (optional) T -- matrix which diagonalizes Q (over it's fraction field)

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [0,1,-1])
        sage: Q
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 0 1 ]
        [ * -1 ]
        sage: Q.rational_diagonal_form()
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ -2 0 ]
        [ * 1/8 ]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.rational_diagonal_form(return_matrix=True)
        (
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]                                                          ,
        <BLANKLINE>
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        )

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
    """
    n = self.dim()
    Q = copy.deepcopy(self)
    Q.__init__(FractionField(self.base_ring()), self.dim(), self.coefficients())
    MS = MatrixSpace(Q.base_ring(), n, n)
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


    ## Return the appropriate output
    if return_matrix:
        return Q, T
    else:
        return Q





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

    ## TO DO: Cache this result?

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

    EXAMPLES:
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
    """
    Computes the Hasse invariant at a prime `p`, as given on p55 of
    Cassels's book.  If Q is diagonal with coefficients `a_i`, then the
    (Cassels) Hasse invariant is given by

    .. math::

        c_p = \prod_{i < j} (a_i, a_j)_p

    where `(a,b)_p` is the Hilbert symbol at `p`.  The underlying
    quadratic form must be non-degenerate over `Q_p` for this to make
    sense.

    WARNING: This is different from the O'Meara Hasse invariant, which
    allows `i <= j` in the product.  That is given by the method
    hasse_invariant__OMeara(p).

    NOTE: We should really rename this hasse_invariant__Cassels(), and
    set hasse_invariant() as a front-end to it.


    INPUT:
        `p` -- a prime number > 0

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
    """
    Computes the O'Meara Hasse invariant at a prime `p`, as given on
    p167 of O'Meara's book.  If Q is diagonal with coefficients `a_i`,
    then the (Cassels) Hasse invariant is given by

    .. math::

        c_p = \prod_{i <= j} (a_i, a_j)_p

    where `(a,b)_p` is the Hilbert symbol at `p`.

    WARNING: This is different from the (Cassels) Hasse invariant, which
    only allows `i < j` in the product.  That is given by the method
    hasse_invariant(p).


    INPUT:
        `p` -- a prime number > 0

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

        sage: Q=DiagonalQuadraticForm(ZZ,[1,-1,-1])
        sage: [Q.hasse_invariant(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]
        sage: [Q.hasse_invariant__OMeara(p) for p in prime_range(20)]
        [-1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: K.<a>=NumberField(x^2-23)
        sage: Q=DiagonalQuadraticForm(K,[-a,a+2])
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
    """
    Checks if the quadratic form is a sum of hyperbolic planes over
    the p-adic numbers Q_p.

    REFERENCES:
        This criteria follows from Cassels's "Rational Quadratic Forms":
            - local invariants for hyperbolic plane (Lemma 2.4, p58)
            - direct sum formulas (Lemma 2.3 on p58)

    INPUT:
        `p` -- a prime number > 0

    OUTPUT:
        boolean

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.is_hyperbolic("infinity")
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
    if self.dim() % 2 != 0:
        return False

    ## True for the zero form
    if self.dim == 0:
        return True

    ## Compare local invariants
    ## (Note: since the dimension is even, the extra powers of 2 in
    ##        self.det() := Det(2*Q) don't affect the answer!)
    m = ZZ(self.dim() // 2)
    if p == "infinity":
        return (self.signature() == 0)

    elif p == 2:
        return QQ(self.det() * (-1)**m).is_padic_square(p) and (self.hasse_invariant(p) == (-1)**m)    ## Actually, this -1 is the Hilbert symbol (-1,-1)_p

    else:
        return QQ(self.det() * (-1)**m).is_padic_square(p) and (self.hasse_invariant(p) == 1)



def is_anisotropic(self, p):
    """
    Checks if the quadratic form is anisotropic over the p-adic numbers `Q_p`.

    INPUT:
        `p` -- a prime number > 0

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
    n = self.dim()
    D = self.det()

    ## TO DO: Should check that p is prime

    if (n >= 5):
        return False;

    if (n == 4):
        return ( QQ(D).is_padic_square(p) and (self.hasse_invariant(p) == - hilbert_symbol(-1,-1,p)) )

    if (n == 3):
        return (self.hasse_invariant(p) != hilbert_symbol(-1, -D, p))

    if (n == 2):
        return (not QQ(-D).is_padic_square(p))

    if (n == 1):
        return (self[0,0] != 0)

    raise NotImplementedError("Oops!  We haven't established a convention for 0-dim'l quadratic forms... =(")


def is_isotropic(self, p):
    """
    Checks if Q is isotropic over the p-adic numbers `Q_p`.

    INPUT:
        `p` -- a prime number > 0

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
    Returns a list with all of the anisotropic primes of the quadratic form.


    INPUT:
        None

    OUTPUT:
        Returns a list of prime numbers >0.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.anisotropic_primes()
        [2]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.anisotropic_primes()
        [2]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.anisotropic_primes()
        []

    """

    ## Look at all prime divisors of 2 * Det(Q) to find the anisotropic primes...
    possible_primes = prime_divisors(2 * self.det())
    AnisoPrimes = []

    ## DIAGNSOTIC
    #print " Possible anisotropic primes are: " + str(possible_primes)

    for p in possible_primes:
        if (self.is_anisotropic(p)):
            AnisoPrimes += [p]

    ## DIAGNSOTIC
    #print " leaving anisotropic_primes..."

    return AnisoPrimes




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
    M = self.matrix()


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

    EXAMPLES:
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
        I = range(r)
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


