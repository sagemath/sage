"""
Conway-Sloane masses
"""
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.all import kronecker_symbol, legendre_symbol, prime_divisors, is_prime, fundamental_discriminant
from sage.symbolic.constants import pi
from sage.misc.misc_c import prod
from sage.quadratic_forms.special_values import gamma__exact, zeta__exact, quadratic_L_function__exact


def parity(self, allow_rescaling_flag=True):
    """
    Return the parity ("even" or "odd") of an integer-valued quadratic
    form over `ZZ`, defined up to similitude/rescaling of the form so that
    its Jordan component of smallest scale is unimodular.  After this
    rescaling, we say a form is even if it only represents even numbers,
    and odd if it represents some odd number.

    If the 'allow_rescaling_flag' is set to False, then we require that
    the quadratic form have a Gram matrix with coefficients in `ZZ`, and
    look at the unimodular Jordan block to determine its parity.  This
    returns an error if the form is not integer-matrix, meaning that it
    has Jordan components at `p=2` which do not have an integer scale.

    We determine the parity by looking for a 1x1 block in the 0-th
    Jordan component, after a possible rescaling.

    INPUT:

        self -- a quadratic form with base_ring `ZZ`, which we may
                require to have integer Gram matrix.

    OUTPUT:

    One of the strings: "even" or "odd"

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [4, -2, 0, 2, 3, 2]); Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 4 -2 0 ]
        [ * 2 3 ]
        [ * * 2 ]
        sage: Q.parity()
        'even'

    ::

        sage: Q = QuadraticForm(ZZ, 3, [4, -2, 0, 2, 3, 1]); Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 4 -2 0 ]
        [ * 2 3 ]
        [ * * 1 ]
        sage: Q.parity()
        'even'

    ::

        sage: Q = QuadraticForm(ZZ, 3, [4, -2, 0, 2, 2, 2]); Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 4 -2 0 ]
        [ * 2 2 ]
        [ * * 2 ]
        sage: Q.parity()
        'even'

    ::

        sage: Q = QuadraticForm(ZZ, 3, [4, -2, 0, 2, 2, 1]); Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 4 -2 0 ]
        [ * 2 2 ]
        [ * * 1 ]
        sage: Q.parity()
        'odd'

    """
    ## Deal with 0-dim'l forms
    if self.dim() == 0:
        return "even"

    ## Identify the correct Jordan component to use.
    Jordan_list = self.jordan_blocks_by_scale_and_unimodular(2)
    scale_pow_list = [J[0]  for J in Jordan_list]
    min_scale_pow = min(scale_pow_list)
    if allow_rescaling_flag:
        ind = scale_pow_list.index(min_scale_pow)
    else:
        if min_scale_pow < 0:
            raise TypeError("Oops!  If rescaling is not allowed, then we require our form to have an integral Gram matrix.")
        ind = scale_pow_list.index(0)


    ## Find the component of scale (power) zero, and then look for an odd dim'l component.
    J0 = Jordan_list[ind]
    Q0 = J0[1]

    ## The lattice is even if there is no component of scale (power) 0
    if J0 is None:
        return "even"

    ## Look for a 1x1 block in the 0-th Jordan component (which by
    ## convention of the local_normal_form routine will appear first).
    if Q0.dim() == 1:
        return "odd"
    elif Q0[0,1] == 0:
        return "odd"
    else:
        return "even"



def is_even(self, allow_rescaling_flag=True):
    """
    Returns true iff after rescaling by some appropriate factor, the
    form represents no odd integers.  For more details, see parity().

    Requires that Q is defined over `ZZ`.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1, 0, 1])
        sage: Q.is_even()
        False
        sage: Q = QuadraticForm(ZZ, 2, [1, 1, 1])
        sage: Q.is_even()
        True

    """
    return self.parity(allow_rescaling_flag) == "even"


def is_odd(self, allow_rescaling_flag=True):
    """
    Returns true iff after rescaling by some appropriate factor, the
    form represents some odd integers.  For more details, see parity().

    Requires that Q is defined over `ZZ`.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1, 0, 1])
        sage: Q.is_odd()
        True
        sage: Q = QuadraticForm(ZZ, 2, [1, 1, 1])
        sage: Q.is_odd()
        False

    """
    return self.parity(allow_rescaling_flag) == "odd"



def conway_species_list_at_odd_prime(self, p):
    """
    Returns an integer called the 'species' which determines the type
    of the orthogonal group over the finite field `F_p`.

    This assumes that the given quadratic form is a unimodular Jordan
    block at an odd prime `p`.  When the dimension is odd then this
    number is always positive, otherwise it may be positive or
    negative (or zero, but that is considered positive by convention).

    Note: The species of a zero dim'l form is always 0+, so we
    interpret the return value of zero as positive here! =)

    INPUT:

        a positive prime number

    OUTPUT:

        a list of integers

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,10))
        sage: Q.conway_species_list_at_odd_prime(3)
        [6, 2, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,8))
        sage: Q.conway_species_list_at_odd_prime(3)
        [5, 2]
        sage: Q.conway_species_list_at_odd_prime(5)
        [-6, 1]

    """
    ## Sanity Check:
    if not ((p>2) and is_prime(p)):
        raise TypeError("Oops!  We are assuming that p is an odd positive prime number.")

    ## Deal with the zero-dim'l form
    if self.dim() == 0:
        return [0]

    ## List the (unscaled/unimodular) Jordan blocks by their scale power
    jordan_list = self.jordan_blocks_in_unimodular_list_by_scale_power(p)

    ## Make a list of species (including the two zero-dim'l forms missing at either end of the list of Jordan blocks)
    species_list = []
    for tmp_Q in jordan_list:

        ## Some useful variables
        n = tmp_Q.dim()
        d = tmp_Q.det()

        ## Determine the species
        if (n % 2 != 0):                            ## Deal with odd dim'l forms
            species = n
        elif (n % 4 == 2) and (p % 4 == 3):         ## Deal with even dim'l forms
            species = (-1) * legendre_symbol(d, p) * n
        else:
            species = legendre_symbol(d, p) * n

        ## Append the species to the list
        species_list.append(species)

    ## Return the species list
    return species_list



def conway_species_list_at_2(self):
    """
    Returns an integer called the 'species' which determines the type
    of the orthogonal group over the finite field `F_p`.

    This assumes that the given quadratic form is a unimodular Jordan
    block at an odd prime `p`.  When the dimension is odd then this
    number is always positive, otherwise it may be positive or
    negative.

    Note: The species of a zero dim'l form is always 0+, so we
    interpret the return value of zero as positive here! =)

    OUTPUT:

        a list of integers

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,10))
        sage: Q.conway_species_list_at_2()
        [1, 5, 1, 1, 1, 1]

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,8))
        sage: Q.conway_species_list_at_2()
        [1, 3, 1, 1, 1]

    """
    ## Some useful variables
    n = self.dim()
    d = self.det()

    ## Deal with the zero-dim'l form
    if n == 0:
        return 0

    ## List the (unscaled/unimodular) Jordan blocks by their scale power
    jordan_list = self.jordan_blocks_in_unimodular_list_by_scale_power(2)

    ## Make a list of species (including the two zero-dim'l forms missing at either end of the list of Jordan blocks)
    species_list = []

    if jordan_list[0].parity() == "odd":        ## Add an entry for the unlisted "-1" Jordan component as well.
        species_list.append(1)

    for i in range(len(jordan_list)):           ## Add an entry for each (listed) Jordan component

        ## Make the number 2*t in the C-S Table 1.
        d = jordan_list[i].dim()
        if jordan_list[i].is_even():
            two_t = d
        else:
            two_t = ZZ(2) * ((d-1) // 2)

        ## Determine if the form is bound
        if len(jordan_list) == 1:
            is_bound = False
        elif i == 0:
            is_bound = jordan_list[i+1].is_odd()
        elif i == len(jordan_list) - 1:
            is_bound = jordan_list[i-1].is_odd()
        else:
            is_bound = jordan_list[i-1].is_odd() or jordan_list[i+1].is_odd()

        ## Determine the species
        octane = jordan_list[i].conway_octane_of_this_unimodular_Jordan_block_at_2()
        if is_bound or (octane == 2) or (octane == 6):
            species = two_t + 1
        elif (octane == 0) or (octane == 1) or (octane == 7):
            species = two_t
        else:
            species = (-1) * two_t

        ## Append the species to the list
        species_list.append(species)


    if jordan_list[-1].is_odd():        ## Add an entry for the unlisted "s_max + 1" Jordan component as well.
        species_list.append(1)

    ## Return the species list
    return species_list




def conway_octane_of_this_unimodular_Jordan_block_at_2(self):
    """
    Determines the 'octane' of this full unimodular Jordan block at
    the prime `p=2`.  This is an invariant defined `(mod 8)`, ad.

    This assumes that the form is given as a block diagonal form with
    unimodular blocks of size <= 2 and the 1x1 blocks are all in the upper
    leftmost position.

    INPUT:

        none

    OUTPUT:

        an integer 0 <= x <= 7

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.conway_octane_of_this_unimodular_Jordan_block_at_2()
        0
        sage: Q = DiagonalQuadraticForm(ZZ, [1,5,13])
        sage: Q.conway_octane_of_this_unimodular_Jordan_block_at_2()
        3
        sage: Q = DiagonalQuadraticForm(ZZ, [3,7,13])
        sage: Q.conway_octane_of_this_unimodular_Jordan_block_at_2()
        7

    """
    ## Deal with 'even' forms
    if self.parity() == "even":
        d = self.Gram_matrix().det()
        if (d % 8 == 1) or (d % 8 == 7):
            return 0
        else:
            return 4

    ## Deal with 'odd' forms by diagonalizing, and then computing the octane.
    n = self.dim()
    u = self[0,0]
    tmp_diag_vec = [None  for i in range(n)]
    tmp_diag_vec[0] = u       ## This should be an odd integer!
    ind = 1                  ## The next index to diagonalize


    ## Use u to diagonalize the form -- WHAT ARE THE POSSIBLE LOCAL NORMAL FORMS?
    while ind < n:

        ## Check for a 1x1 block and diagonalize it
        if (ind == (n-1)) or (self[ind, ind+1] == 0):
            tmp_diag_vec[ind] = self[ind, ind]
            ind += 1

        ## Diagonalize the 2x2 block
        else:
            B = self[ind, ind+1]
            if (B % 2 != 0):
                raise RuntimeError("Oops, we expected the mixed term to be even! ")

            a = self[ind, ind]
            b = ZZ(B / ZZ(2))
            c = self[ind+1, ind+1]
            tmp_disc = b * b - a * c

            ## Perform the diagonalization
            if (tmp_disc % 8 == 1):                ## 2xy
                tmp_diag_vec[ind] = 1
                tmp_diag_vec[ind+1] = -1
                ind += 2
            elif(tmp_disc % 8 == 5):               ## 2x^2 + 2xy + 2y^2
                tmp_diag_vec[0] = 3*u
                tmp_diag_vec[ind] = -u
                tmp_diag_vec[ind+1] = -u
                ind += 2
                u = tmp_diag_vec[0]
            else:
                raise RuntimeError("Oops!  This should not happen -- the odd 2x2 blocks have disc 1 or 5 (mod 8).")

    ## Compute the octane
    octane = 0
    for a in tmp_diag_vec:
        if a % 4 == 1:
            octane += 1
        elif a % 4 == 3:
            octane += -1
        else:
            raise RuntimeError("Oops!  The diagonal elements should all be odd... =(")

    ## Return its value
    return octane % 8



def conway_diagonal_factor(self, p):
    """
    Computes the diagonal factor of Conway's `p`-mass.

    INPUT:

        `p` -- a prime number > 0

    OUTPUT:

        a rational number > 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,6))
        sage: Q.conway_diagonal_factor(3)
        81/256

    """
     ## Get the species list at p
    if p == 2:
        species_list = self.conway_species_list_at_2()
    else:
        species_list = self.conway_species_list_at_odd_prime(p)

    ## Evaluate the diagonal factor
    diag_factor = QQ(1)
    for s in species_list:
        if s == 0:
            pass
        elif s % 2 == 1:                   ## Note: Here always s > 0.
            diag_factor = diag_factor / (2 * prod([1 - QQ(p)**(-i)  for i in range(2, s, 2)]))
        else:
            diag_factor = diag_factor / (2 * prod([1 - QQ(p)**(-i)  for i in range(2, abs(s), 2)]))
            s_sign = ZZ(s / abs(s))
            diag_factor = diag_factor / (ZZ(1) - s_sign * QQ(p) ** ZZ(-abs(s) / ZZ(2)))

    ## Return the diagonal factor
    return diag_factor



def conway_cross_product_doubled_power(self, p):
    """
    Computes twice the power of p which evaluates the 'cross product'
    term in Conway's mass formula.

    INPUT:

        `p` -- a prime number > 0

    OUTPUT:

        a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,8))
        sage: Q.conway_cross_product_doubled_power(2)
        18
        sage: Q.conway_cross_product_doubled_power(3)
        10
        sage: Q.conway_cross_product_doubled_power(5)
        6
        sage: Q.conway_cross_product_doubled_power(7)
        6
        sage: Q.conway_cross_product_doubled_power(11)
        0
        sage: Q.conway_cross_product_doubled_power(13)
        0

    """
    doubled_power = 0
    dim_list = [J.dim()  for J in self.jordan_blocks_in_unimodular_list_by_scale_power(p)]
    for i in range(len(dim_list)):
        for j in range(i):
            doubled_power += (i-j) * dim_list[i] * dim_list[j]

    return doubled_power



def conway_type_factor(self):
    """
    This is a special factor only present in the mass formula when `p=2`.

    INPUT:

        none

    OUTPUT:

        a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1,8))
        sage: Q.conway_type_factor()
        4

    """
    jordan_list = self.jordan_blocks_in_unimodular_list_by_scale_power(2)
    n2 = sum([J.dim()  for J in jordan_list  if J.is_even()])
    n11 = sum([1  for i in range(len(jordan_list) - 1)  if jordan_list[i].is_odd() and jordan_list[i+1].is_odd()])

    return ZZ(2)**(n11 - n2)



def conway_p_mass(self, p):
    """
    Computes Conway's `p`-mass.

    INPUT:

        `p` -- a prime number > 0

    OUTPUT:

        a rational number > 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, range(1, 6))
        sage: Q.conway_p_mass(2)
        16/3
        sage: Q.conway_p_mass(3)
        729/256

    """
    ## Compute the first two factors of the p-mass
    p_mass = self.conway_diagonal_factor(p) * (p ** (self.conway_cross_product_doubled_power(p) / ZZ(2)))

    ## Multiply by the 'type factor' when p = 2
    if p == 2:
        p_mass *= self.conway_type_factor()

    ## Return the result
    return p_mass



def conway_standard_p_mass(self, p):
    """
    Computes the standard (generic) Conway-Sloane `p`-mass.

    INPUT:

        `p` -- a prime number > 0

    OUTPUT:

        a rational number > 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.conway_standard_p_mass(2)
        2/3

    """
    ## Some useful variables
    n = self.dim()
    if n % 2 == 0:
        s = n // 2
    else:
        s = (n+1) // 2

    ## Compute the inverse of the generic p-mass
    p_mass_inv = 2 * prod([1-p**(-i)  for i in range(2, 2*s, 2)])
    if n % 2 == 0:
        D = (-1)**s * self.det() * (2**n)   ##   We should have something like  D = (-1)**s * self.det() / (2**n), but that's not an integer and here we only care about the square-class.
        #d = self.det()   ## Note: No normalizing power of 2 is needed since the power is even.
        #if not ((p == 2) or (d % p == 0)):
        p_mass_inv *= (1 - kronecker_symbol(fundamental_discriminant(D), p) * p**(-s))

    ## Return the standard p-mass
    return ZZ(1) / p_mass_inv



def conway_standard_mass(self):
    """
    Returns the infinite product of the standard mass factors.

    INPUT:

        none

    OUTPUT:

        a rational number > 0

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [2, -2, 0, 3, -5, 4])
        sage: Q.conway_standard_mass()
        1/6

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.conway_standard_mass()
        1/6

    """
    n = self.dim()
    if n % 2 == 0:
        s = n // 2
    else:
        s = (n+1) // 2

    ## DIAGNOSTIC
    #print "n = ", n
    #print "s = ", s
    #print "Gamma Factor = \n", prod([gamma__exact(j / ZZ(2))  for j in range(1, n+1)])
    #print "Zeta Factor = \n", prod([zeta__exact(2*k)  for k in range(1, s)])
    #print "Pi Factor = \n", pi**((-1) * n * (n+1) / ZZ(4))

    generic_mass = 2 * pi**((-1) * n * (n+1) / ZZ(4)) \
            * prod([gamma__exact(j / ZZ(2))  for j in range(1, n+1)]) \
            * prod([zeta__exact(2*k)  for k in range(1, s)])

    if n % 2 == 0:
        D = (-1)**s * self.det() * (2**n)   ##   We should have something like  D = (-1)**s * self.det() / (2**n), but that's not an integer and here we only care about the square-class.
        generic_mass *= quadratic_L_function__exact(s, D)

    return generic_mass



def conway_mass(self):
    """
    Compute the mass by using the Conway-Sloane mass formula.

    OUTPUT:

    a rational number > 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.conway_mass()
        1/48

        sage: Q = DiagonalQuadraticForm(ZZ, [7,1,1])
        sage: Q.conway_mass()
        3/16

        sage: Q = QuadraticForm(ZZ, 3, [7, 2, 2, 2, 0, 2]) + DiagonalQuadraticForm(ZZ, [1])
        sage: Q.conway_mass()
        3/32

        sage: Q = QuadraticForm(Matrix(ZZ,2,[2,1,1,2]))
        sage: Q.conway_mass()
        1/12
    """
    ## Try to use the cached result
    try:
        return self.__conway_mass
    except AttributeError:
        # Double the form so it's integer-matrix
        Q = self.scale_by_factor(2)

        # Compute the standard mass
        mass = Q.conway_standard_mass()

        # Adjust the p-masses when p|2d
        d = self.det()
        for p in prime_divisors(2*d):
            mass *= (Q.conway_p_mass(p) / Q.conway_standard_p_mass(p))

        # Cache and return the (simplified) result
        self.__conway_mass = QQ(mass.canonicalize_radical()).abs()
        return self.__conway_mass


## ========================================================




#def conway_generic_mass(self):
#    """
#    Computes the generic mass given as
#         2 \pi^{-n(n+1)/4} \prod_{j=1}^{n} \Gamma\(\tfrac{j}{2}\)
#        \zeta(2) \cdots \zeta(2s-2) \zeta_{D}(s)
#    where $n = 2s$ or $2s-1$ depending on the parity of $n$,
#    and $D = (-1)^{s} d$.  We interpret the symbol $\(\frac{D}{p}\)$
#    as 0 if $p\mid 2d$.
#    (Conway and Sloane, Mass formula paper, p??)
#
#    This is possibly equal to
#        2^{-td} * \tau(G) *[\prod_{i=1}^{t} \zeta(1-2i) ]* L(1-t, \chi)
#    where $\dim(Q) = n = 2t$ or $2t+1$, and the last factor is omitted
#    when $n$ is odd.
#    (GHY, Prop 7.4 and 7.5, p121)
#    """
#    RR = RealField(200)
#    n = self.dim()
#    if n % 2 == 0:
#        s = n / 2
#    else:
#        s = (n-1) / 2
#
#    ## Form the generic zeta product
#    ans = 2 * RR(pi)^(-n * (n+1) / 4)
#    for j in range(1,n+1):
#        ans *= gamma(RR(j/2))
#    for j in range(2, 2*s, 2):  ## j = 2, ..., 2s-2
#        ans *= zeta(RR(j))
#
#    ## Extra L-factor for even dimensional forms  -- DO THIS!!!
#    raise NotImplementedError, "This routine is not finished yet... =("
#
#    ## Return the answer
#    return ans






#def conway_p_mass_adjustment(self, p):
#    """
#    Computes the adjustment to give the p-mass from the generic mass.
#    """
#    pass


########################################################################
