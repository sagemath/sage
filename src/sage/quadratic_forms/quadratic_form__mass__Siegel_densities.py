"""
Local Masses and Siegel Densities
"""
########################################################################
#  Computes the local masses (rep'n densities of a form by itself) for a quadratic form over ZZ
#      using the papers of Pall [PSPUM VIII (1965), pp95--105]  for p>2, and Watson [Mathematika
#      23, no. 1, (1976), pp 94--106] for p=2.  These formulas will also work for any local field
#          which is unramified at p=2.
##
#  Copyright by Jonathan Hanke 2007 <jonhanke@gmail.com>
########################################################################

from copy import deepcopy

from sage.arith.misc import kronecker, legendre_symbol, prime_divisors
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.functional import squarefree_part
from sage.misc.misc_c import prod
from sage.misc.mrange import mrange
from sage.quadratic_forms.special_values import gamma__exact, zeta__exact, quadratic_L_function__exact
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


def mass__by_Siegel_densities(self, odd_algorithm='Pall', even_algorithm='Watson'):
    """
    Return the mass of transformations (det 1 and -1).

    .. WARNING::

        This is broken right now...

    INPUT:

    - ``odd_algorithm`` -- algorithm to be used when `p>2`; ``'Pall'`` (only
      one choice for now)
    - ``even_algorithm`` -- algorithm to be used when `p=2`; either
      ``'Kitaoka'`` or ``'Watson'`` (the default)

    REFERENCES:

    - Nipp's Book "Tables of Quaternary Quadratic Forms".
    - Papers of Pall (only for `p>2`) and Watson (for `p=2` -- tricky!).
    - Siegel, Milnor-Hussemoller, Conway-Sloane Paper IV, Kitoaka (all of which
      have problems...)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: m = Q.mass__by_Siegel_densities(); m                                      # needs sage.symbolic
        1/384
        sage: m - (2^Q.dim() * factorial(Q.dim()))^(-1)                                 # needs sage.symbolic
        0

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: m = Q.mass__by_Siegel_densities(); m                                      # needs sage.symbolic
        1/48
        sage: m - (2^Q.dim() * factorial(Q.dim()))^(-1)                                 # needs sage.symbolic
        0
    """
    from sage.symbolic.constants import pi

    # Setup
    n = self.dim()
    s = (n - 1) // 2
    if n % 2 != 0:
        char_d = squarefree_part(2 * self.det())   # Accounts for the det as a QF
    else:
        char_d = squarefree_part(self.det())

    # Form the generic zeta product
    generic_prod = ZZ(2) * (pi)**(-ZZ(n) * (n + 1) / 4)
    ##########################################
    generic_prod *= self.det()**(ZZ(n + 1) / 2)  # ***** This uses the Hessian Determinant ********
    ##########################################
    generic_prod *= prod([gamma__exact(ZZ(j) / 2) for j in range(1, n + 1)])
    generic_prod *= prod([zeta__exact(ZZ(j)) for j in range(2, 2 * s + 1, 2)])
    if n % 2 == 0:
        generic_prod *= quadratic_L_function__exact(n // 2, ZZ(-1)**(n // 2) * char_d)
    # Determine the adjustment factors
    adj_prod = ZZ.one()
    for p in prime_divisors(2 * self.det()):
        # Cancel out the generic factors
        p_adjustment = prod([1 - ZZ(p)**(-j) for j in range(2, 2 * s + 1, 2)])
        if n % 2 == 0:
            p_adjustment *= (1 - kronecker((-1)**(n // 2) * char_d, p) * ZZ(p)**(-n // 2))
        # Insert the new mass factors
        if p == 2:
            if even_algorithm == "Kitaoka":
                p_adjustment = p_adjustment / self.Kitaoka_mass_at_2()
            elif even_algorithm == "Watson":
                p_adjustment = p_adjustment / self.Watson_mass_at_2()
            else:
                raise TypeError("There is a problem -- your even_algorithm argument is invalid.  Try again. =(")
        else:
            if odd_algorithm == "Pall":
                p_adjustment = p_adjustment / self.Pall_mass_density_at_odd_prime(p)
            else:
                raise TypeError("There is a problem -- your optional arguments are invalid.  Try again. =(")

        # Put them together (cumulatively)
        adj_prod *= p_adjustment

        # Extra adjustment for the case of a 2-dimensional form.
    # if (n == 2):
    #    generic_prod *= 2

    # Return the mass
    mass = generic_prod * adj_prod
    return mass


def Pall_mass_density_at_odd_prime(self, p):
    r"""
    Return the local representation density of a form (for
    representing itself) defined over `\ZZ`, at some prime `p>2`.

    REFERENCES:

    Pall's article "The Weight of a Genus of Positive n-ary Quadratic Forms"
    appearing in Proc. Symp. Pure Math. VIII (1965), pp95--105.

    INPUT:

    - ``p`` -- a prime number > 2

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1,0,0,1,0,1])
        sage: Q.Pall_mass_density_at_odd_prime(3)
        [(0, Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 0 0 ]
        [ * 1 0 ]
        [ * * 1 ])] [(0, 3, 8)] [8/9] 8/9
        8/9
    """
    # Check that p is a positive prime -- unnecessary since it's done implicitly in the next step. =)
    if p <= 2:
        raise TypeError("we need p to be a prime > 2")

    # Step 1: Obtain a p-adic (diagonal) local normal form, and
    # compute the invariants for each Jordan block.
    jordan_list = self.jordan_blocks_by_scale_and_unimodular(p)
    modified_jordan_list = [(a, Q.dim(), Q.det()) for a, Q in jordan_list]     # List of pairs (scale, det)

    # Step 2: Compute the list of local masses for each Jordan block
    jordan_mass_list = []
    for (s, n, d) in modified_jordan_list:
        generic_factor = prod([1 - p**(-2 * j) for j in range(1, (n - 1) // 2 + 1)])
        if n % 2 == 0:
            m = n // 2
            generic_factor *= (1 + legendre_symbol(((-1)**m) * d, p) * p**(-m))
        jordan_mass_list = jordan_mass_list + [generic_factor]

    # Step 3: Compute the local mass $\al_p$ at p.
        MJL = modified_jordan_list
    s = len(modified_jordan_list)
    M = [sum([MJL[j][1] for j in range(i, s)]) for i in range(s - 1)]    # Note: It's s-1 since we don't need the last M.
    nu = sum([M[i] * MJL[i][0] * MJL[i][1] for i in range(s - 1)]) - ZZ(sum([J[0] * J[1] * (J[1] - 1) for J in MJL])) / ZZ(2)
    p_mass = prod(jordan_mass_list)
    p_mass *= 2**(s - 1) * p**nu

    print(jordan_list, MJL, jordan_mass_list, p_mass)

    # Return the result
    return p_mass


def Watson_mass_at_2(self):
    """
    Return the local mass of the quadratic form when `p=2`, according
    to Watson's Theorem 1 of "The 2-adic density of a quadratic form"
    in Mathematika 23 (1976), pp 94--106.

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.Watson_mass_at_2()  # WARNING:  WE NEED TO CHECK THIS CAREFULLY!        # needs sage.symbolic
        384
    """
    from sage.functions.all import sgn

    # Make a 0-dim'l quadratic form (for initialization purposes)
    Null_Form = deepcopy(self)
    Null_Form.__init__(ZZ, 0)

    # Step 0: Compute Jordan blocks and bounds of the scales to keep track of
    Jordan_Blocks = self.jordan_blocks_by_scale_and_unimodular(2)
    scale_list = [B[0] for B in Jordan_Blocks]
    s_min = min(scale_list)
    s_max = max(scale_list)

    # Step 1: Compute dictionaries of the diagonal block and 2x2 block for each scale
    diag_dict = {i: Null_Form for i in range(s_min - 2, s_max + 4)}     # Initialize with the zero form
    dim2_dict = {i: Null_Form for i in range(s_min, s_max + 4)}       # Initialize with the zero form
    for s, L in Jordan_Blocks:
        i = 0
        while i < L.dim() - 1 and L[i, i + 1] == 0:      # Find where the 2x2 blocks start
            i += 1
        if i < L.dim() - 1:
            diag_dict[s] = L.extract_variables(range(i))                # Diagonal Form
            dim2_dict[s + 1] = L.extract_variables(range(i, L.dim()))     # Non-diagonal Form
        else:
            diag_dict[s] = L

    # Step 2: Compute three dictionaries of invariants (for n_j, m_j, nu_j)
    n_dict = {j: 0 for j in range(s_min + 1, s_max + 2)}
    m_dict = {j: 0 for j in range(s_min, s_max + 4)}
    for s, L in Jordan_Blocks:
        n_dict[s + 1] = L.dim()
        if diag_dict[s].dim() == 0:
            m_dict[s + 1] = ZZ.one() / ZZ(2) * L.dim()
        else:
            m_dict[s + 1] = ZZ(L.dim() - 1) // ZZ(2)

    nu_dict = {j: n_dict[j + 1] - 2 * m_dict[j + 1]
               for j in range(s_min, s_max + 1)}
    nu_dict[s_max + 1] = 0

    # Step 3: Compute the e_j dictionary
    eps_dict = {}
    for j in range(s_min, s_max + 3):
        two_form = (diag_dict[j - 2] + diag_dict[j] + dim2_dict[j]).scale_by_factor(2)
        j_form = (two_form + diag_dict[j - 1]).change_ring(IntegerModRing(4))

        if j_form.dim() == 0:
            eps_dict[j] = 1
        else:
            iter_vec = [4] * j_form.dim()
            alpha = sum(1 for x in mrange(iter_vec) if j_form(x) == 0)
            beta = sum(1 for x in mrange(iter_vec) if j_form(x) == 2)
            if alpha > beta:
                eps_dict[j] = 1
            elif alpha == beta:
                eps_dict[j] = 0
            else:
                eps_dict[j] = -1

    # Step 4: Compute the quantities nu, q, P, E for the local mass at 2
    nu = sum([j * n_dict[j] * (ZZ(n_dict[j] + 1) / ZZ(2) +
                               sum([n_dict[r] for r in range(j + 1, s_max + 2)])) for j in range(s_min + 1, s_max + 2)])
    q = sum([sgn(nu_dict[j - 1] * (n_dict[j] + sgn(nu_dict[j]))) for j in range(s_min + 1, s_max + 2)])
    P = prod([prod([1 - QQ(4)**(-k) for k in range(1, m_dict[j] + 1)]) for j in range(s_min + 1, s_max + 2)])
    E = prod([ZZ(1) / ZZ(2) * (1 + eps_dict[j] * QQ(2)**(-m_dict[j])) for j in range(s_min, s_max + 3)])

    # Step 5: Compute the local mass for the prime 2.
    mass_at_2 = QQ(2)**(nu - q) * P / E
    return mass_at_2


def Kitaoka_mass_at_2(self):
    """
    Return the local mass of the quadratic form when `p=2`, according
    to Theorem 5.6.3 on pp108--9 of Kitaoka's Book "The Arithmetic of
    Quadratic Forms".

    OUTPUT: a rational number > 0

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.Kitaoka_mass_at_2()   # WARNING:  WE NEED TO CHECK THIS CAREFULLY!
        1/2
    """
    # Make a 0-dim'l quadratic form (for initialization purposes)
    Null_Form = deepcopy(self)
    Null_Form.__init__(ZZ, 0)

    # Step 0: Compute Jordan blocks and bounds of the scales to keep track of
    Jordan_Blocks = self.jordan_blocks_by_scale_and_unimodular(2)
    scale_list = [B[0] for B in Jordan_Blocks]
    s_min = min(scale_list)
    s_max = max(scale_list)

    # Step 1: Compute dictionaries of the diagonal block and 2x2 block for each scale
    diag_dict = {i: Null_Form for i in range(s_min - 2, s_max + 4)}   # Initialize with the zero form
    dim2_dict = {i: Null_Form for i in range(s_min, s_max + 4)}       # Initialize with the zero form
    for s, L in Jordan_Blocks:
        i = 0
        while i < L.dim() - 1 and L[i, i + 1] == 0:      # Find where the 2x2 blocks start
            i += 1
        if i < L.dim() - 1:
            diag_dict[s] = L.extract_variables(range(i))                # Diagonal Form
            dim2_dict[s + 1] = L.extract_variables(range(i, L.dim()))   # Non-diagonal Form
        else:
            diag_dict[s] = L

    # ################  START EDITING HERE  ##################

    # Compute q := sum of the q_j
    q = 0
    for j in range(s_min, s_max + 1):
        if diag_dict[j].dim() > 0:               # Check that N_j is odd (i.e. rep'ns an odd #)
            if diag_dict[j + 1].dim() == 0:
                q += Jordan_Blocks[j][1].dim()        # When N_{j+1} is "even", add n_j
            else:
                q += Jordan_Blocks[j][1].dim() + 1    # When N_{j+1} is "odd", add n_j + 1

    # Compute P = product of the P_j
    P = QQ.one()
    for j in range(s_min, s_max + 1):
        tmp_m = dim2_dict[j].dim() // 2
        P *= prod(QQ.one() - QQ(4**(-k)) for k in range(1, tmp_m + 1))

    # Compute the product E := prod_j (1 / E_j)
    E = QQ.one()
    for j in range(s_min - 1, s_max + 2):
        if (diag_dict[j - 1].dim() == 0) and (diag_dict[j + 1].dim() == 0) and \
           ((diag_dict[j].dim() != 2) or (((diag_dict[j][0, 0] - diag_dict[j][1, 1]) % 4) != 0)):

            # Deal with the complicated case:
            tmp_m = dim2_dict[j].dim() // 2
            if dim2_dict[j].is_hyperbolic(2):
                E *= QQ(2) / (1 + 2**(-tmp_m))
            else:
                E *= QQ(2) / (1 - 2**(-tmp_m))

        else:
            E *= 2

    # Compute the exponent w
    w = QQ.zero()
    for j in range(s_min, s_max + 1):
        n_j = Jordan_Blocks[j][1].dim()
        for k in range(j + 1, s_max + 1):
            n_k = Jordan_Blocks[k][1].dim()
            w += j * n_j * (n_k + QQ(n_j + 1) / 2)

    # Step 5: Compute the local mass for the prime 2.
    mass_at_2 = (QQ(2)**(w - q)) * P * E
    return mass_at_2


def mass_at_two_by_counting_mod_power(self, k):
    r"""
    Compute the local mass at `p=2` assuming that it's stable (mod `2^k`).

    .. NOTE::

        This is **way** too slow to be useful, even when `k=1`.

    .. TODO::

        Remove this routine, or try to compile it!

    INPUT:

    - ``k`` -- integer `\geq 1`

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.mass_at_two_by_counting_mod_power(1)
        4
    """
    R = IntegerModRing(2**k)
    Q1 = self.change_ring(R)
    n = self.dim()
    MS = MatrixSpace(R, n)

    ct = sum(1 for x in mrange([2**k] * (n**2)) if Q1(MS(x)) == Q1)   # Count the solutions mod 2^k
    return ZZ.one() / 2 * (ZZ(ct) / ZZ(2)**(k * n * (n - 1) / 2))
