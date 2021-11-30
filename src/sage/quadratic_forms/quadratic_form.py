"""
Quadratic Forms Overview

AUTHORS:

- Jon Hanke (2007-06-19)
- Anna Haensch (2010-07-01): Formatting and ReSTification
- Simon Brandhorst (2019-10-15): :meth:`quadratic_form_from_invariants`
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein and Jonathan Hanke
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from warnings import warn
from copy import deepcopy

from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.element import is_Matrix
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.ring import Ring
from sage.misc.functional import denominator, is_even
from sage.arith.all import GCD, LCM
from sage.rings.all import Ideal, QQ
from sage.rings.ring import is_Ring, PrincipalIdealDomain
from sage.structure.sage_object import SageObject
from sage.structure.element import is_Vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module_element import vector
from sage.quadratic_forms.genera.genus import genera
from sage.quadratic_forms.quadratic_form__evaluate import QFEvaluateVector, QFEvaluateMatrix


def QuadraticForm__constructor(R, n=None, entries=None):
    """
    Wrapper for the QuadraticForm class constructor.  This is meant
    for internal use within the QuadraticForm class code only.  You
    should instead directly call QuadraticForm().

    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form import QuadraticForm__constructor
        sage: QuadraticForm__constructor(ZZ, 3)   # Makes a generic quadratic form over the integers
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 0 0 0 ]
        [ * 0 0 ]
        [ * * 0 ]

    """
    return QuadraticForm(R, n, entries)


def is_QuadraticForm(Q):
    """
    Determine if the object Q is an element of the QuadraticForm class.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
        sage: from sage.quadratic_forms.quadratic_form import is_QuadraticForm
        sage: is_QuadraticForm(Q)  ##random
        True
        sage: is_QuadraticForm(2)  ##random
        False
    """
    return isinstance(Q, QuadraticForm)


def quadratic_form_from_invariants(F, rk, det, P, sminus):
    r"""
    Return a rational quadratic form with given invariants.

    INPUT:

    - ``F`` -- the base field; currently only ``QQ`` is allowed
    - ``rk`` -- integer; the rank
    - ``det`` -- rational; the determinant
    - ``P`` -- a list of primes where Cassel's Hasse invariant
      is negative
    - ``sminus`` -- integer; the number of negative eigenvalues
      of any Gram matrix

    OUTPUT:

    - a quadratic form with the specified invariants

    Let `(a_1, \ldots, a_n)` be the gram marix of a regular quadratic space.
    Then Cassel's Hasse invariant is defined as

    .. MATH::

        \prod_{i<j} (a_i,a_j),

    where `(a_i,a_j)` denotes the Hilbert symbol.

    ALGORITHM:

    We follow [Kir2016]_.

    EXAMPLES::

        sage: P = [3,5]
        sage: q = quadratic_form_from_invariants(QQ,2,-15,P,1)
        sage: q
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 5 0 ]
        [ * -3 ]
        sage: all(q.hasse_invariant(p)==-1 for p in P)
        True

    TESTS:

    This shows that :trac:`28955` is fixed::

        sage: quadratic_form_from_invariants(QQ,3,2,[2],2)
        Quadratic form in 3 variables over Rational Field with coefficients:
        [ -1 0 0 ]
        [ * 1 0 ]
        [ * * -2 ]

        sage: quadratic_form_from_invariants(QQ,4,2,[2],4)
        Traceback (most recent call last):
        ...
        ValueError: invariants do not define a rational quadratic form
    """
    from sage.arith.misc import hilbert_symbol
    # normalize input
    if F!=QQ:
        raise NotImplementedError('base field must be QQ. If you want this over any field, implement weak approximation.')
    P = [ZZ(p) for p in P]
    rk = ZZ(rk)
    d = QQ(det).squarefree_part()
    sminus = ZZ(sminus)
    # check if the invariants define a global quadratic form
    if d.sign() != (-1)**sminus:
        raise ValueError("invariants do not define a rational quadratic form")
    if rk == 1 and len(P) != 0:
        raise ValueError("invariants do not define a rational quadratic form")
    if rk == 2:
        for p in P:
            if QQ(-d).is_padic_square(p):
                raise ValueError("invariants do not define a rational quadratic form")
    f = 0
    if sminus % 4 in (2, 3):
        f = 1
    if (f + len(P)) % 2 == 1:
        raise ValueError("invariants do not define a rational quadratic form")
    D = []
    while rk >= 2:
        if rk >= 4:
            if sminus > 0:
                a = ZZ(-1)
            else:
                a = ZZ(1)
        elif rk == 3:
            Pprime = [p for p in P if hilbert_symbol(-1, -d, p)==1]
            Pprime += [p for p in (2*d).prime_divisors()
                       if hilbert_symbol(-1, -d, p)==-1 and p not in P]
            if sminus > 0:
                a = ZZ(-1)
            else:
                a = ZZ(1)
            for p in Pprime:
                if d.valuation(p) % 2 == 0:
                    a *= p
            assert all((a*d).valuation(p)%2==1 for p in Pprime)
        elif rk == 2:
            S = P
            if sminus == 2:
                S += [-1]
            a = QQ.hilbert_symbol_negative_at_S(S,-d)
            a = ZZ(a)
        P = ([p for p in P if hilbert_symbol(a, -d, p) == 1]
            +[p for p in (2*a*d).prime_divisors()
              if hilbert_symbol(a, -d, p)==-1 and p not in P])
        sminus = max(0, sminus-1)
        rk = rk - 1
        d = a*d
        D.append(a.squarefree_part())
    d = d.squarefree_part()
    D.append(d)
    return DiagonalQuadraticForm(QQ,D)


class QuadraticForm(SageObject):
    r"""
    The ``QuadraticForm`` class represents a quadratic form in n variables with
    coefficients in the ring R.

    INPUT:

    The constructor may be called in any of the following ways.

    #. ``QuadraticForm(R, n, entries)``, where

       - `R` -- ring for which the quadratic form is defined
       - `n` -- an integer >= 0
       - ``entries`` -- a list of `n(n+1)/2` coefficients of the quadratic form
         in `R` (given lexicographically, or equivalently, by rows of the
         matrix)

    #. ``QuadraticForm(R, n)``, where

       - `R` -- a ring
       - `n` -- a symmetric `n \times n` matrix with even diagonal (relative to
         `R`)

    #. ``QuadraticForm(R)``, where

       - `R` -- a symmetric `n \times n` matrix with even diagonal (relative to
         its base ring)

    If the keyword argument ``unsafe_initialize`` is True, then the subsequent
    fields may by used to force the external initialization of various fields
    of the quadratic form. Currently the only fields which can be set are:

    - ``number_of_automorphisms``
    - ``determinant``


    OUTPUT:

    quadratic form

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
        sage: Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]

    ::

        sage: Q = QuadraticForm(QQ, 3, [1,2,3,4/3 ,5,6])
        sage: Q
        Quadratic form in 3 variables over Rational Field with coefficients:
        [ 1 2 3 ]
        [ * 4/3 5 ]
        [ * * 6 ]
        sage: Q[0,0]
        1
        sage: Q[0,0].parent()
        Rational Field

    ::

        sage: Q = QuadraticForm(QQ, 7, range(28))
        sage: Q
        Quadratic form in 7 variables over Rational Field with coefficients:
        [ 0 1 2 3 4 5 6 ]
        [ * 7 8 9 10 11 12 ]
        [ * * 13 14 15 16 17 ]
        [ * * * 18 19 20 21 ]
        [ * * * * 22 23 24 ]
        [ * * * * * 25 26 ]
        [ * * * * * * 27 ]

    ::

        sage: Q = QuadraticForm(QQ, 2, range(1,4))
        sage: A = Matrix(ZZ,2,2,[-1,0,0,1])
        sage: Q(A)
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1 -2 ]
        [ * 3 ]

    ::

        sage: m = matrix(2,2,[1,2,3,4])
        sage: m + m.transpose()
        [2 5]
        [5 8]
        sage: QuadraticForm(m + m.transpose())
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 1 5 ]
        [ * 4 ]

    ::

        sage: QuadraticForm(ZZ, m + m.transpose())
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 1 5 ]
        [ * 4 ]

    ::

        sage: QuadraticForm(QQ, m + m.transpose())
        Quadratic form in 2 variables over Rational Field with coefficients:
        [ 1 5 ]
        [ * 4 ]
    """

    ## Import specialized methods:
    ## ---------------------------

    ## Routines to compute the p-adic local normal form
    from sage.quadratic_forms.quadratic_form__local_normal_form import \
            find_entry_with_minimal_scale_at_prime, \
            local_normal_form, \
            jordan_blocks_by_scale_and_unimodular, \
            jordan_blocks_in_unimodular_list_by_scale_power

    ## Routines to perform elementary variable substitutions
    from sage.quadratic_forms.quadratic_form__variable_substitutions import \
            swap_variables, \
            multiply_variable, \
            divide_variable, \
            scale_by_factor, \
            extract_variables, \
            elementary_substitution, \
            add_symmetric

    ## Routines to compute p-adic field invariants
    from sage.quadratic_forms.quadratic_form__local_field_invariants import \
            rational_diagonal_form, \
            _rational_diagonal_form_and_transformation, \
            signature_vector, \
            signature, \
            hasse_invariant, \
            hasse_invariant__OMeara, \
            is_hyperbolic, \
            is_anisotropic, \
            is_isotropic, \
            anisotropic_primes, \
            compute_definiteness, \
            compute_definiteness_string_by_determinants, \
            is_positive_definite, \
            is_negative_definite, \
            is_indefinite, \
            is_definite

    ## Routines to compute local densities by the reduction procedure
    from sage.quadratic_forms.quadratic_form__local_density_congruence import \
            count_modp_solutions__by_Gauss_sum, \
            local_good_density_congruence_odd, \
            local_good_density_congruence_even, \
            local_good_density_congruence, \
            local_zero_density_congruence, \
            local_badI_density_congruence, \
            local_badII_density_congruence, \
            local_bad_density_congruence, \
            local_density_congruence, \
            local_primitive_density_congruence

    ## Routines to compute local densities by counting solutions of various types
    from sage.quadratic_forms.quadratic_form__count_local_2 import \
            count_congruence_solutions_as_vector, \
            count_congruence_solutions, \
            count_congruence_solutions__good_type, \
            count_congruence_solutions__zero_type, \
            count_congruence_solutions__bad_type, \
            count_congruence_solutions__bad_type_I, \
            count_congruence_solutions__bad_type_II

    ## Routines to be called by the user to compute local densities
    from sage.quadratic_forms.quadratic_form__local_density_interfaces import \
            local_density, \
            local_primitive_density

    ## Routines for computing with ternary forms
    from sage.quadratic_forms.quadratic_form__ternary_Tornaria import \
            disc, \
            content, \
            adjoint, \
            antiadjoint, \
            is_adjoint, \
            reciprocal, \
            omega, \
            delta, \
            level__Tornaria, \
            discrec, \
            hasse_conductor, \
            clifford_invariant, \
            clifford_conductor, \
            basiclemma, \
            basiclemmavec, \
            xi, \
            xi_rec, \
            lll, \
            representation_number_list, \
            representation_vector_list, \
            is_zero, \
            is_zero_nonsingular, \
            is_zero_singular

    ## Routines to compute the theta function
    from sage.quadratic_forms.quadratic_form__theta import \
            theta_series, \
            theta_series_degree_2, \
            theta_by_pari, \
            theta_by_cholesky

    ## Routines to compute the product of all local densities
    from sage.quadratic_forms.quadratic_form__siegel_product import \
            siegel_product

    ## Routines to compute p-neighbors
    from sage.quadratic_forms.quadratic_form__neighbors import \
            find_primitive_p_divisible_vector__random, \
            find_primitive_p_divisible_vector__next, \
            find_p_neighbor_from_vec, \
            neighbor_iteration, \
            orbits_lines_mod_p

    ## Routines to reduce a given quadratic form
    from sage.quadratic_forms.quadratic_form__reduction_theory import \
            reduced_binary_form1, \
            reduced_ternary_form__Dickson, \
            reduced_binary_form, \
            minkowski_reduction, \
            minkowski_reduction_for_4vars__SP
    ## Wrappers for Conway-Sloane genus routines (in ./genera/)
    from sage.quadratic_forms.quadratic_form__genus import \
            global_genus_symbol, \
            local_genus_symbol, \
            CS_genus_symbol_list


    ## Routines to compute local masses for ZZ.
    from sage.quadratic_forms.quadratic_form__mass import \
            shimura_mass__maximal, \
            GHY_mass__maximal
    from sage.quadratic_forms.quadratic_form__mass__Siegel_densities import \
            mass__by_Siegel_densities, \
            Pall_mass_density_at_odd_prime, \
            Watson_mass_at_2, \
            Kitaoka_mass_at_2, \
            mass_at_two_by_counting_mod_power
    from sage.quadratic_forms.quadratic_form__mass__Conway_Sloane_masses import \
            parity, \
            is_even, \
            is_odd, \
            conway_species_list_at_odd_prime, \
            conway_species_list_at_2, \
            conway_octane_of_this_unimodular_Jordan_block_at_2, \
            conway_diagonal_factor, \
            conway_cross_product_doubled_power, \
            conway_type_factor, \
            conway_p_mass, \
            conway_standard_p_mass, \
            conway_standard_mass, \
            conway_mass
#            conway_generic_mass, \
#            conway_p_mass_adjustment

    ## Routines to check local representability of numbers
    from sage.quadratic_forms.quadratic_form__local_representation_conditions import \
            local_representation_conditions, \
            is_locally_universal_at_prime, \
            is_locally_universal_at_all_primes, \
            is_locally_universal_at_all_places, \
            is_locally_represented_number_at_place, \
            is_locally_represented_number

    ## Routines to make a split local covering of the given quadratic form.
    from sage.quadratic_forms.quadratic_form__split_local_covering import \
            cholesky_decomposition, \
            vectors_by_length, \
            complementary_subform_to_vector, \
            split_local_cover

    ## Routines to make automorphisms of the given quadratic form.
    from sage.quadratic_forms.quadratic_form__automorphisms import \
            basis_of_short_vectors, \
            short_vector_list_up_to_length, \
            short_primitive_vector_list_up_to_length, \
            _compute_automorphisms, \
            automorphism_group, \
            automorphisms, \
            number_of_automorphisms, \
            set_number_of_automorphisms

    ## Routines to test the local and global equivalence/isometry of two quadratic forms.
    from sage.quadratic_forms.quadratic_form__equivalence_testing import \
            is_globally_equivalent_to, \
            is_locally_equivalent_to, \
            has_equivalent_Jordan_decomposition_at_prime, \
            is_rationally_isometric

    ## Routines for solving equations of the form Q(x) = c.
    from sage.quadratic_forms.qfsolve import solve



    def __init__(self, R, n=None, entries=None, unsafe_initialization=False, number_of_automorphisms=None, determinant=None):
        """
        EXAMPLES::

            sage: s = QuadraticForm(ZZ, 4, range(10))
            sage: s.dim()
            4

        TESTS::

            sage: s == loads(dumps(s))
            True
            sage: QuadraticForm(ZZ, -1)
            Traceback (most recent call last):
            ...
            ValueError: the size must be a non-negative integer, not -1

            sage: x = polygen(ZZ, 'x')
            sage: QuadraticForm(x**2)
            Traceback (most recent call last):
            ....
            TypeError: wrong input for QuadraticForm
        """
        # Deal with:  QuadraticForm(ring, matrix)
        matrix_init_flag = False
        if isinstance(R, Ring):
            if is_Matrix(n):
                # Test if n is symmetric and has even diagonal
                if not self._is_even_symmetric_matrix_(n, R):
                    raise TypeError("Oops!  The matrix is not a symmetric with even diagonal defined over R.")

                # Rename the matrix and ring
                M = n
                M_ring = R
                matrix_init_flag = True

        elif not is_Matrix(R):
            # first argument, if not a ring, must be a matrix
            raise TypeError('wrong input for QuadraticForm')
        else:
            # Deal with:  QuadraticForm(matrix)
            # Test if R is symmetric and has even diagonal
            if not self._is_even_symmetric_matrix_(R):
                raise TypeError("Oops!  The matrix is not a symmetric with even diagonal.")

            # Rename the matrix and ring
            M = R
            M_ring = R.base_ring()
            matrix_init_flag = True

        ## Perform the quadratic form initialization
        if matrix_init_flag:
            self.__n = ZZ(M.nrows())
            self.__base_ring = M_ring
            self.__coeffs = []
            for i in range(M.nrows()):
                for j in range(i, M.nrows()):
                    if (i == j):
                        self.__coeffs += [ M_ring(M[i,j] / 2) ]
                    else:
                        self.__coeffs += [ M_ring(M[i,j]) ]

            return

        ## -----------------------------------------------------------

        ## Verify the size of the matrix is an integer >= 0
        n = ZZ(n)
        if n < 0:
            raise ValueError("the size must be a non-negative integer, not {}".format(n))

        # Store the relevant variables
        N = n * (n + 1) // 2
        self.__n = n
        self.__base_ring = R
        self.__coeffs = [self.__base_ring.zero() for i in range(N)]

        # Check if entries is a list, tuple or iterator for the
        # current size, and if so, write the upper-triangular matrix
        if entries is not None:
            try:
                entries = list(entries)
            except TypeError:
                raise TypeError('entries must be an iterable')

            if len(entries) == N:
                for i in range(N):
                    self.__coeffs[i] = self.__base_ring(entries[i])
            else:
                raise TypeError("Oops! The entries " + str(entries) + " must be a list of size n(n+1)/2.")

        ## -----------------------------------------------------------

        ## Process possible forced initialization of various fields
        self._external_initialization_list = []
        if unsafe_initialization:

            ## Set the number of automorphisms
            if number_of_automorphisms is not None:
                self.set_number_of_automorphisms(number_of_automorphisms)
                #self.__number_of_automorphisms = number_of_automorphisms
                #self.__external_initialization_list.append('number_of_automorphisms')

            ## Set the determinant
            if determinant is not None:
                self.__det = determinant
                self._external_initialization_list.append('determinant')

    def list_external_initializations(self):
        """
        Return a list of the fields which were set externally at
        creation, and not created through the usual QuadraticForm
        methods.  These fields are as good as the external process
        that made them, and are thus not guaranteed to be correct.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,0,5])
            sage: Q.list_external_initializations()
            []
            sage: T = Q.theta_series()
            sage: Q.list_external_initializations()
            []
            sage: Q = QuadraticForm(ZZ, 2, [1,0,5], unsafe_initialization=False, number_of_automorphisms=3, determinant=0)
            sage: Q.list_external_initializations()
            []

        ::

            sage: Q = QuadraticForm(ZZ, 2, [1,0,5], unsafe_initialization=False, number_of_automorphisms=3, determinant=0)
            sage: Q.list_external_initializations()
            []
            sage: Q = QuadraticForm(ZZ, 2, [1,0,5], unsafe_initialization=True, number_of_automorphisms=3, determinant=0)
            sage: Q.list_external_initializations()
            ['number_of_automorphisms', 'determinant']
        """
        return deepcopy(self._external_initialization_list)


    def __pari__(self):
        """
        Return a PARI-formatted Hessian matrix for Q.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,0,5])
            sage: Q.__pari__()
            [2, 0; 0, 10]

        """
        return self.matrix().__pari__()

    def _pari_init_(self):
        """
        Return a PARI-formatted Hessian matrix for Q, as string.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,0,5])
            sage: Q._pari_init_()
            'Mat([2,0;0,10])'
        """
        return self.matrix()._pari_init_()

    def _repr_(self):
        """
        Give a text representation for the quadratic form given as an upper-triangular matrix of coefficients.

        EXAMPLES::

            sage: QuadraticForm(ZZ, 2, [1,3,5])
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 1 3 ]
            [ * 5 ]

        """
        n = self.dim()
        out_str = "Quadratic form in " + str(n) + " variables over " + str(self.base_ring()) + " with coefficients: \n"
        for i in range(n):
            if i > 0:
                out_str += '\n'
            out_str += "[ "
            for j in range(n):
                if (i > j):
                    out_str += "* "
                else:
                    out_str += str(self[i,j]) + " "
            out_str += "]"
        return out_str


    def _latex_(self):
        """
        Give a LaTeX representation for the quadratic form given as an upper-triangular matrix of coefficients.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [2,3,5])
            sage: Q._latex_()
            'Quadratic form in 2 variables over Integer Ring with coefficients: \\newline\\left[ \\begin{array}{cc}2 & 3 &  * & 5 & \\end{array} \\right]'

        """
        n = self.dim()
        out_str = ""
        out_str += "Quadratic form in " + str(n) + " variables over " + str(self.base_ring())
        out_str += " with coefficients: \\newline"
        out_str += "\\left[ \\begin{array}{" + n * "c" + "}"
        for i in range(n):
            for j in range(n):
                if (i > j):
                    out_str += " * & "
                else:
                    out_str += str(self[i,j]) + " & "
#            if i < (n-1):
#                out_str += "\\"
        out_str += "\\end{array} \\right]"
        return out_str

    def __getitem__(self, ij):
        """
        Return the coefficient `a_{ij}` of `x_i * x_j`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: matrix(ZZ, 3, 3, [Q[i,j]  for i in range(3) for j in range(3)])
            [1 2 3]
            [2 4 5]
            [3 5 6]

        """
        ## Unpack the list of indices
        i, j =  ij
        i = int(i)
        j = int(j)

        ## Ensure we're using upper-triangular coordinates
        if i > j:
            tmp = i
            i = j
            j = tmp

        return self.__coeffs[i*self.__n - i*(i-1)//2 + j - i]

    def __setitem__(self, ij, coeff):
        """
        Set the coefficient `a_{ij}` in front of `x_i * x_j`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
            sage: Q
            Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 5 ]
            [ * * 6 ]
            sage: Q[2,1] = 17
            sage: Q
            Quadratic form in 3 variables over Integer Ring with coefficients:
            [ 1 2 3 ]
            [ * 4 17 ]
            [ * * 6 ]

        """
        ## Unpack the list of indices
        i, j =  ij
        i = int(i)
        j = int(j)

        ## TO DO:  Verify that 0 <= i, j <= (n-1)

        ## Ensure we're using upper-triangular coordinates
        if i > j:
            tmp = i
            i = j
            j = tmp

        ## Set the entry
        try:
            self.__coeffs[i*self.__n - i*(i-1)//2 + j -i] = self.__base_ring(coeff)
        except Exception:
            raise RuntimeError("Oops!  This coefficient can't be coerced to an element of the base ring for the quadratic form.")

    def __hash__(self):
        r"""
        TESTS::

            sage: Q1 = QuadraticForm(QQ, 2, [1,1,1])
            sage: Q2 = QuadraticForm(QQ, 2, [1,1,1])
            sage: Q3 = QuadraticForm(QuadraticField(2), 2, [1,1,1])
            sage: hash(Q1) == hash(Q2)
            True
            sage: hash(Q1) == hash(Q3)
            False
        """
        return hash(self.__base_ring) ^ hash(tuple(self.__coeffs))

    def __eq__(self, right):
        """
        Determines if two quadratic forms are equal.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,4,10])
            sage: Q == Q
            True

            sage: Q1 = QuadraticForm(QQ, 2, [1,4,10])
            sage: Q == Q1
            False

            sage: Q2 = QuadraticForm(ZZ, 2, [1,4,-10])
            sage: Q == Q1
            False
            sage: Q == Q2
            False
            sage: Q1 == Q2
            False

        """
        if not isinstance(right, QuadraticForm):
            return False
        return (self.__base_ring == right.__base_ring) and (self.__coeffs == right.__coeffs)


    def __add__(self, right):
          """
          Return the direct sum of two quadratic forms.

          EXAMPLES::

              sage: Q = QuadraticForm(ZZ, 2, [1,4,10])
              sage: Q
              Quadratic form in 2 variables over Integer Ring with coefficients:
              [ 1 4 ]
              [ * 10 ]
              sage: Q2 = QuadraticForm(ZZ, 2, [1,4,-10])
              sage: Q + Q2
              Quadratic form in 4 variables over Integer Ring with coefficients:
              [ 1 4 0 0 ]
              [ * 10 0 0 ]
              [ * * 1 4 ]
              [ * * * -10 ]
          """
          if not isinstance(right, QuadraticForm):
              raise TypeError("Oops!  Can't add these objects since they're not both quadratic forms. =(")
          elif (self.base_ring() != right.base_ring()):
              raise TypeError("Oops!  Can't add these since the quadratic forms don't have the same base rings... =(")
          else:
              Q = QuadraticForm(self.base_ring(), self.dim() + right.dim())
              n = self.dim()
              m = right.dim()

              for i in range(n):
                  for j in range(i,n):
                      Q[i,j] = self[i,j]

              for i in range(m):
                  for j in range(i,m):
                      Q[n+i,n+j] = right[i,j]

              return Q


    def sum_by_coefficients_with(self, right):
          """
          Return the sum (on coefficients) of two quadratic forms of the same size.

          EXAMPLES::

              sage: Q = QuadraticForm(ZZ, 2, [1,4,10])
              sage: Q
              Quadratic form in 2 variables over Integer Ring with coefficients:
              [ 1 4 ]
              [ * 10 ]
              sage: Q+Q
              Quadratic form in 4 variables over Integer Ring with coefficients:
              [ 1 4 0 0 ]
              [ * 10 0 0 ]
              [ * * 1 4 ]
              [ * * * 10 ]

              sage: Q2 = QuadraticForm(ZZ, 2, [1,4,-10])
              sage: Q.sum_by_coefficients_with(Q2)
              Quadratic form in 2 variables over Integer Ring with coefficients:
              [ 2 8 ]
              [ * 0 ]

          """
          if not isinstance(right, QuadraticForm):
              raise TypeError("Oops!  Can't add these objects since they're not both quadratic forms. =(")
          elif (self.__n != right.__n):
              raise TypeError("Oops!  Can't add these since the quadratic forms don't have the same sizes... =(")
          elif (self.__base_ring != right.__base_ring):
              raise TypeError("Oops!  Can't add these since the quadratic forms don't have the same base rings... =(")
          else:
              return QuadraticForm(self.__base_ring, self.__n, [self.__coeffs[i] + right.__coeffs[i]  for i in range(len(self.__coeffs))])


# ========================  CHANGE THIS TO A TENSOR PRODUCT?!?  Even in Characteristic 2?!?  =======================
#    def __mul__(self, right):
#        """
#        Multiply (on the right) the quadratic form Q by an element of the ring that Q is defined over.
#
#        EXAMPLES::
#
#            sage: Q = QuadraticForm(ZZ, 2, [1,4,10])
#            sage: Q*2
#            Quadratic form in 2 variables over Integer Ring with coefficients:
#            [ 2 8 ]
#            [ * 20 ]
#
#            sage: Q+Q == Q*2
#            True
#
#        """
#        try:
#            c = self.base_ring()(right)
#        except Exception:
#            raise TypeError, "Oh no! The multiplier cannot be coerced into the base ring of the quadratic form. =("
#
#        return QuadraticForm(self.base_ring(), self.dim(), [c * self.__coeffs[i]  for i in range(len(self.__coeffs))])
# =========================================================================================================================



    def __call__(self, v):
        """
        Evaluate this quadratic form Q on a vector or matrix of elements
        coercible to the base ring of the quadratic form.  If a vector
        is given then the output will be the ring element Q(`v`), but if a
        matrix is given then the output will be the quadratic form Q'
        which in matrix notation is given by:

        .. MATH::
                Q' = v^t * Q * v.


        EXAMPLES:

        Evaluate a quadratic form at a vector::

            sage: Q = QuadraticForm(QQ, 3, range(6))
            sage: Q
            Quadratic form in 3 variables over Rational Field with coefficients:
            [ 0 1 2 ]
            [ * 3 4 ]
            [ * * 5 ]
            sage: Q([1,2,3])
            89
            sage: Q([1,0,0])
            0
            sage: Q([1,1,1])
            15

        Evaluate a quadratic form using a column matrix::

            sage: Q = QuadraticForm(QQ, 2, range(1,4))
            sage: A = Matrix(ZZ,2,2,[-1,0,0,1])
            sage: Q(A)
            Quadratic form in 2 variables over Rational Field with coefficients:
            [ 1 -2 ]
            [ * 3 ]
            sage: Q([1,0])
            1
            sage: type(Q([1,0]))
            <... 'sage.rings.rational.Rational'>
            sage: Q = QuadraticForm(QQ, 2, range(1,4))
            sage: Q(matrix(2, [1,0]))
            Quadratic form in 1 variables over Rational Field with coefficients:
            [ 1 ]

        Simple 2x2 change of variables::

            sage: Q = QuadraticForm(ZZ, 2, [1,0,1])
            sage: Q
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 1 0 ]
            [ * 1 ]
            sage: M = Matrix(ZZ, 2, 2, [1,1,0,1])
            sage: M
            [1 1]
            [0 1]
            sage: Q(M)
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 1 2 ]
            [ * 2 ]

        Some more tests::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
            sage: Q([1,2,3])
            14
            sage: v = vector([1,2,3])
            sage: Q(v)
            14
            sage: t = tuple([1,2,3])
            sage: Q(v)
            14
            sage: M = Matrix(ZZ, 3, [1,2,3])
            sage: Q(M)
            Quadratic form in 1 variables over Integer Ring with coefficients:
            [ 14 ]

        """
        ## If we are passed a matrix A, return the quadratic form Q(A(x))
        ## (In matrix notation: A^t * Q * A)
        n = self.dim()

        if is_Matrix(v):
            ## Check that v has the correct number of rows
            if v.nrows() != n:
                raise TypeError("the matrix must have {} rows".format(n))

            ## Create the new quadratic form
            m = v.ncols()
            Q2 = QuadraticForm(self.base_ring(), m)
            return QFEvaluateMatrix(self, v, Q2)

        elif (is_Vector(v) or isinstance(v, (list, tuple))):
            ## Check the vector/tuple/list has the correct length
            if not (len(v) == n):
                raise TypeError("your vector needs to have length {}".format(n))

            ## TO DO:  Check that the elements can be coerced into the base ring of Q -- on first elt.
            if len(v) > 0:
                try:
                    self.base_ring()(v[0])
                except Exception:
                    raise TypeError("your vector is not coercible to the base ring of the quadratic form")

            ## Attempt to evaluate Q[v]
            return QFEvaluateVector(self, v)

        else:
            raise TypeError




# =====================================================================================================

    def _is_even_symmetric_matrix_(self, A, R=None):
        """
        Tests if a matrix is symmetric, defined over R, and has even diagonal in R.

        INPUT:
            A -- matrix

            R -- ring

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [2,3,5])
            sage: A = Q.matrix()
            sage: A
            [ 4  3]
            [ 3 10]
            sage: Q._is_even_symmetric_matrix_(A)
            True
            sage: A[0,0] = 1
            sage: Q._is_even_symmetric_matrix_(A)
            False

        """
        if not is_Matrix(A):
            raise TypeError("A is not a matrix.")

        ring_coerce_test = True
        if R is None:            ## This allows us to omit the ring from the variables, and take it from the matrix
            R = A.base_ring()
            ring_coerce_test = False

        if not isinstance(R, Ring):
            raise TypeError("R is not a ring.")

        if not A.is_square():
            return False

        ## Test that the matrix is symmetric
        n = A.nrows()
        for i in range(n):
            for j in range(i+1, n):
                if A[i,j] != A[j,i]:
                    return False

        ## Test that all entries coerce to R
        if not ((A.base_ring() == R) or ring_coerce_test):
            try:
                for i in range(n):
                    for j in range(i, n):
                        R(A[i,j])
            except Exception:
                return False

        ## Test that the diagonal is even (if 1/2 isn't in R)
        if not R(2).is_unit():
            for i in range(n):
                if not is_even(R(A[i,i])):
                    return False

        return True


# =====================================================================================================

    def matrix(self):
        """
        Return the Hessian matrix A for which Q(X) =  `(1/2) * X^t * A * X`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 3, range(6))
            sage: Q.matrix()
            [ 0  1  2]
            [ 1  6  4]
            [ 2  4 10]

        """
        return self.Hessian_matrix()

    def Hessian_matrix(self):
        """
        Return the Hessian matrix A for which Q(X) = `(1/2) * X^t * A * X`.

        EXAMPLES::

            sage: Q = QuadraticForm(QQ, 2, range(1,4))
            sage: Q
            Quadratic form in 2 variables over Rational Field with coefficients:
            [ 1 2 ]
            [ * 3 ]
            sage: Q.Hessian_matrix()
            [2 2]
            [2 6]
            sage: Q.matrix().base_ring()
            Rational Field

        """
        mat_entries = []
        for i in range(self.dim()):
            for j in range(self.dim()):
                if (i == j):
                    mat_entries += [ 2 * self[i,j] ]
                else:
                    mat_entries += [ self[i,j] ]

        return matrix(self.base_ring(), self.dim(), self.dim(), mat_entries)

    def Gram_matrix_rational(self):
        """
        Return a (symmetric) Gram matrix A for the quadratic form Q,
        meaning that

        .. MATH::

            Q(x) = x^t * A * x,

        defined over the fraction field of the base ring.

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
            sage: A = Q.Gram_matrix_rational(); A
            [1 0 0 0]
            [0 3 0 0]
            [0 0 5 0]
            [0 0 0 7]
            sage: A.base_ring()
            Rational Field

        """
        return (ZZ(1) / ZZ(2)) * self.matrix()

    def Gram_matrix(self):
        """
        Return a (symmetric) Gram matrix A for the quadratic form Q,
        meaning that

        .. MATH::

            Q(x) = x^t * A * x,

        defined over the base ring of Q.  If this is not possible,
        then a TypeError is raised.

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
            sage: A = Q.Gram_matrix(); A
            [1 0 0 0]
            [0 3 0 0]
            [0 0 5 0]
            [0 0 0 7]
            sage: A.base_ring()
            Integer Ring

        """
        A = (ZZ(1) / ZZ(2)) * self.matrix()
        n = self.dim()

        ## Test to see if it has an integral Gram matrix
        Int_flag = True
        for i in range(n):
            for j in range(i,n):
                Int_flag = Int_flag and A[i,j] in self.base_ring()

        ## Return the Gram matrix, or an error
        if Int_flag:
            return MatrixSpace(self.base_ring(), n, n)(A)
        else:
            raise TypeError("Oops!  This form does not have an integral Gram matrix. =(")

    def has_integral_Gram_matrix(self):
        """
        Return whether the quadratic form has an integral Gram matrix (with respect to its base ring).

        A warning is issued if the form is defined over a field, since in that case the return is trivially true.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [7,8,9])
            sage: Q.has_integral_Gram_matrix()
            True

        ::

            sage: Q = QuadraticForm(ZZ, 2, [4,5,6])
            sage: Q.has_integral_Gram_matrix()
            False

        """
        # Warning over fields
        if self.base_ring().is_field():
            warn("Warning -- A quadratic form over a field always has integral Gram matrix.  Do you really want to do this?!?")

        # Determine integrality of the Gram matrix
        flag = True
        try:
            self.Gram_matrix()
        except Exception:
            flag = False

        return flag

    def gcd(self):
        """
        Return the greatest common divisor of the coefficients of the
        quadratic form (as a polynomial).

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 4, range(1, 21, 2))
            sage: Q.gcd()
            1

        ::

            sage: Q = QuadraticForm(ZZ, 4, range(0, 20, 2))
            sage: Q.gcd()
            2
        """
        if self.base_ring() != ZZ:
            raise TypeError("Oops! The given quadratic form must be defined over ZZ.")

        return GCD(self.coefficients())


    def polynomial(self,names='x'):
        r"""
        Return the polynomial in 'n' variables of the quadratic form in the ring 'R[names].'

        INPUT:

            -'self' - a quadratic form over a commutative ring.
            -'names' - the name of the variables. Digits will be appended to the name for each different canonical
            variable e.g x1, x2, x3 etc.

        OUTPUT:

            The polynomial form of the quadratic form.

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(QQ,[1, 3, 5, 7])
            sage: P = Q.polynomial(); P
            x0^2 + 3*x1^2 + 5*x2^2 + 7*x3^2

        ::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: Z = F.ring_of_integers()
            sage: Q = QuadraticForm(Z,3,[2*a, 3*a, 0 , 1 - a, 0, 2*a + 4])
            sage: P = Q.polynomial(names='y'); P
            2*a*y0^2 + 3*a*y0*y1 + (-a + 1)*y1^2 + (2*a + 4)*y2^2
            sage: Q = QuadraticForm(F,4,[a, 3*a, 0, 1 - a, a - 3, 0, 2*a + 4, 4 + a, 0, 1])
            sage: Q.polynomial(names='z')
            a*z0^2 + (3*a)*z0*z1 + (a - 3)*z1^2 + (a + 4)*z2^2 + (-a + 1)*z0*z3 + (2*a + 4)*z1*z3 + z3^2
            sage: B.<i,j,k> = QuaternionAlgebra(F,-1,-1)
            sage: Q = QuadraticForm(B, 3, [2*a, 3*a, i, 1 - a, 0, 2*a + 4])
            sage: Q.polynomial()
            Traceback (most recent call last):
            ...
            ValueError: Can only create polynomial rings over commutative rings.
        """
        B = self.base_ring()
        n = self.dim()
        M = matrix(B, n)
        for i in range(n):
            for j in range(i, n):
                M[i,j] = self[i,j]
        try:
            R = PolynomialRing(self.base_ring(), names, n)
        except Exception:
            raise ValueError('Can only create polynomial rings over commutative rings.')
        V = vector(R.gens())
        P = (V*M).dot_product(V)
        return P




    def is_primitive(self):
        """
        Determines if the given integer-valued form is primitive
        (i.e. not an integer (>1) multiple of another integer-valued
        quadratic form).

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [2,3,4])
            sage: Q.is_primitive()
            True
            sage: Q = QuadraticForm(ZZ, 2, [2,4,8])
            sage: Q.is_primitive()
            False

        """
        return (self.gcd() == 1)

    def primitive(self):
        """
        Return a primitive version of an integer-valued quadratic form, defined over `ZZ`.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [2,3,4])
            sage: Q.primitive()
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 2 3 ]
            [ * 4 ]
            sage: Q = QuadraticForm(ZZ, 2, [2,4,8])
            sage: Q.primitive()
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 1 2 ]
            [ * 4 ]

        """
        if self.base_ring() != ZZ:
            raise TypeError("Oops! The given quadratic form must be defined over ZZ.")

        g = self.gcd()
        return QuadraticForm(self.base_ring(), self.dim(), [ZZ(x/g)  for x in self.coefficients()])

    def adjoint_primitive(self):
        """
        Return the primitive adjoint of the quadratic form, which is
        the smallest discriminant integer-valued quadratic form whose
        matrix is a scalar multiple of the inverse of the matrix of
        the given quadratic form.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.adjoint_primitive()
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 3 -2 ]
            [ *  1 ]

        """
        return QuadraticForm(self.Hessian_matrix().adjoint_classical()).primitive()

    def dim(self):
        """
        Gives the number of variables of the quadratic form.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.dim()
            2
            sage: parent(Q.dim())
            Integer Ring
            sage: Q = QuadraticForm(Q.matrix())
            sage: Q.dim()
            2
            sage: parent(Q.dim())
            Integer Ring
        """
        return self.__n

    def base_ring(self):
        """
        Gives the ring over which the quadratic form is defined.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.base_ring()
            Integer Ring

        """
        return self.__base_ring


    def coefficients(self):
        """
        Gives the matrix of upper triangular coefficients,
        by reading across the rows from the main diagonal.

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.coefficients()
            [1, 2, 3]

        """
        return self.__coeffs


    def det(self):
        """
        Gives the determinant of the Gram matrix of 2*Q, or
        equivalently the determinant of the Hessian matrix of Q.

        (Note: This is always defined over the same ring as the
        quadratic form.)

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.det()
            8

        """
        try:
            return self.__det
        except AttributeError:
            ## Compute the determinant
            if self.dim() == 0:
                new_det = self.base_ring()(1)
            else:
                new_det = self.matrix().det()

            ## Cache and return the determinant
            self.__det = new_det
            return new_det


    def Gram_det(self):
        """
        Gives the determinant of the Gram matrix of Q.

        (Note: This is defined over the fraction field of the ring of
        the quadratic form, but is often not defined over the same
        ring as the quadratic form.)

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, [1,2,3])
            sage: Q.Gram_det()
            2

        """
        return self.det() / ZZ(2**self.dim())


    def base_change_to(self, R):
        """
        Alters the quadratic form to have all coefficients
        defined over the new base_ring R.  Here R must be
        coercible to from the current base ring.

        Note: This is preferable to performing an explicit
        coercion through the base_ring() method, which does
        not affect the individual coefficients.  This is
        particularly useful for performing fast modular
        arithmetic evaluations.

        INPUT:
            R -- a ring

        OUTPUT:
            quadratic form

        EXAMPLES::

            sage: Q = DiagonalQuadraticForm(ZZ,[1,1]); Q
            Quadratic form in 2 variables over Integer Ring with coefficients:
            [ 1 0 ]
            [ * 1 ]

        ::

            sage: Q1 = Q.base_change_to(IntegerModRing(5)); Q1
            Quadratic form in 2 variables over Ring of integers modulo 5 with coefficients:
            [ 1 0 ]
            [ * 1 ]

            sage: Q1([35,11])
            1

        """
        ## Check that a canonical coercion is possible
        if not is_Ring(R):
            raise TypeError("Oops!  R is not a ring. =(")
        if not R.has_coerce_map_from(self.base_ring()):
            raise TypeError("Oops!  There is no canonical coercion from " + str(self.base_ring()) + " to R.")

        ## Return the coerced form
        return QuadraticForm(R, self.dim(), [R(x) for x in self.coefficients()])


    def level(self):
        r"""
        Determines the level of the quadratic form over a PID, which is a
        generator for the smallest ideal `N` of `R` such that N * (the matrix of
        2*Q)^(-1) is in R with diagonal in 2*R.

        Over `\ZZ` this returns a non-negative number.

        (Caveat: This always returns the unit ideal when working over a field!)

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, range(1,4))
            sage: Q.level()
            8

            sage: Q1 = QuadraticForm(QQ, 2, range(1,4))
            sage: Q1.level()      # random
            UserWarning: Warning -- The level of a quadratic form over a field is always 1.  Do you really want to do this?!?
            1

            sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
            sage: Q.level()
            420

        """
        ## Try to return the cached level
        try:
            return self.__level
        except AttributeError:

            ## Check that the base ring is a PID
            if not isinstance(self.base_ring(), PrincipalIdealDomain):
                raise TypeError("Oops!  The level (as a number) is only defined over a Principal Ideal Domain.  Try using level_ideal().")


            ## Warn the user if the form is defined over a field!
            if self.base_ring().is_field():
                warn("Warning -- The level of a quadratic form over a field is always 1.  Do you really want to do this?!?")
                #raise RuntimeError, "Warning -- The level of a quadratic form over a field is always 1.  Do you really want to do this?!?"


            ## Check invertibility and find the inverse
            try:
                mat_inv = self.matrix()**(-1)
            except ZeroDivisionError:
                raise TypeError("Oops!  The quadratic form is degenerate (i.e. det = 0). =(")

            ## Compute the level
            inv_denoms = []
            for i in range(self.dim()):
                for j in range(i, self.dim()):
                    if (i == j):
                        inv_denoms += [denominator(mat_inv[i,j] / 2)]
                    else:
                        inv_denoms += [denominator(mat_inv[i,j])]
            lvl = LCM(inv_denoms)
            lvl = Ideal(self.base_ring()(lvl)).gen()
            ##############################################################
            ## To do this properly, the level should be the inverse of the
            ## fractional ideal (over R) generated by the entries whose
            ## denominators we take above. =)
            ##############################################################

            ## Normalize the result over ZZ
            if self.base_ring() == IntegerRing():
                lvl = abs(lvl)

            ## Cache and return the level
            self.__level = lvl
            return lvl

    def level_ideal(self):
        """
        Determines the level of the quadratic form (over R), which is the
        smallest ideal N of R such that N * (the matrix of 2*Q)^(-1) is
        in R with diagonal in 2*R.
        (Caveat: This always returns the principal ideal when working over a field!)

        WARNING:  THIS ONLY WORKS OVER A PID RING OF INTEGERS FOR NOW!
              (Waiting for Sage fractional ideal support.)

        EXAMPLES::

            sage: Q = QuadraticForm(ZZ, 2, range(1,4))
            sage: Q.level_ideal()
            Principal ideal (8) of Integer Ring

        ::

            sage: Q1 = QuadraticForm(QQ, 2, range(1,4))
            sage: Q1.level_ideal()
            Principal ideal (1) of Rational Field

        ::

            sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
            sage: Q.level_ideal()
            Principal ideal (420) of Integer Ring

        """
        ##############################################################
        ## To do this properly, the level should be the inverse of the
        ## fractional ideal (over R) generated by the entries whose
        ## denominators we take above. =)
        ##############################################################

        return Ideal(self.base_ring()(self.level()))

    def bilinear_map(self, v, w):
        r"""
        Return the value of the associated bilinear map on two vectors

        Given a quadratic form `Q` over some base ring `R` with
        characteristic not equal to 2, this gives the image of two
        vectors with coefficients in `R` under the associated bilinear
        map `B`, given by the relation `2 B(v,w) = Q(v) + Q(w) - Q(v+w)`.

        INPUT:

        `v, w` -- two vectors

        OUTPUT:

        an element of the base ring `R`.

        EXAMPLES:

        First, an example over `\ZZ`::

            sage: Q = QuadraticForm(ZZ,3,[1,4,0,1,4,1])
            sage: v = vector(ZZ,(1,2,0))
            sage: w = vector(ZZ,(0,1,1))
            sage: Q.bilinear_map(v,w)
            8

        This also works over `\QQ`::

            sage: Q = QuadraticForm(QQ,2,[1/2,2,1])
            sage: v = vector(QQ,(1,1))
            sage: w = vector(QQ,(1/2,2))
            sage: Q.bilinear_map(v,w)
            19/4

        The vectors must have the correct length::

            sage: Q = DiagonalQuadraticForm(ZZ,[1,7,7])
            sage: v = vector((1,2))
            sage: w = vector((1,1,1))
            sage: Q.bilinear_map(v,w)
            Traceback (most recent call last):
            ...
            TypeError: vectors must have length 3

        This does not work if the characteristic is 2::

            sage: Q = DiagonalQuadraticForm(GF(2),[1,1,1])
            sage: v = vector((1,1,1))
            sage: w = vector((1,1,1))
            sage: Q.bilinear_map(v,w)
            Traceback (most recent call last):
            ...
            TypeError: not defined for rings of characteristic 2
        """
        if len(v) != self.dim() or len(w) != self.dim():
            raise TypeError("vectors must have length " + str(self.dim()))
        if self.base_ring().characteristic() == 2:
            raise TypeError("not defined for rings of characteristic 2")
        return (self(v+w) - self(v) - self(w))/2

    genera = staticmethod(genera)

# ============================================================================


def DiagonalQuadraticForm(R, diag):
    """
    Return a quadratic form over `R` which is a sum of squares.

    INPUT:

    - `R` -- ring
    - ``diag`` -- list/tuple of elements coercible to R

    OUTPUT:

    quadratic form

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 0 0 ]
        [ * 3 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]
    """
    Q = QuadraticForm(R, len(diag))
    for i in range(len(diag)):
        Q[i, i] = diag[i]
    return Q
