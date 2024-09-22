# distutils: libraries = m NTL_LIBRARIES Lfunction
# distutils: extra_compile_args = NTL_CFLAGS -O3 -ffast-math
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
Rubinstein's lcalc library

This is a wrapper around Michael Rubinstein's lcalc.
See http://oto.math.uwaterloo.ca/~mrubinst/L_function_public/CODE/.

AUTHORS:

- Rishikesh (2010): added compute_rank() and hardy_z_function()
- Yann Laigle-Chapuy (2009): refactored
- Rishikesh (2009): initial version
"""

# ****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport str_to_bytes

from sage.libs.gmp.mpz cimport *
from sage.libs.mpfr cimport *
from sage.rings.integer cimport Integer

from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_mpfr import ComplexField
CCC = ComplexField()

from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr import RealField
RRR = RealField()
pi = RRR.pi()

initialize_globals()

##############################################################################
# Lfunction: base class for L-functions
##############################################################################

cdef class Lfunction:
    # virtual class
    def __init__(self, name, what_type_L, dirichlet_coefficient,
                 period, Q, OMEGA, gamma, lambd, pole, residue):
        """
        Initialization of `L`-function objects.
        See derived class for details, this class is not supposed to be
        instantiated directly.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: Lfunction_from_character(DirichletGroup(5)[1])
            L-function with complex Dirichlet coefficients
        """
        cdef int i              # for indexing loops
        cdef Integer tmpi       # for accessing integer values
        cdef RealNumber tmpr    # for accessing real values
        cdef ComplexNumber tmpc # for accessing complex values

        _name = str_to_bytes(name)
        cdef char *NAME = _name
        cdef int what_type = what_type_L

        tmpi = Integer(period)
        cdef int Period = mpz_get_si(tmpi.value)
        tmpr = RRR(Q)
        cdef double q=mpfr_get_d(tmpr.value, MPFR_RNDN)
        tmpc = CCC(OMEGA)
        cdef c_Complex w=new_Complex(mpfr_get_d(tmpc.__re, MPFR_RNDN), mpfr_get_d(tmpc.__im, MPFR_RNDN))

        cdef int A=len(gamma)
        cdef double *g=new_doubles(A+1)
        cdef c_Complex *l=new_Complexes(A+1)
        for i from 0 <= i < A:
            tmpr = RRR(gamma[i])
            g[i+1] = mpfr_get_d(tmpr.value, MPFR_RNDN)
            tmpc = CCC(lambd[i])
            l[i+1] = new_Complex(mpfr_get_d(tmpc.__re, MPFR_RNDN), mpfr_get_d(tmpc.__im, MPFR_RNDN))

        cdef int n_poles = len(pole)
        cdef c_Complex *p = new_Complexes(n_poles +1)
        cdef c_Complex *r = new_Complexes(n_poles +1)
        for i from 0 <= i < n_poles:
            tmpc=CCC(pole[i])
            p[i+1] = new_Complex(mpfr_get_d(tmpc.__re, MPFR_RNDN), mpfr_get_d(tmpc.__im, MPFR_RNDN))
            tmpc=CCC(residue[i])
            r[i+1] = new_Complex(mpfr_get_d(tmpc.__re, MPFR_RNDN), mpfr_get_d(tmpc.__im, MPFR_RNDN))

        self._init_fun(NAME, what_type, dirichlet_coefficient, Period, q,  w,  A, g, l, n_poles, p, r)

        if name:
            name += ': '

        self._repr = name + 'L-function'

        del_doubles(g)
        del_Complexes(l)
        del_Complexes(p)
        del_Complexes(r)

    def __repr__(self):
        """
        Return string representation of this `L`-function.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: Lfunction_from_character(DirichletGroup(5)[1])
            L-function with complex Dirichlet coefficients

            sage: Lfunction_Zeta()
            The Riemann zeta function
        """
        return self._repr

    def value(self, s, derivative=0):
        """
        Compute the value of the `L`-function at ``s``.

        INPUT:

        - ``s`` -- a complex number
        - ``derivative`` -- integer (default: 0);  the derivative to be evaluated
        - ``rotate`` -- boolean (default: ``False``); if True, this returns the value of the
          Hardy Z-function (sometimes called the Riemann-Siegel Z-function or
          the Siegel Z-function)

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] # This is a quadratic character
            sage: L = Lfunction_from_character(chi, type='int')
            sage: (L.value(0.5) - 0.231750947504016).abs() < 1e-8
            True
            sage: v = L.value(0.2 + 0.4*I)
            sage: (v - (0.102558603193 + 0.190840777924*I)).abs() < 1e-8
            True
            sage: L = Lfunction_from_character(chi, type='double')
            sage: (L.value(0.6) - 0.274633355856345).abs() < 1e-8
            True
            sage: v = L.value(0.6 + I)
            sage: (v - (0.362258705721 + 0.43388825062*I)).abs() < 1e-8
            True

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L = Lfunction_from_character(chi, type='complex')
            sage: v = L.value(0.5)
            sage: (v - (0.763747880117 + 0.21696476751*I)).abs() < 1e-8
            True
            sage: v = L.value(0.6 + 5*I)
            sage: (v - (0.702723260619 - 1.10178575243*I)).abs() < 1e-8
            True

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: L = Lfunction_Zeta()
            sage: (L.value(0.5) + 1.46035450880).abs() < 1e-8
            True
            sage: v = L.value(0.4 + 0.5*I)
            sage: (v - (-0.450728958517 - 0.780511403019*I)).abs() < 1e-8
            True
        """
        cdef ComplexNumber complexified_s = CCC(s)
        cdef c_Complex z = new_Complex(mpfr_get_d(complexified_s.__re, MPFR_RNDN), mpfr_get_d(complexified_s.__im, MPFR_RNDN))
        cdef c_Complex result = self._value(z, derivative)
        return CCC(result.real(),result.imag())

    def hardy_z_function(self, s):
        """
        Compute the Hardy Z-function of the `L`-function at s.

        INPUT:

        - ``s`` -- a complex number with imaginary part between -0.5 and 0.5

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2]  # Quadratic character
            sage: L = Lfunction_from_character(chi, type='int')
            sage: (L.hardy_z_function(0) - 0.231750947504).abs() < 1e-8
            True
            sage: L.hardy_z_function(0.5).imag().abs() < 1e-8
            True

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L = Lfunction_from_character(chi, type='complex')
            sage: (L.hardy_z_function(0) - 0.793967590477).abs() < 1e-8
            True
            sage: L.hardy_z_function(0.5).imag().abs() < 1e-8
            True

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: E = EllipticCurve([-82,0])
            sage: L = Lfunction_from_elliptic_curve(E, number_of_coeffs=40000)
            sage: (L.hardy_z_function(2.1) - (-0.006431791768)).abs() < 1e-8
            True
        """
        #This takes s -> .5 + I*s
        cdef ComplexNumber complexified_s = CCC(0.5)+ CCC(0,1)*CCC(s)
        cdef c_Complex z = new_Complex(mpfr_get_d(complexified_s.__re, MPFR_RNDN), mpfr_get_d(complexified_s.__im, MPFR_RNDN))
        cdef c_Complex result = self._hardy_z_function(z)
        return CCC(result.real(),result.imag())

    def compute_rank(self):
        """
        Compute the analytic rank (the order of vanishing at the center) of
        of the `L`-function.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] # This is a quadratic character
            sage: L = Lfunction_from_character(chi, type='int')
            sage: L.compute_rank()
            0

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: E = EllipticCurve([-82,0])
            sage: L = Lfunction_from_elliptic_curve(E, number_of_coeffs=40000)
            sage: L.compute_rank()
            3
        """
        return self._compute_rank()

    def _N(self, T):
        """
        Compute the number of zeroes upto height `T` using the formula for
        `N(T)` with the error of `S(T)`. Please do not use this. It is only
        for debugging

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L=Lfunction_from_character(chi, type='complex')
            sage: L._N(10) # abs tol 1e-8
            4.0
        """
        cdef RealNumber real_T=RRR(T)
        cdef double double_T = mpfr_get_d(real_T.value, MPFR_RNDN)
        cdef double res_d = self._typedN(double_T)
        return RRR(res_d)

    def find_zeros(self, T1, T2, stepsize):
        """
        Finds zeros on critical line between ``T1`` and ``T2`` using step size
        of stepsize. This function might miss zeros if step size is too
        large. This function computes the zeros of the `L`-function by using
        change in signs of areal valued function whose zeros coincide with
        the zeros of `L`-function.

        Use :meth:`find_zeros_via_N` for slower but more rigorous computation.

        INPUT:

        - ``T1`` -- a real number giving the lower bound
        - ``T2`` -- a real number giving the upper bound
        - ``stepsize`` -- step size to be used for the zero search

        OUTPUT: list of the imaginary parts of the zeros which were found

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] # This is a quadratic character
            sage: L = Lfunction_from_character(chi, type='int')
            sage: L.find_zeros(5,15,.1)
            [6.64845334472..., 9.83144443288..., 11.9588456260...]
            sage: L = Lfunction_from_character(chi, type='double')
            sage: L.find_zeros(1,15,.1)
            [6.64845334472..., 9.83144443288..., 11.9588456260...]

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L = Lfunction_from_character(chi, type='complex')
            sage: L.find_zeros(-8,8,.1)
            [-4.13290370521..., 6.18357819545...]

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: L = Lfunction_Zeta()
            sage: L.find_zeros(10,29.1,.1)
            [14.1347251417..., 21.0220396387..., 25.0108575801...]
        """
        cdef doublevec result
        cdef double myresult
        cdef int i
        cdef RealNumber real_T1 = RRR(T1)
        cdef RealNumber real_T2 = RRR(T2)
        cdef RealNumber real_stepsize = RRR(stepsize)
        sig_on()
        self._find_zeros_v( mpfr_get_d(real_T1.value, MPFR_RNDN), mpfr_get_d(real_T2.value, MPFR_RNDN), mpfr_get_d(real_stepsize.value, MPFR_RNDN),&result)
        sig_off()
        i=result.size()
        returnvalue = []
        for i in range(result.size()):
            returnvalue.append(  RRR(result.ind(i)))
        result.clear()
        return returnvalue

    #The default values are from L.h. See L.h
    def find_zeros_via_N(self, count=0, start=0, max_refine=1025, rank=-1):
        """
        Find ``count`` zeros (in order of increasing magnitude) and output
        their imaginary parts. This function verifies that no zeros
        are missed, and that all values output are indeed zeros.

        If this `L`-function is self-dual (if its Dirichlet coefficients
        are real, up to a tolerance of 1e-6), then only the zeros with
        positive imaginary parts are output. Their conjugates, which
        are also zeros, are not output.

        INPUT:

        - ``count`` -- number of zeros to be found
        - ``start`` -- (default: 0) how many initial zeros to skip
        - ``max_refine`` -- when some zeros are found to be missing, the step
          size used to find zeros is refined. max_refine gives an upper limit
          on when lcalc should give up. Use default value unless you know
          what you are doing.
        - ``rank`` -- integer (default: -1); analytic rank of the `L`-function.
          If -1 is passed, then we attempt to compute it. (Use default if in
          doubt)

        OUTPUT: list of the imaginary parts of the zeros that have been found

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L = Lfunction_from_character(chi, type='int')
            sage: L.find_zeros_via_N(3)
            [6.64845334472..., 9.83144443288..., 11.9588456260...]
            sage: L = Lfunction_from_character(chi, type='double')
            sage: L.find_zeros_via_N(3)
            [6.64845334472..., 9.83144443288..., 11.9588456260...]

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L = Lfunction_from_character(chi, type='complex')
            sage: zeros = L.find_zeros_via_N(3)
            sage: (zeros[0] - (-4.13290370521286)).abs() < 1e-8
            True
            sage: (zeros[1]  - 6.18357819545086).abs() < 1e-8
            True
            sage: (zeros[2]  - 8.45722917442320).abs() < 1e-8
            True

        ::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: L = Lfunction_Zeta()
            sage: L.find_zeros_via_N(3)
            [14.1347251417..., 21.0220396387..., 25.0108575801...]
        """

        # This is the default value for message_stamp, but we have to
        # pass it explicitly since we're passing in the next argument,
        # our &result pointer.
        cdef const char* message_stamp = ""
        cdef doublevec result
        sig_on()
        self._find_zeros(count, start, max_refine, rank, message_stamp, &result)
        sig_off()
        returnvalue = []
        for i in range(result.size()):
            returnvalue.append(  RRR(result.ind(i)))
        result.clear()
        return returnvalue

    # Needs to be overridden
    cdef void _init_fun(self, char *NAME, int what_type, dirichlet_coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r) noexcept:
        raise NotImplementedError

    cdef c_Complex _value(self,c_Complex s,int derivative) noexcept:
        raise NotImplementedError

    cdef c_Complex _hardy_z_function(self,c_Complex s) noexcept:
        raise NotImplementedError

    cdef int _compute_rank(self) noexcept:
        raise NotImplementedError

    cdef double _typedN(self,double T) noexcept:
        raise NotImplementedError

    cdef void _find_zeros_v(self,double T1, double T2, double stepsize, doublevec *result) noexcept:
        raise NotImplementedError

    cdef int _find_zeros(self, long count, long start, double max_refine, int rank, const char* message_stamp, doublevec *result) noexcept:
        raise NotImplementedError

##############################################################################
# Lfunction_I: L-functions with integer Dirichlet Coefficients
##############################################################################

cdef class Lfunction_I(Lfunction):
    r"""
    The ``Lfunction_I`` class is used to represent `L`-functions
    with integer Dirichlet Coefficients. We assume that `L`-functions
    satisfy the following functional equation.

    .. MATH::

        \Lambda(s) = \omega Q^s \overline{\Lambda(1-\bar s)}

    where

    .. MATH::

        \Lambda(s) = Q^s \left( \prod_{j=1}^a \Gamma(\kappa_j s + \gamma_j) \right) L(s)

    See (23) in :arxiv:`math/0412181`

    INPUT:

    - ``what_type_L`` -- integer; this should be set to 1 if the coefficients
      are periodic and 0 otherwise

    - ``dirichlet_coefficient`` -- list of Dirichlet coefficients of the
      `L`-function. Only first `M` coefficients are needed if they are periodic.

    - ``period`` -- if the coefficients are periodic, this should be the
      period of the coefficients

    - ``Q`` -- see above

    - ``OMEGA`` -- see above

    - ``kappa`` -- list of the values of `\kappa_j` in the functional equation

    - ``gamma`` -- list of the values of `\gamma_j` in the functional equation

    - ``pole`` -- list of the poles of `L`-function

    - ``residue`` -- list of the residues of the `L`-function

    .. NOTE::

         If an `L`-function satisfies `\Lambda(s) = \omega Q^s \Lambda(k-s)`,
         by replacing `s` by `s+(k-1)/2`, one can get it in the form we need.
    """

    def __init__(self, name, what_type_L, dirichlet_coefficient,
                 period, Q, OMEGA, gamma, lambd, pole, residue):
        r"""
        Initialize an `L`-function with integer coefficients.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L=Lfunction_from_character(chi, type='int')
            sage: type(L)
            <class 'sage.libs.lcalc.lcalc_Lfunction.Lfunction_I'>
        """
        Lfunction.__init__(self, name, what_type_L, dirichlet_coefficient,  period, Q, OMEGA, gamma,lambd, pole,residue)
        self._repr += " with integer Dirichlet coefficients"

    # override
    cdef void _init_fun(self, char *NAME, int what_type, dirichlet_coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r) noexcept:
        cdef int N = len(dirichlet_coeff)
        cdef Integer tmpi
        cdef int * coeffs = new_ints(N+1)  # lcalc ignores 0the coefficient
        for i from 0 <= i< N by 1:
            tmpi=Integer(dirichlet_coeff[i])
            coeffs[i+1] = mpz_get_si(tmpi.value)
        self.thisptr=new_c_Lfunction_I(NAME, what_type,  N, coeffs, Period, q,  w,  A, g, l, n_poles, p, r)
        del_ints(coeffs)

    cdef inline c_Complex _value(self,c_Complex s,int derivative) noexcept:
        return (<c_Lfunction_I *>(self.thisptr)).value(s, derivative, "pure")

    cdef inline c_Complex _hardy_z_function(self,c_Complex s) noexcept:
        return (<c_Lfunction_I *>(self.thisptr)).value(s, 0, "rotated pure")

    cdef int _compute_rank(self) noexcept:
        return (<c_Lfunction_I *>(self.thisptr)).compute_rank()

    cdef void _find_zeros_v(self, double T1, double T2, double stepsize, doublevec *result) noexcept:
        (<c_Lfunction_I *>self.thisptr).find_zeros_v(T1,T2,stepsize,result[0])

    cdef double _typedN(self, double T) noexcept:
        return (<c_Lfunction_I *>self.thisptr).N(T)

    cdef int _find_zeros(self, long count, long start, double max_refine, int rank, const char* message_stamp, doublevec *result) noexcept:
        (<c_Lfunction_I *>self.thisptr).find_zeros(count, start, max_refine, rank, message_stamp, result)

    # debug tools
    def _print_data_to_standard_output(self):
        """
        This is used in debugging. It prints out information from
        the C++ object behind the scenes. It will use standard output.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L=Lfunction_from_character(chi, type='int')
            sage: L._print_data_to_standard_output() # tol 1e-8
            -----------------------------------------------
            <BLANKLINE>
            Name of L_function:
            number of dirichlet coefficients = 5
            coefficients are periodic
            b[1] = 1
            b[2] = -1
            b[3] = -1
            b[4] = 1
            b[5] = 0
            <BLANKLINE>
            Q = 1.26156626101
            OMEGA = (1,0)
            a = 1 (the quasi degree)
            gamma[1] =0.5    lambda[1] =(0,0)
            <BLANKLINE>
            <BLANKLINE>
            number of poles (of the completed L function) = 0
            -----------------------------------------------
            <BLANKLINE>
        """
        (<c_Lfunction_I *>self.thisptr).print_data_L()

    def __dealloc__(self):
        """
        Deallocate memory used
        """
        del_c_Lfunction_I(<c_Lfunction_I *>(self.thisptr))

##############################################################################
# Lfunction_D: `L`-functions with double (real) Dirichlet Coefficients
##############################################################################

cdef class Lfunction_D(Lfunction):
    r"""
    The ``Lfunction_D`` class is used to represent `L`-functions
    with real Dirichlet coefficients. We assume that `L`-functions
    satisfy the following functional equation.

    .. MATH::

        \Lambda(s) = \omega Q^s \overline{\Lambda(1-\bar s)}

    where

    .. MATH::

        \Lambda(s) = Q^s \left( \prod_{j=1}^a \Gamma(\kappa_j s + \gamma_j) \right) L(s)

    See (23) in :arxiv:`math/0412181`

    INPUT:

    - ``what_type_L`` -- integer; this should be set to 1 if the coefficients are
      periodic and 0 otherwise

    - ``dirichlet_coefficient`` -- list of Dirichlet coefficients of the
      `L`-function. Only first `M` coefficients are needed if they are periodic.

    - ``period`` -- if the coefficients are periodic, this should be the
      period of the coefficients

    - ``Q`` -- see above

    - ``OMEGA`` -- see above

    - ``kappa`` -- list of the values of `\kappa_j` in the functional equation

    - ``gamma`` -- list of the values of `\gamma_j` in the functional equation

    - ``pole`` -- list of the poles of `L`-function

    - ``residue`` -- list of the residues of the `L`-function

    .. NOTE::

        If an `L`-function satisfies `\Lambda(s) = \omega Q^s \Lambda(k-s)`,
        by replacing `s` by `s+(k-1)/2`, one can get it in the form we need.
    """
    def __init__(self, name, what_type_L, dirichlet_coefficient,
                 period, Q, OMEGA, gamma, lambd, pole, residue):
        r"""
        Initialize an `L`-function with real coefficients.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L=Lfunction_from_character(chi, type='double')
            sage: type(L)
            <class 'sage.libs.lcalc.lcalc_Lfunction.Lfunction_D'>
        """
        Lfunction.__init__(self, name, what_type_L, dirichlet_coefficient,  period, Q, OMEGA, gamma,lambd, pole,residue)
        self._repr += " with real Dirichlet coefficients"

    # override
    cdef void _init_fun(self, char *NAME, int what_type, dirichlet_coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r) noexcept:
        cdef int i
        cdef RealNumber tmpr
        cdef int N = len(dirichlet_coeff)
        cdef double * coeffs = new_doubles(N+1)  # lcalc ignores 0th position
        for i in range(N):
            tmpr = RRR(dirichlet_coeff[i])
            coeffs[i+1] = mpfr_get_d(tmpr.value, MPFR_RNDN)
        self.thisptr=new_c_Lfunction_D(NAME, what_type,  N, coeffs, Period, q,  w,  A, g, l, n_poles, p, r)
        del_doubles(coeffs)

    cdef inline c_Complex _value(self,c_Complex s,int derivative) noexcept:
        return (<c_Lfunction_D *>(self.thisptr)).value(s, derivative, "pure")

    cdef inline c_Complex _hardy_z_function(self,c_Complex s) noexcept:
        return (<c_Lfunction_D *>(self.thisptr)).value(s, 0, "rotated pure")

    cdef inline int _compute_rank(self) noexcept:
        return (<c_Lfunction_D *>(self.thisptr)).compute_rank()

    cdef void _find_zeros_v(self, double T1, double T2, double stepsize, doublevec *result) noexcept:
        (<c_Lfunction_D *>self.thisptr).find_zeros_v(T1,T2,stepsize,result[0])

    cdef double _typedN(self, double T) noexcept:
        return (<c_Lfunction_D *>self.thisptr).N(T)

    cdef int _find_zeros(self, long count, long start,double max_refine, int rank, const char* message_stamp, doublevec *result) noexcept:
        (<c_Lfunction_D *>self.thisptr).find_zeros(count, start, max_refine, rank, message_stamp, result)

    # debug tools
    def _print_data_to_standard_output(self):
        """
        This is used in debugging. It prints out information from
        the C++ object behind the scenes. It will use standard output.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[2] #This is a quadratic character
            sage: L=Lfunction_from_character(chi, type='double')
            sage: L._print_data_to_standard_output() # tol 1e-8
            -----------------------------------------------
            <BLANKLINE>
            Name of L_function:
            number of dirichlet coefficients = 5
            coefficients are periodic
            b[1] = 1
            b[2] = -1
            b[3] = -1
            b[4] = 1
            b[5] = 0
            <BLANKLINE>
            Q = 1.26156626101
            OMEGA = (1,0)
            a = 1 (the quasi degree)
            gamma[1] =0.5    lambda[1] =(0,0)
            <BLANKLINE>
            <BLANKLINE>
            number of poles (of the completed L function) = 0
            -----------------------------------------------
            <BLANKLINE>
        """
        (<c_Lfunction_D *>self.thisptr).print_data_L()

    def __dealloc__(self):
        """
        Deallocate memory used
        """
        del_c_Lfunction_D(<c_Lfunction_D *>(self.thisptr))

##############################################################################
# Lfunction_C: L-functions with Complex Dirichlet Coefficients
##############################################################################

cdef class Lfunction_C:
    r"""
    The ``Lfunction_C`` class is used to represent `L`-functions
    with complex Dirichlet Coefficients. We assume that `L`-functions
    satisfy the following functional equation.

    .. MATH::

        \Lambda(s) = \omega Q^s \overline{\Lambda(1-\bar s)}

    where

    .. MATH::

        \Lambda(s) = Q^s \left( \prod_{j=1}^a \Gamma(\kappa_j s + \gamma_j) \right) L(s)

    See (23) in :arxiv:`math/0412181`

    INPUT:

    - ``what_type_L`` -- integer; this should be set to 1 if the coefficients are
      periodic and 0 otherwise

    - ``dirichlet_coefficient`` -- list of Dirichlet coefficients of the
      `L`-function. Only first `M` coefficients are needed if they are periodic.

    - ``period`` -- if the coefficients are periodic, this should be the
      period of the coefficients

    - ``Q`` -- see above

    - ``OMEGA`` -- see above

    - ``kappa`` -- list of the values of `\kappa_j` in the functional equation

    - ``gamma`` -- list of the values of `\gamma_j` in the functional equation

    - ``pole`` -- list of the poles of `L`-function

    - ``residue`` -- list of the residues of the `L`-function

    .. NOTE::

        If an `L`-function satisfies `\Lambda(s) = \omega Q^s \Lambda(k-s)`,
        by replacing `s` by `s+(k-1)/2`, one can get it in the form we need.
    """

    def __init__(self, name, what_type_L, dirichlet_coefficient,
                 period, Q, OMEGA, gamma, lambd, pole, residue):
        r"""
        Initialize an `L`-function with complex coefficients.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L=Lfunction_from_character(chi, type='complex')
            sage: type(L)
            <class 'sage.libs.lcalc.lcalc_Lfunction.Lfunction_C'>
        """
        Lfunction.__init__(self, name, what_type_L, dirichlet_coefficient,  period, Q, OMEGA, gamma,lambd, pole,residue)
        self._repr += " with complex Dirichlet coefficients"

    # override
    cdef void _init_fun(self, char *NAME, int what_type, dirichlet_coeff, long long Period, double q,  c_Complex w, int A, double *g, c_Complex *l, int n_poles, c_Complex *p, c_Complex *r) noexcept:
        cdef int i
        cdef int N = len(dirichlet_coeff)
        cdef ComplexNumber tmpc

        cdef c_Complex * coeffs = new_Complexes(N+1)
        coeffs[0]=new_Complex(0,0)
        for i from 0 <= i< N by 1:
            tmpc=CCC(dirichlet_coeff[i])
            coeffs[i+1] = new_Complex(mpfr_get_d(tmpc.__re, MPFR_RNDN), mpfr_get_d(tmpc.__im, MPFR_RNDN))

        self.thisptr = new_c_Lfunction_C(NAME, what_type,  N, coeffs, Period, q,  w,  A, g, l, n_poles, p, r)

        del_Complexes(coeffs)

    cdef inline c_Complex _value(self,c_Complex s,int derivative) noexcept:
        return (<c_Lfunction_C *>(self.thisptr)).value(s, derivative, "pure")

    cdef inline c_Complex _hardy_z_function(self,c_Complex s) noexcept:
        return (<c_Lfunction_C *>(self.thisptr)).value(s, 0,"rotated pure")

    cdef inline int _compute_rank(self) noexcept:
        return (<c_Lfunction_C *>(self.thisptr)).compute_rank()

    cdef void _find_zeros_v(self, double T1, double T2, double stepsize, doublevec *result) noexcept:
        (<c_Lfunction_C *>self.thisptr).find_zeros_v(T1,T2,stepsize,result[0])

    cdef double _typedN(self, double T) noexcept:
        return (<c_Lfunction_C *>self.thisptr).N(T)

    cdef int _find_zeros(self, long count, long start, double max_refine, int rank, const char* message_stamp, doublevec *result) noexcept:
        (<c_Lfunction_C *>self.thisptr).find_zeros(count, start, max_refine, rank, message_stamp, result)

    # debug tools
    def _print_data_to_standard_output(self):
        """
        This is used in debugging. It prints out information from
        the C++ object behind the scenes. It will use standard output.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: chi = DirichletGroup(5)[1]
            sage: L=Lfunction_from_character(chi, type='complex')
            sage: L._print_data_to_standard_output() # tol 1e-8
            -----------------------------------------------
            <BLANKLINE>
            Name of L_function:
            number of dirichlet coefficients = 5
            coefficients are periodic
            b[1] = (1,0)
            b[2] = (0,1)
            b[3] = (0,-1)
            b[4] = (-1,0)
            b[5] = (0,0)
            <BLANKLINE>
            Q = 1.26156626101
            OMEGA = (0.850650808352,0.525731112119)
            a = 1 (the quasi degree)
            gamma[1] =0.5    lambda[1] =(0.5,0)
            <BLANKLINE>
            <BLANKLINE>
            number of poles (of the completed L function) = 0
            -----------------------------------------------
            <BLANKLINE>
        """
        (<c_Lfunction_C *>self.thisptr).print_data_L()

    def __dealloc__(self):
        """
        Deallocate memory used
        """
        del_c_Lfunction_C(<c_Lfunction_C *>(self.thisptr))


##############################################################################
# Zeta function
##############################################################################

cdef class Lfunction_Zeta(Lfunction):
    r"""
    The ``Lfunction_Zeta`` class is used to generate the Riemann zeta function.
    """
    def __init__(self):
        r"""
        Initialize the Riemann zeta function.

        EXAMPLES::

            sage: from sage.libs.lcalc.lcalc_Lfunction import *
            sage: sage.libs.lcalc.lcalc_Lfunction.Lfunction_Zeta()
            The Riemann zeta function
        """
        self.thisptr = new_c_Lfunction_Zeta()
        self._repr = "The Riemann zeta function"

    cdef inline c_Complex _value(self,c_Complex s,int derivative) noexcept:
        return (<c_Lfunction_Zeta *>(self.thisptr)).value(s, derivative, "pure")

    cdef inline c_Complex _hardy_z_function(self,c_Complex s) noexcept:
        return (<c_Lfunction_Zeta *>(self.thisptr)).value(s, 0, "rotated pure")

    cdef inline int _compute_rank(self) noexcept:
        return (<c_Lfunction_Zeta *>(self.thisptr)).compute_rank()

    cdef void _find_zeros_v(self, double T1, double T2, double stepsize, doublevec *result) noexcept:
        (<c_Lfunction_Zeta *>self.thisptr).find_zeros_v(T1,T2,stepsize,result[0])

    cdef double _typedN(self, double T) noexcept:
        return (<c_Lfunction_Zeta *>self.thisptr).N(T)

    cdef int _find_zeros(self, long count, long start, double max_refine, int rank, const char* message_stamp, doublevec *result) noexcept:
        (<c_Lfunction_Zeta *>self.thisptr).find_zeros(count, start, max_refine, rank, message_stamp, result)

    def __dealloc__(self):
        """
        Deallocate memory used
        """
        del_c_Lfunction_Zeta(<c_Lfunction_Zeta *>(self.thisptr))


##############################################################################
# Tools
##############################################################################

def Lfunction_from_character(chi, type='complex'):
    """
    Given a primitive Dirichlet character, this function returns
    an lcalc `L`-function object for the `L`-function of the character.

    INPUT:

    - ``chi`` -- a Dirichlet character
    - ``use_type`` -- string (default: ``'complex'``); type used for the Dirichlet
      coefficients. This can be ``'int'``, ``'double'`` or ``'complex'``.

    OUTPUT: `L`-function object for ``chi``

    EXAMPLES::

        sage: from sage.libs.lcalc.lcalc_Lfunction import Lfunction_from_character
        sage: Lfunction_from_character(DirichletGroup(5)[1])
        L-function with complex Dirichlet coefficients
        sage: Lfunction_from_character(DirichletGroup(5)[2], type='int')
        L-function with integer Dirichlet coefficients
        sage: Lfunction_from_character(DirichletGroup(5)[2], type='double')
        L-function with real Dirichlet coefficients
        sage: Lfunction_from_character(DirichletGroup(5)[1], type='int')
        Traceback (most recent call last):
        ...
        ValueError: For non quadratic characters you must use type="complex"
    """
    if (not chi.is_primitive()):
        raise TypeError("Dirichlet character is not primitive")

    modulus=chi.modulus()
    if chi.is_even():
        a=0
    else:
        a=1

    Q=(RRR(modulus/pi)).sqrt()
    poles=[]
    residues=[]
    period=modulus
    OMEGA=1.0/ ( CCC(0,1)**a * (CCC(modulus)).sqrt()/chi.gauss_sum() )

    if type == "complex":
        dir_coeffs = [CCC(chi(n)) for n in range(1, modulus + 1)]
        return Lfunction_C("", 1,dir_coeffs, period,Q,OMEGA,[.5],[a/2.],poles,residues)
    if type not in ["double", "int"]:
        raise ValueError("unknown type")
    if chi.order() != 2:
        raise ValueError("For non quadratic characters you must use type=\"complex\"")
    if type == "double":
        dir_coeffs = [RRR(chi(n)) for n in range(1, modulus + 1)]
        return Lfunction_D("", 1,dir_coeffs, period,Q,OMEGA,[.5],[a/2.],poles,residues)
    if type == "int":
        dir_coeffs = [Integer(chi(n)) for n in range(1, modulus + 1)]
        return Lfunction_I("", 1,dir_coeffs, period,Q,OMEGA,[.5],[a/2.],poles,residues)


def Lfunction_from_elliptic_curve(E, number_of_coeffs=10000):
    """
    Given an elliptic curve E, return an `L`-function object for
    the function `L(s, E)`.

    INPUT:

    - ``E`` -- an elliptic curve
    - ``number_of_coeffs`` -- integer (default: 10000); the number of
      coefficients to be used when constructing the `L`-function object. Right
      now this is fixed at object creation time, and is not automatically
      set intelligently.

    OUTPUT: `L`-function object for ``L(s, E)``

    EXAMPLES::

        sage: from sage.libs.lcalc.lcalc_Lfunction import Lfunction_from_elliptic_curve
        sage: L = Lfunction_from_elliptic_curve(EllipticCurve('37'))
        sage: L
        L-function with real Dirichlet coefficients
        sage: L.value(0.5).abs() < 1e-8
        True
        sage: (L.value(0.5, derivative=1) - 0.305999773835200).abs() < 1e-6
        True
    """
    Q = RRR(E.conductor()).sqrt() / RRR(2 * pi)
    poles = []
    residues = []
    dir_coeffs = E.anlist(number_of_coeffs)
    dir_coeffs = [RRR(dir_coeffs[i]) / (RRR(i)).sqrt()
                  for i in range(1, number_of_coeffs)]
    OMEGA = E.root_number()
    return Lfunction_D("", 2, dir_coeffs, 0, Q, OMEGA, [1], [.5],
                       poles, residues)
