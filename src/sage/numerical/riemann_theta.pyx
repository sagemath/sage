#distutils: libraries = pari
r"""
Numerical computation of Riemann theta functions
================================================

This module implements arbitrary precision numerical computation of Riemann
theta functions with characteristics and their derivatives. We consider the following
definitions.

* Let `g` be a positive integer

* Let `\Omega` be a `g\times g` *Riemann matrix*, i.e., a symmetric complex matrix with positive definite imaginary part.

* A *characteristic* of level `N`, where `N` is a positive integer is a `2\times g` matrix with rows `\epsilon/N` and `\delta/N`. One normally only considers *reduced* characteristics, where the entries of `\epsilon,\delta` are in `\{0,\ldots,N-1\}`.

For a row vector `z\in \Bold{C}^g` the *Riemann theta function* of `\Omega` evaluated
at `z` is

.. MATH::
    \theta\begin{bmatrix} \epsilon/N\\\delta/N\end{bmatrix}(z,\Omega)=\sum_{n\in\Bold{Z}^g}e^{\pi i\left((n+\epsilon/N)\Omega(n+\epsilon/N)^T+2(n+\epsilon/N)(z+\delta/N)^T\right)}

In addition, we also consider partial derivatives of Riemann theta functions with respect
to the components of `z`.

See [DHBvHS2004]_ and [AC2019]_ for a description of the basic description of the summation
strategy and the relevant error bounds that allow for efficient computation.
The main features of the present implementation are:

* It allows for multiprecision computations

* It allows for characteristics and derivatives

* The implementation is particularly optimized for computing multiple partial derivatives of a Riemann theta function with given characteristic and evaluation point.

AUTHORS:

 - Nils Bruin, Sohrab Ganjian (2021-09-08): initial version

"""
# ****************************************************************************
#       Copyright (C) 2021 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.modules.free_module_element cimport FreeModuleElement, FreeModuleElement_generic_dense
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.rings.real_mpfr cimport RealNumber, RealField_class, RealField
from sage.libs.mpfr.types cimport mpfr_rnd_t, mpfr_t, mpfr_prec_t
from sage.libs.mpfr cimport *
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from sage.modules.vector_real_double_dense cimport Vector_real_double_dense
from sage.matrix.matrix_real_double_dense cimport Matrix_real_double_dense
from sage.rings.complex_mpc cimport *
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.libs.mpc cimport *
from cypari2.gen cimport Gen
from cypari2.types cimport *
from cypari2.paridecl cimport *

from sage.libs.pari import pari
from math import pi as double_pi
from math import sqrt as double_sqrt
from math import exp as double_exp
from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
import sage.libs.mpmath.all as mpall
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.arith.misc import binomial as binom
from sage.modules.vector_modn_dense import Vector_modn_dense

cdef class Vector_long:
    r"""
    Vector of system "long" integers

    This is only a very thin wrapper to give Cython code efficient access
    to an array of "long" that can be placed in Python data structures,
    so implemented functionality is minimal and almost none of it is
    available outside of cython.
    """
    cdef long n
    cdef long *vec

    def __cinit__(self, long n):
        r"""
        Allocate vector.

        INPUT:

        - ``n`` -- integer. Length of vector to allocate.

        OUTPUT: The allocated (but uninitialized!) vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_long
            sage: V = Vector_long(10)
        """
        self.n = n
        self.vec = <long *> PyMem_Malloc(n * sizeof(long))

    def __dealloc__(self):
        r"""
        Deallocate vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_long
            sage: V = Vector_long(10)
            sage: del V
        """
        PyMem_Free(self.vec)

    def __len__(self):
        r"""
        Return length of vector.

        OUTPUT: Length of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_long
            sage: V = Vector_long(10)
            sage: len(V)
            10

        """
        return self.n

    def __repr__(self):
        r"""
        Return string representation of vector.

        OUTPUT: String representation of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_long
            sage: V = Vector_long(4)
            sage: repr(V) # random
            '<Vector_long [139766604523248, 139774110506560, 1, 8704]>'

        """
        return "<Vector_long {}>".format(list(self))

    def __getitem__(self, i):
        r"""
        Return an entry from vector.

        INPUT:

        - ``i`` -- integer. Index of entry to retrieve.

        OUTPUT: entry value.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_long
            sage: V = Vector_long(4)
            sage: V[0] # random
            139766604523248

        """

        if i < 0 or i >= self.n:
            raise IndexError("index out of range")
        else:
            return self.vec[i]

    cdef assign(self, L):
        r"""
        Assign values from a sequence-type.

        INPUT:

        - ``L`` -- sequence. Must consist of elements convertible to system longs.

        This method is cdef only, so cannot be tested from python.
        """
        if len(L) != self.n:
            raise ValueError("mismatch in length")
        for i in range(self.n):
            self.vec[i] = L[i]

    @staticmethod
    cdef from_list(object L):
        r"""
        Return vector initialized from sequence-type.

        INPUT:

        - ``L`` -- sequence. Must consist of elements convertible to system longs.

        OUTPUT: Vector initialized with values from ``L``.

        Note that this method is a static method: it is supposed to be invoked
        via ``Vector_long.from_list(...)``. It combines allocation and initialization
        from a python data type.

        This method is cdef only, so cannot be tested from python.
        """
        cdef Vector_long v = Vector_long(len(L))
        v.assign(L)
        return v

    cdef assign_scaled_diff(self, long scaling, Vector_long v, Vector_long w):
        r"""
        Store the difference of scaled difference of two vectors

        Stores ``scaling*v-w`` in this vector.

        INPUT:

        - ``scaling`` -- system long. Scaling factor.
        - ``v`` -- Vector_long.
        - ``w`` -- Vector_long.

        This method is cdef only, so cannot be tested from python.
        """
        if self.n != v.n or self.n != w.n:
            raise ValueError("dimension mismatch")
        cdef long i
        for i in range(self.n):
            self.vec[i]=scaling*v.vec[i]-w.vec[i]


cdef class Vector_mpfr:
    r"""
    Vector of mpfr reals.

    This is only a very thin wrapper to give Cython code efficient access
    to an array of "mpfr_t" that can be placed in Python data structures,
    so implemented functionality is minimal and almost none of it is
    available outside of cython.
    """
    cdef long n
    cdef mpfr_t* vec
    cdef mpfr_rnd_t rnd
    cdef RealField_class RR
    cdef mpfr_prec_t prec

    def __cinit__(self, RealField_class RR, long n):
        r"""
        Allocate vector.

        INPUT:

        - ``R`` RealField_class -- Field to inherit precision and rounding from
        - ``n`` integer -- length of vector

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: V
            <Vector_mpfr [NaN, NaN, NaN]>

        """
        cdef long i
        self.RR = RR
        self.rnd = RR.rnd
        self.prec = RR._prec
        self.n = n
        self.vec = <mpfr_t *> PyMem_Malloc(n * sizeof(mpfr_t))
        for i in range(self.n):
            mpfr_init2(self.vec[i], self.prec)

    def __dealloc__(self):
        r"""
        Deallocate vector.

         EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: del V

        """
        cdef long i
        for i in range(self.n):
            mpfr_clear(self.vec[i])
        PyMem_Free(self.vec)

    def __len__(self):
        r"""Return length of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: len(V)
            3

        """
        return self.n

    def __getitem__(self, i):
        r"""
        Return an entry from vector.

        INPUT:

        - ``i`` -- integer. Index of entry to retrieve.

        OUTPUT: entry value.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: V[0]
            NaN

        """

        cdef RealNumber a
        if i < 0 or i >= self.n:
            raise IndexError("index out of range")
        else:
            a = self.RR._new()
            mpfr_set(a.value, self.vec[i], self.rnd)
            return a

    def __repr__(self):
        r"""
        Return string representation of vector

        OUTPUT: string representation

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: repr(V)
            '<Vector_mpfr [NaN, NaN, NaN]>'

        """
        return "<Vector_mpfr {}>".format(list(self))

    cdef assign(self, FreeModuleElement_generic_dense L):
        r"""
        Assign vector from sage vector.

        INPUT:

        - `L` -- sage vector over same RealField as this vector.

        This method is cdef only, so cannot be tested from python.
        """
        cdef long i
        cdef RealNumber a
        if self.RR is not L.base_ring():
            raise ValueError("parent mismatch")
        if self.n != len(L):
            raise ValueError("dimension mismatch")
        for i in range(self.n):
            a = L[i]
            mpfr_set(self.vec[i], a.value, self.rnd)

    @staticmethod
    cdef Vector_mpfr from_vector(FreeModuleElement_generic_dense L):
        r"""
        Return vector initialized from sage vector.

        INPUT:

        - ``L`` -- sage vector over same RealField as this vector.

        OUTPUT: Vector initialized with values from ``L``.

        Note that this method is a static method: it is supposed to be invoked
        via ``Vector_mpfr.from_vector(...)``. It combines allocation and initialization
        from a python data type.

        This method is cdef only, so cannot be tested from python.
        """
        cdef Vector_mpfr v = Vector_mpfr(L.base_ring(), len(L))
        v.assign(L)
        return v

    cdef assign_sum_si(self, Vector_mpfr v, Vector_long w):
        r"""
        Store sum ``v+w`` of a real and integer vector

        INPUT:

        - ``v`` -- Vector_mpfr.
        - ``w`` -- Vector_long.

        This method is cdef only, so cannot be tested from python.
        """
        cdef long i
        if self.n != v.n or self.n != w.n:
            raise ValueError("dimension mismatch")
        for i in range(self.n):
            mpfr_add_si(self.vec[i], v.vec[i], w.vec[i], self.rnd)


cdef class Vector_mpc:
    r"""
    Vector of mpc complex numbers.

    This is only a very thin wrapper to give Cython code efficient access
    to an array of "mpc_t" that can be placed in Python data structures,
    so implemented functionality is minimal and almost none of it is
    available outside of cython.
    """
    cdef mpc_t * vec
    cdef mpfr_rnd_t rnd
    cdef long n
    cdef ComplexNumber z

    def __cinit__(self, CC, long n):
        r"""
        Allocate vector.

        INPUT:

        - ``CC`` -- ComplexField. Field to inherit precision and rounding from.
        - ``n`` -- integer. Length of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpfr
            sage: RR = RealField(40)
            sage: V = Vector_mpfr(RR,3)
            sage: V
            <Vector_mpfr [NaN, NaN, NaN]>
        """
        cdef long i
        self.n = n
        self.z = CC.zero()
        self.rnd = (<RealField_class> CC._real_field()).rnd
        self.vec = <mpc_t *> PyMem_Malloc(n * sizeof(mpc_t))
        for i in range(self.n):
            mpc_init2(self.vec[i], self.z._prec)

    def __dealloc__(self):
        r"""
        Deallocate vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpc
            sage: CC = ComplexField(40)
            sage: V = Vector_mpc(CC,3)
            sage: del V
        """
        cdef long i
        for i in range(self.n):
            mpc_clear(self.vec[i])
        PyMem_Free(self.vec)

    def __len__(self):
        r"""
        Return length of vector

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpc
            sage: CC = ComplexField(40)
            sage: V = Vector_mpc(CC,3)
            sage: len(V)
            3
        """
        return self.n

    def __getitem__(self, long i):
        r"""
        Return an entry from vector.

        INPUT:

        - ``i`` -- integer. Index of entry to retrieve.

        OUTPUT: entry value.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpc
            sage: CC = ComplexField(40)
            sage: V = Vector_mpc(CC,3)
            sage: V[0]
            NaN + NaN*I
        """
        if i < 0 or i >= self.n:
            raise IndexError("index out of bounds")
        cdef ComplexNumber r = self.z._new()
        mpfr_set(r.__re, self.vec[i].re, self.rnd)
        mpfr_set(r.__im, self.vec[i].im, self.rnd)
        return r

    def __repr__(self):
        r"""
        Return string representation of vector

        OUTPUT: string representation

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import Vector_mpc
            sage: CC = ComplexField(40)
            sage: V = Vector_mpc(CC,3)
            sage: repr(V)
            '<Vector_mpc [NaN + NaN*I, NaN + NaN*I, NaN + NaN*I]>'
        """
        return "<Vector_mpc {}>".format(list(self))


cdef class NormCholesky:
    r"""Class for evaluating positive-definite norms.

    For this class, the norm is given by the lower triangular cholesky
    decomposition of the positive-definite Gram matrix of the norm. The class
    is aimed at providing highly optimized action from cython, so functionality
    on python level is limited.
    """
    cdef long n
    cdef RealField_class RR
    cdef mpfr_t r1, r2
    cdef mpfr_prec_t prec
    cdef mpfr_rnd_t rnd
    cdef mpfr_t *Clist

    def __cinit__(self, RealField_class RR, long n):
        r"""
        Allocate object.

        INPUT:

        - ``RR`` -- RealField. Field to inherit precision and rounding from.
        - ``n`` -- integer. Dimension of space on which the norm is defined.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormCholesky
            sage: RR = RealField(40)
            sage: nm = NormCholesky(RR,3)
        """
        cdef long k
        cdef long Clength
        self.n = n
        self.RR = RR
        self.prec = RR._prec
        self.rnd = RR.rnd
        mpfr_init2(self.r1, self.prec)  # initialize two registers r1,r2
        mpfr_init2(self.r2, self.prec)
        Clength = (n*(n+1))//2
        self.Clist = <mpfr_t*> PyMem_Malloc(Clength * sizeof(mpfr_t))
        for k in range(Clength):
            mpfr_init2(self.Clist[k], self.prec)

    def __dealloc__(self):
        r"""
        Deallocate object.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormCholesky
            sage: RR = RealField(40)
            sage: nm = NormCholesky(RR,3)
            sage: del nm
        """
        cdef long k
        mpfr_clear(self.r1)
        mpfr_clear(self.r2)
        for k in range((self.n*(self.n+1))//2):
            mpfr_clear(self.Clist[k])
        PyMem_Free(self.Clist)

    cdef assign(self, Matrix_generic_dense C):
        r"""
        Initialize norm from lower triangular Cholesky decomposition

        INPUT:

        - ``C`` -- real matrix. Matrix is assumed to be lower triangular.

        This method is cdef only, so cannot be tested from python.
        """
        cdef long i, j, k

        if self.RR is not C.base_ring():
            raise ValueError("parent mismatch")
        if C.nrows() != self.n or C.ncols() != self.n:
            raise ValueError("matrix not square or not of right dimension")

        k = 0
        for j in range(self.n):
            for i in range(j, self.n):
                mpfr_set(self.Clist[k], (<RealNumber> C.get_unsafe(i, j)).value, self.rnd)
                k += 1

    @staticmethod
    cdef NormCholesky from_cholesky_matrix(Matrix_generic_dense C):
        r"""
        Allocate and initialize norm from lower triangular Cholesky decomposition.

        Cython-level method.

        INPUT:

        - ``C`` -- real matrix. Matrix is assumed to be lower triangular.

        This method is cdef only, so cannot be tested from python.
        """
        cdef NormCholesky NC = NormCholesky(C.base_ring(), C.nrows())
        NC.assign(C)
        return NC

    @staticmethod
    def init(Matrix_generic_dense C):
        r"""
        Allocate and initialize norm from lower triangular Cholesky decomposition.

        Python-level wrapper.

        INPUT:

        - ``C`` -- real matrix. Matrix is assumed to be lower triangular.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormCholesky
            sage: RR = RealField(40)
            sage: C=matrix(RR,2,2,[1,0,1,1])
            sage: nm=NormCholesky.init(C)

        """
        return NormCholesky.from_cholesky_matrix(C)

    cdef mpfr_norm(self, mpfr_t s, mpfr_t *v):
        r"""
        Compute norm of a vector and place in preallocated mpfr_t.

        INPUT:

        - ``s`` -- mpfr_t. Location for result
        - ``v`` -- mpfr_t*. Vector to compute the norm of.

        This method is cdef only, so cannot be tested from python.
        """
        cdef long i, j, k

        mpfr_set_zero(s, 1)
        k = 0
        for j in range(self.n):
            mpfr_set_zero(self.r2, 1)
            for i in range(j, self.n):
                mpfr_mul(self.r1, v[i], self.Clist[k], self.rnd)
                k += 1
                mpfr_add(self.r2, self.r2, self.r1, self.rnd)
            mpfr_sqr(self.r1, self.r2, self.rnd)
            mpfr_add(s, s, self.r1, self.rnd)

    def __call__(self, FreeModuleElement_generic_dense v):
        r"""
        Return norm of vector.

        INPUT:

        - ``v`` -- vector. Vector to compute norm of.

        OUTPUT: Norm of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormCholesky
            sage: RR = RealField(40)
            sage: C=matrix(RR,2,2,[1,0,1,1])
            sage: nm=NormCholesky.init(C)
            sage: v=vector(RR,2,[2,3])
            sage: nm(v)
            34.000000000
            sage: v*C*C.T*v
            34.000000000
        """
        cdef RealNumber s
        cdef Vector_mpfr w
        s = self.RR._new()
        w = Vector_mpfr.from_vector(v)
        self.mpfr_norm(s.value, w.vec)
        return s


cdef class NormGramInt:
    r"""
    Class for computation of norms of integer vectors given by a real-valued Gram matrix.

    The class is aimed at providing highly optimized action from
    cython, so functionality on python level is limited.
    """

    cdef long n
    cdef RealField_class RR
    cdef mpfr_t r
    cdef mpfr_prec_t prec
    cdef Matrix_generic_dense G
    cdef mpfr_rnd_t rnd
    cdef mpfr_t *Glist

    def __cinit__(self, RealField_class RR, long n):
        r"""
        Allocate object.

        INPUT:

        - ``RR`` -- RealField. Field to inherit precision and rounding from.
        - ``n`` -- integer. Dimension of space on which the norm is defined.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormGramInt
            sage: RR = RealField(40)
            sage: nm=NormGramInt(RR,3)
        """
        cdef long k, Glength
        self.n = n
        self.RR = RR
        self.rnd = RR.rnd
        self.prec = RR._prec
        mpfr_init2(self.r, self.prec)
        Glength = (n*(n+1))//2
        self.Glist = <mpfr_t*> PyMem_Malloc(Glength * sizeof(mpfr_t))
        for k in range(Glength):
            mpfr_init2(self.Glist[k], self.prec)

    def __dealloc__(self):
        r"""
        Deallocate object.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormGramInt
            sage: RR = RealField(40)
            sage: nm=NormGramInt(RR,3)
            sage: del nm
        """
        cdef long k
        mpfr_clear(self.r)
        for k in range((self.n*(self.n+1)) // 2):
            mpfr_clear(self.Glist[k])
        PyMem_Free(self.Glist)

    cdef assign(self, Matrix_generic_dense G):
        r"""
        Initialize norm from Gram matrix.

        INPUT:

        - ``G`` -- real matrix. Matrix is assumed to be symmetric.

        This routine only accesses the lower triangular part of ``G``.

        This method is cdef only, so cannot be tested from python.
        """
        if self.RR is not G.base_ring():
            raise ValueError("base ring mismatch")
        if G.nrows() != self.n or G.ncols() != self.n:
            raise ValueError("dimension mismatch")
        cdef long k = 0
        for i in range(self.n):
            mpfr_set(self.Glist[k], (<RealNumber>G.get_unsafe(i, i)).value, self.rnd)
            k += 1
        for i in range(1, self.n):
            for j in range(i):
                mpfr_set(self.Glist[k], (<RealNumber>G.get_unsafe(i, j)).value, self.rnd)
                # note that we multiply the off-diagonal coefficients by 2
                mpfr_mul_si(self.Glist[k], self.Glist[k], 2, self.rnd)
                k += 1

    @staticmethod
    cdef NormGramInt from_gram_matrix(Matrix_generic_dense G):
        r"""
        Allocate and initialize norm from Gram matrix.

        Cython-level method.

        INPUT:

        - ``G`` -- real matrix. Matrix is assumed to be symmetric.

        This routine only accesses the lower triangular part of ``G``.

        This method is cdef only, so cannot be tested from python.
        """
        cdef NormGramInt NG = NormGramInt(G.base_ring(), G.nrows())
        NG.assign(G)
        return NG

    @staticmethod
    def init(Matrix_generic_dense G):
        r"""
        Allocate and initialize norm from Gram matrix.

        Python-level wrapper.

        INPUT:

        - ``G`` -- real matrix. Matrix is assumed to be lower triangular.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormGramInt
            sage: RR = RealField(40)
            sage: G=matrix(RR,2,2,[1,2,2,1])
            sage: nm=NormGramInt.init(G)

        """
        return NormGramInt.from_gram_matrix(G)

    cdef NormGramInt scaled_by(self, long C):
        r"""
        Return new norm object, scaled by `1/C`.

        INPUT:

        - ``C`` -- integer. Scaling factor

        OUTPUT: New norm object, with Gram matrix divided by ``C``.

        This method is cdef only, so cannot be tested from python.
        """
        cdef NormGramInt NG = NormGramInt(self.RR, self.n)
        cdef long k
        for k in range((self.n*(self.n+1)) // 2):
            mpfr_div_si(NG.Glist[k], self.Glist[k], C, self.rnd)
        return NG

    cdef mpfr_norm(self, mpfr_t result, long* v):
        r"""
        Compute norm of a vector and place in preallocated mpfr_t.

        INPUT:

        - ``s`` -- mpfr_t. Location for result
        - ``v`` -- long*. Vector to compute the norm of.

        This method is cdef only, so cannot be tested from python.
        """
        cdef long i, j, k
        cdef mpfr_t * Glist = self.Glist
        cdef mpfr_rnd_t rnd = self.rnd
        mpfr_set_zero(result, 1)
        k = 0
        for i in range(self.n):
            mpfr_mul_si(self.r, Glist[k], v[i]**2, rnd)
            k += 1
            mpfr_add(result, result, self.r, rnd)
        for i in range(1, self.n):
            for j in range(i):
                mpfr_mul_si(self.r, Glist[k], v[i]*v[j], rnd)
                k += 1
                mpfr_add(result, result, self.r, rnd)

    def __call__(self, w):
        r"""
        Return norm of vector.

        INPUT:

        - ``v`` -- vector. Vector to compute norm of.

        OUTPUT: Norm of vector.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import NormGramInt
            sage: RR = RealField(40)
            sage: G=matrix(RR,2,2,[1,2,2,1])
            sage: nm=NormGramInt.init(G)
            sage: v=vector(ZZ,2,[2,3])
            sage: nm(v)
            37.000000000
            sage: v*G*v
            37.000000000
        """
        if len(w) != self.n:
            raise ValueError("vector dimension mismatch")
        cdef Vector_long v = Vector_long.from_list(w)
        cdef RealNumber s = self.RR._new()  # reserve space for our sum
        self.mpfr_norm(s.value, v.vec)
        return s


def imag_func(a):
    r"""
    Return result of calling ``imag`` method on argument.

    INPUT:

    - ``a`` -- object

    OUTPUT: ``a.imag()``

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import imag_func
        sage: imag_func(CC.0)
        1.00000000000000
    """
    return a.imag()


def real_func(a):
    r"""
    Return result of calling ``real`` method on argument.

    INPUT:

    - ``a`` -- object

    OUTPUT: ``a.real()``

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import real_func
        sage: real_func(CC.0)
        0.000000000000000
    """
    return a.real()


def round_func(a):
    r"""
    Return result of calling ``round`` method on argument.

    INPUT:

    - ``a`` -- object

    OUTPUT: ``a.round()``

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import round_func
        sage: round_func(2.7)
        3

    """
    return a.round()


def Rbound(Y, tol):
    r"""
    Compute radius for Riemann theta function summation.

    See Theorem 2 in [DHBvHS2004]_.

    INPUT:

    - ``Y`` -- real matrix. Positive definite (imaginary part of a Riemann matrix)
    - ``tol`` -- real number. Tolerance for error allowed in theta function summation.

    OUTPUT: radius for summation to restrict error term to be within specified tolerance.

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import Rbound
        sage: Y = matrix(RR,2,2,[1,0,0,1])
        sage: Rbound(Y,10^(-10))
        5.70985786129...

    """
    cdef RealField_class RR
    cdef RealNumber rho, Rmin, R, fR, dfR, fRo, half_g, pi
    cdef long g

    RR = RealField(53)
    pi = RR.pi()
    g = Y.nrows()
    half_g = RR(g)/2
    rho = RR(pari.qfminim(Y, flag=2)[1])
    rho = (pi*rho).sqrt()

    def f(R):
        return (g*2**(g-1)* (RR(pari.incgam(half_g, (R-rho/2)**2))/rho**g)) - tol

    def df(R):
        return -2*(2*R-rho)**(g-1) * g * (-(2*R-rho)**2/4).exp()/rho**g

    Rmin = (RR(2*g).sqrt()+rho) / 2
    R = Rmin
    fR = f(R)
    # because the Gamma function is so highly convex, Newton iteration rather severely underestimates where the
    # zero lies. That means Newton iteration is actually pretty slow. It would help if we could start
    # with a much better initial approximation.
    if fR > 0:
        while True:
            fRo = fR
            dfR = df(R)
            R = R - fR/dfR
            fR = f(R)
            if fR.abs() >= fRo.abs():
                break
        assert R > Rmin
    return R


def Rbound_deriv(Y, N, tol):
    r"""
    Compute radius for Riemann theta function summation with derivatives.

    See Theorem 3.1 in [AC2019]_.

    INPUT:

    - ``Y`` -- real matrix. Positive definite (imaginary part of a Riemann matrix)
    - ``N`` -- integer. Order of derivative.
    - ``tol`` -- real number. Tolerance for error allowed in theta function summation.

    OUTPUT: radius for summation to restrict error term to be within specified tolerance.

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import Rbound_deriv
        sage: Y = matrix(RR,2,2,[1,0,0,1])
        sage: Rbound_deriv(Y,3,10^(-10))
        6.6689474473...
    """
    cdef RealField_class RR
    cdef RealNumber rho, Rmin, R, fR, dfR, fRo, pi
    cdef long g

    RR = RealField(53)
    pi = RR.pi()
    g = Y.nrows()

    rho = RR(pari.qfminim(Y, flag=2)[1])
    rho = (pi*rho).sqrt()
    Rmin = (((RR(g**2+8*N)).sqrt() + g + 2*N).sqrt() + rho)/2
    Yinv_norm = RR(1/min(Y.change_ring(RDF).SVD()[1].diagonal()).sqrt())

    gRR = RR(g)
    g_over_two = gRR/2
    g_sqrt = gRR.sqrt()

    C = [binom(N, i)*pi**(-i/2) * (Yinv_norm)**i * g_sqrt**(N-i) for i in range(N+1)]
    LC = (2*pi)**N * g_over_two * (2/rho)**g

    def f(R):
        return (LC * sum(C[i] * RR(pari.incgam((gRR+i)/2, (R-rho/2)**2)) for i in range(N+1))) - tol

    def df(R):
        return -2* LC * (-(R-rho/2)**2).exp() * sum(C[i] * (R-rho/2)**(g+i-1) for i in range(N+1))

    R = Rmin
    fR = f(R)
    if fR > 0:
        while True:
            fRo = fR
            dfR = df(R)
            R = R - fR/dfR
            fR = f(R)
            if fR.abs() >= fRo.abs():
                break
        assert R > Rmin
    return R


def cholesky_decomposition(G):
    r"""
    Return Cholesky decomposition of a positive definite real matrix.

    The Cholesky decomposition of a real positive definite matrix `G` is
    a lower triangular matrix `G` such that `G = C C^T`.

    This routine wraps the multiprecision implementation in the mp library.

    INPUT:

    - ``G`` -- real matrix. Positive definite.

    OUTPUT: The cholesky decomposition ``C`` such that ``G == C * C.T``

    EXAMPLES::

        sage: from sage.numerical.riemann_theta import cholesky_decomposition
        sage: RR = RealField(100)
        sage: G = matrix(RR, 3,3, [1,1/2,1/4,1/2,1,1/5,1/4,1/5,1])
        sage: C = cholesky_decomposition(G)
        sage: max(abs(a) for a in (C*C.T - G).list())
        0.00000000000000000000000000000
    """
    R = G.parent()
    prec = R.base_ring().prec()
    mpall.mp.prec = prec  # set work precision in "mp library"
    with mpall.workprec(prec):
        Cmp = mpall.matrix([mpall.sage_to_mpmath(list(c), prec) for c in G])
        M = mpall.cholesky(Cmp)
        C = R([mpall.mpmath_to_sage(c, prec) for c in M])
    return C


cdef class RiemannTheta:
    r"""
    Object for numerical computation of Riemann Theta functions with characteristics and derivatives

    INPUT:

    - ``Omega`` -- Complex matrix. The Riemann matrix for which to compute
        Riemann Theta functions. The precision of the base ring determines the default
        tolerance used in computing Riemann Theta function values.

    EXAMPLES:

    We go through a very simple example that illustrates the basic features.
    First we define the Riemann matrix and its RiemannTheta object::

        sage: from sage.numerical.riemann_theta import RiemannTheta
        sage: CC = ComplexField(80)
        sage: Omega = matrix(CC,2,2,[3*I,0,0,5*I])
        sage: RT = RiemannTheta(Omega)

    By default, the object evaluates the theta nullwerte, but we can specify
    the value of `z` at which we want to evaluate. Theta functions are fully
    periodic under translation by integer vectors::

        sage: RT()
        1.000161700487241998...
        sage: RT(z=(1,1))
        1.000161700487241998...

    The Riemann theta function is even, so its first order partial derivatives
    are zero. We also demonstrate we can compute them as a vector in one go.
    Note that partial derivatives are indicated by the index with respect to which
    the partial derivative should be taken::

        sage: RT(derivs=[0])
        0.00000000000000000000000
        sage: RT(derivs=[1])
        0.00000000000000000000000
        sage: RT(derivs=[[0],[1]])
        (0.00000000000000000000000, 0.00000000000000000000000)

    We can also compute higher order derivatives, so the following gives us the
    hessian matrix at ``z=0``::

        sage: matrix(2,2,RT(derivs=[[0,0],[0,1],[1,0],[1,1]]))
        [  -0.0063717804307107830675222     -0.00000000000000000000000]
        [    -0.00000000000000000000000 -0.000011900851943023954077279]

    Characteristics can be given as ``[[eps_1,...,eps_g],[delta_1,...,delta_g],N]]``,
    where ``N`` gives the level and the ``eps_i, delta_j`` are integers specifying
    the characteristic. Alternatively, one can give the characteristic as a
    ``2g``-dimensional vector over ``Z/NZ``::

        sage: c = [[1,0],[1,0],2]
        sage: v = vector(GF(2),[1,0,1,0])
        sage: RT(char=c)  # abs tol = 1e-24
        1.2071813646649472697563e-25 - 2.8725542860550068599620e-26*I
        sage: RT(char=v)  # abs tol = 1e-24
        1.2071813646649472697563e-25 - 2.8725542860550068599620e-26*I
        sage: RT(char=c, derivs=[[0],[1]]) # abs tol = 1e-24
        (-0.59552188399685576910149 - 1.1412196198763623205771e-49*I, -3.5185834728040112058953e-32 + 2.5226254523149252440284e-56*I)

    We check that for the genus 2 curve

    .. MATH::

        C: y^2=(x-2)(x-3)(x-5)(x-7)(x-11)(x-13),

    the gradients of the odd theta characteristics of level 2 are proportional to the
    roots `2, 3, 5, 7, 11, 13`. We give a period matrix for this curve relative to a
    cohomology basis that is defined over `Q`, derive the Riemann matrix, compute the
    gradients of the odd characteristics (with respect to the original cohomology basis!)
    and check that the ratios give us back the roots listed::

        sage: from sage.numerical.riemann_theta import RiemannTheta
        sage: from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
        sage: A = matrix(CC,2,4,[ -0.100985221999616*I, -0.0576741242160428*I, 0.170602500958820, 0.137052375058957, -0.257755342052576*I, -0.684089137456259*I, 0.685128296898840, 1.18843441146637])
        sage: Omega1 = A[:, :2]
        sage: Omega2 = A[:, 2:]
        sage: Omega1i = numerical_inverse(Omega1)
        sage: Omega = Omega1i*Omega2
        sage: odd = [v for v in GF(2)^4 if v[:2]*v[2:] == 1]
        sage: RT = RiemannTheta(Omega)
        sage: gradients = [vector(RT(derivs=[[0],[1]],char=c))*Omega1i for c in odd]
        sage: sorted([(-v[0]/v[1]).algdep(1) for v in gradients])
        [x - 13, x - 11, x - 7, x - 5, x - 3, x - 2]

    """
    cdef long g
    cdef long nvec
    cdef Matrix_generic_dense Y
    cdef Matrix_generic_dense X
    cdef RealField_class RR
    cdef mpfr_prec_t prec
    cdef object CC
    cdef RealNumber pi
    cdef Matrix_generic_dense Yinv
    cdef Matrix_generic_dense T
    cdef NormCholesky Ynorm
    cdef dict Xnorm_dict
    cdef MPComplexField_class MPCC
    cdef mpc_t c1, c2
    cdef mpfr_t r1, r2, r3
    cdef object ZZg, CCg
    cdef dict Rbound_dict

    def __init__(self, Matrix_generic_dense Omega):
        r"""
        Initialize object.

        INPUT:

        - ``Omega`` -- complex matrix. The Riemann matrix.

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import RiemannTheta
            sage: RT = RiemannTheta(matrix(CC,2,2,[2*I,0,0,3*I]))

        """
        self.CC = Omega.base_ring()
        self.pi = self.RR.pi()
        self.Yinv = numerical_inverse(self.Y)

        # be careful: we scale the cholesky decomposition to already have the pi in it.
        self.T = cholesky_decomposition(self.pi*self.Y)
        self.Ynorm = NormCholesky.from_cholesky_matrix(self.T)
        self.Xnorm_dict = {1: NormGramInt.from_gram_matrix(self.X)}
        self.Rbound_dict = {}
        self.MPCC = MPComplexField_class(self.prec)
        self.ZZg = ZZ**(self.g)
        self.CCg = (self.CC)**self.g

    def __cinit__(self, Matrix_generic_dense Omega):
        r"""
        Initialize object.

        INPUT:

        - ``Omega`` -- complex matrix. The Riemann matrix.

        Cython level initialization (mainly allocation)

        This method is cdef only, so cannot be tested from python.
        """
        self.Y = Omega.apply_map(imag_func)
        self.X = Omega.apply_map(real_func)
        self.RR = self.Y.base_ring()
        self.prec = self.RR._prec
        self.g = self.Y.nrows()

        mpfr_init2(self.r1, self.prec)
        mpfr_init2(self.r2, self.prec)
        mpfr_init2(self.r3, self.prec)
        mpc_init2(self.c1, self.prec)
        mpc_init2(self.c2, self.prec)

    def _dealloc__(self):
        r"""
        Deallocate object.

        This method is cdef only, so cannot be tested from python.
        """
        mpfr_clear(self.r1)
        mpfr_clear(self.r2)
        mpfr_clear(self.r3)
        mpc_clear(self.c1)
        mpc_clear(self.c2)

    cdef NormGramInt scaled_Xnorm(self, long N):
        r"""
        Return a scaled X-GramNormInt object.

        When computing characteristics of level ``N`` we end up scaling
        the denominator of our characteristics into the Gram matrices we
        compute, so that we can leave our input as small integers.
        This routine caches the objects for efficiency.

        INPUT:

        - ``N`` -- scaling factor

        OUTPUT: a ``GramNormInt`` object scaled by 1/N^2

        This method is cdef only, so cannot be tested from python.
        """
        try:
            return self.Xnorm_dict[N]
        except KeyError:
            r = (<NormGramInt> self.Xnorm_dict[1]).scaled_by(N*N)
            self.Xnorm_dict[N] = r
            return r

    def __call__(self, z=None, object char=None, derivs=[], tol=None):
        r"""
        Evaluate Riemann theta function with characteristic and derivatives.

        INPUT:

        - ``z`` -- vector; optional (default 0). Point to evaluate theta function at.
        - ``char`` -- list or vector; optional (default 0). Characteristic.
            The characteristic can either be specified as a list
            ``[[eps_1, ...,eps_g],[delta_1, ..., delta_g],N]``, where ``N`` is the level
            of the characteristic and the ``eps_i, delta_j`` are integers describing the
            characteristic, or as a ``2*g`-dimensional vector over ``ZZ/ N*ZZ``. In the latter case,
            the level ``N`` is read off from the base ring and the vector is taken to be
            ``[eps_1, ..., eps_g, delta_1, ..., delta_g]``.
        - ``derivs`` -- list; optional (default []). Derivatives. It can be a list
            of integers ``[i_1,...,i_n]``, in which case it is taken to mean the derivative of
            order ``n``, obtained by taking the partial derivative with respect to
            ``z[i_1], ..., z[i_n]``. Alternatively, it can be a list of lists of integers,
            in which case the values of the derivatives indicated by the members of the list
            are returned as a tuple, in order.
        - ``tol`` -- real number; optional. Tolerance allowed in approximation. The default is
            the tolerance indicated by the precision of the base ring. Note that the tolerance
            controlled is the tolerance in the approximation of the periodic part of the Riemann
            Theta function (see [DHBvHS2004]_). Furthermore, floating point rounding in iterated
            summations may perturb the lower bits.

        OUTPUT: A complex number of a tuple of them; the value(s) of the indicated Riemann
        Theta function(s).

        EXAMPLES::

            sage: from sage.numerical.riemann_theta import RiemannTheta
            sage: RT = RiemannTheta(matrix(CC,2,2,[2*I,0,0,3*I]))
            sage: RT(z=(0,0),char=[[1,0],[0,1],2],derivs=[0,0,0]).abs() # abs tol = 1e-15
            2.88494706892332e-16
        """
        if z is None:
            z = self.CCg.zero()
        else:
            z = self.CCg(z)

        if char is None:
            eps = self.ZZg.zero()
            delta = self.ZZg.zero()
            N = 1
        elif isinstance(char, FreeModuleElement):
            R = char.base_ring()
            rnk = char.parent().rank()
            # note R.characteristic is different
            if R.characteristic() != R.order() or rnk != 2*self.g:
                raise TypeError("invalid characteristic specification")
            eps = self.ZZg(char[:self.g])
            delta = self.ZZg(char[self.g:])
            N = R.order()
        else:
            if len(char) != 3 or len(char[0]) != self.g or len(char[1]) != self.g:
                raise TypeError("invalid characteristic specification")
            eps = self.ZZg(char[0])
            delta = self.ZZg(char[1])
            N = int(char[2])

        if len(derivs) == 0:
            vecresult = False
            derivs = [[]]
        else:
            try:
                len(derivs[0]) > 0
                vecresult = True
            except TypeError:
                vecresult = False
                derivs=[derivs]

        if any(d < 0 or d >= self.g for l in derivs for d in l):
            raise ValueError("invalid value in derivative list")
        derivs = [Vector_long.from_list(d) for d in derivs]

        if tol is None:
            tol = self.RR(2)**(-self.prec)
        result = self._eval_vector_(z, eps, delta, N, derivs, tol)
        if vecresult:
            return tuple(result)
        else:
            return result[0]

    cdef double Rbound(self, maxnderiv, tol):
        r"""
        Return radius bound for summation

        Value is cached for efficiency.

        INPUT:

        - ``maxnderiv`` -- integer. (Maximal) order of derivative for which to compute bound.
        - ``tol`` -- real number. Tolerance allowed.

        OUTPUT: radius for summation.

        This method is cdef only, so cannot be tested from python.
        """
        key = (maxnderiv, tol)
        try:
            return self.Rbound_dict[key]
        except KeyError:
            pass
        if maxnderiv == 0:
            R = Rbound(self.Y, tol)
        else:
            R = Rbound_deriv(self.Y, maxnderiv, tol)
        self.Rbound_dict[key] = R
        return R

    cdef _eval_vector_(self, FreeModuleElement_generic_dense z, Vector_integer_dense eps, Vector_integer_dense delta, long N, list derivs, RealNumber tol):
        r"""
        Return computed evaluation(s)

        This routine implements the actual computation. You should probably call
        it through the ``__call__`` wrapper instead of directly.

        INPUT:

        - ``z`` -- vector over complex field. Evaluation point.
        - ``eps`` -- integer vector. Part of characteristic.
        - ``delta`` -- integer vector. Part of characteristic.
        - ``N`` -- integer. Level of characteristic.
        - ``derivs`` -- list of lists of integers. Derivatives to compute.
        - ``tol`` -- real number. Tolerance permitted in approximations.

        OUTPUT: Tuple of computed theta values.

        This method is cdef only, so cannot be tested from python.

        """
        # allocations and unpacking to local variables
        cdef long g = self.g
        cdef mpfr_rnd_t rnd = self.RR.rnd
        cdef mpc_rnd_t mpc_rnd = self.MPCC.__rnd
        cdef Vector_long nvec = Vector_long(g)
        cdef Vector_mpfr wvec = Vector_mpfr(self.RR, g)
        cdef Vector_long derivvec
        cdef long npoints, i, j, k, flag
        if eps is None:
            eps = self.ZZg.zero()
        else:
            eps = self.ZZg(eps)
        if delta is None:
            delta = self.ZZg.zero()
        else:
            delta = self.ZZg(delta)

        cdef NormGramInt Xnorm = self.scaled_Xnorm(N)

        # first unpacking and set-up computations
        cdef FreeModuleElement_generic_dense x = z.apply_map(real_func)
        cdef FreeModuleElement_generic_dense y = z.apply_map(imag_func)
        cdef FreeModuleElement_generic_dense Yinv_y = (self.Yinv * y)
        cdef Vector_integer_dense roundYinv_y = Yinv_y.apply_map(round_func)
        cdef FreeModuleElement_generic_dense fracYinv_y = Yinv_y - roundYinv_y
        cdef FreeModuleElement_generic_dense c = fracYinv_y + eps/N

        # store some quantities for fast access
        # note that xvec is actually (2/N)*(x+delta). This is the scaling
        # with which it gets used later.
        cdef Vector_mpfr xvec = Vector_mpfr.from_vector((2/N)*(x+delta/N))
        cdef Vector_mpfr cvec = Vector_mpfr.from_vector(c)
        cdef Vector_long Netavec = Vector_long.from_list(N*roundYinv_y-eps)

        # allocate space for return values
        cdef Vector_mpc s = Vector_mpc(self.CC, len(derivs))
        for i in range(s.n):
            mpc_set_si(s.vec[i], 0, mpc_rnd)

        # compute the enumeration radius
        cdef long maxnderiv = max((<Vector_long> l).n for l in derivs)
        cdef double R = self.Rbound(maxnderiv, tol)
        cdef double Rsqr = R**2

        # use pari's Finke-Pohst implementation to get lattice points.
        # we need points enumerated in a ball centered at c, but
        # the implementation only supports balls centered at 0. Hence we
        # increase the radius to make sure we include the ball we want.
        # We later throw out the points we do not need.
        cdef Gen V
        npoints, _, V = pari.qfminim(self.Y, (R+(self.Ynorm(c)/self.pi).sqrt())**2/double_pi, flag=2)
        npoints = npoints // 2

        # note that Pari's Finke-Pohst only returns one representative for each
        # {v,-v} pair and leaves out the 0 vector. So we do the 0-vector separately
        # and keep track of the sign we looked at with a flag that will toggle
        # through each iteration.
        for i in range(g):
            nvec.vec[i] = 0
        j = 0
        flag = 1

        # The following loop does the actual summation. The arithmetic is
        # written out in explicit mpfr-calls for optimized performance
        # (it saves us allocation overhead, because we are using preallocated
        # memory).
        while True:
            # w = n + c, where c = [[Y^(-1)y]] + epsilon
            wvec.assign_sum_si(cvec, nvec)
            # r1 = w*Y*w
            self.Ynorm.mpfr_norm(self.r1, wvec.vec)
            # based on the norm of r1, we can tell
            # if the point lies in our ball.
            if mpfr_get_d(self.r1, rnd) <= Rsqr:
                # r1 = -w*Y*w
                mpfr_neg(self.r1, self.r1, rnd)
                # nvec = N*n - N*eta = N*(n - [Y^-1y] + epsilon)
                nvec.assign_scaled_diff(N, nvec, Netavec)
                # r2 = (n-eta)*X*(n-eta) [note: Xnorm is scaled by 1/N^2]
                Xnorm.mpfr_norm(self.r2, nvec.vec)

                # r2 = (n-eta)*X*(n-eta) + 2*(x+delta)*(n-eta)
                # [note: xvec is scaled for this]
                for i in range(self.g):
                    mpfr_mul_si(self.r3, xvec.vec[i], nvec.vec[i], rnd)
                    mpfr_add(self.r2, self.r2, self.r3, rnd)

                # r3 = pi*((n-eta)*X*(n-eta) + 2*(x+delta)*(n-eta))
                mpfr_mul(self.r3, self.r2, self.pi.value, rnd)

                # c1 = exp(pi*i*((n-eta)*X*(n-eta) + 2*(x+delta)*(n-eta)))
                mpfr_sin_cos(self.c1.im, self.c1.re, self.r3, rnd)

                # r2 = exp(r1) = exp (-w*Y*w)
                mpfr_exp(self.r2, self.r1, rnd)

                # c1 = c1*r2 [note: this is the summation term for no derivatives]
                mpc_mul_fr(self.c1, self.c1, self.r2, mpc_rnd)

                # loop through the derivative descriptions
                for k in range(s.n):
                    derivvec = derivs[k]
                    mpc_set(self.c2, self.c1, mpc_rnd)
                    # multiply with appropriate product.
                    for i in range(derivvec.n):
                        mpc_mul_si(self.c2, self.c2, nvec[derivvec.vec[i]], mpc_rnd)

                    # add c1 to the running sum s
                    mpc_add(s.vec[k], s.vec[k], self.c2, mpc_rnd)
            # go to next iteration:
            # if flag == 1, go to next point; otherwise negate current point.
            if flag:
                j += 1
                if j > npoints:
                    break
                flag = 0
                for i in range(g):
                    nvec.vec[i] = itos(gcoeff(V.g, i+1, j))
            else:
                flag = 1
                for i in range(g):
                    nvec.vec[i] = -itos(gcoeff(V.g, i+1, j))

        # compute inner product y*Yinv_y
        # set r1 = 0
        mpfr_set_si(self.r1, 0, rnd)
        for i in range(g):
            # r1 += y[i]*Yinv_y[i]
            mpfr_fma(self.r1, (<RealNumber>y.get_unsafe(i)).value, (<RealNumber>Yinv_y.get_unsafe(i)).value, self.r1, rnd)
        # r2 = r1*pi = pi*y*Yinv_y
        mpfr_mul(self.r2, self.r1, self.pi.value, rnd)
        # r1 = exp (pi*y*Yinv_y)
        mpfr_exp(self.r1, self.r2, rnd)
        # s *= r1
        for k in range(s.n):
            derivvec = derivs[k]
            mpc_mul_fr(s.vec[k], s.vec[k], self.r1, mpc_rnd)
            # r3 = 2*pi
            mpfr_mul_si(self.r3, self.pi.value, 2, rnd)
            # r2 = (2*pi)/N
            mpfr_div_si(self.r2, self.r3, N, rnd)
            # r3 = (2*pi/N)^nderiv
            mpfr_pow_si(self.r3, self.r2, derivvec.n, rnd)
            # c1 = i^nderiv
            mpc_rootofunity(self.c1, 4, derivvec.n % 4, mpc_rnd)
            # c2 = (2*pi*i/N)^nderiv
            mpc_mul_fr(self.c2, self.c1, self.r3, mpc_rnd)
            # s *= c2
            mpc_mul(s.vec[k], s.vec[k], self.c2, mpc_rnd)
        return s
