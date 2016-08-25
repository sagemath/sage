r"""
fpLLL library

Wrapper for the fpLLL library by Damien Stehle, Xavier Pujol and David Cade
found at http://perso.ens-lyon.fr/damien.stehle/fplll/.

This wrapper provides access to fpLLL's LLL, BKZ and enumeration
implementations.

AUTHORS:

- Martin Albrecht (2007-10) initial release
- Martin Albrecht (2014-03) update to fpLLL 4.0 interface

"""

#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#       Copyright (C) 2014 Martin Albrecht <martinralbrecht@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "cysignals/signals.pxi"

from sage.libs.gmp.mpz cimport *
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.misc.superseded import deprecation

from sage.misc.decorators import rename_keyword

cdef inline int _check_precision(int precision) except -1:
    """
    Check whether the provided precision is within valid bounds. If not raise
    a ``TypeError``.

    INPUT:

    - ``precision`` -- an integer

    """
    if precision < 0:
        raise TypeError("precision must be >= 0")

cdef inline int _check_eta(float eta) except -1:
    r"""
    Check whether the provided parameter `\eta` is within valid bounds. If
    not raise a ``TypeError``.

    INPUT:

    - ``eta`` -- a floating point number

    """
    if eta < 0.5:
        raise TypeError("eta must be >= 0.5")

cdef inline int _check_delta(float delta) except -1:
    """
    Check whether the provided parameter `\delta` is within valid bounds. If
    not raise a ``TypeError``.

    INPUT:

    - ``delta`` -- a floating point number

    """
    if delta <= 0.25:
        raise TypeError("delta must be > 0.25")
    elif delta > 1.0:
        raise TypeError("delta must be <= 1.0")

cdef inline FloatType check_float_type(object float_type):
    cdef FloatType float_type_

    if float_type == "default" or float_type is None:
        float_type_= FT_DEFAULT
    elif float_type == "double":
        float_type_ = FT_DOUBLE
    elif float_type == "long double":
        float_type_ = FT_LONG_DOUBLE
    elif float_type == "dpe":
        float_type_ = FT_DPE
    elif float_type == "mpfr":
        float_type_ = FT_MPFR
    else:
        raise ValueError("Float type '%s' unknown" % float_type)
    return float_type_

cdef class FP_LLL:
    """
    A basic wrapper class to support conversion to/from Sage integer matrices
    and executing LLL/BKZ computations.

    .. note:

      Usually you don't want to create this object yourself but use the LLL
      method of the integer matrices/lattices.

    """


    def __cinit__(self, Matrix_integer_dense A):
        """
        Construct a new fpLLL wrapper for the matrix ``A``.

        INPUT:

        - ``A`` -- a matrix over the integers

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10)
            sage: FP_LLL(A)
            fpLLL wrapper acting on a 10 x 10 matrix

            sage: A = matrix(ZZ, 2, 0)
            sage: FP_LLL(A).fast()
            doctest:...: DeprecationWarning: You can just call LLL() instead
            See http://trac.sagemath.org/15976 for details.

            sage: A = matrix(ZZ, 0, 2)
            sage: FP_LLL(A)
            Traceback (most recent call last):
            ...
            ValueError: fpLLL cannot handle matrices with zero rows.

        """
        cdef unsigned long i,j

        if A._nrows == 0:
            raise ValueError('fpLLL cannot handle matrices with zero rows.')

        self._lattice = new ZZ_mat[mpz_t](A._nrows,A._ncols)
        cdef mpz_t tmp
        mpz_init(tmp)
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < A._ncols:
                A.get_unsafe_mpz(i,j,tmp)
                # mpz_set(self._lattice[0][i][j],tmp)
                self._lattice[0][i][j].set(tmp)
        mpz_clear(tmp)

    def __dealloc__(self):
        """
        Destroy internal data.
        """
        del self._lattice

    def __repr__(self):
        """
        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10)
            sage: FP_LLL(A) # indirect doctest
            fpLLL wrapper acting on a 10 x 10 matrix
        """
        return "fpLLL wrapper acting on a %d x %d matrix" % (
            self._lattice.getRows(), self._lattice.getCols())

    def _sage_(self):
        """
        Return a Sage representation of this matrix.

        OUTPUT:

        A Sage representation of this matrix.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10)
            sage: fpLLL = FP_LLL(A)
            sage: fpLLL._sage_() == A
            True
        """
        return to_sage(self._lattice)

    @rename_keyword(deprecation=19572, method='algorithm')
    def LLL(self, float delta=LLL_DEF_DELTA, float eta=LLL_DEF_ETA,
            algorithm=None, float_type=None, int precision=0,
            verbose=False, siegel=False, early_red=False):
        r"""
        `(\delta,\eta)`-LLL reduce this lattice.

        INPUT:

        - ``delta`` -- (default: ``0.99``) parameter `0.25 < \delta < 1.0`
        - ``eta `` -- (default: ``0.51``) parameter `0.5 \leq \eta <
          \sqrt{\delta}`
        - ``algorithm`` -- (default: ``None``) can be one of the following:

          * ``'wrapper'`` (``None``)
          * ``'proved'``
          * ``'fast'``
          * ``'heuristic'``

        - ``float_type`` -- (default: ``None``) can be one of the following:

          * ``None`` - for automatic choice
          * ``'double'``
          * ``'long double'``
          * ``'dpe'``
          * ``'mpfr'``

        - ``precision`` -- (default: ``0`` for automatic choice) precision
          to use
        - ``verbose`` -- (default: ``False``) be verbose
        - ``siegel`` -- (default: ``False``) use Siegel conditioning
        - ``early_red`` -- (default: ``False``) use early reduction

        OUTPUT:

        Nothing is returned but the internal state is modified.

        EXAMPLES::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.LLL(algorithm="wrapper")
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]
            sage: L.is_LLL_reduced()
            True
            sage: L.hermite_form() == A.hermite_form()
            True

            sage: set_random_seed(0)
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.LLL(algorithm="proved")
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced()
            True
            sage: L.hermite_form() == A.hermite_form()
            True

            sage: A = random_matrix(ZZ,10,10,x=-(10^5),y=10^5)
            sage: f = FP_LLL(A)
            sage: f.LLL(algorithm="fast")
            sage: L = f._sage_()
            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True

            sage: set_random_seed(0)
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.LLL(algorithm="fast", early_red=True)
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True

            sage: set_random_seed(0)
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.LLL(algorithm="heuristic")
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True

            sage: set_random_seed(0)
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.LLL(algorithm="heuristic", early_red=True)
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        _check_delta(delta)
        _check_eta(eta)
        _check_precision(precision)

        cdef LLLMethod method_
        if algorithm == "wrapper" or algorithm is None:
            method_ = LM_WRAPPER
        elif algorithm == "proved":
            method_ = LM_PROVED
        elif algorithm == "heuristic":
            method_ = LM_HEURISTIC
        elif algorithm == "fast":
            method_ = LM_FAST
        else:
            raise ValueError("algorithm '{}' unknown".format(algorithm))

        cdef int flags = LLL_DEFAULT

        if verbose:
            flags |= LLL_VERBOSE
        if early_red:
            flags |= LLL_EARLY_RED
        if siegel:
            flags |= LLL_SIEGEL

        if float_type is None and method_ == LM_FAST:
            float_type = "double"

        if method_ == LM_WRAPPER and check_float_type(float_type) != FT_DEFAULT:
            raise ValueError("fpLLL's LLL wrapper function requires "
                             "float type None")
        if method_ == LM_FAST and \
                check_float_type(float_type) not in (FT_DOUBLE, FT_LONG_DOUBLE):
            raise ValueError("fpLLL's LLL fast function requires "
                             "float type 'double' or 'long double'")

        sig_on()
        cdef int r = lllReduction(self._lattice[0], delta, eta, method_,
                                  check_float_type(float_type), precision, flags)
        sig_off()
        if r:
            raise RuntimeError( str(getRedStatusStr(r)) )

    def BKZ(self, int block_size, double delta=LLL_DEF_DELTA,
            float_type=None, int precision=0, int max_loops=0, int max_time=0,
            verbose=False, no_lll=False, bounded_lll=False, auto_abort=False):
        r"""
        Run BKZ reduction.

        INPUT:

        - ``block_size`` -- an integer from 1 to ``nrows``
        - ``delta`` -- (default: ``0.99``) LLL parameter `0.25 < \delta < 1.0`
        - ``float_type`` -- (default: ``None``) can be one of the following:

          * ``None`` - for automatic choice
          * ``'double'``
          * ``'long double'``
          * ``'dpe'``
          * ``'mpfr'``

        - ``verbose`` -- (default: ``False``) be verbose
        - ``no_lll`` -- (default: ``False``) to use LLL
        - ``bounded_lll`` -- (default: ``False``) bounded LLL
        - ``precision`` -- (default: ``0`` for automatic choice) bit
          precision to use if ``fp`` is ``'rr'``
        - ``max_loops`` -- (default: ``0`` for no restriction) maximum
          number of full loops
        - ``max_time`` -- (default: ``0`` for no restricion) stop after
          time seconds (up to loop completion)
        - ``auto_abort`` -- (default: ``False``) heuristic, stop when the
          average slope of `\log(\lVert b_i^* \rVert)` does not decrease
          fast enough

        OUTPUT:

        Nothing is returned but the internal state is modified.

        EXAMPLES::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=60, q=2^90, seed=42)
            sage: F = FP_LLL(A)

            sage: F.LLL()
            sage: F._sage_()[0].norm().n()
            7.810...

            sage: F.BKZ(10)
            sage: F._sage_()[0].norm().n()
            6.164...


            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=60, q=2^90, seed=42)
            sage: F = FP_LLL(A)
            sage: F.BKZ(10, max_loops=10, verbose=True)
            Entering BKZ:
            ...
            loops limit exceeded in BKZ
            sage: F._sage_()[0].norm().n()
            6.480...
        """
        if block_size <= 0:
            raise ValueError("block size must be > 0")
        if max_loops < 0:
            raise ValueError("maximum number of loops must be >= 0")
        if max_time < 0:
            raise ValueError("maximum time must be >= 0")

        _check_delta(delta)
        _check_precision(precision)

        cdef BKZParam o

        o.delta = delta
        o.blockSize = block_size
        o.flags = BKZ_DEFAULT

        if verbose:
            o.flags |= BKZ_VERBOSE
        if no_lll:
            o.flags |= BKZ_NO_LLL
        if bounded_lll:
            o.flags |= BKZ_BOUNDED_LLL
        if auto_abort:
            o.flags |= BKZ_AUTO_ABORT
        if max_loops:
            o.flags |= BKZ_MAX_LOOPS
            o.maxLoops = max_loops
        if max_time:
            o.flags |= BKZ_MAX_TIME
            o.maxTime = max_time

        sig_on()
        cdef int r = bkzReduction(self._lattice, NULL, o, check_float_type(float_type), precision)
        sig_off()
        if r:
            if r in (RED_BKZ_LOOPS_LIMIT, RED_BKZ_TIME_LIMIT):
                if verbose:
                    print(str(getRedStatusStr(r)))
            else:
                raise RuntimeError(str(getRedStatusStr(r)))


    def HKZ(self):
        """
        Run HKZ reduction.

        OUTPUT:

        Nothing is returned but the internal state is modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=10, q=2^60, seed=42)
            sage: F = FP_LLL(A)
            sage: F.HKZ()
            sage: F._sage_()
            [ -8  27   7  19  10  -5  14  34   4 -18]
            [-22  23   3 -14  11  30 -12  26  17  26]
            [-20   6 -18  33 -26  16   8 -15 -14 -26]
            [ -2  30   9 -30 -28 -19  -7 -28  12 -15]
            [-16   1  25 -23 -11 -21 -39   4 -34 -13]
            [-27  -2 -24 -67  32 -13  -6   0  15  -4]
            [  9 -12   7  31  22  -7 -63  11  27  36]
            [ 14  -4   0 -21 -17  -7  -9  35  79 -22]
            [-17 -16  54  21   0 -17  28 -45  -6  12]
            [ 43  16   6  30  24  17 -39 -46 -18 -22]

        """
        sig_on()
        cdef int r = hkzReduction(self._lattice[0])
        sig_off()
        if r:
            raise RuntimeError( str(getRedStatusStr(r)) )

    @rename_keyword(deprecation=19572, method='algorithm')
    def shortest_vector(self, algorithm=None):
        """
        Return a shortest vector.

        INPUT:

        - ``algorithm`` - (default: ``"proved"``) ``"proved"`` or ``"fast"``

        OUTPUT:

        A shortest non-zero vector for this lattice.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=40, q=2^60, seed=42)
            sage: F = FP_LLL(A)
            sage: F.shortest_vector('proved')  == F.shortest_vector('fast')
            True

        """
        cdef SVPMethod method_
        if algorithm == "proved" or algorithm is None:
            method_ = SVPM_PROVED
        elif algorithm == "fast":
            method_ = SVPM_FAST
        else:
            raise ValueError("algorithm '{}' unknown".format(algorithm))

        cdef int r

        cdef ZZ_mat[mpz_t] u

        r = lllReduction(self._lattice[0], u, LLL_DEF_DELTA, LLL_DEF_ETA,
                         LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT)
        if r:
            raise RuntimeError( str(getRedStatusStr(r)) )

        cdef vector[Z_NR[mpz_t]] solCoord
        cdef vector[Z_NR[mpz_t]] solution

        sig_on()
        r = shortestVector(self._lattice[0], solCoord, method_)
        sig_off()
        if r:
            raise RuntimeError("fpLLL's SVP solver returned an error ({:d})".format(r))

        vectMatrixProduct(solution, solCoord, self._lattice[0])

        VS = ZZ**self._lattice.getCols()
        cdef Vector_integer_dense v = Vector_integer_dense(VS, 0)

        for i in range(len(v)):
            mpz_set(v._entries[i], solution[i].getData())
        return v

    #
    # Deprecated interfaces
    #

    def wrapper(self, int precision=0, float eta=0.51, float delta=0.99):
        r"""
        Perform LLL reduction using fpLLL's \code{wrapper}
        implementation. This implementation invokes a sequence of
        floating point LLL computations such that

        * the computation is reasonably fast (based on an heuristic model)
        * the result is proven to be LLL reduced.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{\delta}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`

        OUTPUT:

        Nothing is returned but the internal state is modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.wrapper()
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]
            sage: L.is_LLL_reduced()
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)
        cdef int ret = 0

        cdef wrapper *w = wrapper_new(self._lattice, 0, eta, delta)
        sig_on()
        ret = w.LLL()
        sig_off()
        wrapper_delete(w)
        if ret != 0:
            raise RuntimeError("fpLLL returned {:d} != 0".format(ret))


    def proved(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's ``proved``
        implementation. This implementation is the only provable
        correct floating point implementation in the free software
        world. Provability is only guaranteed if the ``'mpfr'``
        implementation is chosen.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{δ}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`
        - ``implementation`` -- (default: ``"mpfr"``) which floating point
          implementation to use, can be one of the following:

          * ``"double"``
          * ``"dpe"``
          * ``"mpfr"``

        OUTPUT:

        Nothing is returned but the internal state modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.proved()
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced()
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)

        cdef proved_double *pdouble
        cdef proved_mpfr *pmpfr
        cdef proved_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = proved_double_new(self._lattice, precision, eta, delta)
           sig_on()
           ret = pdouble.LLL()
           sig_off()
           proved_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = proved_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           sig_on()
           ret = pdpe.LLL()
           sig_off()
           proved_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = proved_mpfr_new(self._lattice, precision, eta, delta)
           sig_on()
           ret = pmpfr.LLL()
           sig_off()
           proved_mpfr_delete(pmpfr)

        if ret != 0:
            raise RuntimeError("fpLLL returned {:d} != 0".format(ret))

    def fast(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        r"""
        Perform LLL reduction using fpLLL's fast
        implementation.  This implementation is the fastest floating
        point implementation currently available in the free software world.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{\delta}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`
        - ``implementation`` -- ignored

        OUTPUT:

        Nothing is returned but the internal is state modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10,x=-(10^5),y=10^5)
            sage: f = FP_LLL(A)
            sage: f.fast()
            sage: L = f._sage_()
            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)

        cdef fast_double *pdouble
        cdef int ret = 0

        pdouble = fast_double_new(self._lattice, precision, eta, delta)
        sig_on()
        ret = pdouble.LLL()
        sig_off()
        fast_double_delete(pdouble)

    def fast_early_red(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        r"""
        Perform LLL reduction using fpLLL's fast
        implementation with early reduction.

        This implementation inserts some early reduction steps inside
        the execution of the 'fast' LLL algorithm. This sometimes makes the
        entries of the basis smaller very quickly. It occurs in
        particular for lattice bases built from minimal polynomial or
        integer relation detection problems.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{\delta}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`
        - ``implementation`` -- (default: ``"mpfr"``) which floating point
          implementation to use, can be one of the following:

          * ``"double"``
          * ``"dpe"``
          * ``"mpfr"``

        OUTPUT:

        Nothing is returned but the internal state modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.fast_early_red()
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True
         """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)
        cdef int ret = 0

        cdef fast_early_red_double *pdouble

        pdouble = fast_early_red_double_new(self._lattice, precision, eta, delta)
        sig_on()
        ret = pdouble.LLL()
        sig_off()
        fast_early_red_double_delete(pdouble)

    def heuristic(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        r"""
        Perform LLL reduction using fpLLL's heuristic
        implementation.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{\delta}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`
        - ``implementation`` -- (default: ``"mpfr"``) which floating point
          implementation to use, can be one of the following:

          * ``"double"``
          * ``"dpe"``
          * ``"mpfr"``

        OUTPUT:

        Nothing is returned but the internal state modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.heuristic()
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)

        cdef heuristic_double *pdouble
        cdef heuristic_mpfr *pmpfr
        cdef heuristic_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
            pdouble = heuristic_double_new(self._lattice, precision, eta, delta)
            sig_on()
            ret = pdouble.LLL()
            sig_off()
            heuristic_double_delete(pdouble)
        elif implementation == "dpe":
            pdpe = heuristic_dpe_new(self._lattice, precision, <double>eta, <double>delta)
            sig_on()
            ret = pdpe.LLL()
            sig_off()
            heuristic_dpe_delete(pdpe)
        elif implementation == "mpfr":
            pmpfr = heuristic_mpfr_new(self._lattice, precision, eta, delta)
            sig_on()
            ret = pmpfr.LLL()
            sig_off()
            heuristic_mpfr_delete(pmpfr)

    def heuristic_early_red(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        r"""
        Perform LLL reduction using fpLLL's heuristic
        implementation with early reduction.

        This implementation inserts some early reduction steps inside
        the execution of the 'fast' LLL algorithm. This sometimes
        makes the entries of the basis smaller very quickly. It occurs
        in particular for lattice bases built from minimal polynomial
        or integer relation detection problems.

        INPUT:

        - ``precision`` -- (default: auto) internal precision
        - ``eta`` -- (default: ``0.51``) LLL parameter `\eta` with
          `1/2 \leq \eta < \sqrt{\delta}`
        - ``delta`` -- (default: ``0.99``) LLL parameter `\delta` with
          `1/4 < \delta \leq 1`
        - ``implementation`` -- (default: ``"mpfr"``) which floating point
          implementation to use, can be one of the following:

          * ``"double"``
          * ``"dpe"``
          * ``"mpfr"``

        OUTPUT:

        Nothing is returned but the internal state modified.

        EXAMPLE::

            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10); A
            [   -8     2     0     0     1    -1     2     1   -95    -1]
            [   -2   -12     0     0     1    -1     1    -1    -2    -1]
            [    4    -4    -6     5     0     0    -2     0     1    -4]
            [   -6     1    -1     1     1    -1     1    -1    -3     1]
            [    1     0     0    -3     2    -2     0    -2     1     0]
            [   -1     1     0     0     1    -1     4    -1     1    -1]
            [   14     1    -5     4    -1     0     2     4     1     1]
            [   -2    -1     0     4    -3     1    -5     0    -2    -1]
            [   -9    -1    -1     3     2     1    -1     1    -2     1]
            [   -1     2    -7     1     0     2     3 -1955   -22    -1]

            sage: F = FP_LLL(A)
            sage: F.heuristic_early_red()
            sage: L = F._sage_(); L
            [   1    0    0   -3    2   -2    0   -2    1    0]
            [  -1    1    0    0    1   -1    4   -1    1   -1]
            [  -2    0    0    1    0   -2   -1   -3    0   -2]
            [  -2   -2    0   -1    3    0   -2    0    2    0]
            [   1    1    1    2    3   -2   -2    0    3    1]
            [  -4    1   -1    0    1    1    2    2   -3    3]
            [   1   -3   -7    2    3   -1    0    0   -1   -1]
            [   1   -9    1    3    1   -3    1   -1   -1    0]
            [   8    5   19    3   27    6   -3    8  -25  -22]
            [ 172  -25   57  248  261  793   76 -839  -41  376]

            sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
            True
            sage: L.hermite_form() == A.hermite_form()
            True
        """
        deprecation(15976, 'You can just call LLL() instead')

        _check_precision(precision)
        _check_eta(eta)
        _check_delta(delta)

        cdef heuristic_early_red_double *pdouble
        cdef heuristic_early_red_mpfr *pmpfr
        cdef heuristic_early_red_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = heuristic_early_red_double_new(self._lattice, precision, eta, delta)
           sig_on()
           ret = pdouble.LLL()
           sig_off()
           heuristic_early_red_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = heuristic_early_red_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           sig_on()
           ret = pdpe.LLL()
           sig_off()
           heuristic_early_red_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = heuristic_early_red_mpfr_new(self._lattice, precision, eta, delta)
           sig_on()
           ret = pmpfr.LLL()
           sig_off()
           heuristic_early_red_mpfr_delete(pmpfr)

#
# Lattice basis generators
#

def gen_intrel(int d, int b):
    r"""
    Return a `(d+1 \times d)`-dimensional knapsack-type random lattice,
    where the `x_i`'s are random ``b`` bits integers.

    INPUT:

    - ``d`` -- dimension
    - ``b`` -- bitsize of entries

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_intrel
        sage: A = gen_intrel(10,10); A
        [116   1   0   0   0   0   0   0   0   0   0]
        [331   0   1   0   0   0   0   0   0   0   0]
        [303   0   0   1   0   0   0   0   0   0   0]
        [963   0   0   0   1   0   0   0   0   0   0]
        [456   0   0   0   0   1   0   0   0   0   0]
        [225   0   0   0   0   0   1   0   0   0   0]
        [827   0   0   0   0   0   0   1   0   0   0]
        [381   0   0   0   0   0   0   0   1   0   0]
        [ 99   0   0   0   0   0   0   0   0   1   0]
        [649   0   0   0   0   0   0   0   0   0   1]

        sage: L = A.LLL(); L
        [ 1  1  1  0  0  0  0 -1  1  0  0]
        [ 1  0  1  0  0 -1  1  0  0 -1  0]
        [ 0  0  1  1  0 -1  0 -1  0  0  1]
        [ 0 -1  0 -1 -1  1  0  1  0  1  0]
        [-1 -1  0 -1  0 -1  1  0  0  0  1]
        [ 0  1 -1  0  0 -1  1  1 -1  0  0]
        [ 0  0  0  0 -1  1  1  0  1 -1  0]
        [ 1 -1 -1  0  0 -1 -1  0  1  1  1]
        [-1  0  0 -1 -1  0 -1  1  2 -1  0]
        [-1 -1  0  0  1  0  2  0  0  0 -2]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](d,d+1)
    A.gen_intrel(b)

    B = to_sage(A)
    del A
    return B

def gen_simdioph(int d, int b, int b2):
    """
    Return a ``d``-dimensional simultaneous diophantine approximation random
    lattice, where the ``d`` `x_i`'s are random ``b`` bits integers.

    INPUT:

    - ``d`` -- dimension
    - ``b`` -- bitsize of entries
    - ``b2`` -- bitsize of entries

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_simdioph
        sage: A = gen_simdioph(10,10,3); A
        [   8  395  975  566  213  694  254  629  303  597]
        [   0 1024    0    0    0    0    0    0    0    0]
        [   0    0 1024    0    0    0    0    0    0    0]
        [   0    0    0 1024    0    0    0    0    0    0]
        [   0    0    0    0 1024    0    0    0    0    0]
        [   0    0    0    0    0 1024    0    0    0    0]
        [   0    0    0    0    0    0 1024    0    0    0]
        [   0    0    0    0    0    0    0 1024    0    0]
        [   0    0    0    0    0    0    0    0 1024    0]
        [   0    0    0    0    0    0    0    0    0 1024]

        sage: L = A.LLL(); L
        [ 192  264 -152  272   -8  272  -48 -264  104   -8]
        [-128 -176 -240  160 -336  160   32  176  272 -336]
        [ -24 -161  147  350  385  -34  262  161  115  257]
        [ 520   75 -113  -74 -491   54  126  -75  239 -107]
        [-376 -133  255   22  229  150  350  133   95 -411]
        [-168 -103    5  402 -377 -238 -214  103 -219 -249]
        [-352   28  108 -328 -156  184   88  -28  -20  356]
        [ 120 -219  289  298  123  170 -286  219  449 -261]
        [ 160 -292   44   56  164  568  -40  292  -84 -348]
        [-192  760  152 -272    8 -272   48  264 -104    8]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](d,d)
    A.gen_simdioph(b, b2)

    B = to_sage(A)
    del A
    return B

def gen_uniform(int nr, int nc, int b):
    r"""
    Return a `(nr \times nc)` matrix where the entries are random ``b``
    bits integers.

    INPUT:

    - ``nr`` -- row dimension
    - ``nc`` -- column dimension
    - ``b`` -- bitsize of entries

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_uniform
        sage: A = gen_uniform(10,10,12); A
        [ 980 3534  533 3303 2491 2960 1475 3998  105  162]
        [1766 3683 2782  668 2356 2149 1887 2327  976 1151]
        [1573  438 1480  887 1490  634 3820 3379 4074 2669]
        [ 215 2054 2388 3214 2459  250 2921 1395 3626  434]
        [ 638 4011 3626 1864  633 1072 3651 2339 2926 1004]
        [3731  439 1087 1088 2627 3446 2669 1419  563 2079]
        [1868 3196 3712 4016 1451 2589 3327  712  647 1057]
        [2068 2761 3479 2552  197 1258 1544 1116 3090 3667]
        [1394  529 1683 1781 1779 3032   80 2712  639 3047]
        [3695 3888 3139  851 2111 3375  208 3766 3925 1465]

        sage: L = A.LLL(); L
        [  200 -1144  -365   755  1404  -218  -937   321  -718   790]
        [  623   813   873  -595  -422   604  -207  1265 -1418  1360]
        [ -928  -816   479  1951  -319 -1295   827   333  1232   643]
        [-1802 -1904  -952   425  -141   697   300  1608  -501  -767]
        [ -572 -2010  -734   358 -1981  1101  -870    64   381  1106]
        [  853  -223   767  1382  -529  -780  -500  1507 -2455 -1190]
        [-1016 -1755  1297 -2210  -276  -114   712   -63   370   222]
        [ -430  1471   339  -513  1361  2715  2076  -646 -1406   -60]
        [-3390   748    62   775   935  1697  -306  -618    88  -452]
        [  713 -1115  1887  -563   733  2443   816   972   876 -2074]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](nr,nc)
    A.gen_uniform(b)

    B = to_sage(A)
    del A
    return B

def gen_ntrulike(int d, int b, int q):
    r"""
    Generate a NTRU-like lattice of dimension `(2d \times 2d)`, with the
    coefficients `h_i` chosen as random `b` bits integers and parameter `q`::

        [[ 1 0 ... 0 h0      h1 ... h_{d-1} ]
         [ 0 1 ... 0 h1      h2 ... h0      ]
         [ ................................ ]
         [ 0 0 ... 1 h_{d-1} h0 ... h_{d-1} ]
         [ 0 0 ... 0 q       0  ...  0      ]
         [ 0 0 ... 0 0       q  ...  0      ]
         [ ................................ ]
         [ 0 0 ... 0 0       0  ...  q      ]]

    INPUT:

    - ``d`` -- dimension
    - ``b`` -- bitsize of entries
    - ``q`` -- the `q` above

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_ntrulike
        sage: A = gen_ntrulike(5,10,12); A
        [  1   0   0   0   0 320 351 920 714  66]
        [  0   1   0   0   0 351 920 714  66 320]
        [  0   0   1   0   0 920 714  66 320 351]
        [  0   0   0   1   0 714  66 320 351 920]
        [  0   0   0   0   1  66 320 351 920 714]
        [  0   0   0   0   0  12   0   0   0   0]
        [  0   0   0   0   0   0  12   0   0   0]
        [  0   0   0   0   0   0   0  12   0   0]
        [  0   0   0   0   0   0   0   0  12   0]
        [  0   0   0   0   0   0   0   0   0  12]

        sage: L = A.LLL(); L
        [-1 -1  0  0  0  1  1 -2  0 -2]
        [-1  0  0  0 -1 -2  1  1 -2  0]
        [ 0 -1 -1  0  0  1 -2  0 -2  1]
        [ 0  0  1  1  0  2  0  2 -1 -1]
        [ 0  0  0  1  1  0  2 -1 -1  2]
        [-2 -1 -2  1  1  1  0  1  1  0]
        [-1 -2  1  1 -2  0  1  0  1  1]
        [ 2 -1 -1  2  1 -1  0 -1  0 -1]
        [-1 -1  2  1  2 -1 -1  0 -1  0]
        [ 1 -2 -1 -2  1  0  1  1  0  1]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](2*d,2*d)
    A.gen_ntrulike(b, q)

    B = to_sage(A)
    del A
    return B

def gen_ntrulike2(int d, int b, int q):
    """
    Like :func:`gen_ntrulike` but with the `q` vectors coming first.

    INPUT:

    - ``d`` -- dimension
    - ``b`` -- bitsize of entries
    - ``q`` -- see :func:`gen_ntrulike`

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_ntrulike2
        sage: A = gen_ntrulike2(5,10,12); A
        [ 12   0   0   0   0   0   0   0   0   0]
        [  0  12   0   0   0   0   0   0   0   0]
        [  0   0  12   0   0   0   0   0   0   0]
        [  0   0   0  12   0   0   0   0   0   0]
        [  0   0   0   0  12   0   0   0   0   0]
        [902 947 306  40 908   1   0   0   0   0]
        [947 306  40 908 902   0   1   0   0   0]
        [306  40 908 902 947   0   0   1   0   0]
        [ 40 908 902 947 306   0   0   0   1   0]
        [908 902 947 306  40   0   0   0   0   1]

        sage: L = A.LLL(); L
        [ 1  0  0  2 -3 -2  1  1  0  0]
        [-1  0 -2  1  2  2  1 -2 -1  0]
        [ 0  2 -1 -2  1  0 -2 -1  2  1]
        [ 0  3  0  1  3  1  0 -1  1  0]
        [ 2 -1  0 -2  1  1 -2 -1  0  2]
        [ 0 -1  0 -1 -1  1  4 -1 -1  0]
        [ 2  1  1  1 -1 -3 -2 -1 -1 -1]
        [-1  0 -1  0 -1  4 -1 -1  0  1]
        [ 0  1 -2  1  1 -1  0  1 -3 -2]
        [-2  1  1  0  1 -3 -2 -1  0  1]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](2*d,2*d)
    A.gen_ntrulike2(b,q)

    B = to_sage(A)
    del A
    return B

def gen_ajtai(int d, float alpha):
    r"""
    Return Ajtai-like `(d \times d)`-matrix of floating point parameter
    `\alpha`. The matrix is lower-triangular, `B_{imi}` is
    `~2^{ (d-i+1)^{\alpha} }` and `B_{i,j}` is `~B_{j,j} / 2` for `j < i`.

    INPUT:

    - ``d`` -- dimension
    - ``alpha`` -- the `\alpha` above

    OUTPUT:

    An integer lattice.

    EXAMPLE::

        sage: from sage.libs.fplll.fplll import gen_ajtai
        sage: A = gen_ajtai(10, 0.7); A # random output
        [117   0   0   0   0   0   0   0   0   0]
        [ 11  55   0   0   0   0   0   0   0   0]
        [-47  21 104   0   0   0   0   0   0   0]
        [ -3 -22 -16  95   0   0   0   0   0   0]
        [ -8 -21  -3 -28  55   0   0   0   0   0]
        [-33 -15 -30  37   8  52   0   0   0   0]
        [-35  21  41 -31 -23  10  21   0   0   0]
        [ -9  20 -34 -23 -18 -13  -9  63   0   0]
        [-11  14 -38 -16 -26 -23  -3  11   9   0]
        [ 15  21  35  37  12   6  -2  10   1  17]

        sage: L = A.LLL(); L # random output
        [  4   7  -3  21 -14 -17  -1  -1  -8  17]
        [-20   0  -6   6 -11  -4 -19  10   1  17]
        [-22  -1   8 -21  18 -29   3  11   9   0]
        [ 31   8  20   2 -12  -4 -27 -22 -18   0]
        [ -2   6  -4   7  -8 -10   6  52  -9   0]
        [  3  -7  35 -12 -29  23   3  11   9   0]
        [-16  -6 -16  37   2  11  -1  -9   7 -34]
        [ 11  55   0   0   0   0   0   0   0   0]
        [ 11  14  38  16  26  23   3  11   9   0]
        [ 13 -28  -1   7 -11  11 -12   3  54   0]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

    """
    cdef ZZ_mat[mpz_t] *A = new ZZ_mat[mpz_t](d,d)
    A.gen_ajtai(alpha)

    B = to_sage(A)
    del A
    return B

cdef to_sage(ZZ_mat[mpz_t] *A):
    """
    Return a Sage integer matrix for ``A``. ``A`` is not destroyed.

    INPUT:

    - ``A`` -- ZZ_mat
    """
    cdef int i,j

    cdef Matrix_integer_dense B = <Matrix_integer_dense>matrix(ZZ, A.getRows(), A.getCols())

    for i from 0 <= i < A.getRows():
        for j from 0 <= j < A.getCols():
            B.set_unsafe_mpz(i,j,A[0][i][j].getData())
    return B
