r"""
Dense Matrices over `\mathbb F_q`, with `q<255`

This module is a wrapper for version 2.4.24 of the Aachen
`C-MeatAxe <http://www.math.rwth-aachen.de/homes/MTX/download.html>`_,
improved by an implementation of the Winograd-Strassen multiplication
algorithm. It provides matrices over the finite field `\mathbb F_q`,
where `q\le 255`.

By default, it is only used when `q` is odd and not prime, because other
matrix implementations in SageMath perform better for prime fields or in
characteristic two.

AUTHORS:

- Simon King (2015-09): initial version

"""

#*****************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

## Define an environment variable that enables MeatAxe to find
## its multiplication tables.

from sage.env import DOT_SAGE
import os
cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *s, Py_ssize_t len)
    char* PyString_AsString(object string)
MtxLibDir = PyString_AsString(os.path.join(DOT_SAGE,'meataxe'))

####################
#
# import sage types
#
####################

from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod import IntegerMod_int
from sage.matrix.constructor import random_matrix
from sage.rings.arith import is_prime_power, factor
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.randstate import current_randstate
from sage.misc.cachefunc import cached_method, cached_function
from sage.structure.element cimport Element, ModuleElement, RingElement, Matrix

from libc.stdlib cimport free
from libc.string cimport memset, memcpy

cimport sage.matrix.matrix0

include "cysignals/memory.pxi"

####################
#
# auxiliary functions
#
####################
import sys
from libc.string cimport memcpy

# Fast conversion from field to int and int to field
cdef class FieldConverter_class:
    """
    An auxiliary class, used to convert between <int> and finite field element

    This class is for non-prime fields only. The method
    :meth:`int_to_field` exists for speed. The method
    :meth:`field_to_int` exists in order to have a common interface
    for elements of prime and non-prime fields; see
    :class:`PrimeFieldConverter_class`.

    EXAMPLE::

        sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class  # optional: meataxe
        sage: F.<y> = GF(125)
        sage: C = FieldConverter_class(F)               # optional: meataxe
        sage: C.int_to_field(15)                        # optional: meataxe
        3*y
        sage: F.fetch_int(15)                           # optional: meataxe
        3*y
        sage: %timeit C.int_to_field(15)    # not tested
        625 loops, best of 3: 1.04 µs per loop
        sage: %timeit F.fetch_int(15)       # not tested
        625 loops, best of 3: 3.97 µs per loop
        sage: C.field_to_int(y)                         # optional: meataxe
        5
        sage: y.integer_representation()
        5

    """
    def __init__(self, field):
        """
        INPUT:

        A finite *non-prime* field. This assumption is not tested.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class # optional: meataxe
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)           # optional: meataxe
            sage: C.int_to_field(15)                    # optional: meataxe
            3*y
            sage: F.fetch_int(15)
            3*y
            sage: C.field_to_int(y)                     # optional: meataxe
            5
            sage: y.integer_representation()
            5

        """
        self.field = field._cache.fetch_int
    cpdef object int_to_field(self, int x):
        """
        Fetch a python int into the field.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class  # optional: meataxe
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)           # optional: meataxe
            sage: C.int_to_field(15)                    # optional: meataxe
            3*y
            sage: F.fetch_int(15)
            3*y

        """
        return self.field(x)
    cpdef int field_to_int(self, x):
        """
        Represent a field element by a python int.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class  # optional: meataxe
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)           # optional: meataxe
            sage: C.field_to_int(y)                     # optional: meataxe
            5
            sage: y.integer_representation()
            5

        """
        return x.integer_representation()

cdef class PrimeFieldConverter_class(FieldConverter_class):
    """
    An auxiliary class, used to convert between <int> and finite field element

    This class is for prime fields only. The methods
    :meth:`int_to_field` and :meth:`field_to_int` exist in order to
    have a common interface for elements of prime and non-prime fields;
    see :class:`FieldConverter_class`.

    EXAMPLE::

        sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class # optional: meataxe
        sage: F = GF(5)
        sage: C = PrimeFieldConverter_class(F)      # optional: meataxe
        sage: C.int_to_field(int(2))                # optional: meataxe
        2
        sage: F(2)
        2
        sage: C.field_to_int(F(2))                  # optional: meataxe
        2
        sage: int(F(2))
        2

    """
    def __init__(self, field):
        """
        INPUT:

        A finite *prime* field. This assumption is not tested.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class  # optional: meataxe
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)  # optional: meataxe
            sage: C.int_to_field(int(2))            # optional: meataxe
            2
            sage: F(2)
            2
            sage: C.field_to_int(F(2))              # optional: meataxe
            2
            sage: int(F(2))
            2

        """
        self.field = field
    cpdef object int_to_field(self, int x):
        """
        Fetch a python int into the field.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class  # optional: meataxe
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)  # optional: meataxe
            sage: C.int_to_field(int(2))            # optional: meataxe
            2
            sage: F(2)
            2

        """
        return IntegerMod_int(self.field, x)
    cpdef int field_to_int(self, x):
        """
        Represent a field element by a python int.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class  # optional: meataxe
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)      # optional: meataxe
            sage: C.field_to_int(F(2))                  # optional: meataxe
            2
            sage: int(F(2))
            2

        """
        return int(x)

cdef dict _converter_cache = {}
cdef FieldConverter_class FieldConverter(field):
    """
    Return a :class:`FieldConverter_class` or :class:`PrimeFieldConverter_class` instance,
    depending whether the field is prime or not.

    EXAMPLE::

        sage: MS = MatrixSpace(GF(5^3,'y'),2)
        sage: A = MS.random_element()
        sage: A*2 == A+A    # indirect doctest
        True
        sage: A = MS.random_element()
        sage: A*2 == A+A
        True

    """
    try:
        return _converter_cache[field]
    except KeyError:
        if field.is_prime_field():
            return _converter_cache.setdefault(field, PrimeFieldConverter_class(field))
        return _converter_cache.setdefault(field, FieldConverter_class(field))

######################################
## Error handling for MeatAxe, to prevent immediate exit of the program

cdef dict ErrMsg = {
    "Not enough memory": MemoryError,
    "Time limit exceeded": RuntimeError,
    "Division by zero": ZeroDivisionError,
    "Bad file format": IOError,
    "Bad argument": ValueError,
    "Argument out of range": IndexError,

    "Matrix not in echelon form": ValueError,
    "Matrix not square": ArithmeticError,
    "Incompatible objects": TypeError,

    "Bad syntax, try `-help'": SyntaxError,
    "Bad usage of option, try `-help'": ValueError,
    "Bad number of arguments, try `-help'": ValueError,

    "Not a matrix": TypeError,
    "Not a permutation": TypeError
}

from cpython.exc cimport PyErr_SetObject

cdef void ErrorHandler(MtxErrorRecord_t *err):
    PyErr_SetObject(ErrMsg.get(err.Text, SystemError), "{} in file {} (line {})".format(err.Text, err.FileInfo.BaseName, err.LineNo))

MtxSetErrorHandler(ErrorHandler)

######################################
##
## Wrapper for MeatAxe matrices
##
######################################

cdef class Matrix_gfpn_dense(Matrix_dense):
    r"""
    Dense matrices over `\mathbb F_q`, `q<255` odd and not prime.

    NOTE:

    This class uses a major modification of the Aachen C-MeatAxe
    as backend. In principle, it would also work for prime fields
    and in characteristic two. However, other matrices in Sage,
    relying on linbox, m4ri or m4rie, are more efficient in these
    cases.

    EXAMPLES::

        sage: M = MatrixSpace(GF(25,'z'),2,3)([1,2,3,4,5,6])
        sage: print(M)
        [1 2 3]
        [4 0 1]
        sage: type(M)     # optional: meataxe
        <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>

    The documentation of the ``__init__`` methods shows further
    ways of creating a :class:`Matrix_gfpn_dense` instance.
    However, these should only be of internal use.

    """
##################
## Init, Dealloc, Copy
    def __cinit__(self, parent=None, entries=None, *args, **kwds):
        """
        TESTS::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense  # optional: meataxe
            sage: Matrix_gfpn_dense.__new__(Matrix_gfpn_dense)   # optional: meataxe
            []
            sage: Matrix_gfpn_dense(MatrixSpace(GF(64,'z'),4), None)  # optional: meataxe
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        """
        if parent is None:  # this makes Matrix_gfpn_dense.__new__(Matrix_gfpn_dense) work,
                            # returning a non-initialised matrix
            return
        if isinstance(parent, basestring): # this allows to provide a file when initialising a matrix
            return
        cdef int f = parent.base_ring().order()
        cdef int nrows = parent.nrows()
        cdef int ncols = parent.ncols()
        self.Data = MatAlloc(f, nrows, ncols)

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense  # optional: meataxe
            sage: Matrix_gfpn_dense.__new__(Matrix_gfpn_dense)   # optional: meataxe
            []
            sage: M = None
            sage: M = Matrix_gfpn_dense(MatrixSpace(GF(64,'z'),4), None)  # optional: meataxe
            sage: del M    # indirect doctest
        """
        if self.Data != NULL:
            MatFree(self.Data)
            self.Data = NULL

    def __init__(self, parent, data=None, mutable=True, copy=False, coerce=False):
        """
        Matrix extension class using libmeataxe as backend

        INPUT:

        Instances of this class can be created by providing one of
        the following input data, where ``q<255`` is a prime power,
        ``m,n`` are non-negative integers, and `a_{11},...,a_{mn}`
        can be coerced into ``GF(q)``. Note that a user should
        create these instances via the matrix constructors; what
        we explain here is for internal use only!

        - None => empty matrix over an unspecified field (used for unpickling)
        - a string ``f`` ==> load matrix from the file named ``f``
        - A matrix space of `m\\times n` matrices over GF(q) and either

          - a list `[a_{11},a_{12},...,a_{1n},a_{21},...,a_{m1},...,a_{mn}]`,
            which results in a matrix with the given marks
          - ``None``, which is the fastest way to creata a zero matrix.
          - an element of GF(q), which results in a diagonal matrix with the
            given element on the diagonal.

        If the optional parameter ``mutable`` is ``False`` (by default,
        it is ``True``), the resulting matrix can not be changed, and
        it can be used as key in a Python dictionary.

        The arguments ``copy`` and ``coerce`` are ignored, they are only
        here for a common interface with :class:`~sage.matrix.matrix.Matrix`.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense  # optional: meataxe

        1. Creating an empty matrix::

            sage: Matrix_gfpn_dense(None)  # optional: meataxe
            []

        2. Creating a zero (3x2)-matrix::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(4,'z'),3,2))  # optional: meataxe
            [0 0]
            [0 0]
            [0 0]

        3. Creating a matrix from a list or list of lists::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[1,2,3,4,5,6])  # optional: meataxe
            [1 2 3]
            [4 0 1]
            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[[1,2,3],[4,5,6]])    # optional: meataxe
            [1 2 3]
            [4 0 1]

        4. Creating a diagonal matrix::

            sage: M = Matrix_gfpn_dense(MatrixSpace(GF(7),5),2); M  # optional: meataxe
            [2 0 0 0 0]
            [0 2 0 0 0]
            [0 0 2 0 0]
            [0 0 0 2 0]
            [0 0 0 0 2]

        5. Creating a matrix from a file in MeatAxe format.

           This is not tested.

        TESTS::

            sage: MS = MatrixSpace(GF(125,'y'),2)  # indirect doctest
            sage: A = MS(0)
            sage: A.left_kernel()
            Vector space of degree 2 and dimension 2 over Finite Field in y of size 5^3
            Basis matrix:
            [1 0]
            [0 1]
            sage: A.right_kernel()
            Vector space of degree 2 and dimension 2 over Finite Field in y of size 5^3
            Basis matrix:
            [1 0]
            [0 1]

        """
        if parent is None:
            self._is_immutable = False
            self._ncols = 0
            self._nrows = 0
            self._cache = {}
            return
        if isinstance(parent, basestring): # load from file
            FILE = os.path.realpath(parent)
            try:
                fsock = open(FILE,"rb",0)
                fsock.close()
            except (OSError,IOError):
                return
            self.Data = MatLoad(FILE)
            FfSetField(self.Data.Field)
            B = GF(self.Data.Field, 'z')
            parent = MatrixSpace(B, self.Data.Nor, self.Data.Noc)
            self._is_immutable = False
            self._parent = parent
            self._base_ring = B
            self._converter = FieldConverter(B)
            self._ncols = self.Data.Noc
            self._nrows = self.Data.Nor
            self._cache = {}
            return

        if not self.Data: # should have been initialised by __cinit__
            raise MemoryError("Error allocating memory for MeatAxe matrix")
        Matrix_dense.__init__(self, parent)
        self._is_immutable = not mutable
        B = self._base_ring
        self._converter = FieldConverter(B)
        if data is None:
            return

        cdef int i,j
        cdef FEL f
        cdef PTR x
        if not isinstance(data,list):
            if not data:
                return
            if self._nrows != self._ncols:
                raise ValueError("Cannot initialise non-square matrix from {}".format(data))
            f = FfFromInt(self._converter.field_to_int(self._coerce_element(data)))
            x = self.Data.Data
            for j from 0 <= j < self.Data.Noc:
                FfInsert(x,j,f)
                FfStepPtr(&x)
            return

        x = self.Data.Data
        cdef int nr = self.Data.Nor
        cdef int nc = self.Data.Noc
        assert self._ncols == nc
        assert self._nrows == nr
        if nr==0 or nc==0:
            return
        if len(data)<nr:
            raise ValueError("Expected a list of size at least the number of rows")
        cdef list dt, dt_i
        FfSetField(self.Data.Field)
        FfSetNoc(nc)
        if isinstance(data[0],list):
            # The matrix is given by a list of rows
            dt = data
            for i from 0 <= i < nr:
                idx = 0
                dt_i = dt[i]
                for j from 0 <= j < nc:
                    FfInsert(x, j, FfFromInt(self._converter.field_to_int(self._coerce_element(dt_i[j]))))
                FfStepPtr(&(x))
        else:
            # It is supposed to be a flat list of all entries, sorted by rows
            dtnext = data.__iter__().next
            for i from 0 <= i < nr:
                for j from 0 <= j < nc:
                    bla = self._converter.field_to_int(self._coerce_element(dtnext()))
                    FfInsert(x, j, FfFromInt(bla))
                FfStepPtr(&(x))

    cdef Matrix_gfpn_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        r"""
        Return a new matrix with no entries set.
        """
        cdef Matrix_gfpn_dense res
        res = self.__class__.__new__(self.__class__)

        if nrows == self._nrows and ncols == self._ncols:
            res._parent = self._parent
        else:
            res._parent = self.matrix_space(nrows, ncols)
        res._ncols  = ncols
        res._nrows  = nrows
        res._base_ring = self._base_ring
        res._converter = self._converter
        return res

    def __copy__(self):
        """
        Return a copy of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(GF(25,'x'), 3, 20)([20*[0],20*[0],[1]+19*[0]])
            sage: N = copy(M)   # indirect doctest
            sage: print(N)
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            sage: N== M
            True
            sage: N is M
            False
            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense  # optional: meataxe
            sage: M = Matrix_gfpn_dense('')   # optional: meataxe
            sage: N = copy(M)
            sage: N                         # optional: meataxe
            []
            sage: N == M
            True
            sage: N is M
            False
        """
        cdef Matrix_gfpn_dense retval = self._new(self._nrows, self._ncols)
        retval._is_immutable = False  # a copy of a matrix is mutable!
        retval._cache = dict(self._cache.iteritems()) if self._cache is not None else {}
        if self.Data:
            retval.Data = MatDup(self.Data)
        else:
            retval.Data = NULL
        return retval

    def __reduce__(self):
        """
        TESTS::

            sage: M = MatrixSpace(GF(9,'x'),10,10).random_element()
            sage: M == loads(dumps(M))   # indirect doctest
            True
            sage: M is loads(dumps(M))
            False
        """
        cdef char* d
        cdef int i,NR
        cdef PTR p
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
            return mtx_unpickle, (self._parent, self.Data.Nor, self.Data.Noc,
                        PyString_FromStringAndSize(<char*>self.Data.Data,self.Data.RowSize * self.Data.Nor),
                        not self._is_immutable) # for backward compatibility with the group cohomology package
        else:
            return mtx_unpickle, (0, 0, 0, '', not self._is_immutable)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get an element without checking.

        TEST::

            sage: F.<z> = GF(9)
            sage: M = MatrixSpace(F,3)(sorted(list(F)))
            sage: type(M)               # optional: meataxe
            <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: M                     # indirect doctest
            [      0       1       2]
            [      z   z + 1   z + 2]
            [    2*z 2*z + 1 2*z + 2]

        """
        if self.Data == NULL:
            raise IndexError("Matrix is empty")
        FfSetField(self.Data.Field)
        return self._converter.int_to_field(FfToInt(FfExtract(MatGetPtr(self.Data,i), j)))

    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j):
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        return FfToInt(FfExtract(MatGetPtr(self.Data,i), j))

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set values without bound checking.

        TESTS:

        The following test would have failed in a preliminary version
        of this MeatAxe wrapper::

            sage: K.<x> = GF(125)
            sage: M = MatrixSpace(K,9,9)()
            sage: N = MatrixSpace(GF(9,'x'),20).random_element()
            sage: M[2,2] = x
            sage: M
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 x 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]

        """
        # ASSUMPTION: value's parent is the base ring
        if self.Data == NULL:
            raise IndexError("Matrix is empty")
        FfSetField(self.Data.Field)
        FfInsert(MatGetPtr(self.Data,i), j, FfFromInt(self._converter.field_to_int(value)))

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        FfInsert(FfGetPtr(self.Data.Data,i), j, FfFromInt(value))

    def randomize(self, density=None, nonzero=False, *args, **kwds):
        """
        Fill the matrix with random values.

        INPUT:

        - ``density`` (optional real number between zero and one) --
          the expected density of the resulting matrix
        - ``nonzero`` (optional bool, default ``False``) --
          If true, all inserted marks are non-zero.

        EXAMPLE::

            sage: MS = MatrixSpace(GF(27,'z'),6,6)
            sage: M = MS.random_element()       # indirect doctest
            sage: M                             # optional: meataxe
            [              1           z + 1     z^2 + z + 1             z^2       2*z^2 + z           z + 1]
            [2*z^2 + 2*z + 2   2*z^2 + z + 2         z^2 + 1 2*z^2 + 2*z + 2         z^2 + z   2*z^2 + z + 1]
            [        2*z + 2     z^2 + z + 2           z + 2 2*z^2 + 2*z + 2           2*z^2           2*z^2]
            [  2*z^2 + z + 2             z^2           z + 2         z^2 + z       2*z^2 + 2         z^2 + 2]
            [      2*z^2 + z             2*z 2*z^2 + 2*z + 1       2*z^2 + 1 2*z^2 + 2*z + 1       2*z^2 + z]
            [        2*z + 1         z^2 + z             z^2             z^2     2*z^2 + 2*z           z + 1]
            sage: type(M)                           # optional: meataxe
            <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: MS.random_element(nonzero=True)   # optional: meataxe
            [            2*z               1   z^2 + 2*z + 1   2*z^2 + z + 1             z^2     z^2 + z + 1]
            [    2*z^2 + 2*z   2*z^2 + z + 2         2*z + 1       z^2 + 2*z     2*z^2 + 2*z             z^2]
            [        z^2 + z     z^2 + z + 2 2*z^2 + 2*z + 1         z^2 + 2               1           2*z^2]
            [              z     2*z^2 + 2*z           2*z^2         2*z + 1           z + 2           z + 2]
            [        z^2 + z             z^2           z + 2     2*z^2 + 2*z         2*z + 1         z^2 + z]
            [    z^2 + z + 2       2*z^2 + z             z^2           z + 1     2*z^2 + 2*z   z^2 + 2*z + 1]
            sage: MS.random_element(density=0.5)    # optional: meataxe
            [        z^2 + 2               0   z^2 + 2*z + 2       2*z^2 + z               0     z^2 + z + 2]
            [              0               1               0               0               0               0]
            [  2*z^2 + z + 1   2*z^2 + z + 2               0     z^2 + z + 2               0     z^2 + z + 1]
            [              0               0               0               0               0               0]
            [2*z^2 + 2*z + 2               0               0   2*z^2 + z + 2               0         2*z + 1]
            [              0       2*z^2 + z               0               1               0   2*z^2 + z + 1]

        """
        self.check_mutability()
        cdef int fl = self.Data.Field
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.clear_cache()

        cdef PTR x
        cdef unsigned char *y
        x = self.Data.Data
        cdef int nr = self.Data.Nor
        cdef int nc = self.Data.Noc
        cdef int i, j, k

        FfSetField(fl)
        FfSetNoc(nc)
        cdef int O, MPB, tmp
        randint = current_randstate().c_random
        randdouble = current_randstate().c_rand_double

        if not nonzero:
            if density == 1:
                MPB = 0
                tmp = fl
                while tmp <= 256:
                    MPB += 1
                    tmp *= fl
                O = (fl**MPB)
                sig_on()
                if nc%MPB:
                    for i from 0 <= i < nr:
                        y = <unsigned char*>x
                        for j from 0 <= j < FfCurrentRowSizeIo-1:
                            y[j] = randint()%O
                        y[FfCurrentRowSizeIo-1] = randint()%(fl**(nc%MPB))
                        FfStepPtr(&(x))
                else:
                    for i from 0 <= i < nr:
                        y = <unsigned char*>x
                        for j from 0 <= j < FfCurrentRowSizeIo:
                            y[j] = randint()%O
                        FfStepPtr(&(x))
                sig_off()
            else:
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        if randdouble() < density:
                            FfInsert(x, j, FfFromInt( (randint()%fl) ))
                    FfStepPtr(&(x))
                sig_off()
        else:
            if density == 1:
                fl -= 1
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        FfInsert(x, j, FfFromInt( (randint()%fl)+1 ))
                    FfStepPtr(&(x))
                sig_off()
            else:
                fl -= 1
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        if randdouble() < density:
                            FfInsert(x, j, FfFromInt( (randint()%fl)+1 ))
                    FfStepPtr(&(x))
                sig_off()

## Debugging
#    def show_contents(self, r=None):
#        FfSetField(self.Data.Field)
#        FfSetNoc(self.Data.Noc)
#        cdef PTR p
#        cdef size_t i, j
#        if r is not None:
#            r_min = r
#            r_max = r+1
#        else:
#            r_min = 0
#            r_max = self.Data.Nor
#        for i in range(r_min, r_max):
#            p = FfGetPtr(self.Data.Data, i)
#            for j from 0<=j<self.Data.RowSize:
#                print("%3.3d" % p[j])
#            print

##################
## comparison
    cpdef int _cmp_(left, right) except -2:
        """
        Compare two :class:`Matrix_gfpn_dense` matrices

        Of course, '<' and '>' doesn't make much sense for matrices.

        EXAMPLES::

            sage: M = MatrixSpace(GF(125,'x'),3,20)([20*[0],20*[0],[1]+19*[0]])
            sage: N = copy(M)
            sage: M == N
            True
            sage: M != N
            False
            sage: M < N
            False
            sage: N[2,19] = 1
            sage: M == N
            False
            sage: M != N
            True
        """
        cdef Matrix_gfpn_dense self = left
        cdef Matrix_gfpn_dense N = right
        if self is None or N is None:
            return -1
        cdef char* d1
        cdef char* d2
        if self.Data == NULL:
            if N.Data == NULL:
                return 0
            else:
                return 1
        elif N.Data == NULL:
            return -1
        if self.Data.Field != N.Data.Field:
            if self.Data.Field > N.Data.Field:
                return 1
            return -1
        if self.Data.Noc != N.Data.Noc:
            if self.Data.Noc > N.Data.Noc:
                return 1
            return -1
        if self.Data.Nor != N.Data.Nor:
            if self.Data.Nor > N.Data.Nor:
                return 1
            return -1
        d1 = <char*>(self.Data.Data)
        d2 = <char*>(N.Data.Data)
        cdef str s1 = PyString_FromStringAndSize(d1,self.Data.RowSize * self.Data.Nor)
        cdef str s2 = PyString_FromStringAndSize(d2,N.Data.RowSize * N.Data.Nor)
        if s1 != s2:
            if s1 > s2:
                return 1
            return -1
        return 0

    cdef list _rowlist_(self, i, j=-1):
        "M._rowlist_(i): Return row <i> as a list of python ints"
        cdef int k
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
        else:
            raise ValueError("Matrix is empty")
        if (i<0) or (i>=self.Data.Nor):
            raise IndexError("Index {} out of range 0..{}",format(i,self.Data.Nor-1))
        cdef PTR p
        p = MatGetPtr(self.Data,i)
        L = [FfToInt(FfExtract(p,k)) for k in range(self.Data.Noc)]
        if j!=-1:
            if not(isinstance(j,int) or isinstance(j,Integer)):
                raise TypeError("Second index must be an integer")
            if j >= self.Data.Nor:
                raise IndexError("Index out of range")
            for k from i < k <= j:
                FfStepPtr(&(p)) # This is only called after MatGetPtr, hence, after FfSetNoc.
                L.extend([FfToInt(FfExtract(p,l)) for l in range(self.Data.Noc)])
        return L

    def _list(self):
        """
        Return a flat list of all entries of this matrix.

        The result is cached.

        EXAMPLES::

            sage: MatrixSpace(GF(9,'x'),3)(sorted(list(GF(9,'x')))).list()  # indirect doctest
            [0, 1, 2, x, x + 1, x + 2, 2*x, 2*x + 1, 2*x + 2]

        """
        cdef list x = self.fetch('list')
        if not x is None:
            return x
        x = []
        cdef int i
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
        else:
            raise IndexError("Matrix is empty")
        cdef PTR p
        p = self.Data.Data
        sig_on()
        for i from 1<=i<self.Data.Nor:
            x.extend([self._converter.int_to_field(FfToInt(FfExtract(p,j))) for j in range(self.Data.Noc)])
            FfStepPtr(&(p))
        x.extend([self._converter.int_to_field(FfToInt(FfExtract(p,j))) for j in range(self.Data.Noc)])
        sig_off()
        self.cache('list', x)
        return x

#########################
## Arithmetics
    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col):
        """
        Rescale row number `i` in-place by multiplication with the scalar `s`.

        The argument ``start_col`` is ignored. The scalar `s` is
        converted into the base ring.

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K,5,5)(sorted(list(K)))
            sage: M
            [      0       1       2       3       4]
            [      x   x + 1   x + 2   x + 3   x + 4]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]
            sage: M.rescale_row(1, 3)   # indirect doctest
            sage: M
            [      0       1       2       3       4]
            [    3*x 3*x + 3 3*x + 1 3*x + 4 3*x + 2]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]
            sage: M.rescale_row(4, x)
            sage: M
            [      0       1       2       3       4]
            [    3*x 3*x + 3 3*x + 1 3*x + 4 3*x + 2]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [4*x + 2       2   x + 2 2*x + 2 3*x + 2]

        """
        if start_col != 0 or self.Data == NULL:
            raise ValueError("We can only rescale a full row of a non-empty matrix")
        FfMulRow(MatGetPtr(self.Data, i), FfFromInt(self._converter.field_to_int(self._base_ring(s))))

    cdef add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, multiple, Py_ssize_t start_col):
        """
        Add the ``multiple``-fold of row ``row_from`` in-place to row ``row_to``.

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K,5,5)(sorted(list(K)))
            sage: M
            [      0       1       2       3       4]
            [      x   x + 1   x + 2   x + 3   x + 4]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]
            sage: M.add_multiple_of_row(2, 4, x)  # indirect doctest
            sage: M
            [      0       1       2       3       4]
            [      x   x + 1   x + 2   x + 3   x + 4]
            [  x + 2 2*x + 3 3*x + 4     4*x       1]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]

        """
        if start_col != 0 or self.Data == NULL:
            raise ValueError("We can only rescale a full row of a non-empty matrix")
        FfAddMulRow(MatGetPtr(self.Data, row_to), MatGetPtr(self.Data, row_from), FfFromInt(self._converter.field_to_int(self._base_ring(multiple))))

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        """
        Swap the rows ``row1`` and ``row2`` in-place.

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K,5,5)(sorted(list(K)))
            sage: M
            [      0       1       2       3       4]
            [      x   x + 1   x + 2   x + 3   x + 4]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]
            sage: M.swap_rows(1, 3)    # indirect doctest
            sage: M
            [      0       1       2       3       4]
            [    3*x 3*x + 1 3*x + 2 3*x + 3 3*x + 4]
            [    2*x 2*x + 1 2*x + 2 2*x + 3 2*x + 4]
            [      x   x + 1   x + 2   x + 3   x + 4]
            [    4*x 4*x + 1 4*x + 2 4*x + 3 4*x + 4]

        """
        FfSwapRows(MatGetPtr(self.Data, row1), MatGetPtr(self.Data, row2))

    def trace(self):
        """
        Trace of this matrix, i.e., the sum of diagonal elements.

        EXAMPLES::

            sage: K.<x> = GF(125)
            sage: MatrixSpace(K,7,7)(x).trace()
            2*x

        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        return self._converter.int_to_field(FfToInt(MatTrace(self.Data)))

    def stack(self, Matrix_gfpn_dense other):
        """
        Stack two matrices of the same number of columns.

        EXAMPLES::

            sage: M = MatrixSpace(GF(9,'x'),1,9)(sorted(list(GF(9,'x'))))
            sage: M
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]
            sage: M.stack(M)
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]

        """
        if self._ncols != other._ncols:
            raise TypeError("Both numbers of columns must match.")
        if self._nrows == 0 or self.Data == NULL:
            return other.__copy__()
        if other._nrows == 0 or other.Data == NULL:
            return self.__copy__()
        cdef Matrix_gfpn_dense OUT = self._new(self._nrows+other._nrows, self._ncols)
        OUT.Data = MatAlloc(self.Data.Field, self.Data.Nor+other.Data.Nor, self.Data.Noc)
        memcpy(OUT.Data.Data, self.Data.Data, FfCurrentRowSize*self.Data.Nor)
        memcpy(MatGetPtr(OUT.Data, self.Data.Nor), other.Data.Data, FfCurrentRowSize*other.Data.Nor)
        return OUT

    cpdef _add_(self, right):
        """
        TESTS::

            sage: K.<x> = GF(9)
            sage: M = MatrixSpace(K,3,3)(sorted(list(K)))
            sage: N = MatrixSpace(K,3,3)(2*x)
            sage: M+N           # indirect doctest
            [    2*x       1       2]
            [      x       1   x + 2]
            [    2*x 2*x + 1   x + 2]

        """
        cdef Matrix_gfpn_dense Self = self
        cdef Matrix_gfpn_dense Right = right
        assert Self is not None
        assert Right is not None
        if Self.Data == NULL or Right.Data == NULL:
            raise NotImplementedError("The matrices must not be empty")
        cdef Matrix_gfpn_dense Left = Self.__copy__()
        Left._cache = {}
        MatAdd(Left.Data, Right.Data)
        return Left

    cpdef _sub_(self, right):
        """
        TESTS::

            sage: K.<x> = GF(9)
            sage: M = MatrixSpace(K,3,3)(sorted(list(K)))
            sage: N = MatrixSpace(K,3,3)(2*x)
            sage: M-N    # indirect doctest
            [      x       1       2]
            [      x 2*x + 1   x + 2]
            [    2*x 2*x + 1       2]

        """
        cdef Matrix_gfpn_dense Self = self
        cdef Matrix_gfpn_dense Right = right
        assert Self is not None
        assert Right is not None
        if Self.Data == NULL or Right.Data == NULL:
            raise NotImplementedError("The matrices must not be empty")
        cdef Matrix_gfpn_dense Left = Self.__copy__()
        Left._is_immutable = False
        Left._cache = {}
        MatAddMul(Left.Data, Right.Data, mtx_taddinv[1])
        return Left

    def __neg__(self):
        """
        TESTS::

            sage: M = MatrixSpace(GF(9,'x'),3,3)(sorted(list(GF(9,'x'))))
            sage: -M
            [      0       2       1]
            [    2*x 2*x + 2 2*x + 1]
            [      x   x + 2   x + 1]

        ::

            sage: M = MatrixSpace(GF(125,'x'),10,30).random_element()
            sage: N = MatrixSpace(GF(125,'x'),10,30).random_element()
            sage: M + (-N) == M - N == -(N - M)
            True

        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        return self._lmul_(self._base_ring(-1))

    cpdef _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: M = MatrixSpace(GF(9,'x'),3,3)(sorted(list(GF(9,'x'))))
            sage: K.<x> = GF(9)
            sage: M = MatrixSpace(K,3,3)(sorted(list(K)))
            sage: x*M    # indirect doctest
            [      0       x     2*x]
            [  x + 1 2*x + 1       1]
            [2*x + 2       2   x + 2]
            sage: M*x    # indirect doctest
            [      0       x     2*x]
            [  x + 1 2*x + 1       1]
            [2*x + 2       2   x + 2]
            sage: -M == (-1)*M
            True

        """
        if self.Data == NULL:
            return self.__copy__()
        FfSetField(self.Data.Field)
        cdef Matrix_gfpn_dense OUT = self.__copy__()
        OUT._cache = {}
        MatMulScalar(OUT.Data, FfFromInt(self._converter.field_to_int(right)))
        return OUT

    cdef int _strassen_default_cutoff(self, sage.matrix.matrix0.Matrix right) except -2:
        # Surprisingly, Winograd-Strassen can compete with school book
        # multiplication for smallish matrices, and of course it is
        # asymptotically faster. So, we used it by default.
        return 0

    cpdef Matrix_gfpn_dense _multiply_classical(Matrix_gfpn_dense self, Matrix_gfpn_dense right):
        """
        Multiplication using the cubic school book multiplication algorithm.

        EXAMPLES:

        Since by default the asymptotically faster Strassen-Winograd
        multiplication algorithm is used, the following is a valid
        consistency check::

            sage: M = MatrixSpace(GF(9,'x'),1000,500).random_element()
            sage: N = MatrixSpace(GF(9,'x'),500,2000).random_element()
            sage: M*N == M._multiply_classical(N)                       # optional: meataxe
            True

        """
        "multiply two meataxe matrices by the school book algorithm"
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        cdef Matrix_gfpn_dense OUT = self._new(self._nrows, right._ncols)
        sig_on()
        OUT.Data = MatDup(self.Data)
        MatMul(OUT.Data,right.Data)
        sig_off()
        OUT._is_immutable = False
        OUT._cache = {}
        return OUT

    cpdef Matrix_gfpn_dense _multiply_strassen(Matrix_gfpn_dense self, Matrix_gfpn_dense right, cutoff=0):
        """
        Matrix multiplication using the asymptotically fast Strassen-Winograd algorithm.

        INPUT:

        - ``right`` -- a matrix of dimensions suitable to do multiplication
        - ``cutoff`` (optional integer) -- indicates the minimal size of submatrices
          that will be considered in the divide-and-conquer algorithm. The size is
          *not* expressed by the number of rows/columns, but the rowsize expressed
          in bytes. Depending on the base field, one byte may represent up to eight
          entries in a matrix row. The default is ``sizeof(long)^2/2`` byte.

        EXAMPLES:

        We test that different cutoffs yield the same result::

            sage: M = MatrixSpace(GF(9,'x'),1500,600).random_element()
            sage: N = MatrixSpace(GF(9,'x'),600,1500).random_element()
            sage: M._multiply_strassen(N) == M._multiply_strassen(N,80) == M._multiply_strassen(N,2) # optional: meataxe
            True

        """
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        MS = self.matrix_space(self._nrows, right._ncols, False)
        cdef Matrix_gfpn_dense OUT = Matrix_gfpn_dense(MS, None)
        # Now, OUT.Data is initialised, which is needed for MatMulStrassen to work.
        cutoff = cutoff//sizeof(long)
        StrassenSetCutoff(cutoff)
        sig_on()
        MatMulStrassen(OUT.Data, self.Data, right.Data)
        sig_off()
        return OUT

    cdef _mul_long(self, long n):
        """
        Multiply an MTX matrix with a field element represented by an integer
        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        cdef Matrix_gfpn_dense left
        cdef FEL r
        if n < 0:
            r = mtx_taddinv[FfFromInt(-n)]
        else:
            r = FfFromInt(n)
        left = self.__copy__()
        left._cache = {}
        MatMulScalar(left.Data, r)
        return left

    def __div__(Matrix_gfpn_dense self, p):
        """
        Divide a matrix by a scalar.

        EXAMPLES::

            sage: K.<x> = GF(9)
            sage: M = MatrixSpace(K,3,3)(sorted(list(K)))
            sage: M
            [      0       1       2]
            [      x   x + 1   x + 2]
            [    2*x 2*x + 1 2*x + 2]
            sage: M/2                   # indirect doctest
            [      0       2       1]
            [    2*x 2*x + 2 2*x + 1]
            [      x   x + 2   x + 1]
            sage: M/x
            [      0   x + 2 2*x + 1]
            [      1       x 2*x + 2]
            [      2   x + 1     2*x]

        """
        if self.Data == NULL:
            return self.__copy__()
        if not p:
            raise ZeroDivisionError
        if p not in self._base_ring:
            raise ValueError("{} is not a scalar".format(p))
        p = self._base_ring(p)
        FfSetField(self.Data.Field)
        cdef Matrix_gfpn_dense OUT = self.__copy__()
        OUT._cache = {}
        cdef FEL r = mtx_tmultinv[FfFromInt(self._converter.field_to_int(p))]
        MatMulScalar(OUT.Data, r)
        return OUT

    def __invert__(Matrix_gfpn_dense self):
        """
        Multiplicative inverse of this matrix (if available)

        TESTS::

            sage: MS = MatrixSpace(GF(9,'x'),500)
            sage: while 1:
            ....:     M = MS.random_element()
            ....:     if M.rank() == 500:
            ....:         break
            sage: Minv = ~M    # indirect doctest
            sage: Minv*M == M*Minv == 1
            True

        We use the occasion to demonstrate that errors in MeatAxe are
        correctly handled in Sage::

            sage: MS = MatrixSpace(GF(25,'x'),5)
            sage: while 1:
            ....:     M = MS.random_element(density=0.4)
            ....:     if M.rank() < 5:
            ....:         break
            sage: ~M                    # optional: meataxe
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero in file matinv.c (line 50)

        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        cdef Matrix_gfpn_dense OUT = self._new(self._nrows, self._ncols)
        OUT._is_immutable = False
        OUT._cache = {}
        sig_on()
        try:
            OUT.Data = MatInverse(self.Data)
        except ZeroDivisionError:
            # Attempting to invert singular matrices happens
            # in the tests, and we make the special case here
            # so that the sig_on/off count is fine.
            sig_off()
            raise
        sig_off()
        if OUT.Data != NULL:
            return OUT
        raise ArithmeticError("This matrix is not invertible")

    def transpose(Matrix_gfpn_dense self):
        """
        Return the transposed matrix.

        EXAMPLES::

            sage: K.<x> = GF(9)
            sage: M = MatrixSpace(K, 2,4)(sorted(list(K)[1:]))
            sage: M
            [      1       2       x   x + 1]
            [  x + 2     2*x 2*x + 1 2*x + 2]
            sage: M.transpose()
            [      1   x + 2]
            [      2     2*x]
            [      x 2*x + 1]
            [  x + 1 2*x + 2]

        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        cdef Matrix_gfpn_dense OUT = self._new(self._ncols, self._nrows)
        OUT._is_immutable = False
        OUT._cache = {}
        OUT.Data = MatTransposed(self.Data)
        return OUT

    def order(self):
        """
        Return the multiplicative order of this matrix.

        EXAMPLES::

            sage: K.<x> = GF(27)
            sage: M = MatrixSpace(K, 4)([2*x^2 + 2*x, 2*x^2 + x, 2*x^2 + x + 1,
            ....: x^2 + x + 2, x + 2, x^2, 2*x + 2, 2*x^2 + 2*x, 2*x^2 + 1,
            ....: 1, 2, x^2 + 2*x + 1, x^2 + x + 2, x + 1, 2*x^2 + 2*x, x^2 + x])
            sage: M.order()                 # optional: meataxe
            104
            sage: M^104 == 1
            True
            sage: M^103 == 1
            False

        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if (self.Data.Nor <> self.Data.Noc):
            raise ValueError("only defined for square matrices")
        o = MatOrder(self.Data)
        if o==-1:
            raise ArithmeticError("order too large")
        else:
            return o

###################
## Gauss algorithm

    def left_kernel_matrix(self):
        """
        Return the null space of this matrix, represented as a matrix.

        NOTE:

        - For a matrix `M`, ``M.left_kernel_matrix()*M`` is a null matrix.
        - The command `M.left_kernel()` uses a generic implementation in Sage,
          that relies on computing the echelon form of the transposed
          matrix. This method however uses a MeatAxe function to compute
          the left kernel matrix.

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K, 10)()
            sage: entries = [((0, 2), x), ((0, 4), 3*x + 2),
            ....: ((0, 8), 2*x), ((1, 1), x + 3), ((1, 5), 3*x),
            ....: ((1, 6), x + 4), ((2, 3), 2*x), ((2, 5), 4*x + 1),
            ....: ((2, 6), 4), ((3, 4), x + 4), ((3, 5), x + 1),
            ....: ((5, 5), 3*x), ((5, 7), x + 3), ((6, 1), x),
            ....: ((6, 2), x + 1), ((6, 5), x + 1), ((8, 2), 4),
            ....: ((8, 8), 4), ((8, 9), x + 3), ((9, 8), 4*x + 2)]
            sage: for (i,j),v in entries: M[i,j] = v
            sage: M.left_kernel()
            Vector space of degree 10 and dimension 2 over Finite Field in x of size 5^2
            Basis matrix:
            [0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 0]
            sage: M.left_kernel_matrix()    # optional: meataxe
            [0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 0]

        """
        cdef Matrix_gfpn_dense OUT = self.fetch("left_kernel_matrix")
        if OUT is not None:
            return OUT
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        OUT = type(self).__new__(type(self))
        OUT.Data = MatNullSpace(self.Data)
        OUT._nrows = OUT.Data.Nor
        OUT._ncols = OUT.Data.Noc
        OUT._is_immutable = False
        OUT._parent = self.matrix_space(OUT._nrows, OUT._ncols, False)
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        self.cache("left_kernel_matrix", OUT)
        return OUT

    def _echelon_in_place_classical(self, reduced=True, **kwds):
        """
        Change this matrix into echelon form, using classical Gaussian elimination.

        INPUT:

        - ``reduced`` (optional, default ``True``) -- will result
          in the row-reduced echelon form (otherwise, only a
          semi-echelon form results).

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K, 10)()
            sage: entries = [((0, 2), x), ((0, 4), 3*x + 2),
            ....: ((0, 8), 2*x), ((1, 1), x + 3), ((1, 5), 3*x),
            ....: ((1, 6), x + 4), ((2, 3), 2*x), ((2, 5), 4*x + 1),
            ....: ((2, 6), 4), ((3, 4), x + 4), ((3, 5), x + 1),
            ....: ((5, 5), 3*x), ((5, 7), x + 3), ((6, 1), x),
            ....: ((6, 2), x + 1), ((6, 5), x + 1), ((8, 2), 4),
            ....: ((8, 8), 4), ((8, 9), x + 3), ((9, 8), 4*x + 2)]
            sage: for (i,j),v in entries: M[i,j] = v
            sage: M
            [      0       0       x       0 3*x + 2       0       0       0     2*x       0]
            [      0   x + 3       0       0       0     3*x   x + 4       0       0       0]
            [      0       0       0     2*x       0 4*x + 1       4       0       0       0]
            [      0       0       0       0   x + 4   x + 1       0       0       0       0]
            [      0       0       0       0       0       0       0       0       0       0]
            [      0       0       0       0       0     3*x       0   x + 3       0       0]
            [      0       x   x + 1       0       0   x + 1       0       0       0       0]
            [      0       0       0       0       0       0       0       0       0       0]
            [      0       0       4       0       0       0       0       0       4   x + 3]
            [      0       0       0       0       0       0       0       0 4*x + 2       0]
            sage: M.echelon_form()   # indirect doctest
            [      0       1       0       0       0       0       0       0       0 4*x + 4]
            [      0       0       1       0       0       0       0       0       0 4*x + 2]
            [      0       0       0       1       0       0       0       0       0 3*x + 4]
            [      0       0       0       0       1       0       0       0       0 3*x + 3]
            [      0       0       0       0       0       1       0       0       0 2*x + 3]
            [      0       0       0       0       0       0       1       0       0       x]
            [      0       0       0       0       0       0       0       1       0 2*x + 2]
            [      0       0       0       0       0       0       0       0       1       0]
            [      0       0       0       0       0       0       0       0       0       0]
            [      0       0       0       0       0       0       0       0       0       0]

        A semi-echelon form can be produced by invoking the single-underscore
        method directly::

            sage: N = copy(M)
            sage: N._echelon_in_place_classical(reduced=False)      # optional: meataxe
            sage: N                                                 # optional: meataxe
            [      0       0       x       0 3*x + 2       0       0       0     2*x       0]
            [      0   x + 3       0       0       0     3*x   x + 4       0       0       0]
            [      0       0       0     2*x       0 4*x + 1       4       0       0       0]
            [      0       0       0       0   x + 4   x + 1       0       0       0       0]
            [      0       0       0       0       0     3*x       0   x + 3       0       0]
            [      0       0       0       0       0       0 2*x + 2     4*x 3*x + 3       0]
            [      0       0       0       0       0       0       0   x + 1       1   x + 3]
            [      0       0       0       0       0       0       0       0 4*x + 2       0]
            [      0       0       0       0       0       0       0       0       0       0]
            [      0       0       0       0       0       0       0       0       0       0]

        TESTS:

        We verify that the above echelon form is consistent with Sage's generic
        implementation of dense matrices::

            sage: type(M)                           # optional: meataxe
            <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: MS = M.parent()
            sage: from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            sage: MS._MatrixSpace__matrix_class = Matrix_generic_dense
            sage: X = MS(M._list())
            sage: type(X)
            <type 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: X.echelon_form()
            [      0       1       0       0       0       0       0       0       0 4*x + 4]
            [      0       0       1       0       0       0       0       0       0 4*x + 2]
            [      0       0       0       1       0       0       0       0       0 3*x + 4]
            [      0       0       0       0       1       0       0       0       0 3*x + 3]
            [      0       0       0       0       0       1       0       0       0 2*x + 3]
            [      0       0       0       0       0       0       1       0       0       x]
            [      0       0       0       0       0       0       0       1       0 2*x + 2]
            [      0       0       0       0       0       0       0       0       1       0]
            [      0       0       0       0       0       0       0       0       0       0]
            [      0       0       0       0       0       0       0       0       0       0]

        The following was a problem in a preliminary version of the code::

            sage: K.<a> = GF(25)
            sage: M = MatrixSpace(K, 2, 4)([4, 4, 1, 0, 0, 2*a+1, a+2, 1])
            sage: M
            [      4       4       1       0]
            [      0 2*a + 1   a + 2       1]
            sage: M.echelonize()
            sage: M
            [      1       0 3*a + 4 2*a + 2]
            [      0       1     2*a 3*a + 3]

        """
        if self._nrows == 0 or self._ncols == 0:
            self.cache('in_echelon_form',True)
            self.cache('rank', 0)
            self.cache('pivots', ())
            return self
        MatEchelonize(self.Data)
        self._cache = {}
        # Now, self.Data is in semi-echelon form.
        r = self.Data.Nor
        cdef size_t i, j, pos
        cdef PTR old, dest, src
        cdef FEL piv
        self.cache('rank', r)
        # Next, we do permutations to achieve the reduced echelon form,
        # if requested.
        sig_on()
        if reduced:
            pivs = [(self.Data.PivotTable[i],i) for i in range(r)]
            pivs.sort()
            if pivs != [(self.Data.PivotTable[i],i) for i in range(r)] or self.Data.Nor < self._nrows:
                # We copy the row one by one, sorting their pivot positions
                old = self.Data.Data
                self.Data.Data = FfAlloc(self._nrows)
                for i, (pos,j) in enumerate(pivs):
                    # We have to move row j to row i
                    dest = self.Data.Data+FfCurrentRowSize*i
                    memcpy(dest, old+FfCurrentRowSize*j, FfCurrentRowSize)
                    self.Data.PivotTable[i] = pos
                free(old)
                self.Data.Nor = self._nrows
            # Now, the pivot columns are strictly increasing.
            # We now normalize each row, and annulate everything
            # above the pivot (currently, we only know that the matrix
            # is zero below the pivots).
            for i from 0 <= i < r:
                src = MatGetPtr(self.Data, i)
                piv = FfExtract(src, self.Data.PivotTable[i])
                assert piv!=FF_ZERO
                if piv != FF_ONE:
                    FfMulRow(src, mtx_tmultinv[piv])
                for j from 0 <= j < i:
                    dest = MatGetPtr(self.Data, j)
                    piv = FfExtract(dest, self.Data.PivotTable[i])
                    if piv != FF_ONE:
                        FfAddMulRow(dest, src, mtx_taddinv[piv])
                    else:
                        FfSubRow(dest, src)
        elif self.Data.Nor < self._nrows:
            # Some rows may have vanished. In SageMath, we
            # want that the number of rows does not change,
            # thus, we have to append zero rows.
            self.Data.Data = <PTR>check_realloc(self.Data.Data, FfCurrentRowSize*self._nrows)
            memset(self.Data.Data + FfCurrentRowSize*self.Data.Nor, FF_ZERO, FfCurrentRowSize*(self._nrows-self.Data.Nor))
            self.Data.Nor = self._nrows
        sig_off()
        self.cache('pivots', tuple(self.Data.PivotTable[i] for i in range(r)))
        self.cache('in_echelon_form',True)

def mtx_unpickle(f, int nr, int nc, str Data, bint m):
    """
    Helper function for unpickling.

    TESTS::

        sage: M = MatrixSpace(GF(9,'x'),10,10).random_element()
        sage: M == loads(dumps(M))   # indirect doctest
        True
        sage: M is loads(dumps(M))
        False
    """
    cdef Matrix_gfpn_dense OUT
    OUT = Matrix_gfpn_dense.__new__(Matrix_gfpn_dense)
    if isinstance(f, (int, long)):
        # This is for old pickles created with the group cohomology spkg
        Matrix_dense.__init__(OUT, MatrixSpace(GF(f, 'z'), nr, nc))
    else:
        Matrix_dense.__init__(OUT, f)
        f = OUT._base_ring.order()
    OUT.Data = MatAlloc(f, nr, nc)
    OUT._is_immutable = not m
    OUT._converter = FieldConverter(OUT._base_ring)
    cdef char *x
    if Data:
        x = PyString_AsString(Data)
        memcpy(OUT.Data.Data, x, OUT.Data.RowSize*OUT.Data.Nor)
    return OUT
