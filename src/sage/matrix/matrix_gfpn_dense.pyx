# sage_setup: distribution = sagemath-meataxe
# distutils: libraries = mtx
# sage.doctest: optional - meataxe

r"""
Dense Matrices over `\mathbb F_q`, with `q<255`.

This module is a wrapper for the `Aachen C-MeatAxe library <../spkg/meataxe.html>`_,
improved by an implementation of the Winograd-Strassen multiplication
algorithm. It provides matrices over the finite field `\mathbb F_q`,
where `q\le 255`.

By default, it is only used when `q` is odd and not prime, because other
matrix implementations in SageMath perform better for prime fields or in
characteristic two.

.. NOTE::

    The examples shown here will only work when the `meataxe
    <../spkg/meataxe.html>` package has been installed.

AUTHORS:

- Simon King (2015-09): initial version
"""

# ***************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from cysignals.memory cimport check_realloc, check_malloc, sig_free
from cpython.bytes cimport PyBytes_AsString, PyBytes_FromStringAndSize
from cysignals.signals cimport sig_on, sig_off, sig_check
cimport cython

####################
#
# import sage types
#
####################

from sage.cpython.string cimport str_to_bytes
from sage.cpython.string import FS_ENCODING
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod import IntegerMod_int
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.randstate import current_randstate
from sage.misc.randstate cimport randstate
from sage.structure.element cimport Element, Matrix
from sage.structure.richcmp import rich_to_bool
from sage.matrix.args cimport MatrixArgs_init

from libc.string cimport memset, memcpy

cimport sage.matrix.matrix0

# The following import is just to ensure that meataxe_init() is called.
import sage.libs.meataxe

####################
#
# auxiliary functions
#
####################

# Fast conversion from field to int and int to field
cdef class FieldConverter_class:
    """
    An auxiliary class, used to convert between meataxe ``FEL`` and
    finite field elements.

    This class is for Givaro finite fields only. It exists in order to
    have a common interface for elements of prime and non-prime fields;
    see :class:`PrimeFieldConverter_class`.

    .. NOTE::

        This class is really meant to be used only from Cython.
        We do make the methods available from Python, but that is
        for testing only.

    .. WARNING::

        Before calling the ``fel_to_field`` or ``field_to_fel`` methods,
        one should call the ``FfSetField`` function.

    EXAMPLES::

        sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
        sage: F.<y> = GF(125)
        sage: M = MatrixSpace(F, 2, 2, implementation=Matrix_gfpn_dense).one()
        sage: C = M._converter
        sage: C.fel_to_field(15)
        3*y
        sage: F.from_integer(15)
        3*y
        sage: C.field_to_fel(y)
        5
        sage: y.to_integer()
        5
    """
    def __init__(self, field):
        """
        INPUT:

        - ``field`` -- a finite field with Givaro implementation and at most 251
          elements. These assumptions are not tested.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: FieldConverter_class(GF(13^2))
            <sage.matrix.matrix_gfpn_dense.FieldConverter_class object at ...>
        """
        self.field = field._cache.fetch_int
        self.zero_FEL = self.field_to_fel(field.zero())

    cpdef fel_to_field(self, FEL x) noexcept:
        """
        Fetch a python int into the field.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)
            sage: C.fel_to_field(15)
            3*y
            sage: F.from_integer(15)
            3*y
        """
        return self.field(FfToInt(x))

    cpdef FEL field_to_fel(self, x) except 255:
        """
        Represent a field element by a python int.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)
            sage: C.field_to_fel(y)
            5
            sage: y.to_integer()
            5

        TESTS:

        Test invalid input::

            sage: C.field_to_fel('foo')
            Traceback (most recent call last):
            ...
            AttributeError: 'str' object has no attribute 'to_integer'...
        """
        return FfFromInt(x.to_integer())


cdef class PrimeFieldConverter_class(FieldConverter_class):
    """
    An auxiliary class, used to convert between meataxe ``FEL`` and
    finite field elements.

    This class is for prime fields only. It exists in order to have a
    common interface for elements of prime and non-prime fields; see
    :class:`FieldConverter_class`.

    .. NOTE::

        This class is really meant to be used only from Cython.
        We do make the methods available from Python, but that is
        for testing only.

    .. WARNING::

        Before calling the ``fel_to_field`` or ``field_to_fel`` methods,
        one should call the ``FfSetField`` function.

    EXAMPLES::

        sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
        sage: F = GF(5)
        sage: M = MatrixSpace(F, 2, 2, implementation=Matrix_gfpn_dense).one()
        sage: C = M._converter
        sage: C.fel_to_field(2)
        2
        sage: F(2)
        2
        sage: C.field_to_fel(F(2))
        2
        sage: int(F(2))
        2
    """
    def __init__(self, field):
        """
        INPUT:

        - ``field`` -- a finite *prime* field with at most 251 elements.
          This assumption is not tested.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: PrimeFieldConverter_class(GF(251))
            <sage.matrix.matrix_gfpn_dense.PrimeFieldConverter_class object at ...>
        """
        self.field = field

    cpdef fel_to_field(self, FEL x) noexcept:
        """
        Fetch a python int into the field.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)
            sage: C.fel_to_field(2)
            2
        """
        return IntegerMod_int(self.field, FfToInt(x))

    cpdef FEL field_to_fel(self, x) except 255:
        """
        Represent a field element by a python int.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)
            sage: C.field_to_fel(F(2))
            2

        TESTS:

        Test invalid input::

            sage: C.field_to_fel('foo')
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        return FfFromInt(x)


cdef dict _converter_cache = {}
cdef FieldConverter_class FieldConverter(field) noexcept:
    """
    Return a :class:`FieldConverter_class` or :class:`PrimeFieldConverter_class` instance,
    depending whether the field is prime or not.

    EXAMPLES::

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
##
## Wrapper for MeatAxe matrices
##
######################################

cdef Matrix_gfpn_dense new_mtx(Matrix_t* mat, Matrix_gfpn_dense template) noexcept:
    """
    Create a new ``Matrix_gfpn_dense`` from a meataxe matrix

    INPUT:

    - ``mat`` -- (``Matrix_t*``) a meataxe matrix

    - ``template`` -- an optional matrix with the same base ring as the
      new matrix being created. This can also be ``None``.

    ``template`` is only used for efficiency: various objects are simply
    taken from the template instead of being constructed. The condition
    that the base ring is the same is not checked.
    """
    assert mat is not NULL

    if template is not None:
        B = template._base_ring
        conv = template._converter
        src = template.Data
        if src is not NULL and src.Nor == mat.Nor and src.Noc == mat.Noc:
            # Same dimensions
            MS = template._parent
        else:
            # Different dimensions
            MS = MatrixSpace(B, mat.Nor, mat.Noc, implementation=Matrix_gfpn_dense)
    else:
        B = GF(mat.Field, 'z')
        conv = FieldConverter(B)
        MS = MatrixSpace(B, mat.Nor, mat.Noc, implementation=Matrix_gfpn_dense)

    ret = <Matrix_gfpn_dense>Matrix_gfpn_dense.__new__(Matrix_gfpn_dense, MS)
    ret.Data = mat
    ret._converter = conv
    return ret


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
        sage: type(M)
        <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>

    The documentation of the ``__init__`` methods shows further
    ways of creating a :class:`Matrix_gfpn_dense` instance.
    However, these should only be of internal use.
    """
##################
## Init, Dealloc, Copy

    def __dealloc__(self):
        """
        TESTS::

            sage: MS = MatrixSpace(GF(64), 0)
            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
            sage: Matrix_gfpn_dense.__new__(Matrix_gfpn_dense, MS)
            []
            sage: M = None
            sage: M = Matrix_gfpn_dense(MatrixSpace(GF(64,'z'),4), None)
            sage: del M    # indirect doctest
        """
        if self.Data != NULL:
            MatFree(self.Data)

    def __init__(self, parent, entries=None, copy=None, bint coerce=False, *, bint mutable=True):
        r"""
        Matrix extension class using libmeataxe as backend.

        INPUT:

        - ``parent`` -- a matrix space over ``GF(q)`` with `q < 255`

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- ignored

        - ``mutable`` -- if ``False``, the resulting matrix cannot be
          changed, and it can be used as key in a Python dictionary

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense

        1. Creating a zero (3x2)-matrix::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(4,'z'),3,2))
            [0 0]
            [0 0]
            [0 0]

        2. Creating a matrix from a list or list of lists::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[1,2,3,4,5,6])
            [1 2 3]
            [4 0 1]
            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[[1,2,3],[4,5,6]])
            [1 2 3]
            [4 0 1]

        3. Creating a diagonal matrix::

            sage: M = Matrix_gfpn_dense(MatrixSpace(GF(7),5),2); M
            [2 0 0 0 0]
            [0 2 0 0 0]
            [0 0 2 0 0]
            [0 0 0 2 0]
            [0 0 0 0 2]

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
        ma = MatrixArgs_init(parent, entries)
        cdef long fl = ma.base.order()
        cdef long nr = ma.nrows
        cdef long nc = ma.ncols

        self.Data = MatAlloc(fl, nr, nc)
        self._converter = FieldConverter(ma.base)

        cdef PTR x = self.Data.Data
        assert self.Data.Noc == nc
        assert self.Data.Nor == nr
        FfSetField(fl)
        FfSetNoc(nc)
        it = ma.iter(False)
        cdef long i,j
        for i in range(nr):
            for j in range(nc):
                c = self._converter.field_to_fel(self._coerce_element(next(it)))
                FfInsert(x, j, c)
                sig_check()
            FfStepPtr(&x)

        self._is_immutable = not mutable

    @staticmethod
    def from_filename(filename):
        """
        Create a matrix from a file in MeatAxe format. If the file
        does not exist, an error raised by the MeatAxe library is
        propagated::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
            sage: Matrix_gfpn_dense.from_filename('foobarNONEXISTING_FILE')
            Traceback (most recent call last):
            ...
            OSError: foobarNONEXISTING_FILE: No such file or directory in file os.c (line 255)
            sage: Matrix_gfpn_dense.from_filename('')
            Traceback (most recent call last):
            ...
            ValueError: cannot construct meataxe matrix from empty filename
        """
        if not filename:
            raise ValueError("cannot construct meataxe matrix from empty filename")

        if type(filename) is not bytes:
            filename = str_to_bytes(filename, FS_ENCODING,
                                    'surrogateescape')
        sig_on()
        try:
            mat = MatLoad(filename)
        finally:
            sig_off()
        return new_mtx(mat, None)

    def __copy__(self):
        """
        Return a mutable copy of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(GF(25,'x'), 3, 20)([20*[0],20*[0],[1]+19*[0]])
            sage: N = copy(M)   # indirect doctest
            sage: print(N)
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            sage: N == M
            True
            sage: N is M
            False
            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
            sage: M = Matrix_gfpn_dense.__new__(Matrix_gfpn_dense, parent(M))
            sage: copy(M)
            Traceback (most recent call last):
            ...
            ValueError: cannot copy an uninitialized matrix
        """
        if self.Data is NULL:
            raise ValueError("cannot copy an uninitialized matrix")
        sig_on()
        try:
            mat = MatDup(self.Data)
        finally:
            sig_off()
        retval = new_mtx(mat, self)
        if self._cache:
            retval._cache = dict(self._cache)
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
        cdef char* x
        cdef size_t i
        cdef PTR p
        cdef size_t pickle_size
        cdef bytes pickle_str
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
            pickle_size = FfCurrentRowSizeIo*self.Data.Nor
            d = <char*>check_malloc(pickle_size)
            p = self.Data.Data
            x = d
            for i in range(self.Data.Nor):
                memcpy(x, p, FfCurrentRowSizeIo)
                sig_check()
                x += FfCurrentRowSizeIo
                FfStepPtr(&p)
            pickle_str = PyBytes_FromStringAndSize(d, pickle_size)
            sig_free(d)
            return mtx_unpickle, (self._parent, self.Data.Nor, self.Data.Noc,
                        pickle_str,
                        not self._is_immutable) # for backward compatibility with the group cohomology package
        else:
            return mtx_unpickle, (0, 0, 0, '', not self._is_immutable)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j) noexcept:
        """
        Get an element without checking.

        TESTS::

            sage: F.<z> = GF(9)
            sage: M = MatrixSpace(F,3)(sorted(list(F)))
            sage: type(M)
            <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: M                     # indirect doctest
            [      0       1       2]
            [      z   z + 1   z + 2]
            [    2*z 2*z + 1 2*z + 2]
        """
        if self.Data == NULL:
            raise IndexError("Matrix is empty")
        FfSetField(self.Data.Field)
        return self._converter.fel_to_field(FfExtract(MatGetPtr(self.Data,i), j))

    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j) noexcept:
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        return FfToInt(FfExtract(MatGetPtr(self.Data,i), j))

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: F.<z> = GF(9)
            sage: M = MatrixSpace(F,2,5)(sorted(list(F))+[0])
            sage: M.zero_pattern_matrix()  # indirect doctest
            [1 0 0 0 0]
            [0 0 0 0 1]
        """
        FfSetField(self.Data.Field)
        return FfExtract(MatGetPtr(self.Data,i), j) == self._converter.zero_FEL

    cpdef Matrix_gfpn_dense get_slice(self, Py_ssize_t i, Py_ssize_t j) noexcept:
        """
        Return a horizontal slice of this matrix.

        NOTE:

        ``M[i:j]`` may return a matrix that uses a different backend than
        MeatAxe. This method is useful when the slice has to be of type
        :class:`Matrix_gfpn_dense`.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense as MTX
            sage: M = MTX(MatrixSpace(GF(7), 5, 3), [[0,1,2], [1,2,3], [2,3,4], [3,4,5], [4,5,6]])
            sage: M
            [0 1 2]
            [1 2 3]
            [2 3 4]
            [3 4 5]
            [4 5 6]
            sage: M.get_slice(1,3)
            [1 2 3]
            [2 3 4]
            sage: type(_) is MTX
            True
        """
        if not 0 <= i < j <= self.Data.Nor:
            raise IndexError("Indices i={}, j={} violate the condition 0 < i < j < {}".format(i,j,self.Data.Nor))
        sig_on()
        try:
            mat = MatAlloc(self.Data.Field, j-i, self.Data.Noc)
            memcpy(mat.Data, FfGetPtr(self.Data.Data, i), FfCurrentRowSize*(j-i))
        finally:
            sig_off()
        return new_mtx(mat, self)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value) noexcept:
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
        FfInsert(MatGetPtr(self.Data,i), j, self._converter.field_to_fel(value))

    cdef void set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value) noexcept:
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        FfInsert(FfGetPtr(self.Data.Data,i), j, FfFromInt(value))

    cdef set_slice_unsafe(self, Py_ssize_t i, Matrix_gfpn_dense S) noexcept:
        # Overwrite the self[i:i+S.nrows()] by the contents of S.
        #
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that the dimensions of self and S match!
        sig_on()
        try:
            memcpy(FfGetPtr(self.Data.Data, i), S.Data.Data, FfCurrentRowSize*S.Data.Nor)
        finally:
            sig_off()

    def randomize(self, density=None, nonzero=False, *args, **kwds):
        """
        Fill the matrix with random values.

        INPUT:

        - ``density`` -- (optional) real number between zero and one the
          expected density of the resulting matrix
        - ``nonzero`` -- boolean (default: ``False``); if ``True``, all
          inserted marks are nonzero

        EXAMPLES::

            sage: MS = MatrixSpace(GF(27,'z'),6,6)
            sage: M = MS.random_element()       # indirect doctest
            sage: type(M)
            <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: M = MS.random_element(nonzero=True)
            sage: all(M[i,j] for i in range(6) for j in range(6))
            True
            sage: avg_density = sum(MS.random_element(density=0.5).density()
            ....:                   for _ in range(100))
            sage: avg_density /= 100
            sage: RR(avg_density)  # abs tol 0.05
            0.5

        TESTS:

        The following tests against a bug that was fixed in :issue:`23352`.
        This test could fail for some seed, but it would be highly unlikely::

            sage: MS = MatrixSpace(GF(9,'x'),1,5)
            sage: any(MS.random_element()[0,4] for _ in range(50))
            True
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
        cdef int i, j

        FfSetField(fl)
        FfSetNoc(nc)
        cdef int MPB, tmp
        cdef unsigned char O
        cdef randstate RandState = current_randstate()

        if not nonzero:
            if density == 1:
                MPB = 0
                tmp = fl
                while tmp <= 256:
                    MPB += 1
                    tmp *= fl
                O = (fl**MPB)
                if nc%MPB:
                    for i in range(nr):
                        y = <unsigned char*>x
                        for j in range(FfCurrentRowSizeIo-1):
                            y[j] = RandState.c_random()%O
                            sig_check()
                        for j in range(nc-(nc%MPB), nc):
                            FfInsert(x, j, FfFromInt( (RandState.c_random()%fl) ))
                            sig_check()
                        FfStepPtr(&(x))
                else:
                    for i in range(nr):
                        y = <unsigned char*>x
                        for j in range(FfCurrentRowSizeIo):
                            y[j] = RandState.c_random()%O
                            sig_check()
                        FfStepPtr(&(x))
            else:
                for i in range(nr):
                    for j in range(nc):
                        if RandState.c_rand_double() < density:
                            FfInsert(x, j, FfFromInt( (RandState.c_random()%fl) ))
                            sig_check()
                    FfStepPtr(&(x))
        else:
            if density == 1:
                fl -= 1
                for i in range(nr):
                    for j in range(nc):
                        FfInsert(x, j, FfFromInt( (RandState.c_random()%fl)+1 ))
                        sig_check()
                    FfStepPtr(&(x))
            else:
                fl -= 1
                for i in range(nr):
                    for j in range(nc):
                        if RandState.c_rand_double() < density:
                            FfInsert(x, j, FfFromInt( (RandState.c_random()%fl)+1 ))
                            sig_check()
                    FfStepPtr(&(x))

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
    cpdef _richcmp_(left, right, int op) noexcept:
        """
        Compare two :class:`Matrix_gfpn_dense` matrices.

        Of course, '<' and '>' do not make much sense for matrices.

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
            return rich_to_bool(op, -1)
        if self.Data == NULL:
            if N.Data == NULL:
                return rich_to_bool(op, 0)
            else:
                return rich_to_bool(op, 1)
        elif N.Data == NULL:
            return rich_to_bool(op, -1)
        if self.Data.Field != N.Data.Field:
            if self.Data.Field > N.Data.Field:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)
        if self.Data.Noc != N.Data.Noc:
            if self.Data.Noc > N.Data.Noc:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)
        if self.Data.Nor != N.Data.Nor:
            if self.Data.Nor > N.Data.Nor:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)

        cdef char* d1 = <char*>self.Data.Data
        cdef char* d2 = <char*>N.Data.Data
        cdef Py_ssize_t total_size = self.Data.RowSize
        total_size *= self.Data.Nor
        cdef bytes s1, s2
        s1 = PyBytes_FromStringAndSize(d1, total_size)
        s2 = PyBytes_FromStringAndSize(d2, total_size)
        if s1 != s2:
            if s1 > s2:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)
        return rich_to_bool(op, 0)

    cpdef list _rowlist_(self, i, j=-1) noexcept:
        """
        Return rows as a flat list of python ints.

        INPUT:

        - ``i`` -- index of the first row to be extracted
        - ``j`` -- (default: -1) -1, or index of the last
          row to be extracted

        OUTPUT:

        If `j=-1` then only the `i`-th row is returned as a list.
        Otherwises, rows `i` to `j` (both included) are returned
        as a list of integers.

        EXAMPLES::

            sage: k.<x> = GF(25,'x')
            sage: M = random_matrix(GF(25,'x'), 5,5)
            sage: M = matrix(k,[[      4,     4*x, x + 3, 4*x + 2, 3*x + 4],
            ....:               [  x + 2, 3*x + 1,     3,       0,       3],
            ....:               [    3*x, 2*x + 4,     1,       0,     2*x],
            ....:               [4*x + 4, 2*x + 3,   4*x,       1, 3*x + 1],
            ....:               [3*x + 3,   x + 3, x + 2,   x + 1, 3*x + 2]])
            sage: M._rowlist_(1)
            [7, 16, 3, 0, 3]
            sage: [M[1,i]._int_repr() for i in range(5)]
            ['7', '16', '3', '0', '3']
            sage: M._rowlist_(2,4)
            [15, 14, 1, 0, 10, 24, 13, 20, 1, 16, 18, 8, 7, 6, 17]
            sage: [[M[i,j]._int_repr() for j in range(5)] for i in range(2,5)]
            [['15', '14', '1', '0', '10'],
             ['24', '13', '20', '1', '16'],
             ['18', '8', '7', '6', '17']]
        """
        cdef int k
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
        else:
            raise ValueError("Matrix is empty")
        if (i < 0) or (i >= self.Data.Nor):
            raise IndexError("Index {} out of range 0..{}",format(i,self.Data.Nor-1))
        cdef PTR p
        p = MatGetPtr(self.Data,i)
        L = [FfToInt(FfExtract(p,k)) for k in range(self.Data.Noc)]
        if j!=-1:
            if not(isinstance(j,int) or isinstance(j,Integer)):
                raise TypeError("Second index must be an integer")
            if j >= self.Data.Nor:
                raise IndexError("Index out of range")
            for k in range(i, j):
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

        TESTS:

        This works correctly on matrices with zero rows or columns::

            sage: M = Matrix(GF(9), 0, 3); M._list()
            []
            sage: M = Matrix(GF(9), 3, 0); M._list()
            []
            sage: M = Matrix(GF(9), 0, 0); M._list()
            []
        """
        L = self.fetch('list')
        if L is not None:
            return L
        if self.Data is NULL:
            raise IndexError("matrix is empty")
        FfSetField(self.Data.Field)
        FfSetNoc(self.Data.Noc)
        cdef list x = []
        cdef int i, j
        p = self.Data.Data
        for i in range(self.Data.Nor):
            sig_check()
            for j in range(self.Data.Noc):
                x.append(self._converter.fel_to_field(FfExtract(p,j)))
            FfStepPtr(&p)
        self.cache('list', x)
        return x

#########################
## Arithmetics
    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col) noexcept:
        """
        Rescale row number `i` in-place by multiplication with the scalar `s`.

        The scalar `s` is converted into the base ring.

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

        TESTS:

        The following tests against bugs fixed in :issue:`25490`. It shows
        that the optional argument ``start_col`` is correctly dealt with,
        in an example where several marks are packed into one byte, and
        so that temporarily the current field and row size are changed.
        ::

            sage: L.<y> = GF(9)
            sage: N = MatrixSpace(L,5,5)((3*sorted(list(L)))[:25])
            sage: N
            [      0       1       2       y   y + 1]
            [  y + 2     2*y 2*y + 1 2*y + 2       0]
            [      1       2       y   y + 1   y + 2]
            [    2*y 2*y + 1 2*y + 2       0       1]
            [      2       y   y + 1   y + 2     2*y]
            sage: M = MatrixSpace(K,1,25)(sorted(list(K)))
            sage: N.rescale_row(1,1/y,1)
            sage: N
            [      0       1       2       y   y + 1]
            [  y + 2       2   y + 1     2*y       0]
            [      1       2       y   y + 1   y + 2]
            [    2*y 2*y + 1 2*y + 2       0       1]
            [      2       y   y + 1   y + 2     2*y]
        """
        if self.Data == NULL or start_col >= self.Data.Noc:
            return
        FfSetField(self.Data.Field)
        cdef FEL c = self._converter.field_to_fel(self._base_ring(s))
        cdef ssize_t byte_offset = start_col//MPB     # how many full bytes will remain unchanged?
        cdef ssize_t remains = start_col % MPB        # how many cols have to be treated separately?
        cdef ssize_t noc                              # what bunch of cols will be treated together?
        cdef int j
        cdef PTR row_head = MatGetPtr(self.Data, i) + byte_offset
        if remains:
            for j in range(remains,MPB):
                FfInsert(row_head, j, FfMul(FfExtract(row_head,j), c))
            byte_offset += 1
            row_head    += 1
        noc = self.Data.Noc - byte_offset*MPB
        if noc:
            FfSetNoc(noc)
            FfMulRow(row_head, c)

    cdef add_multiple_of_row_c(self, Py_ssize_t row_to, Py_ssize_t row_from, multiple, Py_ssize_t start_col) noexcept:
        """
        Add the ``multiple``-fold of row ``row_from`` in-place to row ``row_to``, beginning with ``start_col``

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

        TESTS:

        The following tests against bugs fixed at :issue:`25490`. It demonstrates
        that the optional argument ``start_col`` is correctly dealt with,
        in a situation where several marks are packed into one byte
        and the current field and row size are temporarily changed.
        ::

            sage: L.<y> = GF(9)
            sage: N = MatrixSpace(L,5,5)((3*sorted(list(L)))[:25])
            sage: N
            [      0       1       2       y   y + 1]
            [  y + 2     2*y 2*y + 1 2*y + 2       0]
            [      1       2       y   y + 1   y + 2]
            [    2*y 2*y + 1 2*y + 2       0       1]
            [      2       y   y + 1   y + 2     2*y]
            sage: M = MatrixSpace(K,1,25)(sorted(list(K)))
            sage: N.add_multiple_of_row(2,1,-N[2,1]/N[1,1],1)
            sage: N
            [      0       1       2       y   y + 1]
            [  y + 2     2*y 2*y + 1 2*y + 2       0]
            [      1       0       2 2*y + 1   y + 2]
            [    2*y 2*y + 1 2*y + 2       0       1]
            [      2       y   y + 1   y + 2     2*y]
            sage: y-(2*y+1)/y
            2
        """
        if self.Data == NULL or start_col >= self.Data.Noc:
            return
        FfSetField(self.Data.Field)
        cdef FEL c = self._converter.field_to_fel(self._base_ring(multiple))
        cdef ssize_t byte_offset = start_col//MPB     # how many full bytes will remain unchanged?
        cdef ssize_t remains = start_col % MPB        # how many cols have to be treated separately?
        cdef ssize_t noc                              # what bunch of cols will be treated together?
        cdef PTR row_to_head = MatGetPtr(self.Data, row_to)+byte_offset
        cdef PTR row_from_head = MatGetPtr(self.Data, row_from)+byte_offset
        if remains:
            for j in range(remains,MPB):
                FfInsert(row_to_head, j, FfAdd(FfExtract(row_to_head,j), FfMul(FfExtract(row_from_head,j), c)))
            byte_offset += 1
            row_to_head   += 1
            row_from_head += 1
        noc = self.Data.Noc - byte_offset*MPB
        if noc:
            FfSetNoc(noc)
            FfAddMulRow(row_to_head, row_from_head, c)

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2) noexcept:
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
        if not self.Data:
            raise ValueError("This matrix is empty")
        FfSetField(self.Data.Field)
        FfSetNoc(self.Data.Noc)
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
        return self._converter.fel_to_field(MatTrace(self.Data))

    cdef _stack_impl(self, bottom) noexcept:
        r"""
        Stack ``self`` on top of ``bottom``.

        INPUT:

        - ``bottom`` -- a matrix with the same number of columns as ``self``

        EXAMPLES::

            sage: M = MatrixSpace(GF(9,'x'),1,9)(sorted(list(GF(9,'x'))))
            sage: M
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]
            sage: M.stack(M)
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]
            [      0       1       2       x   x + 1   x + 2     2*x 2*x + 1 2*x + 2]

        Check that we can stack a vector (:issue:`31708`)::

            sage: R.<a> = GF(3^3)
            sage: M = matrix(R, [[1,1],[0,a+1]])
            sage: M.stack(vector(R, [a,0]))
            [    1     1]
            [    0 a + 1]
            [    a     0]
        """
        cdef Matrix_gfpn_dense other = <Matrix_gfpn_dense> bottom

        if self._nrows == 0 or self.Data == NULL:
            return other.__copy__()
        if other._nrows == 0 or other.Data == NULL:
            return self.__copy__()
        sig_on()
        try:
            mat = MatAlloc(self.Data.Field, self.Data.Nor+other.Data.Nor, self.Data.Noc)
            memcpy(mat.Data, self.Data.Data, FfCurrentRowSize*self.Data.Nor)
            memcpy(MatGetPtr(mat, self.Data.Nor), other.Data.Data, FfCurrentRowSize*other.Data.Nor)
        finally:
            sig_off()
        return new_mtx(mat, self)

    cpdef _add_(self, right) noexcept:
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
        sig_on()
        try:
            mat = MatDup(Self.Data)
            MatAdd(mat, Right.Data)
        finally:
            sig_off()
        return new_mtx(mat, self)

    cpdef _sub_(self, right) noexcept:
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
        sig_on()
        try:
            mat = MatDup(Self.Data)
            MatAddMul(mat, Right.Data, mtx_taddinv[1])
        finally:
            sig_off()
        return new_mtx(mat, self)

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

    cpdef _lmul_(self, Element right) noexcept:
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

        TESTS:

        Make sure that :issue:`25076` remains fixed::

            sage: M == M*int(4) == int(4)*M
            True
        """
        if self.Data == NULL:
            return self.__copy__()
        FfSetField(self.Data.Field)
        sig_on()
        try:
            mat = MatDup(self.Data)
            MatMulScalar(mat, self._converter.field_to_fel(right))
        finally:
            sig_off()
        return new_mtx(mat, self)

    cdef int _strassen_default_cutoff(self, sage.matrix.matrix0.Matrix right) except -2:
        # Surprisingly, Winograd-Strassen can compete with school book
        # multiplication for smallish matrices, and of course it is
        # asymptotically faster. So, we used it by default.
        return 0

    cpdef Matrix_gfpn_dense _multiply_classical(Matrix_gfpn_dense self, Matrix_gfpn_dense right) noexcept:
        """
        Multiplication using the cubic school book multiplication algorithm.

        EXAMPLES:

        Since by default the asymptotically faster Strassen-Winograd
        multiplication algorithm is used, the following is a valid
        consistency check::

            sage: M = MatrixSpace(GF(9,'x'),1000,500).random_element()
            sage: N = MatrixSpace(GF(9,'x'),500,2000).random_element()
            sage: M*N == M._multiply_classical(N)
            True
        """
        "multiply two meataxe matrices by the school book algorithm"
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        sig_on()
        try:
            mat = MatDup(self.Data)
            MatMul(mat, right.Data)
        finally:
            sig_off()
        return new_mtx(mat, self)

    cpdef Matrix_gfpn_dense _multiply_strassen(Matrix_gfpn_dense self, Matrix_gfpn_dense right, cutoff=0) noexcept:
        """
        Matrix multiplication using the asymptotically fast Strassen-Winograd algorithm.

        INPUT:

        - ``right`` -- a matrix of dimensions suitable to do multiplication
        - ``cutoff`` -- (optional integer) indicates the minimal size of submatrices
          that will be considered in the divide-and-conquer algorithm. The size is
          *not* expressed by the number of rows/columns, but the rowsize expressed
          in bytes. Depending on the base field, one byte may represent up to eight
          entries in a matrix row. The default is ``sizeof(long)^2/2`` byte.

        EXAMPLES:

        We test that different cutoffs yield the same result::

            sage: M = MatrixSpace(GF(9,'x'),1500,600).random_element()
            sage: N = MatrixSpace(GF(9,'x'),600,1500).random_element()
            sage: M._multiply_strassen(N) == M._multiply_strassen(N,80) == M._multiply_strassen(N,2)
            True
        """
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        StrassenSetCutoff(cutoff // sizeof(long))
        sig_on()
        try:
            mat = MatAlloc(self.Data.Field, self._nrows, right._ncols)
            MatMulStrassen(mat, self.Data, right.Data)
        finally:
            sig_off()
        return new_mtx(mat, self)

    cdef _mul_long(self, long n) noexcept:
        """
        Multiply an MTX matrix with a field element represented by an integer.

        TESTS::

            sage: M = random_matrix(GF(9,'x'), 64,51)
            sage: M == M*int(4) == int(4)*M
            True
            sage: M*int(-1)+M == 0
            True
        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        FfSetField(self.Data.Field)
        cdef FEL r
        with cython.cdivision(False):
            r = FfFromInt(n % FfChar)
        sig_on()
        try:
            mat = MatDup(self.Data)
            MatMulScalar(mat, r)
        finally:
            sig_off()
        return new_mtx(mat, self)

    def __truediv__(Matrix_gfpn_dense self, p):
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
        cdef FEL r = mtx_tmultinv[self._converter.field_to_fel(p)]
        sig_on()
        try:
            mat = MatDup(self.Data)
            MatMulScalar(mat, r)
        finally:
            sig_off()
        return new_mtx(mat, self)

    def __invert__(Matrix_gfpn_dense self):
        """
        Multiplicative inverse of this matrix (if available).

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
            sage: ~M
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero in file matinv.c (line 50)
        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        sig_on()
        try:
            mat = MatInverse(self.Data)
        finally:
            # Inverting a singular matrix causes MatInverse to raise a
            # ZeroDivisionError (see sage_meataxe_error_handler).
            # So we need to support exceptions here.
            sig_off()
        return new_mtx(mat, self)

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
        sig_on()
        try:
            mat = MatTransposed(self.Data)
        finally:
            sig_off()
        return new_mtx(mat, self)

    def order(self):
        """
        Return the multiplicative order of this matrix.

        EXAMPLES::

            sage: K.<x> = GF(27)
            sage: M = MatrixSpace(K, 4)([2*x^2 + 2*x, 2*x^2 + x, 2*x^2 + x + 1,
            ....: x^2 + x + 2, x + 2, x^2, 2*x + 2, 2*x^2 + 2*x, 2*x^2 + 1,
            ....: 1, 2, x^2 + 2*x + 1, x^2 + x + 2, x + 1, 2*x^2 + 2*x, x^2 + x])
            sage: M.order()
            104
            sage: M^104 == 1
            True
            sage: M^103 == 1
            False
        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if self.Data.Nor != self.Data.Noc:
            raise ValueError("only defined for square matrices")
        sig_on()
        try:
            o = MatOrder(self.Data)
        finally:
            sig_off()
        if o == -1:
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
            sage: M.left_kernel_matrix()
            [0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 0]
        """
        cdef Matrix_gfpn_dense OUT = self.fetch("left_kernel_matrix")
        if OUT is not None:
            return OUT
        if self.Data is NULL:
            raise ValueError("The matrix must not be empty")
        sig_on()
        try:
            mat = MatNullSpace(self.Data)
        finally:
            sig_off()
        OUT = new_mtx(mat, self)
        self.cache("left_kernel_matrix", OUT)
        return OUT

    cpdef _echelon_in_place(self, str algorithm) noexcept:
        """
        Change this matrix into echelon form, using classical
        Gaussian elimination, and return the pivots.

        .. NOTE::

            Use :meth:`_echelon_in_place_classical`, which can take the
            keyword ``reduced``.

        EXAMPLES::

            sage: K.<x> = GF(25)
            sage: M = MatrixSpace(K, 9, 3)()
            sage: entries = [((0, 2), x), ((0, 4), 3*x + 2),
            ....: ((0, 8), 2*x), ((1, 1), x + 3), ((1, 5), 3*x),
            ....: ((1, 6), x + 4), ((2, 0), 2*x), ((2, 5), 4*x + 1)]
            sage: for (i,j),v in entries: M[j,i] = v
            sage: M
            [      0       0     2*x]
            [      0   x + 3       0]
            [      x       0       0]
            [      0       0       0]
            [3*x + 2       0       0]
            [      0     3*x 4*x + 1]
            [      0   x + 4       0]
            [      0       0       0]
            [    2*x       0       0]
            sage: M._echelon_in_place("classical")
            (0, 1, 2)
            sage: M
            [1 0 0]
            [0 1 0]
            [0 0 1]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        return self._echelon_in_place_classical()

    def _echelon_in_place_classical(self, bint reduced=True, **kwds):
        """
        Change this matrix into echelon form, using classical Gaussian elimination
        and return its pivots.

        INPUT:

        - ``reduced`` -- boolean (default: ``True``); will result
          in the row-reduced echelon form (otherwise, only a
          semi-echelon form results)

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
            sage: N._echelon_in_place_classical(reduced=False)
            (2, 1, 3, 4, 5, 6, 7, 8)
            sage: N
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

            sage: type(M)
            <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: MS = MatrixSpace(M.base_ring(), M.nrows(), M.ncols(), implementation='generic')
            sage: X = MS(M)
            sage: type(X)
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
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
            return ()
        self._cache = None
        sig_on()
        try:
            MatEchelonize(self.Data)
        finally:
            sig_off()
        # Now, self.Data is in semi-echelon form.
        r = self.Data.Nor
        cdef size_t i, j, pos
        cdef PTR old, dest, src
        cdef FEL piv
        self.cache('rank', r)
        # Next, we do permutations to achieve the reduced echelon form,
        # if requested.
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
                    sig_check()
                sig_free(old)
                self.Data.Nor = self._nrows
            # Now, the pivot columns are strictly increasing.
            # We now normalize each row, and annulate everything
            # above the pivot (currently, we only know that the matrix
            # is zero below the pivots).
            for i in range(r):
                src = MatGetPtr(self.Data, i)
                piv = FfExtract(src, self.Data.PivotTable[i])
                assert piv!=FF_ZERO
                if piv != FF_ONE:
                    FfMulRow(src, mtx_tmultinv[piv])
                sig_check()
                for j in range(i):
                    dest = MatGetPtr(self.Data, j)
                    piv = FfExtract(dest, self.Data.PivotTable[i])
                    if piv != FF_ONE:
                        FfAddMulRow(dest, src, mtx_taddinv[piv])
                    else:
                        FfSubRow(dest, src)
                    sig_check()
        elif self.Data.Nor < self._nrows:
            # Some rows may have vanished. In SageMath, we
            # want that the number of rows does not change,
            # thus, we have to append zero rows.
            self.Data.Data = <PTR>check_realloc(self.Data.Data, FfCurrentRowSize*self._nrows)
            memset(self.Data.Data + FfCurrentRowSize*self.Data.Nor, FF_ZERO, FfCurrentRowSize*(self._nrows-self.Data.Nor))
            self.Data.Nor = self._nrows
        self.cache('pivots', tuple(self.Data.PivotTable[i] for i in range(r)))
        self.cache('in_echelon_form',True)
        return self._cache['pivots']


from sage.misc.superseded import deprecation_cython as deprecation


def mtx_unpickle(f, int nr, int nc, data, bint m):
    r"""
    Helper function for unpickling.

    EXAMPLES::

        sage: K.<x> = GF(9)
        sage: M = MatrixSpace(K,10,10).random_element()
        sage: M == loads(dumps(M))   # indirect doctest
        True
        sage: M is loads(dumps(M))
        False

    We also test pickles with zero rows and columns, as they may constitute
    corner cases. Note that in the following case, if ``sizeof(long)==8``,
    two matrix entries are stored in one byte, and therefore the last byte of
    a row is only half filled::

        sage: M = matrix(K,3,5, [x, 0, 1, 0,-1, 0, 0, 0, 0, 0, -1, 0, 1, 0, x])
        sage: loads(dumps(M)) == M
        True
        sage: M = matrix(K,3,5, [0, 1, 0,-1, 0, x, 0, 1, 0, -1, 0, 0, 0, 0, 0])
        sage: loads(dumps(M)) == M
        True

    TESTS:

    We test that a pickle created by one machine can be understood
    by other machines with different architecture (see :issue:`23411`).
    Internally, a row is stored in a memory block of length a multiple
    of ``sizeof(long)``, which may be machine dependent, but in a pickle,
    only the bytes actually containing data of the row are stored, which
    is machine independent. We chose a matrix over the field with `13` elements.
    Since `13^2<255<13^3`, two columns will be stored in one byte. Our matrix
    has five columns, thus, one row will occupy three bytes in the pickle,
    but eight bytes (if ``sizeof(long)==8``) in memory, and the pickle
    string will be six bytes, since we have two rows::

        sage: s = b'Uq\x82\xa7\x8bh'
        sage: len(s)
        6
        sage: from sage.matrix.matrix_gfpn_dense import mtx_unpickle, Matrix_gfpn_dense
        sage: MS = MatrixSpace(GF(13), 2, 5, implementation=Matrix_gfpn_dense)
        sage: N = mtx_unpickle(MS, 2, 5, s, True)
        sage: N
        [ 6  7  8  9 10]
        [12 11 10  9  8]
        sage: type(N)
        <class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>

    We demonstrate that a slightly different pickle format can be understood
    as well, that was at some point used by some optional package::

        sage: N == mtx_unpickle(int(13), 2, 5, s, True)
        True

    In a previous version of this optional module, the whole memory chunk
    used to store the matrix was stored. The result would have been, as
    in the following example, a string of length 16. Unpickling works, but
    results in a warning::

        sage: t = b'Uq\x82\x00\x00\x00\x00\x00\xa7\x8bh\x00\x00\x00\x00\x00'
        sage: len(t)
        16
        sage: N == mtx_unpickle(MS, 2, 5, t, True)
        doctest:warning
        ...
        DeprecationWarning: Reading this pickle may be machine dependent
        See https://github.com/sagemath/sage/issues/23411 for details.
        True

    Unpickling would even work in the case that the machine creating
    the deprecated pickle had ``sizeof(long)==9``::

        sage: t = b'Uq\x82\x00\x00\x00\x00\x00\x00\xa7\x8bh\x00\x00\x00\x00\x00\x00'
        sage: len(t)
        18
        sage: N == mtx_unpickle(MS, 2, 5, t, True)
        True

    The data may be empty, which results in the zero matrix::

        sage: mtx_unpickle(MS, 2, 5, b'', True)
        [0 0 0 0 0]
        [0 0 0 0 0]

    We test further corner cases. A :exc:`ValueError` is raised if the number
    of bytes in the pickle does not comply with either the old or the new
    pickle format (we test several code paths here)::

        sage: t = b'Uq\x82\x00\x00\x00\x00\x00\xa7\x8bh\x00\x00\x00\x00\x00\x00'
        sage: mtx_unpickle(MS, 2, 5, t, True)
        Traceback (most recent call last):
        ...
        ValueError: Expected a pickle with 3*2 bytes, got 17 instead
        sage: t = b'Uq\x82\x00\x00\x00\x00\x00\xa7\x8bh\x00\x00\x00\x00\x00\x00'
        sage: mtx_unpickle(MS, 2, 5, t[:4], True)
        Traceback (most recent call last):
        ...
        ValueError: Expected a pickle with 3*2 bytes, got 2*2 instead
        sage: MS = MatrixSpace(GF(13), 0, 5)
        sage: mtx_unpickle(MS, 0, 5, s, True)
        Traceback (most recent call last):
        ...
        ValueError: This matrix pickle contains data, thus, the number of rows
        and columns must be positive
        sage: MS = MatrixSpace(GF(13), 3, 5)
        sage: mtx_unpickle(MS, 2, 5, s, True)
        Traceback (most recent call last):
        ...
        ValueError: Inconsistent dimensions in this matrix pickle
        sage: mtx_unpickle(MatrixSpace(GF(19),0,5), 0, 5, b'', True)
        []
    """
    # The expected input type is bytes. However, Python-2 legacy pickles do
    # not distinguish between str and bytes. If such pickle is unpickled
    # in Python-3, Sage will receive a str in `latin1` encoding. Therefore,
    # in the following line, we use a helper function that would return bytes,
    # regardless whether the input is bytes or str.
    cdef bytes Data = str_to_bytes(data, encoding='latin1')
    if isinstance(f, int):
        # This is for old pickles created with the group cohomology spkg
        MS = MatrixSpace(GF(f, 'z'), nr, nc, implementation=Matrix_gfpn_dense)
    else:
        MS = f
        if MS.nrows() != nr or MS.ncols() != nc:
            raise ValueError("Inconsistent dimensions in this matrix pickle")
        f = MS.base_ring().order()

    OUT = <Matrix_gfpn_dense>Matrix_gfpn_dense.__new__(Matrix_gfpn_dense, MS)
    sig_on()
    try:
        OUT.Data = MatAlloc(f, nr, nc)
    finally:
        sig_off()
    OUT._is_immutable = not m
    OUT._converter = FieldConverter(OUT._base_ring)
    cdef char *x
    cdef PTR pt
    cdef size_t lenData = len(Data)
    cdef size_t pickled_rowsize
    cdef size_t i
    if Data:
        if nr <= 0 or nc <= 0:
            raise ValueError("This matrix pickle contains data, thus, the number of rows and columns must be positive")
        pickled_rowsize = lenData//nr
        if lenData != pickled_rowsize*nr:
            raise ValueError(f"Expected a pickle with {FfCurrentRowSizeIo}*{nr} bytes, got {lenData} instead")
        x = PyBytes_AsString(Data)
        if pickled_rowsize == FfCurrentRowSizeIo:
            pt = OUT.Data.Data
            for i in range(nr):
                memcpy(pt,x,FfCurrentRowSizeIo)
                x += FfCurrentRowSizeIo
                FfStepPtr(&(pt))
                sig_check()
        elif pickled_rowsize >= FfCurrentRowSizeIo:
            deprecation(23411, "Reading this pickle may be machine dependent")
            if pickled_rowsize == FfCurrentRowSize:
                memcpy(OUT.Data.Data, x, OUT.Data.RowSize*OUT.Data.Nor)
            else:
                pt = OUT.Data.Data
                for i in range(nr):
                    memcpy(pt,x,FfCurrentRowSizeIo)
                    x += pickled_rowsize
                    FfStepPtr(&(pt))
                    sig_check()
        else:
            raise ValueError(f"Expected a pickle with {FfCurrentRowSizeIo}*{nr} bytes, got {pickled_rowsize}*{nr} instead")
    return OUT
