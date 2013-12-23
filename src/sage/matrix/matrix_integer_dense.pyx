"""
Dense matrices over the integer ring

AUTHORS:

- William Stein

- Robert Bradshaw

EXAMPLES::

    sage: a = matrix(ZZ, 3,3, range(9)); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: a.det()
    0
    sage: a[0,0] = 10; a.det()
    -30
    sage: a.charpoly()
    x^3 - 22*x^2 + 102*x + 30
    sage: b = -3*a
    sage: a == b
    False
    sage: b < a
    True

TESTS::

    sage: a = matrix(ZZ,2,range(4), sparse=False)
    sage: TestSuite(a).run()
    sage: Matrix(ZZ,0,0).inverse()
    []
"""

########## *** IMPORTANT ***
# If you're working on this code, we *always* assume that
#   self._matrix[i] = self._entries[i*self._ncols]
# !!!!!!!! Do not break this!
# This is assumed in the _rational_kernel_iml

######################################################################
#       Copyright (C) 2006,2007 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
######################################################################

from sage.modules.vector_integer_dense cimport Vector_integer_dense

from sage.misc.misc import verbose, get_verbose, cputime

from sage.rings.arith import previous_prime
from sage.structure.element import is_Element
from sage.structure.proof.proof import get_flag as get_proof_flag

from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense

#########################################################
# PARI C library
from sage.libs.pari.gen cimport gen
from sage.libs.pari.pari_instance cimport PariInstance

import sage.libs.pari.pari_instance
cdef PariInstance pari = sage.libs.pari.pari_instance.pari

include "sage/libs/pari/decl.pxi"
include "sage/libs/pari/pari_err.pxi"

cdef extern from "convert.h":
    cdef void t_INT_to_ZZ( mpz_t value, long *g )

#########################################################

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/gmp.pxi"
include "sage/ext/random.pxi"

cdef extern from "math.h":
    double log(double x)


from sage.ext.multi_modular import MultiModularBasis
from sage.ext.multi_modular cimport MultiModularBasis

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.integer_ring import ZZ, IntegerRing_class
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence

from matrix_modn_dense_float cimport Matrix_modn_dense_template
from matrix_modn_dense_float cimport Matrix_modn_dense_float
from matrix_modn_dense_double cimport Matrix_modn_dense_double

from matrix_mod2_dense import Matrix_mod2_dense
from matrix_mod2_dense cimport Matrix_mod2_dense

from matrix_modn_dense cimport is_Matrix_modn_dense


from matrix2 import decomp_seq

import sage.modules.free_module
import sage.modules.free_module_element

from matrix cimport Matrix

cimport sage.structure.element

import sage.matrix.matrix_space as matrix_space

################
# Used for modular HNF
from sage.rings.fast_arith cimport arith_int
cdef arith_int ai
ai = arith_int()
################

######### linbox interface ##########
from sage.libs.linbox.linbox cimport Linbox_integer_dense
cdef Linbox_integer_dense linbox
linbox = Linbox_integer_dense()
USE_LINBOX_POLY = True


########## iml -- integer matrix library ###########

cdef extern from "iml.h":

    enum SOLU_POS:
        LeftSolu = 101
        RightSolu = 102

    long nullspaceMP(long n, long m,
                     mpz_t *A, mpz_t * *mp_N_pass)

    void nonsingSolvLlhsMM (SOLU_POS solupos, long n, \
                       long m, mpz_t *mp_A, mpz_t *mp_B, mpz_t *mp_N, \
                       mpz_t mp_D)



cdef class Matrix_integer_dense(matrix_dense.Matrix_dense):   # dense or sparse
    r"""
    Matrix over the integers.

    On a 32-bit machine, they can have at most `2^{32}-1` rows or
    columns.  On a 64-bit machine, matrices can have at most
    `2^{64}-1` rows or columns.

    EXAMPLES::

        sage: a = MatrixSpace(ZZ,3)(2); a
        [2 0 0]
        [0 2 0]
        [0 0 2]
        sage: a = matrix(ZZ,1,3, [1,2,-3]); a
        [ 1  2 -3]
        sage: a = MatrixSpace(ZZ,2,4)(2); a
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square
    """
    ########################################################################
    # LEVEL 1 functionality
    # x * __cinit__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################

    def __cinit__(self, parent, entries, coerce, copy):
        """
        Create and allocate memory for the matrix. Does not actually
        initialize any of the memory.

        INPUT:


        -  ``parent, entries, coerce, copy`` - as for
           __init__.


        EXAMPLES::

            sage: from sage.matrix.matrix_integer_dense import Matrix_integer_dense
            sage: a = Matrix_integer_dense.__new__(Matrix_integer_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        .. warning::

           This is for internal use only, or if you really know what
           you're doing.
        """
        self._parent = parent
        self._base_ring = ZZ
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._pivots = None

        # Allocate an array where all the entries of the matrix are stored.
        sig_on()
        self._entries = <mpz_t *>sage_malloc(sizeof(mpz_t) * (self._nrows * self._ncols))
        sig_off()
        if self._entries == NULL:
            raise MemoryError("out of memory allocating a matrix")

        # Allocate an array of pointers to the rows, which is useful for
        # certain algorithms.
        ##################################
        # *IMPORTANT*: FOR MATRICES OVER ZZ, WE ALWAYS ASSUME THAT
        # THIS ARRAY IS *not* PERMUTED.  This should be OK, since all
        # algorithms are multi-modular.
        ##################################
        self._matrix = <mpz_t **> sage_malloc(sizeof(mpz_t*)*self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError("out of memory allocating a matrix")

        # Set each of the pointers in the array self._matrix to point
        # at the memory for the corresponding row.
        cdef Py_ssize_t i, k
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    cdef _init_linbox(self):
        sig_on()
        linbox.set(self._matrix, self._nrows, self._ncols)
        sig_off()

    def __copy__(self):
        r"""
        Returns a new copy of this matrix.

        EXAMPLES::

            sage: a = matrix(ZZ,1,3, [1,2,-3]); a
            [ 1  2 -3]
            sage: b = a.__copy__(); b
            [ 1  2 -3]
            sage: b is a
            False
            sage: b == a
            True
        """
        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent,
                                         0, 0, 0)
        cdef Py_ssize_t i
        sig_on()
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init_set(A._entries[i], self._entries[i])
        sig_off()
        A._initialized = True
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    def __hash__(self):
        r"""
        Returns hash of self.

        self must be immutable.

        EXAMPLES::

            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        ::

            sage: a.set_immutable()
            sage: hash(a)
            8
        """
        return self._hash()

    def __dealloc__(self):
        """
        Frees all the memory allocated for this matrix. EXAMPLE::

            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: del a
        """
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < (self._nrows * self._ncols):
                mpz_clear(self._entries[i])
            sage_free(self._entries)
            sage_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        r"""
        Initialize a dense matrix over the integers.

        INPUT:


        -  ``parent`` - a matrix space

        -  ``entries`` - list - create the matrix with those
           entries along the rows.

        -  ``other`` - a scalar; entries is coerced to an
           integer and the diagonal entries of this matrix are set to that
           integer.

        -  ``coerce`` - whether need to coerce entries to the
           integers (program may crash if you get this wrong)

        -  ``copy`` - ignored (since integers are immutable)


        EXAMPLES:

        The __init__ function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a `2 \times 2` and a `1\times 4` matrix::

            sage: matrix(ZZ,2,2,range(4))
            [0 1]
            [2 3]
            sage: Matrix(ZZ,1,4,range(4))
            [0 1 2 3]

        If the number of columns isn't given, it is determined from the
        number of elements in the list.

        ::

            sage: matrix(ZZ,2,range(4))
            [0 1]
            [2 3]
            sage: matrix(ZZ,2,range(6))
            [0 1 2]
            [3 4 5]

        Another way to make a matrix is to create the space of matrices and
        coerce lists into it.

        ::

            sage: A = Mat(ZZ,2); A
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: A(range(4))
            [0 1]
            [2 3]

        Actually it is only necessary that the input can be coerced to a
        list, so the following also works::

            sage: v = reversed(range(4)); type(v)
            <type 'listreverseiterator'>
            sage: A(v)
            [3 2]
            [1 0]

        Matrices can have many rows or columns (in fact, on a 64-bit
        machine they could have up to `2^64-1` rows or columns)::

            sage: v = matrix(ZZ,1,10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over Integer Ring
        """
        cdef Py_ssize_t i, j
        cdef bint is_list
        cdef Integer x

        if entries is None:
            x = ZZ(0)
            is_list = 0
        elif isinstance(entries, (int,long)) or is_Element(entries):
            try:
                x = ZZ(entries)
            except TypeError:
                self._initialized = False
                raise TypeError("unable to coerce entry to an integer")
            is_list = 0
        else:
            entries = list(entries)
            is_list = 1

        if is_list:

            # Create the matrix whose entries are in the given entry list.
            if len(entries) != self._nrows * self._ncols:
                sage_free(self._entries)
                sage_free(self._matrix)
                self._entries = NULL
                raise TypeError("entries has the wrong length")
            if coerce:
                for i from 0 <= i < self._nrows * self._ncols:
                    x = ZZ(entries[i])
                    # todo -- see integer.pyx and the TODO there; perhaps this could be
                    # sped up by creating a mpz_init_set_sage function.
                    mpz_init_set(self._entries[i], x.value)
                self._initialized = True
            else:
                for i from 0 <= i < self._nrows * self._ncols:
                    mpz_init_set(self._entries[i], (<Integer> entries[i]).value)
                self._initialized = True
        else:

            # If x is zero, make the zero matrix and be done.
            if mpz_sgn(x.value) == 0:
                self._zero_out_matrix()
                return

            # the matrix must be square:
            if self._nrows != self._ncols:
                sage_free(self._entries)
                sage_free(self._matrix)
                self._entries = NULL
                raise TypeError("nonzero scalar matrix must be square")

            # Now we set all the diagonal entries to x and all other entries to 0.
            self._zero_out_matrix()
            j = 0
            for i from 0 <= i < self._nrows:
                mpz_set(self._entries[j], x.value)
                j = j + self._nrows + 1
            self._initialized = True

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set position i,j of this matrix to x.

        (VERY UNSAFE - value *must* be of type Integer).

        INPUT:


        -  ``i`` - row

        -  ``j`` - column

        -  ``value`` - The value to set self[i,j] to. value
           MUST be of type Integer


        EXAMPLES::

            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        mpz_set(self._matrix[i][j], (<Integer>value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Returns (j, i) entry of self as a new Integer.

        .. warning::

           This is very unsafe; it assumes i and j are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            6
        """
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, self._matrix[i][j])
        return z

    def _pickle(self):
        """
        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1,193,15,-2,3,0])
            sage: a._pickle()
            ('1 61 f -2 3 0', 0)

            sage: S = ModularSymbols(250,4,sign=1).cuspidal_submodule().new_subspace().decomposition() # long time
            sage: S == loads(dumps(S)) # long time
            True
        """
        return self._pickle_version0(), 0

    cdef _pickle_version0(self):
        """
        EXAMPLES::

            sage: matrix(ZZ,1,3,[1,193,15])._pickle()   # indirect doctest
            ('1 61 f', 0)

        """
        return self._export_as_string(32)

    cpdef _export_as_string(self, int base=10):
        """
        Return space separated string of the entries in this matrix, in the
        given base. This is optimized for speed.

        INPUT: base -an integer = 36; (default: 10)

        EXAMPLES::

            sage: m = matrix(ZZ,2,3,[1,2,-3,1,-2,-45])
            sage: m._export_as_string(10)
            '1 2 -3 1 -2 -45'
            sage: m._export_as_string(16)
            '1 2 -3 1 -2 -2d'
        """
        # TODO: *maybe* redo this to use mpz_import and mpz_export
        # from sec 5.14 of the GMP manual. ??
        cdef int i, j, len_so_far, m, n
        cdef char *a
        cdef char *s, *t, *tmp

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*10
            s = <char*> sage_malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            sig_on()
            for i from 0 <= i < self._nrows * self._ncols:
                m = mpz_sizeinbase (self._entries[i], base)
                if len_so_far + m + 2 >= n:
                    # copy to new string with double the size
                    n = 2*n + m + 1
                    tmp = <char*> sage_malloc(n * sizeof(char))
                    strcpy(tmp, s)
                    sage_free(s)
                    s = tmp
                    t = s + len_so_far
                #endif
                mpz_get_str(t, base, self._entries[i])
                m = strlen(t)
                len_so_far = len_so_far + m + 1
                t = t + m
                t[0] = <char>32
                t[1] = <char>0
                t = t + 1
            sig_off()
            data = str(s)[:-1]
            sage_free(s)
        return data

    def _unpickle(self, data, int version):
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError("unknown matrix version (=%s)"%version)

    cdef _unpickle_version0(self, data):
        cdef Py_ssize_t i, n
        data = data.split()
        n = self._nrows * self._ncols
        if len(data) != n:
            raise RuntimeError("invalid pickle data.")
        for i from 0 <= i < n:
            s = data[i]
            if mpz_init_set_str(self._entries[i], s, 32):
                raise RuntimeError("invalid pickle data")
        self._initialized = True


    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    cdef _zero_out_matrix(self):
        """
        Set this matrix to be the zero matrix. This is only for internal
        use. (Note: this matrix must NOT already have initialised
        entries.)
        """
        # TODO: This is about 6-10 slower than MAGMA doing what seems to be the same thing.
        # Moreover, NTL can also do this quickly.  Why?   I think both have specialized
        # small integer classes. (dmharvey: yes, NTL does not allocate any memory when
        # initializing a ZZ to zero.)
        sig_on()
        cdef Py_ssize_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(self._entries[i])
        sig_off()
        self._initialized = True

    cdef _new_unitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the integers with the given number of rows
        and columns. All memory is allocated for this matrix, but its
        entries have not yet been filled in.
        """
        P = self._parent.matrix_space(nrows, ncols)
        return Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)


    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_
    # x * cdef _sub_
    # x * cdef _mul_
    # x * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    # cdef _mul_(self, Matrix right):
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):

    def __nonzero__(self):
        r"""
        Tests whether self is the zero matrix.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ, 2, 3)(range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.__nonzero__()
            True
            sage: (a - a).__nonzero__()
            False

        ::

            sage: a = MatrixSpace(ZZ, 0, 3)()
            sage: a.__nonzero__()
            False
            sage: a = MatrixSpace(ZZ, 3, 0)()
            sage: a.__nonzero__()
            False
            sage: a = MatrixSpace(ZZ, 0, 0)()
            sage: a.__nonzero__()
            False
        """
        cdef mpz_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_sgn(self._entries[i]):
                return True
        return False

    def _multiply_linbox(self, Matrix right):
        """
        Multiply matrices over ZZ using linbox.

        .. warning::

           This is very slow right now, i.e., linbox is very slow.

        EXAMPLES::

            sage: A = matrix(ZZ,2,3,range(6))
            sage: A*A.transpose()
            [ 5 14]
            [14 50]
            sage: A._multiply_linbox(A.transpose())
            [ 5 14]
            [14 50]
        """
        cdef int e
        cdef Matrix_integer_dense ans, B
        B = right
        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        self._init_linbox()
        sig_on()
        linbox.matrix_matrix_multiply(ans._matrix, B._matrix, B._nrows, B._ncols)
        sig_off()
        return ans

    def _multiply_classical(self, Matrix right):
        """
        EXAMPLE::

            sage: n = 3
            sage: a = MatrixSpace(ZZ,n,n)(range(n^2))
            sage: b = MatrixSpace(ZZ,n,n)(range(1, n^2 + 1))
            sage: a._multiply_classical(b)
            [ 18  21  24]
            [ 54  66  78]
            [ 90 111 132]
        """
        if self._ncols != right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")

        cdef Py_ssize_t i, j, k, l, nr, nc, snc
        cdef mpz_t *v
        cdef object parent

        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        if self._nrows == right._nrows:
            # self acts on the space of right
            parent = right.parent()
        if self._ncols == right._ncols:
            # right acts on the space of self
            parent = self.parent()
        else:
            parent = self.matrix_space(nr, nc)

        cdef Matrix_integer_dense M, _right
        _right = right

        M = Matrix_integer_dense.__new__(Matrix_integer_dense, parent, None, None, None)
        Matrix.__init__(M, parent)

        # M has memory allocated but entries are not initialized
        cdef mpz_t *entries
        entries = M._entries

        cdef mpz_t s, z
        mpz_init(s)
        mpz_init(z)

        sig_on()
        l = 0
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpz_set_si(s,0)   # set s = 0
                v = self._matrix[i]
                for k from 0 <= k < snc:
                    mpz_mul(z, v[k], _right._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init_set(entries[l], s)
                l += 1
        sig_off()
        mpz_clear(s)
        mpz_clear(z)
        M._initialized = True
        return M

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):

        #############
        # see the tune_multiplication function below.
        n = max(self._nrows, self._ncols, right._nrows, right._ncols)
        if n <= 20:
            return self._multiply_classical(right)
        a = self.height(); b = right.height()
        if float(max(a,b)) / float(n) >= 0.70:
            return self._multiply_classical(right)
        else:
            return self._multiply_multi_modular(right)

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i
        cdef Integer _x
        _x = Integer(right)
        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(M._entries[i])
            mpz_mul(M._entries[i], self._entries[i], _x.value)
        M._initialized = True
        return M

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two dense matrices over ZZ.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9))
            sage: a+a
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]
            sage: b = MatrixSpace(ZZ,3)(range(9))
            sage: b.swap_rows(1,2)
            sage: a+b
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        sig_on()
        cdef mpz_t *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_integer_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_init(row_ans[j])
                mpz_add(row_ans[j], row_self[j], row_right[j])
        sig_off()
        M._initialized = True
        return M

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract two dense matrices over ZZ.

        EXAMPLES::

            sage: M = Mat(ZZ,3)
            sage: a = M(range(9)); b = M(reversed(range(9)))
            sage: a - b
            [-8 -6 -4]
            [-2  0  2]
            [ 4  6  8]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        sig_on()
        cdef mpz_t *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_integer_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_init(row_ans[j])
                mpz_sub(row_ans[j], row_self[j], row_right[j])
        sig_off()
        M._initialized = True
        return M


    cdef int _cmp_c_impl(self, Element right) except -2:
        r"""
        Compares self with right, examining entries in lexicographic (row
        major) ordering.

        EXAMPLES::

            sage: Matrix(ZZ, [[0, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 10], [20, 30]]))
            0
            sage: Matrix(ZZ, [[0, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 15], [20, 30]]))
            -1
            sage: Matrix(ZZ, [[5, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 15], [20, 30]]))
            1
            sage: Matrix(ZZ, [[5, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 10], [25, 30]]))
            1
        """
        cdef mpz_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k

        for i from 0 <= i < self._nrows:
            a = self._matrix[i]
            b = (<Matrix_integer_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                k = mpz_cmp(a[j], b[j])
                if k:
                    if k < 0:
                        return -1
                    else:
                        return 1
        return 0


    cdef Vector _vector_times_matrix_(self, Vector v):
        """
        Returns the vector times matrix product.

        INPUT:


        -  ``v`` - a free module element.


        OUTPUT: The vector times matrix product v\*A.

        EXAMPLES::

            sage: B = matrix(ZZ,2, [1,2,3,4])
            sage: V = ZZ^2
            sage: w = V([-1,5])
            sage: w*B
            (14, 18)
        """
        cdef Vector_integer_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpz_t x

        M = self._row_ambient_module()
        w = <Vector_integer_dense> v
        ans = M.zero_vector()

        mpz_init(x)
        for i from 0 <= i < self._ncols:
            mpz_set_si(x, 0)
            for j from 0 <= j < self._nrows:
                mpz_addmul(x, w._entries[j], self._matrix[j][i])
            mpz_set(ans._entries[i], x)
        mpz_clear(x)
        return ans


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #    * Specialized echelon form
    ########################################################################

    def _clear_denom(self):
        """
        INPUT:

        -  ``self`` - a matrix

        OUTPUT:  self, 1

        EXAMPLES::

            sage: a = matrix(ZZ,2,[1,2,3,4])
            sage: a._clear_denom()
            (
            [1 2]
            [3 4], 1
            )
        """
        return self, ZZ(1)

    def charpoly(self, var='x', algorithm='linbox'):
        """
        INPUT:


        -  ``var`` - a variable name

        -  ``algorithm`` - 'linbox' (default) 'generic'


        .. note::

           Linbox charpoly disabled on 64-bit machines, since it hangs
           in many cases.

        EXAMPLES::

            sage: A = matrix(ZZ,6, range(36))
            sage: f = A.charpoly(); f
            x^6 - 105*x^5 - 630*x^4
            sage: f(A) == 0
            True
            sage: n=20; A = Mat(ZZ,n)(range(n^2))
            sage: A.charpoly()
            x^20 - 3990*x^19 - 266000*x^18
            sage: A.minpoly()
            x^3 - 3990*x^2 - 266000*x

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:trac:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(ZZ, 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: A.charpoly('x')
            x^2 - 3*x - 2
            sage: A.charpoly('y')
            y^2 - 3*y - 2
            sage: A._cache['charpoly_linbox']
            x^2 - 3*x - 2

        """
        cache_key = 'charpoly_%s' % algorithm
        g = self.fetch(cache_key)
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'linbox' and not USE_LINBOX_POLY:
            algorithm = 'generic'

        if algorithm == 'linbox':
            g = self._charpoly_linbox(var)
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.charpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

        self.cache(cache_key, g)
        return g

    def minpoly(self, var='x', algorithm='linbox'):
        """
        INPUT:


        -  ``var`` - a variable name

        -  ``algorithm`` - 'linbox' (default) 'generic'


        .. note::

           Linbox charpoly disabled on 64-bit machines, since it hangs
           in many cases.

        EXAMPLES::

            sage: A = matrix(ZZ,6, range(36))
            sage: A.minpoly()
            x^3 - 105*x^2 - 630*x
            sage: n=6; A = Mat(ZZ,n)([k^2 for k in range(n^2)])
            sage: A.minpoly()
            x^4 - 2695*x^3 - 257964*x^2 + 1693440*x
        """
        key = 'minpoly_%s_%s'%(algorithm, var)
        x = self.fetch(key)
        if x: return x

        if algorithm == 'linbox' and not USE_LINBOX_POLY:
            algorithm = 'generic'

        if algorithm == 'linbox':
            g = self._minpoly_linbox(var)
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.minpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

        self.cache(key, g)
        return g

    def _minpoly_linbox(self, var='x'):
        return self._poly_linbox(var=var, typ='minpoly')

    def _charpoly_linbox(self, var='x'):
        if self.is_zero():  # program around a bug in linbox on 32-bit linux
            x = self.base_ring()[var].gen()
            return x ** self._nrows
        return self._poly_linbox(var=var, typ='charpoly')

    def _poly_linbox(self, var='x', typ='minpoly'):
        """
        INPUT:


        -  ``var`` - 'x'

        -  ``typ`` - 'minpoly' or 'charpoly'
        """
        time = verbose('computing %s of %s x %s matrix using linbox'%(typ, self._nrows, self._ncols))
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")
        if self._nrows <= 1:
            return matrix_dense.Matrix_dense.charpoly(self, var)
        self._init_linbox()
        if typ == 'minpoly':
            sig_on()
            v = linbox.minpoly()
            sig_off()
        else:
            sig_on()
            v = linbox.charpoly()
            sig_off()
        R = self._base_ring[var]
        verbose('finished computing %s'%typ, time)
        return R(v)


    def height(self):
        """
        Return the height of this matrix, i.e., the max absolute value of
        the entries of the matrix.

        OUTPUT: A nonnegative integer.

        EXAMPLE::

            sage: a = Mat(ZZ,3)(range(9))
            sage: a.height()
            8
            sage: a = Mat(ZZ,2,3)([-17,3,-389,15,-1,0]); a
            [ -17    3 -389]
            [  15   -1    0]
            sage: a.height()
            389
        """
        cdef mpz_t h
        cdef Integer x

        self.mpz_height(h)
        x = Integer()
        x.set_from_mpz(h)
        mpz_clear(h)

        return x

    cdef int mpz_height(self, mpz_t height) except -1:
        """
        Used to compute the height of this matrix.

        INPUT:


        -  ``height`` - a GMP mpz_t (that has not been
           initialized!)


        OUTPUT: sets the value of height to the height of this matrix,
        i.e., the max absolute value of the entries of the matrix.
        """
        cdef mpz_t x, h
        cdef Py_ssize_t i

        sig_on()
        mpz_init_set_si(h, 0)
        mpz_init(x)

        for i from 0 <= i < self._nrows * self._ncols:
            mpz_abs(x, self._entries[i])
            if mpz_cmp(h, x) < 0:
                mpz_set(h, x)

        mpz_init_set(height, h)
        mpz_clear(h)
        mpz_clear(x)
        sig_off()

        return 0   # no error occurred.

    cpdef double _log_avg_sq1(self) except -1.0:
        """
        Return the logarithm of the average of `x^2 + 1`, where `x`
        ranges over the matrix entries.

        This is used to determine which determinant algorithm to use.

        TESTS::

            sage: M = random_matrix(ZZ,100)
            sage: L1 = M._log_avg_sq1()
            sage: L2 = log(RR(sum([i*i+1 for i in M.list()])/10000))
            sage: abs(L1 - L2) < 1e-13
            True
            sage: matrix(ZZ,10)._log_avg_sq1()
            0.0
        """
        cdef Py_ssize_t i
        cdef Py_ssize_t N = self._nrows * self._ncols

        # wsq * 4^wsq_e = sum of entries squared plus number of entries
        cdef double wsq = N
        cdef long wsq_e = 0

        cdef double d
        cdef long e

        sig_on()
        for i in range(N):
            d = mpz_get_d_2exp(&e, self._entries[i])
            if (e > wsq_e):
                wsq = ldexp(wsq, 2*(wsq_e - e))
                wsq_e = e
            elif (e < wsq_e):
                d = ldexp(d, e - wsq_e)
            wsq += d*d
        sig_off()

        # Compute log(wsq * 4^wsq_e / N) =
        # log(wsq * 4^wsq_e / N) = log(wsq/N) + wsq_e * log(4)
        return log(wsq/N) + wsq_e * log(4.0)

    def _multiply_multi_modular(left, Matrix_integer_dense right):
        """
        Multiply this matrix by ``left`` using a multi modular algorithm.

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 3, range(5,11))
            sage: N = Matrix(ZZ, 3, 2, range(15,21))
            sage: M._multiply_multi_modular(N)
            [310 328]
            [463 490]
            sage: M._multiply_multi_modular(-N)
            [-310 -328]
            [-463 -490]
        """
        cdef Integer h
        cdef mod_int *moduli
        cdef int i, n, k

        h = left.height() * right.height() * left.ncols()
        verbose('multiplying matrices of height %s and %s'%(left.height(),right.height()))
        mm = MultiModularBasis(h)
        res = left._reduce(mm)
        res_right = right._reduce(mm)
        k = len(mm)
        for i in range(k):  # yes, I could do this with zip, but to conserve memory...
            t = cputime()
            res[i] *= res_right[i]
            verbose('multiplied matrices modulo a prime (%s/%s)'%(i+1,k), t)
        result = left.new_matrix(left.nrows(), right.ncols())
        _lift_crt(result, res, mm)  # changes result
        return result

    cdef void reduce_entry_unsafe(self, Py_ssize_t i, Py_ssize_t j, Integer modulus):
        # Used for p-adic matrices.
        if mpz_cmp(self._matrix[i][j], modulus.value) >= 0 or mpz_cmp_ui(self._matrix[i][j], 0) < 0:
            mpz_mod(self._matrix[i][j], self._matrix[i][j], modulus.value)

    def _mod_int(self, modulus):
        """
        Reduce the integer matrix modulo a positive integer.

        EXAMPLES::
            sage: matrix(QQ,2,[1,0,0,1]).change_ring(GF(2)) - 1
            [0 0]
            [0 0]

        """
        cdef mod_int c = modulus
        if int(c) != modulus:
            raise OverflowError
        else:
            return self._mod_int_c(modulus)

    cdef _mod_two(self):
        cdef Matrix_mod2_dense res
        res = Matrix_mod2_dense.__new__(Matrix_mod2_dense, matrix_space.MatrixSpace(IntegerModRing(2), self._nrows, self._ncols, sparse=False), None, None, None)
        res.__init__(matrix_space.MatrixSpace(IntegerModRing(2), self._nrows, self._ncols, sparse=False), self.list(), None, None)
        return res

    cdef _mod_int_c(self, mod_int p):
        from matrix_modn_dense_float import MAX_MODULUS as MAX_MODULUS_FLOAT
        from matrix_modn_dense_double import MAX_MODULUS as MAX_MODULUS_DOUBLE

        cdef Py_ssize_t i, j
        cdef mpz_t* self_row

        cdef float* res_row_f
        cdef Matrix_modn_dense_float res_f

        cdef double* res_row_d
        cdef Matrix_modn_dense_double res_d

        if p == 2:
            return self._mod_two()
        elif p < MAX_MODULUS_FLOAT:
            res_f = Matrix_modn_dense_float.__new__(Matrix_modn_dense_float,
                                                    matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None)
            for i from 0 <= i < self._nrows:
                self_row = self._matrix[i]
                res_row_f = res_f._matrix[i]
                for j from 0 <= j < self._ncols:
                    res_row_f[j] = <float>mpz_fdiv_ui(self_row[j], p)
            return res_f

        elif p < MAX_MODULUS_DOUBLE:
            res_d = Matrix_modn_dense_double.__new__(Matrix_modn_dense_double,
                                                     matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None)
            for i from 0 <= i < self._nrows:
                self_row = self._matrix[i]
                res_row_d = res_d._matrix[i]
                for j from 0 <= j < self._ncols:
                    res_row_d[j] = <double>mpz_fdiv_ui(self_row[j], p)
            return res_d
        else:
            raise ValueError("p to big.")

    def _reduce(self, moduli):
        from matrix_modn_dense_float import MAX_MODULUS as MAX_MODULUS_FLOAT
        from matrix_modn_dense_double import MAX_MODULUS as MAX_MODULUS_DOUBLE

        if isinstance(moduli, (int, long, Integer)):
            return self._mod_int(moduli)
        elif isinstance(moduli, list):
            moduli = MultiModularBasis(moduli)

        cdef MultiModularBasis mm
        mm = moduli

        res = []
        for p in mm:
            if p < MAX_MODULUS_FLOAT:
                res.append( Matrix_modn_dense_float.__new__(Matrix_modn_dense_float,
                                                            matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False),
                                                            None, None, None) )
            elif p < MAX_MODULUS_DOUBLE:
                res.append( Matrix_modn_dense_double.__new__(Matrix_modn_dense_double,
                                                             matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False),
                                                             None, None, None) )
            else:
                raise ValueError("p=%d too big."%p)

        cdef size_t i, k, n
        cdef Py_ssize_t nr, nc

        n = len(mm)
        nr = self._nrows
        nc = self._ncols

        cdef mod_int **row_list
        row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
        if row_list == NULL:
            raise MemoryError("out of memory allocating multi-modular coefficient list")
        for k in range(n):
            row_list[k] = <mod_int*>sage_malloc(sizeof(mod_int)*nc)
            if row_list[k] == NULL:
                for i in range(k):
                    sage_free(row_list[i])
                sage_free(row_list)
                raise MemoryError("out of memory allocating multi-modular coefficient list")

        sig_on()
        for i from 0 <= i < nr:
            for k from 0 <= k < n:
                mm.mpz_reduce_vec(self._matrix[i], row_list, nc)
                if isinstance(res[k], Matrix_modn_dense_float):
                    for j in range(nc):
                        (<Matrix_modn_dense_float>res[k])._matrix[i][j] = (<float>row_list[k][j])%(<Matrix_modn_dense_float>res[k]).p
                else:
                    for j in range(nc):
                        (<Matrix_modn_dense_double>res[k])._matrix[i][j] = (<double>row_list[k][j])%(<Matrix_modn_dense_double>res[k]).p
        sig_off()

        for k in range(n):
            sage_free(row_list[k])
        sage_free(row_list)
        return res

    def _echelon_in_place_classical(self):
        cdef Matrix_integer_dense E
        E = self.echelon_form()

        cdef int i
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_set(self._entries[i], E._entries[i])

        self.clear_cache()

    def _echelon_strassen(self):
        raise NotImplementedError

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: m = matrix(ZZ,2,3,[1,2,-3,1,-2,-45])
            sage: m._magma_init_(magma)
            'Matrix(IntegerRing(),2,3,StringToIntegerSequence("1 2 -3 1 -2 -45"))'
            sage: magma(m)                                               # optional - magma
            [  1   2  -3]
            [  1  -2 -45]
        """
        w = self._export_as_string(base=10)
        return 'Matrix(IntegerRing(),%s,%s,StringToIntegerSequence("%s"))'%(
            self.nrows(), self.ncols(), w)

    def symplectic_form(self):
        r"""
            Find a symplectic basis for self if self is an anti-symmetric,
            alternating matrix.

            Returns a pair (F, C) such that the rows of C form a symplectic
            basis for self and F = C \* self \* C.transpose().

            Raises a ValueError if self is not anti-symmetric, or self is not
            alternating.

            Anti-symmetric means that `M = -M^t`. Alternating means
            that the diagonal of `M` is identically zero.

            A symplectic basis is a basis of the form
            `e_1, \ldots, e_j, f_1, \ldots f_j, z_1, \dots, z_k`
            such that

            -  `z_i M v^t` = 0 for all vectors `v`

            -  `e_i M {e_j}^t = 0` for all `i, j`

            -  `f_i M {f_j}^t = 0` for all `i, j`

            -  `e_i M {f_i}^t = 1` for all `i`

            -  `e_i M {f_j}^t = 0` for all `i` not equal
                `j`.

            The ordering for the factors `d_{i} | d_{i+1}` and for
            the placement of zeroes was chosen to agree with the output of
            ``smith_form``.

            See the example for a pictorial description of such a basis.

            EXAMPLES::

                sage: E = matrix(ZZ, 5, 5, [0, 14, 0, -8, -2, -14, 0, -3, -11, 4, 0, 3, 0, 0, 0, 8, 11, 0, 0, 8, 2, -4, 0, -8, 0]); E
                [  0  14   0  -8  -2]
                [-14   0  -3 -11   4]
                [  0   3   0   0   0]
                [  8  11   0   0   8]
                [  2  -4   0  -8   0]
                sage: F, C = E.symplectic_form()
                sage: F
                [ 0  0  1  0  0]
                [ 0  0  0  2  0]
                [-1  0  0  0  0]
                [ 0 -2  0  0  0]
                [ 0  0  0  0  0]
                sage: F == C * E * C.transpose()
                True
                sage: E.smith_form()[0]
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 2 0 0]
                [0 0 0 2 0]
                [0 0 0 0 0]
            """
        import sage.matrix.symplectic_basis
        return sage.matrix.symplectic_basis.symplectic_basis_over_ZZ(self)

    hermite_form = echelon_form

    def echelon_form(self, algorithm="default", proof=None, include_zero_rows=True,
                     transformation=False, D=None):
        r"""
        Return the echelon form of this matrix over the integers, also known
        as the hermit normal form (HNF).

        INPUT:

        - ``algorithm`` -- String. The algorithm to use. Valid options are:

          - ``'default'`` -- Let Sage pick an algorithm (default). Up
            to 10 rows or columns: pari with flag 0; Up to 75 rows or
            columns: pari with flag 1; Larger: use padic algorithm.

          - ``'padic'`` - an asymptotically fast p-adic modular
            algorithm, If your matrix has large coefficients and is
            small, you may also want to try this.

          - ``'pari'`` - use PARI with flag 1

          - ``'pari0'`` - use PARI with flag 0

          - ``'pari4'`` - use PARI with flag 4 (use heuristic LLL)

          - ``'ntl'`` - use NTL (only works for square matrices of
            full rank!)

        -  ``proof`` - (default: True); if proof=False certain
           determinants are computed using a randomized hybrid p-adic
           multimodular strategy until it stabilizes twice (instead of up to
           the Hadamard bound). It is *incredibly* unlikely that one would
           ever get an incorrect result with proof=False.

        -  ``include_zero_rows`` - (default: True) if False,
           don't include zero rows

        -  ``transformation`` - if given, also compute
           transformation matrix; only valid for padic algorithm

        -  ``D`` - (default: None) if given and the algorithm
           is 'ntl', then D must be a multiple of the determinant and this
           function will use that fact.

        OUTPUT:

        The Hermite normal form (=echelon form over `\ZZ`) of self.

        EXAMPLES::

            sage: A = MatrixSpace(ZZ,2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]
            sage: A = MatrixSpace(ZZ,5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]

        Getting a transformation matrix in the nonsquare case::

            sage: A = matrix(ZZ,5,3,[1..15])
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: H
            [1 2 3]
            [0 3 6]
            sage: U
            [  0   0   0   4  -3]
            [  0   0   0  13 -10]
            sage: U*A == H
            True

        TESTS: Make sure the zero matrices are handled correctly::

            sage: m = matrix(ZZ,3,3,[0]*9)
            sage: m.echelon_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: m = matrix(ZZ,3,1,[0]*3)
            sage: m.echelon_form()
            [0]
            [0]
            [0]
            sage: m = matrix(ZZ,1,3,[0]*3)
            sage: m.echelon_form()
            [0 0 0]

        The ultimate border case!

        ::

            sage: m = matrix(ZZ,0,0,[])
            sage: m.echelon_form()
            []

        .. note::

           If 'ntl' is chosen for a non square matrix this function
           raises a ValueError.

        Special cases: 0 or 1 rows::

            sage: a = matrix(ZZ, 1,2,[0,-1])
            sage: a.hermite_form()
            [0 1]
            sage: a.pivots()
            (1,)
            sage: a = matrix(ZZ, 1,2,[0,0])
            sage: a.hermite_form()
            [0 0]
            sage: a.pivots()
            ()
            sage: a = matrix(ZZ,1,3); a
            [0 0 0]
            sage: a.echelon_form(include_zero_rows=False)
            []
            sage: a.echelon_form(include_zero_rows=True)
            [0 0 0]

        Illustrate using various algorithms.::

            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari0')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari4')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='padic')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='default')
            [1 2 3]
            [0 3 6]
            [0 0 0]

        The 'ntl' algorithm doesn't work on matrices that do not have full rank.::

            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='ntl')
            Traceback (most recent call last):
            ...
            ValueError: ntl only computes HNF for square matrices of full rank.
            sage: matrix(ZZ,3,[0] +[2..9]).hermite_form(algorithm='ntl')
            [1 0 0]
            [0 1 0]
            [0 0 3]

        TESTS:

        This example illustrated trac 2398::

            sage: a = matrix([(0, 0, 3), (0, -2, 2), (0, 1, 2), (0, -2, 5)])
            sage: a.hermite_form()
            [0 1 2]
            [0 0 3]
            [0 0 0]
            [0 0 0]

        Check that #12280 is fixed::

            sage: m = matrix([(-2, 1, 9, 2, -8, 1, -3, -1, -4, -1),
            ...               (5, -2, 0, 1, 0, 4, -1, 1, -2, 0),
            ...               (-11, 3, 1, 0, -3, -2, -1, -11, 2, -2),
            ...               (-1, 1, -1, -2, 1, -1, -1, -1, -1, 7),
            ...               (-2, -1, -1, 1, 1, -2, 1, 0, 2, -4)]).stack(
            ...               200 * identity_matrix(ZZ, 10))
            sage: matrix(ZZ,m).hermite_form(algorithm='pari', include_zero_rows=False)
            [  1   0   2   0  13   5   1 166  72  69]
            [  0   1   1   0  20   4  15 195  65 190]
            [  0   0   4   0  24   5  23  22  51 123]
            [  0   0   0   1  23   7  20 105  60 151]
            [  0   0   0   0  40   4   0  80  36  68]
            [  0   0   0   0   0  10   0 100 190 170]
            [  0   0   0   0   0   0  25   0 100 150]
            [  0   0   0   0   0   0   0 200   0   0]
            [  0   0   0   0   0   0   0   0 200   0]
            [  0   0   0   0   0   0   0   0   0 200]
            sage: matrix(ZZ,m).hermite_form(algorithm='padic', include_zero_rows=False)
            [  1   0   2   0  13   5   1 166  72  69]
            [  0   1   1   0  20   4  15 195  65 190]
            [  0   0   4   0  24   5  23  22  51 123]
            [  0   0   0   1  23   7  20 105  60 151]
            [  0   0   0   0  40   4   0  80  36  68]
            [  0   0   0   0   0  10   0 100 190 170]
            [  0   0   0   0   0   0  25   0 100 150]
            [  0   0   0   0   0   0   0 200   0   0]
            [  0   0   0   0   0   0   0   0 200   0]
            [  0   0   0   0   0   0   0   0   0 200]
        """
        if self._nrows == 0 or self._ncols == 0:
            self.cache('pivots', ())
            self.cache('rank', 0)
            if transformation:
                return self, self
            return self

        key = 'hnf-%s-%s'%(include_zero_rows,transformation)
        ans = self.fetch(key)
        if ans is not None: return ans

        cdef Py_ssize_t nr, nc, n, i, j
        nr = self._nrows
        nc = self._ncols
        n = nr if nr >= nc else nc
        if algorithm is 'default':
            if transformation: algorithm = 'padic'
            else:
                if n <= 10: algorithm = 'pari0'
                elif n <= 75: algorithm = 'pari'
                else: algorithm = 'padic'

        cdef bint pari_big = 0
        if algorithm.startswith('pari'):
            if self.height().ndigits() > 10000 or n >= 50:
                pari_big = 1

        cdef Matrix_integer_dense H_m

        proof = get_proof_flag(proof, "linear_algebra")
        pivots = None
        rank = None

        if algorithm == "padic":
            import matrix_integer_dense_hnf
            if transformation:
                H_m, U, pivots = matrix_integer_dense_hnf.hnf_with_transformation(self, proof=proof)
                if not include_zero_rows:
                    r = H_m.rank()
                    H_m = H_m[:r]
                    U = U[:r]
            else:
                H_m, pivots = matrix_integer_dense_hnf.hnf(self,
                                   include_zero_rows=include_zero_rows, proof=proof)
            self.cache('pivots', tuple(pivots))
            self.cache('rank', len(pivots))


        elif algorithm == 'pari0':
            if transformation:
                raise ValueError("transformation matrix only available with p-adic algorithm")
            if pari_big:
                H_m = self._hnf_pari_big(0, include_zero_rows=include_zero_rows)
            else:
                H_m = self._hnf_pari(0, include_zero_rows=include_zero_rows)

        elif algorithm == 'pari':
            if transformation:
                raise ValueError("transformation matrix only available with p-adic algorithm")
            if pari_big:
                H_m = self._hnf_pari_big(1, include_zero_rows=include_zero_rows)
            else:
                H_m = self._hnf_pari(1, include_zero_rows=include_zero_rows)

        elif algorithm == 'pari4':
            if transformation:
                raise ValueError("transformation matrix only available with p-adic algorithm")
            if pari_big:
                H_m = self._hnf_pari_big(4, include_zero_rows=include_zero_rows)
            else:
                H_m = self._hnf_pari(4, include_zero_rows=include_zero_rows)

        elif algorithm == 'ntl':
            if transformation:
                raise ValueError("transformation matrix only available with p-adic algorithm")

            if nr != nc:
                raise ValueError("ntl only computes HNF for square matrices of full rank.")

            import sage.libs.ntl.ntl_mat_ZZ
            v =  sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols)
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    v[i,j] = self.get_unsafe(nr-i-1,nc-j-1)

            try:
                w = v.HNF(D=D)
            except RuntimeError: # HNF may fail if a nxm matrix has rank < m
                raise ValueError("ntl only computes HNF for square matrices of full rank.")

            rank = w.nrows()

            if include_zero_rows:
                H_m = self.new_matrix()
            else:
                H_m = self.new_matrix(nrows=w.nrows())

            nr = w.nrows()
            nc = w.ncols()

            for i from 0 <= i < w.nrows():
                for j from 0 <= j < w.ncols():
                    H_m[i,j] = w[nr-i-1,nc-j-1]

        else:
            raise TypeError("algorithm '%s' not understood"%(algorithm))

        H_m.set_immutable()
        if pivots is None:
            from matrix_integer_dense_hnf import pivots_of_hnf_matrix
            pivots = tuple(pivots_of_hnf_matrix(H_m))
            rank = len(pivots)
        else:
            pivots = tuple(pivots)

        H_m.cache('pivots', pivots)
        self.cache('pivots', pivots)

        H_m.cache('rank', rank)
        self.cache('rank',rank)

        H_m.cache('in_echelon_form', True)


        if transformation:
            U.set_immutable()
            ans = H_m, U
        else:
            ans = H_m
        self.cache(key, ans)
        return ans

    def saturation(self, p=0, proof=None, max_dets=5):
        r"""
        Return a saturation matrix of self, which is a matrix whose rows
        span the saturation of the row span of self. This is not unique.

        The saturation of a `\ZZ` module `M`
        embedded in `\ZZ^n` is the a module `S` that
        contains `M` with finite index such that
        `\ZZ^n/S` is torsion free. This function takes the
        row span `M` of self, and finds another matrix of full rank
        with row span the saturation of `M`.

        INPUT:


        -  ``p`` - (default: 0); if nonzero given, saturate
           only at the prime `p`, i.e., return a matrix whose row span
           is a `\ZZ`-module `S` that contains self and
           such that the index of `S` in its saturation is coprime to
           `p`. If `p` is None, return full saturation of
           self.

        -  ``proof`` - (default: use proof.linear_algebra());
           if False, the determinant calculations are done with proof=False.

        -  ``max_dets`` - (default: 5); technical parameter -
           max number of determinant to compute when bounding prime divisor of
           self in its saturation.


        OUTPUT:


        -  ``matrix`` - a matrix over ZZ


        .. note::

           The result is *not* cached.

        ALGORITHM: 1. Replace input by a matrix of full rank got from a
        subset of the rows. 2. Divide out any common factors from rows. 3.
        Check max_dets random dets of submatrices to see if their GCD
        (with p) is 1 - if so matrix is saturated and we're done. 4.
        Finally, use that if A is a matrix of full rank, then
        `hnf(transpose(A))^{-1}*A` is a saturation of A.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, 5, [-51, -1509, -71, -109, -593, -19, -341, 4, 86, 98, 0, -246, -11, 65, 217])
            sage: A.echelon_form()
            [      1       5    2262   20364   56576]
            [      0       6   35653  320873  891313]
            [      0       0   42993  386937 1074825]
            sage: S = A.saturation(); S
            [  -51 -1509   -71  -109  -593]
            [  -19  -341     4    86    98]
            [   35   994    43    51   347]

        Notice that the saturation spans a different module than A.

        ::

            sage: S.echelon_form()
            [ 1  2  0  8 32]
            [ 0  3  0 -2 -6]
            [ 0  0  1  9 25]
            sage: V = A.row_space(); W = S.row_space()
            sage: V.is_submodule(W)
            True
            sage: V.index_in(W)
            85986
            sage: V.index_in_saturation()
            85986

        We illustrate each option::

            sage: S = A.saturation(p=2)
            sage: S = A.saturation(proof=False)
            sage: S = A.saturation(max_dets=2)
        """
        proof = get_proof_flag(proof, "linear_algebra")
        import matrix_integer_dense_saturation
        return matrix_integer_dense_saturation.saturation(self, p=p, proof=proof, max_dets=max_dets)

    def index_in_saturation(self, proof=None):
        """
        Return the index of self in its saturation.

        INPUT:


        -  ``proof`` - (default: use proof.linear_algebra());
           if False, the determinant calculations are done with proof=False.


        OUTPUT:


        -  ``positive integer`` - the index of the row span of
           this matrix in its saturation


        ALGORITHM: Use Hermite normal form twice to find an invertible
        matrix whose inverse transforms a matrix with the same row span as
        self to its saturation, then compute the determinant of that
        matrix.

        EXAMPLES::

            sage: A = matrix(ZZ, 2,3, [1..6]); A
            [1 2 3]
            [4 5 6]
            sage: A.index_in_saturation()
            3
            sage: A.saturation()
            [1 2 3]
            [1 1 1]
        """
        proof = get_proof_flag(proof, "linear_algebra")
        import matrix_integer_dense_saturation
        return matrix_integer_dense_saturation.index_in_saturation(self, proof=proof)

    def pivots(self):
        """
        Return the pivot column positions of this matrix.

        OUTPUT: a tuple of Python integers: the position of the
        first nonzero entry in each row of the echelon form.

        EXAMPLES::

            sage: n = 3; A = matrix(ZZ,n,range(n^2)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.pivots()
            (0, 1)
            sage: A.echelon_form()
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]
        """
        p = self.fetch('pivots')
        if not p is None: return tuple(p)

        cdef Matrix_integer_dense E
        E = self.echelon_form()

        # Now we determine the pivots from the matrix E as quickly as we can.
        # For each row, we find the first nonzero position in that row -- it is the pivot.
        cdef Py_ssize_t i, j, k
        cdef mpz_t *row
        p = []
        k = 0
        for i from 0 <= i < E._nrows:
            row = E._matrix[i]   # pointer to ith row
            for j from k <= j < E._ncols:
                if mpz_cmp_si(row[j], 0) != 0:  # nonzero position
                    p.append(j)
                    k = j+1  # so start at next position next time
                    break
        p = tuple(p)
        self.cache('pivots', p)
        return p

    #### Elementary divisors

    def elementary_divisors(self, algorithm='pari'):
        """
        Return the elementary divisors of self, in order.


        .. warning::

           This is MUCH faster than the smith_form function.

        The elementary divisors are the invariants of the finite abelian
        group that is the cokernel of *left* multiplication of this matrix.
        They are ordered in reverse by divisibility.

        INPUT:


        -  ``self`` - matrix

        -  ``algorithm`` - (default: 'pari')

           - ``'pari'``: works robustly, but is slower.

           - ``'linbox'`` - use linbox (currently off, broken)


        OUTPUT: list of integers


        .. note::

           These are the invariants of the cokernel of *left* multiplication::

               sage: M = Matrix([[3,0,1],[0,1,0]])
               sage: M
               [3 0 1]
               [0 1 0]
               sage: M.elementary_divisors()
               [1, 1]
               sage: M.transpose().elementary_divisors()
               [1, 1, 0]

        EXAMPLES::

            sage: matrix(3, range(9)).elementary_divisors()
            [1, 3, 0]
            sage: matrix(3, range(9)).elementary_divisors(algorithm='pari')
            [1, 3, 0]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [1, 1, 1, 687]

        ::

            sage: M = matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2])
            sage: M.elementary_divisors()
            [1, 1, 6]

        This returns a copy, which is safe to change::

            sage: edivs = M.elementary_divisors()
            sage: edivs.pop()
            6
            sage: M.elementary_divisors()
            [1, 1, 6]

        .. seealso::

           :meth:`smith_form`
        """
        d = self.fetch('elementary_divisors')
        if not d is None:
            return d[:]
        if self._nrows == 0 or self._ncols == 0:
            d = []
        else:
            if algorithm == 'linbox':
                raise ValueError("linbox too broken -- currently Linbox SNF is disabled.")
                # This fails in linbox: a = matrix(ZZ,2,[1, 1, -1, 0, 0, 0, 0, 1])
                try:
                    d = self._elementary_divisors_linbox()
                except RuntimeError:
                    import sys
                    sys.stderr.write("DONT PANIC -- switching to using PARI (which will work fine)\n")
                    algorithm = 'pari'
            if algorithm == 'pari':
                d = self._pari_().matsnf(0).python()
                i = d.count(0)
                d.sort()
                if i > 0:
                    d = d[i:] + [d[0]]*i
            elif not (algorithm in ['pari', 'linbox']):
                raise ValueError("algorithm (='%s') unknown"%algorithm)
        self.cache('elementary_divisors', d)
        return d[:]

    def _elementary_divisors_linbox(self):
        self._init_linbox()
        sig_on()
        d = linbox.smithform()
        sig_off()
        return d

    def smith_form(self):
        r"""
        Returns matrices S, U, and V such that S = U\*self\*V, and S is in
        Smith normal form. Thus S is diagonal with diagonal entries the
        ordered elementary divisors of S.

        .. warning::

           The elementary_divisors function, which returns the
           diagonal entries of S, is VASTLY faster than this function.

        The elementary divisors are the invariants of the finite abelian
        group that is the cokernel of this matrix. They are ordered in
        reverse by divisibility.

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0 0]
            [0 3 0]
            [0 0 0]
            sage: U
            [ 0  1  0]
            [ 0 -1  1]
            [-1  2 -1]
            sage: V
            [-1  4  1]
            [ 1 -3 -2]
            [ 0  0  1]
            sage: U*A*V
            [1 0 0]
            [0 3 0]
            [0 0 0]

        It also makes sense for nonsquare matrices::

            sage: A = Matrix(ZZ,3,2,range(6))
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0]
            [0 2]
            [0 0]
            sage: U
            [ 0  1  0]
            [ 0 -1  1]
            [-1  2 -1]
            sage: V
            [-1  3]
            [ 1 -2]
            sage: U * A * V
            [1 0]
            [0 2]
            [0 0]

        Empty matrices are handled sensibly (see trac #3068)::

            sage: m = MatrixSpace(ZZ, 2,0)(0); d,u,v = m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(ZZ, 0,2)(0); d,u,v = m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(ZZ, 0,0)(0); d,u,v = m.smith_form(); u*m*v == d
            True

        .. seealso::

           :meth:`elementary_divisors`
        """
        v = self._pari_().matsnf(1).python()
        if self._ncols == 0: v[0] = self.matrix_space(ncols = self._nrows)(1)
        if self._nrows == 0: v[1] = self.matrix_space(nrows = self._ncols)(1)
        # need to reverse order of rows of U, columns of V, and both of D.
        D = self.matrix_space()([v[2][i,j] for i in xrange(self._nrows-1,-1,-1) for j in xrange(self._ncols-1,-1,-1)])

        if self._ncols == 0:
            # silly special cases for matrices with 0 columns (PARI has a unique empty matrix)
            U = self.matrix_space(ncols = self._nrows)(1)
        else:
            U = self.matrix_space(ncols = self._nrows)([v[0][i,j] for i in xrange(self._nrows-1,-1,-1) for j in xrange(self._nrows)])

        if self._nrows == 0:
            # silly special cases for matrices with 0 rows (PARI has a unique empty matrix)
            V = self.matrix_space(nrows = self._ncols)(1)
        else:
            V = self.matrix_space(nrows = self._ncols)([v[1][i,j] for i in xrange(self._ncols) for j in xrange(self._ncols-1,-1,-1)])

        return D, U, V

    def frobenius(self, flag=0, var='x'):
        """
        Return the Frobenius form (rational canonical form) of this
        matrix.

        INPUT:

        -  ``flag`` -- 0 (default), 1 or 2 as follows:

            -  ``0`` -- (default) return the Frobenius form of this
               matrix.

            -  ``1`` -- return only the elementary divisor
               polynomials, as polynomials in var.

            -  ``2`` -- return a two-components vector [F,B] where F
               is the Frobenius form and B is the basis change so that
               `M=B^{-1}FB`.

        -  ``var`` -- a string (default: 'x')

        ALGORITHM: uses PARI's matfrobenius()

        EXAMPLES::

            sage: A = MatrixSpace(ZZ, 3)(range(9))
            sage: A.frobenius(0)
            [ 0  0  0]
            [ 1  0 18]
            [ 0  1 12]
            sage: A.frobenius(1)
            [x^3 - 12*x^2 - 18*x]
            sage: A.frobenius(1, var='y')
            [y^3 - 12*y^2 - 18*y]
            sage: F, B = A.frobenius(2)
            sage: A == B^(-1)*F*B
            True
            sage: a=matrix([])
            sage: a.frobenius(2)
            ([], [])
            sage: a.frobenius(0)
            []
            sage: a.frobenius(1)
            []
            sage: B = random_matrix(ZZ,2,3)
            sage: B.frobenius()
            Traceback (most recent call last):
            ...
            ArithmeticError: frobenius matrix of non-square matrix not defined.

        AUTHORS:

        - Martin Albrect (2006-04-02)

        TODO: - move this to work for more general matrices than just over
        Z. This will require fixing how PARI polynomials are coerced to
        Sage polynomials.
        """
        if not self.is_square():
            raise ArithmeticError("frobenius matrix of non-square matrix not defined.")

        v = self._pari_().matfrobenius(flag)
        if flag==0:
            return self.matrix_space()(v.python())
        elif flag==1:
            r = PolynomialRing(self.base_ring(), names=var)
            retr = []
            for f in v:
                retr.append(eval(str(f).replace("^","**"), {'x':r.gen()}, r.gens_dict()))
            return retr
        elif flag==2:
            F = matrix_space.MatrixSpace(QQ, self.nrows())(v[0].python())
            B = matrix_space.MatrixSpace(QQ, self.nrows())(v[1].python())
            return F, B

    def _right_kernel_matrix(self, **kwds):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` - determines which algorithm to use, options are:

          - 'pari' - use the ``matkerint()`` function from the PARI library
          - 'padic' - use the p-adic algorithm from the IML library
          - 'default' - use a heuristic to decide which of the two above
            routines is fastest.  This is the default value.

        - ``proof`` - this is passed to the p-adic IML algorithm.
          If not specified, the global flag for linear algebra will be used.

        OUTPUT:

        Returns a pair.  First item is the string is either
        'computed-pari-int' or 'computed-iml-int', which identifies
        the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integers, as computed by either the IML or PARI libraries.

        EXAMPLES::

            sage: A = matrix(ZZ, [[4, 7, 9, 7, 5, 0],
            ...                   [1, 0, 5, 8, 9, 1],
            ...                   [0, 1, 0, 1, 9, 7],
            ...                   [4, 7, 6, 5, 1, 4]])

            sage: result = A._right_kernel_matrix(algorithm='pari')
            sage: result[0]
            'computed-pari-int'
            sage: X = result[1]; X
            [-26  31 -30  21   2 -10]
            [-47 -13  48 -14 -11  18]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 2)
            True

            sage: result = A._right_kernel_matrix(algorithm='padic')
            sage: result[0]
            'computed-iml-int'
            sage: X = result[1]; X
            [-469  214  -30  119  -37    0]
            [ 370 -165   18  -91   30   -2]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 2)
            True

            sage: result = A._right_kernel_matrix(algorithm='default')
            sage: result[0]
            'computed-pari-int'
            sage: result[1]
            [-26  31 -30  21   2 -10]
            [-47 -13  48 -14 -11  18]

            sage: result = A._right_kernel_matrix()
            sage: result[0]
            'computed-pari-int'
            sage: result[1]
            [-26  31 -30  21   2 -10]
            [-47 -13  48 -14 -11  18]

        With the 'default' given as the algorithm, several heuristics are
        used to determine if PARI or IML ('padic') is used.  The code has
        exact details, but roughly speaking, relevant factors are: the
        absolute size of the matrix, or the relative dimensions, or the
        magnitude of the entries. ::

            sage: A = random_matrix(ZZ, 18, 11)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'
            sage: A = random_matrix(ZZ, 18, 11, x = 10^200)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 60)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 55)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'

        TESTS:

        We test three trivial cases. PARI is used for small matrices,
        but we let the heuristic decide that.  ::

            sage: A = matrix(ZZ, 0, 2)
            sage: A._right_kernel_matrix()[1]
            [1 0]
            [0 1]
            sage: A = matrix(ZZ, 2, 0)
            sage: A._right_kernel_matrix()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
            sage: A = zero_matrix(ZZ, 4, 3)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over the integers for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)

        algorithm = kwds.pop('algorithm', None)
        if algorithm == None:
          algorithm = 'default'

        if algorithm == 'default':
            # The heuristic here could be auto-tuned, stored for
            # different architecture, etc.  What I've done below here
            # I just got by playing around with examples.  This is
            # *dramatically* better than doing absolutely nothing
            # (i.e., always choosing 'padic'), but is of course
            # far from optimal.   -- William Stein
            if max(self._nrows, self._ncols) <= 10:
                # pari much better for very small matrices, as long as entries aren't huge.
                algorithm = 'pari'
            if max(self._nrows, self._ncols) <= 50:
                # when entries are huge, padic relatively good.
                h = self.height().ndigits()
                if h < 100:
                    algorithm = 'pari'
                else:
                    algorithm = 'padic'
            elif self._nrows <= self._ncols + 3:
                # the padic algorithm is much better for bigger
                # matrices if there are nearly more columns than rows
                # (that is its forte)
                algorithm = 'padic'
            else:
                algorithm = 'pari'

        if algorithm == 'pari':
            K = self._pari_().matkerint().mattranspose().python()
            format = 'computed-pari-int'
        if algorithm == 'padic':
            proof = kwds.pop('proof', None)
            proof = get_proof_flag(proof, "linear_algebra")
            K = self._rational_kernel_iml().transpose().saturation(proof=proof)
            format = 'computed-iml-int'
        tm = verbose("done computing right kernel matrix over the integers for %sx%s matrix" % (self.nrows(), self.ncols()),level=1, t=tm)
        return (format, K)

    def _adjoint(self):
        """
        Return the adjoint of this matrix.

        Assumes self is a square matrix (checked in adjoint).

        EXAMPLES::

            sage: m = matrix(ZZ,3,[1..9])
            sage: m.adjoint()
            [ -3   6  -3]
            [  6 -12   6]
            [ -3   6  -3]
        """
        return self.parent()(self._pari_().matadjoint().python())

    def _ntl_(self):
        r"""
        ntl.mat_ZZ representation of self.

        EXAMPLE::

            sage: a = MatrixSpace(ZZ,200).random_element(x=-2, y=2)    # -2 to 2
            sage: A = a._ntl_()

        .. note::

           NTL only knows dense matrices, so if you provide a sparse
           matrix NTL will allocate memory for every zero entry.
        """
        import sage.libs.ntl.ntl_mat_ZZ
        return sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols, self.list())


    ####################################################################################
    # LLL
    ####################################################################################
    def LLL_gram(self):
        """
        LLL reduction of the lattice whose gram matrix is self.

        INPUT:


        -  ``M`` - gram matrix of a definite quadratic form


        OUTPUT:


        -  ``U`` - unimodular transformation matrix such that
           U.transpose() \* M \* U  is LLL-reduced.

        ALGORITHM: Use PARI

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 2, [5,3,3,2]) ; M
            [5 3]
            [3 2]
            sage: U = M.LLL_gram(); U
            [-1  1]
            [ 1 -2]
            sage: U.transpose() * M * U
            [1 0]
            [0 1]

        Semidefinite and indefinite forms no longer raise a ValueError::

            sage: Matrix(ZZ,2,2,[2,6,6,3]).LLL_gram()
            [-3 -1]
            [ 1  0]
            sage: Matrix(ZZ,2,2,[1,0,0,-1]).LLL_gram()
            [ 0 -1]
            [ 1  0]

        """
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        n = self.nrows()
        # maybe should be /unimodular/ matrices ?
        P = self._pari_()
        try:
            U = P.lllgramint()
        except (RuntimeError, ArithmeticError), msg:
            raise ValueError("not a definite matrix")
        MS = matrix_space.MatrixSpace(ZZ,n)
        U = MS(U.python())
        # Fix last column so that det = +1
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U

    def BKZ(self, delta=None, fp="rr", block_size=10, prune=0, use_givens=False):
        """
        Block Korkin-Zolotarev reduction.

        INPUT:


        -  ``fp`` - 'fp' - double precision: NTL's FP or
           fpLLL's double

        -  ``'qd'`` - quad doubles: NTL's QP

        -  ``'qd1'`` - quad doubles: uses quad_float precision
           to compute Gram-Schmidt, but uses double precision in the search
           phase of the block reduction algorithm. This seems adequate for
           most purposes, and is faster than 'qd', which uses quad_float
           precision uniformly throughout.

        -  ``'xd'`` - extended exponent: NTL's XD

        -  ``'rr'`` - arbitrary precision (default)

        -  ``block_size`` - specifies the size of the blocks
           in the reduction. High values yield shorter vectors, but the
           running time increases exponentially with
           ``block_size``. ``block_size`` should be
           between 2 and the number of rows of ``self`` (default:
           10)

        -  ``prune`` - The optional parameter
           ``prune`` can be set to any positive number to invoke
           the Volume Heuristic from [Schnorr and Horner, Eurocrypt '95]. This
           can significantly reduce the running time, and hence allow much
           bigger block size, but the quality of the reduction is of course
           not as good in general. Higher values of ``prune`` mean
           better quality, and slower running time. When ``prune``
           == 0, pruning is disabled. Recommended usage: for
           ``block_size`` = 30, set 10 = ``prune`` =
           15.

        -  ``use_givens`` - use Given's orthogonalization.
           This is a bit slower, but generally much more stable, and is really
           the preferred orthogonalization strategy. For a nice description of
           this, see Chapter 5 of [G. Golub and C. van Loan, Matrix
           Computations, 3rd edition, Johns Hopkins Univ. Press, 1996].


        EXAMPLE::

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ()
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]
            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ(use_givens=True)
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

        ::

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ(fp="fp")
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]
        """
        if delta is None:
            delta = 0.99
        elif delta <= 0.25:
            raise TypeError("delta must be > 0.25")
        elif delta > 1:
            raise TypeError("delta must be <= 1")
        delta = float(delta)

        if fp is None:
            fp = "rr"

        if fp == "fp":
            algorithm = "BKZ_FP"
        elif fp == "qd":
            algorithm = "BKZ_QP"
        elif fp == "qd1":
            algorithm = "BKZ_QP1"
        elif fp == "xd":
            algorithm = "BKZ_XD"
        elif fp == "rr":
            algorithm = "BKZ_RR"
        else:
            raise TypeError("fp parameter not understood.")

        block_size = int(block_size)

        if prune < 0:
            raise TypeError("prune must be >= 0")
        prune = int(prune)

        if get_verbose() >= 2:
            verbose = True
        else:
            verbose = False

        A = self._ntl_()

        if algorithm == "BKZ_FP":
            if not use_givens:
                r = A.BKZ_FP(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)
            else:
                r = A.G_BKZ_FP(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)

        elif algorithm == "BKZ_QP":
            if not use_givens:
                r = A.BKZ_QP(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)
            else:
                r = A.G_BKZ_QP(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)

        elif algorithm == "BKZ_QP1":
            if not use_givens:
                r = A.BKZ_QP1(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)
            else:
                r = A.G_BKZ_QP1(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)

        elif algorithm == "BKZ_XD":
            if not use_givens:
                r = A.BKZ_XD(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)
            else:
                r = A.G_BKZ_XD(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)

        elif algorithm == "BKZ_RR":
            if not use_givens:
                r = A.BKZ_RR(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)
            else:
                r = A.G_BKZ_RR(U=None, delta=delta, BlockSize=block_size, prune=prune, verbose=verbose)

        self.cache("rank",ZZ(r))
        R = <Matrix_integer_dense>self.new_matrix(entries=map(ZZ,A.list()))
        return R

    def LLL(self, delta=None, eta=None, algorithm="fpLLL:wrapper", fp=None, prec=0, early_red = False, use_givens = False):
        r"""
        Returns LLL reduced or approximated LLL reduced lattice R for this
        matrix interpreted as a lattice.

        A lattice `(b_1, b_2, ..., b_d)` is
        `(\delta, \eta)` -LLL-reduced if the two following
        conditions hold:

        -  For any `i>j`, we have `|mu_{i, j}| <= \eta`,
        -  For any `i<d`, we have
           `\delta |b_i^*|^2 <= |b_{i + 1}^* + mu_{i + 1, i} b_i^* |^2`,

        where `mu_{i,j} = <b_i, b_j^*>/<b_j^*,b_j^*>` and
        `b_i^*` is the `i`-th vector of the Gram-Schmidt
        orthogonalisation of `(b_1, b_2, ..., b_d)`.

        The default reduction parameters are `\delta=3/4` and
        `eta=0.501`. The parameters `\delta` and
        `\eta` must satisfy: `0.25 < \delta <= 1.0` and
        `0.5 <= \eta < sqrt(\delta)`. Polynomial time
        complexity is only guaranteed for `\delta < 1`.

        The lattice is returned as a matrix. Also the rank (and the
        determinant) of self are cached if those are computed during the
        reduction. Note that in general this only happens when self.rank()
        == self.ncols() and the exact algorithm is used.

        INPUT:


        -  ``delta`` - parameter as described above (default:
           3/4)

        -  ``eta`` - parameter as described above (default:
           0.501), ignored by NTL

        -  ``algorithm`` - string (default: "fpLLL:wrapper")
           one of the algorithms mentioned below

        -  ``fp``

            -  None - NTL's exact reduction or fpLLL's
               wrapper

            -  ``'fp'`` - double precision: NTL's FP or fpLLL's
               double

            -  ``'qd'`` - quad doubles: NTL's QP

            -  ``'xd'`` - extended exponent: NTL's XD or fpLLL's
               dpe

            -  ``'rr'`` - arbitrary precision: NTL'RR or fpLLL's
               MPFR

        -  ``prec`` - precision, ignored by NTL (default: auto
           choose)

        -  ``early_red`` - perform early reduction, ignored by
           NTL (default: False)

        -  ``use_givens`` - use Givens orthogonalization
           (default: False) only applicable to approximate reductions and NTL.
           This is more stable but slower.

           Also, if the verbose level is = 2, some more verbose output is
           printed during the calculation if NTL is used.

           AVAILABLE ALGORITHMS:

        -  ``NTL:LLL`` - NTL's LLL + fp

        -  ``fpLLL:heuristic`` - fpLLL's heuristic + fp

        -  ``fpLLL:fast`` - fpLLL's fast

        -  ``fpLLL:wrapper`` - fpLLL's automatic choice
           (default)


        OUTPUT: a matrix over the integers

        EXAMPLE::

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.LLL()
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

        We compute the extended GCD of a list of integers using LLL, this
        example is from the Magma handbook::

            sage: Q = [ 67015143, 248934363018, 109210, 25590011055, 74631449, \
                        10230248, 709487, 68965012139, 972065, 864972271 ]
            sage: n = len(Q)
            sage: S = 100
            sage: X = Matrix(ZZ, n, n + 1)
            sage: for i in xrange(n):
            ...       X[i,i + 1] = 1
            sage: for i in xrange(n):
            ...       X[i,0] = S*Q[i]
            sage: L = X.LLL()
            sage: M = L.row(n-1).list()[1:]
            sage: M
            [-3, -1, 13, -1, -4, 2, 3, 4, 5, -1]
            sage: add([Q[i]*M[i] for i in range(n)])
            -1

        TESTS::

            sage: matrix(ZZ, 0, 0).LLL()
            []
            sage: matrix(ZZ, 3, 0).LLL()
            []
            sage: matrix(ZZ, 0, 3).LLL()
            []

            sage: M = matrix(ZZ, [[1,2,3],[31,41,51],[101,201,301]])
            sage: A = M.LLL()
            sage: A
            [ 0  0  0]
            [-1  0  1]
            [ 1  1  1]
            sage: B = M.LLL(algorithm='NTL:LLL')
            sage: C = M.LLL(algorithm='NTL:LLL', fp=None)
            sage: D = M.LLL(algorithm='NTL:LLL', fp='fp')
            sage: F = M.LLL(algorithm='NTL:LLL', fp='xd')
            sage: G = M.LLL(algorithm='NTL:LLL', fp='rr')
            sage: A == B == C == D == F == G
            True
            sage: H = M.LLL(algorithm='NTL:LLL', fp='qd')
            Traceback (most recent call last):
            ...
            TypeError: algorithm NTL:LLL_QD not supported

        ALGORITHM: Uses the NTL library by Victor Shoup or fpLLL library by
        Damien Stehle depending on the chosen algorithm.

        REFERENCES:

        - ``ntl.mat_ZZ`` or ``sage.libs.fplll.fplll`` for details on
          the used algorithms.
        """
        if self.ncols()==0 or self.nrows()==0:
            verbose("Trivial matrix, nothing to do")
            return self

        tm = verbose("LLL of %sx%s matrix (algorithm %s)"%(self.nrows(), self.ncols(), algorithm))
        import sage.libs.ntl.all
        ntl_ZZ = sage.libs.ntl.all.ZZ

        from sage.libs.fplll.fplll import FP_LLL

        if get_verbose() >= 2: verb = True
        else: verb = False

        # auto choice

        # FP choice
        if algorithm == 'NTL:LLL':
            if fp == None:
                algorithm = 'NTL:LLL_FP'
            elif fp == 'fp':
                algorithm = 'NTL:LLL_FP'
            elif fp == 'qd':
                algorithm = 'NTL:LLL_QD'
            elif fp == 'xd':
                algorithm = 'NTL:LLL_XD'
            elif fp == 'rr':
                algorithm = 'NTL:LLL_RR'
        elif algorithm == 'fpLLL:heuristic':
            if fp == None:
                raise TypeError("if 'fpLLL:heuristic' is chosen, a floating point number implementation must be chosen")
            elif fp == 'fp':
                fp = 'double'
            elif fp == 'qd':
                raise TypeError("fpLLL does not support quad doubles.")
            elif fp == 'xd':
                fp = 'dpe'
            elif fp == 'rr':
                fp = 'mpfr'

        if algorithm == "NTL:LLL":
            if delta is None:
                delta = ZZ(3)/ZZ(4)
            elif delta <= ZZ(1)/ZZ(4):
                raise TypeError("delta must be > 1/4")
            elif delta > 1:
                raise TypeError("delta must be <= 1")
            delta = QQ(delta)
            a = delta.numer()
            b = delta.denom()

        else:
            if delta is None:
                delta = 0.99
            elif delta <= 0.25:
                raise TypeError("delta must be > 0.25")
            elif delta > 1:
                raise TypeError("delta must be <= 1")
            delta = float(delta)

            if eta is None:
                eta = 0.501
            elif eta < 0.5:
                raise TypeError("eta must be >= 0.5")

        if prec < 0:
            raise TypeError("precision prec must be >= 0")
        int(prec)

        if algorithm.startswith('NTL:'):
            A = sage.libs.ntl.all.mat_ZZ(self.nrows(),self.ncols(),map(ntl_ZZ,self.list()))

            if algorithm == "NTL:LLL":
                r, det2 = A.LLL(a,b, verbose=verb)
                det2 = ZZ(det2)
                try:
                    det = ZZ(det2.sqrt_approx())
                    self.cache("det", det)
                except TypeError:
                    pass
            elif algorithm == "NTL:LLL_FP":
                if use_givens:
                    r = A.G_LLL_FP(delta, verbose=verb)
                else:
                    r = A.LLL_FP(delta, verbose=verb)
            elif algorithm == "NTL:LLL_QP":
                if use_givens:
                    r = A.G_LLL_QP(delta, verbose=verb)
                else:
                    r = A.LLL_QP(delta, verbose=verb)
            elif algorithm == "NTL:LLL_XD":
                if use_givens:
                    r = A.G_LLL_XD(delta, verbose=verb)
                else:
                    r = A.LLL_XD(delta, verbose=verb)
            elif algorithm == "NTL:LLL_RR":
                if use_givens:
                    r = A.G_LLL_RR(delta, verbose=verb)
                else:
                    r = A.LLL_XD(delta, verbose=verb)
            else:
                raise TypeError("algorithm %s not supported"%algorithm)

            r = ZZ(r)

            R = <Matrix_integer_dense>self.new_matrix(entries=map(ZZ,A.list()))
            self.cache("rank",r)

        elif algorithm.startswith('fpLLL:'):

            A = sage.libs.fplll.fplll.FP_LLL(self)
            if algorithm == 'fpLLL:wrapper':
                A.wrapper(prec, eta, delta)
            elif algorithm == 'fpLLL:heuristic':
                if early_red:
                    A.heuristic_early_red(prec,eta,delta,fp)
                else:
                    A.heuristic(prec,eta,delta,fp)
            elif algorithm == 'fpLLL:fast':
                if early_red:
                    A.fast_early_red(prec,eta,delta)
                else:
                    A.fast(prec,eta,delta)
            elif algorithm == 'fpLLL:proved':
                A.proved(prec,eta,delta)
            else:
                raise TypeError("algorithm %s not supported"%algorithm)
            R = A._sage_()
        else:
            raise TypeError("algorithm %s not supported"%algorithm)

        verbose("LLL finished", tm)
        return R

    def is_LLL_reduced(self, delta=None, eta=None):
        r"""
        Return ``True`` if this lattice is
        `(\delta, \eta)`-LLL reduced. See ``self.LLL``
        for a definition of LLL reduction.

        INPUT:


        -  ``delta`` - parameter as described above (default:
           3/4)

        -  ``eta`` - parameter as described above (default:
           0.501)


        EXAMPLE::

            sage: A = random_matrix(ZZ, 10, 10)
            sage: L = A.LLL()
            sage: A.is_LLL_reduced()
            False
            sage: L.is_LLL_reduced()
            True
        """
        if eta is None:
            eta = 0.501
        if delta is None:
            delta = ZZ(3)/ZZ(4)

        if delta <= ZZ(1)/ZZ(4):
            raise TypeError("delta must be > 1/4")
        elif delta > 1:
            raise TypeError("delta must be <= 1")

        if eta < 0.5:
            raise TypeError("eta must be >= 0.5")

        # this is pretty slow
        import sage.modules.misc
        G, mu = sage.modules.misc.gram_schmidt(self.rows())
        #For any $i>j$, we have $|mu_{i, j}| <= \eta$
        for e in mu.list():
            if e.abs() > eta:
                return False

        #For any $i<d$, we have $\delta |b_i^*|^2 <= |b_{i+1}^* + mu_{i+1, i} b_i^* |^2$
        norms = [G[i].norm()**2 for i in range(len(G))]
        for i in xrange(1,self.nrows()):
            if norms[i] < (delta - mu[i,i-1]**2) * norms[i-1]:
                return False
        return True

    def prod_of_row_sums(self, cols):
        """
        Return the product of the sums of the entries in the submatrix of
        self with given columns.

        INPUT:


        -  ``cols`` - a list (or set) of integers representing
           columns of self.


        OUTPUT: an integer

        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.prod_of_row_sums([0,2])
            40
            sage: (1+3)*(4+6)
            40
            sage: a.prod_of_row_sums(set([0,2]))
            40
        """
        cdef Py_ssize_t c, row
        cdef mpz_t s, pr
        mpz_init(s)
        mpz_init(pr)

        mpz_set_si(pr, 1)
        for row from 0 <= row < self._nrows:
            mpz_set_si(s, 0)
            for c in cols:
                if c<0 or c >= self._ncols:
                    raise IndexError("matrix column index out of range")
                mpz_add(s, s, self._matrix[row][c])
            mpz_mul(pr, pr, s)
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, pr)
        mpz_clear(s)
        mpz_clear(pr)
        return z

    def _linbox_sparse(self):
        cdef Py_ssize_t i, j
        v = ['%s %s M'%(self._nrows, self._ncols)]
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if mpz_cmp_si(self._matrix[i][j], 0):
                    v.append('%s %s %s'%(i+1,j+1,self.get_unsafe(i,j)))
        v.append('0 0 0\n')
        return '\n'.join(v)

    def _linbox_dense(self):
        cdef Py_ssize_t i, j
        v = ['%s %s x'%(self._nrows, self._ncols)]
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                v.append(str(self.get_unsafe(i,j)))
        return ' '.join(v)

    def rational_reconstruction(self, N):
        """
        Use rational reconstruction to lift self to a matrix over the
        rational numbers (if possible), where we view self as a matrix
        modulo N.

        INPUT:


        -  ``N`` - an integer


        OUTPUT:


        -  ``matrix`` - over QQ or raise a ValueError


        EXAMPLES: We create a random 4x4 matrix over ZZ.

        ::

            sage: A = matrix(ZZ, 4, [4, -4, 7, 1, -1, 1, -1, -12, -1, -1, 1, -1, -3, 1, 5, -1])

        There isn't a unique rational reconstruction of it::

            sage: A.rational_reconstruction(11)
            Traceback (most recent call last):
            ...
            ValueError: Rational reconstruction of 4 (mod 11) does not exist.

        We throw in a denominator and reduce the matrix modulo 389 - it
        does rationally reconstruct.

        ::

            sage: B = (A/3 % 389).change_ring(ZZ)
            sage: B.rational_reconstruction(389) == A/3
            True

        TEST:

        Check that ticket #9345 is fixed::

            sage: A = random_matrix(ZZ, 3, 3)
            sage: A.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
        """
        import misc
        return misc.matrix_integer_dense_rational_reconstruction(self, N)

    def randomize(self, density=1, x=None, y=None, distribution=None, \
                  nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        The parameters are the same as the ones for the integer ring's
        ``random_element`` function.

        If ``x`` and ``y`` are given, randomized entries of this matrix have
        to be between ``x`` and ``y`` and have density 1.

        INPUT:

        -  ``self`` - a mutable matrix over ZZ

        -  ``density`` - a float between 0 and 1

        -  ``x, y`` - if not ``None``, these are passed to the
           ``ZZ.random_element`` function as the upper and lower endpoints in
           the  uniform distribution

        -  ``distribution`` - would also be passed into ``ZZ.random_element``
           if given

        -  ``nonzero`` - bool (default: ``False``); whether the new entries
           are guaranteed to be zero

        OUTPUT:

        -  None, the matrix is modified in-place

        EXAMPLES::

            sage: A = matrix(ZZ, 2,3, [1..6]); A
            [1 2 3]
            [4 5 6]
            sage: A.randomize()
            sage: A
            [-8  2  0]
            [ 0  1 -1]
            sage: A.randomize(x=-30,y=30)
            sage: A
            [  5 -19  24]
            [ 24  23  -9]
        """
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.check_mutability()
        self.clear_cache()

        cdef randstate rstate = current_randstate()

        cdef Py_ssize_t i, j, k, nc, num_per_row
        global state, ZZ

        cdef IntegerRing_class the_integer_ring = ZZ

        if not nonzero:
            # Original code, before adding the ``nonzero`` option.
            sig_on()
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    the_integer_ring._randomize_mpz(self._entries[i], x, y, \
                                                    distribution)
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random()%nc
                        the_integer_ring._randomize_mpz(self._matrix[i][k], \
                                                        x, y, distribution)
            sig_off()
        else:
            # New code, to implement the ``nonzero`` option.  Note that this
            # code is almost the same as above, the only difference being that
            # each entry is set until it's non-zero.
            sig_on()
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    while mpz_sgn(self._entries[i]) == 0:
                        the_integer_ring._randomize_mpz(self._entries[i], \
                            x, y, distribution)
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        while mpz_sgn(self._matrix[i][k]) == 0:
                            the_integer_ring._randomize_mpz(self._matrix[i][k],\
                                                            x, y, distribution)
            sig_off()

    #### Rank

    def rank(self):
        """
        Return the rank of this matrix.

        OUTPUT:


        -  ``nonnegative integer`` - the rank


        .. note::

           The rank is cached.

        ALGORITHM: First check if the matrix has maxim possible rank by
        working modulo one random prime. If not call LinBox's rank
        function.

        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.rank()
            2
            sage: a = matrix(ZZ,3,3,[1..9]); a
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: a.rank()
            2

        Here's a bigger example - the rank is of course still 2::

            sage: a = matrix(ZZ,100,[1..100^2]); a.rank()
            2
        """
        r = self.fetch('rank')
        if not r is None: return r

        if self._nrows <= 6 and self._ncols <= 6 and self.height().ndigits() <= 10000:
            r = self._rank_pari()
        # Can very quickly detect full rank by working modulo p.
        r = self._rank_modp()
        if r == self._nrows or r == self._ncols:
            self.cache('rank', r)
            return r
        # Detecting full rank didn't work -- use LinBox's general algorithm.
        r = self._rank_linbox()
        self.cache('rank', r)
        return r

    def _rank_linbox(self):
        """
        Compute the rank of this matrix using Linbox.
        """
        self._init_linbox()
        cdef unsigned long r
        sig_on()
        r = linbox.rank()
        sig_off()
        return Integer(r)

    def _rank_modp(self, p=46337):
        A = self._mod_int_c(p)
        return A.rank()

    #### Determinant

    def determinant(self, algorithm='default', proof=None, stabilize=2):
        r"""
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm``

          - ``'default'`` -- automatically determine which algorithm
                             to use depending on the matrix.

          - ``'padic'`` -  uses a p-adic / multimodular
            algorithm that relies on code in IML and linbox

          - ``'linbox'`` - calls linbox det (you *must* set
            proof=False to use this!)

          - ``'ntl'`` - calls NTL's det function

          - ``'pari'`` - uses PARI

        -  ``proof`` - bool or None; if None use
           proof.linear_algebra(); only relevant for the padic algorithm.

           .. note::

              It would be *VERY VERY* hard for det to fail even with
              proof=False.

        -  ``stabilize`` - if proof is False, require det to be
           the same for this many CRT primes in a row. Ignored if proof is
           True.


        ALGORITHM: The p-adic algorithm works by first finding a random
        vector v, then solving A\*x = v and taking the denominator
        `d`. This gives a divisor of the determinant. Then we
        compute `\det(A)/d` using a multimodular algorithm and the
        Hadamard bound, skipping primes that divide `d`.

        TIMINGS: This is perhaps the fastest implementation of determinants
        in the world. E.g., for a 500x500 random matrix with 32-bit entries
        on a core2 duo 2.6Ghz running OS X, Sage takes 4.12 seconds,
        whereas Magma takes 62.87 seconds (both with proof False). With
        proof=True on the same problem Sage takes 5.73 seconds. For another
        example, a 200x200 random matrix with 1-digit entries takes 4.18
        seconds in pari, 0.18 in Sage with proof True, 0.11 in Sage with
        proof False, and 0.21 seconds in Magma with proof True and 0.18 in
        Magma with proof False.

        EXAMPLES::

            sage: A = matrix(ZZ,8,8,[3..66])
            sage: A.determinant()
            0

        ::

            sage: A = random_matrix(ZZ,20,20)
            sage: D1 = A.determinant()
            sage: A._clear_cache()
            sage: D2 = A.determinant(algorithm='ntl')
            sage: D1 == D2
            True

        We have a special-case algorithm for 4 x 4 determinants::

            sage: A = matrix(ZZ,4,[1,2,3,4,4,3,2,1,0,5,0,1,9,1,2,3])
            sage: A.determinant()
            270

        Next we try the Linbox det. Note that we must have proof=False.

        ::

            sage: A = matrix(ZZ,5,[1,2,3,4,5,4,6,3,2,1,7,9,7,5,2,1,4,6,7,8,3,2,4,6,7])
            sage: A.determinant(algorithm='linbox')
            Traceback (most recent call last):
            ...
            RuntimeError: you must pass the proof=False option to the determinant command to use LinBox's det algorithm
            sage: A.determinant(algorithm='linbox',proof=False)
            -21
            sage: A._clear_cache()
            sage: A.determinant()
            -21

        A bigger example::

            sage: A = random_matrix(ZZ,30)
            sage: d = A.determinant()
            sage: A._clear_cache()
            sage: A.determinant(algorithm='linbox',proof=False) == d
            True

        TESTS:

        This shows that we can compute determinants for all sizes up to
        80.  The check that the determinant of a squared matrix is a
        square is a sanity check that the result is probably correct::

            sage: for s in [1..80]:  # long time (6s on sage.math, 2013)
            ...       M = random_matrix(ZZ, s)
            ...       d = (M*M).determinant()
            ...       assert d.is_square()
        """
        d = self.fetch('det')
        if d is not None:
            return d
        if not self.is_square():
            raise ValueError("self must be a square matrix")

        cdef Py_ssize_t n = self.nrows()

        cdef Integer det4x4
        if n <= 4:
            if n <= 3:
                # Use generic special cased code.
                return matrix_dense.Matrix_dense.determinant(self)
            else:
                det4x4 = ZZ(0)
                four_dim_det(det4x4.value, self._entries)
                return det4x4

        proof = get_proof_flag(proof, "linear_algebra")

        cdef double difficulty
        if algorithm == 'default':
            # These heuristics are by Jeroen Demeyer (#14007).  There
            # is no mathematics behind this, it was experimentally
            # observed to work well.  I tried various estimates for
            # the "difficulty" of a matrix, and this one worked best.
            # I tested matrices with entries uniformly distributed in
            # [0,n] as well as random_matrix(ZZ,s).
            #
            # linbox works sometimes better for large matrices with
            # mostly small entries, but it is never much faster than
            # padic (and it only works with proof=False), so we never
            # default to using linbox.
            difficulty = (self._log_avg_sq1() + 2.0) * (n * n)
            if difficulty <= 800:
                algorithm = 'pari'
            elif n <= 48 or (proof and n <= 72) or (proof and n <= 400 and difficulty <= 600000):
                algorithm = 'ntl'
            else:
                algorithm = 'padic'

        if algorithm == 'padic':
            import matrix_integer_dense_hnf
            d = matrix_integer_dense_hnf.det_padic(self, proof=proof, stabilize=stabilize)
        elif algorithm == 'linbox':
            if proof:
                raise RuntimeError("you must pass the proof=False option to the determinant command to use LinBox's det algorithm")
            d = self._det_linbox()
        elif algorithm == 'pari':
            d = self._det_pari()
        elif algorithm == 'ntl':
            d = self._det_ntl()
        else:
            raise TypeError("algorithm '%s' not known"%(algorithm))

        self.cache('det', d)
        return d

    def _det_linbox(self):
        """
        Compute the determinant of this matrix using Linbox.
        """
        self._init_linbox()
        sig_on()
        d = linbox.det()
        sig_off()
        return Integer(d)

    def _det_ntl(self):
        """
        Compute the determinant of this matrix using NTL.
        """
        sig_on()
        d = self._ntl_().determinant()
        sig_off()
        return Integer(d)

    #### Rational kernel, via IML
    def _rational_kernel_iml(self):
        """
        IML: Return the rational kernel of this matrix (acting from the
        left), considered as a matrix over QQ. I.e., returns a matrix K
        such that self\*K = 0, and the number of columns of K equals the
        nullity of self.

        AUTHORS:

        - William Stein
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.matrix_space(self._ncols, 0).zero_matrix()

        cdef long dim
        cdef mpz_t *mp_N
        time = verbose('computing null space of %s x %s matrix using IML'%(self._nrows, self._ncols))
        sig_on()
        dim = nullspaceMP (self._nrows, self._ncols, self._entries, &mp_N)
        sig_off()
        P = self.matrix_space(self._ncols, dim)

        # Now read the answer as a matrix.
        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        for i from 0 <= i < dim*self._ncols:
            mpz_init_set(M._entries[i], mp_N[i])
            mpz_clear(mp_N[i])
        sage_free(mp_N)

        verbose("finished computing null space", time)
        M._initialized = True
        return M

    def _invert_iml(self, use_nullspace=False, check_invertible=True):
        """
        Invert this matrix using IML. The output matrix is an integer
        matrix and a denominator.

        INPUT:


        -  ``self`` - an invertible matrix

        -  ``use_nullspace`` - (default: False): whether to
           use nullspace algorithm, which is slower, but doesn't require
           checking that the matrix is invertible as a precondition.

        -  ``check_invertible`` - (default: True) whether to
           check that the matrix is invertible.


        OUTPUT: A, d such that A\*self = d


        -  ``A`` - a matrix over ZZ

        -  ``d`` - an integer


        ALGORITHM: Uses IML's p-adic nullspace function.

        EXAMPLES::

            sage: a = matrix(ZZ,3,[1,2,5, 3,7,8, 2,2,1])
            sage: b, d = a._invert_iml(); b,d
            (
            [  9  -8  19]
            [-13   9  -7]
            [  8  -2  -1], 23
            )
            sage: a*b
            [23  0  0]
            [ 0 23  0]
            [ 0  0 23]

        AUTHORS:

        - William Stein
        """
        if self._nrows != self._ncols:
            raise TypeError("self must be a square matrix.")

        P = self.parent()
        time = verbose('computing inverse of %s x %s matrix using IML'%(self._nrows, self._ncols))
        if use_nullspace:
            A = self.augment(P.identity_matrix())
            K = A._rational_kernel_iml()
            d = -K[self._nrows,0]
            if K.ncols() != self._ncols or d == 0:
                raise ZeroDivisionError("input matrix must be nonsingular")
            B = K[:self._nrows]
            verbose("finished computing inverse using IML", time)
            return B, d
        else:
            if check_invertible and self.rank() != self._nrows:
                raise ZeroDivisionError("input matrix must be nonsingular")
            return self._solve_iml(P.identity_matrix(), right=True)

    def _solve_right_nonsingular_square(self, B, check_rank=True):
        r"""
        If self is a matrix `A` of full rank, then this function
        returns a vector or matrix `X` such that `A X = B`.
        If `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix. The base
        ring of `X` is the integers unless a denominator is needed
        in which case the base ring is the rational numbers.

        .. note::

           In Sage one can also write ``A  B`` for
           ``A.solve_right(B)``, i.e., Sage implements the "the
           MATLAB/Octave backslash operator".

        .. note::

           This is currently only implemented when A is square.

        INPUT:


        -  ``B`` - a matrix or vector

        -  ``check_rank`` - bool (default: True); if True
           verify that in fact the rank is full.


        OUTPUT: a matrix or vector over `\QQ`

        EXAMPLES::

            sage: a = matrix(ZZ, 2, [0, -1, 1, 0])
            sage: v = vector(ZZ, [2, 3])
            sage: a \ v
            (3, -2)

        Note that the output vector or matrix is always over
        `\QQ`.

        ::

            sage: parent(a\v)
            Vector space of dimension 2 over Rational Field

        We solve a bigger system where the answer is over the rationals.

        ::

            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: v = vector(ZZ, [1,2,3])
            sage: w = a \ v; w
            (2/15, -4/15, 7/15)
            sage: parent(w)
            Vector space of dimension 3 over Rational Field
            sage: a * w
            (1, 2, 3)

        We solve a system where the right hand matrix has multiple
        columns.

        ::

            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: b = matrix(ZZ, 3, 2, [1,5, 2, -3, 3, 0])
            sage: w = a \ b; w
            [ 2/15 -19/5]
            [-4/15 -27/5]
            [ 7/15 98/15]
            sage: a * w
            [ 1  5]
            [ 2 -3]
            [ 3  0]

        TESTS: We create a random 100x100 matrix and solve the
        corresponding system, then verify that the result is correct.
        (Note that this test is very risky without having a seeded
        random number generator!)

        ::

            sage: n = 100
            sage: a = random_matrix(ZZ,n)
            sage: v = vector(ZZ,n,range(n))
            sage: x = a \ v
            sage: a * x == v
            True
        """
        t = verbose('starting IML solve_right...')
        # It would probably be much better to rewrite linbox so it
        # throws an error instead of ** going into an infinite loop **
        # in the non-full rank case.  In any case, we do this for now,
        # since rank is very fast and infinite loops are evil.
        if check_rank and self.rank() < self.nrows():
            raise ValueError("self must be of full rank.")

        if not self.is_square():
            raise NotImplementedError("the input matrix must be square.")

        if is_Vector(B):
            if self.nrows() != B.degree():
                raise ValueError("number of rows of self must equal degree of B.")
        elif self.nrows() != B.nrows():
                raise ValueError("number of rows of self must equal number of rows of B.")

        if self.nrows() == 0:
            return B

        matrix = True
        C = B
        if not isinstance(B, Matrix_integer_dense):
            if is_Vector(B):
                matrix = False
                C = self.matrix_space(self.nrows(), 1)(B.list())
            else:
                raise NotImplementedError

        if C.ncols() >= 2*self.ncols():
            # likely much better to just invert then multiply
            X = self**(-1)*C
            verbose('finished solve_right (via inverse)', t)
            return X

        X, d = self._solve_iml(C, right=True)
        if d != 1:
            X = (1/d) * X
        if not matrix:
            # Convert back to a vector
            X = (X.base_ring() ** X.nrows())(X.list())
        verbose('finished solve_right via IML', t)
        return X

    def _solve_iml(self, Matrix_integer_dense B, right=True):
        """
        Let A equal self be a square matrix. Given B return an integer
        matrix C and an integer d such that self C\*A == d\*B if right is
        False or A\*C == d\*B if right is True.

        OUTPUT:


        -  ``C`` - integer matrix

        -  ``d`` - integer denominator


        EXAMPLES::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,3,4, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B,right=False); C
            [  6 -18 -15  27]
            [  0  24  24 -36]
            [  4 -12  -6  -2]

        ::

            sage: d
            12

        ::

            sage: C*A == d*B
            True

        ::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,4,3, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B)
            sage: C
            [ 12  40  28]
            [-12  -4  -4]
            [ -6 -25 -16]
            [ 12  34  16]

        ::

            sage: d
            12

        ::

            sage: A*C == d*B
            True

        Test wrong dimensions::

            sage: A = random_matrix(ZZ, 4, 4)
            sage: B = random_matrix(ZZ, 2, 3)
            sage: B._solve_iml(A)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
            sage: A._solve_iml(B, right=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of columns must match self's number of rows
            sage: A._solve_iml(B, right=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of rows must match self's number of columns

        Check that this can be interrupted properly (:trac:`15453`)::

            sage: A = random_matrix(ZZ, 2000, 2000)
            sage: B = random_matrix(ZZ, 2000, 2000)
            sage: t0 = walltime()
            sage: alarm(2); A._solve_iml(B)  # long time
            Traceback (most recent call last):
            ...
            AlarmInterrupt
            sage: t = walltime(t0)
            sage: t < 10 or t
            True

        ALGORITHM: Uses IML.

        AUTHORS:

        - Martin Albrecht
        """
        cdef int i
        cdef mpz_t *mp_N, mp_D
        cdef Matrix_integer_dense M
        cdef Integer D

        if self._nrows != self._ncols:
            # This is *required* by the IML function we call below.
            raise ValueError("self must be a square matrix")

        if self.nrows() == 1:
            return B, self[0,0]

        cdef SOLU_POS solu_pos

        if right:
            if self._ncols != B._nrows:
                raise ArithmeticError("B's number of rows must match self's number of columns")

            n = self._ncols
            m = B._ncols

            P = self.matrix_space(n, m)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)

            solu_pos = RightSolu

        else: # left
            if self._nrows != B._ncols:
                raise ArithmeticError("B's number of columns must match self's number of rows")

            n = self._ncols
            m = B._nrows

            P = self.matrix_space(m, n)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = m, ncols = n), Integer(1)

            solu_pos = LeftSolu

        sig_check()

        mp_N = <mpz_t *> sage_malloc( n * m * sizeof(mpz_t) )
        for i from 0 <= i < n * m:
            mpz_init(mp_N[i])
        mpz_init(mp_D)

        try:
            sig_on()
            nonsingSolvLlhsMM(solu_pos, n, m, self._entries, B._entries, mp_N, mp_D)
            sig_off()

            M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
            for i from 0 <= i < n*m:
                mpz_init_set(M._entries[i], mp_N[i])
            M._initialized = True

            D = PY_NEW(Integer)
            mpz_set(D.value, mp_D)

            return M, D
        finally:
            mpz_clear(mp_D)
            for i from 0 <= i < n*m:
                mpz_clear(mp_N[i])
            sage_free(mp_N)

    def _rational_echelon_via_solve(self):
        r"""
        Computes information that gives the reduced row echelon form (over
        QQ!) of a matrix with integer entries.

        INPUT:


        -  ``self`` - a matrix over the integers.


        OUTPUT:


        -  ``pivots`` - ordered list of integers that give the
           pivot column positions

        -  ``nonpivots`` - ordered list of the nonpivot column
           positions

        -  ``X`` - matrix with integer entries

        -  ``d`` - integer


        If you put standard basis vectors in order at the pivot columns,
        and put the matrix (1/d)\*X everywhere else, then you get the
        reduced row echelon form of self, without zero rows at the bottom.

        .. note::

           IML is the actual underlying `p`-adic solver that we
           use.

        AUTHORS:

        - William Stein

        ALGORITHM: I came up with this algorithm from scratch. As far as I
        know it is new. It's pretty simple, but it is ... (fast?!).

        Let A be the input matrix.


        #. Compute r = rank(A).

        #. Compute the pivot columns of the transpose `A^t` of
           `A`. One way to do this is to choose a random prime
           `p` and compute the row echelon form of `A^t`
           modulo `p` (if the number of pivot columns is not equal to
           `r`, pick another prime).

        #. Let `B` be the submatrix of `A` consisting of
           the rows corresponding to the pivot columns found in the previous
           step. Note that, aside from zero rows at the bottom, `B`
           and `A` have the same reduced row echelon form.

        #. Compute the pivot columns of `B`, again possibly by
           choosing a random prime `p` as in [2] and computing the
           Echelon form modulo `p`. If the number of pivot columns is
           not `r`, repeat with a different prime. Note - in this step
           it is possible that we mistakenly choose a bad prime `p`
           such that there are the right number of pivot columns modulo
           `p`, but they are at the wrong positions - e.g., imagine
           the augmented matrix `[pI|I]` - modulo `p` we would
           miss all the pivot columns. This is OK, since in the next step we
           would detect this, as the matrix we would obtain would not be in
           echelon form.

        #. Let `C` be the submatrix of `B` of pivot
           columns. Let `D` be the complementary submatrix of
           `B` of all all non-pivot columns. Use a `p`-adic
           solver to find the matrix `X` and integer `d` such
           that `C (1/d) X=D`. I.e., solve a bunch of linear systems
           of the form `Cx = v`, where the columns of `X` are
           the solutions.

        #. Verify that we had chosen the correct pivot columns. Inspect the
           matrix `X` obtained in step 5. If when used to construct
           the echelon form of `B`, `X` indeed gives a matrix
           in reduced row echelon form, then we are done - output the pivot
           columns, `X`, and `d`. To tell if `X` is
           correct, do the following: For each column of `X`
           (corresponding to non-pivot column `i` of `B`),
           read up the column of `X` until finding the first nonzero
           position `j`; then verify that `j` is strictly less
           than the number of pivot columns of `B` that are strictly
           less than `i`. Otherwise, we got the pivots of `B`
           wrong - try again starting at step 4, but with a different random
           prime.
        """
        if self._nrows == 0:
            pivots = []
            nonpivots = range(self._ncols)
            X = self.__copy__()
            d = Integer(1)
            return pivots, nonpivots, X, d
        elif self._ncols == 0:
            pivots = []
            nonpivots = []
            X = self.__copy__()
            d = Integer(1)
            return pivots, nonpivots, X, d

        from matrix_modn_dense_double import MAX_MODULUS
        A = self
        # Step 1: Compute the rank

        t = verbose('computing rank', level=2, caller_name='p-adic echelon')
        r = A.rank()
        verbose('done computing rank', level=2, t=t, caller_name='p-adic echelon')

        cdef randstate rstate = current_randstate()

        if r == self._nrows:
            # The input matrix already has full rank.
            B = A
        else:
            # Steps 2 and 3: Extract out a submatrix of full rank.
            i = 0
            while True:
                p = previous_prime(rstate.c_random() % (MAX_MODULUS-15000) + 10000)
                P = A._mod_int(p).transpose().pivots()
                if len(P) == r:
                    B = A.matrix_from_rows(P)
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError("Bug in _rational_echelon_via_solve in finding linearly independent rows.")

        _i = 0
        while True:
            _i += 1
            if _i == 50:
                raise RuntimeError("Bug in _rational_echelon_via_solve -- pivot columns keep being wrong.")

            # Step 4: Now we instead worry about computing the reduced row echelon form of B.
            i = 0
            while True:
                p = previous_prime(rstate.c_random() % (MAX_MODULUS-15000) + 10000)
                pivots = B._mod_int(p).pivots()
                if len(pivots) == r:
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError("Bug in _rational_echelon_via_solve in finding pivot columns.")

            # Step 5: Apply p-adic solver
            C = B.matrix_from_columns(pivots)
            pivots_ = set(pivots)
            non_pivots = [i for i in range(B.ncols()) if not i in pivots_]
            D = B.matrix_from_columns(non_pivots)
            t = verbose('calling IML solver', level=2, caller_name='p-adic echelon')
            X, d = C._solve_iml(D, right=True)
            t = verbose('finished IML solver', level=2, caller_name='p-adic echelon', t=t)

            # Step 6: Verify that we had chosen the correct pivot columns.
            pivots_are_right = True
            for z in range(X.ncols()):
                if not pivots_are_right:
                    break
                i = non_pivots[z]
                np = len([k for k in pivots if k < i])
                for j in reversed(range(X.nrows())):
                    if X[j,z] != 0:
                        if j < np:
                            break # we're good -- go on to next column of X
                        else:
                            pivots_are_right = False
                            break
            if pivots_are_right:
                break

        #end while


        # Finally, return the answer.
        return pivots, non_pivots, X, d

    def decomposition(self, **kwds):
        """
        Returns the decomposition of the free module on which this matrix A
        acts from the right (i.e., the action is x goes to x A), along with
        whether this matrix acts irreducibly on each factor. The factors
        are guaranteed to be sorted in the same way as the corresponding
        factors of the characteristic polynomial, and are saturated as ZZ
        modules.

        INPUT:


        -  ``self`` - a matrix over the integers

        -  ``**kwds`` - these are passed onto to the
           decomposition over QQ command.


        EXAMPLES::

            sage: t = ModularSymbols(11,sign=1).hecke_matrix(2)
            sage: w = t.change_ring(ZZ)
            sage: w.list()
            [3, -1, 0, -2]
        """
        F = self.charpoly().factor()
        if len(F) == 1:
            V = self.base_ring()**self.nrows()
            return decomp_seq([(V, F[0][1]==1)])

        A = self.change_ring(QQ)
        X = A.decomposition(**kwds)
        V = ZZ**self.nrows()
        if isinstance(X, tuple):
            D, E = X
            D = [(W.intersection(V), t) for W, t in D]
            E = [(W.intersection(V), t) for W, t in E]
            return decomp_seq(D), decomp_seq(E)
        else:
            return decomp_seq([(W.intersection(V), t) for W, t in X])

    def _add_row_and_maintain_echelon_form(self, row, pivots):
        """
        Assuming self is a full rank n x m matrix in reduced row Echelon
        form over ZZ and row is a vector of degree m, this function creates
        a new matrix that is the echelon form of self with row appended to
        the bottom.

        .. warning::

           It is assumed that self is in echelon form.

        INPUT:


        -  ``row`` - a vector of degree m over ZZ

        -  ``pivots`` - a list of integers that are the pivot
           columns of self.


        OUTPUT:


        -  ``matrix`` - a matrix of in reduced row echelon form
           over ZZ

        -  ``pivots`` - list of integers


        ALGORITHM: For each pivot column of self, we use the extended
        Euclidean algorithm to clear the column. The result is a new matrix
        B whose row span is the same as self.stack(row), and whose last row
        is 0 if and only if row is in the QQ-span of the rows of self. If
        row is not in the QQ-span of the rows of self, then row is nonzero
        and suitable to be inserted into the top n rows of A to form a new
        matrix that is in reduced row echelon form. We then clear that
        corresponding new pivot column.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, [1, 0, 110, 0, 3, 112, 0, 0, 221]); a
            [  1   0 110]
            [  0   3 112]
            [  0   0 221]
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1,2])
            (
            [1 0 0]
            [0 1 0]
            [0 0 1], [0, 1, 2]
            )
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[0,0,0]),[0,1,2])
            (
            [  1   0 110]
            [  0   3 112]
            [  0   0 221], [0, 1, 2]
            )
            sage: a = matrix(ZZ, 2, [1, 0, 110, 0, 3, 112])
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1])
            (
            [  1   0 110]
            [  0   1 219]
            [  0   0 545], [0, 1, 2]
            )
        """
        from sage.all import get_memory_usage
        cdef Py_ssize_t i, j, piv, n = self._nrows, m = self._ncols

        import constructor

        # 0. Base case
        if self.nrows() == 0:
            pos = row.nonzero_positions()
            if len(pos) > 0:
                pivots = [pos[0]]
                if row[pivots[0]] < 0:
                    row *= -1
            else:
                pivots = []
            return constructor.matrix([row]), pivots


        if row == 0:
            return self, pivots
        # 1. Create a new matrix that has row as the last row.
        row_mat = constructor.matrix(row)
        A = self.stack(row_mat)
        # 2. Working from the left, clear each column to put
        #    the resulting matrix back in echelon form.
        for i, p in enumerate(pivots):
            # p is the i-th pivot
            b = A[n,p]
            if not b:
                continue

            # (a). Take xgcd of pivot positions in last row and in ith
            # row.

            # TODO (optimize) -- change to use direct call to gmp and
            # no bounds checking!
            a = A[i,p]
            if not a:
                raise ZeroDivisionError("claimed pivot is not a pivot")
            if b % a == 0:
                # (b) Subtract a multiple of row i from row n.
                c = b // a
                if c:
                    for j in range(m):
                        A[n,j] -= c * A[i,j]
            else:
                # (b). More elaborate.
                #  Replace the ith row by s*A[i] + t*A[n], which will
                # have g in the i,p position, and replace the last row by
                # (b//g)*A[i] - (a//g)*A[n], which will have 0 in the i,p
                # position.
                g, s, t = a.xgcd(b)
                if not g:
                    raise ZeroDivisionError("claimed pivot is not a pivot (got a 0 gcd)")

                row_i = A.row(i)
                row_n = A.row(n)

                ag = a//g; bg = b//g

                new_top = s*row_i  +  t*row_n
                new_bot = bg*row_i - ag*row_n


                # OK -- now we have to make sure the top part of the matrix
                # but with row i replaced by
                #     r = s*row_i[j]  +  t*row_n[j]
                # is put in rref.  We do this by recursively calling this
                # function with the top part of A (all but last row) and the
                # row r.

                zz = range(A.nrows()-1)
                del zz[i]
                top_mat = A.matrix_from_rows(zz)
                new_pivots = list(pivots)
                del new_pivots[i]

                top_mat, pivots = top_mat._add_row_and_maintain_echelon_form(new_top, new_pivots)
                w = top_mat._add_row_and_maintain_echelon_form(new_bot, pivots)
                return w
        # 3. If it turns out that the last row is nonzero,
        #    insert last row in A sliding other rows down.
        v = A.row(n)
        new_pivots = list(pivots)
        if v != 0:
            i = v.nonzero_positions()[0]
            assert not (i in pivots), 'WARNING: bug in add_row -- i (=%s) should not be a pivot'%i

            # If pivot entry is negative negate this row.
            if v[i] < 0:
                v = -v

            # Determine where the last row should be inserted.
            new_pivots.append(i)
            new_pivots.sort()
            import bisect
            j = bisect.bisect(pivots, i)
            # The new row should go *before* row j, so it becomes row j
            A = A.insert_row(j, v)
        try:
            _clear_columns(A, new_pivots, A.nrows())
        except RuntimeError:
            raise ZeroDivisionError("mistake in claimed pivots")
        if A.row(A.nrows() - 1) == 0:
            A = A.matrix_from_rows(range(A.nrows()-1))
        return A, new_pivots

    #####################################################################################
    # Hermite form modulo D
    # This code below is by E. Burcin.  Thanks!
    #####################################################################################
    cdef _new_uninitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the integers with the given number of rows
        and columns. All memory is allocated for this matrix, but its
        entries have not yet been filled in.
        """
        P = self._parent.matrix_space(nrows, ncols)
        return Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)

    def _hnf_mod(self, D):
        """
        INPUT:


        -  ``D`` - a small integer that is assumed to be a
           multiple of 2\*det(self)


        OUTPUT:


        -  ``matrix`` - the Hermite normal form of self.
        """
        t = verbose('hermite mod %s'%D, caller_name='matrix_integer_dense')
        cdef Matrix_integer_dense res = self._new_uninitialized_matrix(self._nrows, self._ncols)
        self._hnf_modn(res, D)
        verbose('finished hnf mod', t, caller_name='matrix_integer_dense')
        return res

    cdef int _hnf_modn(Matrix_integer_dense self, Matrix_integer_dense res,
            mod_int det) except -1:
        """
        Puts self into HNT form modulo det. Changes self in place.
        """
        cdef long long *res_l
        cdef Py_ssize_t i,j,k
        res_l = self._hnf_modn_impl(det, self._nrows, self._ncols)
        k = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpz_init_set_si(res._matrix[i][j], res_l[k])
                k += 1
        res._initialized = True
        sage_free(res_l)


    cdef long long* _hnf_modn_impl(Matrix_integer_dense self, mod_int det,
            Py_ssize_t nrows, Py_ssize_t ncols) except NULL:
        cdef long long *res, *T_ent, **res_rows, **T_rows, *B
        cdef Py_ssize_t i, j, k
        cdef long long R, mod, T_i_i, T_j_i, c1, c2, q, t
        cdef int u, v, d
        cdef mpz_t m

        # allocate memory for result matrix
        res = <long long*> sage_malloc(sizeof(long long)*ncols*nrows)
        if res == NULL:
            raise MemoryError("out of memory allocating a matrix")
        res_rows = <long long**> sage_malloc(sizeof(long long*)*nrows)
        if res_rows == NULL:
            sage_free(res)
            raise MemoryError("out of memory allocating a matrix")

        # allocate memory for temporary matrix
        T_ent = <long long*> sage_malloc(sizeof(long long)*ncols*nrows)
        if T_ent == NULL:
            sage_free(res)
            sage_free(res_rows)
            raise MemoryError("out of memory allocating a matrix")
        T_rows = <long long**> sage_malloc(sizeof(long long*)*nrows)
        if T_rows == NULL:
            sage_free(res)
            sage_free(res_rows)
            sage_free(T_ent)
            raise MemoryError("out of memory allocating a matrix")

        # allocate memory for temporary row vector
        B = <long long*>sage_malloc(sizeof(long long)*nrows)
        if B == NULL:
            sage_free(res)
            sage_free(res_rows)
            sage_free(T_ent)
            sage_free(T_rows)
            raise MemoryError("out of memory allocating a matrix")

        # initialize the row pointers
        k = 0
        for i from 0 <= i < nrows:
            res_rows[i] = res + k
            T_rows[i] = T_ent + k
            k += nrows


        mpz_init(m)
        # copy entries from self to temporary matrix
        k = 0
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                mpz_mod_ui(m, self._matrix[i][j], det)
                T_ent[k] = mpz_get_si(m)
                k += 1
        mpz_clear(m)


        # initialize variables
        i = 0
        j = 0
        R = det

        while 1:
            if j == nrows-1:
                T_i_i = T_rows[i][i]
                d = ai.c_xgcd_int(T_i_i, R, &u, &v)
                for k from 0 <= k < i:
                    res_rows[i][k] = 0
                for k from i <= k < ncols:
                    t = (u*T_rows[i][k])%R
                    if t < 0:
                        t += R
                    res_rows[i][k] = t
                if res_rows[i][i] == 0:
                    res_rows[i][i] = R
                d = res_rows[i][i]
                for j from 0 <= j < i:
                    q = res_rows[j][i]/d
                    for k from i <= k < ncols:
                        u = (res_rows[j][k] - q*res_rows[i][k])%R
                        if u < 0:
                            u += R
                        res_rows[j][k] = u

                R = R/d
                i += 1
                j = i
                if i == nrows :
                    break # return res
                if T_rows[i][i] == 0:
                    T_rows[i][i] = R
                continue


            j += 1
            if T_rows[j][i] == 0:
                continue

            T_i_i = T_rows[i][i]
            T_j_i = T_rows[j][i]
            d = ai.c_xgcd_int(T_i_i , T_j_i, &u, &v)
            if d != T_i_i:
                for k from i <= k < ncols:
                    B[k] = (u*T_rows[i][k] + v*T_rows[j][k])
            c1 = T_i_i/d
            c2 = -T_j_i/d
            for k from i <= k < ncols:
                T_rows[j][k] = (c1*T_rows[j][k] + c2*T_rows[i][k])%R
            if d != T_i_i:
                for k from i <= k < ncols:
                    T_rows[i][k] = B[k]%R

        sage_free(B)
        sage_free(res_rows)
        sage_free(T_ent)
        sage_free(T_rows)
        return res


    #####################################################################################
    # operations with matrices
    #####################################################################################
    def stack(self, bottom, subdivide=False):
        r"""
        Return the matrix self on top of bottom: [ self ] [ bottom ]

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 3, range(6))
            sage: N = Matrix(ZZ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]

        A vector may be stacked below a matrix. ::

            sage: A = matrix(ZZ, 2, 4, range(8))
            sage: v = vector(ZZ, 4, range(4))
            sage: A.stack(v)
            [0 1 2 3]
            [4 5 6 7]
            [0 1 2 3]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``bottom``.  For more details about how subdivisions
        are managed when stacking, see
        :meth:`sage.matrix.matrix1.Matrix.stack`.  ::

            sage: A = matrix(ZZ, 3, 4, range(12))
            sage: B = matrix(ZZ, 2, 4, range(8))
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [-----------]
            [ 0  1  2  3]
            [ 4  5  6  7]

        TESTS:

        Stacking a dense matrix atop a sparse one should work::

            sage: M = Matrix(ZZ, 2, 3, range(6))
            sage: M.is_sparse()
            False
            sage: N = diagonal_matrix([10,11,12], sparse=True)
            sage: N.is_sparse()
            True
            sage: P = M.stack(N); P
            [ 0  1  2]
            [ 3  4  5]
            [10  0  0]
            [ 0 11  0]
            [ 0  0 12]
            sage: P.is_sparse()
            False
        """
        if hasattr(bottom, '_vector_'):
            bottom = bottom.row()
        if self._ncols != bottom.ncols():
            raise TypeError("number of columns must be the same")
        if not (self._base_ring is bottom.base_ring()):
            bottom = bottom.change_ring(self._base_ring)
        cdef Matrix_integer_dense other = bottom.dense_matrix()
        cdef Matrix_integer_dense M
        M = self.new_matrix(nrows = self._nrows + other._nrows, ncols = self.ncols())
        cdef Py_ssize_t i, k
        k = self._nrows * self._ncols
        for i from 0 <= i < k:
            mpz_set(M._entries[i], self._entries[i])
        for i from 0 <= i < other._nrows * other._ncols:
            mpz_set(M._entries[k + i], other._entries[i])
        if subdivide:
            M._subdivide_on_stack(self, other)
        return M

    def augment(self, right, subdivide=False):
        r"""
        Returns a new matrix formed by appending the matrix
        (or vector) ``right`` on the right side of ``self``.

        INPUT:

        - ``right`` - a matrix, vector or free module element, whose
          dimensions are compatible with ``self``.

        - ``subdivide`` - default: ``False`` - request the resulting
          matrix to have a new subdivision, separating ``self`` from ``right``.

        OUTPUT:

        A new matrix formed by appending ``right`` onto the right side of ``self``.
        If ``right`` is a vector (or free module element) then in this context
        it is appropriate to consider it as a column vector.  (The code first
        converts a vector to a 1-column matrix.)

        EXAMPLES::

            sage: A = matrix(ZZ, 4, 5, range(20))
            sage: B = matrix(ZZ, 4, 3, range(12))
            sage: A.augment(B)
            [ 0  1  2  3  4  0  1  2]
            [ 5  6  7  8  9  3  4  5]
            [10 11 12 13 14  6  7  8]
            [15 16 17 18 19  9 10 11]

        A vector may be augmented to a matrix. ::

            sage: A = matrix(ZZ, 3, 5, range(15))
            sage: v = vector(ZZ, 3, range(3))
            sage: A.augment(v)
            [ 0  1  2  3  4  0]
            [ 5  6  7  8  9  1]
            [10 11 12 13 14  2]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``right``.  For more details about how subdivisions
        are managed when augmenting, see
        :meth:`sage.matrix.matrix1.Matrix.augment`.  ::

            sage: A = matrix(ZZ, 3, 5, range(15))
            sage: B = matrix(ZZ, 3, 3, range(9))
            sage: A.augment(B, subdivide=True)
            [ 0  1  2  3  4| 0  1  2]
            [ 5  6  7  8  9| 3  4  5]
            [10 11 12 13 14| 6  7  8]

        Errors are raised if the sizes are incompatible. ::

            sage: A = matrix(ZZ, [[1, 2],[3, 4]])
            sage: B = matrix(ZZ, [[10, 20], [30, 40], [50, 60]])
            sage: A.augment(B)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, not 2 != 3
        """
        if hasattr(right, '_vector_'):
            right = right.column()
        if self._nrows != right.nrows():
            raise TypeError('number of rows must be the same, not {0} != {1}'.format(self._nrows, right.nrows()))
        if not (self._base_ring is right.base_ring()):
            right = right.change_ring(self._base_ring)

        cdef Matrix_integer_dense other = right.dense_matrix()
        m = self._nrows
        ns, na = self._ncols, other._ncols
        n = ns + na

        cdef Matrix_integer_dense Z
        Z = self.new_matrix(nrows = m, ncols = n)
        cdef Py_ssize_t i, j, p, qs, qa
        p, qs, qa = 0, 0, 0
        for i from 0 <= i < m:
          for j from 0 <= j < ns:
            mpz_set(Z._entries[p], self._entries[qs])
            p = p + 1
            qs = qs + 1
          for j from 0 <= j < na:
            mpz_set(Z._entries[p], other._entries[qa])
            p = p + 1
            qa = qa + 1
        if subdivide:
          Z._subdivide_on_augment(self, other)
        return Z

    def insert_row(self, Py_ssize_t index, row):
        """
        Create a new matrix from self with.

        INPUT:

        - ``index`` - integer

        - ``row`` - a vector

        EXAMPLES::

            sage: X = matrix(ZZ,3,range(9)); X
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: X.insert_row(1, [1,5,-10])
            [  0   1   2]
            [  1   5 -10]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(0, [1,5,-10])
            [  1   5 -10]
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(3, [1,5,-10])
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            [  1   5 -10]
        """
        cdef Matrix_integer_dense res = self._new_uninitialized_matrix(self._nrows + 1, self._ncols)
        cdef Py_ssize_t j, k
        cdef Integer z
        if index < 0:
            raise ValueError("index must be nonnegative")
        if index > self._nrows:
            raise ValueError("index must be less than number of rows")
        for j from 0 <= j < self._ncols * index:
            mpz_init_set(res._entries[j], self._entries[j])

        k = 0
        for j from self._ncols * index <= j < self._ncols * (index+1):
            z = row[k]
            mpz_init_set(res._entries[j], z.value)
            k += 1

        for j from self._ncols * (index+1) <= j < (self._nrows + 1)*self._ncols:
            mpz_init_set(res._entries[j], self._entries[j - self._ncols])

        res._initialized = True
        return res

    def _delete_zero_columns(self):
        """
        Return matrix obtained from self by deleting all zero columns along
        with the positions of those columns.

        OUTPUT: matrix list of integers

        EXAMPLES::

            sage: a = matrix(ZZ, 2,3, [1,0,3,-1,0,5]); a
            [ 1  0  3]
            [-1  0  5]
            sage: a._delete_zero_columns()
            (
            [ 1  3]
            [-1  5], [1]
            )
        """
        C = self.columns()
        zero_cols = [i for i,v in enumerate(self.columns()) if v.is_zero()]
        s = set(zero_cols)
        nonzero_cols = [i for i in range(self.ncols()) if not (i in s)]
        return self.matrix_from_columns(nonzero_cols), zero_cols

    def _insert_zero_columns(self, cols):
        """
        Return matrix obtained by self by inserting zero columns so that
        the columns with positions specified in cols are all 0.

        INPUT:

        -  ``cols`` - list of nonnegative integers

        OUTPUT: matrix

        EXAMPLES::

            sage: a = matrix(ZZ, 2,3, [1,0,3,-1,0,5]); a
            [ 1  0  3]
            [-1  0  5]
            sage: b, cols = a._delete_zero_columns()
            sage: b
            [ 1  3]
            [-1  5]
            sage: cols
            [1]
            sage: b._insert_zero_columns(cols)
            [ 1  0  3]
            [-1  0  5]
        """
        if len(cols) == 0:
            return self
        cdef Py_ssize_t i, c, r, nc = max(self._ncols + len(cols), max(cols)+1)
        cdef Matrix_integer_dense A = self.new_matrix(self._nrows, nc)
        # Now fill in the entries of A that come from self.
        cols_set = set(cols)
        cols_ins = [j for j in range(nc) if j not in cols_set]
        for r from 0 <= r < self._nrows:
            i = 0
            for c in cols_ins:
                # The following does this quickly: A[r,c] = self[r,i]
                mpz_set(A._matrix[r][c], self._matrix[r][i])
                i += 1
        return A

    def _factor_out_common_factors_from_each_row(self):
        """
        Very very quickly modifies self so that the gcd of the entries in
        each row is 1 by dividing each row by the common gcd.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, [-9, 3, -3, -36, 18, -5, -40, -5, -5, -20, -45, 15, 30, -15, 180])
            sage: a
            [ -9   3  -3 -36  18]
            [ -5 -40  -5  -5 -20]
            [-45  15  30 -15 180]
            sage: a._factor_out_common_factors_from_each_row()
            sage: a
            [ -3   1  -1 -12   6]
            [ -1  -8  -1  -1  -4]
            [ -3   1   2  -1  12]
        """
        self.check_mutability()

        cdef mpz_t g
        mpz_init(g)
        cdef Py_ssize_t i, j
        cdef mpz_t* row

        for i from 0 <= i < self._nrows:
            mpz_set_ui(g, 0)
            row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_gcd(g, g, row[j])
                if mpz_cmp_ui(g, 1) == 0:
                    break
            if mpz_cmp_ui(g, 1) != 0:
                # divide through row
                for j from 0 <= j < self._ncols:
                    mpz_divexact(row[j], row[j], g)
        mpz_clear(g)

    def gcd(self):
        """
        Return the gcd of all entries of self; very fast.

        EXAMPLES::

            sage: a = matrix(ZZ,2, [6,15,-6,150])
            sage: a.gcd()
            3
        """
        cdef Integer g = Integer(0)
        cdef Py_ssize_t i, j
        cdef mpz_t* row

        for i from 0 <= i < self._nrows:
            row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_gcd(g.value, g.value, row[j])
                if mpz_cmp_ui(g.value, 1) == 0:
                    return g
        return g

    def _change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this matrix
        into the given ring.

        EXAMPLES::

            sage: a = matrix(ZZ,2,[1,-7,3,5])
            sage: a._change_ring(RDF)
            [ 1.0 -7.0]
            [ 3.0  5.0]
        """
        if ring == RDF:
            import change_ring
            return change_ring.integer_to_real_double_dense(self)
        else:
            raise NotImplementedError

    def _singular_(self, singular=None):
        r"""
        Return Singular representation of this integer matrix.

        INPUT:


        -  ``singular`` - Singular interface instance (default:
           None)


        EXAMPLE::

            sage: A = random_matrix(ZZ,3,3)
            sage: As = singular(A); As
            -8     2     0
            0     1    -1
            2     1   -95
            sage: As.type()
            'intmat'
        """
        if singular is None:
            from sage.interfaces.singular import singular as singular_default
            singular = singular_default

        name = singular._next_var_name()
        values = str(self.list())[1:-1]
        singular.eval("intmat %s[%d][%d] = %s"%(name, self.nrows(), self.ncols(), values))

        from sage.interfaces.singular import SingularElement
        return SingularElement(singular, 'foobar', name, True)

    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:

        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.

        ::

            sage: A = matrix(ZZ,2,3,xrange(6))
            sage: type(A)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: B = A.transpose()
            sage: print B
            [0 3]
            [1 4]
            [2 5]
            sage: print A
            [0 1 2]
            [3 4 5]

        ``.T`` is a convenient shortcut for the transpose::

            sage: A.T
            [0 3]
            [1 4]
            [2 5]

        ::

            sage: A.subdivide(None, 1); A
            [0|1 2]
            [3|4 5]
            sage: A.transpose()
            [0 3]
            [---]
            [1 4]
            [2 5]
        """
        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense,
                                          self._parent.matrix_space(self._ncols,self._nrows),
                                          0,False,False)
        cdef Py_ssize_t i,j
        sig_on()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpz_init_set(A._matrix[j][i], self._matrix[i][j])
        sig_off()
        A._initialized = True
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide(col_divs, row_divs)
        return A

    def antitranspose(self):
        """
        Returns the antitranspose of self, without changing self.

        EXAMPLES::

            sage: A = matrix(2,3,range(6))
            sage: type(A)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: A.antitranspose()
            [5 2]
            [4 1]
            [3 0]
            sage: A
            [0 1 2]
            [3 4 5]

            sage: A.subdivide(1,2); A
            [0 1|2]
            [---+-]
            [3 4|5]
            sage: A.antitranspose()
            [5|2]
            [-+-]
            [4|1]
            [3|0]
        """
        nr , nc = (self._nrows, self._ncols)

        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense,
                                          self._parent.matrix_space(nc,nr),
                                          0,False,False)

        cdef Py_ssize_t i,j
        cdef Py_ssize_t ri,rj # reversed i and j
        sig_on()
        ri = nr
        for i from 0 <= i < nr:
            rj = nc
            ri =  ri-1
            for j from 0 <= j < nc:
                rj = rj-1
                mpz_init_set(A._matrix[rj][ri], self._matrix[i][j])
        sig_off()
        A._initialized = True

        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide([nc - t for t in reversed(col_divs)],
                        [nr - t for t in reversed(row_divs)])
        return A

    def _pari_(self):
        """
        Return PARI C-library version of this matrix.

        EXAMPLES::

            sage: a = matrix(ZZ,2,2,[1,2,3,4])
            sage: a._pari_()
            [1, 2; 3, 4]
            sage: pari(a)
            [1, 2; 3, 4]
            sage: type(pari(a))
            <type 'sage.libs.pari.gen.gen'>
        """
        return pari.integer_matrix(self._matrix, self._nrows, self._ncols, 0)

    def _det_pari(self, int flag=0):
        """
        Determinant of this matrix using Gauss-Bareiss. If (optional)
        flag is set to 1, use classical Gaussian elimination.

        For efficiency purposes, this det is computed entirely on the
        PARI stack then the PARI stack is cleared.  This function is
        most useful for very small matrices.

        EXAMPLES::

            sage: matrix(ZZ,3,[1..9])._det_pari()
            0
            sage: matrix(ZZ,3,[1..9])._det_pari(1)
            0
        """
        pari_catch_sig_on()
        cdef GEN d = det0(pari_GEN(self), flag)
        # now convert d to a Sage integer e
        cdef Integer e = Integer()
        t_INT_to_ZZ(e.value, d)
        pari.clear_stack()
        return e

    def _rank_pari(self):
        """
        Rank of this matrix, computed using PARI.  The computation is
        done entirely on the PARI stack, then the PARI stack is
        cleared.  This function is most useful for very small
        matrices.

        EXAMPLES::
            sage: matrix(ZZ,3,[1..9])._rank_pari()
            2
        """
        pari_catch_sig_on()
        cdef long r = rank(pari_GEN(self))
        pari.clear_stack()
        return r

    def _hnf_pari(self, int flag=0, bint include_zero_rows=True):
        """
        Hermite form of this matrix, computed using PARI.  The
        computation is done entirely on the PARI stack, then the PARI
        stack is cleared.  This function is only useful for small
        matrices, and can crash on large matrices (e.g., if the PARI
        stack overflows).

        INPUT:

        - ``flag`` -- 0 (default), 1, 3 or 4 (see docstring for
          gp.mathnf).

        - ``include_zero_rows`` -- boolean. if False, do not include
          any of the zero rows at the bottom of the matrix in the
          output.

        .. NOTE::

            In no cases is the transformation matrix returned by this
            function.

        EXAMPLES::

            sage: matrix(ZZ,3,[1..9])._hnf_pari()
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(1)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(3)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(4)
            [1 2 3]
            [0 3 6]
            [0 0 0]

        Check that ``include_zero_rows=False`` works correctly::

            sage: matrix(ZZ,3,[1..9])._hnf_pari(0, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(1, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(3, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(4, include_zero_rows=False)
            [1 2 3]
            [0 3 6]

        Check that :trac:`12346` is fixed::

            sage: pari('mathnf(Mat([0,1]), 4)')
            [Mat(1), [1, 0; 0, 1]]
        """
        cdef GEN A
        pari_catch_sig_on()
        A = pari._new_GEN_from_mpz_t_matrix_rotate90(self._matrix, self._nrows, self._ncols)
        cdef GEN H = mathnf0(A, flag)
        B = self.extract_hnf_from_pari_matrix(H, flag, include_zero_rows)
        pari.clear_stack()  # This calls pari_catch_sig_off()
        return B


    def _hnf_pari_big(self, int flag=0, bint include_zero_rows=True):
        """
        Hermite form of this matrix, computed using PARI.

        INPUT:

        - ``flag`` -- 0 (default), 1, 3 or 4 (see docstring for
          gp.mathnf).

        - ``include_zero_rows`` -- boolean. if False, do not include
          any of the zero rows at the bottom of the matrix in the
          output.

        .. NOTE::

            In no cases is the transformation matrix returned by this
            function.

        EXAMPLES::

            sage: a = matrix(ZZ,3,3,[1..9])
            sage: a._hnf_pari_big(flag=0, include_zero_rows=True)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: a._hnf_pari_big(flag=1, include_zero_rows=True)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: a._hnf_pari_big(flag=3, include_zero_rows=True)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: a._hnf_pari_big(flag=4, include_zero_rows=True)
            [1 2 3]
            [0 3 6]
            [0 0 0]

        Check that ``include_zero_rows=False`` works correctly::

            sage: matrix(ZZ,3,[1..9])._hnf_pari_big(0, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari_big(1, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari_big(3, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari_big(4, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
        """
        cdef gen H = pari.integer_matrix(self._matrix, self._nrows, self._ncols, 1)
        H = H.mathnf(flag)
        pari_catch_sig_on()
        B = self.extract_hnf_from_pari_matrix(H.g, flag, include_zero_rows)
        pari.clear_stack()  # This calls pari_catch_sig_off()
        return B

    cdef extract_hnf_from_pari_matrix(self, GEN H, int flag, bint include_zero_rows):
        # Throw away the transformation matrix (yes, we should later
        # code this to keep track of it).
        if flag > 0:
            H = gel(H,1)

        # Figure out how many columns we got back.
        cdef Py_ssize_t H_nc = glength(H)  # number of columns
        # Now get the resulting Hermite form matrix back to Sage, suitably re-arranged.
        cdef Matrix_integer_dense B
        if include_zero_rows:
            B = self.new_matrix()
        else:
            B = self.new_matrix(nrows=H_nc)
        for i in range(self._ncols):
            for j in range(H_nc):
                t_INT_to_ZZ(B._matrix[j][self._ncols-i-1], gcoeff(H, i+1, H_nc-j))
        return B

cdef inline GEN pari_GEN(Matrix_integer_dense B):
    r"""
    Create the PARI GEN object on the stack defined by the integer
    matrix B. This is used internally by the function for conversion
    of matrices to PARI.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN A = pari._new_GEN_from_mpz_t_matrix(B._matrix, B._nrows, B._ncols)
    return A


    #####################################################################################

cdef _clear_columns(Matrix_integer_dense A, pivots, Py_ssize_t n):
    # Clear all columns
    cdef Py_ssize_t i, k, p, l, m = A._ncols
    cdef mpz_t** matrix = A._matrix
    cdef mpz_t c, t
    sig_on()
    mpz_init(c)
    mpz_init(t)
    for i from 0 <= i < len(pivots):
        p = pivots[i]
        for k from 0 <= k < n:
            if k != i:
                if mpz_cmp_si(matrix[k][p],0):
                    mpz_fdiv_q(c, matrix[k][p], matrix[i][p])
                    # subtract off c*v from row k; resulting A[k,i] entry will be < b, hence in Echelon form.
                    for l from 0 <= l < m:
                        mpz_mul(t, c, matrix[i][l])
                        mpz_sub(matrix[k][l], matrix[k][l], t)
    mpz_clear(c)
    mpz_clear(t)
    sig_off()

###############################################################





cpdef _lift_crt(Matrix_integer_dense M, residues, moduli=None):
    """
    INPUT:

    - ``M`` -- A ``Matrix_integer_dense``. Will be modified to hold
      the output.

    - ``residues`` -- a list of ``Matrix_modn_dense_template``. The
      matrix to reconstruct modulo primes.

    OUTPUT:

    The matrix whose reductions modulo primes are the input
    ``residues``.

    TESTS::

        sage: from sage.matrix.matrix_integer_dense import _lift_crt
        sage: T1 = Matrix(Zmod(5), 4, 4, [1, 4, 4, 0, 2, 0, 1, 4, 2, 0, 4, 1, 1, 4, 0, 3])
        sage: T2 = Matrix(Zmod(7), 4, 4, [1, 4, 6, 0, 2, 0, 1, 2, 4, 0, 6, 6, 1, 6, 0, 5])
        sage: T3 = Matrix(Zmod(11), 4, 4, [1, 4, 10, 0, 2, 0, 1, 9, 8, 0, 10, 6, 1, 10, 0, 9])
        sage: _lift_crt(Matrix(ZZ, 4, 4), [T1, T2, T3])
        [ 1  4 -1  0]
        [ 2  0  1  9]
        [-3  0 -1  6]
        [ 1 -1  0 -2]

        sage: from sage.ext.multi_modular import MultiModularBasis
        sage: mm = MultiModularBasis([5,7,11])
        sage: _lift_crt(Matrix(ZZ, 4, 4), [T1, T2, T3], mm)
        [ 1  4 -1  0]
        [ 2  0  1  9]
        [-3  0 -1  6]
        [ 1 -1  0 -2]

    The modulus must be smaller than the maximum for the multi-modular
    reconstruction (using ``mod_int``) and also smaller than the limit
    for ``Matrix_modn_dense_double`` to be able to represent the
    ``residues`` ::

        sage: from sage.ext.multi_modular import MAX_MODULUS as MAX_multi_modular
        sage: from sage.matrix.matrix_modn_dense_double import MAX_MODULUS as MAX_modn_dense_double
        sage: MAX_MODULUS = min(MAX_multi_modular, MAX_modn_dense_double)
        sage: p0 = previous_prime(MAX_MODULUS)
        sage: p1 = previous_prime(p0)
        sage: mmod = [matrix(GF(p0), [[-1, 0, 1, 0, 0, 1, 1, 0, 0, 0, p0-1, p0-2]]),
        ....:         matrix(GF(p1), [[-1, 0, 1, 0, 0, 1, 1, 0, 0, 0, p1-1, p1-2]])]
        sage: _lift_crt(Matrix(ZZ, 1, 12), mmod)
        [-1  0  1  0  0  1  1  0  0  0 -1 -2]
    """

    cdef size_t n, i, j, k
    cdef Py_ssize_t nr, nc

    n = len(residues)
    if n == 0:   # special case: obviously residues[0] wouldn't make sense here.
        return M
    nr = residues[0].nrows()
    nc = residues[0].ncols()

    if moduli is None:
        moduli = MultiModularBasis([m.base_ring().order() for m in residues])
    else:
        if len(residues) != len(moduli):
            raise IndexError("Number of residues (%s) does not match number of moduli (%s)"%(len(residues), len(moduli)))

    cdef MultiModularBasis mm
    mm = moduli

    for b in residues:
        if not is_Matrix_modn_dense(b):
            raise TypeError("Can only perform CRT on list of matrices mod n.")

    cdef mod_int **row_list
    row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
    if row_list == NULL:
        raise MemoryError("out of memory allocating multi-modular coefficient list")

    sig_on()
    for k in range(n):
        row_list[k] = <mod_int *>sage_malloc(sizeof(mod_int) * nc)
        if row_list[k] == NULL:
            raise MemoryError("out of memory allocating multi-modular coefficient list")

    for i in range(nr):
        for k in range(n):
            (<Matrix_modn_dense_template>residues[k])._copy_row_to_mod_int_array(row_list[k],i)
        mm.mpz_crt_vec(M._matrix[i], row_list, nc)

    for k in range(n):
        sage_free(row_list[k])
    sage_free(row_list)

    sig_off()
    return M

#######################################################

# Conclusions:
#  On OS X Intel, at least:
#    - if log_2(height) >= 0.70 * nrows, use classical

def tune_multiplication(k, nmin=10, nmax=200, bitmin=2,bitmax=64):
    """
    Compare various multiplication algorithms.

    INPUT:

    - ``k`` - integer; affects numbers of trials

    - ``nmin`` - integer; smallest matrix to use

    - ``nmax`` - integer; largest matrix to use

    - ``bitmin`` - integer; smallest bitsize

    - ``bitmax`` - integer; largest bitsize

    OUTPUT:

    -  ``prints what doing then who wins`` - multimodular
       or classical

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense import tune_multiplication
        sage: tune_multiplication(2, nmin=10, nmax=60, bitmin=2,bitmax=8)
        10 2 0.2
        ...
    """
    from constructor import random_matrix
    from sage.rings.integer_ring import ZZ
    for n in range(nmin,nmax,10):
        for i in range(bitmin, bitmax, 4):
            A = random_matrix(ZZ, n, n, x = 2**i)
            B = random_matrix(ZZ, n, n, x = 2**i)
            t0 = cputime()
            for j in range(k//n + 1):
                C = A._multiply_classical(B)
            t0 = cputime(t0)
            t1 = cputime()
            for j in range(k//n+1):
                C = A._multiply_multi_modular(B)
            t1 = cputime(t1)
            print n, i, float(i)/float(n)
            if t0 < t1:
                print 'classical'
            else:
                print 'multimod'


##############################################################
# Some people really really really want to make sure their
# 4x4 determinant is really really really fast.
##############################################################

cdef int four_dim_det(mpz_t r,mpz_t *x) except -1:
    """
    Internal function used in computing determinants of 4x4 matrices.

    TESTS::

        sage: A = matrix(ZZ,4,[1,0,3,0,4,3,2,1,0,5,0,0,9,1,2,3])
        sage: A.determinant()     # indirect doctest
        25
    """
    cdef mpz_t a,b

    sig_on()
    mpz_init(a)
    mpz_init(b)

    mpz_mul(a,x[3], x[6] ); mpz_submul(a,x[2], x[7] )
    mpz_mul(b,x[9], x[12]); mpz_submul(b,x[8], x[13])
    mpz_mul(r,a,b)
    mpz_mul(a,x[1], x[7] ); mpz_submul(a,x[3], x[5] )
    mpz_mul(b,x[10],x[12]); mpz_submul(b,x[8], x[14])
    mpz_addmul(r,a,b)
    mpz_mul(a,x[2], x[5] ); mpz_submul(a,x[1], x[6] )
    mpz_mul(b,x[11],x[12]); mpz_submul(b,x[8], x[15])
    mpz_addmul(r,a,b)
    mpz_mul(a,x[3], x[4] ); mpz_submul(a,x[0], x[7] )
    mpz_mul(b,x[10],x[13]); mpz_submul(b,x[9], x[14])
    mpz_addmul(r,a,b)
    mpz_mul(a,x[0], x[6] ); mpz_submul(a,x[2], x[4] )
    mpz_mul(b,x[11],x[13]); mpz_submul(b,x[9], x[15])
    mpz_addmul(r,a,b)
    mpz_mul(a,x[1], x[4] ); mpz_submul(a,x[0], x[5] )
    mpz_mul(b,x[11],x[14]); mpz_submul(b,x[10],x[15])
    mpz_addmul(r,a,b)

    mpz_clear(a)
    mpz_clear(b)
    sig_off()
    return 0

#The above was generated by the following Sage code:
#def idx(x):
#    return [i for i in range(16) if x[i]]
#
#def Det4Maker():
#    Q = PolynomialRing(QQ,['x%s'%ZZ(i).str(16) for i in range(16)])
#    M = Matrix(Q,4,4,Q.gens())
#    D = M.det()
#    Dm = [(m.exponents()[0],D[m]) for m in D.monomials()]
#    Dm.sort(reverse=True)
#    MM = [[Dm[i],Dm[i+1]] for i in range(0,len(Dm),2)]
#
#    S = []
#    for j,((t,s),(r,o)) in enumerate(MM):
#        a,b,c,d = [ZZ(i).str(16) for i in range(16) if t[i]]
#        e,f,g,h = [ZZ(i).str(16) for i in range(16) if r[i]]
#        S.append((c+d+g+h,j))
#    S.sort(lambda x,y: cmp(x[0],y[0]))
#    MM = [MM[s[1]] for s in S]
#    MMM = [[MM[i],MM[i+1]] for i in range(0,len(MM),2)]
#
#    N = 0
#    Cstrux = ""
#    Pstrux = ""
#    for ((t1,s1),(t2,s2)),((t3,s3),(t4,s4)) in MMM:
#        a,b,c,d = idx(t1)
#        e,f = idx(t2)[2:]
#        h,i = idx(t3)[:2]
#
#        if s1 == 1:
#            a,b,h,i = h,i,a,b
#
#        Pstrux+= 'a = x[%s]*x[%s]; '%(a,b)
#        Cstrux+= 'mpz_mul(a,x[%s],x[%s]); '%(a,b)
#
#        Pstrux+= 'a-=x[%s]*x[%s]; '%(h,i)
#        Cstrux+= 'mpz_submul(a,x[%s],x[%s])\n'%(h,i)
#
#        if s4*s1 == 1:
#            c,d,e,f = e,f,c,d
#
#        Pstrux+= 'b = x[%s]*x[%s]; '%(c,d)
#        Cstrux+= 'mpz_mul(b,x[%s],x[%s]); '%(c,d)
#
#        Pstrux+= 'b-=x[%s]*x[%s]; '%(e,f)
#        Cstrux+= 'mpz_submul(b,x[%s],x[%s])\n'%(e,f)
#
#
#        if N:
#            Pstrux+= 'r+= a*b\n'
#            Cstrux+= 'mpz_addmul(r,a,b)\n'
#        else:
#            Pstrux+= 'r = a*b\n'
#            Cstrux+= 'mpz_mul(r,a,b)\n'
#            N = 1
#    return Cstrux,Pstrux
#
#print Det4Maker()[0]

#Note, one can prove correctness of the above by setting
#    sage: Q = PolynomialRing(QQ,['x%s'%ZZ(i).str(16) for i in range(16)])
#    sage: M = Matrix(Q,4,4,Q.gens())
#    sage: D = M.det()
#    sage: x = Q.gens()
#and then evaluating the contents of Det4Maker()[1]...
#    sage: a = x[3]*x[6]; a-=x[2]*x[7]; b = x[9]*x[12]; b-=x[8]*x[13]; r = a*b
#    sage: a = x[1]*x[7]; a-=x[3]*x[5]; b = x[10]*x[12]; b-=x[8]*x[14]; r+= a*b
#    sage: a = x[2]*x[5]; a-=x[1]*x[6]; b = x[11]*x[12]; b-=x[8]*x[15]; r+= a*b
#    sage: a = x[3]*x[4]; a-=x[0]*x[7]; b = x[10]*x[13]; b-=x[9]*x[14]; r+= a*b
#    sage: a = x[0]*x[6]; a-=x[2]*x[4]; b = x[11]*x[13]; b-=x[9]*x[15]; r+= a*b
#    sage: a = x[1]*x[4]; a-=x[0]*x[5]; b = x[11]*x[14]; b-=x[10]*x[15]; r+= a*b
#and comparing results:
#    sage: r-D
#    0
#bingo!

