"""
Dense matrices over GF(2) using the M4RI library.

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

EXAMPLES::

    sage: a = matrix(GF(2),3,range(9),sparse=False); a
    [0 1 0]
    [1 0 1]
    [0 1 0]
    sage: a.rank()
    2
    sage: type(a)
    <type 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>
    sage: a[0,0] = 1
    sage: a.rank()
    3
    sage: parent(a)
    Full MatrixSpace of 3 by 3 dense matrices over Finite Field of size 2

    sage: a^2
    [0 1 1]
    [1 0 0]
    [1 0 1]
    sage: a+a
    [0 0 0]
    [0 0 0]
    [0 0 0]

    sage: b = a.new_matrix(2,3,range(6)); b
    [0 1 0]
    [1 0 1]

    sage: a*b
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 3 by 3 dense matrices over Finite Field of size 2' and 'Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 2'
    sage: b*a
    [1 0 1]
    [1 0 0]

    sage: TestSuite(a).run()
    sage: TestSuite(b).run()

    sage: a.echelonize(); a
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: b.echelonize(); b
    [1 0 1]
    [0 1 0]

TESTS::

    sage: FF = FiniteField(2)
    sage: V = VectorSpace(FF,2)
    sage: v = V([0,1]); v
    (0, 1)
    sage: W = V.subspace([v])
    sage: W
    Vector space of degree 2 and dimension 1 over Finite Field of size 2
    Basis matrix:
    [0 1]
    sage: v in W
    True

    sage: M = Matrix(GF(2), [[1,1,0],[0,1,0]])
    sage: M.row_space()
    Vector space of degree 3 and dimension 2 over Finite Field of size 2
    Basis matrix:
    [1 0 0]
    [0 1 0]

    sage: M = Matrix(GF(2), [[1,1,0],[0,0,1]])
    sage: M.row_space()
    Vector space of degree 3 and dimension 2 over Finite Field of size 2
    Basis matrix:
    [1 1 0]
    [0 0 1]

TODO:

- make LinBox frontend and use it

    - charpoly ?
    - minpoly ?

- make Matrix_modn_frontend and use it (?)
"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2007,2008,2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
include "cysignals/memory.pxi"

cimport matrix_dense
from libc.stdio cimport *
from sage.structure.element cimport (Matrix, Vector, parent_c,
                                     ModuleElement, Element)
from sage.modules.free_module_element cimport FreeModuleElement
from sage.libs.gmp.random cimport *
from sage.misc.functional import log
from sage.misc.randstate cimport randstate, current_randstate
from sage.misc.misc import verbose, get_verbose, cputime
from sage.modules.free_module import VectorSpace
from sage.modules.vector_mod2_dense cimport Vector_mod2_dense

cdef extern from "gd.h":
    ctypedef struct gdImagePtr "gdImagePtr":
        pass

    gdImagePtr gdImageCreateFromPng(FILE *f)
    gdImagePtr gdImageCreateFromPngPtr(int size, void *s)
    gdImagePtr gdImageCreate(int x, int y)
    void gdImagePng(gdImagePtr im, FILE *out)
    void *gdImagePngPtr(gdImagePtr im, int *size)
    void gdImageDestroy(gdImagePtr im)
    int gdImageSX(gdImagePtr im)
    int gdImageSY(gdImagePtr im)
    int gdImageGetPixel(gdImagePtr im, int x, int y)
    void gdImageSetPixel(gdImagePtr im, int x, int y, int value)
    int gdImageColorAllocate(gdImagePtr im, int r, int g, int b)
    void gdImageFilledRectangle(gdImagePtr im, int x1, int y1, int x2, int y2, int color)
    void gdFree(void *m)

## from sage.libs.linbox.linbox cimport Linbox_mod2_dense
## cdef Linbox_mod2_dense linbox
## linbox = Linbox_mod2_dense()

cdef object called

cdef void init_m4ri():
    global called
    if called is None:
        m4ri_build_all_codes()
        called = True

init_m4ri()

def free_m4ri():
    """
    Free global Gray code tables.
    """
    m4ri_destroy_all_codes()



cdef class Matrix_mod2_dense(matrix_dense.Matrix_dense):   # dense or sparse
    """
    Dense matrix over GF(2).
    """
    ########################################################################
    # LEVEL 1 functionality
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce, alloc=True):
        """
        Dense matrix over GF(2) constructor.

        INPUT:

        - ``parent`` - MatrixSpace.
        - ``entries`` - may be list or 0 or 1
        - ``copy`` - ignored, elements are always copied
        - ``coerce`` - ignored, elements are always coerced to ints % 2

        EXAMPLES::

            sage: type(random_matrix(GF(2),2,2))
            <type 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>

            sage: Matrix(GF(2),3,3,1) # indirect doctest
            [1 0 0]
            [0 1 0]
            [0 0 1]

        See :trac:`10858`::

            sage: matrix(GF(2),0,[]) * vector(GF(2),0,[])
            ()
        """
        matrix_dense.Matrix_dense.__init__(self, parent)

        if alloc:
            self._entries = mzd_init(self._nrows, self._ncols)

        # cache elements
        self._zero = self._base_ring(0)
        self._one = self._base_ring(1)

    def __dealloc__(self):
        if self._entries:
            mzd_free(self._entries)
            self._entries = NULL

    def __init__(self, parent, entries, copy, coerce):
        """
        Dense matrix over GF(2) constructor.

        INPUT:

        - ``parent`` - MatrixSpace.
        - ``entries`` - may be list or 0 or 1
        - ``copy`` - ignored, elements are always copied
        - ``coerce`` - ignored, elements are always coerced to ints % 2

        EXAMPLES::

            sage: type(random_matrix(GF(2),2,2))
            <type 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>

            sage: Matrix(GF(2),3,3,1)
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: Matrix(GF(2),2,2,[1,1,1,0])
            [1 1]
            [1 0]

            sage: Matrix(GF(2),2,2,4)
            [0 0]
            [0 0]

            sage: Matrix(GF(2),1,1, 1/3)
            [1]
            sage: Matrix(GF(2),1,1, [1/3])
            [1]

        TESTS::

            sage: Matrix(GF(2),0,0)
            []
            sage: Matrix(GF(2),2,0)
            []
            sage: Matrix(GF(2),0,2)
            []
        """
        cdef int i,j,e

        if entries is None:
            return

        R = self.base_ring()

        # scalar ?
        if not isinstance(entries, list):
            if self._nrows and self._ncols and R(entries) == 1:
                mzd_set_ui(self._entries, 1)
            return

        # all entries are given as a long list
        if len(entries) != self._nrows * self._ncols:
            raise IndexError("The vector of entries has the wrong length.")

        k = 0

        for i from 0 <= i < self._nrows:
            sig_check()
            for j from 0 <= j < self._ncols:
                mzd_write_bit(self._entries,i,j, R(entries[k]))
                k = k + 1

    def __hash__(self):
        r"""
        The has of a matrix is computed as `\oplus i*a_i` where the
        `a_i` are the flattened entries in a matrix (by row, then by
        column).

        EXAMPLE::

            sage: B = random_matrix(GF(2),3,3)
            sage: B.set_immutable()
            sage: {B:0} # indirect doctest
            {[0 1 0]
            [0 1 1]
            [0 0 0]: 0}
            sage: M = random_matrix(GF(2), 123, 321)
            sage: M.set_immutable()
            sage: MZ = M.change_ring(ZZ)
            sage: MZ.set_immutable()
            sage: hash(M) == hash(MZ)
            True
            sage: MS = M.sparse_matrix()
            sage: MS.set_immutable()
            sage: hash(M) == hash(MS)
            True

        TEST::

            sage: A = matrix(GF(2),2,0)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            0

        """
        if not self._is_immutable:
            raise TypeError("mutable matrices are unhashable")

        x = self.fetch('hash')
        if not x is None:
            return x

        if self._nrows == 0 or self._ncols == 0:
            return 0

        cdef unsigned long i, j, truerow
        cdef unsigned long start, shift
        cdef m4ri_word row_xor
        cdef m4ri_word end_mask = __M4RI_LEFT_BITMASK(self._ncols%m4ri_radix)
        cdef m4ri_word top_mask, bot_mask
        cdef m4ri_word cur
        cdef m4ri_word* row

        # running_xor is the xor of all words in the matrix, as if the rows
        # in the matrix were written out consecutively, without regard to
        # word boundaries.
        cdef m4ri_word running_xor = 0
        # running_parity is the number of extra words that must be xor'd.
        cdef unsigned long running_parity = 0


        for i from 0 <= i < self._entries.nrows:

            # All this shifting and masking is because the
            # rows are word-aligned.
            row = self._entries.rows[i]
            start = (i*self._entries.ncols) >> 6
            shift = (i*self._entries.ncols) & 0x3F
            bot_mask = __M4RI_LEFT_BITMASK(m4ri_radix - shift)
            top_mask = ~bot_mask

            if self._entries.width > 1:
                row_xor = row[0]
                running_parity ^= start & parity_mask(row[0] & bot_mask)

                for j from 1 <= j < self._entries.width - 1:
                    row_xor ^= row[j]
                    cur = ((row[j-1] >> (63-shift)) >> 1) ^ (row[j] << shift)
                    running_parity ^= (start+j) & parity_mask(cur)

                running_parity ^= (start+j) & parity_mask(row[j-1] & top_mask)

            else:
                j = 0
                row_xor = 0

            cur = row[j] & end_mask
            row_xor ^= cur
            running_parity ^= (start+j) & parity_mask(cur & bot_mask)
            running_parity ^= (start+j+1) & parity_mask(cur & top_mask)

            running_xor ^= (row_xor << shift) ^ ((row_xor >> (63-shift)) >> 1)

        cdef unsigned long bit_is_set
        cdef unsigned long h

        # Now we assemble the running_parity and running_xor to get the hash.
        # Viewing the flattened matrix as a list of a_i, the hash is the xor
        # of the i for which a_i is non-zero. We split i into the lower m4ri_radix
        # bits and the rest, so i = i1 << m4ri_radix + i0. Now two matching i0
        # would cancel, so we only need the parity of how many of each
        # possible i0 occur. This is stored in the bits of running_xor.
        # Similarly, running_parity is the xor of the i1 needed. It's called
        # parity because i1 is constant across a word, and for each word
        # the number of i1 to add is equal to the number of set bits in that
        # word (but because two cancel, we only need keep track of the
        # parity.

        h = m4ri_radix * running_parity
        for i from 0 <= i < m4ri_radix:
            bit_is_set = (running_xor >> i) & 1
            h ^= (m4ri_radix-1) & ~(bit_is_set-1) & i

        if h == -1:
            h = -2

        self.cache('hash', h)
        return h

    # this exists for compatibility with matrix_modn_dense
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        """
        Set the (i,j) entry of self to the int value.
        """
        mzd_write_bit(self._entries, i, j, int(value))

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        mzd_write_bit(self._entries, i, j, int(value))

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        if mzd_read_bit(self._entries, i, j):
            return self._one
        else:
            return self._zero


    def str(self, rep_mapping=None, zero=None, plus_one=None, minus_one=None):
        r"""
        Return a nice string representation of the matrix.

        INPUT:

        - ``rep_mapping`` - a dictionary or callable used to override
          the usual representation of elements.  For a dictionary,
          keys should be elements of the base ring and values the
          desired string representation.

        - ``zero`` - string (default: ``None``); if not ``None`` use
          the value of ``zero`` as the representation of the zero
          element.

        - ``plus_one`` - string (default: ``None``); if not ``None``
          use the value of ``plus_one`` as the representation of the
          one element.

        - ``minus_one`` - Ignored.  Only for compatibility with
          generic matrices.

        EXAMPLE::

            sage: B = random_matrix(GF(2),3,3)
            sage: B # indirect doctest
            [0 1 0]
            [0 1 1]
            [0 0 0]
            sage: block_matrix([[B, 1], [0, B]])
            [0 1 0|1 0 0]
            [0 1 1|0 1 0]
            [0 0 0|0 0 1]
            [-----+-----]
            [0 0 0|0 1 0]
            [0 0 0|0 1 1]
            [0 0 0|0 0 0]
            sage: B.str(zero='.')
            '[. 1 .]\n[. 1 1]\n[. . .]'
        """
        if self._nrows ==0 or self._ncols == 0:
            return "[]"

        cdef int i,j, last_i
        cdef list s = []
        empty_row = " "*(self._ncols*2-1)
        cdef char *row_s
        cdef char *div_s

        # Set the mapping based on keyword arguments
        # We ignore minus_one (it's only there for compatibility with Matrix)
        if rep_mapping is not None or zero is not None or plus_one is not None:
        # Shunt mappings off to the generic code since they might not be single characters
            return matrix_dense.Matrix_dense.str(self, rep_mapping=rep_mapping, zero=zero, plus_one=plus_one)

        cdef list row_div, col_div
        if self._subdivisions is not None:
            row_s = empty_row
            div_s = row_divider = b"[%s]" % ("-" * (self._ncols*2-1))
            row_div, col_div = self.subdivisions()
            last_i = 0
            for i in col_div:
                if i == last_i or i == self._ncols:
                    # Adjacent column divisions messy, use generic code
                    return matrix_dense.Matrix_dense.str(self, rep_mapping)
                row_s[2*i-1] = '|'
                div_s[2*i] = '+'
                last_i = i

        for i from 0 <= i < self._nrows:
            row_s = row = b"[%s]" % empty_row
            for j from 0 <= j < self._ncols:
                row_s[1+2*j] = c'0' + mzd_read_bit(self._entries,i,j)
            s.append(row)

        if self._subdivisions is not None:
            for i in reversed(row_div):
                s.insert(i, row_divider)

        return "\n".join(s)

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the ``i``'th row of this matrix as a vector.

        This row is a dense vector if and only if the matrix is a dense
        matrix.

        INPUT:

        - ``i`` - integer

        - ``from_list`` - bool (default: ``False``); if ``True``,
          returns the ``i``'th element of ``self.rows()`` (see
          :func:`rows`), which may be faster, but requires building a
          list of all rows the first time it is called after an entry
          of the matrix is changed.

        EXAMPLES::

            sage: A = random_matrix(GF(2),10,10); A
            [0 1 0 1 1 0 0 0 1 1]
            [0 1 1 1 0 1 1 0 0 1]
            [0 0 0 1 0 1 0 0 1 0]
            [0 1 1 0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1 0 1 1]
            [0 0 1 1 1 1 0 0 0 0]
            [1 1 1 1 0 1 0 1 1 0]
            [0 0 0 1 1 0 0 0 1 1]
            [1 0 0 0 1 1 1 0 1 1]
            [1 0 0 1 1 0 1 0 0 0]

            sage: A.row(0)
            (0, 1, 0, 1, 1, 0, 0, 0, 1, 1)

            sage: A.row(-1)
            (1, 0, 0, 1, 1, 0, 1, 0, 0, 0)

            sage: A.row(2,from_list=True)
            (0, 0, 0, 1, 0, 1, 0, 0, 1, 0)

            sage: A = Matrix(GF(2),1,0)
            sage: A.row(0)
            ()
        """
        if self._nrows == 0:
            raise IndexError("matrix has no rows")
        if i >= self._nrows or i < -self._nrows:
            raise IndexError("row index out of range")
        if i < 0:
            i = i + self._nrows
        if from_list:
            return self.rows(copy=False)[i]
        cdef Py_ssize_t j
        cdef Vector_mod2_dense z = Vector_mod2_dense.__new__(Vector_mod2_dense)
        z._init(self._ncols, VectorSpace(self.base_ring(),self._ncols))
        if self._ncols:
            mzd_submatrix(z._entries, self._entries, i, 0, i+1, self._ncols)
        return z

    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _mul_
    #   * cpdef _cmp_
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0

    cpdef _add_(self, right):
        """
        Matrix addition.

        INPUT:

        - right -- matrix of dimension self.nrows() x self.ncols()

        EXAMPLES::

            sage: A = random_matrix(GF(2),10,10)
            sage: A + A == Matrix(GF(2),10,10,0)
            True

            sage: A = random_matrix(GF(2),257,253)
            sage: A + A == Matrix(GF(2),257,253,0) # indirect doctest
            True

        TESTS::

            sage: A = matrix(GF(2),2,0)
            sage: A+A
            []
            sage: A = matrix(GF(2),0,2)
            sage: A+A
            []
            sage: A = matrix(GF(2),0,0)
            sage: A+A
            []
        """
        cdef Matrix_mod2_dense A
        A = Matrix_mod2_dense.__new__(Matrix_mod2_dense, self._parent, 0, 0, 0, alloc=False)
        if self._nrows == 0 or self._ncols == 0:
            return A
        A._entries = mzd_add(NULL, self._entries,(<Matrix_mod2_dense>right)._entries)

        return A

    cpdef _sub_(self, right):
        """
        Matrix addition.

        INPUT:

        - right -- matrix of dimension self.nrows() x self.ncols()

        EXAMPLES::

            sage: A = random_matrix(GF(2),10,10)
            sage: A - A == Matrix(GF(2),10,10,0)  # indirect doctest
            True
        """
        return self._add_(right)

    cdef _matrix_times_vector_(self, Vector v):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),10^4,10^4)
            sage: v0 = random_matrix(GF(2),10^4,1)
            sage: v1 = v0.column(0)
            sage: r0 = A*v0
            sage: r1 = A*v1
            sage: r0.column(0) == r1
            True

        TESTS:

        Check that :trac:`19378` is fixed::

            sage: m = matrix(GF(2), 11, 0)
            sage: v = vector(GF(2), 0)
            sage: m * v
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        cdef mzd_t *tmp
        if not isinstance(v, Vector_mod2_dense):
            M = VectorSpace(self._base_ring, self._nrows)
            v = M(v)
        if self.ncols() != v.degree():
            raise ArithmeticError("number of columns of matrix must equal degree of vector")

        VS = VectorSpace(self._base_ring, self._nrows)
        # If the vector is 0-dimensional, the result will be the 0-vector
        if not self.ncols():
            return VS.zero()
        cdef Vector_mod2_dense c = Vector_mod2_dense.__new__(Vector_mod2_dense)
        c._init(self._nrows, VS)
        c._entries = mzd_init(1, self._nrows)
        if c._entries.nrows and c._entries.ncols:
            tmp = mzd_init(self._nrows, 1)
            _mzd_mul_naive(tmp, self._entries, (<Vector_mod2_dense>v)._entries, 0)
            mzd_transpose(c._entries, tmp)
            mzd_free(tmp)
        return c

    cdef _matrix_times_matrix_(self, Matrix right):
        """
        Matrix multiplication.

        ALGORITHM: Uses the 'Method of the Four Russians
        Multiplication', see :func:`_multiply_m4rm`.
        """
        if get_verbose() >= 2:
            verbose('matrix multiply of %s x %s matrix by %s x %s matrix'%(
                self._nrows, self._ncols, right._nrows, right._ncols))

        return self._multiply_strassen(right, 0)

    cpdef Matrix_mod2_dense _multiply_m4rm(Matrix_mod2_dense self, Matrix_mod2_dense right, int k):
        """
        Multiply matrices using the 'Method of the Four Russians
        Multiplication' (M4RM) or Konrod's method.

        The algorithm is based on an algorithm by Arlazarov, Dinic,
        Kronrod, and Faradzev [ADKF70] and appeared in [AHU]. This
        implementation is based on a description given in Gregory
        Bard's 'Method of the Four Russians Inversion' paper [B06].

        INPUT:

        - right -- Matrix
        - k -- parameter `k` for the Gray Code table size. If `k=0` a suitable
          value is chosen by the function. (`0<= k <= 16`, default: 0)

        EXAMPLE::

              sage: A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              sage: B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              sage: A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              sage: B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              sage: A._multiply_m4rm(B, 0)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]

        TESTS::

            sage: A = random_matrix(GF(2),0,0)
            sage: B = random_matrix(GF(2),0,0)
            sage: A._multiply_m4rm(B, 0)
            []
            sage: A = random_matrix(GF(2),3,0)
            sage: B = random_matrix(GF(2),0,3)
            sage: A._multiply_m4rm(B, 0)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: A = random_matrix(GF(2),0,3)
            sage: B = random_matrix(GF(2),3,0)
            sage: A._multiply_m4rm(B, 0)
            []

        ALGORITHM: Uses the 'Method of the Four Russians'
        multiplication as implemented in the M4RI library.

        REFERENCES:

        ..  [AHU] \A. Aho, J. Hopcroft, and J. Ullman. 'Chapter 6:
                     Matrix Multiplication and Related Operations.'
                     The Design and Analysis of Computer
                     Algorithms. Addison-Wesley, 1974.

        ..  [ADKF70] \V. Arlazarov, E. Dinic, M. Kronrod, and
                     I. Faradzev. 'On Economical Construction of the
                     Transitive Closure of a Directed Graph.'
                     Dokl. Akad. Nauk. SSSR No. 194 (in Russian),
                     English Translation in Soviet Math Dokl. No. 11,
                     1970.

        ..  [Bard06] \G. Bard. 'Accelerating Cryptanalysis with the
                     Method of Four Russians'. Cryptography E-Print
                     Archive (http://eprint.iacr.org/2006/251.pdf),
                     2006.
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        if get_verbose() >= 2:
            verbose('m4rm multiply of %s x %s matrix by %s x %s matrix'%(
                self._nrows, self._ncols, right._nrows, right._ncols))

        cdef Matrix_mod2_dense ans

        ans = self.new_matrix(nrows = self._nrows, ncols = right._ncols)
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans
        sig_on()
        ans._entries = mzd_mul_m4rm(ans._entries, self._entries, right._entries, k)
        sig_off()
        return ans

    def _multiply_classical(Matrix_mod2_dense self, Matrix_mod2_dense right):
        """
        Classical `O(n^3)` multiplication.

        This can be quite fast for matrix vector multiplication but
        the other routines fall back to this implementation in that
        case anyway.

        EXAMPLE::

              sage: A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              sage: B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              sage: A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              sage: B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              sage: A._multiply_classical(B)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]

        TESTS::

            sage: A = random_matrix(GF(2),0,0)
            sage: B = random_matrix(GF(2),0,0)
            sage: A._multiply_classical(B)
            []
            sage: A = random_matrix(GF(2),3,0)
            sage: B = random_matrix(GF(2),0,3)
            sage: A._multiply_classical(B)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: A = random_matrix(GF(2),0,3)
            sage: B = random_matrix(GF(2),3,0)
            sage: A._multiply_classical(B)
            []
        """
        cdef Matrix_mod2_dense A
        A = self.new_matrix(nrows = self._nrows, ncols = right._ncols)
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return A
        A._entries = mzd_mul_naive(A._entries, self._entries,(<Matrix_mod2_dense>right)._entries)
        return A

    cpdef Matrix_mod2_dense _multiply_strassen(Matrix_mod2_dense self, Matrix_mod2_dense right, int cutoff):
        r"""
        Strassen-Winograd `O(n^{2.807})` multiplication [Str69].

        This implementation in M4RI is inspired by Sage's generic
        Strassen implementation [BHS08] but uses a more memory
        efficient operation schedule [DP08].

        The performance of this routine depends on the parameter
        cutoff. On many modern machines 2048 should give acceptable
        performance, a good rule of thumb for calculating the optimal
        cutoff would that two matrices of the cutoff size should fit
        in L2 cache, so: `cutoff = \sqrt{L2 * 8 * 1024^2 / 2}` where
        `L2` is the size of the L2 cache in MB.

        INPUT:

        - ``right`` - a matrix of matching dimensions.
        - ``cutoff`` - matrix dimension where M4RM should be used
          instead of Strassen (default: let M4RI decide)

        EXAMPLE::

              sage: A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              sage: B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              sage: A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              sage: B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              sage: A._multiply_strassen(B, 0)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]
              sage: A = random_matrix(GF(2),2701,3000)
              sage: B = random_matrix(GF(2),3000,3172)
              sage: A._multiply_strassen(B, 256) == A._multiply_m4rm(B, 0)
              True

        TESTS::

            sage: A = random_matrix(GF(2),0,0)
            sage: B = random_matrix(GF(2),0,0)
            sage: A._multiply_strassen(B, 0)
            []
            sage: A = random_matrix(GF(2),3,0)
            sage: B = random_matrix(GF(2),0,3)
            sage: A._multiply_strassen(B, 0)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: A = random_matrix(GF(2),0,3)
            sage: B = random_matrix(GF(2),3,0)
            sage: A._multiply_strassen(B, 0)
            []

        ALGORITHM: Uses Strassen-Winograd matrix multiplication with
        M4RM as base case as implemented in the M4RI library.

        REFERENCES:

        ..  [Str69] Volker Strassen. Gaussian elimination is not
                    optimal. Numerische Mathematik, 13:354-356, 1969.

        ..  [BHS08] Robert Bradshaw, David Harvey and William
                    Stein. strassen_window_multiply_c. strassen.pyx,
                    Sage 3.0, 2008. http://www.sagemath.org

        ..  [DP08] Jean-Guillaume Dumas and Clement Pernet. Memory
                   efficient scheduling of Strassen-Winograd's matrix
                   multiplication algorithm. arXiv:0707.2347v1, 2008.
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_mod2_dense ans
        #ans = self.new_matrix(nrows = self._nrows, ncols = right._ncols)
        # The following is a little faster:
        ans = self.matrix_space(self._nrows, right._ncols, sparse=False).zero_matrix().__copy__()
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans

        sig_on()
        ans._entries = mzd_mul(ans._entries, self._entries, right._entries, cutoff)
        sig_off()
        return ans

    def __neg__(self):
        """
        EXAMPLES::

            sage: A = random_matrix(GF(2),100,100)
            sage: A - A == A - -A
            True
        """
        return self.__copy__()

    def __invert__(self):
        """
        Inverts self using the 'Method of the Four Russians'
        inversion.

        If ``self`` is not invertible a ``ZeroDivisionError`` is
        raised.

        EXAMPLE::

            sage: A = Matrix(GF(2),3,3, [0, 0, 1, 0, 1, 1, 1, 0, 1])
            sage: MS = A.parent()
            sage: A
            [0 0 1]
            [0 1 1]
            [1 0 1]
            sage: ~A
            [1 0 1]
            [1 1 0]
            [1 0 0]
            sage: A * ~A == ~A * A == MS(1)
            True

        TESTS::

            sage: A = matrix(GF(2),0,0)
            sage: A^(-1)
            []
        """
        cdef int k = 0
        cdef mzd_t *I
        cdef Matrix_mod2_dense A

        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        if self._ncols == 0:
            return self.__copy__()

        if self.rank() != self._nrows:
            raise ZeroDivisionError("Matrix does not have full rank.")

        A = Matrix_mod2_dense.__new__(Matrix_mod2_dense, self._parent, 0, 0, 0, alloc = False)
        sig_on()
        A._entries = mzd_inv_m4ri(NULL, self._entries, 0)
        sig_off()

        if A._entries==NULL:
            raise ZeroDivisionError("input matrix must be nonsingular")
        else:
            return A

    def __copy__(self):
        """
        Returns a copy of ``self``.

        EXAMPLES::

             sage: MS = MatrixSpace(GF(2),3,3)
             sage: A = MS(1)
             sage: A.__copy__() == A
             True
             sage: A.__copy__() is A
             False

             sage: A = random_matrix(GF(2),100,100)
             sage: A.__copy__() == A
             True
             sage: A.__copy__() is A
             False

             sage: A.echelonize()
             sage: A.__copy__() == A
             True

        """
        cdef Matrix_mod2_dense A
        A = Matrix_mod2_dense.__new__(Matrix_mod2_dense, self._parent, 0, 0, 0)

        if self._nrows and self._ncols:
            mzd_copy(A._entries, self._entries)

        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())

        return A

    def _list(self):
        """
        Returns list of the elements of ``self`` in row major
        ordering.

        EXAMPLE::

            sage: A = Matrix(GF(2),2,2,[1,0,1,1])
            sage: A
            [1 0]
            [1 1]
            sage: A.list() #indirect doctest
            [1, 0, 1, 1]

        TESTS::

            sage: A = Matrix(GF(2),3,0)
            sage: A.list()
            []
        """
        cdef int i,j
        l = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if mzd_read_bit(self._entries,i,j):
                    l.append(self._one)
                else:
                    l.append(self._zero)
        return l

    # def _dict(self):

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * __deepcopy__
    #    * Matrix windows -- only if you need strassen for that base
    ########################################################################

    def echelonize(self, algorithm='heuristic', cutoff=0, reduced=True, **kwds):
        """
        Puts self in (reduced) row echelon form.

        INPUT:

        - self -- a mutable matrix
        - algorithm

            - 'heuristic' -- uses M4RI and PLUQ (default)
            - 'm4ri' -- uses M4RI
            - 'pluq' -- uses PLUQ factorization
            - 'classical' -- uses classical Gaussian elimination

        - k --  the parameter 'k' of the M4RI algorithm. It MUST be between 1
          and 16 (inclusive). If it is not specified it will be calculated as
          3/4 * log_2( min(nrows, ncols) ) as suggested in the M4RI paper.
        - reduced -- return reduced row echelon form (default:True)

        EXAMPLE::

             sage: A = random_matrix(GF(2), 10, 10)
             sage: B = A.__copy__(); B.echelonize() # fastest
             sage: C = A.__copy__(); C.echelonize(k=2) # force k
             sage: E = A.__copy__(); E.echelonize(algorithm='classical') # force Gaussian elimination
             sage: B == C == E
             True

        TESTS::

             sage: VF2 = VectorSpace(GF(2),2)
             sage: WF2 = VF2.submodule([VF2([1,1])])
             sage: WF2
             Vector space of degree 2 and dimension 1 over Finite Field of size 2
             Basis matrix:
             [1 1]

             sage: A2 = matrix(GF(2),2,[1,0,0,1])
             sage: A2.kernel()
             Vector space of degree 2 and dimension 0 over Finite Field of size 2
             Basis matrix:
             []

        ALGORITHM:

        Uses M4RI library

        REFERENCES:

        .. [Bard06] \G. Bard. 'Accelerating Cryptanalysis with the Method of
           Four Russians'. Cryptography E-Print Archive
           (http://eprint.iacr.org/2006/251.pdf), 2006.
        """
        if self._nrows == 0 or self._ncols == 0:
            self.cache('in_echelon_form',True)
            self.cache('rank', 0)
            self.cache('pivots', ())
            return self
        cdef int k, n, full

        full = int(reduced)

        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form

        if algorithm == 'heuristic':

            self.check_mutability()
            self.clear_cache()

            sig_on()
            r =  mzd_echelonize(self._entries, full)
            sig_off()

            self.cache('in_echelon_form',True)
            self.cache('rank', r)
            self.cache('pivots', tuple(self._pivots()))

        elif algorithm == 'm4ri':

            self.check_mutability()
            self.clear_cache()

            if 'k' in kwds:
                k = int(kwds['k'])

                if k<1 or k>16:
                    raise RuntimeError("k must be between 1 and 16")
                k = round(k)
            else:
                k = 0

            sig_on()
            r =  mzd_echelonize_m4ri(self._entries, full, k)
            sig_off()

            self.cache('in_echelon_form',True)
            self.cache('rank', r)
            self.cache('pivots', tuple(self._pivots()))


        elif algorithm == 'pluq':

            self.check_mutability()
            self.clear_cache()

            sig_on()
            r =  mzd_echelonize_pluq(self._entries, full)
            sig_off()

            self.cache('in_echelon_form',True)
            self.cache('rank', r)
            self.cache('pivots', tuple(self._pivots()))

        elif algorithm == 'linbox':

            #self._echelonize_linbox()
            raise NotImplementedError

        elif algorithm == 'classical':

            # for debugging purposes only, it is slow
            self._echelon_in_place_classical()
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

    def _pivots(self):
        """
        Returns the pivot columns of ``self`` if ``self`` is in
        row echelon form.

        EXAMPLE::

            sage: A = matrix(GF(2),5,5,[0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1])
            sage: E = A.echelon_form()
            sage: E
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            [0 0 0 0 0]
            sage: E._pivots()
            [1, 2, 3, 4]
        """
        if not self.fetch('in_echelon_form'):
            raise RuntimeError("self must be in reduced row echelon form first.")
        pivots = []
        cdef Py_ssize_t i, j, nc
        nc = self._ncols
        i = 0
        while i < self._nrows:
            for j from i <= j < nc:
                if mzd_read_bit(self._entries, i, j):
                    pivots.append(j)
                    i += 1
                    break
            if j == nc:
                break
        return pivots

    def randomize(self, density=1, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        INPUT:

        -  ``density`` - float; proportion (roughly) to be considered for
           changes
        -  ``nonzero`` - Bool (default: ``False``); whether the new entries
           are forced to be non-zero

        OUTPUT:

        -  None, the matrix is modified in-space

        EXAMPLES::

            sage: A = matrix(GF(2), 5, 5, 0)
            sage: A.randomize(0.5); A
            [0 0 0 1 1]
            [0 1 0 0 1]
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 0 1 0]
            sage: A.randomize(); A
            [0 0 1 1 0]
            [1 1 0 0 1]
            [1 1 1 1 0]
            [1 1 1 1 1]
            [0 0 1 1 0]

        TESTS:

        With the libc random number generator random(), we had problems
        where the ranks of all of these matrices would be the same
        (and they would all be far too low).  This verifies that the
        problem is gone, with Mersenne Twister::

            sage: MS2 = MatrixSpace(GF(2), 1000)
            sage: [MS2.random_element().rank() for i in range(5)]
            [999, 998, 1000, 999, 999]

        Testing corner case::

            sage: A = random_matrix(GF(2),3,0)
            sage: A
            []
        """
        if self._ncols == 0 or self._nrows == 0:
            return

        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.check_mutability()
        self.clear_cache()

        cdef randstate rstate = current_randstate()

        cdef int i, j, k
        cdef int nc
        cdef unsigned int low, high
        cdef m4ri_word mask = 0

        # Original code, before adding the ``nonzero`` option.
        if not nonzero:
            if density == 1:
                assert(sizeof(m4ri_word) == 8)
                mask = __M4RI_LEFT_BITMASK(self._entries.ncols % m4ri_radix)
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._entries.width:
                        # for portability we get 32-bit twice rather than 64-bit once
                        low = gmp_urandomb_ui(rstate.gmp_state, 32)
                        high = gmp_urandomb_ui(rstate.gmp_state, 32)
                        self._entries.rows[i][j] = m4ri_swap_bits( ((<unsigned long long>high)<<32) | (<unsigned long long>low) )
                    self._entries.rows[i][self._entries.width - 1] &= mask
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                sig_on()
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random()%nc
                        mzd_write_bit(self._entries, i, k, rstate.c_random() % 2)
                sig_off()

        # New code for the case when ``nonzero`` is ``True``.
        else:
            sig_on()
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    if rstate.c_rand_double() <= density:
                        mzd_write_bit(self._entries, i, j, 1)
            sig_off()

    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        """
        EXAMPLE::

            sage: A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]

            sage: A.rescale_row(0,0,0); A
            [0 0 0]
            [0 1 1]
            [0 0 0]
        """
        if (int(multiple)%2) == 0:
            mzd_row_clear_offset(self._entries, row, start_col);

    cdef add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, multiple,
                               Py_ssize_t start_col):
        """
        EXAMPLE::

            sage: A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            sage: A.add_multiple_of_row(0,1,1,0); A
            [0 0 1]
            [0 1 1]
            [0 0 0]
        """
        if (int(multiple)%2) != 0:
            mzd_row_add_offset(self._entries, row_to, row_from, start_col)

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        """
        EXAMPLE::

            sage: A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            sage: A.swap_rows(0,1); A
            [0 1 1]
            [0 1 0]
            [0 0 0]
        """
        mzd_row_swap(self._entries, row1, row2)

    cdef swap_columns_c(self, Py_ssize_t col1, Py_ssize_t col2):
        """
        EXAMPLE::

            sage: A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]

            sage: A.swap_columns(0,1); A
            [1 0 0]
            [1 0 1]
            [0 0 0]

            sage: A = random_matrix(GF(2),3,65)

            sage: B = A.__copy__()
            sage: B.swap_columns(0,1)
            sage: B.swap_columns(0,1)
            sage: A == B
            True

            sage: A.swap_columns(0,64)
            sage: A.column(0) == B.column(64)
            True
            sage: A.swap_columns(63,64)
            sage: A.column(63) == B.column(0)
            True
        """
        mzd_col_swap(self._entries, col1, col2)



    def _magma_init_(self, magma):
        """
        Returns a string of self in ``Magma`` form. Does not return
        ``Magma`` object but string.

        EXAMPLE::

            sage: A = random_matrix(GF(2),3,3)
            sage: A._magma_init_(magma)                             # optional - magma
            'Matrix(GF(2),3,3,StringToIntegerSequence("0 1 0 0 1 1 0 0 0"))'
            sage: A = random_matrix(GF(2),100,100)
            sage: B = random_matrix(GF(2),100,100)
            sage: magma(A*B) == magma(A) * magma(B)                 # optional - magma
            True

        TESTS::

            sage: A = random_matrix(GF(2),0,3)
            sage: magma(A)                          # optional - magma
            Matrix with 0 rows and 3 columns
            sage: A = matrix(GF(2),2,3,[0,1,1,1,0,0])
            sage: A._magma_init_(magma)             # optional - magma
            'Matrix(GF(2),2,3,StringToIntegerSequence("0 1 1 1 0 0"))'
            sage: magma(A)                          # optional - magma
            [0 1 1]
            [1 0 0]
        """
        s = self.base_ring()._magma_init_(magma)
        return 'Matrix(%s,%s,%s,StringToIntegerSequence("%s"))'%(
            s, self._nrows, self._ncols, self._export_as_string())

    def determinant(self):
        """
        Return the determinant of this matrix over GF(2).

        EXAMPLES::

            sage: matrix(GF(2),2,[1,1,0,1]).determinant()
            1
            sage: matrix(GF(2),2,[1,1,1,1]).determinant()
            0
        """
        if not self.is_square():
            raise ValueError("self must be a square matrix")
        return self.base_ring()(1 if self.rank() == self.nrows() else 0)

    def transpose(self):
        """
        Returns transpose of self and leaves self untouched.

        EXAMPLE::

            sage: A = Matrix(GF(2),3,5,[1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0])
            sage: A
            [1 0 1 0 0]
            [0 1 1 0 0]
            [1 1 0 1 0]
            sage: B = A.transpose(); B
            [1 0 1]
            [0 1 1]
            [1 1 0]
            [0 0 1]
            [0 0 0]
            sage: B.transpose() == A
            True

        ``.T`` is a convenient shortcut for the transpose::

            sage: A.T
            [1 0 1]
            [0 1 1]
            [1 1 0]
            [0 0 1]
            [0 0 0]

        TESTS::

            sage: A = random_matrix(GF(2),0,40)
            sage: A.transpose()
            40 x 0 dense matrix over Finite Field of size 2 (use the '.str()' method to see the entries)

            sage: A = Matrix(GF(2), [1,0])
            sage: B = A.transpose()
            sage: A[0,0] = 0
            sage: B[0,0]
            1
        """
        cdef Matrix_mod2_dense A = self.new_matrix(ncols = self._nrows,  nrows = self._ncols)
        if self._nrows == 0 or self._ncols == 0:
            return A

        A._entries = mzd_transpose(A._entries, self._entries)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    cpdef int _cmp_(self, right) except -2:
        """
        Compares ``self`` with ``right``. While equality and
        inequality are clearly defined, ``<`` and ``>`` are not.  For
        those first the matrix dimensions of ``self`` and ``right``
        are compared. If these match then ``<`` means that there is a
        position smallest (i,j) in ``self`` where ``self[i,j]`` is
        zero but ``right[i,j]`` is one. This (i,j) is smaller than the
        (i,j) if ``self`` and ``right`` are exchanged for the
        comparison.

        EXAMPLE::

            sage: A = MatrixSpace(GF(2),3,3).one()
            sage: B = copy(MatrixSpace(GF(2),3,3).one())
            sage: B[0,1] = 1
            sage: A < B
            True

        TESTS::

            sage: A = matrix(GF(2),2,0)
            sage: B = matrix(GF(2),2,0)
            sage: A < B
            False
        """
        if self._nrows == 0 or self._ncols == 0:
            return 0
        return mzd_cmp(self._entries, (<Matrix_mod2_dense>right)._entries)


    def augment(self, right, subdivide=False):
        r"""
        Augments ``self`` with ``right``.

        EXAMPLE::

            sage: MS = MatrixSpace(GF(2),3,3)
            sage: A = MS([0, 1, 0, 1, 1, 0, 1, 1, 1]); A
            [0 1 0]
            [1 1 0]
            [1 1 1]
            sage: B = A.augment(MS(1)); B
            [0 1 0 1 0 0]
            [1 1 0 0 1 0]
            [1 1 1 0 0 1]
            sage: B.echelonize(); B
            [1 0 0 1 1 0]
            [0 1 0 1 0 0]
            [0 0 1 0 1 1]
            sage: C = B.matrix_from_columns([3,4,5]); C
            [1 1 0]
            [1 0 0]
            [0 1 1]
            sage: C == ~A
            True
            sage: C*A == MS(1)
            True

        A vector may be augmented to a matrix. ::

            sage: A = matrix(GF(2), 3, 4, range(12))
            sage: v = vector(GF(2), 3, range(3))
            sage: A.augment(v)
            [0 1 0 1 0]
            [0 1 0 1 1]
            [0 1 0 1 0]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``right``.  For more details about how subdivisions
        are managed when augmenting, see
        :meth:`sage.matrix.matrix1.Matrix.augment`.  ::

            sage: A = matrix(GF(2), 3, 5, range(15))
            sage: B = matrix(GF(2), 3, 3, range(9))
            sage: A.augment(B, subdivide=True)
            [0 1 0 1 0|0 1 0]
            [1 0 1 0 1|1 0 1]
            [0 1 0 1 0|0 1 0]

        TESTS::

            sage: A = random_matrix(GF(2),2,3)
            sage: B = random_matrix(GF(2),2,0)
            sage: A.augment(B)
            [0 1 0]
            [0 1 1]

            sage: B.augment(A)
            [0 1 0]
            [0 1 1]

            sage: M = Matrix(GF(2), 0, 0, 0)
            sage: N = Matrix(GF(2), 0, 19, 0)
            sage: W = M.augment(N)
            sage: W.ncols()
            19
            sage: M = Matrix(GF(2), 0, 1, 0)
            sage: N = Matrix(GF(2), 0, 1, 0)
            sage: M.augment(N)
            []

        Check that :trac:`19165` is solved::

            sage: m = matrix(GF(2), 2, range(4))
            sage: m.augment(matrix(GF(2), 2, range(4), sparse=True))
            [0 1 0 1]
            [0 1 0 1]

            sage: m.augment(1)
            Traceback (most recent call last):
            ...
            TypeError: right must either be a matrix or a vector. Not
            <type 'sage.rings.integer.Integer'>
        """
        cdef Matrix_mod2_dense other

        if isinstance(right, FreeModuleElement):
            right = right.column()

        if isinstance(right, Matrix_mod2_dense):
            other = <Matrix_mod2_dense> right
        elif isinstance(right, Matrix):
            from sage.matrix.constructor import matrix
            other = matrix(self.base_ring(),
                           right.nrows(),
                           right.ncols(),
                           right.list(),
                           sparse=False)
        else:
            raise TypeError("right must either be a matrix or a vector. Not {}".format(type(right)))

        if self._nrows != other._nrows:
            raise TypeError("Both numbers of rows must match.")

        if self._ncols == 0:
            return other.__copy__()
        if other._ncols == 0:
            return self.__copy__()

        cdef Matrix_mod2_dense Z
        Z = self.new_matrix(ncols = self._ncols + other._ncols)
        if self._nrows == 0:
            return Z
        Z._entries = mzd_concat(Z._entries, self._entries, other._entries)
        if subdivide:
            Z._subdivide_on_augment(self, other)
        return Z

    cdef _stack_impl(self, bottom):
        r"""
        Stack ``self`` on top of ``bottom``.

        EXAMPLE::

            sage: A = matrix(GF(2),2,2,[1,0,0,1])
            sage: B = matrix(GF(2),2,2,[0,1,1,0])
            sage: A.stack(B)
            [1 0]
            [0 1]
            [0 1]
            [1 0]
            sage: B.stack(A)
            [0 1]
            [1 0]
            [1 0]
            [0 1]

        A vector may be stacked below a matrix. ::

            sage: A = matrix(GF(2), 2, 5, range(10))
            sage: v = vector(GF(2), 5, range(5))
            sage: A.stack(v)
            [0 1 0 1 0]
            [1 0 1 0 1]
            [0 1 0 1 0]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``bottom``.  For more details about how subdivisions
        are managed when stacking, see
        :meth:`sage.matrix.matrix1.Matrix.stack`.  ::

            sage: A = matrix(GF(2), 3, 5, range(15))
            sage: B = matrix(GF(2), 1, 5, range(5))
            sage: A.stack(B, subdivide=True)
            [0 1 0 1 0]
            [1 0 1 0 1]
            [0 1 0 1 0]
            [---------]
            [0 1 0 1 0]

        TESTS::

            sage: A = random_matrix(GF(2),0,3)
            sage: B = random_matrix(GF(2),3,3)
            sage: A.stack(B)
            [0 1 0]
            [0 1 1]
            [0 0 0]

            sage: B.stack(A)
            [0 1 0]
            [0 1 1]
            [0 0 0]

            sage: M = Matrix(GF(2), 0, 0, 0)
            sage: N = Matrix(GF(2), 19, 0, 0)
            sage: W = M.stack(N)
            sage: W.nrows()
            19
            sage: M = Matrix(GF(2), 1, 0, 0)
            sage: N = Matrix(GF(2), 1, 0, 0)
            sage: M.stack(N)
            []
        """
        cdef Matrix_mod2_dense other = <Matrix_mod2_dense>bottom
        cdef Matrix_mod2_dense Z
        Z = self.new_matrix(nrows=self._nrows + other._nrows, ncols=self._ncols)
        if self._ncols > 0:
            Z._entries = mzd_stack(Z._entries, self._entries, other._entries)
        return Z

    def submatrix(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return submatrix from the index row, col (inclusive) with
        dimension nrows x ncols.

        INPUT:

        - row -- index of start row
        - col -- index of start column
        - nrows -- number of rows of submatrix
        - ncols -- number of columns of submatrix

        EXAMPLES::

             sage: A = random_matrix(GF(2),200,200)
             sage: A[0:2,0:2] == A.submatrix(0,0,2,2)
             True
             sage: A[0:100,0:100] == A.submatrix(0,0,100,100)
             True
             sage: A == A.submatrix(0,0,200,200)
             True

             sage: A[1:3,1:3] == A.submatrix(1,1,2,2)
             True
             sage: A[1:100,1:100] == A.submatrix(1,1,99,99)
             True
             sage: A[1:200,1:200] == A.submatrix(1,1,199,199)
             True

        TESTS for handling of default arguments (:trac:`18761`)::

             sage: A.submatrix(17,15) == A.submatrix(17,15,183,185)
             True
             sage: A.submatrix(row=100,col=37,nrows=1,ncols=3) == A.submatrix(100,37,1,3)
             True
        """
        cdef Matrix_mod2_dense A

        cdef int highr, highc

        if nrows < 0:
            nrows = self._nrows - row
            if nrows < 0:
                nrows = 0

        if ncols < 0:
            ncols = self._ncols - col
            if ncols < 0:
                ncols = 0

        highr = row + nrows
        highc = col + ncols

        if row < 0:
            raise TypeError("Expected row >= 0, but got %d instead."%row)

        if col < 0:
            raise TypeError("Expected col >= 0, but got %d instead."%col)

        if highc > self._entries.ncols:
            raise TypeError("Expected highc <= self.ncols(), but got %d > %d instead."%(highc, self._entries.ncols))

        if highr > self._entries.nrows:
            raise TypeError("Expected highr <= self.nrows(), but got %d > %d instead."%(highr, self._entries.nrows))

        A = self.new_matrix(nrows = nrows, ncols = ncols)
        if ncols == 0 or nrows == 0:
            return A
        A._entries = mzd_submatrix(A._entries, self._entries, row, col, highr, highc)
        return A

    def __reduce__(self):
        """
        Serialize ``self``.

        EXAMPLE::

            sage: A = random_matrix(GF(2),10,10)
            sage: f,s = A.__reduce__()
            sage: f(*s) == A
            True
        """
        cdef int i,j, r,c, size

        r, c = self.nrows(), self.ncols()
        if r == 0 or c == 0:
            return unpickle_matrix_mod2_dense_v1, (r, c, None, 0)

        sig_on()
        cdef gdImagePtr im = gdImageCreate(c, r)
        sig_off()
        cdef int black = gdImageColorAllocate(im, 0, 0, 0)
        cdef int white = gdImageColorAllocate(im, 255, 255, 255)
        gdImageFilledRectangle(im, 0, 0, c-1, r-1, white)
        for i from 0 <= i < r:
            for j from 0 <= j < c:
                if mzd_read_bit(self._entries, i, j):
                    gdImageSetPixel(im, j, i, black )

        cdef signed char *buf = <signed char*>gdImagePngPtr(im, &size)

        data = [buf[i] for i in range(size)]
        gdFree(buf)
        gdImageDestroy(im)
        return unpickle_matrix_mod2_dense_v1, (r,c, data, size)

    cpdef _export_as_string(self):
        """
        Return space separated string of the entries in this matrix.

        EXAMPLES::

            sage: w = matrix(GF(2),2,3,[1,0,1,1,1,0])
            sage: w._export_as_string()
            '1 0 1 1 1 0'
        """
        cdef Py_ssize_t i, j, k, n
        cdef char *s
        cdef char *t

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*2 + 2
            s = <char*> sig_malloc(n * sizeof(char))
            k = 0
            sig_on()
            for i in range(self._nrows):
                for j in range(self._ncols):
                    s[k] = <char>(48 + (1 if mzd_read_bit(self._entries,i,j) else 0))  # "0" or "1"
                    k += 1
                    s[k] = <char>32  # space
                    k += 1
            sig_off()
            s[k-1] = <char>0
            data = str(s)
            sig_free(s)
        return data

    def density(self, approx=False):
        """
        Return the density of this matrix.

        By density we understand the ration of the number of nonzero
        positions and the self.nrows() * self.ncols(), i.e. the number
        of possible nonzero positions.

        INPUT:

        - approx -- return floating point approximation (default: False)

        EXAMPLE::

            sage: A = random_matrix(GF(2),1000,1000)
            sage: d = A.density(); d
            62483/125000

            sage: float(d)
            0.499864

            sage: A.density(approx=True)
            0.499864000...

            sage: float(len(A.nonzero_positions())/1000^2)
            0.499864
        """
        if approx:
            from sage.rings.real_mpfr import create_RealNumber
            return create_RealNumber(mzd_density(self._entries, 1))
        else:
            return matrix_dense.Matrix_dense.density(self)

    def rank(self, algorithm='ple'):
        """
        Return the rank of this matrix.

        On average 'ple' should be faster than 'm4ri' and hence it is
        the default choice. However, for small - i.e. quite few
        thousand rows & columns - and sparse matrices 'm4ri' might be
        a better choice.

        INPUT:

        - ``algorithm`` - either "ple" or "m4ri"

        EXAMPLE::

            sage: A = random_matrix(GF(2), 1000, 1000)
            sage: A.rank()
            999

            sage: A = matrix(GF(2),10, 0)
            sage: A.rank()
            0
        """
        x = self.fetch('rank')
        if not x is None:
            return x
        if self._nrows == 0 or self._ncols == 0:
            return 0
        cdef mzd_t *A = mzd_copy(NULL, self._entries)
        cdef mzp_t *P
        cdef mzp_t *Q

        if algorithm == 'ple':
            P = mzp_init(self._entries.nrows)
            Q = mzp_init(self._entries.ncols)
            r = mzd_ple(A, P, Q, 0)
            mzp_free(P)
            mzp_free(Q)
        elif algorithm == 'm4ri':
            r = mzd_echelonize_m4ri(A, 0, 0)
        else:
            raise ValueError("Algorithm '%s' unknown."%algorithm)
        mzd_free(A)
        self.cache('rank', r)
        return r

    def _right_kernel_matrix(self, **kwds):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``kwds`` - these are provided for consistency with other versions
          of this method.  Here they are ignored as there is no optional
          behavior available.

        OUTPUT:

        Returns a pair.  First item is the string 'computed-pluq'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integers mod 2, as computed by the M4RI library
        using PLUQ matrix decomposition.

        EXAMPLES::

            sage: A = matrix(GF(2), [
            ...                      [0, 1, 0, 0, 1, 0, 1, 1],
            ...                      [1, 0, 1, 0, 0, 1, 1, 0],
            ...                      [0, 0, 1, 0, 0, 1, 0, 1],
            ...                      [1, 0, 1, 1, 0, 1, 1, 0],
            ...                      [0, 0, 1, 0, 0, 1, 0, 1],
            ...                      [1, 1, 0, 1, 1, 0, 0, 0]])
            sage: A
            [0 1 0 0 1 0 1 1]
            [1 0 1 0 0 1 1 0]
            [0 0 1 0 0 1 0 1]
            [1 0 1 1 0 1 1 0]
            [0 0 1 0 0 1 0 1]
            [1 1 0 1 1 0 0 0]
            sage: result = A._right_kernel_matrix()
            sage: result[0]
            'computed-pluq'
            sage: result[1]
            [0 1 0 0 1 0 0 0]
            [0 0 1 0 0 1 0 0]
            [1 1 0 0 0 0 1 0]
            [1 1 1 0 0 0 0 1]
            sage: X = result[1].transpose()
            sage: A*X == zero_matrix(GF(2), 6, 4)
            True

        TESTS:

        We test the three trivial cases.  Matrices with no rows or columns will
        cause segfaults in the M4RI code, so we protect against that instance.
        Computing a kernel or a right kernel matrix should never pass these
        problem matrices here. ::

            sage: A = matrix(GF(2), 0, 2)
            sage: A._right_kernel_matrix()
            Traceback (most recent call last):
            ...
            ValueError: kernels of matrices mod 2 with zero rows or zero columns cannot be computed
            sage: A = matrix(GF(2), 2, 0)
            sage: A._right_kernel_matrix()
            Traceback (most recent call last):
            ...
            ValueError: kernels of matrices mod 2 with zero rows or zero columns cannot be computed
            sage: A = zero_matrix(GF(2), 4, 3)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over integers mod 2 for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        if self.nrows()==0 or self.ncols()==0:
            raise ValueError("kernels of matrices mod 2 with zero rows or zero columns cannot be computed")
        cdef Matrix_mod2_dense M
        cdef mzd_t *A = mzd_copy(NULL, self._entries)
        # Despite the name, this next call returns X such that M*X = 0
        cdef mzd_t *k = mzd_kernel_left_pluq(A, 0)
        mzd_free(A)
        if k != NULL:
            M = self.new_matrix(nrows = k.ncols, ncols = k.nrows)
            mzd_transpose(M._entries, k)
            mzd_free(k)
        else:
            M = self.new_matrix(nrows = 0, ncols = self._ncols)
        verbose("done computing right kernel matrix over integers mod 2 for %sx%s matrix" % (self.nrows(), self.ncols()),level=1, t=tm)
        return 'computed-pluq', M

# Used for hashing
cdef int i, k
cdef unsigned long parity_table[256]
for i from 0 <= i < 256:
    parity_table[i] = 1 & ((i) ^ (i>>1) ^ (i>>2) ^ (i>>3) ^
                           (i>>4) ^ (i>>5) ^ (i>>6) ^ (i>>7))

# gmp's ULONG_PARITY may use special
# assembly instructions, could be faster
cpdef inline unsigned long parity(m4ri_word a):
    """
    Returns the parity of the number of bits in a.

    EXAMPLES::

        sage: from sage.matrix.matrix_mod2_dense import parity
        sage: parity(1)
        1L
        sage: parity(3)
        0L
        sage: parity(0x10000101011)
        1L
    """
    if sizeof(m4ri_word) == 8:
        a ^= a >> 32
    a ^= a >> 16
    a ^= a >> 8
    return parity_table[a & 0xFF]

cdef inline unsigned long parity_mask(m4ri_word a):
    return -parity(a)


def unpickle_matrix_mod2_dense_v1(r, c, data, size):
    """
    Deserialize a matrix encoded in the string ``s``.

    INPUT:

    - r -- number of rows of matrix
    - c -- number of columns of matrix
    - s -- a string
    - size -- length of the string s

    EXAMPLE::

        sage: A = random_matrix(GF(2),100,101)
        sage: _,(r,c,s,s2) = A.__reduce__()
        sage: from sage.matrix.matrix_mod2_dense import unpickle_matrix_mod2_dense_v1
        sage: unpickle_matrix_mod2_dense_v1(r,c,s,s2) == A
        True
        sage: loads(dumps(A)) == A
        True
    """
    from sage.matrix.constructor import Matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    cdef int i, j
    cdef Matrix_mod2_dense A

    A = <Matrix_mod2_dense>Matrix(GF(2),r,c)
    if r == 0 or c == 0:
        return A

    cdef signed char *buf = <signed char*>sig_malloc(size)
    for i from 0 <= i < size:
        buf[i] = data[i]

    sig_on()
    cdef gdImagePtr im = gdImageCreateFromPngPtr(size, buf)
    sig_off()

    sig_free(buf)

    if gdImageSX(im) != c or gdImageSY(im) != r:
        raise TypeError("Pickled data dimension doesn't match.")


    for i from 0 <= i < r:
        for j from 0 <= j < c:
            mzd_write_bit(A._entries, i, j, 1-gdImageGetPixel(im, j, i))
    gdImageDestroy(im)
    return A

def from_png(filename):
    """
    Returns a dense matrix over GF(2) from a 1-bit PNG image read from
    ``filename``. No attempt is made to verify that the filename string
    actually points to a PNG image.

    INPUT:

    - filename -- a string

    EXAMPLE::

        sage: from sage.matrix.matrix_mod2_dense import from_png, to_png
        sage: A = random_matrix(GF(2),10,10)
        sage: fn = tmp_filename()
        sage: to_png(A, fn)
        sage: B = from_png(fn)
        sage: A == B
        True
    """
    from sage.matrix.constructor import Matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    cdef int i,j,r,c
    cdef Matrix_mod2_dense A

    fn = open(filename,"r") # check filename
    fn.close()

    cdef FILE *f = fopen(filename, "rb")
    sig_on()
    cdef gdImagePtr im = gdImageCreateFromPng(f)
    sig_off()

    c, r = gdImageSX(im), gdImageSY(im)

    A = <Matrix_mod2_dense>Matrix(GF(2),r,c)

    for i from 0 <= i < r:
        for j from 0 <= j < c:
            mzd_write_bit(A._entries, i, j, 1-gdImageGetPixel(im, j, i))
    fclose(f)
    gdImageDestroy(im)
    return A

def to_png(Matrix_mod2_dense A, filename):
    """
    Saves the matrix ``A`` to filename as a 1-bit PNG image.

    INPUT:

    - ``A`` - a matrix over GF(2)
    - ``filename`` - a string for a file in a writable position

    EXAMPLE::

        sage: from sage.matrix.matrix_mod2_dense import from_png, to_png
        sage: A = random_matrix(GF(2),10,10)
        sage: fn = tmp_filename()
        sage: to_png(A, fn)
        sage: B = from_png(fn)
        sage: A == B
        True
    """
    cdef int i,j, r,c
    r, c = A.nrows(), A.ncols()
    if r == 0 or c == 0:
        raise TypeError("Cannot write image with dimensions %d x %d"%(c,r))
    fn = open(filename,"w") # check filename
    fn.close()
    cdef gdImagePtr im = gdImageCreate(c, r)
    cdef FILE * out = fopen(filename, "wb")
    cdef int black = gdImageColorAllocate(im, 0, 0, 0)
    cdef int white = gdImageColorAllocate(im, 255, 255, 255)
    gdImageFilledRectangle(im, 0, 0, c-1, r-1, white)
    for i from 0 <= i < r:
        for j from 0 <= j < c:
            if mzd_read_bit(A._entries, i, j):
                gdImageSetPixel(im, j, i, black )

    gdImagePng(im, out)
    gdImageDestroy(im)
    fclose(out)

def pluq(Matrix_mod2_dense A, algorithm="standard", int param=0):
    """
    Return PLUQ factorization of A.

    INPUT:

    - A -- matrix
    - algorithm

      * 'standard' asymptotically fast (default)
      * 'mmpf' M4RI inspired
      * 'naive' naive cubic

    - param -- either k for 'mmpf' is chosen or matrix multiplication cutoff
      for 'standard' (default: 0)

    EXAMPLE::

        sage: from sage.matrix.matrix_mod2_dense import pluq
        sage: A = random_matrix(GF(2),4,4); A
        [0 1 0 1]
        [0 1 1 1]
        [0 0 0 1]
        [0 1 1 0]

        sage: LU, P, Q = pluq(A)
        sage: LU
        [1 0 1 0]
        [1 1 0 0]
        [0 0 1 0]
        [1 1 1 0]

        sage: P
        [0, 1, 2, 3]

        sage: Q
        [1, 2, 3, 3]
    """
    cdef Matrix_mod2_dense B = A.__copy__()
    cdef mzp_t *p = mzp_init(A._entries.nrows)
    cdef mzp_t *q = mzp_init(A._entries.ncols)

    if algorithm == "standard":
        sig_on()
        mzd_pluq(B._entries, p, q, param)
        sig_off()
    elif algorithm == "mmpf":
        sig_on()
        _mzd_pluq_russian(B._entries, p, q, param)
        sig_off()
    elif algorithm == "naive":
        sig_on()
        _mzd_pluq_naive(B._entries, p, q)
        sig_off()
    else:
        raise ValueError("Algorithm '%s' unknown."%algorithm)

    P = [p.values[i] for i in range(A.nrows())]
    Q = [q.values[i] for i in range(A.ncols())]
    mzp_free(p)
    mzp_free(q)
    return B,P,Q

def ple(Matrix_mod2_dense A, algorithm="standard", int param=0):
    """
    Return PLE factorization of A.

    INPUT:

    - A -- matrix
    - algorithm

      - 'standard' asymptotically fast (default)
      - 'russian' M4RI inspired
      - 'naive' naive cubic

    - param -- either k for 'mmpf' is chosen or matrix multiplication
      cutoff for 'standard' (default: 0)

    EXAMPLE::

        sage: from sage.matrix.matrix_mod2_dense import ple
        sage: A = random_matrix(GF(2),4,4); A
        [0 1 0 1]
        [0 1 1 1]
        [0 0 0 1]
        [0 1 1 0]

        sage: LU, P, Q = ple(A)
        sage: LU
        [1 0 0 1]
        [1 1 0 0]
        [0 0 1 0]
        [1 1 1 0]

        sage: P
        [0, 1, 2, 3]

        sage: Q
        [1, 2, 3, 3]

        sage: A = random_matrix(GF(2),1000,1000)
        sage: ple(A) == ple(A,'russian') == ple(A,'naive')
        True
    """
    cdef Matrix_mod2_dense B = A.__copy__()
    cdef mzp_t *p = mzp_init(A._entries.nrows)
    cdef mzp_t *q = mzp_init(A._entries.ncols)

    if algorithm == 'standard':
        sig_on()
        mzd_ple(B._entries, p, q, param)
        sig_off()
    elif algorithm == "russian":
        sig_on()
        _mzd_ple_russian(B._entries, p, q, param)
        sig_off()
    elif algorithm == "naive":
        sig_on()
        _mzd_ple_naive(B._entries, p, q)
        sig_off()
    else:
        raise ValueError("Algorithm '%s' unknown."%algorithm)

    P = [p.values[i] for i in range(A.nrows())]
    Q = [q.values[i] for i in range(A.ncols())]
    mzp_free(p)
    mzp_free(q)
    return B,P,Q
