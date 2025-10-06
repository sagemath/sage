from cpython.object cimport PyObject
from sage.structure.element cimport Element, Matrix
from sage.structure.parent cimport Parent


cdef enum entries_type:
    # flags appearing in the types below
    MA_FLAG_ASSUME_SQUARE = 0x01_00  # Assume a square matrix if only 1 dimension given
    MA_FLAG_DIM_REQUIRED  = 0x02_00  # Dimensions must be given
    MA_FLAG_BASE_REQUIRED = 0x04_00  # Base must be given or is trivial to determine
    MA_FLAG_FINAL         = 0x10_00  # Can occur after finalize()
    MA_FLAG_SPARSE        = 0x20_00  # Sparse by default

    # types of input entries
    MA_ENTRIES_UNKNOWN    = 0        # anything
    MA_ENTRIES_ZERO       = 0x17_01  # zero matrix
    MA_ENTRIES_SCALAR     = 0x17_02  # single scalar value
    MA_ENTRIES_SEQ_SEQ    = 0x10_03  # list of lists
    MA_ENTRIES_SEQ_FLAT   = 0x10_04  # flat list
    MA_ENTRIES_SEQ_SPARSE = 0x37_05  # list of SparseEntry instances
    MA_ENTRIES_CALLABLE   = 0x13_06  # function for the entries
    MA_ENTRIES_MATRIX     = 0x14_07  # Sage matrix
    MA_ENTRIES_MAPPING    = 0x01_08  # dict
    MA_ENTRIES_METHOD     = 0x00_09  # function returning the matrix
    MA_ENTRIES_NDARRAY    = 0x00_0A  # numpy array

    # value for exception raised
    MA_EXCEPT = -1


cdef class SparseEntry:
    cdef public long i, j
    cdef public object entry


cdef inline SparseEntry make_SparseEntry(long i, long j, entry):
    e = <SparseEntry>SparseEntry.__new__(SparseEntry)
    e.i = i
    e.j = j
    e.entry = entry
    return e


cdef class MatrixArgs:
    # For integers, -1 means unknown
    # For Python objects, None means unknown
    cdef public Parent space  # parent of matrix
    cdef public Parent base   # parent of entries
    cdef public long nrows, ncols
    cdef public object row_keys, column_keys
    cdef public object entries
    cdef entries_type typ
    cdef public bint sparse
    cdef public dict kwds     # **kwds for MatrixSpace()
    cdef bint is_finalized

    cpdef Matrix matrix(self, bint convert=?)
    cpdef element(self, bint immutable=?)
    cpdef list list(self, bint convert=?)
    cpdef dict dict(self, bint convert=?)

    cdef inline bint ref_safe(self) noexcept:
        """
        Can we safely return self.entries without making a copy?
        A refcount of 1 means that self.entries is the only reference.
        """
        return (<PyObject*>self.entries).ob_refcnt == 1

    cdef inline bint need_to_convert(self, x) noexcept:
        """Is ``x`` not an element of ``self.base``?"""
        if not isinstance(x, Element):
            return True
        return (<Element>x)._parent is not self.base

    cdef int set_base_from_entries(self, entries) except -1

    cdef inline int setdefault_base(self, B) except -1:
        """
        Set the base ring if not previously set.
        """
        if self.base is not None:
            return 0
        self.base = <Parent?>B

    cdef inline int set_ncols(self, long n) except -1:
        """
        Set the number of columns with consistency checking: if the
        value was previously set, it must remain the same.
        """
        if n < 0:
            raise ArithmeticError("number of columns must be nonnegative")
        cdef long p = self.ncols
        if p != -1 and p != n:
            raise ValueError(f"inconsistent number of columns: should be {p} "
                             f"but got {n}")
        if self.column_keys is not None and n != len(self.column_keys):
            raise ValueError(f"inconsistent number of columns: should be cardinality of {self.column_keys} "
                             f"but got {n}")
        self.ncols = n

    cdef inline int set_nrows(self, long n) except -1:
        """
        Set the number of rows with consistency checking: if the
        value was previously set, it must remain the same.
        """
        if n < 0:
            raise ArithmeticError("number of rows must be nonnegative")
        cdef long p = self.nrows
        if p != -1 and p != n:
            raise ValueError(f"inconsistent number of rows: should be {p} but got {n}")
        if self.row_keys is not None and n != len(self.row_keys):
            raise ValueError(f"inconsistent number of rows: should be cardinality of {self.row_keys} "
                             f"but got {n}")
        self.nrows = n

    cdef inline int _ensure_nrows_ncols(self) except -1:
        r"""
        Make sure that the number of rows and columns is set.
        If ``row_keys`` or ``column_keys`` is not finite, this can raise an exception.
        """
        if self.nrows == -1:
            self.nrows = len(self.row_keys)
        if self.ncols == -1:
            self.ncols = len(self.column_keys)

    cpdef int set_column_keys(self, column_keys) except -1
    cpdef int set_row_keys(self, row_keys) except -1
    cpdef int set_space(self, space) except -1

    cdef int finalize(self) except -1
    cdef int process_mapping(self) except -1
    cdef int process_method(self) except -1
    cdef int process_ndarray(self) except -1
    cdef int finalize_seq_seq(self) except -1
    cdef int finalize_seq_scalar(self) except -1
    cdef int finalize_callable(self) except -1
    cdef entries_type get_type(self) except MA_EXCEPT
    cdef entries_type sequence_type(self) except MA_EXCEPT

    cdef int set_seq_flat(self, entries) except -1


cpdef MatrixArgs MatrixArgs_init(space, entries)
