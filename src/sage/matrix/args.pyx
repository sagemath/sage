# cython: wraparound=False
# cython: boundscheck=False
"""
Helpers for creating matrices
"""

#*****************************************************************************
#       Copyright (C) 2018 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport cython
from cpython.sequence cimport PySequence_Fast
from cysignals.signals cimport sig_check

MatrixSpace = None

from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.structure.coerce cimport (coercion_model,
        is_numpy_type, py_scalar_parent)
from sage.structure.element cimport Element, RingElement, Vector
from sage.arith.long cimport pyobject_to_long
from sage.misc.misc_c import sized_iter
from sage.categories import monoids
from sage.misc.superseded import deprecation_cython


try:
    from cypari2.gen import Gen
except ImportError:
    Gen = ()


CommutativeMonoids = monoids.Monoids().Commutative()


cdef inline bint element_is_scalar(Element x) noexcept:
    """
    Should this element be considered a scalar (as opposed to a vector)?
    """
    # This is mostly guesswork, but after some experimentation the
    # checks below turned out to work well in practice...
    if isinstance(x, RingElement):
        return True
    if x._parent in CommutativeMonoids:
        return True
    return False


# We disable cyclic garbage collection and subclassing for efficiency.
# This class is mainly meant to be used internally by MatrixArgs and
# typically not in user code.
@cython.final
@cython.no_gc
cdef class SparseEntry:
    """
    Specialized class for dealing with sparse input in
    :class:`MatrixArgs`. An instance of ``SparseEntry`` represents
    one position in a matrix to be constructed. To construct a sparse
    matrix, one would typically make a list of such.

    Previous versions of Sage used a ``dict`` as data structure for
    sparse input, but that is not so suitable because the keys are not
    guaranteed to be of the correct format. There is also the
    performance cost of creating tuples of Python integers.

    Users of this class are expected to know what they are doing, so
    the indices are not checked when constructing a matrix.

    INPUT:

    - ``i``, ``j`` -- row and column index

    - ``entry`` -- value to be put at position `(i,j)`

    EXAMPLES::

        sage: from sage.matrix.args import SparseEntry
        sage: SparseEntry(123, 456, "abc")
        SparseEntry(123, 456, 'abc')
        sage: SparseEntry(1/3, 2/3, x)                                                  # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: unable to convert rational 1/3 to an integer
    """

    def __init__(self, i, j, entry):
        self.i = pyobject_to_long(i)
        self.j = pyobject_to_long(j)
        self.entry = entry

    def __iter__(self):
        """
        Iteration for convenience, such as converting to a tuple.

        EXAMPLES::

            sage: from sage.matrix.args import SparseEntry
            sage: e = SparseEntry(123, 456, "abc")
            sage: tuple(e)
            (123, 456, 'abc')
        """
        yield self.i
        yield self.j
        yield self.entry

    def __repr__(self):
        return f"{type(self).__name__}({self.i}, {self.j}, {self.entry!r})"


# We disable cyclic garbage collection for efficiency. This is not a
# problem because instances should exist only for a short time.
@cython.no_gc
cdef class MatrixArgs:
    """
    Collect arguments for constructing a matrix.

    This class is meant to pass around arguments, for example from the
    global :func:`matrix` constructor to the matrix space or to the
    element class constructor.

    A typical use case is first creating a ``MatrixArgs`` instance,
    possibly adjusting the attributes. This instance can then be passed
    around and a matrix can be constructed from it using the
    :meth:`matrix` method. Also, a flat list can be constructed using
    :meth:`list` or a sparse dict using :meth:`dict`. It is safe to
    construct multiple objects (of the same or a different kind) from
    the same ``MatrixArgs`` instance.

    ``MatrixArgs`` also supports iteration using the :meth:`iter`
    method. This is a more low-level interface.

    When ``MatrixArgs`` produces output, it is first *finalized*. This
    means that all missing attributes are derived or guessed. After
    finalization, you should no longer change the attributes or it will
    end up in an inconsistent state. You can also finalize explicitly by
    calling the :meth:`finalized` method.

    A ``MatrixArgs`` can contain invalid input. This is not checked when
    constructing the ``MatrixArgs`` instance, but it is checked either
    when finalizing or when constructing an object from it.

    .. WARNING::

        Instances of this class should only be temporary, they are not
        meant to be stored long-term.

    EXAMPLES::

        sage: from sage.matrix.args import MatrixArgs
        sage: ma = MatrixArgs(2, 2, (x for x in range(4))); ma
        <MatrixArgs for None; typ=UNKNOWN; entries=<generator ...>>
        sage: ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_FLAT; entries=[0, 1, 2, 3]>

    Many types of input are possible::

        sage: ma = MatrixArgs(2, 2, entries=None); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=ZERO; entries=None>
        sage: ma.matrix()
        [0 0]
        [0 0]
        sage: ma = MatrixArgs(2, 2, entries={}); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 sparse matrices
         over Integer Ring; typ=SEQ_SPARSE; entries=[]>
        sage: ma.matrix()
        [0 0]
        [0 0]
        sage: ma = MatrixArgs(2, 2, entries=[1,2,3,4]); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_FLAT; entries=[1, 2, 3, 4]>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(2, 2, entries=math.pi); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Real Double Field; typ=SCALAR; entries=3.141592653589793>
        sage: ma.matrix()
        [3.141592653589793               0.0]
        [              0.0 3.141592653589793]
        sage: ma = MatrixArgs(2, 2, entries=pi); ma.finalized()                         # needs sage.symbolic
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Symbolic Ring; typ=SCALAR; entries=pi>
        sage: ma.matrix()                                                               # needs sage.symbolic
        [pi  0]
        [ 0 pi]
        sage: ma = MatrixArgs(ZZ, 2, 2, entries={(0,0): 7}); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 sparse matrices
         over Integer Ring; typ=SEQ_SPARSE; entries=[SparseEntry(0, 0, 7)]>
        sage: ma.matrix()
        [7 0]
        [0 0]
        sage: ma = MatrixArgs(ZZ, 2, 2, entries=((1,2), (3,4))); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_SEQ; entries=((1, 2), (3, 4))>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(ZZ, 2, 2, entries=(1,2,3,4)); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_FLAT; entries=(1, 2, 3, 4)>
        sage: ma.matrix()
        [1 2]
        [3 4]

        sage: # needs sage.libs.pari
        sage: ma = MatrixArgs(QQ, entries=pari("[1,2;3,4]")); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Rational Field; typ=SEQ_FLAT; entries=[1, 2, 3, 4]>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(QQ, 2, 2, entries=pari("[1,2,3,4]")); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Rational Field; typ=SEQ_FLAT; entries=[1, 2, 3, 4]>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(QQ, 2, 2, entries=pari("3/5")); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Rational Field; typ=SCALAR; entries=3/5>
        sage: ma.matrix()
        [3/5   0]
        [  0 3/5]

        sage: ma = MatrixArgs(entries=matrix(2,2)); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=MATRIX; entries=[0 0]
                                                [0 0]>
        sage: ma.matrix()
        [0 0]
        [0 0]
        sage: ma = MatrixArgs(2, 2, entries=lambda i,j: 1+2*i+j); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_FLAT; entries=[1, 2, 3, 4]>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(ZZ, 2, 2, entries=lambda i,j: 1+2*i+j); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=CALLABLE; entries=<function ...>>
        sage: ma.matrix()
        [1 2]
        [3 4]

        sage: # needs numpy
        sage: from numpy import array
        sage: ma = MatrixArgs(array([[1,2],[3,4]])); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Integer Ring; typ=SEQ_SEQ; entries=array([[1, 2], [3, 4]])>
        sage: ma.matrix()
        [1 2]
        [3 4]
        sage: ma = MatrixArgs(array([[1.,2.],[3.,4.]])); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices
         over Real Double Field; typ=MATRIX; entries=[1.0 2.0]
                                                     [3.0 4.0]>
        sage: ma.matrix()
        [1.0 2.0]
        [3.0 4.0]
        sage: ma = MatrixArgs(RealField(20), array([[1.,2.],[3.,4.]])); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices over Real Field
         with 20 bits of precision; typ=MATRIX; entries=[1.0 2.0]
                                                        [3.0 4.0]>
        sage: ma.matrix()
        [1.0000 2.0000]
        [3.0000 4.0000]

        sage: # needs sage.graphs
        sage: ma = MatrixArgs(graphs.CycleGraph(3)); ma.finalized()
        <MatrixArgs for Full MatrixSpace of 3 by 3 dense matrices
         over Integer Ring; typ=MATRIX; entries=[0 1 1]
                                                [1 0 1]
                                                [1 1 0]>
        sage: ma.matrix()
        [0 1 1]
        [1 0 1]
        [1 1 0]

        sage: ma = MatrixArgs([vector([0,1], sparse=True),
        ....:                  vector([0,0], sparse=True)], sparse=True)
        sage: ma.finalized()
        <MatrixArgs for Full MatrixSpace of 2 by 2 sparse matrices over
         Integer Ring; typ=SEQ_SPARSE; entries=[SparseEntry(0, 1, 1)]>
        sage: ma.matrix()
        [0 1]
        [0 0]

    Test invalid input::

        sage: MatrixArgs(ZZ, 2, 2, entries='abcd').finalized()
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'abcd' to a matrix
        sage: MatrixArgs(ZZ, 2, 2, entries=MatrixArgs()).finalized()
        Traceback (most recent call last):
        ...
        TypeError: unable to convert <MatrixArgs for None; typ=UNKNOWN; entries=None> to a matrix
    """

    def __cinit__(self):
        self.nrows = -1
        self.ncols = -1
        self.sparse = -1
        self.kwds = {}

    def __init__(self, *args, base_ring=None, nrows=None, ncols=None, entries=None,
                 sparse=None, row_keys=None, column_keys=None, space=None, **kwds):
        """
        Parse arguments for creating a new matrix.

        See :func:`matrix` for documentation.

        This typically does not raise errors for invalid input, except
        when the arguments cannot be parsed.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: MatrixArgs().finalized()
            <MatrixArgs for Full MatrixSpace of 0 by 0 dense matrices over
             Integer Ring; typ=ZERO; entries=None>
            sage: MatrixArgs(1).finalized()
            <MatrixArgs for Full MatrixSpace of 1 by 1 dense matrices over
             Integer Ring; typ=ZERO; entries=None>
            sage: MatrixArgs(1, 1, 3).finalized()
            <MatrixArgs for Full MatrixSpace of 1 by 1 dense matrices over
             Integer Ring; typ=SCALAR; entries=3>
            sage: MatrixArgs(1, 1, 1, 1).finalized()
            Traceback (most recent call last):
            ...
            TypeError: too many arguments in matrix constructor
            sage: MatrixArgs(3, nrows=1, ncols=1).finalized()
            <MatrixArgs for Full MatrixSpace of 1 by 1 dense matrices over
             Integer Ring; typ=SCALAR; entries=3>
            sage: MatrixArgs(3, nrows=1).finalized()
            <MatrixArgs for Full MatrixSpace of 1 by 1 dense matrices over
             Integer Ring; typ=SCALAR; entries=3>
            sage: MatrixArgs(3, ncols=1).finalized()
            <MatrixArgs for Full MatrixSpace of 1 by 1 dense matrices over
             Integer Ring; typ=SCALAR; entries=3>
        """
        if "ring" in kwds.keys():
            deprecation_cython(issue_number=33380, message="ring is deprecated (keyword will be removed in the future). Use base_ring instead", stacklevel=3)
            self.base = kwds.pop("ring")
        else:
            self.base = base_ring

        if nrows is not None:
            self.set_nrows(pyobject_to_long(nrows))
        if ncols is not None:
            self.set_ncols(pyobject_to_long(ncols))
        self.entries = entries
        if sparse is not None:
            self.sparse = sparse
        if row_keys is not None:
            self.set_row_keys(row_keys)
        if column_keys is not None:
            self.set_column_keys(column_keys)
        if space is not None:
            self.set_space(space)
        self.kwds.update(kwds)

        # Parse positional arguments (base, nrows, ncols, entries)
        # where each of them is optional.
        cdef Py_ssize_t argi = 0, argc = len(args)
        if argi == argc: return

        # fast check for certain types of entries which cannot be
        # confused with a base ring or a number.
        arg = args[argc - 1]
        if self.entries is None and isinstance(arg, (list, tuple, dict)):
            self.entries = arg
            argc -= 1
            if argi == argc:
                return

        # check for base ring argument
        if self.base is None and isinstance(args[argi], Parent):
            self.base = args[argi]
            argi += 1
            if argi == argc:
                return

        # check positional nrows and ncols argument
        # even if redundant with row_keys, column_keys given as keywords;
        # but do not check for positional row_keys, column_keys arguments
        # -- we do not allow those, as they would be too easy to
        # confuse with entries
        cdef int k
        cdef long v
        if self.nrows == -1 and self.ncols == -1:
            for k in range(2):
                arg = args[argi]
                if is_numpy_type(type(arg)):
                    import numpy
                    if isinstance(arg, numpy.ndarray):
                        break
                try:
                    v = pyobject_to_long(arg)
                except TypeError:
                    break
                else:
                    if k == 0:
                        self.set_nrows(v)
                    else:
                        self.set_ncols(v)
                    argi += 1
                    if argi == argc:
                        return

        # check for entries argument
        if self.entries is None:
            self.entries = args[argi]
            argi += 1
            if argi == argc: return

        raise TypeError("too many arguments in matrix constructor")

    def __repr__(self):
        """Print representation for debugging"""
        t = "(invalid)"
        if self.typ == MA_ENTRIES_UNKNOWN:
            t = "UNKNOWN"
        elif self.typ == MA_ENTRIES_ZERO:
            t = "ZERO"
        elif self.typ == MA_ENTRIES_SCALAR:
            t = "SCALAR"
        elif self.typ == MA_ENTRIES_SEQ_SEQ:
            t = "SEQ_SEQ"
        elif self.typ == MA_ENTRIES_SEQ_FLAT:
            t = "SEQ_FLAT"
        elif self.typ == MA_ENTRIES_SEQ_SPARSE:
            t = "SEQ_SPARSE"
        elif self.typ == MA_ENTRIES_CALLABLE:
            t = "CALLABLE"
        elif self.typ == MA_ENTRIES_MATRIX:
            t = "MATRIX"
        elif self.typ == MA_ENTRIES_MAPPING:
            t = "MAPPING"
        elif self.typ == MA_ENTRIES_METHOD:
            t = "METHOD"
        elif self.typ == MA_ENTRIES_NDARRAY:
            t = "NDARRAY"
        return f"<{type(self).__name__} for {self.space}; typ={t}; entries={self.entries!r}>"

    def __reduce__(self):
        """
        We intentionally do not support pickling because there should
        not be any use case for it.

        TESTS::

            sage: from sage.matrix.args import MatrixArgs
            sage: dumps(MatrixArgs())
            Traceback (most recent call last):
            ...
            RuntimeError: pickling MatrixArgs instances is not allowed
            sage: copy(MatrixArgs())
            Traceback (most recent call last):
            ...
            RuntimeError: pickling MatrixArgs instances is not allowed
        """
        raise RuntimeError(f"pickling {type(self).__name__} instances is not allowed")

    def __iter__(self):
        """
        Default iteration (dense with conversion).

        EXAMPLES::

            sage: from sage.matrix.args import SparseEntry, MatrixArgs
            sage: ma = MatrixArgs(ZZ, 2, 3, iter(range(6)))
            sage: list(ma)
            [0, 1, 2, 3, 4, 5]
        """
        return self.iter()

    def iter(self, bint convert=True, bint sparse=False):
        """
        Iteration over the entries in the matrix.

        INPUT:

        - ``convert`` -- if ``True``, the entries are converted to the
          base right; if ``False``, the entries are returned as given

        - ``sparse`` -- see OUTPUT below

        OUTPUT: iterator

        - If ``sparse`` is False: yield all entries of the matrix in
          the following order::

            [1 2 3]
            [4 5 6]

        - If ``sparse`` is True: yield instances of
          :class:`SparseEntry`. The indices ``(i, j)`` are guaranteed to
          lie within the matrix. Zero entries in the input are *not*
          skipped.

        .. WARNING::

            If an iterator is given as input to :class:`MatrixArgs`, it
            may be exhausted breaking any further usage. Otherwise, it
            is safe to iterate multiple times.

        EXAMPLES::

            sage: from sage.matrix.args import SparseEntry, MatrixArgs
            sage: ma = MatrixArgs(ZZ, 2, 3, iter(range(6)))
            sage: list(ma.iter())
            [0, 1, 2, 3, 4, 5]
            sage: ma = MatrixArgs(ZZ, 3, 3, [SparseEntry(0, 0, 0)])
            sage: list(ma.iter())
            Traceback (most recent call last):
            ...
            TypeError: dense iteration is not supported for sparse input

        Sparse examples::

            sage: ma = MatrixArgs(3, 3, pi)                                             # needs sage.symbolic
            sage: list(ma.iter(sparse=True))                                            # needs sage.symbolic
            [SparseEntry(0, 0, pi), SparseEntry(1, 1, pi), SparseEntry(2, 2, pi)]
            sage: ma = MatrixArgs(2, 3)
            sage: list(ma.iter(sparse=True))
            []
            sage: ma = MatrixArgs(2, 2, lambda i, j: i > j)
            sage: list(ma.iter(convert=False, sparse=True))
            [SparseEntry(0, 0, False),
             SparseEntry(0, 1, False),
             SparseEntry(1, 0, True),
             SparseEntry(1, 1, False)]
            sage: ma = MatrixArgs(2, 2, {(1,0):88, (0,1):89})
            sage: sorted(tuple(x) for x in ma.iter(sparse=True))
            [(0, 1, 89), (1, 0, 88)]
            sage: ma = MatrixArgs(QQ, 2, 1, {(1,0):88, (0,1):89})
            sage: ma.finalized()
            Traceback (most recent call last):
            ...
            IndexError: invalid column index 1
            sage: ma = MatrixArgs(QQ, 1, 2, {(1,0):88, (0,1):89})
            sage: ma.finalized()
            Traceback (most recent call last):
            ...
            IndexError: invalid row index 1
        """
        self.finalize()

        cdef long i, j
        cdef SparseEntry se
        if self.typ == MA_ENTRIES_ZERO:
            if sparse:
                pass
            else:
                self._ensure_nrows_ncols()
                zero = self.base.zero()
                for i in range(self.nrows):
                    for j in range(self.ncols):
                        sig_check()
                        yield zero
        elif self.typ == MA_ENTRIES_SCALAR:
            self._ensure_nrows_ncols()
            diag = self.entries
            if convert and self.need_to_convert(diag):
                diag = self.base(diag)
            if sparse:
                for i in range(self.nrows):
                    sig_check()
                    yield make_SparseEntry(i, i, diag)
            else:
                zero = self.base.zero()
                for i in range(self.nrows):
                    for j in range(self.ncols):
                        sig_check()
                        yield diag if (i == j) else zero
        elif self.typ == MA_ENTRIES_SEQ_SEQ:
            self._ensure_nrows_ncols()
            row_iter = sized_iter(self.entries, self.nrows)
            for i in range(self.nrows):
                row = sized_iter(next(row_iter), self.ncols)
                for j in range(self.ncols):
                    sig_check()
                    x = next(row)
                    if convert and self.need_to_convert(x):
                        x = self.base(x)
                    if sparse:
                        yield make_SparseEntry(i, j, x)
                    else:
                        yield x
        elif self.typ == MA_ENTRIES_SEQ_FLAT:
            self._ensure_nrows_ncols()
            it = sized_iter(self.entries, self.nrows * self.ncols)
            for i in range(self.nrows):
                for j in range(self.ncols):
                    sig_check()
                    x = next(it)
                    if convert and self.need_to_convert(x):
                        x = self.base(x)
                    if sparse:
                        yield make_SparseEntry(i, j, x)
                    else:
                        yield x
        elif self.typ == MA_ENTRIES_SEQ_SPARSE:
            if sparse:
                for t in self.entries:
                    sig_check()
                    se = <SparseEntry>t
                    if convert and self.need_to_convert(se.entry):
                        se = make_SparseEntry(se.i, se.j, self.base(se.entry))
                    yield se
            else:
                raise TypeError("dense iteration is not supported for sparse input")
        elif self.typ == MA_ENTRIES_CALLABLE:
            f = self.entries
            row_keys = self.row_keys
            if row_keys is None:
                row_keys = range(self.nrows)
            column_keys = self.column_keys
            if column_keys is None:
                column_keys = range(self.ncols)
            for i in row_keys:
                for j in column_keys:
                    sig_check()
                    x = f(i, j)
                    if convert and self.need_to_convert(x):
                        x = self.base(x)
                    if sparse:
                        yield make_SparseEntry(i, j, x)
                    else:
                        yield x
        elif self.typ == MA_ENTRIES_MATRIX:
            m = self.entries
            for i in range(self.nrows):
                for j in range(self.ncols):
                    sig_check()
                    x = m[i, j]
                    if convert and self.need_to_convert(x):
                        x = self.base(x)
                    if sparse:
                        yield make_SparseEntry(i, j, x)
                    else:
                        yield x
        else:
            raise AssertionError(f"unknown MatrixArgs type {self.typ}")

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: len(MatrixArgs(3, 14))
            42
        """
        self.finalize()
        self._ensure_nrows_ncols()
        return self.nrows * self.ncols

    cpdef Matrix matrix(self, bint convert=True):
        """
        Return the entries of the matrix as a Sage Matrix.

        If possible, this returns a direct reference to
        ``self.entries`` without copying.

        INPUT:

        - ``convert`` -- if ``True``, the matrix is guaranteed to have
          the correct parent matrix space. If ``False``, the input matrix
          may be returned even if it lies in the wrong space.

        .. NOTE::

            This changes ``self.entries`` to the returned matrix. This
            means that it is unsafe to access the ``self.entries``
            attribute after calling this method.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: M = matrix(2, 3, range(6), sparse=True)

        ::

            sage: ma = MatrixArgs(M); ma.finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 3 sparse matrices
             over Integer Ring; typ=MATRIX; entries=[0 1 2]
            [3 4 5]>
            sage: ma.matrix()
            [0 1 2]
            [3 4 5]

        ::

            sage: ma = MatrixArgs(M, sparse=False); ma.finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 3 dense matrices
             over Integer Ring; typ=MATRIX; entries=[0 1 2]
                                                    [3 4 5]>
            sage: ma.matrix()
            [0 1 2]
            [3 4 5]

        ::

            sage: ma = MatrixArgs(RDF, M); ma.finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 3 sparse matrices
             over Real Double Field; typ=MATRIX; entries=[0 1 2]
                                                         [3 4 5]>
            sage: ma.matrix(convert=False)
            [0 1 2]
            [3 4 5]
            sage: ma.matrix()
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]

        If we remove our reference to the input ``M``, no copy is made::

            sage: idM = id(M)
            sage: ma = MatrixArgs(M)
            sage: del M
            sage: M = ma.matrix()
            sage: id(M) == idM
            True

        Immutable matrices are never copied::

            sage: M.set_immutable()
            sage: matrix(M) is M
            True
        """
        self.finalize()

        cdef Matrix M
        for _ in range(1):
            if self.typ == MA_ENTRIES_MATRIX:
                safe = self.ref_safe()  # Check this before assigning M
                M = <Matrix>self.entries
                if M._parent is self.space or not convert:
                    if not (safe or M.is_immutable()):
                        M = M.__copy__()
                    break
        else:
            space = self.space
            if not isinstance(space, MatrixSpace):
                space = space.zero().matrix(side='left').parent()
            M = space(self, coerce=convert)

        # Also store the matrix to support multiple calls of matrix()
        self.entries = M
        self.typ = MA_ENTRIES_MATRIX
        return M

    cpdef element(self, bint immutable=False):
        r"""
        Return the matrix or morphism.

        INPUT:

        - ``immutable`` -- boolean; if ``True``, the result will be immutable

        OUTPUT: an element of ``self.space``

        .. NOTE::

            This may change ``self.entries``, making it unsafe to access the
            ``self.entries`` attribute after calling this method.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: M = matrix(2, 3, range(6), sparse=True)
            sage: ma = MatrixArgs(M); ma.finalized()
            <MatrixArgs for
              Full MatrixSpace of 2 by 3 sparse matrices over Integer Ring;
             typ=MATRIX; entries=[0 1 2]
                                 [3 4 5]>
            sage: M2 = ma.element(immutable=True); M2.parent()
            Full MatrixSpace of 2 by 3 sparse matrices over Integer Ring
            sage: M2.is_immutable()
            True

            sage: ma = MatrixArgs(M, row_keys=['u','v'], column_keys=['a','b','c'])
            sage: ma.finalized()
            <MatrixArgs for
              Set of Morphisms
               from Free module generated by {'a', 'b', 'c'} over Integer Ring
                 to Free module generated by {'u', 'v'} over Integer Ring
                 in Category of finite dimensional modules with basis over Integer Ring;
             typ=MATRIX; entries=[0 1 2]
                                 [3 4 5]>
            sage: phi = ma.element(); phi
            Generic morphism:
            From: Free module generated by {'a', 'b', 'c'} over Integer Ring
            To:   Free module generated by {'u', 'v'} over Integer Ring
        """
        self.finalize()
        cdef Matrix M = self.matrix(convert=True)
        if immutable:
            M.set_immutable()
        if isinstance(self.space, MatrixSpace):
            return M
        return self.space(matrix=M, side='left')

    cpdef list list(self, bint convert=True):
        """
        Return the entries of the matrix as a flat list of scalars.

        If possible and ``convert=False``, this returns a direct
        reference to ``self.entries`` without copying.

        INPUT:

        - ``convert`` -- if ``True``, the entries are converted to the base
          ring; otherwise, the entries are returned as given

        .. NOTE::

            This changes ``self.entries`` to the returned list. This
            means that it is unsafe to access the ``self.entries``
            attribute after calling this method.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: L = list(range(6))
            sage: MatrixArgs(2, 3, L).list()
            [0, 1, 2, 3, 4, 5]

        ::

            sage: ma = MatrixArgs(RDF, 2, 3, L)
            sage: ma.list(convert=False)
            [0, 1, 2, 3, 4, 5]
            sage: ma.list()
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

        If we remove our reference to the input ``L`` and
        ``convert=False``, no copy is made::

            sage: idL = id(L)
            sage: ma = MatrixArgs(2, 3, L)
            sage: del L
            sage: L = ma.list(convert=False)
            sage: id(L) == idL
            True
        """
        self.finalize()

        cdef list L
        if self.typ == MA_ENTRIES_SEQ_FLAT and not convert:
            # Try to re-use existing list
            if type(self.entries) is not list:
                L = list(self.entries)
            else:
                safe = self.ref_safe()  # Check this before assigning L
                L = <list>self.entries
                if not safe:
                    L = L[:]
        else:
            L = list(self.iter(convert))

        cdef long N = self.nrows * self.ncols
        if len(L) != N:
            raise ValueError(f"entries has the wrong length (expected {N}, got {len(L)})")

        # Also store the list to support multiple calls of list()
        self.entries = L
        self.typ = MA_ENTRIES_SEQ_FLAT
        return L

    cpdef dict dict(self, bint convert=True):
        """
        Return the entries of the matrix as a :class:`dict`.

        The keys of this :class:`dict` are the nonzero positions ``(i,j)``. The
        corresponding value is the entry at that position. Zero values are skipped.

        If ``convert`` is ``True``, the entries are converted to the base
        ring. Otherwise, the entries are returned as given.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: L = list(range(6))
            sage: MatrixArgs(2, 3, L).dict()
            {(0, 1): 1, (0, 2): 2, (1, 0): 3, (1, 1): 4, (1, 2): 5}

        ::

            sage: ma = MatrixArgs(GF(2), 2, 3, L)
            sage: ma.dict(convert=False)
            {(0, 1): 1, (0, 2): 2, (1, 0): 3, (1, 1): 4, (1, 2): 5}
            sage: ma.dict()
            {(0, 1): 1, (1, 0): 1, (1, 2): 1}
        """
        self.finalize()

        cdef dict D = {}
        for t in self.iter(convert, True):
            se = <SparseEntry>t
            x = se.entry
            if x:
                D[se.i, se.j] = x
        return D

    cpdef int set_column_keys(self, column_keys) except -1:
        """
        Set the column keys with consistency checking.

        If the value was previously set, it must remain the same.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: ma = MatrixArgs(2, 4)
            sage: ma.set_column_keys('xyz')
            Traceback (most recent call last):
            ...
            ValueError: inconsistent column keys: should be of cardinality 4 but got xyz
            sage: ma.set_column_keys('abcd')
            0
            sage: ma.finalized()
            <MatrixArgs for
              Set of Morphisms
               from Free module generated by {'a', 'b', 'c', 'd'} over Integer Ring
                 to Ambient free module of rank 2 over the principal ideal domain
                    Integer Ring
                 in Category of finite dimensional modules with basis over
                    (Dedekind domains and euclidean domains
                     and noetherian rings and infinite enumerated sets
                     and metric spaces);
             typ=ZERO; entries=None>
        """
        if self.column_keys is not None and self.column_keys != column_keys:
            raise ValueError(f"inconsistent column keys: should be {self.column_keys} "
                             f"but got {column_keys}")
        cdef long p = self.ncols
        if p != -1 and p != len(column_keys):
            raise ValueError(f"inconsistent column keys: should be of cardinality {self.ncols} "
                             f"but got {column_keys}")
        self.column_keys = column_keys

    cpdef int set_row_keys(self, row_keys) except -1:
        """
        Set the row keys with consistency checking.

        If the value was previously set, it must remain the same.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: ma = MatrixArgs(2, 4)
            sage: ma.set_row_keys('xyz')
            Traceback (most recent call last):
            ...
            ValueError: inconsistent row keys: should be of cardinality 2 but got xyz
            sage: ma.set_row_keys(['u', 'v'])
            0
            sage: ma.finalized()
            <MatrixArgs for
              Set of Morphisms
               from Ambient free module of rank 4 over the principal ideal domain
                    Integer Ring
                 to Free module generated by {'u', 'v'} over Integer Ring
                 in Category of finite dimensional modules with basis over
                    (Dedekind domains and euclidean domains
                     and noetherian rings and infinite enumerated sets
                     and metric spaces);
             typ=ZERO; entries=None>
        """
        if self.row_keys is not None and self.row_keys != row_keys:
            raise ValueError(f"inconsistent row keys: should be {self.row_keys} "
                             f"but got {row_keys}")
        if self.nrows != -1 and self.nrows != len(row_keys):
            raise ValueError(f"inconsistent row keys: should be of cardinality {self.nrows} "
                             f"but got {row_keys}")
        self.row_keys = row_keys

    cpdef int set_space(self, space) except -1:
        """
        Set inputs from a given matrix space.

        INPUT:

        - ``space`` -- a :class:`MatrixSpace` or a homset of modules with basis

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: ma = MatrixArgs()
            sage: S = MatrixSpace(QQ, 3, 2, sparse=True)
            sage: _ = ma.set_space(S)
            sage: ma.finalized()
            <MatrixArgs for Full MatrixSpace of 3 by 2 sparse matrices
             over Rational Field; typ=ZERO; entries=None>
            sage: M = ma.matrix(); M
            [0 0]
            [0 0]
            [0 0]
            sage: M.parent() is S
            True

        From a homset::

            sage: C = CombinatorialFreeModule(ZZ, ['a', 'b', 'c'])
            sage: R = CombinatorialFreeModule(ZZ, ['u', 'v'])
            sage: S = Hom(C, R); S
            Set of Morphisms
             from Free module generated by {'a', 'b', 'c'} over Integer Ring
               to Free module generated by {'u', 'v'} over Integer Ring
               in Category of finite dimensional modules with basis over Integer Ring
            sage: ma = MatrixArgs()
            sage: _ = ma.set_space(S)
            sage: ma.finalized()
            <MatrixArgs for Set of Morphisms
             from Free module generated by {'a', 'b', 'c'} over Integer Ring
               to Free module generated by {'u', 'v'} over Integer Ring
               in Category of finite dimensional modules with basis over Integer Ring;
             typ=ZERO; entries=None>
        """
        if self.space is not None:
            return 0  # TODO: ??????
        self.space = <Parent?>space
        try:
            self.set_nrows(space.nrows())
            self.set_ncols(space.ncols())
            self.base = space._base
            self.sparse = space.is_sparse()
        except AttributeError:
            self.set_row_keys(space.codomain().basis().keys())
            self.set_column_keys(space.domain().basis().keys())
            self.base = space.base_ring()

    def finalized(self):
        """
        Determine all missing values.

        Depending on the input, this might be a non-trivial operation.
        In some cases, this will do most of the work of constructing the
        matrix. That is actually not a problem since we eventually want
        to construct the matrix anyway. However, care is taken to avoid
        double work or needless checking or copying.

        OUTPUT: ``self``

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: MatrixArgs(pi).finalized()                                            # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: the dimensions of the matrix must be specified
            sage: MatrixArgs(RR, pi).finalized()                                        # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: the dimensions of the matrix must be specified
            sage: MatrixArgs(2, 3, 0.0).finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 3 dense matrices over Real
             Field with 53 bits of precision; typ=ZERO; entries=0.000000000000000>
            sage: MatrixArgs(RR, 2, 3, 1.0).finalized()
            Traceback (most recent call last):
            ...
            TypeError: nonzero scalar matrix must be square

        Check :issue:`19134`::

            sage: matrix(2, 3, [])
            Traceback (most recent call last):
            ...
            ValueError: sequence too short (expected length 6, got 0)
            sage: matrix(ZZ, 2, 3, [])
            Traceback (most recent call last):
            ...
            ValueError: sequence too short (expected length 6, got 0)
            sage: matrix(2, 3, [1])
            Traceback (most recent call last):
            ...
            ValueError: sequence too short (expected length 6, got 1)

        Check github issue #36065:

            sage: # needs sage.rings.number_field
            sage: class MyAlgebraicNumber(sage.rings.qqbar.AlgebraicNumber):
            ....:     def __bool__(self):
            ....:         raise ValueError
            sage: matrix(1, 1, MyAlgebraicNumber(0))
            [0]
            sage: matrix(1, 1, MyAlgebraicNumber(3))
            [3]
            sage: matrix(1, 2, MyAlgebraicNumber(0))
            Traceback (most recent call last):
            ...
            TypeError: scalar matrix must be square if the value cannot be determined to be zero
            sage: matrix(1, 2, MyAlgebraicNumber(3))
            Traceback (most recent call last):
            ...
            TypeError: scalar matrix must be square if the value cannot be determined to be zero
        """
        self.finalize()
        return self

    cdef int finalize(self) except -1:
        # See finalized() for doc
        if self.is_finalized:
            return 0

        if self.typ == MA_ENTRIES_UNKNOWN:
            self.typ = self.get_type()
            if self.typ == MA_ENTRIES_UNKNOWN:
                raise TypeError(f"unable to convert {self.entries!r} to a matrix")

        # Can we assume a square matrix?
        if self.typ & MA_FLAG_ASSUME_SQUARE:
            # TODO: Handle column_keys/row_keys
            if self.ncols == -1:
                if self.nrows != -1:
                    self.ncols = self.nrows
                elif self.typ == MA_ENTRIES_ZERO:
                    # Special case for the zero matrix:
                    # assume zero rows and columns
                    self.nrows = self.ncols = 0
            elif self.nrows == -1:
                self.nrows = self.ncols

        # Process some types to other types
        if not self.typ & MA_FLAG_FINAL:
            if self.typ == MA_ENTRIES_MAPPING:
                self.process_mapping()
            elif self.typ == MA_ENTRIES_METHOD:
                self.process_method()
            elif self.typ == MA_ENTRIES_NDARRAY:
                self.process_ndarray()
            else:
                raise AssertionError(f"unhandled MatrixArgs type {self.typ}")

        # Handle trivial MATRIX type
        if self.typ == MA_ENTRIES_MATRIX:
            m = <Matrix>self.entries
            self.set_nrows(m._nrows)
            self.set_ncols(m._ncols)
            self.setdefault_base(m._parent._base)
            if self.sparse == -1:
                self.sparse = m.is_sparse()

        # Error if size is required
        if self.typ & MA_FLAG_DIM_REQUIRED:
            if (self.nrows == -1 and self.row_keys is None) or (self.ncols == -1 and self.column_keys is None):
                raise TypeError("the dimensions of the matrix must be specified")

        # Determine base in easy cases
        if self.base is None and self.typ & MA_FLAG_BASE_REQUIRED:
            if self.typ == MA_ENTRIES_ZERO:
                self.base = ZZ
            elif self.typ == MA_ENTRIES_SCALAR:
                if isinstance(self.entries, Element):
                    self.base = (<Element>self.entries)._parent
                else:
                    # Can be None if input is bogus
                    self.base = py_scalar_parent(type(self.entries))
            if self.base is None:
                raise TypeError(f"unable to determine base of {self.entries!r}")

        if ((self.nrows == -1 and self.row_keys is None)
                or (self.ncols == -1 and self.column_keys is None)
                or self.base is None):
            # Determine dimensions or base in the cases where we
            # really need to look at the entries.
            if self.typ == MA_ENTRIES_SEQ_SEQ:
                self.finalize_seq_seq()
            elif self.typ == MA_ENTRIES_SEQ_FLAT:
                self.finalize_seq_scalar()
            elif self.typ == MA_ENTRIES_CALLABLE:
                self.finalize_callable()
            else:
                raise AssertionError(f"nrows={self.nrows}  ncols={self.ncols}  base={self.base}  type={self.typ}")

        # Non-zero scalar matrices must be square
        # also ensure type is MA_ENTRIES_ZERO for scalar zero matrices
        if self.typ == MA_ENTRIES_SCALAR:
            try:
                if not self.entries:
                    self.typ = MA_ENTRIES_ZERO
            except Exception:
                # "not self.entries" has failed, self.entries cannot be determined to be zero
                if self.nrows != self.ncols or self.row_keys != self.column_keys:
                    raise TypeError("scalar matrix must be square if the value cannot be determined to be zero")
            if self.typ == MA_ENTRIES_SCALAR and (self.nrows != self.ncols or self.row_keys != self.column_keys):
                # self.typ is still SCALAR -> "not self.entries" has successfully evaluated, to False
                raise TypeError("nonzero scalar matrix must be square")

        if self.sparse == -1:
            self.sparse = (self.typ & MA_FLAG_SPARSE) != 0

        if self.space is None:
            global MatrixSpace
            if MatrixSpace is None:
                from sage.matrix.matrix_space import MatrixSpace
            nrows = self.nrows
            if nrows == -1:
                nrows = None
            ncols = self.ncols
            if ncols == -1:
                ncols = None
            self.space = MatrixSpace(self.base, nrows, ncols,
                                     sparse=self.sparse,
                                     row_keys=self.row_keys,
                                     column_keys=self.column_keys,
                                     **self.kwds)

        self.is_finalized = True

    cdef int set_base_from_entries(self, entries) except -1:
        """
        Set base ring from the given list of entries.
        """
        if entries:
            B = coercion_model.common_parent(*entries)
            if isinstance(B, type):
                B = py_scalar_parent(B)
            self.base = <Parent?>B
        else:
            self.base = ZZ

    cdef int set_seq_flat(self, entries) except -1:
        """
        Set ``self.entries`` to ``entries``, which should be a flat
        list of entries. Also determine base ring if unknown.
        """
        if self.base is None:
            self.set_base_from_entries(entries)
        self.entries = entries
        self.typ = MA_ENTRIES_SEQ_FLAT

    cdef int process_mapping(self) except -1:
        """
        Convert type ``MAPPING`` to ``SEQ_SPARSE`` or ``SEQ_FLAT``
        (the latter if ``sparse`` is ``False``). This also does the rest
        of finalization, such as determining dimensions and base ring.

        EXAMPLES::

            sage: from sage.matrix.args import MatrixArgs
            sage: ma = MatrixArgs({(2,5): 1/2, (4,-3): 1/3})
            sage: ma = MatrixArgs(2, 2, {(-1,0): 2, (0,-1): 1}, sparse=True)
            sage: ma.finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 2 sparse matrices over
             Integer Ring; typ=SEQ_SPARSE; entries=[SparseEntry(...), SparseEntry(...)]>
            sage: ma = MatrixArgs(2, 2, {(-1,0): 2, (0,-1): 1}, sparse=False)
            sage: ma.finalized()
            <MatrixArgs for Full MatrixSpace of 2 by 2 dense matrices over
             Integer Ring; typ=SEQ_FLAT; entries=[0, 1, 2, 0]>
            sage: ma = MatrixArgs(2, 1, {(1,0): 88, (0,1): 89})
            sage: ma.finalized()
            Traceback (most recent call last):
            ...
            IndexError: invalid column index 1
            sage: ma = MatrixArgs(1, 2, {(1,0): 88, (0,1): 89})
            sage: ma.finalized()
            Traceback (most recent call last):
            ...
            IndexError: invalid row index 1
        """
        cdef list seqsparse = []
        cdef list values = []
        cdef long maxrow = -1, maxcol = -1  # Maximum occurring value for indices
        cdef long i, j
        for (i0, j0), x in self.entries.items():
            sig_check()
            i = pyobject_to_long(i0)
            j = pyobject_to_long(j0)
            if i > maxrow:
                maxrow = i
            elif i < 0:  # Negative indices to index from the end
                if self.nrows > 0:
                    i += self.nrows
                if i < 0:
                    raise IndexError(f"invalid row index {i0}")
            if j > maxcol:
                maxcol = j
            elif j < 0:  # Negative indices to index from the end
                if self.ncols > 0:
                    j += self.ncols
                if j < 0:
                    raise IndexError(f"invalid column index {j0}")
            seqsparse.append(make_SparseEntry(i, j, x))
            values.append(x)
        if self.nrows == -1:
            self.nrows = maxrow + 1
        elif maxrow >= self.nrows:
            raise IndexError(f"invalid row index {maxrow}")
        if self.ncols == -1:
            self.ncols = maxcol + 1
        elif maxcol >= self.ncols:
            raise IndexError(f"invalid column index {maxcol}")

        if self.base is None:
            self.set_base_from_entries(values)

        if self.sparse is False:
            # If we actually want a dense result, convert to MA_ENTRIES_SEQ_FLAT
            self.entries = [self.base.zero()] * (self.nrows * self.ncols)
            self.typ = MA_ENTRIES_SEQ_FLAT
            for t in seqsparse:
                se = <SparseEntry>t
                self.entries[se.i*self.ncols + se.j] = se.entry
        else:
            self.entries = seqsparse
            self.typ = MA_ENTRIES_SEQ_SPARSE

    cdef int process_method(self) except -1:
        """
        Convert type ``METHOD`` to ``MATRIX``.
        """
        args = ()
        if self.base is not None:
            args += (self.base,)
        m = self.entries(*args)
        if not isinstance(m, Matrix):
            raise TypeError(f"{self.entries!r} did not return a matrix")
        self.entries = m
        self.typ = MA_ENTRIES_MATRIX

    cdef int process_ndarray(self) except -1:
        """
        Convert type ``NDARRAY`` to ``MATRIX`` or ``SEQ_SEQ``.
        """
        # TODO: this is old code which should be cleaned up and made
        # more efficient
        e = self.entries
        cdef long inrows, incols
        (inrows, incols) = e.shape
        self.set_nrows(inrows)
        self.set_ncols(incols)
        str_dtype = str(e.dtype)

        if not (e.flags.c_contiguous is True or e.flags.f_contiguous is True):
            raise TypeError('numpy matrix must be either c_contiguous or f_contiguous')

        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        from sage.rings.complex_double import CDF

        if 'float32' in str_dtype:
            m = matrix(RDF, inrows, incols, 0)
            m._replace_self_with_numpy32(e)
            self.typ = MA_ENTRIES_MATRIX
            self.setdefault_base(RDF)
        elif 'float64' in str_dtype:
            m = matrix(RDF, inrows, incols, 0)
            m._replace_self_with_numpy(e)
            self.typ = MA_ENTRIES_MATRIX
            self.setdefault_base(RDF)
        elif 'complex64' in str_dtype:
            m = matrix(CDF, inrows, incols, 0)
            m._replace_self_with_numpy32(e)
            self.typ = MA_ENTRIES_MATRIX
            self.setdefault_base(CDF)
        elif 'complex128' in str_dtype:
            m = matrix(CDF, inrows, incols, 0)
            m._replace_self_with_numpy(e)
            self.typ = MA_ENTRIES_MATRIX
            self.setdefault_base(CDF)
        else:
            if 'int' in str_dtype:
                self.setdefault_base(ZZ)
            m = self.entries
            self.typ = MA_ENTRIES_SEQ_SEQ
        self.entries = m

    cdef int finalize_seq_seq(self) except -1:
        """
        Determine missing input (dimensions, base ring) for type
        ``SEQ_SEQ`` and convert to ``SEQ_FLAT`` in the process.
        """
        e = PySequence_Fast(self.entries, "not a sequence")
        self.set_nrows(len(e))
        if self.nrows == 0:
            if self.ncols == -1 and self.column_keys is None:
                self.set_ncols(0)
            self.setdefault_base(ZZ)
            return 0
        elif (self.ncols != -1 or self.column_keys is not None) and self.base is not None:
            # Everything known => OK
            return 0

        # When sparse and given a list of sparse vectors, convert to a sparse sequence
        cdef long i, j
        if isinstance(e[0], Vector) and all(vec.is_sparse() for vec in e):
            if self.base is None:
                self.base = coercion_model.common_parent(*[(<Vector>vec)._parent._base
                                                           for vec in e])
            self.entries = []
            for i, row in enumerate(e):
                for j, val in (<Vector?>row).iteritems():
                    self.entries.append(make_SparseEntry(i, j, val))

            self.set_ncols(max((<Vector>vec)._parent.ambient_module().rank() for vec in e))
            self.typ = MA_ENTRIES_SEQ_SPARSE
            return 0

        # Process everything and convert to SEQ_FLAT
        cdef list entries = []
        cdef long c
        for row in e:
            c = 0
            for entry in row:
                sig_check()
                c += 1
                entries.append(entry)
            self.set_ncols(c)

        self.set_seq_flat(entries)

    cdef int finalize_seq_scalar(self) except -1:
        """
        Determine missing input (dimensions, base ring) for type
        ``SEQ_FLAT``.
        """
        entries = PySequence_Fast(self.entries, "not a sequence")
        cdef long N = len(entries)

        if self.ncols == 0 or self.nrows == 0 or N == 0:
            # If there are no entries, fill in missing dimensions as 0
            if self.nrows == -1:
                self.nrows = 0
            if self.ncols == -1:
                self.ncols = 0
        elif self.ncols == -1 and self.column_keys is None:
            if self.nrows == -1:
                if self.row_keys is None:
                    # Assume row matrix
                    self.nrows = 1
                    self.ncols = N
                else:
                    self.nrows = len(self.row_keys)
                    self.ncols = N // self.nrows
            else:
                self.ncols = N // self.nrows
        elif self.nrows == -1 and self.row_keys is None:
            if self.ncols == -1:
                self.ncols = len(self.column_keys)
            self.nrows = N // self.ncols

        self.set_seq_flat(entries)

    cdef int finalize_callable(self) except -1:
        """
        Determine base ring for type ``CALLABLE`` and convert to
        ``SEQ_FLAT`` in the process.
        """
        # Dimensions are required, so we must determine the base ring.
        # We do this by converting to the SEQ_FLAT type
        f = self.entries
        cdef list entries = []
        cdef long i, j
        for i in range(self.nrows):
            row = <object>i
            for j in range(self.ncols):
                sig_check()
                entries.append(f(row, j))

        self.set_seq_flat(entries)

    cdef inline entries_type get_type(self) except MA_EXCEPT:
        """
        Return the type of ``self.entries``. In some cases, this might
        change ``self.entries``.

        If the entries are invalid, return ``MA_ENTRIES_UNKNOWN``.

        TESTS:

        Check that :issue:`26655` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(9)
            sage: M = MatrixSpace(F, 2, 2)
            sage: A = M([[1, a], [0, 1]])
            sage: M(pari(A))
            [1 a]
            [0 1]

        Constructing a matrix from a PARI ``t_VEC`` or ``t_COL`` with
        ``t_VEC`` or ``t_COL`` elements is currently not supported::

            sage: M(pari([1, a, 0, 1]))                                                 # needs sage.libs.pari sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NameError: name 'a' is not defined
            sage: M(pari([[1, a], [0, 1]]))                                             # needs sage.libs.pari sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NameError: name 'a' is not defined

        Check that :issue:`38221` is fixed::

            sage: # needs sage.groups
            sage: G = CyclicPermutationGroup(7)
            sage: R = GF(2)
            sage: A = G.algebra(R)
            sage: matrix(A, 3, 3, A.zero())
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: matrix(A, 3, 3, A.one())
            [()  0  0]
            [ 0 ()  0]
            [ 0  0 ()]
        """
        # Check basic Python types. This is very fast, so it doesn't
        # hurt to do these first.
        if self.entries is None:
            return MA_ENTRIES_ZERO
        if isinstance(self.entries, (int, float, complex, Integer)):
            if self.entries:
                return MA_ENTRIES_SCALAR
            return MA_ENTRIES_ZERO
        if isinstance(self.entries, (list, tuple)):
            return self.sequence_type()
        if isinstance(self.entries, dict):
            return MA_ENTRIES_MAPPING

        # Note: some objects are callable, iterable and act like a
        # scalar, e.g. polynomials. So the order of these checks
        # matters a lot. Also the more efficient checks should be
        # done first.

        # We check for Element first because that speeds up some
        # other checks.
        cdef bint is_elt = isinstance(self.entries, Element)
        if is_elt and isinstance(self.entries, Matrix):
            return MA_ENTRIES_MATRIX
        if is_elt and self.base is not None and self.entries.parent() == self.base:
            return MA_ENTRIES_SCALAR
        t = type(self.entries)
        try:
            f = t._matrix_
        except AttributeError:
            pass
        else:
            self.entries = f.__get__(self.entries, t)
            return MA_ENTRIES_METHOD
        if is_elt and element_is_scalar(<Element>self.entries):
            return MA_ENTRIES_SCALAR
        if not is_elt and is_numpy_type(t):
            import numpy
            if isinstance(self.entries, numpy.ndarray):
                # Convert to a numpy array if it was a matrix.
                if t is not numpy.ndarray:
                    self.entries = numpy.array(self.entries)
                return MA_ENTRIES_NDARRAY
            return MA_ENTRIES_SCALAR
        if isinstance(self.entries, Gen):  # PARI object
            from sage.libs.pari.convert_sage_matrix import pari_typ_to_entries_type
            return pari_typ_to_entries_type(self)
        if isinstance(self.entries, MatrixArgs):
            # Prevent recursion
            return MA_ENTRIES_UNKNOWN
        if isinstance(self.entries, str):
            # Blacklist strings, we don't want them to be considered a sequence
            return MA_ENTRIES_UNKNOWN
        try:
            self.entries = list(self.entries)
        except TypeError:
            pass
        else:
            return self.sequence_type()
        if callable(self.entries):
            return MA_ENTRIES_CALLABLE
        if is_elt:  # Last resort
            return MA_ENTRIES_SCALAR
        return MA_ENTRIES_UNKNOWN

    cdef entries_type sequence_type(self) except MA_EXCEPT:
        """
        Return the type of ``self.entries``, where ``self.entries``
        is a sequence.

        If the entries are invalid, return ``MA_ENTRIES_UNKNOWN``.
        """
        if not self.entries:
            return MA_ENTRIES_SEQ_FLAT
        x = self.entries[0]
        if isinstance(x, (list, tuple, Vector)):
            return MA_ENTRIES_SEQ_SEQ
        if type(x) is SparseEntry:
            return MA_ENTRIES_SEQ_SPARSE
        if self.nrows != -1 and self.ncols != -1 and self.ncols != 1:
            # Determine type from the number of entries. Unfortunately,
            # this only works if the number of columns is not 1.
            if len(self.entries) == self.nrows:
                return MA_ENTRIES_SEQ_SEQ
            else:
                return MA_ENTRIES_SEQ_FLAT
        if isinstance(x, (int, float, complex)):
            return MA_ENTRIES_SEQ_FLAT
        if isinstance(x, Element) and element_is_scalar(<Element>x):
            return MA_ENTRIES_SEQ_FLAT
        if isinstance(x, str):
            # Blacklist strings, we don't want them to be considered a sequence
            return MA_ENTRIES_UNKNOWN
        try:
            iter(x)
        except TypeError:
            return MA_ENTRIES_SEQ_FLAT
        else:
            return MA_ENTRIES_SEQ_SEQ


cpdef MatrixArgs MatrixArgs_init(space, entries):
    """
    Construct a :class:`MatrixArgs` object from a matrix space and
    entries. This is the typical use in a matrix constructor.

    If the given entries is already an instance of :class:`MatrixArgs`,
    then just set the space and return the same object.

    EXAMPLES::

        sage: from sage.matrix.args import MatrixArgs_init
        sage: S = MatrixSpace(GF(2), 2, 4)
        sage: ma = MatrixArgs_init(S, {(1, 3): 7})
        sage: M = ma.matrix(); M
        [0 0 0 0]
        [0 0 0 1]
        sage: parent(M) is S
        True
    """
    cdef MatrixArgs ret
    if isinstance(entries, MatrixArgs):
        ret = <MatrixArgs>entries
    else:
        ret = MatrixArgs.__new__(MatrixArgs)
        ret.entries = entries
    ret.set_space(space)
    ret.finalize()
    return ret
