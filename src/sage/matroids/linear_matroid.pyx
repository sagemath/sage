r"""
Linear matroids

When `A` is an `r` times `E` matrix, the linear matroid `M[A]` has groundset
`E` and, for independent sets, all `F` subset of `E` such that the columns of
`M[A]` indexed by `F` are linearly independent.

Construction
============

The recommended way to create a linear matroid is by using the
:func:`Matroid() <sage.matroids.constructor.Matroid>` function, with a
representation matrix `A` as input. This function will intelligently choose
one of the dedicated classes :class:`BinaryMatroid`, :class:`TernaryMatroid`,
:class:`QuaternaryMatroid`, :class:`RegularMatroid` when appropriate. However,
invoking the classes directly is possible too. To get access to them, type::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`. In both cases, it is possible to
provide a reduced matrix `B`, to create the matroid induced by `A = [ I B ]`::

    sage: from sage.matroids.advanced import *
    sage: A = Matrix(GF(2), [[1, 0, 0, 1, 1, 0, 1], [0, 1, 0, 1, 0, 1, 1],
    ....:                    [0, 0, 1, 0, 1, 1, 1]])
    sage: B = Matrix(GF(2), [[1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 1]])
    sage: M1 = Matroid(A)
    sage: M2 = LinearMatroid(A)
    sage: M3 = BinaryMatroid(A)
    sage: M4 = Matroid(reduced_matrix=B)
    sage: M5 = LinearMatroid(reduced_matrix=B)
    sage: isinstance(M1, BinaryMatroid)
    True
    sage: M1.equals(M2)
    True
    sage: M1.equals(M3)
    True
    sage: M1 == M4
    True
    sage: M1.is_field_isomorphic(M5)
    True
    sage: M2 == M3  # comparing LinearMatroid and BinaryMatroid always yields False
    False

Class methods
=============

The ``LinearMatroid`` class and its derivatives inherit all methods from the
:mod:`Matroid <sage.matroids.matroid>` and
:mod:`BasisExchangeMatroid <sage.matroids.basis_exchange_matroid>` classes.
See the documentation for these classes for an overview. In addition, the
following methods are available:

- :class:`LinearMatroid`

    - :func:`base_ring() <sage.matroids.linear_matroid.LinearMatroid.base_ring>`
    - :func:`characteristic() <sage.matroids.linear_matroid.LinearMatroid.characteristic>`
    - :func:`representation() <sage.matroids.linear_matroid.LinearMatroid.representation>`
    - :func:`representation_vectors() <sage.matroids.linear_matroid.LinearMatroid.representation_vectors>`
    - :func:`is_field_equivalent() <sage.matroids.linear_matroid.LinearMatroid.is_field_equivalent>`
    - :func:`is_field_isomorphism() <sage.matroids.linear_matroid.LinearMatroid.is_field_isomorphism>`
    - :func:`has_field_minor() <sage.matroids.linear_matroid.LinearMatroid.has_field_minor>`
    - :func:`fundamental_cycle() <sage.matroids.linear_matroid.LinearMatroid.fundamental_cycle>`
    - :func:`fundamental_cocycle() <sage.matroids.linear_matroid.LinearMatroid.fundamental_cocycle>`
    - :func:`cross_ratios() <sage.matroids.linear_matroid.LinearMatroid.cross_ratios>`
    - :func:`cross_ratio() <sage.matroids.linear_matroid.LinearMatroid.cross_ratio>`

    - :func:`linear_extension() <sage.matroids.linear_matroid.LinearMatroid.linear_extension>`
    - :func:`linear_coextension() <sage.matroids.linear_matroid.LinearMatroid.linear_coextension>`
    - :func:`linear_extension_chains() <sage.matroids.linear_matroid.LinearMatroid.linear_extension_chains>`
    - :func:`linear_coextension_cochains() <sage.matroids.linear_matroid.LinearMatroid.linear_coextension_cochains>`
    - :func:`linear_extensions() <sage.matroids.linear_matroid.LinearMatroid.linear_extensions>`
    - :func:`linear_coextensions() <sage.matroids.linear_matroid.LinearMatroid.linear_coextensions>`

    - :meth:`orlik_terao_algebra() <sage.matroids.linear_matroid.LinearMatroid.orlik_terao_algebra>`

- :class:`BinaryMatroid` has all of the :class:`LinearMatroid` ones, and

    - :func:`bicycle_dimension() <sage.matroids.linear_matroid.BinaryMatroid.bicycle_dimension>`
    - :func:`brown_invariant() <sage.matroids.linear_matroid.BinaryMatroid.brown_invariant>`
    - :func:`is_graphic() <sage.matroids.linear_matroid.BinaryMatroid.is_graphic>`

- :class:`TernaryMatroid` has all of the :class:`LinearMatroid` ones, and

    - :func:`bicycle_dimension() <sage.matroids.linear_matroid.TernaryMatroid.bicycle_dimension>`
    - :func:`character() <sage.matroids.linear_matroid.TernaryMatroid.character>`

- :class:`QuaternaryMatroid` has all of the :class:`LinearMatroid` ones, and

    - :func:`bicycle_dimension() <sage.matroids.linear_matroid.QuaternaryMatroid.bicycle_dimension>`

- :class:`RegularMatroid` has all of the :class:`LinearMatroid` ones, and

    - :func:`is_graphic() <sage.matroids.linear_matroid.RegularMatroid.is_graphic>`

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version

Methods
=======
"""
# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from itertools import product

from cpython.object cimport Py_EQ, Py_NE

from sage.data_structures.bitset_base cimport *
from sage.matrix.constructor import matrix
from sage.matrix.matrix2 cimport Matrix
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.lean_matrix cimport (LeanMatrix, GenericMatrix, BinaryMatrix,
                                        TernaryMatrix, QuaternaryMatrix, PlusMinusOneMatrix,
                                        RationalMatrix)
from sage.matroids.matroid cimport Matroid
from sage.matroids.utilities import newlabel, spanning_stars, spanning_forest, lift_cross_ratios
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.richcmp cimport rich_to_bool

cdef bint GF2_not_defined = True
cdef GF2, GF2_one, GF2_zero

cdef bint GF3_not_defined = True
cdef GF3, GF3_one, GF3_zero, GF3_minus_one

# Implementation note: internally we use data structures from lean_matrix
# instead of Sage's standard Matrix datatypes. This was done so we can use
# highly optimized methods in critical places. Our hope is to do away with
# this in the future. To do so, the following needs to be done:
# 1. Modify the ``__init__`` methods (the lean_matrix constructors are
#                                incompatible with Sage's matrix constructors)
# 2. Look for all lines saying ``# Not a Sage matrix operation`` and provide
#    alternative implementations
# 3. Look for all lines saying ``# Deprecated Sage matrix operation`` and
#    provide alternative implementations
# Below is some code, commented out currently, to get you going.

cdef inline gauss_jordan_reduce(LeanMatrix A, columns):
    return A.gauss_jordan_reduce(columns)   # Not a Sage matrix operation

cdef inline characteristic(LeanMatrix A):
    return A.characteristic()   # Not a Sage matrix operation

# Implementation using default Sage matrices

# cdef gauss_jordan_reduce(Matrix A, columns):
#     """
#     Row-reduce so the lexicographically first basis indexes an identity submatrix.
#     """
#     cdef long r = 0
#     cdef list P = []
#     cdef long c, p, row
#     for c in columns:
#         is_pivot = False
#         for row in range(r, A.nrows()):
#             if A.get_unsafe(row, c) != 0:
#                 is_pivot = True
#                 p = row
#                 break
#         if is_pivot:
#             A.swap_rows_c(p, r)
#             A.rescale_row_c(r, A.get_unsafe(r, c) ** (-1), 0)
#             for row in range(A.nrows()):
#                 if row != r and A.get_unsafe(row, c) != 0:
#                     A.add_multiple_of_row_c(row, r, -A.get_unsafe(row, c), 0)
#             P.append(c)
#             r += 1
#         if r == A.nrows():
#             break
#     return P
#
# cdef inline characteristic(LeanMatrix A):
#     # TODO: use caching for increased speed
#     return A.base_ring().characteristic()

cdef class LinearMatroid(BasisExchangeMatroid):
    r"""
    Linear matroids.

    When `A` is an `r` times `E` matrix, the linear matroid `M[A]` has ground
    set `E` and set of independent sets

        `I(A) =\{F \subseteq E :` the columns of `A` indexed by `F` are linearly independent `\}`

    The simplest way to create a LinearMatroid is by giving only a matrix `A`.
    Then, the groundset defaults to ``range(A.ncols())``. Any iterable object
    ``E`` can be given as a groundset. If ``E`` is a list, then ``E[i]`` will
    label the `i`-th column of `A`. Another possibility is to specify a
    *reduced* matrix `B`, to create the matroid induced by `A = [ I B ]`.

    INPUT:

    - ``matrix`` -- (default: ``None``) a matrix whose column vectors
      represent the matroid
    - ``reduced_matrix`` -- (default: ``None``) a matrix `B` such that
      `[I\ \ B]` represents the matroid, where `I` is an identity matrix with
      the same number of rows as `B`. Only one of ``matrix`` and
      ``reduced_matrix`` should be provided.
    - ``groundset`` -- (default: ``None``) an iterable containing the element
      labels. When provided, must have the correct number of elements: the
      number of columns of ``matrix`` or the number of rows plus the number
      of columns of ``reduced_matrix``.
    - ``ring`` -- (default: ``None``) the desired base ring of the matrix. If
      the base ring is different, an attempt will be made to create a new
      matrix with the correct base ring.
    - ``keep_initial_representation`` -- boolean (default: ``True``); whether
      or not an internal copy of the input matrix should be preserved. This
      can help to see the structure of the matroid (e.g. in the case of
      graphic matroids), and makes it easier to look at extensions. However,
      the input matrix may have redundant rows, and sometimes it is desirable
      to store only a row-reduced copy.

    OUTPUT: a ``LinearMatroid`` instance based on the data above

    .. NOTE::

        The recommended way to generate a linear matroid is through the
        :func:`Matroid() <sage.matroids.constructor.Matroid>` function. It
        will automatically choose more optimized classes when present
        (currently :class:`BinaryMatroid`, :class:`TernaryMatroid`,
        :class:`QuaternaryMatroid`, :class:`RegularMatroid`). For direct
        access to the ``LinearMatroid`` constructor, run::

            sage: from sage.matroids.advanced import *

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: A = Matrix(GF(3), 2, 4, [[1, 0, 1, 1], [0, 1, 1, 2]])
        sage: M = LinearMatroid(A)
        sage: M
        Linear matroid of rank 2 on 4 elements represented over the Finite
        Field of size 3
        sage: sorted(M.groundset())
        [0, 1, 2, 3]
        sage: Matrix(M)
        [1 0 1 1]
        [0 1 1 2]
        sage: M = LinearMatroid(A, 'abcd')
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd']
        sage: B = Matrix(GF(3), 2, 2, [[1, 1], [1, 2]])
        sage: N = LinearMatroid(reduced_matrix=B, groundset='abcd')
        sage: M == N
        True
    """
    def __init__(self, matrix=None, groundset=None, reduced_matrix=None, ring=None, keep_initial_representation=True):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: LinearMatroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],        # indirect doctest
            ....:                                     [0, 1, 1, 2, 3]]))
            Linear matroid of rank 2 on 5 elements represented over the Finite
            Field of size 5
        """
        basis = self._setup_internal_representation(matrix, reduced_matrix, ring, keep_initial_representation)
        if groundset is None:
            groundset = list(range(self._A.nrows() + self._A.ncols()))
        else:
            groundset = list(groundset)
            if len(groundset) != self._A.nrows() + self._A.ncols():
                raise ValueError("size of groundset does not match size of matrix")
        BasisExchangeMatroid.__init__(self, groundset, [groundset[i] for i in basis])
        self._zero = self._A.base_ring()(0)
        self._one = self._A.base_ring()(1)

    def __dealloc__(self):
        """
        Deallocate the memory.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = LinearMatroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],    # indirect doctest
            ....:                                         [0, 1, 1, 2, 3]]))
            sage: M = None
        """
        if self._prow is not NULL:
            sig_free(self._prow)
            self._prow = NULL

    cdef list _setup_internal_representation(self, matrix, reduced_matrix, ring, keep_initial_representation):
        """
        Setup the internal representation matrix ``self._A`` and the array of row- and column indices ``self._prow``.

        Return the displayed basis.
        """
        cdef LeanMatrix A
        cdef long r, c
        cdef list P
        if matrix is not None:
            reduced = False
            if not isinstance(matrix, LeanMatrix):
                if matrix.base_ring() is QQ:
                    A = RationalMatrix(matrix.nrows(), matrix.ncols(), M=matrix)
                else:
                    A = GenericMatrix(matrix.nrows(), matrix.ncols(), M=matrix, ring=ring)
            else:
                A = (<LeanMatrix>matrix).copy()   # Deprecated Sage matrix operation
            if keep_initial_representation:
                self._representation = A.copy()   # Deprecated Sage matrix operation
            P = gauss_jordan_reduce(A, range(A.ncols()))
            self._A = A.matrix_from_rows_and_columns(range(len(P)), [c for c in range(matrix.ncols()) if c not in P])
        else:
            reduced = True
            if not isinstance(reduced_matrix, LeanMatrix):
                if reduced_matrix.base_ring() is QQ:
                    self._A = RationalMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix)
                else:
                    self._A = GenericMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix, ring=ring)
            else:
                self._A = (<LeanMatrix>reduced_matrix).copy()   # Deprecated Sage matrix operation
            P = list(range(self._A.nrows()))
        self._prow = <long* > sig_malloc((self._A.nrows() + self._A.ncols()) * sizeof(long))
        if matrix is not None:
            for r in range(len(P)):
                self._prow[P[r]] = r
            r = 0
            for c in range(A.ncols()):
                if c not in P:
                    self._prow[c] = r
                    r += 1
        else:
            for r in range(self._A.nrows()):
                self._prow[r] = r
            for r in range(self._A.ncols()):
                self._prow[self._A.nrows() + r] = r
        return P

    cpdef _forget(self):
        """
        Remove the internal representation matrix.

        When calling ``Matrix(M)`` after this, the lexicographically first
        basis will be used for the identity matrix.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = LinearMatroid(matrix=Matrix(GF(5), [[1, 1, 0, 1, 1],
            ....:                                         [0, 1, 1, 2, 3]]))
            sage: A = Matrix(M)
            sage: M._forget()
            sage: A == Matrix(M)
            False
        """
        self._representation = None

    cpdef base_ring(self):
        """
        Return the base ring of the matrix representing the matroid.

        EXAMPLES::

            sage: M = Matroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],
            ....:                                   [0, 1, 1, 2, 3]]))
            sage: M.base_ring()
            Finite Field of size 5
        """
        return self._A.base_ring()

    cpdef characteristic(self):
        """
        Return the characteristic of the base ring of the matrix representing
        the matroid.

        EXAMPLES::

            sage: M = Matroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],
            ....:                                   [0, 1, 1, 2, 3]]))
            sage: M.characteristic()
            5
        """
        return characteristic(self._A)

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check if ``self.basis() - x + y`` is again a basis. Internal method.
        """
        return self._A.is_nonzero(self._prow[x], self._prow[y])

    cdef int _exchange(self, long x, long y) except -1:
        """
        Put element indexed by ``x`` into basis, taking out element ``y``.
        Assumptions are that this is a valid basis exchange.

        .. NOTE::

            Safe for noncommutative rings.
        """
        cdef long px, py, r
        px = self._prow[x]
        py = self._prow[y]
        piv = self._A.get_unsafe(px, py)
        pivi = piv ** (-1)
        self._A.rescale_row_c(px, pivi, 0)
        self._A.set_unsafe(px, py, pivi + self._one)       # pivoting without column scaling. Add extra so column does not need adjusting
        for r in range(self._A.nrows()):            # if A and A' are the matrices before and after pivoting, then
            a = self._A.get_unsafe(r, py)       # ker[I A] equals ker[I A'] except for the labelling of the columns
            if a and r != px:
                self._A.add_multiple_of_row_c(r, px, -a, 0)
        self._A.set_unsafe(px, py, pivi)
        self._prow[y] = px
        self._prow[x] = py
        BasisExchangeMatroid._exchange(self, x, y)

    cdef  _exchange_value_internal(self, long x, long y):
        r"""
        Return the (x, y) entry of the current representation.
        """
        return self._A.get_unsafe(self._prow[x], self._prow[y])

    # Sage functions

    def _matrix_(self):
        """
        Return a matrix representation of ``self``.

        OUTPUT:

        A matrix. Either this matrix is equal to the one originally supplied
        by the user, or its displayed basis is the lexicographically least
        basis of the matroid.

        EXAMPLES::

            sage: M = Matroid(matrix=Matrix(GF(5), [[1, 1, 0, 1, 1],
            ....:                                   [0, 1, 1, 2, 3]]))
            sage: M._matrix_()
            [1 1 0 1 1]
            [0 1 1 2 3]
            sage: M._forget()
            sage: M._matrix_()
            [1 0 4 4 3]
            [0 1 1 2 3]
        """
        return self.representation()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = Matroid(matrix=Matrix(GF(5), [[1, 1, 0, 1, 1],
            ....:                                   [0, 1, 1, 2, 3]]))
            sage: repr(M)  # indirect doctest
            'Linear matroid of rank 2 on 5 elements represented over the
            Finite Field of size 5'
        """
        return f'Linear matroid of rank {self.rank()} on {self.size()} elements represented over the {self.base_ring()!r}'

    # representations

    cpdef representation(self, B=None, reduced=False, labels=None, order=None, lift_map=None):
        r"""
        Return a matrix representing the matroid.

        Let `M` be a matroid on `n` elements with rank `r`. Let `E` be an
        ordering of the groundset, as output by
        :func:`M.groundset_list() <sage.matroids.basis_exchange_matroid.BasisExchangeMatroid.groundset_list>`.
        A *representation* of the matroid is an `r \times n` matrix with the
        following property. Consider column `i` to be labeled by `E[i]`,
        and denote by `A[F]` the submatrix formed by the columns labeled by
        the subset `F \subseteq E`. Then for all `F \subseteq E`, the columns
        of `A[F]` are linearly independent if and only if `F` is an
        independent set in the matroid.

        A *reduced representation* is a matrix `D` such that `[I\ \ D]` is a
        representation of the matroid, where `I` is an `r \times r` identity
        matrix. In this case, the rows of `D` are considered to be labeled by
        the first `r` elements of the list ``E``, and the columns by the
        remaining `n - r` elements.

        INPUT:

        - ``B`` -- (default: ``None``) a subset of elements. When provided,
          the representation is such that a basis `B'` that maximally
          intersects `B` is an identity matrix.

        - ``reduced`` -- boolean (default: ``False``); when ``True``, return a
          reduced matrix `D` (so `[I\ \  D]` is a representation of the
          matroid). Otherwise return a full representation matrix.

        - ``labels`` -- (default: ``None``) when ``True``, return additionally
          a list of column labels (if ``reduced=False``) or a list of row
          labels and a list of column labels (if ``reduced=True``).
          The default setting, ``None``, will not return the labels for a full
          matrix, but will return the labels for a reduced matrix.

        - ``order`` -- sequence or ``None`` or ``True`` (default: ``None``)

          - when a sequence, it should be an ordering of the groundset
            elements, and the columns (and, in case of a reduced
            representation, rows) will be presented in the given order,
          - when ``None``, use the same ordering that :meth:`groundset_list`
            uses,
          - when ``True``, return a morphism of free modules instead of a matrix.

        - ``lift_map`` -- (default: ``None``) a dictionary containing the cross
          ratios of the representing matrix in its domain. If provided, the
          representation will be transformed by mapping its cross ratios according
          to ``lift_map``.

        OUTPUT:

        - ``A`` -- a full or reduced representation matrix of ``self``; or
        - ``(A, E)`` -- a full representation matrix ``A`` and a list ``E``
          of column labels; or
        - ``(A, R, C)`` -- a reduced representation matrix and a list ``R`` of
          row labels and a list ``C`` of column labels

        If ``B == None`` and ``reduced == False`` and ``order == None`` then
        this method will always output the same matrix (except when
        ``M._forget()`` is called): either the matrix used as input to create
        the matroid, or a matrix in which the lexicographically least basis
        corresponds to an identity. If only ``order`` is not ``None``, the
        columns of this matrix will be permuted accordingly.

        If a ``lift_map`` is provided, then the resulting matrix will be lifted
        using the method
        :func:`lift_cross_ratios() <sage.matroids.utilities.lift_cross_ratios>`
        See the docstring of this method for further details.

        .. NOTE::

            A shortcut for ``M.representation()`` is ``Matrix(M)``.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.representation()
            [1 0 0 0 1 1 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 0 1]
            sage: Matrix(M) == M.representation()
            True
            sage: M.representation(labels=True)
            (
            [1 0 0 0 1 1 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 0 1], ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            )
            sage: M.representation(B='efg')
            [1 1 0 1 1 0 0]
            [1 0 1 1 0 1 0]
            [1 1 1 0 0 0 1]
            sage: M.representation(B='efg', order='efgabcd')
            [1 0 0 1 1 0 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 1 0]
            sage: M.representation(B='abc', reduced=True)
            (
            [0 1 1 1]
            [1 0 1 1]
            [1 1 0 1], ['a', 'b', 'c'], ['d', 'e', 'f', 'g']
            )
            sage: M.representation(B='efg', reduced=True, labels=False,
            ....:                  order='gfeabcd')
            [1 1 1 0]
            [1 0 1 1]
            [1 1 0 1]

            sage: from sage.matroids.advanced import lift_cross_ratios, lift_map, LinearMatroid
            sage: R = GF(7)
            sage: A = Matrix(R, [[1, 0, 6, 1, 2],[6, 1, 0, 0, 1],[0, 6, 3, 6, 0]])
            sage: M = LinearMatroid(reduced_matrix=A)
            sage: M.representation(lift_map=lift_map('sru'))                            # needs sage.rings.finite_rings
            [     1      0      0      1      0      1      1      1]
            [     0      1      0 -z + 1      1      0      0      1]
            [     0      0      1      0      1      -1 z - 1      0]

        As morphisms::

            sage: M = matroids.catalog.Fano()
            sage: A = M.representation(order=True); A
            Generic morphism:
             From: Free module generated by {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
                   over Finite Field of size 2
             To:   Free module generated by {0, 1, 2} over Finite Field of size 2
            sage: print(A._unicode_art_matrix())
              a b c d e f g
            0⎛1 0 0 0 1 1 1⎞
            1⎜0 1 0 1 0 1 1⎟
            2⎝0 0 1 1 1 0 1⎠
            sage: A = M.representation(B='efg', order=True); A
            Generic morphism:
             From: Free module generated by {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
                   over Finite Field of size 2
             To:   Free module generated by {0, 1, 2} over Finite Field of size 2
            sage: print(A._unicode_art_matrix())
              a b c d e f g
            0⎛1 1 0 1 1 0 0⎞
            1⎜1 0 1 1 0 1 0⎟
            2⎝1 1 1 0 0 0 1⎠
            sage: A = M.representation(B='abc', order=True, reduced=True); A
            Generic morphism:
             From: Free module generated by {'d', 'e', 'f', 'g'}
                   over Finite Field of size 2
             To:   Free module generated by {'a', 'b', 'c'} over Finite Field of size 2
            sage: print(A._unicode_art_matrix())
              d e f g
            a⎛0 1 1 1⎞
            b⎜1 0 1 1⎟
            c⎝1 1 0 1⎠
        """
        cdef LeanMatrix A
        column_keys = None
        if order is True:
            order = self.groundset_list()
            column_keys = tuple(order)
        elif order is None:
            order = self.groundset_list()
        else:
            if not frozenset(order) == self.groundset():
                raise ValueError("elements in argument ``order`` do not correspond to groundset of matroid.")
        order_idx = [self._idx[e] for e in order]
        if not reduced:
            if B is None:
                if self._representation is None:
                    B = set()
                    E = self.groundset_list()
                    i = 0
                    C = self.closure(B)
                    while i < len(E):
                        e = E[i]
                        if e in C:
                            i += 1
                        else:
                            B.add(e)
                            C = self.closure(B)
                            i += 1
                    self._representation = self._basic_representation(B)
                A = self._representation
            else:
                B = self._subset_internal(B)
                A = self._basic_representation(B)
            A = A.matrix_from_rows_and_columns(range(A.nrows()), order_idx)
            Am = A._matrix_()
            if lift_map is not None:
                Am = lift_cross_ratios(Am, lift_map)
            if column_keys is not None:
                Am = matrix(Am, row_keys=range(A.nrows()), column_keys=column_keys)
            if labels:
                return Am, order
            return Am
        else:
            if B is None:
                B = frozenset(self.basis())
            else:
                B = self._subset_internal(B)
            A = self._reduced_representation(B)
            R, C = self._current_rows_cols()
            Ri = []
            Ci = []
            Rl = []
            Cl = []
            for e in order:
                try:
                    i = R.index(e)
                    Ri.append(i)
                    Rl.append(e)
                except ValueError:
                    Ci.append(C.index(e))
                    Cl.append(e)
            A = A.matrix_from_rows_and_columns(Ri, Ci)
            Am = A._matrix_()
            if lift_map is not None:
                Am = lift_cross_ratios(Am, lift_map)
            if column_keys is not None:
                Am = matrix(Am, row_keys=tuple(Rl), column_keys=tuple(Cl))
            if labels or (labels is None and column_keys is None):
                return Am, Rl, Cl
            return Am

    cpdef _current_rows_cols(self, B=None):
        """
        Return the current row and column labels of a reduced matrix.

        INPUT:

        - ``B`` -- (default: ``None``) if provided, first find a basis having
          maximal intersection with ``B``

        OUTPUT:

        - ``R`` -- list of row indices; corresponds to the currently used
          internal basis
        - ``C`` -- list of column indices; corresponds to the complement of
          the current internal basis

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: A = M._reduced_representation('efg')
            sage: R, C = M._current_rows_cols()
            sage: (sorted(R), sorted(C))
            (['e', 'f', 'g'], ['a', 'b', 'c', 'd'])
            sage: R, C = M._current_rows_cols(B='abg')
            sage: (sorted(R), sorted(C))
            (['a', 'b', 'g'], ['c', 'd', 'e', 'f'])
        """
        if B is not None:
            self._move_current_basis(B, set())
        basis = self.basis()
        rows = [0] * self.full_rank()
        for e in basis:
            rows[self._prow[self._idx[e]]] = e
        cols = [0] * self.full_corank()
        for e in self.groundset() - basis:
            cols[self._prow[self._idx[e]]] = e
        return rows, cols

    cpdef LeanMatrix _basic_representation(self, B=None):
        """
        Return a basic matrix representation of the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `M` representing the matroid, where `M[B'] = I` for a basis
        `B'` that maximally intersects the given set `B`.
        If not provided, the current basis used internally is chosen for
        `B'`. For a stable representation, use ``self.representation()``.

        .. NOTE::

            The method self.groundset_list() gives the labelling of the
            columns by the elements of the matroid. The matrix returned
            is a LeanMatrix subclass, which is intended for internal use only.
            Use the ``representation()`` method to get a Sage matrix.

        EXAMPLES::

            sage: M = Matroid(reduced_matrix=Matrix(GF(7), [[1, 1, 1],
            ....:                                           [1, 2, 3]]))
            sage: M._basic_representation()
            LeanMatrix instance with 2 rows and 5 columns over Finite Field of
            size 7
            sage: matrix(M._basic_representation([3, 4]))
            [3 6 2 1 0]
            [5 1 6 0 1]
        """
        cdef LeanMatrix A
        cdef long i
        if B is not None:
            self._move_current_basis(B, set())
        basis = self.basis()
        A = type(self._A)(self.full_rank(), self.size(), ring=self._A.base_ring())
        i = 0
        for e in self._E:
            if e in basis:
                C = self.fundamental_cocycle(basis, e)
                for f in C:
                    A.set_unsafe(i, self._idx[f], C[f])
                i += 1
        return A

    cpdef representation_vectors(self):
        """
        Return a dictionary that associates a column vector with each element
        of the matroid.

        .. SEEALSO::

            :meth:`M.representation() <LinearMatroid.representation>`

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: E = M.groundset_list()
            sage: [M.representation_vectors()[e] for e in E]
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 0),
             (1, 1, 1)]
        """
        R = self._matrix_().columns()
        return {e: R[self._idx[e]] for e in self.groundset()}

    cpdef LeanMatrix _reduced_representation(self, B=None):
        r"""
        Return a reduced representation of the matroid, i.e. a matrix `R` such
        that `[I\ \ R]` represents the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `R` forming a reduced representation of the matroid, with
        rows labeled by a basis `B'` that maximally intersects the given set
        `B`. If not provided, the current basis used internally labels the
        rows.

        .. NOTE::

            The matrix returned is a LeanMatrix subclass, which is intended
            for internal use only. Use the ``representation()`` method to get
            a Sage matrix.

        EXAMPLES::

            sage: M = Matroid(reduced_matrix=Matrix(GF(7), [[1, 1, 1],
            ....:                                           [1, 2, 3]]))
            sage: M._reduced_representation()
            LeanMatrix instance with 2 rows and 3 columns over Finite Field of
            size 7
            sage: matrix(M._reduced_representation([3, 4]))
            [2 3 6]
            [6 5 1]
        """
        if B is not None:
            self._move_current_basis(B, set())
        return self._A.copy()   # Deprecated Sage matrix operation

    # (field) isomorphism

    cpdef bint _is_field_isomorphism(self, LinearMatroid other, morphism) noexcept:  # not safe if self == other
        r"""
        Version of :meth:`<LinearMatroid.is_field_isomorphism>` that does no
        type checking.

        INPUT:

        - ``other`` -- matroid; assumed to have the same base ring as ``self``
        - ``morphism`` -- dictionary mapping the groundset of ``self`` to
          the groundset of ``other``

        OUTPUT: boolean

        .. WARNING::

            This method is not safe if ``self == other``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.Fano().delete(['g'])
            sage: N = BinaryMatroid(Matrix(matroids.Wheel(3)))
            sage: morphism = {'a':0, 'b':1, 'c': 2, 'd':4, 'e':5, 'f':3}
            sage: M._is_field_isomorphism(N, morphism)
            True
        """
        # TODO: ensure this is safe for noncommutative rings
        B = self.basis()
        N = self.groundset() - B
        Bo = frozenset([morphism[e] for e in B])
        No = other.groundset() - Bo
        if not other._is_independent(Bo):
            return False

        C = {}
        for e in B:
            C[e] = self._cocircuit(N | set([e]))
            if other._cocircuit(No | set([morphism[e]])) != frozenset([morphism[f] for f in C[e]]):
                return False

        global GF2, GF2_zero, GF2_one, GF2_not_defined
        try:
            if GF2_not_defined:
                GF2 = GF(2)
                GF2_zero = GF2(0)
                GF2_one = GF2(1)
                GF2_not_defined = False
        except ImportError:
            pass
        else:
            if self.base_ring() == GF2:
                return True

        self._set_current_basis(B)
        other._set_current_basis(Bo)
        normalization = {}
        B = set([b for b in B if len(C[b]) > 1])  # coloops are boring
        N = set(N)
        while B:
            found = False
            for e in B:
                Ce = set(C[e])
                Ce.discard(e)
                N2 = set(Ce - N)
                if N2:
                    found = True
                    f = N2.pop()
                    normalization[e] = self._exchange_value(e, f) * normalization[f] / other._exchange_value(morphism[e], morphism[f])
                    B.discard(e)
                    for f in N2:
                        if self._exchange_value(e, f) * normalization[f] != normalization[e] * other._exchange_value(morphism[e], morphism[f]):
                            return False
                    for f in Ce & N:
                        normalization[f] = (self._one / self._exchange_value(e, f)) * normalization[e] * other._exchange_value(morphism[e], morphism[f])
                        N.discard(f)
                    break
            if not found and N:
                normalization[N.pop()] = self._one
        return True

    cpdef is_field_equivalent(self, other):
        """
        Test for matroid representation equality.

        Two linear matroids `M` and `N` with representation matrices `A` and
        `B` are *field equivalent* if they have the same groundset, and the
        identity map between the groundsets is an isomorphism between the
        representations `A` and `B`. That is, one can be turned into the other
        using only row operations and column scaling.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: boolean

        .. SEEALSO::

            :meth:`M.equals() <sage.matroids.matroid.Matroid.equals>`,
            :meth:`M.is_field_isomorphism() <LinearMatroid.is_field_isomorphism>`,
            :meth:`M.is_field_isomorphic() <LinearMatroid.is_field_isomorphic>`

        EXAMPLES:

        A :class:`BinaryMatroid` and
        :class:`LinearMatroid` use different
        representations of the matroid internally, so `` == ``
        yields ``False``, even if the matroids are equal::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.Fano()
            sage: M1 = LinearMatroid(Matrix(M), groundset=M.groundset_list())
            sage: M2 = Matroid(groundset='abcdefg',
            ....:              reduced_matrix=[[0, 1, 1, 1],
            ....:                              [1, 0, 1, 1],
            ....:                              [1, 1, 0, 1]], field=GF(2))
            sage: M.equals(M1)
            True
            sage: M.equals(M2)
            True
            sage: M.is_field_equivalent(M1)
            True
            sage: M.is_field_equivalent(M2)
            True
            sage: M == M1
            False
            sage: M == M2
            True

        ``LinearMatroid`` instances ``M`` and ``N`` satisfy ``M == N`` if the
        representations are equivalent up to row operations and column
        scaling::

            sage: M1 = Matroid(groundset='abcd',
            ....:          matrix=Matrix(GF(7), [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = Matroid(groundset='abcd',
            ....:          matrix=Matrix(GF(7), [[1, 0, 1, 1], [0, 1, 1, 3]]))
            sage: M3 = Matroid(groundset='abcd',
            ....:          matrix=Matrix(GF(7), [[2, 6, 1, 0], [6, 1, 0, 1]]))
            sage: M1.equals(M2)
            True
            sage: M1.equals(M3)
            True
            sage: M1 == M2
            False
            sage: M1 == M3
            True
            sage: M1.is_field_equivalent(M2)
            False
            sage: M1.is_field_equivalent(M3)
            True
            sage: M1.is_field_equivalent(M1)
            True
        """
        if self is other:
            return True
        if self.base_ring() != other.base_ring():
            return False
        if self.groundset() != other.groundset():
            return False
        if self.full_rank() != other.full_rank():
            return False
        morphism = {e: e for e in self.groundset()}
        return self._is_field_isomorphism(other, morphism)

    cpdef is_field_isomorphism(self, other, morphism):
        r"""
        Test if a provided morphism induces a bijection between represented
        matroids.

        Two represented matroids are *field isomorphic* if the bijection
        ``morphism`` between them induces a field equivalence between their
        representation matrices: the matrices are equal up to row operations
        and column scaling. This implies that the matroids are isomorphic, but
        the converse is false: two isomorphic matroids can be represented by
        matrices that are not field equivalent.

        INPUT:

        - ``other`` -- matroid
        - ``morphism`` -- a map from the groundset of ``self`` to the
          groundset of ``other``. See documentation of the
          :meth:`M.is_isomorphism() <sage.matroids.matroid.Matroid.is_isomorphism>`
          method for more on what is accepted as input.

        OUTPUT: boolean

        .. SEEALSO::

            :meth:`M.is_isomorphism() <sage.matroids.matroid.Matroid.is_isomorphism>`,
            :meth:`M.is_field_equivalent() <LinearMatroid.is_field_equivalent>`,
            :meth:`M.is_field_isomorphic() <LinearMatroid.is_field_isomorphic>`

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: N = matroids.catalog.NonFano()
            sage: N.is_field_isomorphism(M, {e:e for e in M.groundset()})
            False

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.Fano().delete(['g'])
            sage: N = LinearMatroid(reduced_matrix=Matrix(GF(2),
            ....:                       [[-1, 0, 1], [1, -1, 0], [0, 1, -1]]))
            sage: morphism = {'a':0, 'b':1, 'c': 2, 'd':4, 'e':5, 'f':3}
            sage: M.is_field_isomorphism(N, morphism)
            True

            sage: M1 = Matroid(groundset=[0, 1, 2, 3], matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = Matroid(groundset=[0, 1, 2, 3], matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 2, 1]]))
            sage: mf1 = {0:0, 1:1, 2:2, 3:3}
            sage: mf2 = {0:0, 1:1, 2:3, 3:2}
            sage: M1.is_field_isomorphism(M2, mf1)
            False
            sage: M1.is_field_isomorphism(M2, mf2)
            True
        """
        from copy import copy
        if self.base_ring() != other.base_ring():
            return False
        if self.full_rank() != other.full_rank():
            return False
        if self.full_corank() != other.full_corank():
            return False
        if not isinstance(morphism, dict):
            mf = {}
            try:
                for e in self.groundset():
                    mf[e] = morphism[e]
            except (IndexError, TypeError, ValueError):
                try:
                    for e in self.groundset():
                        mf[e] = morphism(e)
                except (TypeError, ValueError):
                    raise TypeError("the morphism argument does not seem to be an isomorphism.")
        else:
            mf = morphism
        if self.groundset().difference(mf.keys()):
            raise ValueError("domain of morphism does not contain groundset of this matroid.")
        if other.groundset().difference([mf[e] for e in self.groundset()]):
            raise ValueError("range of morphism does not contain groundset of other matroid.")

        if self != other:
            return self._is_field_isomorphism(other, mf)
        else:
            return self._is_field_isomorphism(copy(other), mf)

    cpdef _fast_isom_test(self, other):
        """
        Fast (field) isomorphism test for some subclasses.

        INPUT:

        - ``other`` -- a ``LinearMatroid`` instance, of the same subclass as
          ``self``

        OUTPUT:

        - ``None`` -- if the test is inconclusive;
        - ``True`` -- if the matroids were found to be field-isomorphic
        - ``False`` -- if the matroids were found to be non-field-isomorphic

        .. NOTE::

            Intended for internal usage, in particular by the
            ``is_field_isomorphic`` method. Matroids are assumed to be in the
            same subclass.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = BinaryMatroid(reduced_matrix=Matrix(GF(2),
            ....:                 [[1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 1]]))
            sage: M2 = LinearMatroid(reduced_matrix=Matrix(GF(2),
            ....:                 [[1, 1, 0, 1], [1, 0, 1, 1], [1, 1, 0, 1]]))
            sage: M3 = BinaryMatroid(reduced_matrix=Matrix(GF(2),
            ....:                 [[1, 1, 0, 1], [1, 0, 1, 1], [1, 1, 1, 0]]))
            sage: M2._fast_isom_test(M1) is None
            True
            sage: M1._fast_isom_test(M2)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matroids.linear_matroid.LinearMatroid'
            object has no attribute '_invariant'
            sage: M1._fast_isom_test(M3) is None
            True
            sage: Matroid(graphs.WheelGraph(6), regular=True)._fast_isom_test(          # needs sage.graphs
            ....:                                           matroids.Wheel(5))
            True
        """
        pass

    def is_field_isomorphic(self, other):
        """
        Test isomorphism between matroid representations.

        Two represented matroids are *field isomorphic* if there is a
        bijection between their groundsets that induces a field equivalence
        between their representation matrices: the matrices are equal up to
        row operations and column scaling. This implies that the matroids are
        isomorphic, but the converse is false: two isomorphic matroids can be
        represented by matrices that are not field equivalent.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: boolean

        .. SEEALSO::

            :meth:`M.is_isomorphic() <sage.matroids.matroid.Matroid.is_isomorphic>`,
            :meth:`M.is_field_isomorphism() <LinearMatroid.is_field_isomorphism>`,
            :meth:`M.is_field_equivalent() <LinearMatroid.is_field_equivalent>`

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = Matroid(graphs.CompleteGraph(4), regular=True)                   # needs sage.graphs
            sage: M1.is_field_isomorphic(M2)                                            # needs sage.graphs
            True
            sage: M3 = Matroid(bases=M1.bases())
            sage: M1.is_field_isomorphic(M3)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matroids.basis_matroid.BasisMatroid' object
            has no attribute 'base_ring'...
            sage: from sage.matroids.advanced import *
            sage: M4 = BinaryMatroid(Matrix(M1))
            sage: M5 = LinearMatroid(reduced_matrix=Matrix(GF(2), [[-1, 0, 1],
            ....:                                    [1, -1, 0], [0, 1, -1]]))
            sage: M4.is_field_isomorphic(M5)
            True

            sage: M1 = Matroid(groundset=[0, 1, 2, 3], matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = Matroid(groundset=[0, 1, 2, 3], matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 2, 1]]))
            sage: M1.is_field_isomorphic(M2)
            True
            sage: M1.is_field_equivalent(M2)
            False
        """
        if self is other:
            return True
        if self.base_ring() != other.base_ring():
            return False
        if len(self) != len(other):
            return False
        if self.full_rank() != other.full_rank():
            return False
        if self.full_rank() == 0 or self.full_corank() == 0:
            return True
        if self.full_rank() == 1:
            return len(self.loops()) == len(other.loops())
        if self.full_corank() == 1:
            return len(self.coloops()) == len(other.coloops())
        if type(self) is type(other):
            T = self._fast_isom_test(other)
            if T is not None:
                return T

        if self._weak_invariant() != other._weak_invariant():
            return False
        PS = self._weak_partition()
        PO = other._weak_partition()
        if len(PS) != len(PO):
            return False
        if len(PS) == len(self):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self._is_field_isomorphism(other, morphism)

        if self._strong_invariant() != other._strong_invariant():
            return False
        PS = self._strong_partition()
        PO = other._strong_partition()
        if len(PS) != len(PO):
            return False
        if len(PS) == len(self):
            morphism = {}
            for i in range(len(self)):
                morphism[min(PS[i])] = min(PO[i])
            return self._is_field_isomorphism(other, morphism)

        return self.nonbases()._equivalence(lambda sf, ot, morph: self._is_field_isomorphism(other, morph), other.nonbases(), PS, PO) is not None

    def __richcmp__(left, right, op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For linear matroids, in particular, this means field
        equivalence.

        .. TODO::

            In a user guide, write about "pitfalls": testing something like
            ``M in S`` could yield ``False``, even if ``N.equals(M)`` is ``True`` for some
            `N` in `S`.

        .. WARNING::

            This method is linked to __hash__. If you override one, you MUST override the other!

        .. SEEALSO::

            :meth:`<LinearMatroid.is_field_equivalent>`

        EXAMPLES:

        See docstring for :meth:`LinearMatroid.equals>` for more::

            sage: M1 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 3]]))
            sage: M3 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[2, 6, 1, 0], [6, 1, 0, 1]]))
            sage: M1.equals(M2)
            True
            sage: M1.equals(M3)
            True
            sage: M1 != M2  # indirect doctest
            True
            sage: M1 == M3  # indirect doctest
            True
        """
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if type(left) is not type(right):
            return NotImplemented
        if left.is_field_equivalent(right):
            return rich_to_bool(op, 0)
        else:
            return rich_to_bool(op, 1)

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: M1 = Matroid(groundset='abcde', matrix=Matrix(GF(7),
            ....:                         [[1, 0, 1, 1, 1], [0, 1, 1, 2, 3]]))
            sage: M2 = Matroid(groundset='abcde', matrix=Matrix(GF(7),
            ....:                         [[0, 1, 1, 2, 3], [1, 0, 1, 1, 1]]))
            sage: hash(M1) == hash(M2)
            True
            sage: M2 = M1.dual()
            sage: hash(M1) == hash(M2)
            False
        """
        return hash((self.groundset(), self.full_rank(), self._weak_invariant()))

    # minors, dual

    cpdef _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``
        - ``deletions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT: matroid

        EXAMPLES::

            sage: M = Matroid(groundset='abcdefgh', ring=GF(5),
            ....: reduced_matrix=[[2, 1, 1, 0],
            ....:                 [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 2]])
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))
            Linear matroid of rank 3 on 5 elements represented over the Finite
            Field of size 5
        """
        cdef LeanMatrix M
        self._move_current_basis(contractions, deletions)
        rows = list(self.basis() - contractions)
        cols = list(self.cobasis() - deletions)
        M = type(self._A)(len(rows), len(cols), ring=self.base_ring())
        for i in range(len(rows)):
            for j in range(len(cols)):
                M.set_unsafe(i, j, self._exchange_value(rows[i], cols[j]))
        return type(self)(reduced_matrix=M, groundset=rows + cols)

    cpdef dual(self):
        r"""
        Return the dual of the matroid.

        Let `M` be a matroid with groundset `E`. If `B` is the set of bases
        of `M`, then the set `\{E - b : b \in B\}` is the set of bases of
        another matroid, the *dual* of `M`.

        If the matroid is represented by `[I_1\ \ A]`, then the dual is
        represented by `[-A^T\ \ I_2]` for appropriately sized identity
        matrices `I_1, I_2`.

        OUTPUT: the dual matroid

        EXAMPLES::

            sage: A = Matrix(GF(7), [[1, 1, 0, 1],
            ....:                    [1, 0, 1, 1],
            ....:                    [0, 1, 1, 1]])
            sage: B = - A.transpose()
            sage: Matroid(reduced_matrix=A).dual() == Matroid(
            ....:                             reduced_matrix=B,
            ....:                             groundset=[3, 4, 5, 6, 0, 1, 2])
            True
        """
        cdef LeanMatrix R = -self._reduced_representation().transpose()
        rows, cols = self._current_rows_cols()
        return type(self)(reduced_matrix=R, groundset=cols + rows)

    cpdef has_line_minor(self, k, hyperlines=None, certificate=False):
        r"""
        Test if the matroid has a `U_{2, k}`-minor.

        The matroid `U_{2, k}` is a matroid on `k` elements in which every
        subset of at most 2 elements is independent, and every subset of more
        than two elements is dependent.

        The optional argument ``hyperlines`` restricts the search space: this
        method returns ``True`` if `si(M/F)` is isomorphic to `U_{2, l}` with
        `l \geq k` for some `F` in ``hyperlines``, and ``False`` otherwise.

        INPUT:

        - ``k`` -- the length of the line minor
        - ``hyperlines`` -- (default: ``None``) a set of flats of codimension
          2. Defaults to the set of all flats of codimension 2.
        - ``certificate`` -- (default: ``False``) if ``True`` returns
          ``True, F``, where ``F`` is a flat and ``self.minor(contractions=F)``
          has a `U_{2,k}` restriction or ``False, None``

        OUTPUT: boolean or tuple

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: M.has_line_minor(4)
            True
            sage: M.has_line_minor(5)
            False
            sage: M.has_line_minor(k=4, hyperlines=[['a', 'b', 'c']])
            False
            sage: M.has_line_minor(k=4, hyperlines=[['a', 'b', 'c'],
            ....:                                   ['a', 'b', 'd' ]])
            True
            sage: M.has_line_minor(4, certificate=True)
            (True, frozenset({'a', 'b', 'd'}))
            sage: M.has_line_minor(5, certificate=True)
            (False, None)
            sage: M.has_line_minor(k=4, hyperlines=[['a', 'b', 'c'],
            ....:                                   ['a', 'b', 'd' ]], certificate=True)
            (True, frozenset({'a', 'b', 'd'}))
        """
        try:
            if k > len(self.base_ring()) + 1:
                if certificate:
                    return False, None
                return False
        except TypeError:
            pass
        return Matroid.has_line_minor(self, k, hyperlines, certificate)

    cpdef has_field_minor(self, N):
        """
        Check if ``self`` has a minor field isomorphic to ``N``.

        INPUT:

        - ``N`` -- matroid

        OUTPUT: boolean

        .. SEEALSO::

            :meth:`M.minor() <sage.matroids.matroid.Matroid.minor>`,
            :meth:`M.is_field_isomorphic() <LinearMatroid.is_field_isomorphic>`

        .. TODO::

            This important method can (and should) be optimized considerably.
            See [Hli2006]_ p.1219 for hints to that end.

        EXAMPLES::

            sage: M = matroids.Whirl(3)
            sage: matroids.catalog.Fano().has_field_minor(M)
            False
            sage: matroids.catalog.NonFano().has_field_minor(M)
            True
        """
        if self is N:
            return True
        if not isinstance(N, Matroid):
            raise ValueError("N must be a matroid.")
        if self.base_ring() != N.base_ring():
            return False
        rd = self.full_rank() - N.full_rank()
        cd = self.full_corank() - N.full_corank()
        if rd < 0 or cd < 0:
            return False

        R = self._reduced_representation()
        M = type(self)(reduced_matrix=R)

        F = M.flats(rd)
        G = M.dual().flats(cd)
        a = len(M.loops())
        b = len(M.coloops())
        for X in F:
            XB = M.max_independent(X)
            for Y in G:
                YB = M.max_coindependent(Y - XB)
                if len(YB) == cd and len((X - XB) - YB) <= a and len((Y - YB) - XB) <= b and N.is_field_isomorphic(M._minor(contractions=XB, deletions=YB)):
                    return True
        return False

    # cycles, cocycles and cross ratios

    cpdef _exchange_value(self, e, f):
        """
        Return the matrix entry indexed by row `e` and column `f`.

        INPUT:

        - ``e`` -- an element of the groundset
        - ``f`` -- an element of the groundset

        ``e`` should be in the currently active basis, and ``f`` in the
        currently active cobasis.

        OUTPUT:

        The (internal) matrix entry indexed by row `e` and column `f`.

        .. WARNING::

            Intended for internal use. This method does no checks of any kind.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 1, 1], [0, 1, 1, 4]]))
            sage: M._exchange_value(1, 3)
            4
        """
        return self._exchange_value_internal(self._idx[e], self._idx[f])

    cpdef fundamental_cycle(self, B, e):
        """
        Return the fundamental cycle, relative to ``B``, containing element
        ``e``.

        This is the
        :meth:`fundamental circuit <sage.matroids.matroid.Matroid.fundamental_circuit>`
        together with an appropriate signing from the field, such that `Av=0`,
        where `A` is the representation matrix, and `v` the vector
        corresponding to the output.

        INPUT:

        - ``B`` -- a basis of the matroid
        - ``e`` -- an element outside the basis

        OUTPUT:

        A dictionary mapping elements of ``M.fundamental_circuit(B, e)`` to
        elements of the ring.

        .. SEEALSO::

            :meth:`M.fundamental_circuit() <sage.matroids.matroid.Matroid.fundamental_circuit>`

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 1, 1], [0, 1, 1, 4]]))
            sage: v = M.fundamental_cycle([0, 1], 3)
            sage: [v[0], v[1], v[3]]
            [6, 3, 1]
            sage: frozenset(v.keys()) == M.fundamental_circuit([0, 1], 3)
            True
        """
        B = frozenset(B)
        if e in B or not self._is_basis(B):
            return None
        self._move_current_basis(B, set())
        self._move_current_basis(B, set())
        chain = {}
        chain[e] = self._one
        for f in B:
            x = self._exchange_value(f, e)
            if x != 0:
                chain[f] = -x
        return chain

    cpdef fundamental_cocycle(self, B, e):
        """
        Return the fundamental cycle, relative to ``B``, containing element
        ``e``.

        This is the
        :meth:`fundamental cocircuit <sage.matroids.matroid.Matroid.fundamental_cocircuit>`
        together with an appropriate signing from the field, such that `Av=0`,
        where `A` is a representation matrix of the dual, and `v` the vector
        corresponding to the output.

        INPUT:

        - ``B`` -- a basis of the matroid
        - ``e`` -- an element of the basis

        OUTPUT:

        A dictionary mapping elements of ``M.fundamental_cocircuit(B, e)`` to
        elements of the ring.

        .. SEEALSO::

            :meth:`M.fundamental_cocircuit() <sage.matroids.matroid.Matroid.fundamental_cocircuit>`

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 1, 1], [0, 1, 1, 4]]))
            sage: v = M.fundamental_cocycle([0, 1], 0)
            sage: [v[0], v[2], v[3]]
            [1, 1, 1]
            sage: frozenset(v.keys()) == M.fundamental_cocircuit([0, 1], 0)
            True
        """
        B = frozenset(B)
        if e not in B or not self._is_basis(B):
            return None
        self._move_current_basis(B, set())
        cochain = {}
        cochain[e] = self._one
        for f in self.groundset() - set(B):
            x = self._exchange_value(e, f)
            if x != 0:
                cochain[f] = x
        return cochain

    cpdef _line_ratios(self, F):
        """
        Return the set of nonzero ratios of column entries after contracting
        a rank-`r-2` flat ``F``.

        .. WARNING::

            Intended for internal use. Does no checking.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1, 1],
            ....:                    [0, 1, 0, 1, 2, 4], [0, 0, 1, 3, 2, 5]]))
            sage: sorted(M._line_ratios(set([2])))
            [1, 2, 4]
            sage: sorted(M._line_ratios([0]))
            [1, 5]
        """
        self._move_current_basis(F, set())
        X = self.basis().difference(F)
        a = min(X)
        b = max(X)
        rat = set()
        for c in self.cobasis():
            s = self._exchange_value(a, c)
            if s != 0:
                t = self._exchange_value(b, c)
                if t != 0:
                    rat.add(s * (t ** (-1)))
        return rat

    cpdef _line_length(self, F):
        """
        Return ``len(M.contract(F).simplify())``, where ``F`` is assumed to be
        a flat of rank 2 less than the matroid.

        .. WARNING::

            Intended for internal use. Does no checking.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1, 1],
            ....:                    [0, 1, 0, 1, 2, 4], [0, 0, 1, 3, 2, 5]]))
            sage: M._line_length([2])
            5
            sage: M._line_length([0])
            4
        """
        return 2 + len(self._line_ratios(F))

    cpdef _line_cross_ratios(self, F):
        """
        Return the set of cross ratios of column entries after contracting a
        rank-`r-2` flat ``F``.

        Note that these are only the ordered cross ratios!

        .. WARNING::

            Intended for internal use. Does no checking.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1, 1],
            ....:                    [0, 1, 0, 1, 2, 4], [0, 0, 1, 3, 2, 5]]))
            sage: sorted(M._line_cross_ratios(set([2])))
            [2, 4]
            sage: sorted(M._line_cross_ratios([0]))
            [5]
        """
        cr = set()
        rat = self._line_ratios(F)
        while rat:
            r1 = rat.pop()
            for r2 in rat:
                cr.add(r2 / r1)
        return cr

    cpdef cross_ratios(self, hyperlines=None):
        r"""
        Return the set of cross ratios that occur in this linear matroid.

        Consider the following matrix with columns labeled by
        `\{a, b, c, d\}`.

        .. MATH::

            \begin{matrix}
              1 & 0 & 1 & 1\\
              0 & 1 & x & 1
            \end{matrix}

        The cross ratio of the ordered tuple `(a, b, c, d)` equals `x`. The
        set of all cross ratios of a matroid is the set of cross ratios of all
        such minors.

        INPUT:

        - ``hyperlines`` -- (optional) a set of flats of the matroid, of rank
          `r - 2`, where `r` is the rank of the matroid. If not given, then
          ``hyperlines`` defaults to all such flats.

        OUTPUT:

        A list of all cross ratios of this linearly represented matroid that
        occur in rank-2 minors that arise by contracting a flat ``F`` in
        ``hyperlines`` (so by default, those are all cross ratios).

        .. SEEALSO::

            :meth:`M.cross_ratio() <LinearMatroid.cross_ratio>`

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1, 1],
            ....:                            [0, 1, 0, 1, 2, 4],
            ....:                            [0, 0, 1, 3, 2, 5]]))
            sage: sorted(M.cross_ratios())
            [2, 3, 4, 5, 6]
            sage: M = Matroid(graphs.CompleteGraph(5), regular=True)                    # needs sage.graphs
            sage: M.cross_ratios()                                                      # needs sage.graphs
            set()
        """
        if hyperlines is None:
            hyperlines = self.flats(self.full_rank() - 2)
        CR = set()
        for F in hyperlines:
            CR |= self._line_cross_ratios(F)
        CR2 = set()
        while CR:
            cr = CR.pop()
            asc = set([cr, cr ** (-1), -cr + 1, (-cr + 1) ** (-1), cr / (cr - 1), (cr - 1) / cr])
            CR2.update(asc)
            CR.difference_update(asc)
        return CR2

    cpdef cross_ratio(self, F, a, b, c, d):
        r"""
        Return the cross ratio of the four ordered points ``a, b, c, d``
        after contracting a flat ``F`` of codimension 2.

        Consider the following matrix with columns labeled by
        `\{a, b, c, d\}`.

        .. MATH::

            \begin{bmatrix}
              1 & 0 & 1 & 1\\
              0 & 1 & x & 1
            \end{bmatrix}

        The cross ratio of the ordered tuple `(a, b, c, d)` equals `x`. This
        method looks at such minors where ``F`` is a flat to be contracted,
        and all remaining elements other than ``a, b, c, d`` are deleted.

        INPUT:

        - ``F`` -- a flat of codimension 2
        - ``a``, ``b``, ``c``, ``d`` -- elements of the groundset

        OUTPUT:

        The cross ratio of the four points on the line obtained by
        contracting ``F``.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1, 1],
            ....:                            [0, 1, 0, 1, 2, 4],
            ....:                            [0, 0, 1, 3, 2, 6]]))
            sage: M.cross_ratio([0], 1, 2, 3, 5)
            4

            sage: M = Matroid(ring=GF(7), matrix=[[1, 0, 1, 1], [0, 1, 1, 1]])
            sage: M.cross_ratio(set(), 0, 1, 2, 3)
            Traceback (most recent call last):
            ...
            ValueError: points a, b, c, d do not form a 4-point line in M/F
        """
        F = frozenset(F)
        if not F.issubset(self.groundset()):
            raise ValueError("set F must be subset of the groundset")
        if not self.groundset().issuperset([a, b, c, d]):
            raise ValueError("variables a, b, c, d need to be elements of the matroid")
        if self._closure(F | set([a, b])) != self.groundset():
            raise ValueError("set F must be a flat; F with a, b must span the matroid")
        self._move_current_basis(set([a, b]), set([c, d]))
        cr1 = self._exchange_value(a, c) * self._exchange_value(b, d)
        cr2 = self._exchange_value(a, d) * self._exchange_value(b, c)
        if cr1 == 0 or cr2 == 0 or cr1 == cr2:
            raise ValueError("points a, b, c, d do not form a 4-point line in M/F")
        return cr1 / cr2

    cpdef _line_cross_ratio_test(self, F, x, fundamentals):
        r"""
        Check whether the cross ratios involving a fixed element in a fixed
        rank-2 minor are in a specified subset.

        INPUT:

        - ``F`` -- a flat of codimension 2
        - ``x`` -- an element outside ``F``
        - ``fundamentals`` -- set of fundamental elements

        OUTPUT:

        ``True`` if all cross ratios involving ``x`` are in
        ``fundamentals``; ``False`` otherwise.

        .. NOTE::

            This method is intended for checking extensions of a matroid, so
            it is assumed that the cross ratios of `(M/F)-x` are all in the
            desired subset. Moreover, the set of cross ratios is closed under
            reordering of the elements, i.e. if `x` is in ``fundamentals``
            then also `1/x, 1-x, 1/(1-x), x/(x-1), (x-1)/x` are in it.

        .. WARNING::

            Intended for internal use. No checks whatsoever on validity of
            input.

        EXAMPLES::

            sage: M = Matroid(ring=QQ, reduced_matrix=[[1, 1, 1, 0],
            ....:      [1, 1, 0, 1], [1, 0, 1, 1]])
            sage: M._line_cross_ratio_test(set([0]), 6, set([1]))
            True
            sage: M._line_cross_ratio_test(set([4]), 6, set([1]))
            False
            sage: M._line_cross_ratio_test(set([4]), 6, set([1, 2, 1/2, -1]))
            True
        """
        self._move_current_basis(F, set([x]))
        X = self.basis() - F
        a = min(X)
        b = max(X)
        s = self._exchange_value(a, x)
        t = self._exchange_value(b, x)
        if s == 0 or t == 0:
            return True
        try:
            r = s / t
            for c in self.cobasis():  # Only need to check 2x2 submatrices relative to a fixed basis, because of our assumptions
                s = self._exchange_value(a, c)
                t = self._exchange_value(b, c)
                if s != 0 and t != 0:
                    if not s / t / r in fundamentals:
                        return False
        except ZeroDivisionError:
            return False
        return True

    cpdef _cross_ratio_test(self, x, fundamentals, hyperlines=None):
        r"""
        Test if the cross ratios using a given element of this linear matroid
        are contained in a given set of fundamentals.

        INPUT:

        - ``x`` -- an element of the groundset
        - ``fundamentals`` -- a subset of the base ring
        - ``hyperlines`` -- (optional) a set of flats of ``rank=full_rank-2``

        OUTPUT: boolean ``True`` if each cross ratio using ``x`` is an element of
        ``fundamentals``. If ``hyperlines`` is specified, then the method
        tests if each cross ratio in a minor that arises by contracting a flat
        ``F`` in ``hyperlines`` and uses ``x`` is in ``fundamentals``. If
        ``hyperlines`` is not specified, all flats of codimension 2 are
        tested.

        .. NOTE::

            This method is intended for checking extensions of a matroid, so
            it is assumed that the cross ratios of `M/F\setminus x` are all in
            the desired subset. Moreover, the set of fundamentals is closed
            under reordering of the elements, i.e. if `x \in`
            ``fundamentals`` then also `1/x, 1-x, 1/(1-x), x/(x-1), (x-1)/x`
            are in it.

        EXAMPLES::

            sage: M = Matroid(ring=QQ, reduced_matrix=[[1, 1, 1, 0],
            ....:                                 [1, 1, 0, 1], [1, 0, 1, 1]])
            sage: M._cross_ratio_test(6, set([1]))
            False
            sage: M._cross_ratio_test(6, set([1, 2, 1/2, -1]))
            True
        """
        if hyperlines is None:
            hyperlines = [F for F in self.flats(self.full_rank() - 2) if self._line_length(F) > 2]
        if self.rank(self.groundset() - set([x])) < self.full_rank():
            return True
        for F in hyperlines:
            if not self._line_cross_ratio_test(F, x, fundamentals):
                return False
        return True

    # linear extension
    cpdef linear_extension(self, element, chain=None, col=None):
        r"""
        Return a linear extension of this matroid.

        A *linear extension* of the represented matroid `M` by element `e` is
        a matroid represented by `[A\ \ b]`, where `A` is a representation
        matrix of `M` and `b` is a new column labeled by `e`.

        INPUT:

        - ``element`` -- the name of the new element
        - ``col`` -- (default: ``None``) a column to be appended to
          ``self.representation()``; can be any iterable
        - ``chain`` -- (default: ``None``) a dictionary that maps elements of
          the groundset to elements of the base ring

        OUTPUT:

        A linear matroid `N = M([A\ \ b])`, where `A` is a matrix such that
        the current matroid is `M[A]`, and `b` is either given by ``col`` or
        is a weighted combination of columns of `A`, the weights being given
        by ``chain``.

        .. SEEALSO::

            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`.

        EXAMPLES::

            sage: M = Matroid(ring=GF(2), matrix=[[1, 1, 0, 1, 0, 0],
            ....:                                 [1, 0, 1, 0, 1, 0],
            ....:                                 [0, 1, 1, 0, 0, 1],
            ....:                                 [0, 0, 0, 1, 1, 1]])
            sage: M.linear_extension(6, {0:1, 5: 1}).representation()
            [1 1 0 1 0 0 1]
            [1 0 1 0 1 0 1]
            [0 1 1 0 0 1 1]
            [0 0 0 1 1 1 1]
            sage: M.linear_extension(6, col=[0, 1, 1, 1]).representation()
            [1 1 0 1 0 0 0]
            [1 0 1 0 1 0 1]
            [0 1 1 0 0 1 1]
            [0 0 0 1 1 1 1]
        """
        cdef LeanMatrix cl
        cdef long i
        if element in self.groundset():
            raise ValueError("extension element is already in groundset")
        if self._representation is not None and col is not None:
            R = self.base_ring()
            cl = type(self._representation)(self._representation.nrows(), 1, ring=R)
            i = 0
            for x in col:
                if i == self._representation.nrows():
                    raise ValueError("provided column has too many entries")
                cl.set_unsafe(i, 0, R(x))
                i += 1
            if i < self._representation.nrows():
                raise ValueError("provided column has too few entries")
            E = self._E + (element,)
            return type(self)(matrix=self._representation.augment(cl), groundset=E)
        elif col is not None:
            raise ValueError("can only specify column relative to fixed representation. Run self._matrix_() first.")
        else:
            if not isinstance(chain, dict):
                raise TypeError("chain argument needs to be a dictionary")
            return self._linear_extensions(element, [chain])[0]

    cpdef linear_coextension(self, element, cochain=None, row=None):
        r"""
        Return a linear coextension of this matroid.

        A *linear coextension* of the represented matroid `M` by element `e`
        is a matroid represented by

        .. MATH::

            \begin{bmatrix}
                A  & 0\\
                -c & 1
            \end{bmatrix},

        where `A` is a representation matrix of `M`, `c` is a new row, and the
        last column is labeled by `e`.

        This is the dual method of
        :meth:`M.linear_extension() <LinearMatroid.linear_extension>`.

        INPUT:

        - ``element`` -- the name of the new element
        - ``row`` -- (default: ``None``) a row to be appended to
          ``self.representation()``; can be any iterable
        - ``cochain`` -- (default: ``None``) a dictionary that maps elements
          of the groundset to elements of the base ring

        OUTPUT:

        A linear matroid `N = M([A 0; -c 1])`, where `A` is a matrix such that
        the current matroid is `M[A]`, and `c` is either given by ``row``
        (relative to ``self.representation()``) or has nonzero entries given
        by ``cochain``.

        .. NOTE::

            The minus sign is to ensure this method commutes with dualizing.
            See the last example.

        .. SEEALSO::

            :meth:`M.coextension() <sage.matroids.matroid.Matroid.coextension>`,
            :meth:`M.linear_extension() <LinearMatroid.linear_extension>`,
            :meth:`M.dual() <LinearMatroid.dual>`

        EXAMPLES::

            sage: M = Matroid(ring=GF(2), matrix=[[1, 1, 0, 1, 0, 0],
            ....:                                 [1, 0, 1, 0, 1, 0],
            ....:                                 [0, 1, 1, 0, 0, 1],
            ....:                                 [0, 0, 0, 1, 1, 1]])
            sage: M.linear_coextension(6, {0:1, 5: 1}).representation()
            [1 1 0 1 0 0 0]
            [1 0 1 0 1 0 0]
            [0 1 1 0 0 1 0]
            [0 0 0 1 1 1 0]
            [1 0 0 0 0 1 1]
            sage: M.linear_coextension(6, row=[0,1,1,1,0,1]).representation()
            [1 1 0 1 0 0 0]
            [1 0 1 0 1 0 0]
            [0 1 1 0 0 1 0]
            [0 0 0 1 1 1 0]
            [0 1 1 1 0 1 1]

        Coextending commutes with dualizing::

            sage: M = matroids.catalog.NonFano()
            sage: chain = {'a': 1, 'b': -1, 'f': 1}
            sage: M1 = M.linear_coextension('x', chain)
            sage: M2 = M.dual().linear_extension('x', chain)
            sage: M1 == M2.dual()
            True
        """
        cdef LeanMatrix col
        cdef LeanMatrix rw
        cdef long i
        if element in self.groundset():
            raise ValueError("extension element is already in groundset")
        if self._representation is not None and row is not None:
            R = self.base_ring()
            rw = type(self._representation)(1, self._representation.ncols(), ring=R)
            i = 0
            for x in row:
                if i == self._representation.ncols():
                    raise ValueError("provided row has too many entries")
                rw.set_unsafe(0, i, -R(x))
                i += 1
            if i < self._representation.ncols():
                raise ValueError("provided row has too few entries")
            E = self._E + (element,)
            col = type(self._representation)(self._representation.nrows() + 1, 1, ring=self.base_ring())
            col.set_unsafe(self._representation.nrows(), 0, self._one)
            return type(self)(matrix=self._representation.stack(rw).augment(col), groundset=E)
        elif row is not None:
            raise ValueError("can only specify row relative to fixed representation. Run self.representation() first.")
        else:
            if not isinstance(cochain, dict):
                raise TypeError("cochain argument needs to be a dictionary")
            return self._linear_coextensions(element, [cochain])[0]

    cpdef _linear_extensions(self, element, chains):
        r"""
        Return the linear extensions of this matroid representation specified
        by the given chains.

        This is an internal method that does no checking on the input.

        INPUT:

        - ``element`` -- the name of the new element
        - ``chains`` -- list of dictionaries, each of which maps elements of
          the groundset to elements of the base ring

        OUTPUT:

        A list of linear matroids `N = M([A b])`, where `A` is a matrix such
        that the current matroid is `M[A]`, and `b` is a weighted combination
        of columns of `A`, the weights being given by the elements of
        ``chains``.

        EXAMPLES::

            sage: M = Matroid(ring=GF(2), matrix=[[1, 1, 0, 1, 0, 0],
            ....: [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 1]])
            sage: M._linear_extensions(6, [{0:1, 5: 1}])[0].representation()
            [1 1 0 1 0 0 1]
            [1 0 1 0 1 0 1]
            [0 1 1 0 0 1 1]
            [0 0 0 1 1 1 1]
        """
        cdef long i
        cdef LeanMatrix M
        ext = []
        if self._representation is None:
            M = type(self._A)(self.full_rank(), self.size() + 1, self._basic_representation())
        else:
            M = type(self._A)(self._representation.nrows(), self.size() + 1, self._representation)
        E = self._E + (element,)
        D = {E[i]: i for i in range(self.size())}
        for chain in chains:
            for i in range(M.nrows()):
                a = self._zero
                for e in chain:
                    a += M.get_unsafe(i, D[e]) * chain[e]
                M.set_unsafe(i, self.size(), a)
            ext.append(type(self)(matrix=M, groundset=E))
        return ext

    cpdef _linear_coextensions(self, element, cochains):
        r"""
        Return the linear coextensions of this matroid representation
        specified by the given cochains.

        Internal method that does no typechecking.

        INPUT:

        - ``element`` -- the name of the new element
        - ``cochains`` -- list of dictionaries, each of which maps elements
          of the groundset to elements of the base ring

        OUTPUT:

        A list of linear matroids `N = M([A 0; -c 1])`, where `A` is a matrix
        such that the current matroid is `M[A]`, and `c` has nonzero entries
        given by ``cochain``.

        EXAMPLES::

            sage: M = Matroid(ring=GF(2), matrix=[[1, 1, 0, 1, 0, 0],
            ....: [1, 0, 1, 0, 1, 0], [0, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 1]])
            sage: M._linear_coextensions(6, [{0:1, 5: 1}])[0].representation()
            [1 1 0 1 0 0 0]
            [1 0 1 0 1 0 0]
            [0 1 1 0 0 1 0]
            [0 0 0 1 1 1 0]
            [1 0 0 0 0 1 1]
        """
        cdef long i
        cdef LeanMatrix M
        coext = []
        if self._representation is None:
            M = type(self._A)(self.full_rank() + 1, self.size() + 1, self._basic_representation())
        else:
            M = type(self._A)(self._representation.nrows() + 1, self.size() + 1, self._representation)
        M.set_unsafe(M.nrows() - 1, M.ncols() - 1, self._one)
        E = self._E + (element,)
        D = {E[i]: i for i in range(self.size())}
        for cochain in cochains:
            for i in range(M.ncols() - 1):
                M.set_unsafe(M.nrows() - 1, i, 0)
            for e in cochain:
                M.set_unsafe(M.nrows() - 1, D[e], -cochain[e])
            coext.append(type(self)(matrix=M, groundset=E))
        return coext

    cdef _extend_chains(self, C, f, fundamentals=None):
        r"""
        Extend a list of chains for ``self / f`` to a list of chains for
        ``self``.

        See :meth:`linear_extension_chains` for full documentation.
        """
        # assumes connected self, non-loop f, chains with basic support
        R = self.base_ring()
        res = []
        for c in C:
            if not set(c.keys()):
                values = [self._one]
            else:
                if fundamentals is None:
                    values = R
                else:   # generate only chains that satisfy shallow test for 'crossratios in fundamentals'
                    T = frozenset(c.keys())
                    if not self._is_independent(T | set([f])):
                        raise ValueError("_extend_chains can only extend chains with basic support")
                    self._move_current_basis(T | set([f]), set())
                    B = self.basis()
                    mult = {f: self._one}
                    mult2 = {}
                    todo = set([f])
                    todo2 = set()
                    while todo or todo2:
                        while todo:
                            r = todo.pop()
                            cocirc = self.fundamental_cocycle(B, r)
                            for s, v in cocirc.iteritems():
                                if s != r and s not in mult2:
                                    mult2[s] = mult[r] * v
                                    todo2.add(s)
                        while todo2:
                            s = todo2.pop()
                            circ = self.fundamental_cycle(B, s)
                            for t, w in circ.iteritems():
                                if t != s and t not in mult:
                                    mult[t] = mult2[s] / w
                                    if t not in T:
                                        todo.add(t)
                    T2 = set(mult.keys()) & T
                    t = T2.pop()
                    m = -mult[t] * c[t]
                    values = set([fund * m for fund in fundamentals])
                    while T2:
                        t = T2.pop()
                        m = -mult[t] * c[t]
                        values &= set([fund * m for fund in fundamentals])
            for x in values:
                if x != 0:
                    cp = c.copy()
                    cp[f] = x
                    res.append(cp)
            res.append(c)
        ne = newlabel(self._E)
        if fundamentals is not None and self.full_rank() > 1:
            hyperlines = [F for F in self.flats(self.full_rank() - 2) if self._line_length(F) > 2]
            res = [chain for chain in res if len(chain) < 2 or self._linear_extensions(ne, [chain])[0]._cross_ratio_test(ne, fundamentals, hyperlines)]
        return res

    cpdef _linear_extension_chains(self, F, fundamentals=None):  # assumes independent F
        r"""
        Create a list of chains that determine single-element extensions of
        this linear matroid representation.

        .. WARNING::

            Intended for internal use; does no input checking.

        INPUT:

        - ``F`` -- an independent set of elements
        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of chains, so each single-element extension of this linear
        matroid, with support contained in ``F``, is given by one of these
        chains.

        EXAMPLES::

            sage: M = Matroid(reduced_matrix=Matrix(GF(2), [[1, 1, 0],
            ....:                                      [1, 0, 1], [0, 1, 1]]))
            sage: len(M._linear_extension_chains(F=frozenset([0, 1, 2])))
            8
            sage: M._linear_extension_chains(F=frozenset())
            [{}]
            sage: M._linear_extension_chains(F=frozenset([1]))
            [{}, {1: 1}]
            sage: len(M._linear_extension_chains(F=frozenset([0, 1])))
            4
            sage: N = Matroid(ring=QQ, reduced_matrix=[[1, 1, 0],
            ....: [1, 0, 1], [0, 1, 1]])
            sage: N._linear_extension_chains(F=frozenset([0, 1]),
            ....:                           fundamentals=set([1, -1, 1/2, 2]))
            [{0: 1}, {}, {0: 1, 1: 1}, {0: -1, 1: 1}, {1: 1}]
        """

        if not F:
            return [{}]
        if len(F) == 1:
            return [{}, {min(F): self._one}]
        C = self.components()
        if len(C) == 1:
            for f in F:
                sf = frozenset([f])
                ff = self._closure(sf)
                M = self._minor(contractions=sf, deletions=ff - sf)
                if M.is_connected():
                    break
            FM = F & M.groundset()
            chains = M._linear_extension_chains(FM, fundamentals)
            chains = self._extend_chains(chains, f, fundamentals)
        else:
            comp_chains = {}          # make chains for each component
            for comp in C:
                FM = F & comp
                A = self._max_independent(self.groundset() - comp)
                B = self.groundset() - (comp | A)
                M = self._minor(deletions=B, contractions=A)
                M._forget()
                comp_chains[comp] = M._linear_extension_chains(FM, fundamentals)

            chains = [{}]             # make Cartesian product of component chains
            for comp in comp_chains:
                new_chains = []
                for c in chains:
                    for d in comp_chains[comp]:
                        c_new = copy(c)
                        c_new.update(d)
                        new_chains.append(c_new)
                chains = new_chains
        return chains

    cpdef linear_extension_chains(self, F=None, simple=False, fundamentals=None):
        r"""
        Create a list of chains that determine the single-element extensions
        of this linear matroid representation.

        A *chain* is a dictionary, mapping elements from the groundset to
        elements of the base ring, indicating a linear combination of columns
        to form the new column. Think of chains as vectors, only independent
        of representation.

        INPUT:

        - ``F`` -- (default: ``self.groundset()``) a subset of the groundset

        - ``simple`` -- boolean (default: ``False``)

        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of chains, so each single-element extension of this linear
        matroid representation is given by one of these chains.

        If one or more of the above inputs is given, the list is restricted to
        chains

        - so that the support of each chain lies in ``F``, if given
        - so that the chain does not generate a parallel extension or loop, if
          ``simple = True``
        - so that in the extension generated by this chain, the cross ratios
          are restricted to ``fundamentals``, if given.

        .. SEEALSO::

            :meth:`M.linear_extension() <LinearMatroid.linear_extension>`,
            :meth:`M.linear_extensions() <LinearMatroid.linear_extensions>`,
            :meth:`M.cross_ratios() <LinearMatroid.cross_ratios>`

        EXAMPLES::

            sage: M = Matroid(reduced_matrix=Matrix(GF(2),
            ....:                          [[1, 1, 0], [1, 0, 1], [0, 1, 1]]))
            sage: len(M.linear_extension_chains())
            8
            sage: len(M.linear_extension_chains(F=[0, 1]))
            4
            sage: len(M.linear_extension_chains(F=[0, 1], simple=True))
            0
            sage: M.linear_extension_chains(F=[0, 1, 2], simple=True)
            [{0: 1, 1: 1, 2: 1}]
            sage: N = Matroid(ring=QQ,
            ....:         reduced_matrix=[[-1, -1, 0], [1, 0, -1], [0, 1, 1]])
            sage: L = N.linear_extension_chains(F=[0, 1], simple=True,
            ....:                           fundamentals=set([1, -1, 1/2, 2]))
            sage: result = [{0: 1, 1: 1}, {0: -1/2, 1: 1}, {0: -2, 1: 1}]
            sage: all(D in L for D in result)
            True
        """
        if F is None:
            FI = self.basis()
        else:
            FI = self.max_independent(F)
        M = self._minor(contractions=set(), deletions=self.loops())
        M._forget()
        # this is necessary for testing the crossratios of the extension
        # --> skips gauss-jordan reduction when taking minors of M
        # TODO: maybe make this an optional argument for _minor?
        # TODO: make sure the _representation isn't accidentally recreated anywhere (currently this only happens when self.representation() is called)

        chains = M._linear_extension_chains(FI, fundamentals)

        if simple:              # filter out chains that produce a parallel element,
            par = []              # test uses that each supp(chain) is independent
            self._move_current_basis(FI, set())
            B = self.basis()
            for e in self.groundset() - B:
                C = self.fundamental_cycle(B, e)
                C.pop(e)
                par.append(C)
            simple_chains = []
            for c in chains:
                if len(c) < 2:
                    continue
                parallel = False
                for p in par:
                    if set(p.keys()) == set(c.keys()):
                        parallel = True
                        e = min(p)
                        ratio = c[e] / p[e]
                        for f, w in p.iteritems():
                            if c[f] / w != ratio:
                                parallel = False
                                break
                    if parallel:
                        break
                if not parallel:
                    simple_chains.append(c)
            chains = simple_chains
        return chains

    cpdef linear_coextension_cochains(self, F=None, cosimple=False, fundamentals=None):
        r"""
        Create a list of cochains that determine the single-element
        coextensions of this linear matroid representation.

        A cochain is a dictionary, mapping elements from the groundset to
        elements of the base ring. If `A` represents the current matroid, then
        the coextension is given by `N = M([A 0; -c 1])`, with the entries of
        `c` given by the cochain. Note that the matroid does not change when
        row operations are carried out on `A`.

        INPUT:

        - ``F`` -- (default: ``self.groundset()``) a subset of the groundset

        - ``cosimple`` -- boolean (default: ``False``)

        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of cochains, so each single-element coextension of this linear
        matroid representation is given by one of these cochains.

        If one or more of the above inputs is given, the list is restricted to
        chains

        - so that the support of each cochain lies in ``F``, if given
        - so that the cochain does not generate a series extension or coloop,
          if ``cosimple = True``
        - so that in the coextension generated by this cochain, the cross
          ratios are restricted to ``fundamentals``, if given.

        .. SEEALSO::

            :meth:`M.linear_coextension() <LinearMatroid.linear_coextension>`,
            :meth:`M.linear_coextensions() <LinearMatroid.linear_coextensions>`,
            :meth:`M.cross_ratios() <LinearMatroid.cross_ratios>`

        EXAMPLES::

            sage: M = Matroid(reduced_matrix=Matrix(GF(2),
            ....:                          [[1, 1, 0], [1, 0, 1], [0, 1, 1]]))
            sage: len(M.linear_coextension_cochains())
            8
            sage: len(M.linear_coextension_cochains(F=[0, 1]))
            4
            sage: len(M.linear_coextension_cochains(F=[0, 1], cosimple=True))
            0
            sage: M.linear_coextension_cochains(F=[3, 4, 5], cosimple=True)
            [{3: 1, 4: 1, 5: 1}]
            sage: N = Matroid(ring=QQ,
            ....:         reduced_matrix=[[-1, -1, 0], [1, 0, -1], [0, 1, 1]])
            sage: N.linear_coextension_cochains(F=[0, 1], cosimple=True,
            ....:                           fundamentals=set([1, -1, 1/2, 2]))
            [{0: 2, 1: 1}, {0: -1, 1: 1}, {0: 1/2, 1: 1}]
        """
        return self.dual().linear_extension_chains(F=F, simple=cosimple, fundamentals=fundamentals)

    cpdef linear_extensions(self, element=None, F=None, simple=False, fundamentals=None):
        r"""
        Create a list of linear matroids represented by rank-preserving single-element
        extensions of this linear matroid representation.

        INPUT:

        - ``element`` -- (default: ``None``) the name of the new element of
          the groundset

        - ``F`` -- (default: ``None``) a subset of the groundset

        - ``simple`` -- boolean (default: ``False``)

        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of linear matroids represented by rank-preserving single-element extensions of
        this linear matroid representation. In particular, the extension by a coloop is not
        generated.

        If one or more of the above inputs is given, the list is restricted to
        matroids

        - so that the new element is spanned by ``F``, if given
        - so that the new element is not a loop or in a parallel pair, if
          ``simple=True``
        - so that in the representation of the extension, the cross ratios are
          restricted to ``fundamentals``, if given. Note that it is assumed
          that the cross ratios of the input matroid already satisfy this
          condition.

        .. SEEALSO::

            :meth:`M.linear_extension() <LinearMatroid.linear_extension>`,
            :meth:`M.linear_extension_chains() <LinearMatroid.linear_extension_chains>`,
            :meth:`M.cross_ratios() <LinearMatroid.cross_ratios>`

        EXAMPLES::

            sage: M = Matroid(ring=GF(2),
            ....:             reduced_matrix=[[-1, 0, 1], [1, -1, 0], [0, 1, -1]])
            sage: len(M.linear_extensions())
            8
            sage: S = M.linear_extensions(simple=True); S
            [Binary matroid of rank 3 on 7 elements, type (3, 0)]
            sage: S[0].is_field_isomorphic(matroids.catalog.Fano())
            True
            sage: M = Matroid(ring=QQ,
            ....:             reduced_matrix=[[1, 0, 1], [1, 1, 0], [0, 1, 1]])
            sage: S = M.linear_extensions(simple=True,
            ....:                         fundamentals=[1, -1, 1/2, 2])
            sage: len(S)
            7
            sage: any(N.is_isomorphic(matroids.catalog.NonFano())
            ....:     for N in S)
            True
            sage: len(M.linear_extensions(simple=True,
            ....:                         fundamentals=[1, -1, 1/2, 2], F=[0, 1]))
            1
        """
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")
        chains = self.linear_extension_chains(F, simple=simple, fundamentals=fundamentals)
        return self._linear_extensions(element, chains)

    cpdef linear_coextensions(self, element=None, F=None, cosimple=False, fundamentals=None):
        r"""
        Create a list of linear matroids represented by corank-preserving single-element
        coextensions of this linear matroid representation.

        INPUT:

        - ``element`` -- (default: ``None``) the name of the new element of
          the groundset

        - ``F`` -- (default: ``None``) a subset of the groundset

        - ``cosimple`` -- boolean (default: ``False``)

        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of linear matroids represented by corank-preserving single-element
        coextensions of this linear matroid representation. In particular, the coextension
        by a loop is not generated.

        If one or more of the above inputs is given, the list is restricted to
        coextensions

        - so that the new element lies in the cospan of ``F``, if given.
        - so that the new element is no coloop and is not in series with
          another element, if ``cosimple = True``.
        - so that in the representation of the coextension, the cross ratios
          are restricted to ``fundamentals``, if given. Note that it is
          assumed that the cross ratios of the input matroid already satisfy
          this condition.

        .. SEEALSO::

            :meth:`M.linear_coextension() <LinearMatroid.linear_coextension>`,
            :meth:`M.linear_coextension_cochains() <LinearMatroid.linear_coextension_cochains>`,
            :meth:`M.cross_ratios() <LinearMatroid.cross_ratios>`

        EXAMPLES::

            sage: M = Matroid(ring=GF(2),
            ....:         reduced_matrix=[[-1, 0, 1], [1, -1, 0], [0, 1, -1]])
            sage: len(M.linear_coextensions())
            8
            sage: S = M.linear_coextensions(cosimple=True)
            sage: S
            [Binary matroid of rank 4 on 7 elements, type (3, 7)]
            sage: F7 = matroids.catalog.Fano()
            sage: S[0].is_field_isomorphic(F7.dual())
            True
            sage: M = Matroid(ring=QQ,
            ....:            reduced_matrix=[[1, 0, 1], [1, 1, 0], [0, 1, 1]])
            sage: S = M.linear_coextensions(cosimple=True,
            ....:                           fundamentals=[1, -1, 1/2, 2])
            sage: len(S)
            7
            sage: NF7 = matroids.catalog.NonFano()
            sage: any(N.is_isomorphic(NF7.dual()) for N in S)
            True
            sage: len(M.linear_coextensions(cosimple=True,
            ....:                           fundamentals=[1, -1, 1/2, 2],
            ....:                           F=[3, 4]))
            1
        """
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")
        cochains = self.linear_coextension_cochains(F, cosimple=cosimple, fundamentals=fundamentals)
        return self._linear_coextensions(element, cochains)

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data represent an actual matroid.

        Since this matroid is linear, we test the representation matrix.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, or (boolean, dictionary)

        The boolean output value is:

        - ``True`` if the matrix is over a field.
        - ``True`` if the matrix is over a ring and all cross ratios are
          invertible.
        - ``False`` otherwise.

        .. NOTE::

            This function does NOT test if the cross ratios are contained in
            the appropriate set of fundamentals. To that end, use

            ``M.cross_ratios().issubset(F)``

            where ``F`` is the set of fundamentals.

        .. SEEALSO::

            :meth:`M.cross_ratios() <LinearMatroid.cross_ratios>`

        EXAMPLES::

            sage: M = Matroid(ring=QQ, reduced_matrix=Matrix(ZZ,
            ....:                          [[1, 0, 1], [1, 1, 0], [0, 1, 1]]))
            sage: M.is_valid()
            True
            sage: from sage.matroids.advanced import *  # LinearMatroid
            sage: M = LinearMatroid(ring=ZZ, reduced_matrix=Matrix(ZZ,
            ....:                          [[1, 0, 1], [1, 1, 0], [0, 1, 1]]))
            sage: M.is_valid(certificate=True)
            (False, {'error': 'not all cross ratios are invertible'})
        """
        if self.base_ring().is_field():
            return True if not certificate else (True, {})
        try:
            CR = self.cross_ratios()
        except (ArithmeticError, TypeError, ValueError):
            return False if not certificate else (False, {"error": "can't compute cross ratios"})
        for x in CR:
            if not x ** (-1) in self.base_ring():
                return False if not certificate else (False, {"error": "not all cross ratios are invertible"})
        return True if not certificate else (True, {})

    # connectivity

    cpdef _is_3connected_shifting(self, certificate=False):
        r"""
        Return ``True`` if the matroid is 4-connected, ``False`` otherwise. It can
        optionally return a separator as a witness.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); if ``True``,
          then return ``True, None`` if the matroid is 3-connected,
          and ``False,`` `X` otherwise, where `X` is a `<3`-separation

        OUTPUT:

        boolean, or a tuple ``(boolean, frozenset)``

        ALGORITHM:

        The shifting algorithm

        EXAMPLES::

            sage: # needs sage.graphs
            sage: matroids.Uniform(2, 3)._is_3connected_shifting()
            True
            sage: M = Matroid(ring=QQ, matrix=[[1, 0, 0, 1, 1, 0],
            ....:                              [0, 1, 0, 1, 2, 0],
            ....:                              [0, 0, 1, 0, 0, 1]])
            sage: M._is_3connected_shifting()
            False
            sage: N = Matroid(circuit_closures={2: ['abc', 'cdef'],
            ....:                               3: ['abcdef']},
            ....:             groundset='abcdef')
            sage: N._is_3connected_shifting()
            False
            sage: matroids.catalog.BetsyRoss()._is_3connected_shifting()
            True
            sage: M = matroids.catalog.R6()
            sage: M._is_3connected_shifting()                                           # needs sage.rings.finite_rings
            False
            sage: B, X = M._is_3connected_shifting(True)                                # needs sage.rings.finite_rings
            sage: M.connectivity(X) < 3                                                 # needs sage.rings.finite_rings
            True
        """
        if not self.is_connected():
            if certificate:
                return False, self.components()[0]
            else:
                return False
        if self.rank() > self.size() - self.rank():
            return self.dual()._is_3connected_shifting(certificate)

        # the partial matrix
        M2 = self._reduced_representation()
        M = M2._matrix_()
        X, Y = self._current_rows_cols()

        # create mapping between elements and columns
        dX = dict(zip(range(len(X)), X))
        dY = dict(zip(range(len(Y)), Y))

        for (x, y) in spanning_forest(M):
            P_rows=[x]
            P_cols=[y]
            Q_rows=[]
            Q_cols=[]
            sol,cert_pair = M2.shifting_all(P_rows, P_cols, Q_rows, Q_cols, 2)
            if sol:
                if certificate:
                    cert = set([])
                    for x in cert_pair[0]:
                        cert.add(dX[x])
                    for y in cert_pair[1]:
                        cert.add(dY[y])
                    return False, cert
                return False
        if certificate:
            return True, None
        return True

    cpdef _is_4connected_shifting(self, certificate=False):
        r"""
        Return ``True`` if the matroid is 4-connected, ``False`` otherwise. It can
        optionally return a separator as a witness.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); if ``True``,
          then return ``True, None`` if the matroid is 4-connected,
          and ``False,`` `X` otherwise, where `X` is a `<4`-separation

        OUTPUT:

        boolean, or a tuple ``(boolean, frozenset)``

        ALGORITHM:

        The shifting algorithm

        EXAMPLES::

            sage: M = matroids.Uniform(2, 6)
            sage: B, X = M._is_4connected_shifting(True)                                # needs sage.graphs
            sage: (B, M.connectivity(X)<=3)                                             # needs sage.graphs
            (False, True)
            sage: matroids.Uniform(4, 8)._is_4connected_shifting()                      # needs sage.graphs
            True
            sage: M = Matroid(field=GF(2), matrix=[[1,0,0,1,0,1,1,0,0,1,1,1],
            ....:                                  [0,1,0,1,0,1,0,1,0,0,0,1],
            ....:                                  [0,0,1,1,0,0,1,1,0,1,0,1],
            ....:                                  [0,0,0,0,1,1,1,1,0,0,1,1],
            ....:                                  [0,0,0,0,0,0,0,0,1,1,1,1]])
            sage: M._is_4connected_shifting()                                           # needs sage.graphs
            True
        """
        if self.rank()>self.size()-self.rank():
            return self.dual()._is_4connected_shifting(certificate)
        if not self._is_3connected_shifting():
            return self._is_3connected_shifting(certificate)

        # the partial matrix
        M2 = self._reduced_representation()
        M = M2._matrix_()
        X, Y = self._current_rows_cols()

        dX = dict(zip(range(len(X)),X))
        dY = dict(zip(range(len(Y)),Y))
        n = len(X)
        m = len(Y)

        sol = False
        T = spanning_stars(M)

        for (x1, y1) in T:
            # The whiting out
            B = copy(M)
            for (x, y) in product(range(n), range(m)):
                if (x1 != x and y1 != y):
                    if (M[x1, y] != 0 and M[x, y1] != 0 and M[x, y] != 0):
                        B[x, y] = 0

            # remove row x1 and y1
            Xp = list(range(n))
            Xp.remove(x1)
            Yp = list(range(m))
            Yp.remove(y1)

            B = B.matrix_from_rows_and_columns(Xp,Yp)

            # produce a spanning forest of B
            for (x,y) in spanning_forest(B):
                if x >= x1:
                    x = x + 1
                if y >= y1:
                    y = y + 1
                # rank 2 matrix and rank 0 matrix
                P_rows = [x,x1]
                P_cols = [y,y1]
                Q_rows = []
                Q_cols = []
                # make sure the matrix has rank 2
                if M.matrix_from_rows_and_columns(P_rows,P_cols).rank() == 2:
                    sol,cert_pair = M2.shifting_all(P_rows, P_cols, Q_rows, Q_cols, 3)
                    if sol:
                        break
                # rank 1 matrix and rank 1 matrix
                P_rows = [x1]
                P_cols = [y1]
                Q_rows = [x]
                Q_cols = [y]
                # both matrix have rank 1
                sol,cert_pair = M2.shifting_all(P_rows, P_cols, Q_rows, Q_cols, 3)
                if sol:
                    break
            if sol:
                if certificate:
                    (certX, certY) = cert_pair
                    cert = set([])
                    for x in certX:
                        cert.add(dX[x])
                    for y in certY:
                        cert.add(dY[y])
                    return False, cert
                return False
        if certificate:
            return True, None
        return True

    def orlik_terao_algebra(self, R=None, ordering=None, **kwargs):
        """
        Return the Orlik-Terao algebra of ``self``.

        INPUT:

        - ``R`` -- (default: the base ring of ``self``) the base ring
        - ``ordering`` -- (optional) an ordering of the groundset

        .. SEEALSO::

            :class:`~sage.algebras.orlik_terao.OrlikTeraoAlgebra`

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_terao_algebra(); OS
            Orlik-Terao algebra of Wheel(3):
             Regular matroid of rank 3 on 6 elements with 16 bases
             over Integer Ring
            sage: OS.base_ring()
            Integer Ring
            sage: M.orlik_terao_algebra(QQ).base_ring()
            Rational Field

            sage: G = SymmetricGroup(3);                                                # needs sage.groups
            sage: OTG = M.orlik_terao_algebra(QQ, invariant=G)                          # needs sage.groups

            sage: # needs sage.groups
            sage: G = SymmetricGroup(4)
            sage: action = lambda g, x: g(x + 1) - 1
            sage: OTG1 = M.orlik_terao_algebra(QQ, invariant=(G, action))
            sage: OTG2 = M.orlik_terao_algebra(QQ, invariant=(action, G))
            sage: OTG1 is OTG2
            True
        """
        if R is None:
            R = self.base_ring()

        if 'invariant' in kwargs:
            G_action = kwargs.pop('invariant')

            from sage.categories.semigroups import Semigroups
            if len(G_action) > 1 and G_action not in Semigroups:
                G, action = G_action
                if action in Semigroups:
                    G, action = action, G
            else:
                G, action = G_action, None # the None action is g.__call__

            from sage.algebras.orlik_terao import OrlikTeraoInvariantAlgebra

            return OrlikTeraoInvariantAlgebra(R, self, G,
                                              action_on_groundset=action,
                                              ordering=ordering,
                                              **kwargs)

        from sage.algebras.orlik_terao import OrlikTeraoAlgebra
        return OrlikTeraoAlgebra(R, self, ordering)

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle_lean_linear_matroid, (version, data))``, where
        ``unpickle_lean_linear_matroid`` is the name of a function that, when
        called with ``(version, data)``, produces a matroid isomorphic to
        ``self``. ``version`` is an integer (currently 0) and ``data`` is a
        tuple ``(A, E, reduced, name)`` where ``A`` is the representation
        matrix, ``E`` is the groundset of the matroid, ``reduced`` is a
        boolean indicating whether ``A`` is a reduced matrix, and ``name`` is
        a custom name.

        .. WARNING::

            Users should never call this function directly.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1], [0, 1, 0, 1, 2],
            ....:                                           [0, 0, 1, 1, 3]]))
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('U35')
            sage: loads(dumps(M))
            U35
            sage: M = Matroid(Matrix(GF(7), [[1, 0, 1], [1, 0, 1]]))
            sage: N = loads(dumps(M))
            sage: N.representation()
            [1 0 1]
            [1 0 1]
        """
        import sage.matroids.unpickling
        cdef LeanMatrix A
        version = 0
        if self._representation is not None:
            A = self._representation
            gs = self._E
            reduced = False
        else:
            A = self._reduced_representation()
            rows, cols = self._current_rows_cols()
            gs = rows + cols
            reduced = True
        data = (A, gs, reduced, self.get_custom_name())
        return sage.matroids.unpickling.unpickle_linear_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: N = M.relabel({'g': 'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'x']

        TESTS::

            sage: M = matroids.catalog.Fano()
            sage: f = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset_list()]
        M = LinearMatroid(groundset=E, matrix=self._matrix_())
        return M

# Binary matroid

cdef class BinaryMatroid(LinearMatroid):
    r"""
    Binary matroids.

    A binary matroid is a linear matroid represented over the finite field
    with two elements. See :class:`LinearMatroid` for a definition.

    The simplest way to create a ``BinaryMatroid`` is by giving only a matrix
    `A`. Then, the groundset defaults to ``range(A.ncols())``. Any iterable
    object `E` can be given as a groundset. If `E` is a list, then ``E[i]``
    will label the `i`-th column of `A`. Another possibility is to specify a
    *reduced* matrix `B`, to create the matroid induced by `A = [ I B ]`.

    INPUT:

    - ``matrix`` -- (default: ``None``) a matrix whose column vectors
      represent the matroid.
    - ``reduced_matrix`` -- (default: ``None``) a matrix `B` such that
      `[I\ \ B]` represents the matroid, where `I` is an identity matrix with
      the same number of rows as `B`. Only one of ``matrix`` and
      ``reduced_matrix`` should be provided.
    - ``groundset`` -- (default: ``None``) an iterable containing the element
      labels. When provided, must have the correct number of elements: the
      number of columns of ``matrix`` or the number of rows plus the number
      of columns of ``reduced_matrix``.
    - ``ring`` -- (default: ``None``) ignored
    - ``keep_initial_representation`` -- boolean (default: ``True``); whether
      or not an internal copy of the input matrix should be preserved. This
      can help to see the structure of the matroid (e.g. in the case of
      graphic matroids), and makes it easier to look at extensions. However,
      the input matrix may have redundant rows, and sometimes it is desirable
      to store only a row-reduced copy.
    - ``basis`` -- (default: ``None``) when provided, this is an ordered
      subset of ``groundset``, such that the submatrix of ``matrix`` indexed
      by ``basis`` is an identity matrix. In this case, no row reduction takes
      place in the initialization phase.

    OUTPUT: a :class:`BinaryMatroid` instance based on the data above

    .. NOTE::

        An indirect way to generate a binary matroid is through the
        :func:`Matroid() <sage.matroids.constructor.Matroid>` function. This
        is usually the preferred way, since it automatically chooses between
        :class:`BinaryMatroid` and other classes. For direct access to the
        ``BinaryMatroid`` constructor, run::

            sage: from sage.matroids.advanced import *

    EXAMPLES::

        sage: A = Matrix(GF(2), 2, 4, [[1, 0, 1, 1], [0, 1, 1, 1]])
        sage: M = Matroid(A)
        sage: M
        Binary matroid of rank 2 on 4 elements, type (0, 6)
        sage: sorted(M.groundset())
        [0, 1, 2, 3]
        sage: Matrix(M)
        [1 0 1 1]
        [0 1 1 1]
        sage: M = Matroid(matrix=A, groundset='abcd')
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd']
        sage: B = Matrix(GF(2), 2, 2, [[1, 1], [1, 1]])
        sage: N = Matroid(reduced_matrix=B, groundset='abcd')
        sage: M == N
        True
    """
    def __init__(self, matrix=None, groundset=None, reduced_matrix=None, ring=None, keep_initial_representation=True, basis=None):
        """
        See class definition for full documentation.

        .. NOTE::

            The extra argument ``basis``, when provided, is an ordered list of
            elements of the groundset, ordered such that they index a standard
            identity matrix within ``matrix``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: BinaryMatroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],      # indirect doctest
            ....:                                     [0, 1, 1, 2, 3]]))
            Binary matroid of rank 2 on 5 elements, type (1, 7)
        """
        cdef BinaryMatrix A
        cdef long r, c
        cdef list P
        global GF2, GF2_zero, GF2_one, GF2_not_defined
        if GF2_not_defined:
            GF2 = GF(2)
            GF2_zero = GF2(0)
            GF2_one = GF2(1)
            GF2_not_defined = False

        # Setup representation; construct displayed basis
        if matrix is not None:
            A = BinaryMatrix(matrix.nrows(), matrix.ncols(), M=matrix)
            if keep_initial_representation:
                self._representation = A.copy()   # Deprecated Sage matrix operation
            if basis is None:
                P = gauss_jordan_reduce(A, range(A.ncols()))
                A.resize(len(P))   # Not a Sage matrix operation
            self._A = A
        else:
            A = BinaryMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix)
            P = list(range(A.nrows()))
            self._A = A.prepend_identity()   # Not a Sage matrix operation

        # Setup groundset, BasisExchangeMatroid data
        if groundset is None:
            groundset = list(range(self._A.ncols()))
        else:
            if len(groundset) != self._A.ncols():
                raise ValueError("size of groundset does not match size of matrix")
        if basis is None:
            bas = [groundset[i] for i in P]
        else:
            bas = basis
        BasisExchangeMatroid.__init__(self, groundset, bas)

        # Setup index of displayed basis
        self._prow = <long* > sig_malloc((self._A.ncols()) * sizeof(long))
        for c in range(self._A.ncols()):
            self._prow[c] = -1
        if matrix is not None:
            if basis is None:
                for r in range(len(P)):
                    self._prow[P[r]] = r
            else:
                for r in range(self._A.nrows()):
                    self._prow[self._idx[basis[r]]] = r
        else:
            for r in range(self._A.nrows()):
                self._prow[r] = r

        self._zero = GF2_zero
        self._one = GF2_one

    cpdef base_ring(self):
        r"""
        Return the base ring of the matrix representing the matroid,
        in this case `\GF{2}`.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.base_ring()
            Finite Field of size 2
        """
        global GF2
        return GF2

    cpdef characteristic(self):
        """
        Return the characteristic of the base ring of the matrix representing
        the matroid, in this case `2`.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.characteristic()
            2
        """
        return 2

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check if ``self.basis() - x + y`` is again a basis. Internal method.
        """
        return (<BinaryMatrix>self._A).is_nonzero(self._prow[x], y)

    cdef int _exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        cdef long p = self._prow[x]
        self._A.pivot(p, y)   # Not a Sage matrix operation
        self._prow[y] = p
        BasisExchangeMatroid._exchange(self, x, y)

    cdef  __fundamental_cocircuit(self, bitset_t C, long x):
        r"""
        Fill bitset `C` with the incidence vector of the `B`-fundamental cocircuit using ``x``. Internal method using packed elements.
        """
        bitset_copy(C, (<BinaryMatrix>self._A)._M[self._prow[x]])

    cdef _coclosure_internal(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``coclosure``.

        This function overrides the internal function BasisExchangeMatroid._coclosure_internal() of the parent class.
        The implementation should be more efficient for BinaryMatroid, due to the fact that in this class,
        __fundamental_cocircuit is much faster than __fundamental_circuit.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        self._move(self._inside, self._outside)

        bitset_copy(R, F)
        bitset_difference(self._inside, self._current_basis, F)
        cdef long y = bitset_first(self._inside)
        while y >= 0:
            self.__fundamental_cocircuit(self._outside, y)
            bitset_discard(self._outside, y)
            if bitset_issubset(self._outside, F):
                bitset_add(R, y)
            y = bitset_next(self._inside, y + 1)

    cdef  _exchange_value_internal(self, long x, long y):
        r"""
        Return the (x, y) entry of the current representation.
        """
        if (<BinaryMatrix>self._A).get(self._prow[x], y):   # Not a Sage matrix operation
            return self._one
        else:
            return self._zero

    def _repr_(self):
        """
        Return a string representation of ``self``.

        The type consists of :meth:`BinaryMatroid.bicycle_dimension` and
        :meth:`BinaryMatroid.brown_invariant`.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.rename()
            sage: repr(M)  # indirect doctest
            'Binary matroid of rank 3 on 7 elements, type (3, 0)'
        """
        return f'Binary matroid of rank {self.rank()} on {self.size()} elements, type ({self.bicycle_dimension()}, {self.brown_invariant()})'

    cpdef _current_rows_cols(self, B=None):
        """
        Return the current row and column labels of a reduced matrix.

        INPUT:

        - ``B`` -- (default: ``None``) if provided, first find a basis having
          maximal intersection with ``B``

        OUTPUT:

        - ``R`` -- list of row indices; corresponds to the currently used
          internal basis
        - ``C`` -- list of column indices; corresponds to the complement of
          the current internal basis

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: A = M._reduced_representation('efg')
            sage: R, C = M._current_rows_cols()
            sage: (sorted(R), sorted(C))
            (['e', 'f', 'g'], ['a', 'b', 'c', 'd'])
            sage: R, C = M._current_rows_cols(B='abg')
            sage: (sorted(R), sorted(C))
            (['a', 'b', 'g'], ['c', 'd', 'e', 'f'])
        """
        if B is not None:
            self._move_current_basis(B, set())
        basis = self.basis()
        rows = [0] * self.full_rank()
        cols = [0] * self.full_corank()
        c = 0
        for e in self._E:
            if e in basis:
                rows[self._prow[self._idx[e]]] = e
            else:
                cols[c] = e
                c += 1
        return rows, cols

    cpdef LeanMatrix _basic_representation(self, B=None):
        """
        Return a basic matrix representation of the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `M` representing the matroid, where `M[B'] = I` for a basis
        `B'` that maximally intersects the given set `B`. If not provided, the
        current basis used internally is chosen for `B'`. For a stable
        representation, use ``self.representation()``.

        .. NOTE::

            The method self.groundset_list() gives the labelling of the
            columns by the elements of the matroid. The matrix returned
            is a LeanMatrix subclass, which is intended for internal use
            only. Use the ``representation()`` method to get a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M._basic_representation()
            3 x 7 BinaryMatrix
            [1000111]
            [0101011]
            [0011101]
            sage: matrix(M._basic_representation('efg'))
            [1 1 0 1 1 0 0]
            [1 0 1 1 0 1 0]
            [1 1 1 0 0 0 1]
        """
        if B is not None:
            self._move_current_basis(B, set())
        return self._A.copy()   # Deprecated Sage matrix operation

    cpdef LeanMatrix _reduced_representation(self, B=None):
        r"""
        Return a reduced representation of the matroid, i.e. a matrix `R` such
        that `[I\ \ R]` represents the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `R` forming a reduced representation of the matroid, with
        rows labeled by a basis `B'` that maximally intersects the given set
        `B`. If not provided, the current basis used internally labels the
        rows.

        .. NOTE::

            The matrix returned is a LeanMatrix subclass, which is intended
            for internal use only. Use the ``representation()`` method to get
            a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M._reduced_representation()
            3 x 4 BinaryMatrix
            [0111]
            [1011]
            [1101]
            sage: matrix(M._reduced_representation('efg'))
            [1 1 0 1]
            [1 0 1 1]
            [1 1 1 0]
        """
        if B is not None:
            self._move_current_basis(B, set())
        _, cols = self._current_rows_cols()
        return self._A.matrix_from_rows_and_columns(range(self.full_rank()), [self._idx[e] for e in cols])

    # isomorphism

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate = True``, a dictionary giving
        the isomorphism or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: M1 = matroids.catalog.Fano()
            sage: M2 = Matroid(ring=GF(2),
            ....:   reduced_matrix=[[1, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
            sage: M1._is_isomorphic(M2)
            True
            sage: M1._is_isomorphic(M2, certificate=True)
            (True, {'a': 0, 'b': 1, 'c': 2, 'd': 4, 'e': 3, 'f': 5, 'g': 6})

            sage: M1 = matroids.catalog.Fano().delete('a')
            sage: M2 = matroids.Whirl(3)
            sage: M1._is_isomorphic(M2)
            False
            sage: M1._is_isomorphic(M2, certificate=True)
            (False, None)
            sage: M1._is_isomorphic(matroids.Wheel(3))
            True
            sage: M1._is_isomorphic(matroids.Wheel(3), certificate=True)
            (True, {'b': 1, 'c': 2, 'd': 4, 'e': 3, 'f': 5, 'g': 0})
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        if isinstance(other, BinaryMatroid):
            return self.is_field_isomorphic(other)
        else:
            return LinearMatroid._is_isomorphic(self, other)

    cpdef _is_isomorphism(self, other, morphism):
        r"""
        Test if a given bijection is an isomorphism.

        Version of ``is_isomorphism()`` that does no type checking
        of ``morphism``.

        INPUT:

        - ``other`` -- matroid
        - ``morphism`` -- dictionary mapping the groundset of ``self`` to
          the groundset of ``other``

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.catalog.Fano().delete(['a'])
            sage: N = matroids.catalog.Fano().delete(['b'])
            sage: morphism = {'b':'a', 'c':'c', 'd':'e', 'e':'d', 'f':'f', 'g':'g'}
            sage: M._is_isomorphism(N, morphism)
            True
        """
        if isinstance(other, BinaryMatroid):
            return self.is_field_isomorphism(other, morphism)
        else:
            return LinearMatroid._is_isomorphism(self, other, morphism)

    # invariants
    cpdef _make_invariant(self):
        """
        Create an invariant.

        Internal method; see ``_invariant`` for full documentation.

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M._invariant()  # indirect doctest
            (3, 0, 7, 0, 0, 0, 0, 0)
        """
        cdef BinaryMatrix B
        cdef long r, d, i, j
        if self._b_invariant is not None:
            return
        B = (<BinaryMatrix>self._A).copy()   # Deprecated Sage matrix operation
        r = B.nrows()
        b = 0
        d = 0
        i = 0
        U = set()
        while i + d < r:
            for j in range(i, r - d):
                if B.row_len(j) % 2 == 1:   # Not a Sage matrix operation
                    B.swap_rows_c(i, j)
                    break
            if B.row_len(i) % 2 == 1:   # Not a Sage matrix operation
                for j in range(i + 1, r - d):
                    if B.row_inner_product(i, j):   # Not a Sage matrix operation
                        B.add_multiple_of_row_c(j, i, 1, 0)
                if B.row_len(i) % 4 == 1:   # Not a Sage matrix operation
                    b += 1
                else:
                    b -= 1
                U.add(i)
                i += 1
            else:
                for j in range(i + 1, r - d):
                    if B.row_inner_product(i, j):   # Not a Sage matrix operation
                        B.swap_rows_c(i + 1, j)
                        break
                if i + 1 < r - d:
                    if B.row_inner_product(i, i + 1):   # Not a Sage matrix operation
                        for j in range(i + 2, r):
                            if B.row_inner_product(i, j):   # Not a Sage matrix operation
                                B.add_multiple_of_row_c(j, i + 1, 1, 0)
                            if B.row_inner_product(i + 1, j):   # Not a Sage matrix operation
                                B.add_multiple_of_row_c(j, i, 1, 0)
                        if B.row_len(i) % 4 == 2 and B.row_len(i + 1) % 4 == 2:   # Not a Sage matrix operation
                            b += 4
                        i += 2
                    else:
                        d += 1
                        B.swap_rows_c(i, r - d)
                else:
                    d += 1

        doubly_even = True
        for i in range(r - d, r):
            if B.row_len(i) % 4 == 2:   # Not a Sage matrix operation
                doubly_even = False
                break
        if doubly_even:
            b2 = b % 8
        else:
            b2 = None

        Fm = B.row_union(range(r - d, r))   # Not a Sage matrix operation
        Fp = [i for i in B.row_sum(U) if i not in Fm]   # Not a Sage matrix operation
        F0 = [i for i in range(len(self)) if i not in (Fm + Fp)]

        BT = B.transpose()
        self._b_projection = BT._matrix_times_matrix_((B._matrix_times_matrix_(BT))._matrix_times_matrix_(B))
        P = [F0, Fp]
        p = []
        for a in range(2):
            for b in range(a + 1):
                x = 0
                for i in P[a]:
                    for j in P[b]:
                        if self._b_projection.get(i, j) != 0:   # Not a Sage matrix operation
                            x += 1
                p.append(x)
        if d > 0:
            F = F0 + Fp
            self._b_projection = self._b_projection.matrix_from_rows_and_columns(F, F)
        self._b_invariant = tuple([d, b2, len(Fm), len(F0), len(Fp), p[0], p[1], p[2]])
        self._b_partition = tuple([Fm, F0, Fp])

    cpdef _invariant(self):
        r"""
        Return a matroid invariant.

        See [Pen2012]_ for more information.

        OUTPUT:

        A tuple ``(d, b, Lm, L0, Lp, p0, p1, p2)``, with the following
        interpretation:

        - ``d`` -- the :meth:`bicycle dimension <BinaryMatroid.bicycle_dimension>`
        - ``b`` -- the :meth:`Brown invariant <BinaryMatroid.brown_invariant>`
        - ``(Lm, L0, Lp)`` is the triple of lengths of the principal tripartition.
        - ``(p0, p1, p2)`` are the counts of edges in a characteristic graph
          of the matroid, whose vertices are the union of ``F_-`` and ``F_0``
          from the principal tripartition.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BinaryMatroid(matroids.AG(2, 5).representation())
            sage: M._invariant()
            (2, 1, 24, 0, 1, 0, 0, 1)
        """
        if self._b_invariant is None:
            self._make_invariant()
        return self._b_invariant

    cpdef bicycle_dimension(self):
        r"""
        Return the bicycle dimension of the binary matroid.

        The *bicycle dimension* of a linear subspace `V` is
        `\dim(V\cap V^\perp)`. The bicycle dimension of a matroid equals the
        bicycle dimension of its cocycle-space, and is an invariant for binary
        matroids. See [Pen2012]_, [GR2001]_ for more information.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.bicycle_dimension()
            3
        """
        if self._b_invariant is None:
            self._make_invariant()
        return self._b_invariant[0]

    cpdef brown_invariant(self):
        r"""
        Return the value of Brown's invariant for the binary matroid.

        For a binary space `V`, consider the sum
        `B(V):=\sum_{v\in V} i^{|v|}`, where `|v|` denotes the number of
        nonzero entries of a binary vector `v`. The value of the Tutte
        Polynomial in the point `(-i, i)` can be expressed in terms of
        `B(V)`, see [Pen2012]_. If `|v|` equals `2` modulo 4 for some
        `v\in V\cap V^\perp`, then `B(V)=0`. In this case, Browns invariant is
        not defined. Otherwise, `B(V)=\sqrt{2}^k \exp(\sigma \pi i/4)` for
        some integers `k, \sigma`. In that case, `k` equals the bycycle
        dimension of `V`, and Browns invariant for `V` is defined as `\sigma`
        modulo `8`.

        The Brown invariant of a binary matroid equals the Brown invariant of
        its cocycle-space.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: M.brown_invariant()
            0
            sage: M = Matroid(Matrix(GF(2), 3, 8, [[1, 0, 0, 1, 1, 1, 1, 1],
            ....:                                  [0, 1, 0, 1, 1, 0, 0, 0],
            ....:                                  [0, 0, 1, 0, 0, 1, 1, 0]]))
            sage: M.brown_invariant() is None
            True
        """
        if self._b_invariant is None:
            self._make_invariant()
        return self._b_invariant[1]

    cpdef _principal_tripartition(self):
        r"""
        Return the principal tripartition of the binary matroid.

        The principal tripartition is a partition `(F_{-1}, F_0, F_{1})` of
        the groundset. A defining property is the following. It is
        straightforward that if the bicycle dimension of a matroid `M` is `k`,
        then the bicycle dimension of `M\setminus e' is one of `k-1, k, k + 1`
        for each element `e` of `M`. Then if `F_i` denotes the set of elements
        such that the bicycle dimension of `M\setminus e` is `k + i`, we
        obtain the principal tripartition `(F_{-1}, F_0, F_{1})` of `M`.
        See [Pen2012]_ and [GR2001]_.

        OUTPUT:

        ``(F_{-1}, F_0, F_{1})``, the principal tripartition of the matroid.

        EXAMPLES::

            sage: M = matroids.catalog.S8()
            sage: for F in M._principal_tripartition(): print(sorted(F))
            ['a', 'b', 'c', 'e', 'f', 'g']
            ['d']
            ['h']
            sage: M.bicycle_dimension()
            2
            sage: for i in [-1, 0, 1]:
            ....:     print(sorted(e for e in M.groundset()
            ....:                    if (M.delete(e)).bicycle_dimension() == 2 + i))
            ['a', 'b', 'c', 'e', 'f', 'g']
            ['d']
            ['h']
        """
        if self._b_invariant is None:
            self._make_invariant()
        P = self._b_partition
        return frozenset([self._E[e] for e in P[0]]), frozenset([self._E[e] for e in P[1]]), frozenset([self._E[e] for e in P[2]])

    cpdef BinaryMatrix _projection(self):
        """
        Return the projection matrix onto the row space.

        This projection is determined modulo the bicycle space. See [Pen2012]_.

        INPUT:

        - Nothing

        OUTPUT:

        A binary matrix `P`, so that the `e`-th column of `P` is the
        incidence vector of a cocycle `C` such that `C-e` is a cycle. Such a
        `C` is determined up to taking the symmetric difference with bicycles.
        We output the restriction of `P` to rows and columns that are not in
        any bicycle.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = BinaryMatroid(matrix(matroids.catalog.R12()))
            sage: M._projection()
            12 x 12 BinaryMatrix
            [001110111000]
            [001101110100]
            [111011100010]
            [110111010001]
            [101100001011]
            [011100000111]
            [111000101110]
            [110100011101]
            [100010110011]
            [010001110011]
            [001011101110]
            [000111011101]
        """
        if self._b_invariant is None:
            self._make_invariant()
        return self._b_projection

    cpdef BinaryMatrix _projection_partition(self):
        """
        Return the equitable partition of the graph whose incidence matrix is
        the projection matrix of this matroid.

        See method ``._projection()``.

        INPUT:

        - Nothing

        OUTPUT: an ordered partition

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.R12()
        sage: N = BinaryMatroid(reduced_matrix=M.representation(reduced=True,
        ....:                         labels=False), groundset='abcdefghijkl')
        sage: Npp = N._projection_partition(); Npp  # random
        2 x 12 BinaryMatrix
        [110011001100]
        [001100110011]
        sage: sorted(Npp._matrix_().rows())                                             # needs sage.rings.finite_rings
        [(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0), (0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)]
        """
        if self._eq_part is None:
            if self._b_invariant is None:
                self._make_invariant()
            self._eq_part = self._b_projection.equitable_partition()   # Not a Sage matrix operation
        return self._eq_part

    cpdef _fast_isom_test(self, other):
        r"""
        Run a quick test to see if two binary matroids are isomorphic.

        The test is based on comparing strong invariants. See [Pen2012]_ for a
        full account of these invariants.

        INPUT:

        - ``other`` -- a binary matroid

        OUTPUT:

        - ``True``, if ``self`` is isomorphic to ``other``;
        - ``False``, if ``self`` is not isomorphic to ``other``;
        - ``None``, if this test is inconclusive

        EXAMPLES::

           sage: M = matroids.catalog.S8()
           sage: N = matroids.catalog.S8()
           sage: M._fast_isom_test(N) is None
           True
        """
        if self._invariant() != other._invariant():
            return False
        q = self._projection().is_isomorphic(other._projection(), self._projection_partition(), other._projection_partition())   # Not a Sage matrix operation
        if self.bicycle_dimension() == 0:
            return q
        if not q:
            return False

    # minors, dual

    cpdef _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``
        - ``deletions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))
            Binary matroid of rank 2 on 4 elements, type (0, 6)
        """
        self._move_current_basis(contractions, deletions)
        bas = list(self.basis() - contractions)
        R = [self._prow[self._idx[b]] for b in bas]
        F = self.groundset() - (deletions | contractions)
        C = [self._idx[f] for f in F]
        A, C2 = (<BinaryMatrix>self._A).matrix_from_rows_and_columns_reordered(R, C)
        return BinaryMatroid(matrix=A,
                             groundset=[self._E[c] for c in C2],
                             basis=bas,
                             keep_initial_representation=False)

    # graphicness test
    cpdef bint is_graphic(self) noexcept:
        """
        Test if the binary matroid is graphic.

        A matroid is *graphic* if there exists a graph whose edge set equals
        the groundset of the matroid, such that a subset of elements of the
        matroid is independent if and only if the corresponding subgraph is
        acyclic.

        OUTPUT: boolean

        EXAMPLES::

            sage: R10 = matroids.catalog.R10()
            sage: M = Matroid(ring=GF(2), reduced_matrix=R10.representation(
            ....:                                 reduced=True, labels=False))
            sage: M.is_graphic()
            False
            sage: K5 = Matroid(graphs.CompleteGraph(5), regular=True)                   # needs sage.graphs
            sage: M = Matroid(ring=GF(2), reduced_matrix=K5.representation(             # needs sage.graphs sage.rings.finite_rings
            ....:                                 reduced=True, labels=False))
            sage: M.is_graphic()                                                        # needs sage.graphs sage.rings.finite_rings
            True
            sage: M.dual().is_graphic()                                                 # needs sage.graphs
            False

        ALGORITHM:

        In a recent paper, Geelen and Gerards [GG2012]_ reduced the problem to
        testing if a system of linear equations has a solution. While not the
        fastest method, and not necessarily constructive (in the presence of
        2-separations especially), it is easy to implement.
        """
        global GF2
        cdef int r, c
        B= self.basis()
        C = [self._fundamental_cocircuit(B,e) for e in B]

        c = 1
        col = {}
        for e in range(len(B)):
            for f in range(len(B)):
                if e is not f:
                    col[e, f] = c
                    c += 1
        M = []
        r = 0
        for e in range(len(B)):
            for f in range(e):
                for g in range(f):
                    if not C[e].issuperset(C[f] & C[g]):
                        M.append([col[e,f], col[e,g]])
                        r += 1
                    if not C[f].issuperset(C[e] & C[g]):
                        M.append([col[f,e], col[f,g]])
                        r += 1
                    if not C[g].issuperset(C[e] & C[f]):
                        M.append([col[g,e], col[g,f]])
                        r += 1
                    if len(C[e] & C[f] & C[g]) > 0:
                        M.append([0,col[e,f], col[e,g], col[f,e], col[f,g], col[g,e], col[g,f]])
                        r += 1
        cdef BinaryMatrix m = BinaryMatrix(r,c)
        for r in range(len(M)):
            for c in M[r]:
                m.set(r,c)
        # now self is graphic iff there is a binary vector x so that M*x = 0 and x_0 = 1, so:
        return BinaryMatroid(m).corank(frozenset([0])) > 0

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data obey the matroid axioms.

        Since this is a linear matroid over the field `\GF{2}`, this is always
        the case.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: ``True``, or ``(True, {})``

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(2), [[]]))
            sage: M.is_valid()
            True
        """
        return True if not certificate else (True, {})

    # representability

    cpdef binary_matroid(self, randomized_tests=1, verify = True):
        r"""
        Return a binary matroid representing ``self``.

        INPUT:

        - ``randomized_tests`` -- ignored
        - ``verify`` -- ignored

        OUTPUT: a binary matroid

        ALGORITHM:

        ``self`` is a binary matroid, so just return ``self``.

        .. SEEALSO::

            :meth:`M.binary_matroid()
            <sage.matroids.matroid.Matroid.binary_matroid>`

        EXAMPLES::

            sage: N = matroids.catalog.Fano()
            sage: N.binary_matroid() is N
            True
        """
        return self

    cpdef is_binary(self, randomized_tests=1):
        r"""
        Decide if ``self`` is a binary matroid.

        INPUT:

        - ``randomized_tests`` -- ignored

        OUTPUT: boolean

        ALGORITHM:

        ``self`` is a binary matroid, so just return ``True``.

        .. SEEALSO::

            :meth:`M.is_binary() <sage.matroids.matroid.Matroid.is_binary>`

        EXAMPLES::

            sage: N = matroids.catalog.Fano()
            sage: N.is_binary()
            True
        """
        return True

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle_binary_matroid, (version, data))``, where
        ``unpickle_binary_matroid`` is the name of a function that, when
        called with ``(version, data)``, produces a matroid isomorphic to
        ``self``. ``version`` is an integer (currently 0) and ``data`` is a
        tuple ``(A, E, B, name)`` where ``A`` is the representation
        matrix, ``E`` is the groundset of the matroid, ``B`` is the currently
        displayed basis, and ``name`` is a custom name.

        .. WARNING::

            Users should never call this function directly.

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 0, 1],
            ....:        [0, 0, 1, 1]]))
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('U34')
            sage: loads(dumps(M))
            U34
            sage: M = Matroid(Matrix(GF(2), [[1, 0, 1], [1, 0, 1]]))
            sage: loads(dumps(M)).representation()
            [1 0 1]
            [1 0 1]

        TESTS:

        Check that :issue:`23437` is fixed::

            sage: M = matroids.catalog.Fano().dual()
            sage: B = list(M.bases())
            sage: N = loads(dumps(M))
            sage: N.closure(frozenset({'d'}))
            frozenset({'d'})
            sage: N.is_isomorphic(M)
            True
        """
        import sage.matroids.unpickling
        version = 0
        gs = self._E
        cdef list basis = [0] * self.full_rank()
        if self._representation is not None:
            A = self._representation
            basis = None
        else:
            A = self._A
            # current basis ordered so matrix cols form identity matrix:
            basis = self._current_rows_cols()[0]
        data = (A, gs, basis, self.get_custom_name())
        return sage.matroids.unpickling.unpickle_binary_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.Fano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: N = M.relabel({'g': 'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'x']

        TESTS::

            sage: M = matroids.catalog.Fano()
            sage: f = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset_list()]
        M = BinaryMatroid(groundset=E, matrix=self._matrix_())
        return M

cdef class TernaryMatroid(LinearMatroid):
    r"""
    Ternary matroids.

    A ternary matroid is a linear matroid represented over the finite field
    with three elements. See :class:`LinearMatroid` for a definition.

    The simplest way to create a ``TernaryMatroid`` is by giving only a
    matrix `A`. Then, the groundset defaults to ``range(A.ncols())``. Any
    iterable object `E` can be given as a groundset. If `E` is a list, then
    ``E[i]`` will label the `i`-th column of `A`. Another possibility is to
    specify a 'reduced' matrix `B`, to create the matroid induced by
    `A = [I\ \ B]`.

    INPUT:

    - ``matrix`` -- (default: ``None``) a matrix whose column vectors
      represent the matroid.
    - ``reduced_matrix`` -- (default: ``None``) a matrix `B` such that
      `[I\ \ B]` represents the matroid, where `I` is an identity matrix with
      the same number of rows as `B`. Only one of ``matrix`` and
      ``reduced_matrix`` should be provided.
    - ``groundset`` -- (default: ``None``) an iterable containing the element
      labels. When provided, must have the correct number of elements: the
      number of columns of ``matrix`` or the number of rows plus the number
      of columns of ``reduced_matrix``.
    - ``ring`` -- (default: ``None``) ignored
    - ``keep_initial_representation`` -- boolean (default: ``True``); decides
      whether or not an internal copy of the input matrix should be preserved.
      This can help to see the structure of the matroid (e.g. in the case of
      graphic matroids), and makes it easier to look at extensions. However,
      the input matrix may have redundant rows, and sometimes it is desirable
      to store only a row-reduced copy.
    - ``basis`` -- (default: ``None``) when provided, this is an ordered
      subset of ``groundset``, such that the submatrix of ``matrix`` indexed
      by ``basis`` is an identity matrix. In this case, no row reduction takes
      place in the initialization phase.

    OUTPUT: a ``TernaryMatroid`` instance based on the data above

    .. NOTE::

        The recommended way to generate a ternary matroid is through the
        :func:`Matroid() <sage.matroids.constructor.Matroid>` function. This
        is usually the preferred way, since it automatically chooses between
        ``TernaryMatroid`` and other classes. For direct access to the
        ``TernaryMatroid`` constructor, run::

            sage: from sage.matroids.advanced import *

    EXAMPLES::

        sage: A = Matrix(GF(3), 2, 4, [[1, 0, 1, 1], [0, 1, 1, 1]])
        sage: M = Matroid(A); M
        Ternary matroid of rank 2 on 4 elements, type 0-
        sage: sorted(M.groundset())
        [0, 1, 2, 3]
        sage: Matrix(M)
        [1 0 1 1]
        [0 1 1 1]
        sage: M = Matroid(matrix=A, groundset='abcd')
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd']
        sage: B = Matrix(GF(2), 2, 2, [[1, 1], [1, 1]])
        sage: N = Matroid(ring=GF(3), reduced_matrix=B, groundset='abcd')               # needs sage.rings.finite_rings
        sage: M == N                                                                    # needs sage.rings.finite_rings
        True
    """
    def __init__(self, matrix=None, groundset=None, reduced_matrix=None, ring=None, keep_initial_representation=True, basis=None):
        """
        See class definition for full documentation.

        .. NOTE::

            The extra argument ``basis``, when provided, is an ordered list of
            elements of the groundset, ordered such that they index a standard
            identity matrix within ``matrix``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: TernaryMatroid(matrix=Matrix(GF(5), [[1, 0, 1, 1, 1],       # indirect doctest
            ....:                                      [0, 1, 1, 2, 3]]))
            Ternary matroid of rank 2 on 5 elements, type 1+
        """
        cdef TernaryMatrix A
        cdef long r, c
        cdef list P
        global GF3, GF3_zero, GF3_one, GF3_minus_one, GF3_not_defined
        if GF3_not_defined:
            GF3 = GF(3)
            GF3_zero = GF3(0)
            GF3_one = GF3(1)
            GF3_minus_one = GF3(2)
            GF3_not_defined = False

        # Setup representation; construct displayed basis
        if matrix is not None:
            A = TernaryMatrix(matrix.nrows(), matrix.ncols(), M=matrix)
            if keep_initial_representation:
                self._representation = A.copy()   # Deprecated Sage matrix operation
            if basis is None:
                P = gauss_jordan_reduce(A, range(A.ncols()))
                A.resize(len(P))   # Not a Sage matrix operation
            self._A = A
        else:
            A = TernaryMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix)
            P = list(range(A.nrows()))
            self._A = A.prepend_identity()   # Not a Sage matrix operation

        # Setup groundset, BasisExchangeMatroid data
        if groundset is None:
            groundset = list(range(self._A.ncols()))
        else:
            if len(groundset) != self._A.ncols():
                raise ValueError("size of groundset does not match size of matrix")
        if basis is None:
            bas = [groundset[i] for i in P]
        else:
            bas = basis
        BasisExchangeMatroid.__init__(self, groundset, bas)

        # Setup index of displayed basis
        self._prow = <long* > sig_malloc((self._A.ncols()) * sizeof(long))
        for c in range(self._A.ncols()):
            self._prow[c] = -1
        if matrix is not None:
            if basis is None:
                for r in range(len(P)):
                    self._prow[P[r]] = r
            else:
                for r in range(self._A.nrows()):
                    self._prow[self._idx[basis[r]]] = r
        else:
            for r in range(self._A.nrows()):
                self._prow[r] = r

        self._zero = GF3_zero
        self._one = GF3_one
        self._two = GF3_minus_one

    cpdef base_ring(self):
        r"""
        Return the base ring of the matrix representing the matroid, in this
        case `\GF{3}`.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M.base_ring()
            Finite Field of size 3
        """
        global GF3
        return GF3

    cpdef characteristic(self):
        """
        Return the characteristic of the base ring of the matrix representing
        the matroid, in this case `3`.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M.characteristic()
            3
        """
        return 3

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check if ``self.basis() - x + y`` is again a basis. Internal method.
        """
        return (<TernaryMatrix>self._A).is_nonzero(self._prow[x], y)

    cdef int _exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        cdef long p = self._prow[x]
        self._A.pivot(p, y)   # Not a Sage matrix operation
        self._prow[y] = p
        BasisExchangeMatroid._exchange(self, x, y)

    cdef  __fundamental_cocircuit(self, bitset_t C, long x):
        r"""
        Fill bitset `C` with the incidence vector of the `B`-fundamental cocircuit using ``x``. Internal method using packed elements.
        """
        bitset_copy(C, (<TernaryMatrix>self._A)._M0[self._prow[x]])

    cdef _coclosure_internal(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``coclosure``.

        This function overrides the internal function BasisExchangeMatroid._coclosure_internal() of the parent class.
        The implementation should be more efficient for TernaryMatroid, due to the fact that in this class,
        __fundamental_cocircuit is much faster than __fundamental_circuit.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        self._move(self._inside, self._outside)

        bitset_copy(R, F)
        bitset_difference(self._inside, self._current_basis, F)
        cdef long y = bitset_first(self._inside)
        while y >= 0:
            self.__fundamental_cocircuit(self._outside, y)
            bitset_discard(self._outside, y)
            if bitset_issubset(self._outside, F):
                bitset_add(R, y)
            y = bitset_next(self._inside, y + 1)

    cdef  _exchange_value_internal(self, long x, long y):
        r"""
        Return the (x, y) entry of the current representation.
        """
        cdef long t = (<TernaryMatrix>self._A).get(self._prow[x], y)   # Not a Sage matrix operation
        if t == 0:
            return self._zero
        if t == 1:
            return self._one
        if t == -1:
            return self._two

    def _repr_(self):
        """
        Return a string representation of ``self``.

        The type consists of the ``bicycle_dimension`` and the ``character``.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M.rename()
            sage: repr(M)  # indirect doctest
            'Ternary matroid of rank 3 on 7 elements, type 0-'
        """
        S = f'Ternary matroid of rank {self.rank()} on {self.size()} elements, type {self.bicycle_dimension()}'
        if self.character() == 1:
            S = S + '+'
        else:
            S = S + '-'
        return S

    cpdef _current_rows_cols(self, B=None):
        """
        Return the current row and column labels of a reduced matrix.

        INPUT:

        - ``B`` -- (default: ``None``) if provided, first find a basis having
          maximal intersection with ``B``

        OUTPUT:

        - ``R`` -- list of row indices; corresponds to the currently used
          internal basis
        - ``C`` -- list of column indices; corresponds to the complement of
          the current internal basis

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: A = M._reduced_representation('efg')
            sage: R, C = M._current_rows_cols()
            sage: (sorted(R), sorted(C))
            (['e', 'f', 'g'], ['a', 'b', 'c', 'd'])
            sage: R, C = M._current_rows_cols(B='abg')
            sage: (sorted(R), sorted(C))
            (['a', 'b', 'g'], ['c', 'd', 'e', 'f'])
        """
        if B is not None:
            self._move_current_basis(B, set())
        basis = self.basis()
        rows = [0] * self.full_rank()
        cols = [0] * self.full_corank()
        c = 0
        for e in self._E:
            if e in basis:
                rows[self._prow[self._idx[e]]] = e
            else:
                cols[c] = e
                c += 1
        return rows, cols

    cpdef LeanMatrix _basic_representation(self, B=None):
        """
        Return a basic matrix representation of the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `M` representing the matroid, where `M[B'] = I` for a basis
        `B'` that maximally intersects the given set `B`. If not provided, the
        current basis used internally is chosen for `B'`. For a stable
        representation, use ``self.representation()``.

        .. NOTE::

            The method self.groundset_list() gives the labelling of the
            columns by the elements of the matroid. The matrix returned is a
            LeanMatrix subclass, which is intended for internal use only. Use
            the ``representation()`` method to get a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M._basic_representation()
            3 x 7 TernaryMatrix
            [+000+++]
            [0+0+0++]
            [00+++0+]
            sage: matrix(M._basic_representation('efg'))
            [1 2 0 2 1 0 0]
            [1 0 2 2 0 1 0]
            [2 1 1 2 0 0 1]
        """
        if B is not None:
            self._move_current_basis(B, set())
        return self._A.copy()   # Deprecated Sage matrix operation

    cpdef LeanMatrix _reduced_representation(self, B=None):
        r"""
        Return a reduced representation of the matroid, i.e. a matrix `R`
        such that `[I\ \ R]` represents the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `R` forming a reduced representation of the matroid, with
        rows labeled by a basis `B'` that maximally intersects the given set
        `B`. If not provided, the current basis used internally labels the
        rows.

        .. NOTE::

            The matrix returned is a LeanMatrix subclass, which is intended
            for internal use only. Use the ``representation()`` method to get
            a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M._reduced_representation()
            3 x 4 TernaryMatrix
            [0+++]
            [+0++]
            [++0+]
            sage: matrix(M._reduced_representation('efg'))
            [1 2 0 2]
            [1 0 2 2]
            [2 1 1 2]
        """
        if B is not None:
            self._move_current_basis(B, set())
        _, cols = self._current_rows_cols()
        return self._A.matrix_from_rows_and_columns(range(self.full_rank()), [self._idx[e] for e in cols])

    # isomorphism

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``. Internal version that
        performs no checks on input.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate = True``, a dictionary giving
        the isomorphism or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: M1 = matroids.catalog.NonFano().delete('a')
            sage: M2 = matroids.Whirl(3)
            sage: M1._is_isomorphic(M2)
            True

            sage: M2 = matroids.Wheel(3)
            sage: M1._is_isomorphic(M2)
            False
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        if isinstance(other, TernaryMatroid):
            return self.is_field_isomorphic(other)
        else:
            return LinearMatroid._is_isomorphic(self, other)

    # invariants

    cpdef _make_invariant(self):
        """
        Create an invariant.

        Internal method; see ``_invariant`` for full documentation.

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M._invariant()  # indirect doctest
            (0, 2, 0, 4, 3, 0, 12, 12, 3, 0, 0, 0)
        """
        cdef TernaryMatrix T
        cdef long i, j, d, r, x, y
        global GF3

        if self._t_invariant is not None:
            return
        T = (<TernaryMatrix>self._A).copy()   # Deprecated Sage matrix operation
        r = T.nrows()
        d = 0
        c = self._one
        i = 0
        while i < r - d:
            for j in range(i, r - d):
                if T.row_inner_product(j, j) != 0:   # Not a Sage matrix operation
                    if j > i:
                        T.swap_rows_c(i, j)
                    break
                if T.row_inner_product(i, j) != 0:   # Not a Sage matrix operation
                    if j > i:
                        T.add_multiple_of_row_c(i, j, 1, 0)
                    break
            x = T.row_inner_product(i, i)   # Not a Sage matrix operation
            if x == 0:
                d += 1
                T.swap_rows_c(i, r - d)
            else:
                c = c * GF3(x)
                for j in range(i + 1, r - d):
                    y = T.row_inner_product(i, j)   # Not a Sage matrix operation
                    if y == 0:
                        continue
                    if x == y:
                        T.row_subs(j, i)   # Not a Sage matrix operation
                    else:
                        T.add_multiple_of_row_c(j, i, 1, 0)
                i += 1

        TT = T.transpose()
        self._t_projection = TT._matrix_times_matrix_((T._matrix_times_matrix_(TT))._matrix_times_matrix_(T))
        F = frozenset()
        for i in range(r - d, r):
            F = F | frozenset(T.nonzero_positions_in_row(i))
        Fa = frozenset([j for j in range(len(self)) if self._t_projection.get(j, j) == 0]) - F   # Not a Sage matrix operation
        Fb = frozenset([j for j in range(len(self)) if self._t_projection.get(j, j) == 1]) - F   # Not a Sage matrix operation
        Fc = frozenset([j for j in range(len(self)) if self._t_projection.get(j, j) == -1]) - F   # Not a Sage matrix operation

        P = [Fa, Fb, Fc]
        p = []
        for a in range(3):
            for b in range(a + 1):
                x = 0
                for i in P[a]:
                    for j in P[b]:
                        if self._t_projection.get(i, j) != 0:   # Not a Sage matrix operation
                            x += 1
                p.append(x)

        self._t_partition = tuple([F, Fa, Fb, Fc])
        self._t_invariant = tuple([d, c, len(F), len(Fa), len(Fb), len(Fc), p[0], p[1], p[2], p[3], p[4], p[5]])

    cpdef _invariant(self):
        r"""
        Return a matroid invariant.

        See [Pen2012]_ for more information.

        OUTPUT:

        A tuple ``(d, c, L, La, Lb, Lc, p0, p1, p2, p3, p4, p5)``, with the
        following interpretation:

        - ``d`` is the bicycle dimension
        - ``c`` is the character
        - ``(L, La, Lb, Lc)`` is the triple of lengths of the principal
          quadripartition
        - ``(p0, ..., p5)`` counts of edges in a characteristic graph of the
          matroid whose vertex set is the groundset of the matroid,
          restricted to the sets in the principal quadripartition

        EXAMPLES::

           sage: M = matroids.catalog.NonFano()
           sage: M._invariant()
           (0, 2, 0, 4, 3, 0, 12, 12, 3, 0, 0, 0)
        """
        if self._t_invariant is None:
            self._make_invariant()
        return self._t_invariant

    cpdef bicycle_dimension(self):
        r"""
        Return the bicycle dimension of the ternary matroid.

        The bicycle dimension of a linear subspace `V` is
        `\dim(V\cap V^\perp)`. The bicycle dimension of a matroid equals the
        bicycle dimension of its rowspace, and is a matroid invariant.
        See [Pen2012]_.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M.bicycle_dimension()
            0
        """
        if self._t_invariant is None:
            self._make_invariant()
        return self._t_invariant[0]

    cpdef character(self):
        r"""
        Return the character of the ternary matroid.

        For a linear subspace `V` over `GF(3)` with orthogonal basis
        `q_1, \ldots, q_k` the character equals the product of `|q_i|`
        modulo 3, where the product ranges over the `i` such that `|q_i|`
        is not divisible by 3. The character does not depend on the choice of
        the orthogonal basis. The character of a ternary matroid equals the
        character of its cocycle-space, and is an invariant for ternary
        matroids. See [Pen2012]_.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: M.character()
            2
        """
        if self._t_invariant is None:
            self._make_invariant()
        return self._t_invariant[1]

    cpdef _principal_quadripartition(self):
        r"""
        Return an ordered partition of the groundset.

        The partition groups each element `e` of the groundset
        according to the bicycle dimension and the character of `M/e`, where
        `M` is the present matroid.

        EXAMPLES::

            sage: M = matroids.catalog.N1()
            sage: print(M)
            N1: Ternary matroid of rank 5 on 10 elements, type 0+
            sage: P = M._principal_quadripartition()
            sage: for e in sorted(P[0]): print("{} {}".format(e, M/e))
            sage: for e in sorted(P[1]): print("{} {}".format(e, M/e))
            a Ternary matroid of rank 4 on 9 elements, type 1-
            b Ternary matroid of rank 4 on 9 elements, type 1-
            e Ternary matroid of rank 4 on 9 elements, type 1-
            f Ternary matroid of rank 4 on 9 elements, type 1-
            sage: for e in sorted(P[2]): print("{} {}".format(e, M/e))
            d Ternary matroid of rank 4 on 9 elements, type 0-
            i Ternary matroid of rank 4 on 9 elements, type 0-
            sage: for e in sorted(P[3]): print("{} {}".format(e, M/e))
            c Ternary matroid of rank 4 on 9 elements, type 0+
            g Ternary matroid of rank 4 on 9 elements, type 0+
            h Ternary matroid of rank 4 on 9 elements, type 0+
            j Ternary matroid of rank 4 on 9 elements, type 0+
        """
        if self._t_invariant is None:
            self._make_invariant()
        return tuple([[self._E[j] for j in self._t_partition[0]], [self._E[j] for j in self._t_partition[1]], [self._E[j] for j in self._t_partition[2]], [self._E[j] for j in self._t_partition[3]]])

    cpdef TernaryMatrix _projection(self):
        """
        Return the projection matrix onto the row space.

        This projection is determined modulo the bicycle space. See [Pen2012]_.

        INPUT:

        - Nothing

        OUTPUT:

        A ternary matrix `P`, so that the `e`-th column of `P` is the signed
        incidence vector of a cocycle `C` such that `C-e` is a cycle.
        The bicycles `e`-th column is thus determined up to the bicycle
        space. We output the restriction of `P` to rows and columns that are
        not in any bicycle.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = TernaryMatroid(matrix(matroids.catalog.R12()))
            sage: M._projection()
            12 x 12 TernaryMatrix
            [++00-0--0+++]
            [+-+000+0+-+0]
            [0+-0-00+-+0+]
            [000000000000]
            [-0-0-0+-+00+]
            [000000000000]
            [-+00+0000+--]
            [-0+0-00-+0-+]
            [0+-0+00++++-]
            [+-+000+0+-+0]
            [++0000--++00]
            [+0+0+0-+-00-]
        """
        if self._t_invariant is None:
            self._make_invariant()
        return self._t_projection

    cpdef _fast_isom_test(self, other):
        r"""
           Run a quick test to see if two ternary matroids are isomorphic.

           The test is based on comparing strong invariants, including bicycle
           dimension, character, and the principal quadripartition.
           See also [Pen2012]_ .

           INPUT:

           - ``other`` -- a ternary matroid

           OUTPUT:

           - ``True``, if ``self`` is isomorphic to ``other``;
           - ``False``, if ``self`` is not isomorphic to ``other``;
           - ``None``, if the test is inconclusive

           EXAMPLES::

               sage: M = matroids.catalog.T8()
               sage: N = matroids.catalog.P8()
               sage: M._fast_isom_test(N)
               False
           """
        if self._invariant() != other._invariant():
            return False
        return None

    # minors, dual

    cpdef _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``
        - ``deletions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.P8()
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))
            Ternary matroid of rank 3 on 5 elements, type 0-
        """
        self._move_current_basis(contractions, deletions)
        bas = list(self.basis() - contractions)
        R = [self._prow[self._idx[b]] for b in bas]
        F = self.groundset() - (deletions | contractions)
        C = [self._idx[f] for f in F]
        A, C2 = (<TernaryMatrix>self._A).matrix_from_rows_and_columns_reordered(R, C)
        return TernaryMatroid(matrix=A,
                             groundset=[self._E[c] for c in C2],
                             basis=bas,
                             keep_initial_representation=False)

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data obey the matroid axioms.

        Since this is a linear matroid over the field `\GF{3}`, this is always
        the case.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: ``True``, or ``(True, {})``

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(3), [[]]))
            sage: M.is_valid()
            True
        """
        return True if not certificate else (True, {})

    # representability

    cpdef ternary_matroid(self, randomized_tests=1, verify = True):
        r"""
        Return a ternary matroid representing ``self``.

        INPUT:

        - ``randomized_tests`` -- ignored
        - ``verify`` -- ignored

        OUTPUT: a binary matroid

        ALGORITHM:

        ``self`` is a ternary matroid, so just return ``self``.

        .. SEEALSO::

            :meth:`M.ternary_matroid()
            <sage.matroids.matroid.Matroid.ternary_matroid>`

        EXAMPLES::

            sage: N = matroids.catalog.NonFano()
            sage: N.ternary_matroid() is N
            True
        """
        return self

    cpdef is_ternary(self, randomized_tests=1):
        r"""
        Decide if ``self`` is a binary matroid.

        INPUT:

        - ``randomized_tests`` -- ignored

        OUTPUT: boolean

        ALGORITHM:

        ``self`` is a ternary matroid, so just return ``True``.

        .. SEEALSO::

            :meth:`M.is_ternary() <sage.matroids.matroid.Matroid.is_ternary>`

        EXAMPLES::

            sage: N = matroids.catalog.NonFano()
            sage: N.is_ternary()
            True
        """
        return True

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle_ternary_matroid, (version, data))``, where
        ``unpickle_ternary_matroid`` is the name of a function that, when
        called with ``(version, data)``, produces a matroid isomorphic to
        ``self``. ``version`` is an integer (currently 0) and ``data`` is a
        tuple ``(A, E, B, name)`` where ``A`` is the representation
        matrix, ``E`` is the groundset of the matroid, ``B`` is the currently
        displayed basis, and ``name`` is a custom name.

        .. WARNING::

            Users should never call this function directly.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = TernaryMatroid(Matrix(GF(3), [[1, 0, 0, 1],
            ....:              [0, 1, 0, 1], [0, 0, 1, 1]]))
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('U34')
            sage: loads(dumps(M))
            U34
            sage: M = TernaryMatroid(Matrix(GF(3), [[1, 0, 1], [1, 0, 1]]))
            sage: loads(dumps(M)).representation()
            [1 0 1]
            [1 0 1]

        TESTS:

        Check that :issue:`23437` is fixed::

            sage: from sage.matroids.advanced import *
            sage: X_bin = matroids.catalog.Fano().representation()
            sage: X = Matrix(GF(3), X_bin)                                              # needs sage.rings.finite_rings
            sage: M = TernaryMatroid(matrix=X).dual()                                   # needs sage.rings.finite_rings
            sage: B = list(M.bases())
            sage: N = loads(dumps(M))
            sage: N.closure(frozenset({3}))                                             # needs sage.rings.finite_rings
            frozenset({3})
            sage: N.is_isomorphic(M)
            True
        """
        import sage.matroids.unpickling
        version = 0
        gs = self._E
        cdef list basis = [0] * self.full_rank()
        if self._representation is not None:
            A = self._representation
            basis = None
        else:
            A = self._A
            # current basis ordered so matrix cols form identity matrix:
            basis = self._current_rows_cols()[0]
        data = (A, gs, basis, self.get_custom_name())
        return sage.matroids.unpickling.unpickle_ternary_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.NonFano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: N = M.relabel({'g': 'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'x']

        TESTS::

            sage: M = matroids.catalog.NonFano()
            sage: f = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset_list()]
        M = TernaryMatroid(groundset=E, matrix=self._matrix_())
        return M

# Quaternary Matroids

cdef class QuaternaryMatroid(LinearMatroid):
    r"""
    Quaternary matroids.

    A quaternary matroid is a linear matroid represented over the finite field
    with four elements. See :class:`LinearMatroid` for a definition.

    The simplest way to create a ``QuaternaryMatroid`` is by giving only a
    matrix `A`. Then, the groundset defaults to ``range(A.ncols())``. Any
    iterable object `E` can be given as a groundset. If `E` is a list, then
    ``E[i]`` will label the `i`-th column of `A`. Another possibility is to
    specify a 'reduced' matrix `B`, to create the matroid induced by
    `A = [I\ \ B]`.

    INPUT:

    - ``matrix`` -- (default: ``None``) a matrix whose column vectors
      represent the matroid.
    - ``reduced_matrix`` -- (default: ``None``) a matrix `B` such that
      `[I\ \ B]` represents the matroid, where `I` is an identity matrix with
      the same number of rows as `B`. Only one of ``matrix`` and
      ``reduced_matrix`` should be provided.
    - ``groundset`` -- (default: ``None``) an iterable containing the element
      labels. When provided, must have the correct number of elements:
      the number of columns of ``matrix`` or the number of rows plus the
      number of columns of ``reduced_matrix``.
    - ``ring`` -- (default: ``None``) must be a copy of `\GF{4}`
    - ``keep_initial_representation`` -- boolean (default: ``True``); decides
      whether or not an internal copy of the input matrix should be preserved.
      This can help to see the structure of the matroid (e.g. in the case of
      graphic matroids), and makes it easier to look at extensions. However,
      the input matrix may have redundant rows, and sometimes it is desirable
      to store only a row-reduced copy.
    - ``basis`` -- (default: ``None``) when provided, this is an ordered
      subset of ``groundset``, such that the submatrix of ``matrix`` indexed
      by ``basis`` is an identity matrix. In this case, no row reduction takes
      place in the initialization phase.

    OUTPUT: a ``QuaternaryMatroid`` instance based on the data above

    .. NOTE::

        The recommended way to generate a quaternary matroid is through the
        :func:`Matroid() <sage.matroids.constructor.Matroid>` function. This
        is usually the preferred way, since it automatically chooses between
        ``QuaternaryMatroid`` and other classes. For direct access to the
        ``QuaternaryMatroid`` constructor, run::

            sage: from sage.matroids.advanced import *

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: GF4 = GF(4, 'x')
        sage: x = GF4.gens()[0]
        sage: A = Matrix(GF4, 2, 4, [[1, 0, 1, 1], [0, 1, 1, x]])
        sage: M = Matroid(A)
        sage: M
        Quaternary matroid of rank 2 on 4 elements
        sage: sorted(M.groundset())
        [0, 1, 2, 3]
        sage: Matrix(M)
        [1 0 1 1]
        [0 1 1 x]
        sage: M = Matroid(matrix=A, groundset='abcd')
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd']
        sage: GF4p = GF(4, 'y')
        sage: y = GF4p.gens()[0]
        sage: B = Matrix(GF4p, 2, 2, [[1, 1], [1, y]])
        sage: N = Matroid(reduced_matrix=B, groundset='abcd')
        sage: M == N
        False
    """
    def __init__(self, matrix=None, groundset=None, reduced_matrix=None, ring=None, keep_initial_representation=True, basis=None):
        """
        See class definition for full documentation.

        .. NOTE::

            The extra argument ``basis``, when provided, is an ordered list of
            elements of the groundset, ordered such that they index a standard
            identity matrix within ``matrix``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: QuaternaryMatroid(matrix=Matrix(GF(4, 'x'),         # indirect doctest            # needs sage.rings.finite_rings
            ....:                                 [[1, 0, 1, 1, 1],
            ....:                                  [0, 1, 1, 1, 1]]))
            Quaternary matroid of rank 2 on 5 elements
        """
        cdef QuaternaryMatrix A
        cdef long r, c
        cdef list P

        # Setup representation; construct displayed basis
        if matrix is not None:
            A = QuaternaryMatrix(matrix.nrows(), matrix.ncols(), M=matrix, ring=ring)
            if keep_initial_representation:
                self._representation = A.copy()   # Deprecated Sage matrix operation
            if basis is None:
                P = gauss_jordan_reduce(A, range(A.ncols()))
                A.resize(len(P))   # Not a Sage matrix operation
            self._A = A
        else:
            A = QuaternaryMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix, ring=ring)
            P = list(range(A.nrows()))
            self._A = A.prepend_identity()   # Not a Sage matrix operation

        # Setup groundset, BasisExchangeMatroid data
        if groundset is None:
            groundset = list(range(self._A.ncols()))
        else:
            if len(groundset) != self._A.ncols():
                raise ValueError("size of groundset does not match size of matrix")
        if basis is None:
            bas = [groundset[i] for i in P]
        else:
            bas = basis
        BasisExchangeMatroid.__init__(self, groundset, bas)

        # Setup index of displayed basis
        self._prow = <long* > sig_malloc((self._A.ncols()) * sizeof(long))
        for c in range(self._A.ncols()):
            self._prow[c] = -1
        if matrix is not None:
            if basis is None:
                for r in range(len(P)):
                    self._prow[P[r]] = r
            else:
                for r in range(self._A.nrows()):
                    self._prow[self._idx[basis[r]]] = r
        else:
            for r in range(self._A.nrows()):
                self._prow[r] = r

        self._zero = (<QuaternaryMatrix>self._A)._zero
        self._one = (<QuaternaryMatrix>self._A)._one
        self._x_zero = (<QuaternaryMatrix>self._A)._x_zero
        self._x_one = (<QuaternaryMatrix>self._A)._x_one

    cpdef base_ring(self):
        r"""
        Return the base ring of the matrix representing the matroid, in this
        case `\GF{4}`.

        EXAMPLES::

            sage: M = Matroid(ring=GF(4, 'y'), reduced_matrix=[[1, 0, 1],               # needs sage.rings.finite_rings
            ....:                                              [0, 1, 1]])
            sage: M.base_ring()                                                         # needs sage.rings.finite_rings
            Finite Field in y of size 2^2
        """
        return (<QuaternaryMatrix>self._A).base_ring()

    cpdef characteristic(self):
        """
        Return the characteristic of the base ring of the matrix representing
        the matroid, in this case `2`.

        EXAMPLES::

            sage: M = Matroid(ring=GF(4, 'y'), reduced_matrix=[[1, 0, 1],               # needs sage.rings.finite_rings
            ....:                                              [0, 1, 1]])
            sage: M.characteristic()                                                    # needs sage.rings.finite_rings
            2
        """
        return 2

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check if ``self.basis() - x + y`` is again a basis. Internal method.
        """
        return (<QuaternaryMatrix>self._A).is_nonzero(self._prow[x], y)

    cdef int _exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        cdef long p = self._prow[x]
        self._A.pivot(p, y)   # Not a Sage matrix operation
        self._prow[y] = p
        BasisExchangeMatroid._exchange(self, x, y)

    cdef  __fundamental_cocircuit(self, bitset_t C, long x):
        r"""
        Fill bitset `C` with the incidence vector of the `B`-fundamental cocircuit using ``x``. Internal method using packed elements.
        """
        bitset_union(C, (<QuaternaryMatrix>self._A)._M0[self._prow[x]], (<QuaternaryMatrix>self._A)._M1[self._prow[x]])

    cdef _coclosure_internal(self, bitset_t R, bitset_t F):
        """
        Bitpacked version of ``coclosure``.

        This function overrides the internal function BasisExchangeMatroid._coclosure_internal() of the parent class.
        The implementation should be more efficient for QuaternaryMatroid, due to the fact that in this class,
        __fundamental_cocircuit is much faster than __fundamental_circuit.
        """
        bitset_complement(R, F)
        bitset_difference(self._inside, self._current_basis, R)
        bitset_difference(self._outside, R, self._current_basis)
        self._move(self._inside, self._outside)

        bitset_copy(R, F)
        bitset_difference(self._inside, self._current_basis, F)
        cdef long y = bitset_first(self._inside)
        while y >= 0:
            self.__fundamental_cocircuit(self._outside, y)
            bitset_discard(self._outside, y)
            if bitset_issubset(self._outside, F):
                bitset_add(R, y)
            y = bitset_next(self._inside, y + 1)

    cdef  _exchange_value_internal(self, long x, long y):
        r"""
        Return the (x, y) entry of the current representation.
        """
        return (<QuaternaryMatrix>self._A).get(self._prow[x], y)   # Not a Sage matrix operation

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = Matroid(ring=GF(4, 'x'), matrix=[[1, 0, 1], [0, 1, 1]])           # needs sage.rings.finite_rings
            sage: M.rename()                                                            # needs sage.rings.finite_rings
            sage: repr(M)  # indirect doctest                                           # needs sage.rings.finite_rings
            'Quaternary matroid of rank 2 on 3 elements'
        """
        return f'Quaternary matroid of rank {self.rank()} on {self.size()} elements'

    cpdef _current_rows_cols(self, B=None):
        """
        Return the current row and column labels of a reduced matrix.

        INPUT:

        - ``B`` -- (default: ``None``) if provided, first find a basis having
          maximal intersection with ``B``

        OUTPUT:

        - ``R`` -- list of row indices; corresponds to the currently used
          internal basis
        - ``C`` -- list of column indices; corresponds to the complement of
          the current internal basis

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: M = matroids.catalog.Q10()
            sage: A = M._reduced_representation('efghi')
            sage: R, C = M._current_rows_cols()
            sage: (sorted(R), sorted(C))
            (['e', 'f', 'g', 'h', 'i'], ['a', 'b', 'c', 'd', 'j'])
            sage: R, C = M._current_rows_cols(B='abcde')
            sage: (sorted(R), sorted(C))
            (['a', 'b', 'c', 'd', 'e'], ['f', 'g', 'h', 'i', 'j'])
        """
        if B is not None:
            self._move_current_basis(B, set())
        basis = self.basis()
        rows = [0] * self.full_rank()
        cols = [0] * self.full_corank()
        c = 0
        for e in self._E:
            if e in basis:
                rows[self._prow[self._idx[e]]] = e
            else:
                cols[c] = e
                c += 1
        return rows, cols

    cpdef LeanMatrix _basic_representation(self, B=None):
        """
        Return a basic matrix representation of the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `M` representing the matroid, where `M[B'] = I` for a basis
        `B'` that maximally intersects the given set `B`. If not provided, the
        current basis used internally is chosen for `B'`. For a stable
        representation, use ``self.representation()``.

        .. NOTE::

            The method self.groundset_list() gives the labelling of the
            columns by the elements of the matroid. The matrix returned
            is a LeanMatrix subclass, which is intended for internal use only.
            Use the ``representation()`` method to get a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: M._basic_representation()                                             # needs sage.rings.finite_rings
            5 x 10 QuaternaryMatrix
            [100001x00y]
            [01000y1x00]
            [001000y1x0]
            [0001000y1x]
            [00001x00y1]
            sage: matrix(M._basic_representation('efghi'))                              # needs sage.rings.finite_rings
            [    1     0 x + 1     1     0     1     0     0     0     1]
            [    x x + 1     0     0     0     0     0     1     0     1]
            [    0     0     x x + 1     0     0     1     0     0     1]
            [    1     x     0     1     0     0     0     0     1     1]
            [    1     1     1     1     1     0     0     0     0     0]
        """
        if B is not None:
            self._move_current_basis(B, set())
        return self._A.copy()   # Deprecated Sage matrix operation

    cpdef LeanMatrix _reduced_representation(self, B=None):
        r"""
        Return a reduced representation of the matroid, i.e. a matrix `R` such
        that `[I\ \ R]` represents the matroid.

        INPUT:

        - ``B`` -- (default: ``None``) a set of elements of the groundset

        OUTPUT:

        A matrix `R` forming a reduced representation of the matroid, with
        rows labeled by a basis `B'` that maximally intersects the given set
        `B`. If not provided, the current basis used internally labels the
        rows.

        .. NOTE::

            The matrix returned is a LeanMatrix subclass, which is intended
            for internal use only. Use the ``representation()`` method to get
            a Sage matrix.

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: M._reduced_representation()                                           # needs sage.rings.finite_rings
            5 x 5 QuaternaryMatrix
            [1x00y]
            [y1x00]
            [0y1x0]
            [00y1x]
            [x00y1]
            sage: matrix(M._reduced_representation('efghi'))                            # needs sage.rings.finite_rings
            [    1     0 x + 1     1     1]
            [    x x + 1     0     0     1]
            [    0     0     x x + 1     1]
            [    1     x     0     1     1]
            [    1     1     1     1     0]
        """
        if B is not None:
            self._move_current_basis(B, set())
        _, cols = self._current_rows_cols()
        return self._A.matrix_from_rows_and_columns(range(self.full_rank()), [self._idx[e] for e in cols])

    cpdef _make_invariant(self):
        """
        Create an invariant.

        Internal method; see ``_invariant`` for full documentation.

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: M._invariant()  # indirect doctest                                    # needs sage.rings.finite_rings
            (0, 0, 5, 5, 20, 10, 25)
        """
        cdef QuaternaryMatrix Q, QT
        cdef long i, j, d, r

        if self._q_invariant is not None:
            return
        Q = (<QuaternaryMatrix>self._A).copy()   # Deprecated Sage matrix operation
        r = Q.nrows()
        d = 0
        i = 0
        while i < r - d:
            for j in range(i, r - d):
                if Q.row_inner_product(j, j) != 0:   # Not a Sage matrix operation
                    if j > i:
                        Q.swap_rows_c(i, j)
                    break
                y = Q.row_inner_product(i, j)   # Not a Sage matrix operation
                if y != 0:
                    if j > i:
                        Q.add_multiple_of_row_c(i, j, Q._x_zero * y, 0)
                    break
            x = Q.row_inner_product(i, i)   # Not a Sage matrix operation
            if x == 0:
                d += 1
                Q.swap_rows_c(i, r - d)
            else:
                for j in range(i + 1, r - d):
                    y = Q.row_inner_product(j, i)   # Not a Sage matrix operation
                    if y == 0:
                        continue
                    Q.add_multiple_of_row_c(j, i, y, 0)
                i += 1

        QT = Q.transpose()
        QT.conjugate()   # Not a Sage matrix operation
        self._q_projection = QT._matrix_times_matrix_((Q._matrix_times_matrix_(QT))._matrix_times_matrix_(Q))
        F = frozenset()
        for i in range(r - d, r):
            F = F | frozenset(Q.nonzero_positions_in_row(i))
        Fa = frozenset([j for j in range(len(self)) if self._q_projection.get(j, j) == 0]) - F   # Not a Sage matrix operation
        Fb = frozenset([j for j in range(len(self)) if self._q_projection.get(j, j) == 1]) - F   # Not a Sage matrix operation

        P = [Fa, Fb]
        p = []
        for a in range(2):
            for b in range(a + 1):
                x = 0
                for i in P[a]:
                    for j in P[b]:
                        if self._q_projection.get(i, j) != 0:   # Not a Sage matrix operation
                            x += 1
                p.append(x)

        self._q_partition = tuple([F, Fa, Fb])
        self._q_invariant = tuple([d, len(F), len(Fa), len(Fb), p[0], p[1], p[2]])

    cpdef _invariant(self):
        r"""
        Return a matroid invariant.

        See [Pen2012]_ for more information.

        OUTPUT:

        A tuple ``(d, Lm, L0, Lp, p0, p1, p2)``, with the following
        interpretation:

        - ``d`` is the bicycle dimension
        - ``(Lm, L0, Lp)`` is the triple of lengths of the principal
          tripartition
        - ``(p0, p1, p2)`` counts of edges in a characteristic graph of the
          matroid, whose vertices are the union of ``F_-`` and ``F_0`` from
          the principal tripartition

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: M._invariant()                                                        # needs sage.rings.finite_rings
            (0, 0, 5, 5, 20, 10, 25)
        """
        if self._q_invariant is None:
            self._make_invariant()
        return self._q_invariant

    cpdef bicycle_dimension(self):
        r"""
        Return the bicycle dimension of the quaternary matroid.

        The bicycle dimension of a linear subspace `V` is
        `\dim(V\cap V^\perp)`. We use the inner product
        `< v, w >=v_1 w_1^* + ... + v_n w_n^*`, where `w_i^*` is obtained from
        `w_i` by applying the unique nontrivial field automorphism of
        `\GF{4}`.

        The bicycle dimension of a matroid equals the bicycle dimension of its
        rowspace, and is a matroid invariant. See [Pen2012]_.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: M.bicycle_dimension()                                                 # needs sage.rings.finite_rings
            0
        """
        if self._q_invariant is None:
            self._make_invariant()
        return self._q_invariant[0]

    cpdef _principal_tripartition(self):
        r"""
        Return the principal tripartition of the quaternary matroid.

        The principal tripartition is a partition `(F_{-1}, F_0, F_{1})` of
        the groundset. A defining property is the following. It is
        straightforward that if the bicycle dimension of a matroid `M` is `k`,
        then the bicycle dimension of `M\setminus e' is one of `k-1, k, k + 1`
        for each element `e` of `M`. Then if `F_i` denotes the set of elements
        such that the bicycle dimension of `M\setminus e` is `k + i`, we
        obtain the principal tripartition `(F_{-1}, F_0, F_{1})` of `M`.
        See [Pen2012]_, [GR2001]_.

        OUTPUT:

        ``(F_{-1}, F_0, F_{1})``, the principal tripartition of the matroid.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: M = matroids.catalog.Q10().delete('a')
            sage: for F in M._principal_tripartition(): print(sorted(F))
            ['b', 'c', 'd', 'e', 'h', 'i']
            ['f', 'g', 'j']
            []
            sage: M.bicycle_dimension()
            1
            sage: for i in [-1, 0, 1]:
            ....:     print(sorted(e for e in M.groundset()
            ....:                  if (M.delete(e)).bicycle_dimension() == 1 + i))
            ['b', 'c', 'd', 'e', 'h', 'i']
            ['f', 'g', 'j']
            []
        """
        if self._q_invariant is None:
            self._make_invariant()
        P = self._q_partition
        return frozenset([self._E[e] for e in P[0]]), frozenset([self._E[e] for e in P[1]]), frozenset([self._E[e] for e in P[2]])

    cpdef _fast_isom_test(self, other):
        r"""
        Run a quick test to see if two quaternary matroids are isomorphic.

        The test is based on comparing the invariants returned by
        self._invariant().

        INPUT:

        - ``other`` -- a quaternary matroid

        OUTPUT:

        - ``True``, if ``self`` is isomorphic to ``other``;
        - ``False``, if ``self`` is not isomorphic to ``other``;
        - ``None``, if this test is inconclusive

        EXAMPLES::

           sage: M = matroids.catalog.Q10().delete('a')                                 # needs sage.rings.finite_rings
           sage: N = matroids.catalog.Q10().delete('b')                                 # needs sage.rings.finite_rings
           sage: M._fast_isom_test(N) is None                                           # needs sage.rings.finite_rings
           True
        """
        if self._invariant() != other._invariant():
            return False

    # minors, dual

    cpdef _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``
        - ``deletions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.Q10()                                            # needs sage.rings.finite_rings
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))              # needs sage.rings.finite_rings
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))             # needs sage.rings.finite_rings
            Quaternary matroid of rank 4 on 7 elements
        """
        self._move_current_basis(contractions, deletions)
        bas = list(self.basis() - contractions)
        R = [self._prow[self._idx[b]] for b in bas]
        F = self.groundset() - (deletions | contractions)
        C = [self._idx[f] for f in F]
        A, C2 = (<QuaternaryMatrix>self._A).matrix_from_rows_and_columns_reordered(R, C)
        return QuaternaryMatroid(matrix=A,
                             groundset=[self._E[c] for c in C2],
                             basis=bas,
                             keep_initial_representation=False)

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data obey the matroid axioms.

        Since this is a linear matroid over the field `\GF{4}`, this is always
        the case.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: ``True``, or ``(True, {})``

        EXAMPLES::

            sage: M = Matroid(Matrix(GF(4, 'x'), [[]]))                                 # needs sage.rings.finite_rings
            sage: M.is_valid()                                                          # needs sage.rings.finite_rings
            True
        """
        return True if not certificate else (True, {})

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle_quaternary_matroid, (version, data))``, where
        ``unpickle_quaternary_matroid`` is the name of a function that,
        when called with ``(version, data)``, produces a matroid isomorphic to
        ``self``. ``version`` is an integer (currently 0) and ``data`` is a
        tuple ``(A, E, B, name)`` where ``A`` is the representation
        matrix, ``E`` is the groundset of the matroid, ``B`` is the currently
        displayed basis, and ``name`` is a custom name.

        .. WARNING::

            Users should never call this function directly.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: M = Matroid(Matrix(GF(4, 'x'), [[1, 0, 0, 1], [0, 1, 0, 1],
            ....:            [0, 0, 1, 1]]))
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('U34')
            sage: loads(dumps(M))
            U34

        TESTS:

        Check that :issue:`23437` is fixed::

            sage: from sage.matroids.advanced import QuaternaryMatroid
            sage: X_bin = matroids.catalog.Fano().representation()
            sage: X = Matrix(GF(4), X_bin)                                              # needs sage.rings.finite_rings
            sage: M = QuaternaryMatroid(matrix=X).dual()                                # needs sage.rings.finite_rings
            sage: B = list(M.bases())                                                   # needs sage.rings.finite_rings
            sage: N = loads(dumps(M))                                                   # needs sage.rings.finite_rings
            sage: N.closure(frozenset({3}))                                             # needs sage.rings.finite_rings
            frozenset({3})
            sage: N.is_isomorphic(M)                                                    # needs sage.rings.finite_rings
            True
        """
        import sage.matroids.unpickling
        version = 0
        gs = self._E
        cdef list basis = [0] * self.full_rank()
        if self._representation is not None:
            A = self._representation
            basis = None
        else:
            A = self._A
            # current basis ordered so matrix cols form identity matrix:
            basis = self._current_rows_cols()[0]
        data = (A, gs, basis, self.get_custom_name())
        return sage.matroids.unpickling.unpickle_quaternary_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.RelaxedNonFano("abcdefg")
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: N = M.relabel({'g':'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'x']

        TESTS::

            sage: M = matroids.catalog.RelaxedNonFano("abcdefg")
            sage: f = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset_list()]
        M = QuaternaryMatroid(groundset=E, matrix=self._matrix_())
        return M

# Regular Matroids

cdef class RegularMatroid(LinearMatroid):
    r"""
    Regular matroids.

    A regular matroid is a linear matroid represented over the integers by a
    totally unimodular matrix.

    The simplest way to create a ``RegularMatroid`` is by giving only a matrix
    `A`. Then, the groundset defaults to ``range(A.ncols())``.
    Any iterable object `E` can be given as a groundset. If `E` is a list, then ``E[i]`` will label the `i`-th column of `A`.
    Another possibility is to specify a 'reduced' matrix `B`, to create the matroid induced by `A = [ I B ]`.

    INPUT:

    - ``matrix`` -- (default: ``None``) a matrix whose column vectors
      represent the matroid.
    - ``reduced_matrix`` -- (default: ``None``) a matrix `B` such that
      `[I\ \ B]` represents the matroid, where `I` is an identity matrix with
      the same number of rows as `B`. Only one of ``matrix`` and
      ``reduced_matrix`` should be provided.
    - ``groundset`` -- (default: ``None``) an iterable containing the element
      labels. When provided, must have the correct number of elements: the
      number of columns of ``matrix`` or the number of rows plus the number of
      columns of ``reduced_matrix``.
    - ``ring`` -- (default: ``None``) ignored
    - ``keep_initial_representation`` -- boolean (default: ``True``); decides
      whether or not an internal copy of the input matrix should be preserved.
      This can help to see the structure of the matroid (e.g. in the case of
      graphic matroids), and makes it easier to look at extensions. However,
      the input matrix may have redundant rows, and sometimes it is desirable
      to store only a row-reduced copy.
    - ``basis`` -- (default: ``None``) when provided, this is an ordered
      subset of ``groundset``, such that the submatrix of ``matrix`` indexed
      by ``basis`` is an identity matrix. In this case, no row reduction takes
      place in the initialization phase.

    OUTPUT: a ``RegularMatroid`` instance based on the data above

    .. NOTE::

        The recommended way to generate a regular matroid is through the
        :func:`Matroid() <sage.matroids.constructor.Matroid>` function. This
        is usually the preferred way, since it automatically chooses between
        ``RegularMatroid`` and other classes. Moreover, it will test whether
        the input actually yields a regular matroid, unlike this class.
        For direct access to the ``RegularMatroid`` constructor, run::

            sage: from sage.matroids.advanced import *

    .. WARNING::

        No checks are performed to ensure the input data form an actual regular
        matroid! If not, the behavior is unpredictable, and the internal
        representation can get corrupted. If in doubt, run
        :meth:`self.is_valid() <RegularMatroid.is_valid>` to ensure the data
        are as desired.

    EXAMPLES::

        sage: A = Matrix(ZZ, 2, 4, [[1, 0, 1, 1], [0, 1, 1, 1]])
        sage: M = Matroid(A, regular=True); M                                           # needs sage.graphs
        Regular matroid of rank 2 on 4 elements with 5 bases
        sage: sorted(M.groundset())                                                     # needs sage.graphs
        [0, 1, 2, 3]
        sage: Matrix(M)                                                                 # needs sage.graphs
        [1 0 1 1]
        [0 1 1 1]
        sage: M = Matroid(matrix=A, groundset='abcd', regular=True)                     # needs sage.graphs
        sage: sorted(M.groundset())                                                     # needs sage.graphs
        ['a', 'b', 'c', 'd']
    """
    def __init__(self, matrix=None, groundset=None, reduced_matrix=None, ring=None, keep_initial_representation=True):
        """
        See class definition for full documentation.

        .. NOTE::

            The extra argument ``basis``, when provided, is an ordered list of
            elements of the groundset, ordered such that they index a standard
            identity matrix within ``matrix``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: RegularMatroid(matrix=Matrix(ZZ, [[1, 0, 1, 1, 1],    # indirect doctest          # needs sage.graphs
            ....:                                   [0, 1, 1, 1, 1]]))
            Regular matroid of rank 2 on 5 elements with 7 bases
        """
        LinearMatroid.__init__(self, matrix, groundset, reduced_matrix, ring=ZZ, keep_initial_representation=keep_initial_representation)

    cdef list _setup_internal_representation(self, matrix, reduced_matrix, ring, keep_initial_representation):
        """
        Setup the internal representation matrix ``self._A`` and the array of
        row- and column indices ``self._prow``.

        Return the displayed basis.
        """
        cdef PlusMinusOneMatrix A
        cdef long r, c
        cdef list P
        if matrix is not None:
            reduced = False
            if not isinstance(matrix, PlusMinusOneMatrix):
                A = PlusMinusOneMatrix(matrix.nrows(), matrix.ncols(), M=matrix)
            else:
                A = (<PlusMinusOneMatrix>matrix).copy()   # Deprecated Sage matrix operation
            if keep_initial_representation:
                self._representation = A.copy()   # Deprecated Sage matrix operation
            P = gauss_jordan_reduce(A, range(A.ncols()))
            self._A = A.matrix_from_rows_and_columns(range(len(P)), [c for c in range(matrix.ncols()) if c not in P])
        else:
            reduced = True
            if not isinstance(reduced_matrix, PlusMinusOneMatrix):
                self._A = PlusMinusOneMatrix(reduced_matrix.nrows(), reduced_matrix.ncols(), M=reduced_matrix)
            else:
                self._A = (<PlusMinusOneMatrix>reduced_matrix).copy()   # Deprecated Sage matrix operation
            P = list(range(self._A.nrows()))
        self._prow = <long* > sig_malloc((self._A.nrows() + self._A.ncols()) * sizeof(long))
        if matrix is not None:
            for r in range(len(P)):
                self._prow[P[r]] = r
            r = 0
            for c in range(A.ncols()):
                if c not in P:
                    self._prow[c] = r
                    r += 1
        else:
            for r in range(self._A.nrows()):
                self._prow[r] = r
            for r in range(self._A.ncols()):
                self._prow[self._A.nrows() + r] = r
        return P

    cpdef base_ring(self):
        r"""
        Return the base ring of the matrix representing the matroid, in this
        case `\ZZ`.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.base_ring()
            Integer Ring
        """
        return ZZ

    cpdef characteristic(self):
        """
        Return the characteristic of the base ring of the matrix representing
        the matroid, in this case `0`.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.characteristic()
            0
        """
        return 0

    cdef bint _is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check if ``self.basis() - x + y`` is again a basis. Internal method.
        """
        return (<PlusMinusOneMatrix>self._A).is_nonzero(self._prow[x], self._prow[y])

    cdef int _exchange(self, long x, long y) except -1:
        """
        Put element indexed by ``x`` into basis, taking out element ``y``. Assumptions are that this is a valid basis exchange.

        .. NOTE::

            Safe for noncommutative rings.
        """
        cdef long px, py, r
        cdef int a, piv, pivi
        px = self._prow[x]
        py = self._prow[y]
        piv = (<PlusMinusOneMatrix>self._A).get(px, py)   # Not a Sage matrix operation
        pivi = piv  # NOTE: 1 and -1 are their own inverses.
        (<PlusMinusOneMatrix>self._A).rescale_row_c(px, pivi, 0)
        (<PlusMinusOneMatrix>self._A).set(px, py, pivi + 1)       # pivoting without column scaling. Add extra so column does not need adjusting   # Not a Sage matrix operation
        for r in range(self._A.nrows()):                 # if A and A' are the matrices before and after pivoting, then
            a = (<PlusMinusOneMatrix>self._A).get(r, py)       # ker[I A] equals ker[I A'] except for the labelling of the columns   # Not a Sage matrix operation
            if a and r != px:
                (<PlusMinusOneMatrix>self._A).add_multiple_of_row_c(r, px, -a, 0)
        (<PlusMinusOneMatrix>self._A).set(px, py, pivi)   # Not a Sage matrix operation
        self._prow[y] = px
        self._prow[x] = py
        BasisExchangeMatroid._exchange(self, x, y)

    cdef  _exchange_value_internal(self, long x, long y):
        r"""
        Return the (x, y) entry of the current representation.

        .. NOTE::

            This uses get_unsafe(), which returns a Sage ``Integer`` instance,
            rather than the ``int`` returned by ``get``. The
            advantage is that cross ratio tests will return rational numbers
            rather than unwarranted zeroes.
        """
        return (<PlusMinusOneMatrix>self._A).get_unsafe(self._prow[x], self._prow[y])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.rename()
            sage: repr(M)  # indirect doctest
            'Regular matroid of rank 5 on 10 elements with 162 bases'
        """
        return f'Regular matroid of rank {self.rank()} on {self.size()} elements with {self.bases_count()} bases'

    cpdef bases_count(self):
        """
        Count the number of bases.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(5), regular=True)                    # needs sage.graphs
            sage: M.bases_count()                                                       # needs sage.graphs
            125

        ALGORITHM:

        Since the matroid is regular, we use Kirchhoff's Matrix-Tree Theorem.
        See also :wikipedia:`Kirchhoff%27s_theorem`.
        """
        if self._bases_count is None:
            R = self._basic_representation()._matrix_()
            self._bases_count = (R * R.transpose()).det()
        return self._bases_count

    cpdef _projection(self):
        """
        Return the projection matrix onto the row space.

        INPUT:

        - Nothing

        OUTPUT:

        A matrix `P`, defined as follows. If `A` is a representation matrix
        of the matroid, then `Q = A^T (A A^T)^{-1} A`. Finally, `P` is equal
        to `Q` multiplied by the number of bases of the matroid.

        The matrices `P` and `Q` are independent of the choice of `A`, except for column
        scaling. The vector `Qx` is the orthogonal projection of
        the vector `x` onto the row space of `A`. For regular matroids,
        there is an extended Matrix Tree theorem that derives the fraction of
        bases containing a subset by computing the determinant of the
        principal submatrix of `Q` corresponding to that subset. See [Lyo2003]_ .
        Due to the scaling, the entries of `P` are integers.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = RegularMatroid(reduced_matrix=Matrix([[-1, 0, 1],
            ....:                                    [-1, 1, 0], [0, 1, -1]]))
            sage: M._projection()
            [ 8 -4  4 -4  0  4]
            [-4  8 -4 -4  4  0]
            [ 4 -4  8  0  4 -4]
            [-4 -4  0  8 -4 -4]
            [ 0  4  4 -4  8 -4]
            [ 4  0 -4 -4 -4  8]
        """
        if self._r_projection is None:
            R = self._basic_representation()._matrix_()
            self._r_projection = R.transpose() * (R * R.transpose()).adjugate() * R
        return self._r_projection

    cpdef _invariant(self):
        """
        Compute a regular matroid invariant.

        OUTPUT:

        The hash value of a list of pairs `(w, A[w])` and `(w, B[w])` and a number `N`,
        derived form the projection matrix `P` as obtained from ``self._projection()``
        as follows: `A[w]` counts the number of `i` such that `|P[i, i]|=w`, `B[w]` counts
        the number of pairs `(i, j)` such that `|P[i, j]|=w`, and `N` counts the number
        of triples `(i,j,k)` so that `P[i,j]*P'j,k]*P[k,i]` is negative.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: N = matroids.catalog.R10().dual()
            sage: O = matroids.catalog.R12()
            sage: M._invariant() == N._invariant()
            True
            sage: M._invariant() == O._invariant()
            False
        """
        from sage.matroids.utilities import cmp_elements_key
        # TODO: this currently uses Sage matrices internally. Perhaps dependence on those can be eliminated for further speed gains.
        if self._r_invariant is not None:
            return self._r_invariant
        cdef Matrix P
        P = self._projection()
        A = {}
        B = {}
        N = 0
        for i in range(P.nrows()):
            w = P.get_unsafe(i, i)
            if w != 0:
                if w in A:
                    A[w] += 1
                else:
                    A[w] = 1
            for j in range(i):
                w = abs(P.get_unsafe(i, j))
                if w != 0:
                    if w in B:
                        B[w] += 1
                    else:
                        B[w] = 1
                    for k in range(j):
                        if P.get_unsafe(i, j)*P.get_unsafe(j, k)*P.get_unsafe(k, i)<0:
                            N=N+1
        self._r_invariant = hash(tuple([tuple([(w, A[w]) for w in sorted(A, key=cmp_elements_key)]), tuple([(w, B[w]) for w in sorted(B, key=cmp_elements_key)]), N]))
        return self._r_invariant

    cpdef _hypergraph(self):
        """
        Create a bipartite digraph and a vertex partition.

        INPUT:

        - Nothing.

        OUTPUT:

        - ``PV`` -- a partition of the vertices of ``G``
        - ``tups`` -- list of pairs ``(x, y)``, where ``x`` denotes the
          color class of a part and ``y`` the number of elements in that part
        - ``G`` -- a graph

        All are derived from the entries of the projection matrix `P`. The
        partition ``PV`` groups vertices of the form `i` by the value of
        `P[i, i]`. Whenever `P[i, j]` is nonzero, there are edges `i - (i, j)`
        and `j - (i, j)`. Finally, the vertices `(i, j)` are grouped by value
        of `P[i, j]`.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: PV, tups, G = M._hypergraph()                                         # needs sage.graphs
            sage: G                                                                     # needs sage.graphs
            Digraph on 55 vertices
        """
        # NEW VERSION, Uses Sage'S Graph Isomorphism
        from sage.graphs.digraph import DiGraph
        from sage.matroids.utilities import cmp_elements_key
        if self._r_hypergraph is not None:
            return (self._hypergraph_vertex_partition, self._hypergraph_tuples, self._r_hypergraph)
        cdef Matrix P = self._projection()
        A = {}
        B = {}
        V = []
        E = []
        for i in range(P.nrows()):
            e = self._E[i]
            w = P.get_unsafe(i, i)
            if w != 0:
                if w in A:
                    A[w].append(e)
                else:
                    A[w] = [e]
                V.append(e)
            for j in range(i):
                f = self._E[j]
                w = abs(P.get_unsafe(i, j))
                if w != 0:
                    if w in B:
                        B[w].append(frozenset([e, f]))
                    else:
                        B[w] = [frozenset([e, f])]
                    E.append(frozenset([e, f]))
        self._hypergraph_vertex_partition = [[str(x) for x in A[w]] for w in sorted(A, key=cmp_elements_key)] + [['**' + str(x) for x in B[w]] for w in sorted(B, key=cmp_elements_key)]
        self._hypergraph_tuples = [(w, len(A[w])) for w in sorted(A, key=cmp_elements_key)] + [(w, len(B[w])) for w in sorted(B, key=cmp_elements_key)]
        G = DiGraph()
        G.add_vertices([str(x) for x in V] + ['**' + str(x) for x in E])
        # Note that Sage's Graph object attempts a sort on calling G.vertices(), which means vertices have to be well-behaved.
        for X in E:
            Y = list(X)
            G.add_edge(str(Y[0]), '**' + str(X))
            G.add_edge(str(Y[1]), '**' + str(X))

        self._r_hypergraph = G
        return self._hypergraph_vertex_partition, self._hypergraph_tuples, self._r_hypergraph
        # Remnant Of The Old Code That Was Not Yet Translated
        # To Sage'S Graph Isomorphism. Potential Speedup?
        # C = []
        # if len(A) < 5:
        #     for i in range(P.nrows()):
        #         for j in range(i):
        #             if P.get_unsafe(i, j) == 0:
        #                 continue
        #             for k in range(j):
        #                 w = P.get_unsafe(i, j)*P.get_unsafe(j, k)*P.get_unsafe(k, i)
        #                 if w < 0:
        #                     C.append(frozenset([self._E[i], self._E[j], self._E[k]]))
        #     self._r_hypergraph.add_color_class(C)
        #     self._r_hypergraph = self._r_hypergraph.max_refined()
        # return self._r_hypergraph

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate = True``, a dictionary giving
        the isomorphism or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = Matroid(groundset=list(range(6)),                                # needs sage.graphs
            ....:              graph=graphs.CompleteGraph(4), regular=True)
            sage: M1._is_isomorphic(M2)                                                 # needs sage.graphs
            True
            sage: M1._is_isomorphic(M2, certificate=True)                               # needs sage.graphs
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 4})

            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.catalog.Fano()
            sage: M1._is_isomorphic(M2)
            False
            sage: M1._is_isomorphic(M2.delete('a'))
            True
            sage: M1._is_isomorphic(M2.delete('a'), certificate=True)
            (True, {0: 'g', 1: 'b', 2: 'c', 3: 'e', 4: 'd', 5: 'f'})

        Check that :issue:`17316` was fixed::

            sage: from sage.matroids.advanced import *
            sage: Mnew = RegularMatroid(groundset=range(12), matrix=Matrix(ZZ,
            ....: [[ 1, 0, 0, 0, 1, 0, 0,-1,-1, 0, 1, 0],
            ....:  [ 0, 1, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0],
            ....:  [ 0, 0, 1, 0, 0,-1, 1, 0, 1, 0,-1, 0],
            ....:  [ 0, 0, 0, 1, 0, 0,-1, 1, 0, 0, 0, 0],
            ....:  [ 0, 0, 0, 0, 0, 1,-1, 0, 0, 1, 1, 0],
            ....:  [ 0, 0, 0, 0, 1, 0, 0,-1,-1, 0, 0, 1]]))
            sage: Nnew = RegularMatroid(groundset=range(12), matrix=Matrix(ZZ,
            ....: [[1,0,0,0,0,0,1,1,0,0,1,0],
            ....:  [0,1,0,0,0,0,1,0,1,0,1,0],
            ....:  [0,0,1,0,0,0,1,0,0,1,1,0],
            ....:  [0,0,0,1,0,0,1,1,0,0,0,1],
            ....:  [0,0,0,0,1,0,1,0,1,0,0,1],
            ....:  [0,0,0,0,0,1,1,0,0,1,0,1]]))
            sage: Mnew.is_isomorphic(Nnew)
            False
            sage: len(Mnew.circuits()) == len(Nnew.circuits())
            False
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        if isinstance(other, RegularMatroid):
            return self.is_field_isomorphic(other)
        else:
            return LinearMatroid._is_isomorphic(self, other)

    cpdef _fast_isom_test(self, other):
        r"""
        Run a quick test to see if two regular matroids are isomorphic.

        The test is based on:

        * A comparison of the number of bases, which may be computed
          efficiently through the matrix-tree lemma (see self.bases_count()).
        * A comparison of the orthogonal projection matrices of both matroids
          (see self._invariant()).
        * A isomorphism test which makes use of a hypergraph derived from the
          orthogonal projection matrix (see self._hypertest()).

        INPUT:

        - ``other`` -- a regular matroid

        OUTPUT:

        - ``True``, if ``self`` is isomorphic to ``other``;
        - ``False``, if ``self`` is not isomorphic to ``other``;
        - ``None``, if the test is inconclusive.

        EXAMPLES::

           sage: M = matroids.catalog.R10().delete('a')
           sage: N = matroids.catalog.R10().delete('b')
           sage: M._fast_isom_test(N)                                                   # needs sage.graphs
           True
        """
        if self.bases_count() != other.bases_count():
            return False
        if self._invariant() != other._invariant():
            return False
        if self.size() > 8:  # TODO: Optimize the cutoff. _hypertest() is slow for small matroids, and can be fast for larger ones.
            m = self._hypertest(other)
            if m is None:
                return False
            if self._is_field_isomorphism(other, m):
                return True

    cdef _hypertest(self, other):
        """
        Test if the hypergraphs associated with ``self`` and ``other`` are
        isomorphic, and if so return an isomorphism.

        INPUT:

        - ``other`` -- a ``RegularMatroid`` instance

        OUTPUT: dictionary, if the hypergraphs are isomorphic; ``None`` otherwise

        TESTS:

        Check that :issue:`22263` was fixed::

            sage: m1 = Matroid(graph='H?ABC~}')                                         # needs sage.graphs
            sage: m2 = Matroid(graph='H?ACNr}')                                         # needs sage.graphs
            sage: m1.is_isomorphic(m2)                                                  # needs sage.graphs
            False
        """
        from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic
        HS = self._hypergraph()
        HO = other._hypergraph()
        VO = []
        for X in HO[0]:
            VO.extend(X)
        m = isomorphic(HS[2], HO[2], HS[0], VO, 1, 1)
        if m:
            idx={str(f):f for f in other.groundset()}
            return {e:idx[m[str(e)]] for e in self.groundset() if str(e) in m}

    cpdef has_line_minor(self, k, hyperlines=None, certificate=False):
        r"""
        Test if the matroid has a `U_{2, k}`-minor.

        The matroid `U_{2, k}` is a matroid on `k` elements in which every
        subset of at most 2 elements is independent, and every subset of more
        than two elements is dependent.

        The optional argument ``hyperlines`` restricts the search space: this
        method returns ``True`` if `si(M/F)` is isomorphic to `U_{2, l}` with
        `l \geq k` for some `F` in ``hyperlines``, and ``False`` otherwise.

        INPUT:

        - ``k`` -- the length of the line minor
        - ``hyperlines`` -- (default: ``None``) a set of flats of codimension
          2. Defaults to the set of all flats of codimension 2.
        - ``certificate`` -- (default: ``False``) if ``True`` returns
          ``True, F``, where ``F`` is a flat and ``self.minor(contractions=F)``
          has a `U_{2,k}` restriction or ``False, None``

        OUTPUT: boolean or tuple

        .. SEEALSO::

            :meth:`Matroid.has_minor() <sage.matroids.matroid.Matroid.has_minor>`

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.has_line_minor(4)
            False
            sage: M.has_line_minor(4, certificate=True)
            (False, None)
            sage: M.has_line_minor(3)
            True
            sage: M.has_line_minor(3, certificate=True)
            (True, frozenset({'a', 'b', 'c', 'g'}))
            sage: M.has_line_minor(k=3, hyperlines=[['a', 'b', 'c'],
            ....:                                   ['a', 'b', 'd' ]])
            True
            sage: M.has_line_minor(k=3, hyperlines=[['a', 'b', 'c'],
            ....:                                   ['a', 'b', 'd' ]], certificate=True)
            (True, frozenset({'a', 'b', 'c'}))
        """
        if k > 3:
            if certificate:
                return False, None
            return False
        return Matroid.has_line_minor(self, k, hyperlines, certificate)

    cpdef _linear_extension_chains(self, F, fundamentals=None):
        r"""
        Create a list of chains that determine single-element extensions of
        this linear matroid representation.

        .. WARNING::

            Intended for internal use; does no input checking.

        INPUT:

        - ``F`` -- an independent set of elements
        - ``fundamentals`` -- (default: ``None``) a set elements of the base
          ring

        OUTPUT:

        A list of chains, so each single-element regular extension of this
        linear matroid, with support contained in ``F``, is
        given by one of these chains.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: len(M._linear_extension_chains(F=frozenset([0, 1, 2])))
            7
            sage: M._linear_extension_chains(F=frozenset())
            [{}]
            sage: M._linear_extension_chains(F=frozenset([1]))
            [{}, {1: 1}]
            sage: len(M._linear_extension_chains(F=frozenset([0, 1])))
            4
        """
        if fundamentals is None:
            fundamentals = set([1])
        return LinearMatroid._linear_extension_chains(self, F, fundamentals)

    cpdef bint is_graphic(self) noexcept:
        """
        Test if the regular matroid is graphic.

        A matroid is *graphic* if there exists a graph whose edge set equals
        the groundset of the matroid, such that a subset of elements of the
        matroid is independent if and only if the corresponding subgraph is
        acyclic.

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.is_graphic()
            False
            sage: M = Matroid(graphs.CompleteGraph(5), regular=True)                    # needs sage.graphs
            sage: M.is_graphic()                                                        # needs sage.graphs sage.rings.finite_rings
            True
            sage: M.dual().is_graphic()                                                 # needs sage.graphs
            False

        ALGORITHM:

        In a recent paper, Geelen and Gerards [GG2012]_ reduced the problem to
        testing if a system of linear equations has a solution. While not the
        fastest method, and not necessarily constructive (in the presence of
        2-separations especially), it is easy to implement.
        """
        return BinaryMatroid(reduced_matrix=self._reduced_representation()).is_graphic()

    cpdef is_valid(self, certificate=False):
        r"""
        Test if the data obey the matroid axioms.

        Since this is a regular matroid, this function tests if the
        representation matrix is *totally unimodular*, i.e. if all square
        submatrices have determinant in `\{-1, 0, 1\}`.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, or (boolean, dictionary)

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = Matroid(Matrix(ZZ, [[1, 0, 0, 1, 1, 0, 1],
            ....:                         [0, 1, 0, 1, 0, 1, 1],
            ....:                         [0, 0, 1, 0, 1, 1, 1]]),
            ....:             regular=True, check=False)
            sage: M.is_valid(certificate=True)
            (False, {'error': 'the representation matrix is not totally unimodular'})
            sage: M = Matroid(graphs.PetersenGraph())
            sage: M.is_valid()
            True
        """
        M = LinearMatroid(ring=QQ, reduced_matrix=self.representation(self.basis(), True, False))
        CR = M.cross_ratios()
        if CR.issubset(set([1])):
            return True if not certificate else (True, {})
        return False if not certificate else (False, {"error": "the representation matrix is not totally unimodular"})

    # representation

    cpdef bint is_regular(self) noexcept:
        r"""
        Return if ``self`` is regular.

        This is trivially ``True`` for a :class:`RegularMatroid`.

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: M.is_regular()
            True
        """
        return True

    cpdef binary_matroid(self, randomized_tests=1, verify = True):
        r"""
        Return a binary matroid representing ``self``.

        INPUT:

        - ``randomized_tests`` -- ignored
        - ``verify`` -- ignored

        OUTPUT: a binary matroid

        ALGORITHM:

        ``self`` is a regular matroid, so just cast ``self`` to a BinaryMatroid.

        .. SEEALSO::

            :meth:`M.binary_matroid()
            <sage.matroids.matroid.Matroid.binary_matroid>`

        EXAMPLES::

            sage: N = matroids.catalog.R10()
            sage: N.binary_matroid()
            Binary matroid of rank 5 on 10 elements, type (1, None)
        """
        A, E = self.representation(B = self.basis(), reduced = False, labels = True)
        return BinaryMatroid(matrix = A, groundset = E)

    cpdef is_binary(self, randomized_tests=1):
        r"""
        Decide if ``self`` is a binary matroid.

        INPUT:

        - ``randomized_tests`` -- ignored

        OUTPUT: boolean

        ALGORITHM:

        ``self`` is a regular matroid, so just return ``True``.

        .. SEEALSO::

            :meth:`M.is_binary() <sage.matroids.matroid.Matroid.is_binary>`

        EXAMPLES::

            sage: N = matroids.catalog.R10()
            sage: N.is_binary()
            True
        """
        return True

    cpdef ternary_matroid(self, randomized_tests=1, verify = True):
        r"""
        Return a ternary matroid representing ``self``.

        INPUT:

        - ``randomized_tests`` -- ignored
        - ``verify`` -- ignored

        OUTPUT: a ternary matroid

        ALGORITHM:

        ``self`` is a regular matroid, so just cast ``self`` to a TernaryMatroid.

        .. SEEALSO::

            :meth:`M.ternary_matroid()
            <sage.matroids.matroid.Matroid.ternary_matroid>`

        EXAMPLES::

            sage: N = matroids.catalog.R10()
            sage: N.ternary_matroid()
            Ternary matroid of rank 5 on 10 elements, type 4+
        """
        A, E = self.representation(B = self.basis(), reduced = False, labels = True)
        return TernaryMatroid(matrix = A, groundset = E)

    cpdef is_ternary(self, randomized_tests=1):
        r"""
        Decide if ``self`` is a ternary matroid.

        INPUT:

        - ``randomized_tests`` -- ignored

        OUTPUT: boolean

        ALGORITHM:

        ``self`` is a regular matroid, so just return ``True``.

        .. SEEALSO::

            :meth:`M.is_ternary() <sage.matroids.matroid.Matroid.is_ternary>`

        EXAMPLES::

            sage: N = matroids.catalog.R10()
            sage: N.is_ternary()
            True
        """
        return True

    # Copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle_regular_matroid, (version, data))``, where
        ``unpickle_regular_matroid`` is the name of a function that, when
        called with ``(version, data)``, produces a matroid isomorphic to
        ``self``. ``version`` is an integer (currently 0) and ``data`` is a
        tuple ``(A, E, reduced, name)`` where ``A`` is the representation
        matrix, ``E`` is the groundset of the matroid, ``reduced`` is a
        boolean indicating whether ``A`` is a reduced matrix, and ``name`` is
        a custom name.

        .. WARNING::

            Users should never call this function directly.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.R12()
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.rename('R_{12}')
            sage: loads(dumps(M))
            R_{12}
            sage: M = RegularMatroid(Matrix(QQ, [[1, 0, 1], [1, 0, 1]]))
            sage: N = loads(dumps(M))
            sage: N.representation()
            [1 0 1]
            [1 0 1]
        """
        import sage.matroids.unpickling
        cdef LeanMatrix A
        version = 0
        if self._representation is not None:
            A = self._representation
            gs = self._E
            reduced = False
        else:
            A = self._reduced_representation()
            rows, cols = self._current_rows_cols()
            gs = rows + cols
            reduced = True
        data = (A, gs, reduced, self.get_custom_name())
        return sage.matroids.unpickling.unpickle_regular_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.catalog.R10()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
            sage: N = M.relabel({'g': 'x'})
            sage: sorted(N.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'h', 'i', 'j', 'x']

        TESTS::

            sage: M = matroids.catalog.R10()
            sage: f = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([M._relabel_map(f)[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset_list()]
        M = RegularMatroid(groundset=E, matrix=self._matrix_())
        return M
