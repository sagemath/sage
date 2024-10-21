# sage_setup: distribution = sagemath-modules
r"""
Matrices over an arbitrary ring

AUTHORS:

- William Stein

- Martin Albrecht: conversion to Pyrex

- Jaap Spies: various functions

- Gary Zablackis: fixed a sign bug in generic determinant.

- William Stein and Robert Bradshaw - complete restructuring.

- Rob Beezer - refactor kernel functions.

Elements of matrix spaces are of class ``Matrix`` (or a
class derived from Matrix). They can be either sparse or dense, and
can be defined over any base ring.

EXAMPLES:

We create the `2\times 3` matrix

.. MATH::

        \left(\begin{matrix} 1&2&3\\4&5&6 \end{matrix}\right)


as an element of a matrix space over `\QQ`::

    sage: M = MatrixSpace(QQ,2,3)
    sage: A = M([1,2,3, 4,5,6]); A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field

Alternatively, we could create A more directly as follows (which
would completely avoid having to create the matrix space)::

    sage: A = matrix(QQ, 2, [1,2,3, 4,5,6]); A
    [1 2 3]
    [4 5 6]

We next change the top-right entry of `A`. Note that matrix
indexing is `0`-based in Sage, so the top right entry is
`(0,2)`, which should be thought of as "row number
`0`, column number `2`".

::

    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]

Also notice how matrices print. All columns have the same width and
entries in a given column are right justified. Next we compute the
reduced row echelon form of `A`.

::

    sage: A.rref()
    [      1       0 -1933/3]
    [      0       1  1550/3]


Indexing
========

Sage has quite flexible ways of extracting elements or submatrices
from a matrix::

    sage: m=[(1, -2, -1, -1,9), (1, 8, 6, 2,2), (1, 1, -1, 1,4), (-1, 2, -2, -1,4)] ; M = matrix(m)
    sage: M
    [ 1 -2 -1 -1  9]
    [ 1  8  6  2  2]
    [ 1  1 -1  1  4]
    [-1  2 -2 -1  4]

Get the 2 x 2 submatrix of M, starting at row index and column index 1::

    sage: M[1:3,1:3]
    [ 8  6]
    [ 1 -1]

Get the 2 x 3 submatrix of M starting at row index and column index 1::

    sage: M[1:3,[1..3]]
    [ 8  6  2]
    [ 1 -1  1]

Get the second column of M::

    sage: M[:,1]
    [-2]
    [ 8]
    [ 1]
    [ 2]

Get the first row of M::

    sage: M[0,:]
    [ 1 -2 -1 -1  9]

Get the last row of M (negative numbers count from the end)::

    sage: M[-1,:]
    [-1  2 -2 -1  4]

More examples::

    sage: M[range(2),:]
    [ 1 -2 -1 -1  9]
    [ 1  8  6  2  2]
    sage: M[range(2),4]
    [9]
    [2]
    sage: M[range(3),range(5)]
    [ 1 -2 -1 -1  9]
    [ 1  8  6  2  2]
    [ 1  1 -1  1  4]

    sage: M[3,range(5)]
    [-1  2 -2 -1  4]
    sage: M[3,:]
    [-1  2 -2 -1  4]
    sage: M[3,4]
    4

    sage: M[-1,:]
    [-1  2 -2 -1  4]

    sage: A = matrix(ZZ,3,4, [3, 2, -5, 0, 1, -1, 1, -4, 1, 0, 1, -3]); A
    [ 3  2 -5  0]
    [ 1 -1  1 -4]
    [ 1  0  1 -3]

A series of three numbers, separated by colons, like ``n:m:s``, means
numbers from ``n`` up to (but not including) ``m``, in steps of ``s``.
So ``0:5:2`` means the sequence ``[0,2,4]``::

    sage: A[:,0:4:2]
    [ 3 -5]
    [ 1  1]
    [ 1  1]

    sage: A[1:,0:4:2]
    [1 1]
    [1 1]

    sage: A[2::-1,:]
    [ 1  0  1 -3]
    [ 1 -1  1 -4]
    [ 3  2 -5  0]

    sage: A[1:,3::-1]
    [-4  1 -1  1]
    [-3  1  0  1]

    sage: A[1:,3::-2]
    [-4 -1]
    [-3  0]

    sage: A[2::-1,3:1:-1]
    [-3  1]
    [-4  1]
    [ 0 -5]

We can also change submatrices using these indexing features::

    sage: M=matrix([(1, -2, -1, -1,9), (1, 8, 6, 2,2), (1, 1, -1, 1,4), (-1, 2, -2, -1,4)]); M
    [ 1 -2 -1 -1  9]
    [ 1  8  6  2  2]
    [ 1  1 -1  1  4]
    [-1  2 -2 -1  4]

Set the 2 x 2 submatrix of M, starting at row index and column index 1::

    sage: M[1:3,1:3] = [[1,0],[0,1]]; M
    [ 1 -2 -1 -1  9]
    [ 1  1  0  2  2]
    [ 1  0  1  1  4]
    [-1  2 -2 -1  4]

Set the 2 x 3 submatrix of M starting at row index and column index 1::

    sage: M[1:3,[1..3]] = M[2:4,0:3]; M
    [ 1 -2 -1 -1  9]
    [ 1  1  0  1  2]
    [ 1 -1  2 -2  4]
    [-1  2 -2 -1  4]

Set part of the first column of M::

    sage: M[1:,0]=[[2],[3],[4]]; M
    [ 1 -2 -1 -1  9]
    [ 2  1  0  1  2]
    [ 3 -1  2 -2  4]
    [ 4  2 -2 -1  4]

Or do a similar thing with a vector::

    sage: M[1:,0]=vector([-2,-3,-4]); M
    [ 1 -2 -1 -1  9]
    [-2  1  0  1  2]
    [-3 -1  2 -2  4]
    [-4  2 -2 -1  4]

Or a constant::

    sage: M[1:,0]=30; M
    [ 1 -2 -1 -1  9]
    [30  1  0  1  2]
    [30 -1  2 -2  4]
    [30  2 -2 -1  4]


Set the first row of M::

    sage: M[0,:]=[[20,21,22,23,24]]; M
    [20 21 22 23 24]
    [30  1  0  1  2]
    [30 -1  2 -2  4]
    [30  2 -2 -1  4]
    sage: M[0,:]=vector([0,1,2,3,4]); M
    [ 0  1  2  3  4]
    [30  1  0  1  2]
    [30 -1  2 -2  4]
    [30  2 -2 -1  4]
    sage: M[0,:]=-3; M
    [-3 -3 -3 -3 -3]
    [30  1  0  1  2]
    [30 -1  2 -2  4]
    [30  2 -2 -1  4]


    sage: A = matrix(ZZ,3,4, [3, 2, -5, 0, 1, -1, 1, -4, 1, 0, 1, -3]); A
    [ 3  2 -5  0]
    [ 1 -1  1 -4]
    [ 1  0  1 -3]

We can use the step feature of slices to set every other column::

    sage: A[:,0:3:2] = 5; A
    [ 5  2  5  0]
    [ 5 -1  5 -4]
    [ 5  0  5 -3]

    sage: A[1:,0:4:2] = [[100,200],[300,400]]; A
    [  5   2   5   0]
    [100  -1 200  -4]
    [300   0 400  -3]

We can also count backwards to flip the matrix upside down::

    sage: A[::-1,:]=A; A
    [300   0 400  -3]
    [100  -1 200  -4]
    [  5   2   5   0]


    sage: A[1:,3::-1]=[[2,3,0,1],[9,8,7,6]]; A
    [300   0 400  -3]
    [  1   0   3   2]
    [  6   7   8   9]

    sage: A[1:,::-2] = A[1:,::2]; A
    [300   0 400  -3]
    [  1   3   3   1]
    [  6   8   8   6]

    sage: A[::-1,3:1:-1] = [[4,3],[1,2],[-1,-2]]; A
    [300   0  -2  -1]
    [  1   3   2   1]
    [  6   8   3   4]


We save and load a matrix::

    sage: A = matrix(Integers(8),3,range(9))
    sage: loads(dumps(A)) == A
    True

MUTABILITY: Matrices are either immutable or not. When initially
created, matrices are typically mutable, so one can change their
entries. Once a matrix `A` is made immutable using
``A.set_immutable()`` the entries of `A`
cannot be changed, and `A` can never be made mutable again.
However, properties of `A` such as its rank, characteristic
polynomial, etc., are all cached so computations involving
`A` may be more efficient. Once `A` is made
immutable it cannot be changed back. However, one can obtain a
mutable copy of `A` using ``copy(A)``.

EXAMPLES::

    sage: A = matrix(RR,2,[1,10,3.5,2])
    sage: A.set_immutable()
    sage: copy(A) is A
    False

The echelon form method always returns immutable matrices with
known rank.

EXAMPLES::

    sage: A = matrix(Integers(8),3,range(9))
    sage: A.determinant()
    0
    sage: A[0,0] = 5
    sage: A.determinant()
    1
    sage: A.set_immutable()
    sage: A[0,0] = 5
    Traceback (most recent call last):
    ...
    ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).

Implementation and Design
-------------------------

Class Diagram (an x means that class is currently supported)::

    x Matrix
    x   Matrix_sparse
    x     Matrix_generic_sparse
    x     Matrix_integer_sparse
    x     Matrix_rational_sparse
          Matrix_cyclo_sparse
    x     Matrix_modn_sparse
          Matrix_RR_sparse
          Matrix_CC_sparse
          Matrix_RDF_sparse
          Matrix_CDF_sparse

    x  Matrix_dense
    x     Matrix_generic_dense
    x     Matrix_integer_dense
    x     Matrix_rational_dense
          Matrix_cyclo_dense    -- idea: restrict scalars to QQ, compute charpoly there, then factor
    x     Matrix_modn_dense
          Matrix_RR_dense
          Matrix_CC_dense
    x     Matrix_real_double_dense
    x     Matrix_complex_double_dense
    x     Matrix_complex_ball_dense

The corresponding files in the sage/matrix library code directory
are named

::

              [matrix] [base ring] [dense or sparse].

::

    New matrices types can only be implemented in Cython.

    *********** LEVEL 1  **********
    NON-OPTIONAL
    For each base field it is *absolutely* essential to completely
    implement the following functionality for that base ring:

       * __cinit__     -- should use check_allocarray from cysignals.memory
                          (only needed if allocate memory)
       * __init__      -- this signature: 'def __init__(self, parent, entries, copy, coerce)'
       * __dealloc__   -- use sig_free (only needed if allocate memory)
       * set_unsafe(self, size_t i, size_t j, x) -- doesn't do bounds or any other checks; assumes x is in self._base_ring
       * get_unsafe(self, size_t i, size_t j) -- doesn't do checks
       * __richcmp__    -- always the same (I don't know why its needed -- bug in PYREX)

    Note that the __init__ function must construct the all zero matrix if ``entries == None``.

    *********** LEVEL 2  **********

    IMPORTANT (and *highly* recommended):

    After getting the special class with all level 1 functionality to
    work, implement all of the following (they should not change
    functionality, except speed (always faster!) in any way):

       * def _pickle(self):
              return data, version
       * def _unpickle(self, data, int version)
              reconstruct matrix from given data and version; may assume _parent, _nrows, and _ncols are set.
              Use version numbers >= 0 so if you change the pickle strategy then
              old objects still unpickle.
       * cdef _list -- list of underlying elements (need not be a copy)
       * cdef _dict -- sparse dictionary of underlying elements
       * cdef _add_ -- add two matrices with identical parents
       * _matrix_times_matrix_c_impl -- multiply two matrices with compatible dimensions and
                                        identical base rings (both sparse or both dense)
       * cpdef _richcmp_ -- compare two matrices with identical parents
       * cdef _lmul_c_impl -- multiply this matrix on the right by a scalar, i.e., self * scalar
       * cdef _rmul_c_impl -- multiply this matrix on the left by a scalar, i.e., scalar * self
       * __copy__
       * __neg__

    The list and dict returned by _list and _dict will *not* be changed
    by any internal algorithms and are not accessible to the user.

    *********** LEVEL 3  **********
    OPTIONAL:

       * cdef _sub_
       * __invert__
       * _multiply_classical
       * __deepcopy__

    Further special support:
       * Matrix windows -- to support Strassen multiplication for a given base ring.
       * Other functions, e.g., transpose, for which knowing the
         specific representation can be helpful.

    .. NOTE::

       - For caching, use self.fetch and self.cache.

       - Any method that can change the matrix should call
         ``check_mutability()`` first.  There are also many fast cdef'd bounds checking methods.

       - Kernels of matrices
         Implement only a left_kernel() or right_kernel() method, whichever requires
         the least overhead (usually meaning little or no transposing).  Let the
         methods in the matrix2 class handle left, right, generic kernel distinctions.
"""
