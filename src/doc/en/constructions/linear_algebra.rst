.. _chapter-linear_algebra:

**************
Linear algebra
**************

.. index:
   pair: vector space; basis
   pair: vector space; subspace

.. _section-vector_space:

Vector spaces
=============

The ``VectorSpace`` command creates a vector space class, from which
one can create a subspace. Note the basis computed by Sage is
"row reduced".

::

    sage: V = VectorSpace(GF(2),8)
    sage: S = V.subspace([V([1,1,0,0,0,0,0,0]),V([1,0,0,0,0,1,1,0])])
    sage: S.basis()
    [(1, 0, 0, 0, 0, 1, 1, 0), (0, 1, 0, 0, 0, 1, 1, 0)]
    sage: S.dimension()
    2


.. index:
   pair: matrix; powers

.. _section-matrixpower:

Matrix powers
=============

How do I compute matrix powers in Sage? The syntax is illustrated by
the example below.

::

    sage: R = IntegerModRing(51)
    sage: M = MatrixSpace(R,3,3)
    sage: A = M([1,2,3, 4,5,6, 7,8,9])
    sage: A^1000*A^1007
    <BLANKLINE>
    [ 3  3  3]
    [18  0 33]
    [33 48 12]
    sage: A^2007
    <BLANKLINE>
    [ 3  3  3]
    [18  0 33]
    [33 48 12]

.. index:
   pair: matrix; kernel
   single: kernel; nullspace

.. _section-kernel:

Kernels
=======

The kernel is computed by applying the kernel method
to the matrix object. The following examples illustrate the syntax.

::

    sage: M = MatrixSpace(IntegerRing(),4,2)(range(8))
    sage: M.kernel()
    Free module of degree 4 and rank 2 over Integer Ring
    Echelon basis matrix:
    [ 1  0 -3  2]
    [ 0  1 -2  1]

A kernel of dimension one over :math:`\QQ`:

::

    sage: A = MatrixSpace(RationalField(),3)(range(9))
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

A trivial kernel:

::

    sage: A = MatrixSpace(RationalField(),2)([1,2,3,4])
    sage: A.kernel()
    Vector space of degree 2 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: M = MatrixSpace(RationalField(),0,2)(0)
    sage: M
    []
    sage: M.kernel()
    Vector space of degree 0 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: M = MatrixSpace(RationalField(),2,0)(0)
    sage: M.kernel()
    Vector space of degree 2 and dimension 2 over Rational Field
    Basis matrix:
    [1 0]
    [0 1]

Kernel of a zero matrix:

::

    sage: A = MatrixSpace(RationalField(),2)(0)
    sage: A.kernel()
    Vector space of degree 2 and dimension 2 over Rational Field
    Basis matrix:
    [1 0]
    [0 1]

Kernel of a non-square matrix:

::

    sage: A = MatrixSpace(RationalField(),3,2)(range(6))
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

The 2-dimensional kernel of a matrix over a cyclotomic field:

::

    sage: K = CyclotomicField(12); a = K.gen()
    sage: M = MatrixSpace(K,4,2)([1,-1, 0,-2, 0,-a^2-1, 0,a^2-1])
    sage: M
    [             1            -1]
    [             0            -2]
    [             0 -zeta12^2 - 1]
    [             0  zeta12^2 - 1]
    sage: M.kernel()
    Vector space of degree 4 and dimension 2 over Cyclotomic Field of order 12
     and degree 4
    Basis matrix:
    [               0                1                0     -2*zeta12^2]
    [               0                0                1 -2*zeta12^2 + 1]

A nontrivial kernel over a complicated base field.

::

    sage: K = FractionField(PolynomialRing(RationalField(),2,'x'))
    sage: M = MatrixSpace(K, 2)([[K.gen(1),K.gen(0)], [K.gen(1), K.gen(0)]])
    sage: M
    [x1 x0]
    [x1 x0]
    sage: M.kernel()
    Vector space of degree 2 and dimension 1 over Fraction Field of Multivariate
    Polynomial Ring in x0, x1 over Rational Field
    Basis matrix:
     [ 1 -1]

.. index:: Smith normal form, Hermite normal form, Frobenius normal form, rational canonical form

Other methods for integer matrices are ``elementary_divisors``,
``smith_form`` (for the Smith normal form), ``echelon_form``
for the Hermite normal form, ``frobenius`` for the
Frobenius normal form (rational canonical form).


There are many methods for matrices over a field such as
:math:`\QQ` or a finite field: ``row_span``, ``nullity``,
``transpose``, ``swap_rows``, ``matrix_from_columns``,
``matrix_from_rows``, among many others.

See the file ``matrix.py`` for further details.

.. index:: eigenvalues, eigenvectors

.. _section-eigen:

Eigenvectors and eigenvalues
============================

How do you compute eigenvalues and eigenvectors using Sage?

Sage has a full range of functions for computing eigenvalues and both
left and right eigenvectors and eigenspaces.  If our matrix is :math:`A`,
then the ``eigenmatrix_right`` (resp. ``eightmatrix_left``) command also
gives matrices :math:`D` and :math:`P` such that :math:`AP=PD` (resp.
:math:`PA=DP`.)

::

    sage: A = matrix(QQ, [[1,1,0],[0,2,0],[0,0,3]])
    sage: A
    [1 1 0]
    [0 2 0]
    [0 0 3]
    sage: A.eigenvalues()
    [3, 2, 1]
    sage: A.eigenvectors_right()
    [(3, [(0, 0, 1)], 1), (2, [(1, 1, 0)], 1), (1, [(1, 0, 0)], 1)]

    sage: A.eigenspaces_right()
    [(3,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [0 0 1]),
     (2,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [1 1 0]),
     (1,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [1 0 0])]

    sage: D, P = A.eigenmatrix_right()
    sage: D
    [3 0 0]
    [0 2 0]
    [0 0 1]
    sage: P
    [0 1 1]
    [0 1 0]
    [1 0 0]
    sage: A*P == P*D
    True

For eigenvalues outside the fraction field of the base ring of the matrix,
you can choose to have all the eigenspaces output when the algebraic closure
of the field is implemented, such as the algebraic numbers, ``QQbar``.  Or you
may request just a single eigenspace for each irreducible factor of the
characteristic polynomial, since the others may be formed through Galois
conjugation.  The eigenvalues of the matrix below are $\pm\sqrt{3}$ and
we exhibit each possible output.

Also, currently Sage does not implement multiprecision numerical eigenvalues
and eigenvectors, so calling the eigen functions on a matrix from
``CC`` or ``RR`` will probably give inaccurate and nonsensical results
(a warning is also printed).  Eigenvalues and eigenvectors of matrices with
floating point entries (over ``CDF`` and ``RDF``) can be obtained with the
"eigenmatrix" commands.  ::

    sage: MS = MatrixSpace(QQ, 2, 2)
    sage: A = MS([1,-4,1, -1])
    sage: A.eigenspaces_left(format='all')
    [(-1.732050807568878?*I,
      Vector space of degree 2 and dimension 1 over Algebraic Field
      User basis matrix:
      [                        1 -1 - 1.732050807568878?*I]),
     (1.732050807568878?*I,
      Vector space of degree 2 and dimension 1 over Algebraic Field
      User basis matrix:
      [                        1 -1 + 1.732050807568878?*I])]
    sage: A.eigenspaces_left(format='galois')
    [(a0,
      Vector space of degree 2 and dimension 1 over Number Field in a0 with defining polynomial x^2 + 3
      User basis matrix:
      [     1 a0 - 1])]

Another approach is to use the interface with Maxima:

::

    sage: A = maxima("matrix ([1, -4], [1, -1])")
    sage: eig = A.eigenvectors()
    sage: eig.sage()
    [[[-I*sqrt(3), I*sqrt(3)], [1, 1]], [[[1, 1/4*I*sqrt(3) + 1/4]], [[1, -1/4*I*sqrt(3) + 1/4]]]]

This tells us that :math:`\vec{v}_1 = [1,(\sqrt{3}i + 1)/4]` is
an eigenvector of :math:`\lambda_1 = - \sqrt{3}i` (which occurs
with multiplicity one) and
:math:`\vec{v}_2 = [1,(-\sqrt{3}i + 1)/4]` is an eigenvector of
:math:`\lambda_2 =  \sqrt{3}i` (which also occurs with
multiplicity one).

Here are two more examples:

::

    sage: A = maxima("matrix ([11, 0, 0], [1, 11, 0], [1, 3, 2])")
    sage: A.eigenvectors()
    [[[2,11],[1,2]],[[[0,0,1]],[[0,1,1/3]]]]
    sage: A = maxima("matrix ([-1, 0, 0], [1, -1, 0], [1, 3, 2])")
    sage: A.eigenvectors()
    [[[-1,2],[2,1]],[[[0,1,-1]],[[0,0,1]]]]

Warning: Notice how the ordering of the output is reversed, though
the matrices are almost the same.

Finally, you can use Sage's GAP interface as well to compute
"rational" eigenvalues and eigenvectors:

::

    sage: A = libgap([[1,2,3],[4,5,6],[7,8,9]]); A
    [ [ 1, 2, 3 ], [ 4, 5, 6 ], [ 7, 8, 9 ] ]
    sage: libgap(QQ).Eigenvectors(A)
    [ [ 1, -2, 1 ] ]
    sage: libgap(QQ).Eigenvalues(A)
    [ 0 ]

.. _section-rref:

Row reduction
=============

The row reduced echelon form of a matrix is computed as
in the following example.

::

    sage: M = MatrixSpace(RationalField(),2,3)
    sage: A = M([1,2,3, 4,5,6])
    sage: A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field
    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]
    sage: A.echelon_form()
    [      1       0 -1933/3]
    [      0       1  1550/3]

.. index::
   pair: matrix; characteristic polynomial

.. _section-characteristic:

Characteristic polynomial
=========================

The characteristic polynomial is a Sage method
for square matrices.

First a matrix over :math:`\ZZ`:

::

    sage: A = MatrixSpace(IntegerRing(),2)( [[1,2], [3,4]] )
    sage: f = A.charpoly()
    sage: f
    x^2 - 5*x - 2
    sage: f.parent()
    Univariate Polynomial Ring in x over Integer Ring

We compute the characteristic polynomial of a matrix over the
polynomial ring :math:`\ZZ[a]`:

::

    sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
    sage: M = MatrixSpace(R,2)([[a,1], [a,a+1]])
    sage: M
    [    a     1]
    [    a a + 1]
    sage: f = M.charpoly()
    sage: f
    x^2 + (-2*a - 1)*x + a^2
    sage: f.parent()
    Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over
    Integer Ring

    sage: M.trace()
    2*a + 1
    sage: M.determinant()
    a^2

We compute the characteristic polynomial of a matrix over the
multi-variate polynomial ring :math:`\ZZ[u,v]`:

::

    sage: R.<u,v> = PolynomialRing(ZZ,2)
    sage: A = MatrixSpace(R,2)([u,v,u^2,v^2])
    sage: f = A.charpoly(); f
    x^2 + (-v^2 - u)*x - u^2*v + u*v^2

It's a little difficult to distinguish the variables. To fix this,
we might want to rename the indeterminate "Z", which we can easily
do as follows:

.. link

::

    sage: f = A.charpoly('Z'); f
    Z^2 + (-v^2 - u)*Z - u^2*v + u*v^2

.. index::
   pair: solve; linear equations

Solving systems of linear equations
===================================

Using maxima, you can easily solve linear equations:

::

    sage: var('a,b,c')
    (a, b, c)
    sage: eqn = [a+b*c==1, b-a*c==0, a+b==5]
    sage: s = solve(eqn, a,b,c); s
    [[a == -1/4*I*sqrt(79) + 11/4, b == 1/4*I*sqrt(79) + 9/4, c == 1/10*I*sqrt(79) + 1/10], [a == 1/4*I*sqrt(79) + 11/4, b == -1/4*I*sqrt(79) + 9/4, c == -1/10*I*sqrt(79) + 1/10]]

You can even nicely typeset the solution in LaTeX:

::

    sage.: print(latex(s))
    ...

To have the above appear onscreen via xdvi, type ``view(s)``.

You can also solve linear equations symbolically using the
``solve`` command::

    sage: var('x,y,z,a')
    (x, y, z, a)
    sage: eqns = [x + z == y, 2*a*x - y == 2*a^2, y - 2*z == 2]
    sage: solve(eqns, x, y, z)
    [[x == a + 1, y == 2*a, z == a - 1]]

Here is a numerical Numpy example::

    sage: from numpy import arange, eye, linalg
    sage: A = eye(10)       ##   the 10x10 identity matrix
    sage: b = arange(1,11)
    sage: x = linalg.solve(A,b)

Another way to solve a system numerically is to use Sage's octave
interface::

    sage: M33 = MatrixSpace(QQ,3,3)
    sage: A   = M33([1,2,3,4,5,6,7,8,0])
    sage: V3  = VectorSpace(QQ,3)
    sage: b   = V3([1,2,3])
    sage: octave.solve_linear_system(A,b)    # optional - octave
    [-0.333333, 0.666667, 0]
