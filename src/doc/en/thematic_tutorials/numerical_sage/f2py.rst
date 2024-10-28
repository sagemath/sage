f2py
====

F2py is a very nice package that automatically wraps fortran code
and makes it callable from Python. The Fibonacci examples are taken
from the f2py webpage http://cens.ioc.ee/projects/f2py2e/.

From the notebook the magic %fortran will automatically compile any
fortran code in a cell and all the subroutines will become callable
functions (though the names will be converted to lowercase.) As an
example paste the following into a cell. It is important that the
spacing is correct as by default the code is treated as fixed
format fortran and the compiler will complain if things are not in
the correct column. To avoid this, you can write fortran 90 code
instead by making your first line !f90. There will be an example of
this later.

.. code-block:: fortran

    %fortran
    C FILE: FIB1.F
          SUBROUTINE FIB(A,N)
    C
    C     CALCULATE FIRST N FIBONACCI NUMBERS
    C
          INTEGER N
          REAL*8 A(N)
          DO I=1,N
             IF (I.EQ.1) THEN
                A(I) = 0.0D0
             ELSEIF (I.EQ.2) THEN
                A(I) = 1.0D0
             ELSE
                A(I) = A(I-1) + A(I-2)
             ENDIF
          ENDDO
          END
    C END FILE FIB1.F

Now evaluate it. It will be automatically compiled and imported
into Sage (though the name of imported function will be lowercase).
Now we want to try to call it, we need to somehow pass it an array
:math:`A`, and the length of the array :math:`N`. The way it
works is that numpy arrays will be automatically converted to
fortran arrays, and Python scalars converted to fortran scalars. So
to call fib we do the following.

.. CODE-BLOCK:: python

    import numpy
    m=numpy.array([0]*10,dtype=float)
    print(m)
    fib(m,10)
    print(m)

Note that fortran is a function that can be called on any string.
So if you have a fortran program in a file my
prog.f. Then you could do the following

.. CODE-BLOCK:: python

    f=open('my_prog.f','r')
    s=f.read()
    fortran(s)

Now all the functions in my
prog.f are callable.

It is possible to call external libraries in your fortran code. You
simply need to tell f2py to link them in. For example suppose we
wish to write a program to solve a linear equation using lapack (a
linear algebra library). The function we want to use is called
dgesv and it has the following signature.

.. code-block:: fortran

    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )


    *  N       (input) INTEGER
    *          The number of linear equations, i.e., the order of the
    *          matrix A.  N >= 0.
    *
    *  NRHS    (input) INTEGER
    *          The number of right hand sides, i.e., the number of columns
    *          of the matrix B.  NRHS >= 0.
    *
    *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    *          On entry, the N-by-N coefficient matrix A.
    *          On exit, the factors L and U from the factorization
    *          A = P*L*U; the unit diagonal elements of L are not stored.
    *
    *  LDA     (input) INTEGER
    *          The leading dimension of the array A.  LDA >= max(1,N).
    *
    *  IPIV    (output) INTEGER array, dimension (N)
    *          The pivot indices that define the permutation matrix P;
    *          row i of the matrix was interchanged with row IPIV(i).
    *
    *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
    *          On entry, the N-by-NRHS matrix of right hand side matrix B.
    *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
    *
    *  LDB     (input) INTEGER
    *          The leading dimension of the array B.  LDB >= max(1,N).
    *
    *  INFO    (output) INTEGER
    *          = 0:  successful exit
    *          < 0:  if INFO = -i, the i-th argument had an illegal value
    *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    *                has been completed, but the factor U is exactly
    *                singular, so the solution could not be computed.

we could do the following. Note that the order that library are in
the list actually matters as it is the order in which they are
passed to gcc. Also fortran.libraries is simply a list of names of
libraries that are linked in. You can just directly set this list.
So that fortran.libraries=['lapack','blas']is equivalent to the
following.

.. CODE-BLOCK:: python

    fortran.add_library('lapack')
    fortran.add_library('blas')

Now

.. CODE-BLOCK:: fortran

    %fortran
    !f90
    Subroutine LinearEquations(A,b,n)
    Integer n
    Real*8 A(n,n), b(n)
    Integer i, j, pivot(n), ok
    call DGESV(n, 1, A, n, pivot, b, n, ok)
    end

There are a couple things to note about this. As we remarked
earlier, if the first line of the code is !f90, then it will be
treated as fortran 90 code and does not need to be in fixed format.
To use the above try

.. CODE-BLOCK:: python

    a=numpy.random.randn(10,10)
    b=numpy.array(range(10),dtype=float)
    x=b.copy()
    linearequations(a,x,10)
    numpy.dot(a,x)

This will solve the linear system ax=b and store the result in b.
If your library is not in Sage's local/lib or in your path you can
add it to the search path using

.. CODE-BLOCK:: python

    fortran.add_library_path('path').

You can also directly set fortran.library
paths by assignment. It should be a list of paths (strings) to be
passed to gcc. To give you an idea of some more things you can do
with f2py, note that using intent statements you can control the
way the resulting Python function behaves a bit bitter. For example
consider the following modification of our original fibonacci
code.

.. CODE-BLOCK:: fortran

    C FILE: FIB3.F
          SUBROUTINE FIB(A,N)
    C
    C     CALCULATE FIRST N FIBONACCI NUMBERS
    C
          INTEGER N
          REAL*8 A(N)
    Cf2py intent(in) n
    Cf2py intent(out) a
    Cf2py depend(n) a
          DO I=1,N
             IF (I.EQ.1) THEN
                A(I) = 0.0D0
             ELSEIF (I.EQ.2) THEN
                A(I) = 1.0D0
             ELSE
                A(I) = A(I-1) + A(I-2)
             ENDIF
          ENDDO
          END
    C END FILE FIB3.F

Note the comments with the intent statements. This tells f2py that
:math:`n` is an input parameter and :math:`a` is the output.
This is called as

.. CODE-BLOCK:: python

    a=fib(10)

In general you will pass everything declared intent(in) to the
fortran function and everything declared intent(out) will be
returned in a tuple. Note that declaring something intent(in) means
you only care about its value before the function is called not
afterwards. So in the above n tells us how many fiboncci numbers to
compute we need to specify this as an input, however we don't need
to get n back as it doesn't contain anything new. Similarly A is
intent(out) so we don't need A to have a specific value
beforehand, we just care about the contents afterwards. F2py
generates a Python function so you only pass those declared
intent(in) and supplies empty workspaces for the remaining
arguments and it only returns those that are intent(out). All
arguments are intent(in) by default.

Consider now the following

.. CODE-BLOCK:: fortran

    %fortran
            Subroutine Rescale(a,b,n)
            Implicit none
            Integer n,i,j
            Real*8 a(n,n), b
            do i = 1,n
               do j=1,n
                 a(i,j)=b*a(i,j)
               end do
            end do
            end

You might be expecting Rescale(a,n) to rescale a numpy matrix a.
Alas this doesn't work. Anything you pass in is unchanged
afterwards. Note that in the fibonacci example above, the one
dimensional array was changed by the fortran code, similarly the
one dimensional vector b was replaced by its solution in the
example where we called lapack while the matrix A was not changed
even then dgesv says it modifies the input matrix. Why does this
not happen with the two dimensional array. Understanding this
requires that you are aware of the difference between how fortran
and C store arrays. Fortran stores a matrices using column ordering
while C stores them using row ordering. That is the matrix

.. math::

   \left(
   \begin{array}{ccc}
   0 & 1 &2\\
   3 & 4 & 5\\
   \end{array}
   \right)


is stored as

    :math:`(0\, 1\, 2\, 3\, 4\, 5\,) \,\,\,\, \text{ in C}`


    :math:`(0\, 3\,1\, 4\, 2\, 5) \,\,\,\, \text{ in Fortran}`


One dimensional arrays are stored the same in C and Fortran.
Because of this f2py allows the fortran code to operate on one
dimensional vectors in place, so your fortran code will change one
dimensional numpy arrays passed to it. However, since two
dimensional arrays are different by default f2py copies the numpy
array (which is stored in C format) into a second array that is in
the fortran format (i.e. takes the transpose) and that is what is
passed to the fortran function. We will see a way to get around
this copying later. First let us point one way of writing the
rescale function.

.. CODE-BLOCK:: fortran

    %fortran

            Subroutine Rescale(a,b,n)
            Implicit none
            Integer n,i,j
            Real*8 a(n,n), b
    Cf2py intent(in,out) a
            do i = 1,n
               do j=1,n
                 a(i,j)=b*a(i,j)
               end do
            end do
            end

Note that to call this you would use

.. CODE-BLOCK:: python

    b=rescale(a,2.0).

Note here I am not passing in :math:`n` which is the dimension of
:math:`a`. Often f2py can figure this out. This is a good time to
mention that f2py automatically generates some documentation for
the Python version of the function so you can check what you need
to pass to it and what it will return. To use this try

.. CODE-BLOCK:: ipycon

    rescale?

The intent(in,out) directives tells f2py to take the contents of
:math:`a` at the end of the subroutine and return them in a numpy
array. This still may not be what you want. The original
:math:`a` that you pass in is unmodified. If you want to modify
the original :math:`a` that you passed in use intent(inout). This
essentially lets your fortran code work with the data inplace.

.. CODE-BLOCK:: fortran

    %fortran

            Subroutine Rescale(a,b,n)
            Implicit none
            Integer n,i,j
            Real*8 a(n,n), b
    Cf2py intent(inout) a
            do i = 1,n
               do j=1,n
                 a(i,j)=b*a(i,j)
               end do
            end do
            end

If you wish to have fortran code work with numpy arrays in place,
you should make sure that your numpy arrays are stored in fortran's
format. You can ensure this by using the order='FORTRAN' keyword
when creating the arrays, as follows.

.. CODE-BLOCK:: python

    a=numpy.array([[1,2],[3,4]],dtype=float,order='FORTRAN')
    rescale(a,2.0)

After this executes, a will have the rescaled version of itself.
There is one final version which combines the previous two.

.. CODE-BLOCK:: fortran

    %fortran

            Subroutine Rescale(a,b,n)
            Implicit none
            Integer n,i,j
            Real*8 a(n,n), b
    Cf2py intent(in,out,overwrite) a
            do i = 1,n
               do j=1,n
                 a(i,j)=b*a(i,j)
               end do
            end do
            end

The (in,out,overwrite) intent says that if :math:`a` is in FORTRAN
ordering we work in place, however if its not we copy it and return
the contents afterwards. This is sort of the best of both worlds.
Note that if you are repeatedly passing large numpy arrays to
fortran code, it is very important to avoiding copying the array by
using (inout) or (in,out,overwrite). Remember though that your
numpy array must use Fortran ordering to avoid the copying.

For more examples and more advanced usage of F2py you should refer
to the f2py webpage http://cens.ioc.ee/projects/f2py2e/. The
command line f2py tool which is referred to in the f2py
documentation can be called from the Sage shell using

.. CODE-BLOCK:: ipycon

    !f2py
