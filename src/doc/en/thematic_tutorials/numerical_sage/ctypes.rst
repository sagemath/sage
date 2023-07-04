Ctypes
======

Ctypes is a very interesting python package which lets you import
shared object libraries into python and call them directly. I
should say that even though this is called ctypes, it can be used
just as well to call functions from libraries written in fortran.
The only complication is you need to know what a fortran function
looks like to C. This is simple everything is a pointer, so if your
fortran function would be called as foo(A,N) where A is an array
and N is its length, then to call it from C it takes a pointer to
an array of doubles and a pointer to an int. The other thing to be
aware of is that from C, fortran functions usually have an
underscore appended. That is, a fortran function foo would appear
as foo
from C (this is usually the case but is compiler dependent). Having
said this, the following examples are in C.

As an example suppose you write the following simple C program

.. code-block:: c

    #include <stdio.h>

    int sum(double *x,int n)
    {
      int i;
      double counter;
      counter = 0;
      for(i=0;i<n;i++)
        {
          counter=counter+x[i];

        }
      return counter;
    }

which you want to call from python. First make a shared object
library by doing (at the command line)

.. CODE-BLOCK:: shell-session

    $ gcc -c sum.c
    $ gcc -shared -o sum.so sum.o

Note that on OSX -shared should be replaced by -dynamiclib and
sum.so should be called sum.dylib Then you can do

.. CODE-BLOCK:: python

    from ctypes import *
    my_sum=CDLL('sum.so')
    a=numpy.array(range(10),dtype=float)
    my_sum.sum(a.ctypes.data_as(c_void_p),int(10))

Note here that ``a.ctypes.data_as(c_void_p)`` returns a ctypes
object that is void pointer to the underlying
array of a. Note that even though sum takes a double\*, as long as
we have a pointer to the correct data it doesn't matter what its
type is since it will be automatically cast.

Note that actually there are other ways to pass in the required
array of doubles. For example

.. CODE-BLOCK:: python

    a=(c_double*10)()
    for i in range(10):
       a[i]=i
    my_sum.sum(a,int(10))

This example only uses ctypes. Ctypes has wrappers for C data types
so for example

.. CODE-BLOCK:: python

    a=c_double(10.4)

would create a ctypes double object which could be passed to a C
function. Note that there is a byref function that lets you pass
parameters by reference. This is used in the next example. c
double\*10, is a python object that represents an array of 10
doubles and

.. CODE-BLOCK:: python

    a=(c_double*10)()

sets a equal to an array of 10 doubles. I find this method is
usually less useful than using numpy arrays when the data is
mathematical as numpy arrays are more well integrated into python
and sage.

Here is an example of using ctypes to directly call lapack. Note
that this will only work if you have a lapack shared object library
on your system. Also on linux the file would be liblapack.so and
you will probably use dgesv
(OSX use CLAPACK hence the lack of the underscore).

.. CODE-BLOCK:: python

    from ctypes import *
    def ctypes_solve(m,b,n):
        a=CDLL('/usr/lib/liblapack.dylib')
        import numpy
        p=(c_int*n)()
        size=c_int(n)
        ones=c_int(1)
        ok=c_int(0)
        a.dgesv(byref(size),byref(ones),m.ctypes.data_as(c_void_p),
                byref(size),p,b.ctypes.data_as(c_void_p),byref(size),byref(ok))

For completeness, let us consider a way to solve the laplace
equation using C types. Suppose you have written a simple solver in
C and you want to call it from python so you can easily test
different boundary conditions. Your C program might look like
this.

.. code-block:: c

    #include <math.h>
    #include <stdio.h>

    double timestep(double *u,int nx,int ny,double dx,double dy)
    {
      double tmp, err, diff,dx2,dy2,dnr_inv;
      dx2=dx*dx;
      dy2=dy*dy;
      dnr_inv=0.5/(dx2+dy2);
      err = 0.0;
      int i,j;

    for (i=1; i<nx-1; ++i) {
      for (j=1; j<ny-1; ++j) {
        tmp = u[i*nx+j];
        u[i*nx+j] = ((u[(i-1)*nx+j] + u[(i+1)*nx+j])*dy2 +
              (u[i*nx+j-1] + u[i*nx+j+1])*dx2)*dnr_inv;
        diff = u[i*nx+j] - tmp;
        err += diff*diff;
      }
    }

     return sqrt(err);
    }

    double solve_in_C(double *u,int nx,int ny,double dx,double dy)
    {
      double err;
      int iter;
      iter = 0;
      err = 1;
        while(iter <10000 && err > 1e-6)
          {
        err=timestep(u,nx,ny,dx,dy);
        iter++;
          }

      return err;
    }

We can compile it by running at the command line

.. CODE-BLOCK:: shell-session

     $ gcc -c laplace.c
     $ gcc -shared -o laplace.so laplace.o

Now in sage (notebook or command line) execute

.. CODE-BLOCK:: python

    from ctypes import *
    laplace=CDLL('/home/jkantor/laplace.so')
    laplace.timestep.restype=c_double
    laplace.solve_in_C.restype=c_double
    import numpy
    u=numpy.zeros((51,51),dtype=float)
    pi_c=float(pi)
    x=numpy.arange(0,pi_c+pi_c/50,pi_c/50,dtype=float)
    u[0,:]=numpy.sin(x)
    u[50,:]=numpy.sin(x)

    def solve(u):
      iter =0
      err = 2
      n=c_int(int(51))
      pi_c=float(pi/50)
      dx=c_double(pi_c)
      while(iter <5000 and err>1e-6):
         err=laplace.timestep(u.ctypes.data_as(c_void_p),n,n,dx,dx)
         iter+=1
         if(iter %50==0):
            print((err,iter))
      return (u,err,iter)

Note the line laplace.timestep.restype=c
double. By default ctypes assumes the return values are ints. If
they are not you need to tell it by setting restype to the correct
return type. If you execute the above code, then solve(u) will
solve the system. It is comparable to the fortran solution taking 
around .2 seconds. Alternatively you could do

.. CODE-BLOCK:: python

    n=c_int(int(51))
    dx=c_double(float(pi/50))
    laplace.solve_in_C(n.ctypes.data_as(c_void_p),n,n,dx,dx)

which computes the solution entirely in C. This is very fast.
Admittedly we could have had our fortran routines do the
entire solution at the Fortran level and we would have the same
speed.

As I said earlier you can just as easily call a shared object
library that is written in Fortran using ctypes. The key point is
it must be a shared object library and all fortran arguments are
passed by reference, that is as pointers or using byref. Also even
though we used very simple data types, it is possible to deal with
more complicated C structures. For this and more about ctypes see
http://python.net/crew/theller/ctypes/
