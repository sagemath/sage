More Interesting Examples with f2py
===================================

Let us now look at some more interesting examples using f2py. We
will implement a simple iterative method for solving laplace's
equation in a square. Actually, this implementation is taken from
http://www.scipy.org/PerformancePython?highlight=\%28performance\%29
page on the scipy website. It has lots of information on
implementing numerical algorithms in python.

The following fortran code implements a single iteration of a
relaxation method for solving Laplace's equation in a square.

.. CODE-BLOCK:: fortran

    %fortran
          subroutine timestep(u,n,m,dx,dy,error)
          double precision u(n,m)
          double precision dx,dy,dx2,dy2,dnr_inv,tmp,diff
          integer n,m,i,j
    cf2py intent(in) :: dx,dy
    cf2py intent(in,out) :: u
    cf2py intent(out) :: error
    cf2py intent(hide) :: n,m
          dx2 = dx*dx
          dy2 = dy*dy
          dnr_inv = 0.5d0 / (dx2+dy2)
          error = 0d0
          do 200,j=2,m-1
             do 100,i=2,n-1
                tmp = u(i,j)
                u(i,j) = ((u(i-1,j) + u(i+1,j))*dy2+
         &           (u(i,j-1) + u(i,j+1))*dx2)*dnr_inv
                diff = u(i,j) - tmp
                error = error + diff*diff
     100     continue
     200  continue
          error = sqrt(error)
          end

If you do

.. CODE-BLOCK:: ipycon

    timestep?

You find that you need pass timestep a numpy array u, and the grid
spacing dx,dy and it will return the updated u, together with an
error estimate. To apply this to actually solve a problem use this
driver code

.. CODE-BLOCK:: python

    import numpy
    j=complex(0,1)
    num_points=50
    u=numpy.zeros((num_points,num_points),dtype=float)
    pi_c=float(pi)
    x=numpy.r_[0.0:pi_c:num_points*j]
    u[0,:]=numpy.sin(x)
    u[num_points-1,:]=numpy.sin(x)
    def solve_laplace(u,dx,dy):
       iter =0
       err = 2
       while(iter <10000 and err>1e-6):
          (u,err)=timestep(u,dx,dy)
          iter+=1
       return (u,err,iter)

Now call the routine using

.. CODE-BLOCK:: python

    (sol,err,iter)=solve_laplace(u,pi_c/(num_points-1),pi_c/(num_points-1))

This solves the equation with boundary conditions that the right
and left edge of the square are half an oscillation of the sine
function. With a 51x51 grid that we are using I find that it takes
around .2 s to solve this requiring 2750 iterations. If you have
the gnuplot package installed (use optional
packages() to find its name and install
package to install it), then you can plot this using

.. CODE-BLOCK:: python

    import Gnuplot
    g=Gnuplot.Gnuplot(persist=1)
    g('set parametric')
    g('set data style lines')
    g('set hidden')
    g('set contour base')
    g('set zrange [-.2:1.2]')
    data=Gnuplot.GridData(sol,x,x,binary=0)
    g.splot(data)

To see what we have gained by using f2py let us compare the same
algorithm in pure python and a vectorized version using numpy
arrays.

.. CODE-BLOCK:: python

    def slowTimeStep(u,dx,dy):
        """Takes a time step using straight forward Python loops."""
        nx, ny = u.shape
        dx2, dy2 = dx**2, dy**2
        dnr_inv = 0.5/(dx2 + dy2)


        err = 0.0
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                tmp = u[i,j]
                u[i,j] = ((u[i-1, j] + u[i+1, j])*dy2 +
                          (u[i, j-1] + u[i, j+1])*dx2)*dnr_inv
                diff = u[i,j] - tmp
                err += diff*diff

        return u,numpy.sqrt(err)

    def numpyTimeStep(u,dx,dy):
        dx2, dy2 = dx**2, dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u_old=u.copy()
        # The actual iteration
        u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 +
                         (u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv
        v = (u - u_old).flat
        return u,numpy.sqrt(numpy.dot(v,v))

You can try these out by changing the timestep function used in our
driver routine. The python version is slow even on a 50x50 grid. It
takes 70 seconds to solve the system in 3000 iterations. It takes
the numpy routine 2 seconds to reach the error tolerance in around
5000 iterations. In contrast it takes the f2py routine around .2
seconds to reach the error tolerance using 3000 iterations. I
should point out that the numpy routine is not quite the same
algorithm since it is a jacobi iteration while the f2py one is
gauss-seidel. This is why the numpy version requires more
iterations. Even accounting for this you can see the f2py version
appears to be around 5 times faster than the numpy version.
Actually if you try this on a 500x500 grid I find that it takes the
numpy routine 30 seconds to do 500 iterations while it only takes
about 2 seconds for the f2py to do this. So the f2py version is
really about 15 times faster. On smaller grids each actual
iteration is relatively cheap and so the overhead of calling f2py
is more evident, on larger examples where the iteration is
expensive, the advantage of f2py is clear. Even on the small
example it is still very fast. Note that a 500x500 grid in python
would take around half an hour to do 500 iterations.

To my knowledge the fastest that you could implement this algorithm
in matlab would be to vectorize it exactly like the numpy routine
we have. Vector addition in matlab and numpy are comparable. So
unless there is some trick I don't know about, using f2py you can
interactively write code 15 times faster than anything you could
write in matlab (Please correct me if I'm wrong). You can actually
make the f2py version a little bit faster by using
intent(in,out,overwrite) and creating the initial numpy array using
order='FORTRAN'. This eliminates the unnecessary copying that is
occurring in memory.
