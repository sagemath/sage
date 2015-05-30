Plotting
========

We will plot a surface two ways. First we will use easyviz.
Consider the following code::

    import numpy
    from scitools import easyviz
    x = numpy.arange(-8,8,.2)
    xx,yy = numpy.meshgrid(x,x)
    r = numpy.sqrt(xx**2+yy**2) + 0.01
    zz = numpy.sin(r)/r
    easyviz.surfc(x,x,zz)

The function surfc takes a list of x coordinates, and y coordinates
and a numpy array z. Its plots a surface that has height z[i,j] at
the point (x[i],y[i]). Note the use of meshgrid, and vectorized
numpy functions that let us evaluate
:math:`\frac{\sin(\sqrt{x^2+y^2})+1}{\sqrt{x^2+y^2}+1}` over the
grid very easily. We discussed meshgrid at the beginning when we
were talking about numpy. Note that you can drag the plot around
with your mouse and look at it from different angles.

We can make this plot look a bit nicer by adding some shading and
nicer coloring and some labels as follows.

::

    import numpy
    RealNumber=float
    Integer =int
    from scitools import easyviz
    x = numpy.arange(-8,8,.2)
    xx,yy = numpy.meshgrid(x,x)
    r = numpy.sqrt(xx**2+yy**2) + 0.01
    zz = numpy.sin(r)/r
    l = easyviz.Light(lightpos=(-10,-10,5), lightcolor=(1,1,1))
    easyviz.surfc(x,x,zz,shading='interp',colormap=easyviz.jet(),
              zmin=-0.5,zmax=1,clevels=10,
              title='r=sqrt(x**2+y**2)+eps\nsin(r)/r',
              light=l,
              legend='sin',
              )

Let us now try to plot some vector fields. Consider the following
code

::

    import numpy
    from scitools import easyviz
    RealNumber=float
    Integer=int
    j=numpy.complex(0,1)
    w=numpy.zeros((5,5,5))
    u=w+1.0
    xx,yy,zz=numpy.mgrid[-1.0:1.0:5*j,-1:1:5*j,-1:1:5*j]
    easyviz.quiver3(xx,yy,zz,w,w,u)

This should plot a vector field that points up everywhere. The
arguments to quiver3 are 6, :math:`n\times n\times n` arrays. The
first three arrays are the location of the vectors, that is there
will be a vector at :math:`(xx[i,j,k],yy[i,j,k],zz[i,j,k])` for
:math:`0\le i,j,k < n`. The second three arrays are the
directions, i.e., the vector at
:math:`(xx[i,j,k],yy[i,j,k],zz[i,j,k])` points in the direction
:math:`(w[i,j,k],w[i,j,k],u[i,j,k])`.

Now let us give some examples with MayaVi. First lets see how to
plot a function like we did with easyviz.

::

    import numpy
    from mayavi.tools import imv
    x=numpy.arange(-8,8,.2)
    def f(x,y):
        r=numpy.sqrt(x**2+y**2)+.01
        return numpy.sin(r)/r
    imv.surf(x,x,f)

This will open mayavi, and display the plot of the function. The
first two arguments to surf are arrays :math:`x` and :math:`y`,
s.t. the function will be evaluated at :math:`(x[i],y[j])`. The
last argument is the function to graph. It probably looks a bit
different than the easyviz example. Lets try to make it look
similar to the easyviz example. First note that on the left there
is a list of filters and modules. Double-click the warpscalars
button in the filters menu, and change the scale factor from
:math:`1` to say :math:`5`. This should redraw the graph
similar to how easyviz drew it. There are quite a few other options
you can play around with. For example, next click on the module
surfacemap, and you will see you can make the graph transparent by
changing the opacity. You can also change it to a wireframe or make
it plot contours.

TODO: More examples
