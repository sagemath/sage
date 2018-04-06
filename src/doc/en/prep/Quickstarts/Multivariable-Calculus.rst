.. -*- coding: utf-8 -*-

.. linkall

.. _prep-quickstart-multivariate-calculus:

Sage Quickstart for Multivariable Calculus
==========================================

This `Sage <http://www.sagemath.org>`_ quickstart tutorial was developed
for the MAA PREP Workshop "Sage: Using Open\-Source Mathematics Software
with Undergraduates" (funding provided by NSF DUE 0817071).  It is
licensed under the Creative Commons Attribution\-ShareAlike 3.0 license
(`CC BY\-SA <http://creativecommons.org/licenses/by-sa/3.0/>`_).

Because much related material was covered in the calculus tutorial, this
quickstart is designed to:

- Give concrete examples of some things not covered there, without huge
  amounts of commentary, and

- Remind the reader of a few of the things in that tutorial.

Vector Calculus
----------------

In Sage, vectors are primarily linear algebra objects, but they are
slowly becoming simultaneously analytic continuous functions.

Dot Product, Cross Product
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sage: v = vector(RR, [1.2, 3.5, 4.6])
    sage: w = vector(RR, [1.7,-2.3,5.2])
    sage: v*w
    17.9100000000000

::

    sage: v.cross_product(w)
    (28.7800000000000, 1.58000000000000, -8.71000000000000)

Lines and Planes
~~~~~~~~~~~~~~~~

Intersect :math:`x-2y-7z=6` with
:math:`\frac{x-3}{2}=\frac{y+4}{-3}=\frac{z-1}{1}`.  One way to plot the
equation of a line in three dimensions is with ``parametric_plot3d``.

::

    sage: # designed with intersection at t = 2, i.e. (7, -10, 3)
    sage: var('t, x, y')
    (t, x, y)
    sage: line = parametric_plot3d([2*t+3, -3*t-4, t+1], (t, 0, 4),color='red')
    sage: plane = plot3d((1/5)*(-12+x-2*y), (x, 4, 10), (y, -13,-7), opacity=0.5)
    sage: intersect=point3d([7,-10,3],color='black',size=30)
    sage: line+plane+intersect
    Graphics3d Object

Vector\-Valued Functions
~~~~~~~~~~~~~~~~~~~~~~~~

We can make vector-valued functions and do the usual analysis with them.

::

    sage: var('t')
    t
    sage: r=vector((2*t-4, t^2, (1/4)*t^3))
    sage: r
    (2*t - 4, t^2, 1/4*t^3)

::

    sage: r(t=5)
    (6, 25, 125/4)

The following makes the derivative also a vector\-valued expression.

::

    sage: velocity = r.diff(t) # velocity(t) = list(r.diff(t)) also would work
    sage: velocity
    (2, 2*t, 3/4*t^2)

Currently, this expression does *not* function as a function, so we need to
substitute explicitly.

::

    sage: velocity(t=1)
    (2, 2, 3/4)

::

    sage: T=velocity/velocity.norm()

::

    sage: T(t=1).n()
    (0.683486126173409, 0.683486126173409, 0.256307297315028)

Here we compute the arclength between :math:`t=0` and :math:`t=1` by
integrating the normalized derivative. As pointed out in the
:doc:`calculus tutorial <../Calculus>`, the syntax for
``numerical_integral`` is slightly nonstandard -- we just put in the
endpoints, not the variable.

::

    sage: arc_length = numerical_integral(velocity.norm(), 0,1)
    sage: arc_length
    (2.3169847271197814, 2.572369791753217e-14)

We can also plot vector fields, even in three dimensions.

::

    sage: x,y,z=var('x y z')
    sage: plot_vector_field3d((x*cos(z),-y*cos(z),sin(z)), (x,0,pi), (y,0,pi), (z,0,pi),colors=['red','green','blue'])
    Graphics3d Object

If we know a little vector calculus, we can also do line integrals.
Here, based on an example by Ben Woodruff of BYU-Idaho, we
compute a number of quantities in a physical three-dimensional setting.

Try to read through the entire code.  We make an auxiliary function that
will calculate :math:`\int (\text{integrand})\, dt`. We'll use our
formulas relating various differentials to :math:`dt` to easily specify
the integrands.

::

    sage: density(x,y,z)=x^2+z
    sage: r=vector([3*cos(t), 3*sin(t),4*t])
    sage: tstart=0
    sage: tend=2*pi
    sage: t = var('t')
    sage: t_range=(t,tstart,tend)
    sage: def line_integral(integrand):
    ....:     return RR(numerical_integral((integrand).subs(x=r[0], y=r[1],z=r[2]), tstart, tend)[0])
    sage: _ = var('x,y,z,t')
    sage: r_prime = diff(r,t)
    sage: ds=diff(r,t).norm()
    sage: s=line_integral(ds)
    sage: centroid_x=line_integral(x*ds)/s
    sage: centroid_y=line_integral(y*ds)/s
    sage: centroid_z=line_integral(z*ds)/s
    sage: dm=density(x,y,z)*ds
    sage: m=line_integral(dm)
    sage: avg_density = m/s
    sage: moment_about_yz_plane=line_integral(x*dm)
    sage: moment_about_xz_plane=line_integral(y*dm)
    sage: moment_about_xy_plane=line_integral(z*dm)
    sage: center_mass_x = moment_about_yz_plane/m
    sage: center_mass_y = moment_about_xz_plane/m
    sage: center_mass_z = moment_about_xy_plane/m
    sage: Ix=line_integral((y^2+z^2)*dm)
    sage: Iy=line_integral((x^2+z^2)*dm)
    sage: Iz=line_integral((x^2+y^2)*dm)
    sage: Rx = sqrt(Ix/m)
    sage: Ry = sqrt(Iy/m)
    sage: Rz = sqrt(Iz/m)

Finally, we can display everything in a nice :ref:`table <Tables>`.
Recall that we use the ``r"stuff"`` syntax to indicate "raw" strings
so that backslashes from LaTeX won't cause trouble.

.. skip

::

    sage: html.table([
    ....:     [r"Density $\delta(x,y)$", density],
    ....:     [r"Curve $\vec r(t)$",r],
    ....:     [r"$t$ range", t_range],
    ....:     [r"$\vec r'(t)$", r_prime],
    ....:     [r"$ds$, a little bit of arclength", ds],
    ....:     [r"$s$ - arclength", s],
    ....:     [r"Centroid (constant density) $\left(\frac{1}{m}\int x\,ds,\frac{1}{m}\int y\,ds, \frac{1}{m}\int z\,ds\right)$", (centroid_x,centroid_y,centroid_z)],
    ....:     [r"$dm=\delta ds$ - a little bit of mass", dm],
    ....:     [r"$m=\int \delta ds$ - mass", m],
    ....:     [r"average density $\frac{1}{m}\int ds$" , avg_density.n()],
    ....:     [r"$M_{yz}=\int x dm$ - moment about $yz$ plane", moment_about_yz_plane],
    ....:     [r"$M_{xz}=\int y dm$ - moment about $xz$ plane", moment_about_xz_plane],
    ....:     [r"$M_{xy}=\int z dm$ - moment about $xy$ plane", moment_about_xy_plane],
    ....:     [r"Center of mass $\left(\frac1m \int xdm, \frac1m \int ydm, \frac1m \int z dm\right)$", (center_mass_x, center_mass_y, center_mass_z)],
    ....:     [r"$I_x = \int (y^2+z^2) dm$", Ix],[r"$I_y=\int (x^2+z^2) dm$", Iy],[mp(r"$I_z=\int (x^2+y^2)dm$"), Iz],
    ....:     [r"$R_x=\sqrt{I_x/m}$", Rx],[mp(r"$R_y=\sqrt{I_y/m}"), Ry],[mp(r"$R_z=\sqrt{I_z/m}"),Rz]
    ....:     ])

Functions of Several Variables
-------------------------------

This connects directly to other issues of multivariable functions.

How to view these was mostly addressed in the various plotting
tutorials.  Here is a reminder of what can be done.

::

    sage: # import matplotlib.cm; matplotlib.cm.datad.keys()
    sage: # 'Spectral', 'summer', 'blues'
    sage: g(x,y)=e^-x*sin(y)
    sage: contour_plot(g, (x, -2, 2), (y, -4*pi, 4*pi), cmap = 'Blues', contours=10, colorbar=True)
    Graphics object consisting of 1 graphics primitive

Partial Differentiation
~~~~~~~~~~~~~~~~~~~~~~~

The following exercise is from Hass, Weir, and Thomas, University
Calculus, Exercise 12.7.35.  This function has a local minimum at :math:`(4,-2)`.

::

    sage: f(x, y) = x^2 + x*y + y^2 - 6*x + 2

Quiz: Why did we  *not*  need to declare the variables in this case?

::

    sage: fx(x,y)= f.diff(x)
    sage: fy(x,y) = f.diff(y)
    sage: fx; fy
    (x, y) |--> 2*x + y - 6
    (x, y) |--> x + 2*y

::

    sage: f.gradient()
    (x, y) |--> (2*x + y - 6, x + 2*y)

::

    sage: solve([fx==0, fy==0], (x, y))
    [[x == 4, y == -2]]

::

    sage: H = f.hessian()
    sage: H(x,y)
    [2 1]
    [1 2]

And of course if the Hessian has positive determinant and :math:`f_{xx}`
is positive, we have a local minimum.

.. skip

::

    sage: html("$f_{xx}=%s$"%H(4,-2)[0,0])
    sage: html("$D=%s$"%H(4,-2).det())

Notice how we were able to use many things we've done up to now to solve
this.

- Matrices

- Symbolic functions

- Solving

- Differential calculus

- Special formatting commands

- And, below, plotting!

::

    sage: plot3d(f,(x,-5,5),(y,-5,5))+point((4,-2,f(4,-2)),color='red',size=20)
    Graphics3d Object

Multiple Integrals and More
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Naturally, there is lots more that one can do.

::

    sage: f(x,y)=x^2*y
    sage: # integrate in the order dy dx
    sage: f(x,y).integrate(y,0,4*x).integrate(x,0,3)
    1944/5
    sage: # another way to integrate, and in the opposite order too
    sage: integrate( integrate(f(x,y), (x, y/4, 3)), (y, 0, 12) )
    1944/5

::

    sage: var('u v')
    (u, v)
    sage: surface = plot3d(f(x,y), (x, 0, 3.2), (y, 0, 12.3), color = 'blue', opacity=0.3)
    sage: domain = parametric_plot3d([3*u, 4*(3*u)*v,0], (u, 0, 1), (v, 0,1), color = 'green', opacity = 0.75)
    sage: image = parametric_plot3d([3*u, 4*(3*u)*v, f(3*u, 12*u*v)], (u, 0, 1), (v, 0,1), color = 'green', opacity = 1.00)
    sage: surface+domain+image
    Graphics3d Object

Quiz: why did we need to declare variables this time?
