.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

.. _vector_calc_curvilinear:

How to perform vector calculus in curvilinear coordinates
=========================================================

This tutorial introduces some vector calculus capabilities of SageMath within
the 3-dimensional Euclidean space. The corresponding tools have been developed
via the `SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/VectorCalculus/vector_calc_curvilinear.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/VectorCalculus/vector_calc_curvilinear.ipynb>`__ (``binder``).

Using spherical coordinates
---------------------------

To use spherical coordinates `(r,\theta,\phi)`  within the Euclidean 3-space
`\mathbb{E}^3`, it suffices to  declare the latter with the keyword
``coordinates='spherical'``::

    sage: E.<r,th,ph> = EuclideanSpace(coordinates='spherical')
    sage: E
    Euclidean space E^3

Thanks to the notation ``<r,th,ph>`` in the above declaration, the coordinates
`(r,\theta,\phi)` are immediately available as three symbolic variables ``r``,
``th`` and ``ph`` (there is no need to declare them via :func:`var`, i.e. to
type ``r, th, ph = var('r th ph')``)::

    sage: r is E.spherical_coordinates()[1]
    True
    sage: (r, th, ph) == E.spherical_coordinates()[:]
    True
    sage: type(r)
    <class 'sage.symbolic.expression.Expression'>

Moreover, the coordinate LaTeX symbols are already set::

    sage: latex(th)
    {\theta}

The coordinate ranges are::

    sage: E.spherical_coordinates().coord_range()
    r: (0, +oo); th: (0, pi); ph: [0, 2*pi] (periodic)

`\mathbb{E}^3` is endowed with the *orthonormal* vector frame
`(e_r, e_\theta, e_\phi)` associated with spherical coordinates::

    sage: E.frames()
    [Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph)),
     Vector frame (E^3, (e_r,e_th,e_ph))]

In the above output, ``(∂/∂r,∂/∂th,∂/∂ph)`` =
`\left(\frac{\partial}{\partial r}, \frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)`
is the *coordinate* frame associated with `(r,\theta,\phi)`; it is
not an orthonormal frame and will not be used below. The default frame is
the orthonormal one::

    sage: E.default_frame()
    Vector frame (E^3, (e_r,e_th,e_ph))


Defining a vector field
-----------------------

We define a vector field on `\mathbb{E}^3` from its components in
the orthonormal vector frame `(e_r,e_\theta,e_\phi)`::

    sage: v = E.vector_field(r*sin(2*ph)*sin(th)^2 + r,
    ....:                    r*sin(2*ph)*sin(th)*cos(th),
    ....:                    2*r*cos(ph)^2*sin(th), name='v')
    sage: v.display()
    v = (r*sin(2*ph)*sin(th)^2 + r) e_r + r*cos(th)*sin(2*ph)*sin(th) e_th
     + 2*r*cos(ph)^2*sin(th) e_ph

We can access to the components of `v` via the square bracket operator::

    sage: v[1]
    r*sin(2*ph)*sin(th)^2 + r
    sage: v[:]
    [r*sin(2*ph)*sin(th)^2 + r, r*cos(th)*sin(2*ph)*sin(th), 2*r*cos(ph)^2*sin(th)]

A vector field can evaluated at any point of `\mathbb{E}^3`::

    sage: p = E((1, pi/2, pi), name='p')
    sage: p
    Point p on the Euclidean space E^3
    sage: p.coordinates()
    (1, 1/2*pi, pi)
    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = e_r + 2 e_ph

We may define a vector field with generic components::

    sage: u = E.vector_field(function('u_r')(r,th,ph),
    ....:                    function('u_theta')(r,th,ph),
    ....:                    function('u_phi')(r,th,ph),
    ....:                    name='u')
    sage: u.display()
    u = u_r(r, th, ph) e_r + u_theta(r, th, ph) e_th + u_phi(r, th, ph) e_ph
    sage: u[:]
    [u_r(r, th, ph), u_theta(r, th, ph), u_phi(r, th, ph)]

Its value at the point `p` is then::

    sage: up = u.at(p)
    sage: up.display()
    u = u_r(1, 1/2*pi, pi) e_r + u_theta(1, 1/2*pi, pi) e_th
     + u_phi(1, 1/2*pi, pi) e_ph

Differential operators in spherical coordinates
-----------------------------------------------

The standard operators `\mathrm{grad}`, `\mathrm{div}`, `\mathrm{curl}`, etc.
involved in vector calculus are accessible as methods on scalar fields and
vector fields (e.g. ``v.div()``). However, to allow for standard mathematical
notations (e.g. ``div(v)``), let us import the functions
:func:`~sage.manifolds.operators.grad`, :func:`~sage.manifolds.operators.div`,
:func:`~sage.manifolds.operators.curl` and
:func:`~sage.manifolds.operators.laplacian`::

    sage: from sage.manifolds.operators import *


Gradient
~~~~~~~~

We first introduce a scalar field, via its expression in terms of Cartesian
coordinates; in this example, we consider a unspecified function of
`(r,\theta,\phi)`::

    sage: F = E.scalar_field(function('f')(r,th,ph), name='F')
    sage: F.display()
    F: E^3 → ℝ
       (r, th, ph) ↦ f(r, th, ph)

The value of `F` at a point::

    sage: F(p)
    f(1, 1/2*pi, pi)

The gradient of `F`::

    sage: grad(F)
    Vector field grad(F) on the Euclidean space E^3
    sage: grad(F).display()
    grad(F) = d(f)/dr e_r + d(f)/dth/r e_th + d(f)/dph/(r*sin(th)) e_ph
    sage: norm(grad(F)).display()
    |grad(F)|: E^3 → ℝ
       (r, th, ph) ↦ sqrt((r^2*(d(f)/dr)^2 + (d(f)/dth)^2)*sin(th)^2
        + (d(f)/dph)^2)/(r*sin(th))


Divergence
~~~~~~~~~~

The divergence of a vector field::

    sage: s = div(u)
    sage: s.display()
    div(u): E^3 → ℝ
       (r, th, ph) ↦ ((r*d(u_r)/dr + 2*u_r(r, th, ph)
        + d(u_theta)/dth)*sin(th) + cos(th)*u_theta(r, th, ph)
        + d(u_phi)/dph)/(r*sin(th))
    sage: s.expr().expand()
    2*u_r(r, th, ph)/r + cos(th)*u_theta(r, th, ph)/(r*sin(th))
     + diff(u_theta(r, th, ph), th)/r + diff(u_phi(r, th, ph), ph)/(r*sin(th))
     + diff(u_r(r, th, ph), r)

For `v`, we have::

    sage: div(v).expr()
    3

Curl
~~~~

The curl of a vector field::

    sage: s = curl(u)
    sage: s
    Vector field curl(u) on the Euclidean space E^3

::

    sage: s.display()
    curl(u) = (cos(th)*u_phi(r, th, ph) + sin(th)*d(u_phi)/dth
     - d(u_theta)/dph)/(r*sin(th)) e_r - ((r*d(u_phi)/dr + u_phi(r, th, ph))*sin(th)
     - d(u_r)/dph)/(r*sin(th)) e_th + (r*d(u_theta)/dr + u_theta(r, th, ph)
     - d(u_r)/dth)/r e_ph

For `v`, we have::

    sage: curl(v).display()
    curl(v) = 2*cos(th) e_r - 2*sin(th) e_th

The curl of a gradient is always zero::

    sage: curl(grad(F)).display()
    curl(grad(F)) = 0

The divergence of a curl is always zero::

    sage: div(curl(u)).display()
    div(curl(u)): E^3 → ℝ
       (r, th, ph) ↦ 0


Laplacian
~~~~~~~~~

The Laplacian of a scalar field::

    sage: s = laplacian(F)
    sage: s.display()
    Delta(F): E^3 → ℝ
       (r, th, ph) ↦ ((r^2*d^2(f)/dr^2 + 2*r*d(f)/dr
        + d^2(f)/dth^2)*sin(th)^2 + cos(th)*sin(th)*d(f)/dth
        + d^2(f)/dph^2)/(r^2*sin(th)^2)
    sage: s.expr().expand()
    2*diff(f(r, th, ph), r)/r + cos(th)*diff(f(r, th, ph), th)/(r^2*sin(th))
     + diff(f(r, th, ph), th, th)/r^2 + diff(f(r, th, ph), ph, ph)/(r^2*sin(th)^2)
     + diff(f(r, th, ph), r, r)

The Laplacian of a vector field::

    sage: Du = laplacian(u)
    sage: Du.display()
    Delta(u) = ((r^2*d^2(u_r)/dr^2 + 2*r*d(u_r)/dr - 2*u_r(r, th, ph)
     + d^2(u_r)/dth^2 - 2*d(u_theta)/dth)*sin(th)^2 - ((2*u_theta(r, th, ph)
     - d(u_r)/dth)*cos(th) + 2*d(u_phi)/dph)*sin(th) + d^2(u_r)/dph^2)/(r^2*sin(th)^2) e_r
     + ((r^2*d^2(u_theta)/dr^2 + 2*r*d(u_theta)/dr + 2*d(u_r)/dth + d^2(u_theta)/dth^2)*sin(th)^2
     + cos(th)*sin(th)*d(u_theta)/dth - 2*cos(th)*d(u_phi)/dph - u_theta(r, th, ph)
     + d^2(u_theta)/dph^2)/(r^2*sin(th)^2) e_th
     + ((r^2*d^2(u_phi)/dr^2 + 2*r*d(u_phi)/dr
     + d^2(u_phi)/dth^2)*sin(th)^2 + (cos(th)*d(u_phi)/dth + 2*d(u_r)/dph)*sin(th)
     + 2*cos(th)*d(u_theta)/dph - u_phi(r, th, ph) + d^2(u_phi)/dph^2)/(r^2*sin(th)^2) e_ph

Since this expression is quite lengthy, we may ask for a display component by
component::

    sage: Du.display_comp()
    Delta(u)^1 = ((r^2*d^2(u_r)/dr^2 + 2*r*d(u_r)/dr - 2*u_r(r, th, ph) + d^2(u_r)/dth^2
     - 2*d(u_theta)/dth)*sin(th)^2 - ((2*u_theta(r, th, ph) - d(u_r)/dth)*cos(th)
     + 2*d(u_phi)/dph)*sin(th) + d^2(u_r)/dph^2)/(r^2*sin(th)^2)
    Delta(u)^2 = ((r^2*d^2(u_theta)/dr^2 + 2*r*d(u_theta)/dr + 2*d(u_r)/dth
     + d^2(u_theta)/dth^2)*sin(th)^2 + cos(th)*sin(th)*d(u_theta)/dth
     - 2*cos(th)*d(u_phi)/dph - u_theta(r, th, ph) + d^2(u_theta)/dph^2)/(r^2*sin(th)^2)
    Delta(u)^3 = ((r^2*d^2(u_phi)/dr^2 + 2*r*d(u_phi)/dr + d^2(u_phi)/dth^2)*sin(th)^2
     + (cos(th)*d(u_phi)/dth + 2*d(u_r)/dph)*sin(th) + 2*cos(th)*d(u_theta)/dph
     - u_phi(r, th, ph) + d^2(u_phi)/dph^2)/(r^2*sin(th)^2)

We may expand each component::

    sage: for i in E.irange():
    ....:     s = Du[i].expand()
    sage: Du.display_comp()
    Delta(u)^1 = 2*d(u_r)/dr/r - 2*u_r(r, th, ph)/r^2
     - 2*cos(th)*u_theta(r, th, ph)/(r^2*sin(th)) + cos(th)*d(u_r)/dth/(r^2*sin(th))
     + d^2(u_r)/dth^2/r^2 - 2*d(u_theta)/dth/r^2 - 2*d(u_phi)/dph/(r^2*sin(th))
     + d^2(u_r)/dph^2/(r^2*sin(th)^2) + d^2(u_r)/dr^2
    Delta(u)^2 = 2*d(u_theta)/dr/r + 2*d(u_r)/dth/r^2 + cos(th)*d(u_theta)/dth/(r^2*sin(th))
     + d^2(u_theta)/dth^2/r^2 - 2*cos(th)*d(u_phi)/dph/(r^2*sin(th)^2)
     - u_theta(r, th, ph)/(r^2*sin(th)^2) + d^2(u_theta)/dph^2/(r^2*sin(th)^2)
     + d^2(u_theta)/dr^2
    Delta(u)^3 = 2*d(u_phi)/dr/r + cos(th)*d(u_phi)/dth/(r^2*sin(th))
     + d^2(u_phi)/dth^2/r^2 + 2*d(u_r)/dph/(r^2*sin(th))
     + 2*cos(th)*d(u_theta)/dph/(r^2*sin(th)^2) - u_phi(r, th, ph)/(r^2*sin(th)^2)
     + d^2(u_phi)/dph^2/(r^2*sin(th)^2) + d^2(u_phi)/dr^2

As a test, we may check that these formulas coincide with those of Wikipedia's
article :wikipedia:`Del in cylindrical and spherical coordinates <Del_in_cylindrical_and_spherical_coordinates#Del_formula>`.

Using cylindrical coordinates
-----------------------------

The use of cylindrical coordinates `(\rho,\phi,z)` in the Euclidean space
`\mathbb{E}^3` is on the same footing as that of spherical coordinates. To
start with, one has simply to declare::

    sage: E.<rh,ph,z> = EuclideanSpace(coordinates='cylindrical')

The coordinate ranges are then::

    sage: E.cylindrical_coordinates().coord_range()
    rh: (0, +oo); ph: [0, 2*pi] (periodic); z: (-oo, +oo)

The default vector frame is the orthonormal frame `(e_\rho,e_\phi,e_z)`
associated with cylindrical coordinates::

    sage: E.default_frame()
    Vector frame (E^3, (e_rh,e_ph,e_z))

and one may define vector fields from their components in that frame::

    sage: v = E.vector_field(rh*(1+sin(2*ph)), 2*rh*cos(ph)^2, z,
    ....:                    name='v')
    sage: v.display()
    v = rh*(sin(2*ph) + 1) e_rh + 2*rh*cos(ph)^2 e_ph + z e_z
    sage: v[:]
    [rh*(sin(2*ph) + 1), 2*rh*cos(ph)^2, z]

::

    sage: u = E.vector_field(function('u_rho')(rh,ph,z),
    ....:                    function('u_phi')(rh,ph,z),
    ....:                    function('u_z')(rh,ph,z),
    ....:                    name='u')
    sage: u.display()
    u = u_rho(rh, ph, z) e_rh + u_phi(rh, ph, z) e_ph + u_z(rh, ph, z) e_z
    sage: u[:]
    [u_rho(rh, ph, z), u_phi(rh, ph, z), u_z(rh, ph, z)]


Differential operators in cylindrical coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sage: from sage.manifolds.operators import *

The gradient::

    sage: F = E.scalar_field(function('f')(rh,ph,z), name='F')
    sage: F.display()
    F: E^3 → ℝ
       (rh, ph, z) ↦ f(rh, ph, z)
    sage: grad(F)
    Vector field grad(F) on the Euclidean space E^3
    sage: grad(F).display()
    grad(F) = d(f)/drh e_rh + d(f)/dph/rh e_ph + d(f)/dz e_z

The divergence::

    sage: s = div(u)
    sage: s.display()
    div(u): E^3 → ℝ
       (rh, ph, z) ↦ (rh*d(u_rho)/drh + rh*d(u_z)/dz + u_rho(rh, ph, z) + d(u_phi)/dph)/rh
    sage: s.expr().expand()
    u_rho(rh, ph, z)/rh + diff(u_phi(rh, ph, z), ph)/rh + diff(u_rho(rh, ph, z), rh)
     + diff(u_z(rh, ph, z), z)

The curl::

    sage: s = curl(u)
    sage: s
    Vector field curl(u) on the Euclidean space E^3
    sage: s.display()
    curl(u) = -(rh*d(u_phi)/dz - d(u_z)/dph)/rh e_rh + (d(u_rho)/dz - d(u_z)/drh) e_ph
     + (rh*d(u_phi)/drh + u_phi(rh, ph, z) - d(u_rho)/dph)/rh e_z

The Laplacian of a scalar field::

    sage: s = laplacian(F)
    sage: s.display()
    Delta(F): E^3 → ℝ
       (rh, ph, z) ↦ (rh^2*d^2(f)/drh^2 + rh^2*d^2(f)/dz^2 + rh*d(f)/drh
        + d^2(f)/dph^2)/rh^2
    sage: s.expr().expand()
    diff(f(rh, ph, z), rh)/rh + diff(f(rh, ph, z), ph, ph)/rh^2
     + diff(f(rh, ph, z), rh, rh) + diff(f(rh, ph, z), z, z)

The Laplacian of a vector field::

    sage: Du = laplacian(u)
    sage: Du.display()
    Delta(u) = (rh^2*d^2(u_rho)/drh^2 + rh^2*d^2(u_rho)/dz^2 + rh*d(u_rho)/drh
     - u_rho(rh, ph, z) - 2*d(u_phi)/dph + d^2(u_rho)/dph^2)/rh^2 e_rh
     + (rh^2*d^2(u_phi)/drh^2 + rh^2*d^2(u_phi)/dz^2 + rh*d(u_phi)/drh
     - u_phi(rh, ph, z) + d^2(u_phi)/dph^2 + 2*d(u_rho)/dph)/rh^2 e_ph
     + (rh^2*d^2(u_z)/drh^2 + rh^2*d^2(u_z)/dz^2 + rh*d(u_z)/drh
     + d^2(u_z)/dph^2)/rh^2 e_z

::

    sage: for i in E.irange():
    ....:     s = Du[i].expand()
    sage: Du.display_comp()
    Delta(u)^1 = d(u_rho)/drh/rh - u_rho(rh, ph, z)/rh^2 - 2*d(u_phi)/dph/rh^2
     + d^2(u_rho)/dph^2/rh^2 + d^2(u_rho)/drh^2 + d^2(u_rho)/dz^2
    Delta(u)^2 = d(u_phi)/drh/rh - u_phi(rh, ph, z)/rh^2 + d^2(u_phi)/dph^2/rh^2
     + 2*d(u_rho)/dph/rh^2 + d^2(u_phi)/drh^2 + d^2(u_phi)/dz^2
    Delta(u)^3 = d(u_z)/drh/rh + d^2(u_z)/dph^2/rh^2 + d^2(u_z)/drh^2 + d^2(u_z)/dz^2

Again, we may check that the above formulas coincide with those of Wikipedia's
article :wikipedia:`Del in cylindrical and spherical coordinates <Del_in_cylindrical_and_spherical_coordinates#Del_formula>`.

Changing coordinates
--------------------

Given the expression of a vector field in a given coordinate system, SageMath
can compute its expression in another coordinate system, see the tutorial
:ref:`change_coord_euclidean`
