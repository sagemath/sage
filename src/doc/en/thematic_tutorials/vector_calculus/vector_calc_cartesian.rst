.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

How to compute a gradient, a divergence or a curl
=================================================

This tutorial introduces some vector calculus capabilities of SageMath within
the 3-dimensional Euclidean space. The corresponding tools have been developed
via the `SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/VectorCalculus/vector_calc_cartesian.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/VectorCalculus/vector_calc_cartesian.ipynb>`__ (``binder``).


First stage: introduce the Euclidean 3-space
--------------------------------------------

Before evaluating some vector-field operators, one needs to define the arena
in which vector fields live, namely the *3-dimensional Euclidean space*
`\mathbb{E}^3`. In SageMath, we declare it, along with the standard *Cartesian
coordinates* `(x,y,z)`, via :class:`EuclideanSpace`::

    sage: E.<x,y,z> = EuclideanSpace()
    sage: E
    Euclidean space E^3

Thanks to the notation ``<x,y,z>`` in the above declaration, the coordinates
`(x,y,z)` are immediately available as three symbolic variables ``x``,
``y`` and ``z`` (there is no need to declare them via :func:`var`, i.e. to type
``x, y, z = var('x y z')``)::

    sage: x is E.cartesian_coordinates()[1]
    True
    sage: y is E.cartesian_coordinates()[2]
    True
    sage: z is E.cartesian_coordinates()[3]
    True
    sage: type(z)
    <class 'sage.symbolic.expression.Expression'>

Besides, `\mathbb{E}^3` is endowed with the *orthonormal vector frame*
`(e_x, e_y, e_z)` associated with Cartesian coordinates::

    sage: E.frames()
    [Coordinate frame (E^3, (e_x,e_y,e_z))]

At this stage, this is the *default* vector frame on `\mathbb{E}^3`
(being the only vector frame introduced so far)::

    sage: E.default_frame()
    Coordinate frame (E^3, (e_x,e_y,e_z))


Defining a vector field
-----------------------

We define a vector field on `\mathbb{E}^3` from its components in
the vector frame `(e_x,e_y,e_z)`::

    sage: v = E.vector_field(-y, x, sin(x*y*z), name='v')
    sage: v.display()
    v = -y e_x + x e_y + sin(x*y*z) e_z

We can access to the components of `v` via the square bracket operator::

    sage: v[1]
    -y
    sage: v[:]
    [-y, x, sin(x*y*z)]

A vector field can evaluated at any point of `\mathbb{E}^3`::

    sage: p = E((3,-2,1), name='p')
    sage: p
    Point p on the Euclidean space E^3
    sage: p.coordinates()
    (3, -2, 1)
    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = 2 e_x + 3 e_y - sin(6) e_z

Vector fields can be plotted::

    sage: v.plot(max_range=1.5, scale=0.5)
    Graphics3d Object


.. PLOT::

    E = EuclideanSpace(3)
    x, y, z = E.default_chart()[:]
    v = E.vector_field(-y, x, sin(x*y*z), name='v')
    sphinx_plot(v.plot(max_range=1.5, scale=0.5))

For customizing the plot, see the list of options in the documentation of
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.plot`.
For instance, to get a view of the orthogonal projection of `v` in the plane
`y=1`, do::

    sage: v.plot(fixed_coords={y: 1}, ambient_coords=(x,z), max_range=1.5,
    ....:        scale=0.25)
    Graphics object consisting of 81 graphics primitives

.. PLOT::

    E = EuclideanSpace(3)
    x, y, z = E.default_chart()[:]
    v = E.vector_field(-y, x, sin(x*y*z), name='v')
    g = v.plot(fixed_coords={y: 1}, ambient_coords=(x,z), max_range=1.5,
               scale=0.25)
    sphinx_plot(g)

We may define a vector field `u` with generic components `(u_x, u_y, y_z)`::

    sage: u = E.vector_field(function('u_x')(x,y,z),
    ....:                    function('u_y')(x,y,z),
    ....:                    function('u_z')(x,y,z),
    ....:                    name='u')
    sage: u.display()
    u = u_x(x, y, z) e_x + u_y(x, y, z) e_y + u_z(x, y, z) e_z
    sage: u[:]
    [u_x(x, y, z), u_y(x, y, z), u_z(x, y, z)]

Its value at the point `p` is then::

    sage: up = u.at(p)
    sage: up.display()
    u = u_x(3, -2, 1) e_x + u_y(3, -2, 1) e_y + u_z(3, -2, 1) e_z


How to compute various vector products
--------------------------------------

Dot product
~~~~~~~~~~~

The dot (or scalar) product `u\cdot v` of the vector fields `u`
and `v` is obtained by the method
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.dot_product`,
which admits ``dot()`` as a shortcut alias::

    sage: u.dot(v) == u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
    True

`s= u\cdot v` is a *scalar field*, i.e. a map `\mathbb{E}^3 \to \mathbb{R}`::

    sage: s = u.dot(v)
    sage: s
    Scalar field u.v on the Euclidean space E^3
    sage: s.display()
    u.v: E^3 → ℝ
       (x, y, z) ↦ -y*u_x(x, y, z) + x*u_y(x, y, z) + sin(x*y*z)*u_z(x, y, z)

It maps points of `\mathbb{E}^3` to real numbers::

    sage: s(p)
    -sin(6)*u_z(3, -2, 1) + 2*u_x(3, -2, 1) + 3*u_y(3, -2, 1)

Its coordinate expression is::

    sage: s.expr()
    -y*u_x(x, y, z) + x*u_y(x, y, z) + sin(x*y*z)*u_z(x, y, z)


Norm
~~~~

The norm `\|u\|` of the vector field `u` is defined in terms of the dot
product by `\|u\| = \sqrt{u\cdot u}`::

    sage: norm(u) == sqrt(u.dot(u))
    True

It is a scalar field on `\mathbb{E}^3`::

    sage: s = norm(u)
    sage: s
    Scalar field |u| on the Euclidean space E^3
    sage: s.display()
    |u|: E^3 → ℝ
       (x, y, z) ↦ sqrt(u_x(x, y, z)^2 + u_y(x, y, z)^2 + u_z(x, y, z)^2)
    sage: s.expr()
    sqrt(u_x(x, y, z)^2 + u_y(x, y, z)^2 + u_z(x, y, z)^2)

For `v`, we have::

    sage: norm(v).expr()
    sqrt(x^2 + y^2 + sin(x*y*z)^2)


Cross product
~~~~~~~~~~~~~

The cross product `u\times v` is obtained by the method
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.cross_product`,
which admits ``cross()`` as a shortcut alias::

    sage: s = u.cross(v)
    sage: s
    Vector field u x v on the Euclidean space E^3
    sage: s.display()
    u x v = (sin(x*y*z)*u_y(x, y, z) - x*u_z(x, y, z)) e_x
     + (-sin(x*y*z)*u_x(x, y, z) - y*u_z(x, y, z)) e_y
     + (x*u_x(x, y, z) + y*u_y(x, y, z)) e_z

We can check the standard formulas expressing the cross product in terms of
the components::

    sage: all([s[1] == u[2]*v[3] - u[3]*v[2],
    ....:      s[2] == u[3]*v[1] - u[1]*v[3],
    ....:      s[3] == u[1]*v[2] - u[2]*v[1]])
    True


Scalar triple product
~~~~~~~~~~~~~~~~~~~~~

Let us introduce a third vector field, `w` say. As a example, we do not pass
the components as arguments of
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.vector_field`,
as we did for `u` and `v`; instead, we set them in a second stage, via the
square bracket operator, any unset component being assumed to be zero::

    sage: w = E.vector_field(name='w')
    sage: w[1] = x*z
    sage: w[2] = y*z
    sage: w.display()
    w = x*z e_x + y*z e_y

The scalar triple product of the vector fields `u`, `v` and `w` is obtained as
follows::

    sage: triple_product = E.scalar_triple_product()
    sage: s = triple_product(u, v, w)
    sage: s
    Scalar field epsilon(u,v,w) on the Euclidean space E^3
    sage: s.expr()
    -(y*u_x(x, y, z) - x*u_y(x, y, z))*z*sin(x*y*z) - (x^2*u_z(x, y, z)
     + y^2*u_z(x, y, z))*z

Let us check that the scalar triple product of `u`, `v` and `w` is
`u\cdot(v\times w)`::

    sage: s == u.dot(v.cross(w))
    True


How to evaluate the standard differential operators
---------------------------------------------------

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

We first introduce a scalar field, via its expression in terms of
Cartesian coordinates; in this example, we consider some unspecified
function of `(x,y,z)`::

    sage: F = E.scalar_field(function('f')(x,y,z), name='F')
    sage: F.display()
    F: E^3 → ℝ
       (x, y, z) ↦ f(x, y, z)

The value of `F` at a point::

    sage: F(p)
    f(3, -2, 1)

The gradient of `F`::

    sage: grad(F)
    Vector field grad(F) on the Euclidean space E^3
    sage: grad(F).display()
    grad(F) = d(f)/dx e_x + d(f)/dy e_y + d(f)/dz e_z
    sage: norm(grad(F)).display()
    |grad(F)|: E^3 → ℝ
       (x, y, z) ↦ sqrt((d(f)/dx)^2 + (d(f)/dy)^2 + (d(f)/dz)^2)


Divergence
~~~~~~~~~~

The divergence of the vector field `u`::

    sage: s = div(u)
    sage: s.display()
    div(u): E^3 → ℝ
       (x, y, z) ↦ d(u_x)/dx + d(u_y)/dy + d(u_z)/dz

For `v` and `w`, we have::

    sage: div(v).expr()
    x*y*cos(x*y*z)
    sage: div(w).expr()
    2*z

An identity valid for any scalar field `F` and any vector field `u`::

    sage: div(F*u) == F*div(u) + u.dot(grad(F))
    True


Curl
~~~~

The curl of the vector field `u`::

    sage: s = curl(u)
    sage: s
    Vector field curl(u) on the Euclidean space E^3
    sage: s.display()
    curl(u) = (-d(u_y)/dz + d(u_z)/dy) e_x + (d(u_x)/dz - d(u_z)/dx) e_y
     + (-d(u_x)/dy + d(u_y)/dx) e_z

To use the notation ``rot`` instead of ``curl``, simply do::

    sage: rot = curl

An alternative is::

    sage: from sage.manifolds.operators import curl as rot

We have then::

    sage: rot(u).display()
    curl(u) = (-d(u_y)/dz + d(u_z)/dy) e_x + (d(u_x)/dz - d(u_z)/dx) e_y
     + (-d(u_x)/dy + d(u_y)/dx) e_z
    sage: rot(u) == curl(u)
    True

For `v` and `w`, we have::

    sage: curl(v).display()
    curl(v) = x*z*cos(x*y*z) e_x - y*z*cos(x*y*z) e_y + 2 e_z

::

    sage: curl(w).display()
    curl(w) = -y e_x + x e_y

The curl of a gradient is always zero::

    sage: curl(grad(F)).display()
    curl(grad(F)) = 0

The divergence of a curl is always zero::

    sage: div(curl(u)).display()
    div(curl(u)): E^3 → ℝ
       (x, y, z) ↦ 0

An identity valid for any scalar field `F` and any vector field `u` is

.. MATH::

    \mathrm{curl}(Fu) = \mathrm{grad}\, F\times u + F\,  \mathrm{curl}\, u,

as we can check::

    sage: curl(F*u) == grad(F).cross(u) + F*curl(u)
    True


Laplacian
~~~~~~~~~

The Laplacian `\Delta F` of a scalar field `F` is another scalar field::

    sage: s = laplacian(F)
    sage: s.display()
    Delta(F): E^3 → ℝ
       (x, y, z) ↦ d^2(f)/dx^2 + d^2(f)/dy^2 + d^2(f)/dz^2

The following identity holds:

.. MATH::

    \Delta F = \mathrm{div}\left(\mathrm{grad}\, F\right),

as we can check::

    sage: laplacian(F) == div(grad(F))
    True

The Laplacian `\Delta u` of a vector field `u` is another vector field::

    sage: Du = laplacian(u)
    sage: Du
    Vector field Delta(u) on the Euclidean space E^3

whose components are::

    sage: Du.display()
    Delta(u) = (d^2(u_x)/dx^2 + d^2(u_x)/dy^2 + d^2(u_x)/dz^2) e_x
     + (d^2(u_y)/dx^2 + d^2(u_y)/dy^2 + d^2(u_y)/dz^2) e_y
     + (d^2(u_z)/dx^2 + d^2(u_z)/dy^2 + d^2(u_z)/dz^2) e_z

In the Cartesian frame, the components of `\Delta u` are nothing but the
(scalar) Laplacians of the components of `u`, as we can check::

    sage: e = E.cartesian_frame()
    sage: Du == sum(laplacian(u[[i]])*e[i] for i in E.irange())
    True

In the above formula, ``u[[i]]`` return the `i`-th component of `u` as a
scalar field, while ``u[i]`` would have returned the coordinate expression of
this scalar field; besides, ``e`` is the Cartesian frame::

    sage: e[:]
    (Vector field e_x on the Euclidean space E^3,
     Vector field e_y on the Euclidean space E^3,
     Vector field e_z on the Euclidean space E^3)

For the vector fields `v` and `w`, we have::

    sage: laplacian(v).display()
    Delta(v) = -(x^2*y^2 + (x^2 + y^2)*z^2)*sin(x*y*z) e_z
    sage: laplacian(w).display()
    Delta(w) = 0

We have::

    sage: curl(curl(u)).display()
    curl(curl(u)) = (-d^2(u_x)/dy^2 - d^2(u_x)/dz^2 + d^2(u_y)/dxdy
     + d^2(u_z)/dxdz) e_x + (d^2(u_x)/dxdy - d^2(u_y)/dx^2 - d^2(u_y)/dz^2
     + d^2(u_z)/dydz) e_y + (d^2(u_x)/dxdz + d^2(u_y)/dydz - d^2(u_z)/dx^2
     - d^2(u_z)/dy^2) e_z
    sage: grad(div(u)).display()
    grad(div(u)) = (d^2(u_x)/dx^2 + d^2(u_y)/dxdy + d^2(u_z)/dxdz) e_x
     + (d^2(u_x)/dxdy + d^2(u_y)/dy^2 + d^2(u_z)/dydz) e_y
     + (d^2(u_x)/dxdz + d^2(u_y)/dydz + d^2(u_z)/dz^2) e_z

A famous identity is

.. MATH::

    \mathrm{curl}\left(\mathrm{curl}\, u\right) =
    \mathrm{grad}\left(\mathrm{div}\, u\right) - \Delta u .

Let us check it::

    sage: curl(curl(u)) == grad(div(u)) - laplacian(u)
    True


How to customize various symbols
--------------------------------

Customizing the symbols of the orthonormal frame vectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the vectors of the orthonormal frame associated with Cartesian
coordinates are denoted `(e_x,e_y,e_z)`::

    sage: frame = E.cartesian_frame()
    sage: frame
    Coordinate frame (E^3, (e_x,e_y,e_z))

But this can be changed, thanks to the method
:meth:`~sage.manifolds.differentiable.vectorframe.VectorFrame.set_name`::

    sage: frame.set_name('a', indices=('x', 'y', 'z'))
    sage: frame
    Coordinate frame (E^3, (a_x,a_y,a_z))
    sage: v.display()
    v = -y a_x + x a_y + sin(x*y*z) a_z

::

    sage: frame.set_name(('hx', 'hy', 'hz'),
    ....:                latex_symbol=(r'\mathbf{\hat{x}}', r'\mathbf{\hat{y}}',
    ....:                              r'\mathbf{\hat{z}}'))
    sage: frame
    Coordinate frame (E^3, (hx,hy,hz))
    sage: v.display()
    v = -y hx + x hy + sin(x*y*z) hz


Customizing the coordinate symbols
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coordinates symbols are defined within the angle brackets ``<...>`` at the
construction of the Euclidean space. Above we did::

    sage: E.<x,y,z> = EuclideanSpace()

which resulted in the coordinate symbols `(x,y,z)` and in the corresponding
Python variables ``x``, ``y`` and ``z`` (SageMath symbolic expressions). To
use other symbols, for instance `(X,Y,Z)`, it suffices to create ``E`` as::

    sage: E.<X,Y,Z> = EuclideanSpace()

We have then::

    sage: E.atlas()
    [Chart (E^3, (X, Y, Z))]
    sage: E.cartesian_frame()
    Coordinate frame (E^3, (e_X,e_Y,e_Z))
    sage: v = E.vector_field(-Y, X, sin(X*Y*Z), name='v')
    sage: v.display()
    v = -Y e_X + X e_Y + sin(X*Y*Z) e_Z

By default the LaTeX symbols of the coordinate coincide with the letters given
within the angle brackets. But this can be adjusted through the optional
argument ``symbols`` of the function :class:`EuclideanSpace`, which has to be
a string, usually prefixed by *r* (for *raw* string, in order to allow for the
backslash character of LaTeX expressions). This string contains the coordinate
fields separated by a blank space; each field contains the coordinate’s text
symbol and possibly the coordinate’s LaTeX symbol (when the latter is
different from the text symbol), both symbols being separated by a colon
(``:``)::

    sage: E.<xi,et,ze> = EuclideanSpace(symbols=r"xi:\xi et:\eta ze:\zeta")
    sage: E.atlas()
    [Chart (E^3, (xi, et, ze))]
    sage: v = E.vector_field(-et, xi, sin(xi*et*ze), name='v')
    sage: v.display()
    v = -et e_xi + xi e_et + sin(et*xi*ze) e_ze
