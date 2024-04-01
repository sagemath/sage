.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

Advanced aspects: the Euclidean space as a Riemannian manifold
==============================================================

This tutorial introduces some vector calculus capabilities of SageMath within
the 3-dimensional Euclidean space. The corresponding tools have been developed
via the `SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/VectorCalculus/vector_calc_advanced.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/VectorCalculus/vector_calc_advanced.ipynb>`__ (``binder``).


The Euclidean 3-space
---------------------

Let us consider the 3-dimensional Euclidean space `\mathbb{E}^3`, with
Cartesian coordinates `(x,y,z)`::

    sage: E.<x,y,z> = EuclideanSpace()
    sage: E
    Euclidean space E^3

`\mathbb{E}^3` is actually a *Riemannian manifold* (see
:mod:`~sage.manifolds.differentiable.pseudo_riemannian`), i.e. a smooth real
manifold endowed with a positive definite metric tensor::

    sage: E.category()
    Join of
     Category of smooth manifolds over Real Field with 53 bits of precision and
     Category of connected manifolds over Real Field with 53 bits of precision and
     Category of complete metric spaces
    sage: E.base_field() is RR
    True
    sage: E.metric()
    Riemannian metric g on the Euclidean space E^3

Actually ``RR`` is used here as a proxy for the real field (this should be
replaced in the future, see the discussion at :issue:`24456`) and the 53 bits of precision play
of course no role for the symbolic computations.

Let us introduce spherical and cylindrical coordinates on
`\mathbb{E}^3`::

    sage: spherical.<r,th,ph> = E.spherical_coordinates()
    sage: cylindrical.<rh,ph,z> = E.cylindrical_coordinates()

The user atlas of `\mathbb{E}^3` has then three charts::

    sage: E.atlas()
    [Chart (E^3, (x, y, z)), Chart (E^3, (r, th, ph)), Chart (E^3, (rh, ph, z))]

while there are five vector frames defined on `\mathbb{E}^3`::

    sage: E.frames()
    [Coordinate frame (E^3, (e_x,e_y,e_z)),
     Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph)),
     Vector frame (E^3, (e_r,e_th,e_ph)),
     Coordinate frame (E^3, (∂/∂rh,∂/∂ph,∂/∂z)),
     Vector frame (E^3, (e_rh,e_ph,e_z))]

Indeed, there are two frames associated with each of the three coordinate
systems: the coordinate frame (denoted with partial derivatives above) and an
orthonormal frame (denoted by ``e_*`` above), but for Cartesian coordinates,
both frames coincide.

We get the orthonormal spherical and cylindrical frames by::

    sage: spherical_frame = E.spherical_frame()
    sage: spherical_frame
    Vector frame (E^3, (e_r,e_th,e_ph))
    sage: cylindrical_frame = E.cylindrical_frame()
    sage: cylindrical_frame
    Vector frame (E^3, (e_rh,e_ph,e_z))

On the other side, the coordinate frames `\left(\frac{\partial}{\partial r},
\frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)` and
`\left(\frac{\partial}{\partial \rho}, \frac{\partial}{\partial\phi},
\frac{\partial}{\partial z}\right)` are returned by the method
:meth:`~sage.manifolds.differentiable.chart.DiffChart.frame` acting on the
coordinate charts::

    sage: spherical.frame()
    Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph))
    sage: cylindrical.frame()
    Coordinate frame (E^3, (∂/∂rh,∂/∂ph,∂/∂z))


Charts as maps `\mathbb{E}^3 \rightarrow \mathbb{R}^3`
------------------------------------------------------------

The chart of Cartesian coordinates has been constructed at the
declaration of ``E``; let us denote it by ``cartesian``::

    sage: cartesian = E.cartesian_coordinates()
    sage: cartesian
    Chart (E^3, (x, y, z))

Let us consider a point `p\in \mathbb{E}^3`, defined by its
Cartesian coordinates::

    sage: p = E((-1, 1,0), chart=cartesian, name='p')
    sage: p
    Point p on the Euclidean space E^3
    sage: p.parent() is E
    True

The coordinates of `p` in a given coordinate chart are obtained by
letting the corresponding chart act on `p`::

    sage: cartesian(p)
    (-1, 1, 0)
    sage: spherical(p)
    (sqrt(2), 1/2*pi, 3/4*pi)
    sage: cylindrical(p)
    (sqrt(2), 3/4*pi, 0)

Riemannian metric
-----------------

The default metric tensor of `\mathbb{E}^3` is::

    sage: g = E.metric()
    sage: g
    Riemannian metric g on the Euclidean space E^3
    sage: g.display()
    g = dx⊗dx + dy⊗dy + dz⊗dz
    sage: g[:]
    [1 0 0]
    [0 1 0]
    [0 0 1]

The above display in performed in the default frame, which is the
Cartesian one. Of course, we may ask for display with respect to other
frames::

    sage: g.display(spherical_frame)
    g = e^r⊗e^r + e^th⊗e^th + e^ph⊗e^ph
    sage: g[spherical_frame, :]
    [1 0 0]
    [0 1 0]
    [0 0 1]

In the above display, ``e^r`` = `e^r`, ``e^th`` = `e^\theta` and
``e^ph`` = `e^\phi` are the 1-forms defining the coframe dual to the
orthonormal spherical frame `(e_r,e_\theta,e_\phi)`::

    sage: spherical_frame.coframe()
    Coframe (E^3, (e^r,e^th,e^ph))

The fact that the above metric components are either 0 or 1 reflect the
orthonormality of the vector frame `(e_r,e_\theta,e_\phi)`. On the
contrary, in the coordinate frame
`\left(\frac{\partial}{\partial r}, \frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)`,
which is not orthonormal, some components differ from 0 or 1::

    sage: g.display(spherical.frame())
    g = dr⊗dr + (x^2 + y^2 + z^2) dth⊗dth + (x^2 + y^2) dph⊗dph

Note that the components are expressed in terms of the default chart, namely
the Cartesian one. To have them displayed in terms of the spherical chart, we
have to provide the latter as the second argument of the method
``display()``::

    sage: g.display(spherical.frame(), spherical)
    g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph

Since SageMath 8.8, a shortcut is::

    sage: g.display(spherical)
    g = dr⊗dr + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph

The matrix view of the components is obtained via the square bracket operator::

    sage: g[spherical.frame(), :, spherical]
    [            1             0             0]
    [            0           r^2             0]
    [            0             0 r^2*sin(th)^2]

Similarly, for cylindrical coordinates, we have::

    sage: g.display(cylindrical_frame)
    g = e^rh⊗e^rh + e^ph⊗e^ph + e^z⊗e^z
    sage: g.display(cylindrical)
    g = drh⊗drh + rh^2 dph⊗dph + dz⊗dz
    sage: g[cylindrical.frame(), :, cylindrical]
    [   1    0    0]
    [   0 rh^2    0]
    [   0    0    1]

The metric `g` is a *flat*: its Riemann curvature tensor
(see
:meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.riemann`)
is zero::

    sage: g.riemann()
    Tensor field Riem(g) of type (1,3) on the Euclidean space E^3
    sage: g.riemann().display()
    Riem(g) = 0

The metric `g` defines the dot product on `\mathbb{E}^3`::

    sage: u = E.vector_field(x*y, y*z, z*x)
    sage: u.display()
    x*y e_x + y*z e_y + x*z e_z
    sage: v = E.vector_field(-y, x, z^2, name='v')
    sage: v.display()
    v = -y e_x + x e_y + z^2 e_z
    sage: u.dot(v) == g(u,v)
    True

Consequently::

    sage: norm(u) == sqrt(g(u,u))
    True


The Levi-Civita tensor
~~~~~~~~~~~~~~~~~~~~~~

The scalar triple product of `\mathbb{E}^3` is provided by the
Levi-Civita tensor (also called *volume form*) associated with `g`
(and chosen such that `(e_x,e_y,e_z)` is right-handed)::

    sage: epsilon = E.scalar_triple_product()
    sage: epsilon
    3-form epsilon on the Euclidean space E^3
    sage: epsilon is E.volume_form()
    True
    sage: epsilon.display()
    epsilon = dx∧dy∧dz
    sage: epsilon.display(spherical)
    epsilon = r^2*sin(th) dr∧dth∧dph
    sage: epsilon.display(cylindrical)
    epsilon = rh drh∧dph∧dz

Checking that all orthonormal frames introduced above are right-handed::

    sage: ex, ey, ez = E.cartesian_frame()[:]
    sage: epsilon(ex, ey, ez).display()
    epsilon(e_x,e_y,e_z): E^3 → ℝ
       (x, y, z) ↦ 1
       (r, th, ph) ↦ 1
       (rh, ph, z) ↦ 1

::

    sage: epsilon(*spherical_frame)
    Scalar field epsilon(e_r,e_th,e_ph) on the Euclidean space E^3
    sage: epsilon(*spherical_frame).display()
    epsilon(e_r,e_th,e_ph): E^3 → ℝ
       (x, y, z) ↦ 1
       (r, th, ph) ↦ 1
       (rh, ph, z) ↦ 1

::

    sage: epsilon(*cylindrical_frame).display()
    epsilon(e_rh,e_ph,e_z): E^3 → ℝ
       (x, y, z) ↦ 1
       (r, th, ph) ↦ 1
       (rh, ph, z) ↦ 1


Vector fields as derivations
----------------------------

Let `f` be a scalar field on `\mathbb{E}^3`::

    sage: f = E.scalar_field(x^2+y^2 - z^2, name='f')
    sage: f.display()
    f: E^3 → ℝ
       (x, y, z) ↦ x^2 + y^2 - z^2
       (r, th, ph) ↦ -2*r^2*cos(th)^2 + r^2
       (rh, ph, z) ↦ rh^2 - z^2

Vector fields act as derivations on scalar fields::

    sage: v(f)
    Scalar field v(f) on the Euclidean space E^3
    sage: v(f).display()
    v(f): E^3 → ℝ
       (x, y, z) ↦ -2*z^3
       (r, th, ph) ↦ -2*r^3*cos(th)^3
       (rh, ph, z) ↦ -2*z^3
    sage: v(f) == v.dot(f.gradient())
    True

::

    sage: df = f.differential()
    sage: df
    1-form df on the Euclidean space E^3
    sage: df.display()
    df = 2*x dx + 2*y dy - 2*z dz
    sage: v(f) == df(v)
    True

The algebra of scalar fields
----------------------------

The set `C^\infty(\mathbb{E}^3)` of all smooth scalar fields on
`\mathbb{E}^3` forms a commutative algebra over
`\mathbb{R}`::

    sage: CE = E.scalar_field_algebra()
    sage: CE
    Algebra of differentiable scalar fields on the Euclidean space E^3
    sage: CE.category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces
    sage: f in CE
    True

In SageMath terminology, `C^\infty(\mathbb{E}^3)` is the parent of scalar
fields::

    sage: f.parent() is CE
    True


The free module of vector fields
--------------------------------

The set `\mathfrak{X}(\mathbb{E}^3)` of all vector fields on `\mathbb{E}^3` is
a free module of rank 3 over the commutative algebra `C^\infty(\mathbb{E}^3)`::

    sage: XE = v.parent()
    sage: XE
    Free module X(E^3) of vector fields on the Euclidean space E^3
    sage: XE.category()
    Category of finite dimensional modules over Algebra of differentiable
     scalar fields on the Euclidean space E^3
    sage: XE.base_ring()
    Algebra of differentiable scalar fields on the Euclidean space E^3
    sage: XE.base_ring() is CE
    True
    sage: rank(XE)
    3

The bases of the free module `\mathfrak{X}(\mathbb{E}^3)` are nothing but the
vector frames defined on `\mathbb{E}^3`::

    sage: XE.bases()
    [Coordinate frame (E^3, (e_x,e_y,e_z)),
     Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph)),
     Vector frame (E^3, (e_r,e_th,e_ph)),
     Coordinate frame (E^3, (∂/∂rh,∂/∂ph,∂/∂z)),
     Vector frame (E^3, (e_rh,e_ph,e_z))]


Tangent spaces
--------------

A vector field evaluated at a point $p$ is a vector in the tangent space
`T_p\mathbb{E}^3`::

    sage: p
    Point p on the Euclidean space E^3
    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = -e_x - e_y

::

    sage: Tp = vp.parent()
    sage: Tp
    Tangent space at Point p on the Euclidean space E^3
    sage: Tp is E.tangent_space(p)
    True
    sage: Tp.category()
    Category of finite dimensional vector spaces over Symbolic Ring
    sage: dim(Tp)
    3
    sage: isinstance(Tp, FiniteRankFreeModule)
    True

The bases on `T_p\mathbb{E}^3` are inherited from the vector frames of
`\mathbb{E}^3`::

    sage: Tp.bases()
    [Basis (e_x,e_y,e_z) on the Tangent space at Point p on the Euclidean space E^3,
     Basis (∂/∂r,∂/∂th,∂/∂ph) on the Tangent space at Point p on the Euclidean space E^3,
     Basis (e_r,e_th,e_ph) on the Tangent space at Point p on the Euclidean space E^3,
     Basis (∂/∂rh,∂/∂ph,∂/∂z) on the Tangent space at Point p on the Euclidean space E^3,
     Basis (e_rh,e_ph,e_z) on the Tangent space at Point p on the Euclidean space E^3]

For instance, we have::

    sage: spherical_frame.at(p)
    Basis (e_r,e_th,e_ph) on the Tangent space at Point p on the
     Euclidean space E^3
    sage: spherical_frame.at(p) in Tp.bases()
    True


Levi-Civita connection
----------------------

The Levi-Civita connection associated to the Euclidean metric `g` is::

    sage: nabla = g.connection()
    sage: nabla
    Levi-Civita connection nabla_g associated with the Riemannian metric g
     on the Euclidean space E^3

The corresponding Christoffel symbols with respect to Cartesian coordinates
are identically zero: none of them appear in the output of
:meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.christoffel_symbols_display`,
which by default displays only nonzero Christoffel symbols::

    sage: g.christoffel_symbols_display(cartesian)

On the contrary, some of the Christoffel symbols with respect to
spherical coordinates differ from zero::

    sage: g.christoffel_symbols_display(spherical)
    Gam^r_th,th = -r
    Gam^r_ph,ph = -r*sin(th)^2
    Gam^th_r,th = 1/r
    Gam^th_ph,ph = -cos(th)*sin(th)
    Gam^ph_r,ph = 1/r
    Gam^ph_th,ph = cos(th)/sin(th)

By default, only nonzero and nonredundant values are displayed (for instance
`\Gamma^\phi_{\ \, \phi r}` is skipped, since it can be deduced from
`\Gamma^\phi_{\ \, r \phi}` by symmetry on the last two indices).

Similarly, the nonzero Christoffel symbols with respect to cylindrical
coordinates are::

    sage: g.christoffel_symbols_display(cylindrical)
    Gam^rh_ph,ph = -rh
    Gam^ph_rh,ph = 1/rh

The Christoffel symbols are nothing but the connection coefficients in the
corresponding coordinate frame::

    sage: nabla.display(cylindrical.frame(), cylindrical, only_nonredundant=True)
    Gam^rh_ph,ph = -rh
    Gam^ph_rh,ph = 1/rh

The connection coefficients with respect to the orthonormal
(non-coordinate) frames are (again only nonzero values are displayed)::

    sage: nabla.display(spherical_frame, spherical)
    Gam^1_22 = -1/r
    Gam^1_33 = -1/r
    Gam^2_12 = 1/r
    Gam^2_33 = -cos(th)/(r*sin(th))
    Gam^3_13 = 1/r
    Gam^3_23 = cos(th)/(r*sin(th))
    sage: nabla.display(cylindrical_frame, cylindrical)
    Gam^1_22 = -1/rh
    Gam^2_12 = 1/rh

The Levi-Civita connection `\nabla_g` is the connection involved in
the standard differential operators::

    sage: from sage.manifolds.operators import *
    sage: grad(f) == nabla(f).up(g)
    True
    sage: nabla(f) == grad(f).down(g)
    True
    sage: div(u) == nabla(u).trace()
    True
    sage: div(v) == nabla(v).trace()
    True
    sage: laplacian(f) == nabla(nabla(f).up(g)).trace()
    True
    sage: laplacian(u) == nabla(nabla(u).up(g)).trace(1,2)
    True
    sage: laplacian(v) == nabla(nabla(v).up(g)).trace(1,2)
    True
