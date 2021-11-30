r"""
Spheres smoothly embedded in Euclidean Space

Let `E^{n+1}` be a Euclidean space of dimension `n+1` and `c \in E^{n+1}`. An
`n`-sphere with radius `r` and centered at `c`, usually denoted by
`\mathbb{S}^n_r(c)`, smoothly embedded in the Euclidean space `E^{n+1}` is an
`n`-dimensional smooth manifold together with a smooth embedding

.. MATH::

    \iota \colon \mathbb{S}^n_r \to E^{n+1}

whose image consists of all points having the same Euclidean distance to the
fixed point `c`. If we choose Cartesian coordinates `(x_1, \ldots, x_{n+1})` on
`E^{n+1}` with `x(c)=0` then the above translates to

.. MATH::

    \iota(\mathbb{S}^n_r(c)) = \left\{ p \in E^{n+1} : \lVert x(p) \rVert = r \right\}.

This corresponds to the standard `n`-sphere of radius `r` centered at `c`.

AUTHORS:

- Michael Jung (2020): initial version

REFERENCES:

- \M. Berger: *Geometry I&II* [Ber1987]_, [Ber1987a]_
- \J. Lee: *Introduction to Smooth Manifolds* [Lee2013]_

EXAMPLES:

We start by defining a 2-sphere of unspecified radius `r`::

    sage: r = var('r')
    sage: S2_r = manifolds.Sphere(2, radius=r); S2_r
    2-sphere S^2_r of radius r smoothly embedded in the Euclidean space E^3

The embedding `\iota` is constructed from scratch and can be returned by the
following command::

    sage: i = S2_r.embedding(); i
    Differentiable map iota from the 2-sphere S^2_r of radius r smoothly
     embedded in the Euclidean space E^3 to the Euclidean space E^3
    sage: i.display()
    iota: S^2_r → E^3
    on A: (theta, phi) ↦ (x, y, z) = (r*cos(phi)*sin(theta),
                                         r*sin(phi)*sin(theta),
                                         r*cos(theta))

As a submanifold of a Riemannian manifold, namely the Euclidean space,
the 2-sphere admits an induced metric::

    sage: g = S2_r.induced_metric()
    sage: g.display()
    g = r^2 dtheta⊗dtheta + r^2*sin(theta)^2 dphi⊗dphi

The induced metric is also known as the *first fundamental form* (see
:meth:`~sage.manifolds.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.first_fundamental_form`)::

    sage: g is S2_r.first_fundamental_form()
    True

The *second fundamental form* encodes the extrinsic curvature of the
2-sphere as hypersurface of Euclidean space (see
:meth:`~sage.manifolds.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.second_fundamental_form`)::

    sage: K = S2_r.second_fundamental_form(); K
    Field of symmetric bilinear forms K on the 2-sphere S^2_r of radius r
     smoothly embedded in the Euclidean space E^3
    sage: K.display()
    K = r dtheta⊗dtheta + r*sin(theta)^2 dphi⊗dphi

One quantity that can be derived from the second fundamental form is the
Gaussian curvature::

    sage: K = S2_r.gauss_curvature()
    sage: K.display()
    S^2_r → ℝ
    on A: (theta, phi) ↦ r^(-2)

As we have seen, spherical coordinates are initialized by default. To
initialize stereographic coordinates retrospectively, we can use the following
command::

    sage: S2_r.stereographic_coordinates()
    Chart (S^2_r-{NP}, (y1, y2))

To get all charts corresponding to stereographic coordinates, we can use the
:meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.coordinate_charts`::

    sage: stereoN, stereoS = S2_r.coordinate_charts('stereographic')
    sage: stereoN, stereoS
    (Chart (S^2_r-{NP}, (y1, y2)), Chart (S^2_r-{SP}, (yp1, yp2)))

.. SEEALSO::

    See :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.stereographic_coordinates`
    and :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.spherical_coordinates`
    for details.

.. NOTE::

    Notice that the derived quantities such as the embedding as well as the
    first and second fundamental forms must be computed from scratch again
    when new coordinates have been initialized. That makes the usage of
    previously declared objects obsolete.

Consider now a 1-sphere with barycenter `(1,0)` in Cartesian coordinates::

    sage: E2 = EuclideanSpace(2)
    sage: c = E2.point((1,0), name='c')
    sage: S1c.<chi> = E2.sphere(center=c); S1c
    1-sphere S^1(c) of radius 1 smoothly embedded in the Euclidean plane
     E^2 centered at the Point c
    sage: S1c.spherical_coordinates()
    Chart (A, (chi,))

Get stereographic coordinates::

    sage: stereoN, stereoS = S1c.coordinate_charts('stereographic')
    sage: stereoN, stereoS
    (Chart (S^1(c)-{NP}, (y1,)), Chart (S^1(c)-{SP}, (yp1,)))

The embedding takes now the following form in all coordinates::

    sage: S1c.embedding().display()
    iota: S^1(c) → E^2
    on A: chi ↦ (x, y) = (cos(chi) + 1, sin(chi))
    on S^1(c)-{NP}: y1 ↦ (x, y) = (2*y1/(y1^2 + 1) + 1, (y1^2 - 1)/(y1^2 + 1))
    on S^1(c)-{SP}: yp1 ↦ (x, y) = (2*yp1/(yp1^2 + 1) + 1, -(yp1^2 - 1)/(yp1^2 + 1))

Since the sphere is a hypersurface, we can get a normal vector field by using
``normal``::

    sage: n = S1c.normal(); n
    Vector field n along the 1-sphere S^1(c) of radius 1 smoothly embedded in
     the Euclidean plane E^2 centered at the Point c with values on the
     Euclidean plane E^2
    sage: n.display()
    n = -cos(chi) e_x - sin(chi) e_y

Notice that this is just *one* normal field with arbitrary direction,
in this particular case `n` points inwards whereas `-n` points outwards.
However, the vector field `n` is indeed non-vanishing and hence the sphere
admits an orientation (as all spheres do)::

    sage: orient = S1c.orientation(); orient
    [Coordinate frame (S^1(c)-{SP}, (∂/∂yp1)), Vector frame (S^1(c)-{NP}, (f_1))]
    sage: f = orient[1]
    sage: f[1].display()
    f_1 = -∂/∂y1

Notice that the orientation is chosen is such a way that `(\iota_*(f_1), -n)`
is oriented in the ambient Euclidean space, i.e. the last entry is the normal
vector field pointing outwards. Henceforth, the manifold admits
a volume form::

    sage: g = S1c.induced_metric()
    sage: g.display()
    g = dchi⊗dchi
    sage: eps = g.volume_form()
    sage: eps.display()
    eps_g = -dchi

"""

from sage.manifolds.differentiable.pseudo_riemannian_submanifold import \
                                                     PseudoRiemannianSubmanifold
from sage.categories.metric_spaces import MetricSpaces
from sage.categories.manifolds import Manifolds
from sage.categories.topological_spaces import TopologicalSpaces
from sage.rings.real_mpfr import RR
from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace

class Sphere(PseudoRiemannianSubmanifold):
    r"""
    Sphere smoothly embedded in Euclidean Space.

    An `n`-sphere of radius `r`smoothly embedded in a Euclidean space `E^{n+1}`
    is a smooth `n`-dimensional manifold smoothly embedded into `E^{n+1}`,
    such that the embedding constitutes a standard `n`-sphere of radius `r`
    in that Euclidean space (possibly shifted by a point).

    - ``n`` -- positive integer representing dimension of the sphere
    - ``radius`` -- (default: ``1``) positive number that states the radius
      of the sphere
    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      sphere; if ``None``, the name will be set according to the input
      (see convention above)
    - ``ambient_space`` -- (default: ``None``) Euclidean space in which the
      sphere should be embedded; if ``None``, a new instance of Euclidean
      space is created
    - ``center`` -- (default: ``None``) the barycenter of the sphere as point of
      the ambient Euclidean space; if ``None`` the barycenter is set to the
      origin of the ambient space's standard Cartesian coordinates
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the space; if ``None``, it will be set according to the input
      (see convention above)
    - ``coordinates`` -- (default: ``'spherical'``) string describing the
      type of coordinates to be initialized at the sphere's creation; allowed
      values are

          - ``'spherical'`` spherical coordinates (see
            :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.spherical_coordinates`))
          - ``'stereographic'`` stereographic coordinates given by the
            stereographic projection (see
            :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.stereographic_coordinates`)

    - ``names`` -- (default: ``None``) must be a tuple containing
      the coordinate symbols (this guarantees the shortcut operator
      ``<,>`` to function); if ``None``, the usual conventions are used (see
      examples below for details)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    A 2-sphere embedded in Euclidean space::

        sage: S2 = manifolds.Sphere(2); S2
        2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3
        sage: latex(S2)
        \mathbb{S}^{2}

    The ambient Euclidean space is constructed incidentally::

        sage: S2.ambient()
        Euclidean space E^3

    Another call creates another sphere and hence another Euclidean space::

        sage: S2 is manifolds.Sphere(2)
        False
        sage: S2.ambient() is manifolds.Sphere(2).ambient()
        False

    By default, the barycenter is set to the coordinate origin of the
    standard Cartesian coordinates in the ambient Euclidean space::

        sage: c = S2.center(); c
        Point on the Euclidean space E^3
        sage: c.coord()
        (0, 0, 0)

    Each `n`-sphere is a compact manifold and a complete metric space::

        sage: S2.category()
        Join of Category of compact topological spaces and Category of smooth
         manifolds over Real Field with 53 bits of precision and Category of
         connected manifolds over Real Field with 53 bits of precision and
         Category of complete metric spaces

    If not stated otherwise, each `n`-sphere is automatically endowed with
    spherical coordinates::

        sage: S2.atlas()
        [Chart (A, (theta, phi))]
        sage: S2.default_chart()
        Chart (A, (theta, phi))
        sage: spher = S2.spherical_coordinates()
        sage: spher is S2.default_chart()
        True

    Notice that the spherical coordinates do not cover the whole sphere. To
    cover the entire sphere with charts, use stereographic coordinates instead::

        sage: stereoN, stereoS = S2.coordinate_charts('stereographic')
        sage: stereoN, stereoS
        (Chart (S^2-{NP}, (y1, y2)), Chart (S^2-{SP}, (yp1, yp2)))
        sage: list(S2.open_covers())
        [Set {S^2} of open subsets of the 2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3,
         Set {S^2-{NP}, S^2-{SP}} of open subsets of the 2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3]

    .. NOTE::

        Keep in mind that the initialization process of stereographic
        coordinates and their transition maps is computational complex in
        higher dimensions. Henceforth, high computation times are expected with
        increasing dimension.

    """
    @staticmethod
    def __classcall_private__(cls, n=None, radius=1, ambient_space=None,
                              center=None, name=None, latex_name=None,
                              coordinates='spherical', names=None,
                              unique_tag=None):
        r"""
        Determine the correct class to return based upon the input.

        TESTS:

        Each call gives a new instance::

            sage: S2 = manifolds.Sphere(2)
            sage: S2 is manifolds.Sphere(2)
            False

        The dimension can be determined using the ``<...>`` operator::

            sage: S.<x,y> = manifolds.Sphere(coordinates='stereographic'); S
            2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3
            sage: S._first_ngens(2)
            (x, y)

        """
        if n is None:
            if names is None:
                raise ValueError("either n or names must be specified")
            n = len(names)

        # Technical bit for UniqueRepresentation
        from sage.misc.prandom import getrandbits
        from time import time
        if unique_tag is None:
            unique_tag = getrandbits(128)*time()

        return super(cls, Sphere).__classcall__(cls, n, radius=radius,
                                           ambient_space=ambient_space,
                                           center=center,
                                           name=name, latex_name=latex_name,
                                           coordinates=coordinates, names=names,
                                           unique_tag=unique_tag)

    def __init__(self, n, radius=1, ambient_space=None, center=None, name=None,
                 latex_name=None, coordinates='spherical', names=None,
                 category=None, init_coord_methods=None, unique_tag=None):
        r"""
        Construct sphere smoothly embedded in Euclidean space.

        TESTS::

            sage: S2 = manifolds.Sphere(2); S2
            2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3
            sage: S2.metric()
            Riemannian metric g on the 2-sphere S^2 of radius 1 smoothly
             embedded in the Euclidean space E^3
            sage: TestSuite(S2).run()

        """
        # radius
        if radius <= 0:
            raise ValueError('radius must be greater than zero')
        # ambient space
        if ambient_space is None:
            ambient_space = EuclideanSpace(n+1)
        elif not isinstance(ambient_space, EuclideanSpace):
            raise TypeError("the argument 'ambient_space' must be a Euclidean space")
        elif ambient_space._dim != n+1:
            raise ValueError("Euclidean space must have dimension {}".format(n+1))
        if center is None:
            cart = ambient_space.cartesian_coordinates()
            c_coords = [0]*(n+1)
            center = ambient_space.point(c_coords, chart=cart)
        elif center not in ambient_space:
            raise ValueError('{} must be an element of {}'.format(center, ambient_space))
        if name is None:
            name = 'S^{}'.format(n)
            if radius != 1:
                name += r'_{}'.format(radius)
            if center._name:
                name += r'({})'.format(center._name)
            if latex_name is None:
                latex_name = r'\mathbb{S}^{' + str(n) + r'}'
                if radius != 1:
                    latex_name += r'_{{{}}}'.format(radius)
                if center._latex_name:
                    latex_name += r'({})'.format(center._latex_name)
        if category is None:
            category = Manifolds(RR).Smooth() & MetricSpaces().Complete() & \
                       TopologicalSpaces().Compact().Connected()
        # initialize
        PseudoRiemannianSubmanifold.__init__(self, n, name,
                                             ambient=ambient_space,
                                             signature=n, latex_name=latex_name,
                                             metric_name='g', start_index=1,
                                             category=category)
        # set attributes
        self._radius = radius
        self._center = center
        self._coordinates = {}  # established coordinates; values are lists
        self._init_coordinates = {'spherical': self._init_spherical,
                                  'stereographic': self._init_stereographic}
                                 # predefined coordinates
        if init_coord_methods:
            self._init_coordinates.update(init_coord_methods)
        if coordinates not in self._init_coordinates:
            raise ValueError('{} coordinates not available'.format(coordinates))
        # up here, the actual initialization begins:
        self._init_chart_domains()
        self._init_embedding()
        self._init_coordinates[coordinates](names)

    def _init_embedding(self):
        r"""
        Initialize the embedding into Euclidean space.

        TESTS::

            sage: S2 = manifolds.Sphere(2)
            sage: i = S2.embedding(); i
            Differentiable map iota from the 2-sphere S^2 of radius 1 smoothly
             embedded in the Euclidean space E^3 to the Euclidean space E^3
            sage: i.display()
            iota: S^2 → E^3
             on A: (theta, phi) ↦ (x, y, z) = (cos(phi)*sin(theta),
             sin(phi)*sin(theta), cos(theta))

        """
        name = 'iota'
        latex_name = r'\iota'
        iota = self.diff_map(self._ambient, name=name, latex_name=latex_name)
        self.set_embedding(iota)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: S2_3 = manifolds.Sphere(2, radius=3)
            sage: S2_3._repr_()
            '2-sphere S^2_3 of radius 3 smoothly embedded in the Euclidean space E^3'
            sage: S2_3  # indirect doctest
            2-sphere S^2_3 of radius 3 smoothly embedded in the Euclidean space E^3

        """
        s = "{}-sphere {} of radius {} smoothly embedded in " \
            "the {}".format(self._dim, self._name, self._radius, self._ambient)
        if self._center._name:
            s += ' centered at the Point {}'.format(self._center._name)
        return s

    def coordinate_charts(self, coord_name, names=None):
        r"""
        Return a list of all charts belonging to the coordinates ``coord_name``.

        INPUT:

        - ``coord_name`` --  string describing the type of coordinates
        - ``names`` -- (default: ``None``) must be a tuple containing
          the coordinate symbols for the first chart in the list; if
          ``None``, the standard convention is used

        EXAMPLES:

        Spherical coordinates on `S^1`::

            sage: S1 = manifolds.Sphere(1)
            sage: S1.coordinate_charts('spherical')
            [Chart (A, (phi,))]

        Stereographic coordinates on `S^1`::

            sage: stereo_charts = S1.coordinate_charts('stereographic', names=['a'])
            sage: stereo_charts
            [Chart (S^1-{NP}, (a,)), Chart (S^1-{SP}, (ap,))]

        """
        if coord_name not in self._coordinates:
            if coord_name not in self._init_coordinates:
                raise ValueError('{} coordinates not available'.format(coord_name))
            self._init_coordinates[coord_name](names)
        return list(self._coordinates[coord_name])

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates of the default chart.

        This is useful only for the use of Sage preparser::

            sage: preparse("S3.<x,y,z> = manifolds.Sphere(coordinates='stereographic')")
            "S3 = manifolds.Sphere(coordinates='stereographic', names=('x', 'y', 'z',));
             (x, y, z,) = S3._first_ngens(3)"

        TESTS::

            sage: S2 = manifolds.Sphere(2)
            sage: S2._first_ngens(2)
            (theta, phi)
            sage: S2.<u,v> = manifolds.Sphere(2)
            sage: S2._first_ngens(2)
            (u, v)

        """
        return self._def_chart[:]

    def _init_chart_domains(self):
        r"""
        Construct the chart domains on ``self``.

        TESTS::

            sage: S2 = manifolds.Sphere(2)
            sage: list(S2.open_covers())
            [Set {S^2} of open subsets of the 2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3,
             Set {S^2-{NP}, S^2-{SP}} of open subsets of the 2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3]
            sage: frozenset(S2.subsets())  # random
            frozenset({Euclidean 2-sphere S^2 of radius 1,
             Open subset A of the Euclidean 2-sphere S^2 of radius 1,
             Open subset S^2-{NP,SP} of the Euclidean 2-sphere S^2 of radius 1,
             Open subset S^2-{NP} of the Euclidean 2-sphere S^2 of radius 1,
             Open subset S^2-{SP} of the Euclidean 2-sphere S^2 of radius 1})

        """
        # without north pole:
        name = self._name + '-{NP}'
        latex_name = self._latex_name + r'\setminus\{\mathrm{NP}\}'
        self._stereoN_dom = self.open_subset(name, latex_name)
        # without south pole:
        name = self._name + '-{SP}'
        latex_name = self._latex_name + r'\setminus\{\mathrm{SP}\}'
        self._stereoS_dom = self.open_subset(name, latex_name)
        # intersection:
        int = self._stereoN_dom.intersection(self._stereoS_dom)
        int._name = self._name + '-{NP,SP}'
        int._latex_name = self._latex_name + \
                          r'\setminus\{\mathrm{NP}, \mathrm{SP}\}'
        # without half circle:
        self._spher_dom = int.open_subset('A')
        # declare union:
        self.declare_union(self._stereoN_dom, self._stereoS_dom)

    def _shift_coords(self, coordfunc, s='+'):
        r"""
        Shift the coordinates ``coordfunc`` given in Cartesian coordinates to
        the barycenter.

        INPUT:

        - ``coordfunc`` -- coordinate expression given in the standard
          Cartesian coordinates of the ambient Euclidean space
        - ``s`` -- either ``'+'`` or ``'-'`` depending on whether ``coordfunc``
          should be shifted from the coordinate origin to the center or from
          the center to the coordinate origin

        TESTS::

            sage: E3 = EuclideanSpace(3)
            sage: c = E3.point((1,2,3), name='c')
            sage: S2c = manifolds.Sphere(2, ambient_space=E3, center=c)
            sage: spher = S2c.spherical_coordinates()
            sage: cart = E3.cartesian_coordinates()
            sage: coordfunc = S2c.embedding().expr(spher, cart); coordfunc
            (cos(phi)*sin(theta) + 1, sin(phi)*sin(theta) + 2, cos(theta) + 3)
            sage: S2c._shift_coords(coordfunc, s='-')
            (cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta))

        """
        cart = self._ambient.cartesian_coordinates()
        c_coords = cart(self._center)
        if s == '+':
            res = [x + c for x, c in zip(coordfunc, c_coords)]
        elif s == '-':
            res = [x - c for x, c in zip(coordfunc, c_coords)]
        return tuple(res)

    def _init_spherical(self, names=None):
        r"""
        Construct the chart of spherical coordinates.

        TESTS:

        Spherical coordinates on a 2-sphere::

            sage: S2 = manifolds.Sphere(2)
            sage: S2.spherical_coordinates()
            Chart (A, (theta, phi))

        Spherical coordinates on a 1-sphere::

            sage: S1 = manifolds.Sphere(1, coordinates='stereographic')
            sage: spher = S1.spherical_coordinates(); spher  # create coords
            Chart (A, (phi,))
            sage: S1.atlas()
            [Chart (S^1-{NP}, (y1,)),
             Chart (S^1-{SP}, (yp1,)),
             Chart (S^1-{NP,SP}, (y1,)),
             Chart (S^1-{NP,SP}, (yp1,)),
             Chart (A, (phi,)),
             Chart (A, (y1,)),
             Chart (A, (yp1,))]

        """
        # speed-up via simplification method...
        self.set_simplify_function(lambda expr: expr.simplify_trig())

        # get domain...
        A = self._spher_dom

        # initialize coordinates...
        n = self._dim
        if names:
            # add interval:
            names = tuple([x + ':(0,pi)' for x in names[:-1]] +
                          [names[-1] + ':(-pi,pi):periodic'])
        else:
            if n == 1:
                names = ('phi:(-pi,pi):periodic',)
            elif n == 2:
                names = ('theta:(0,pi)', 'phi:(-pi,pi):periodic')
            elif n == 3:
                names = ('chi:(0,pi)', 'theta:(0,pi)', 'phi:(-pi,pi):periodic')
            else:
                names = tuple(["phi_{}:(0,pi)".format(i) for i in range(1,n)] +
                              ["phi_{}:(-pi,pi):periodic".format(n)])
        spher = A.chart(names=names)
        coord = spher[:]

        # manage embedding...
        from sage.misc.misc_c import prod
        from sage.functions.trig import cos, sin

        R = self._radius

        coordfunc = [R*cos(coord[n-1])*prod(sin(coord[i]) for i in range(n-1))]
        coordfunc += [R*prod(sin(coord[i]) for i in range(n))]
        for k in reversed(range(n-1)):
            c = R*cos(coord[k])*prod(sin(coord[i]) for i in range(k))
            coordfunc.append(c)
        cart = self._ambient.cartesian_coordinates()
        # shift coordinates to barycenter:
        coordfunc = self._shift_coords(coordfunc, s='+')
        # add expression to embedding:
        self._immersion.add_expr(spher, cart, coordfunc)
        self.clear_cache()  # clear cache of extrinsic information

        # finish process...
        self._coordinates['spherical'] = [spher]

        # adapt other coordinates...
        if 'stereographic' in self._coordinates:
            self._transition_spher_stereo()

        # reset simplification method...
        self.set_simplify_function('default')

    def stereographic_coordinates(self, pole='north', names=None):
        r"""
        Return stereographic coordinates given by the stereographic
        projection of ``self`` w.r.t. to a given pole.

        INPUT:

        - ``pole`` -- (default: ``'north'``) the pole determining the
          stereographic projection; possible options are ``'north'`` and
          ``'south'``
        - ``names`` -- (default: ``None``) must be a tuple containing
          the coordinate symbols (this guarantees the usage of the shortcut
          operator ``<,>``)

        OUTPUT:

        - the chart of stereographic coordinates w.r.t. to the given pole,
          as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        Let `\mathbb{S}^n_r(c)` be an `n`-sphere of radius `r` smoothly
        embedded in the Euclidean space `E^{n+1}` centered at `c \in E^{n+1}`.
        We denote the north pole of `\mathbb{S}^n_r(c)` by `\mathrm{NP}` and the
        south pole by `\mathrm{SP}`. These poles are uniquely determined by
        the requirement

        .. MATH::

            x(\iota(\mathrm{NP})) &= (0, \ldots, 0, r) + x(c), \\
            x(\iota(\mathrm{SP})) &= (0, \ldots, 0, -r) + x(c).

        The coordinates `(y_1, \ldots, y_n)` (`(y'_1, \ldots, y'_n)`
        respectively) define *stereographic coordinates* on `\mathbb{S}^n_r(c)`
        for the Cartesian coordinates `(x_1, \ldots, x_{n+1})` on `E^{n+1}`
        if they arise from the stereographic projection from
        `\iota(\mathrm{NP})` (`\iota(\mathrm{SP})`) to the hypersurface
        `x_n = x_n(c)`. In concrete formulas, this means:

        .. MATH::

            \left. x \circ \iota \right|_{\mathbb{S}^n_r(c) \setminus \{
            \mathrm{NP}\}}
            &= \left( \frac{2y_1r^2}{r^2+\sum^n_{i=1} y^2_i},
            \ldots, \frac{2y_nr^2}{r^2+\sum^n_{i=1} y^2_i},
            \frac{r\sum^n_{i=1} y^2_i-r^3}{r^2+\sum^n_{i=1} y^2_i} \right) +
            x(c), \\
            \left. x \circ \iota \right|_{\mathbb{S}^n_r(c) \setminus \{
            \mathrm{SP}\}} &= \left( \frac{2y'_1r^2}{r^2+\sum^n_{i=1} y'^2_i},
            \ldots, \frac{2y'_nr^2}{r^2+\sum^n_{i=1} y'^2_i},
            \frac{r^3 - r\sum^n_{i=1} y'^2_i}{r^2+\sum^n_{i=1} y'^2_i} \right) +
            x(c).

        EXAMPLES:

        Initialize a 1-sphere centered at `(1,0)` in the Euclidean plane
        using the shortcut operator::

            sage: E2 = EuclideanSpace(2)
            sage: c = E2.point((1,0), name='c')
            sage: S1.<a> = E2.sphere(center=c, coordinates='stereographic'); S1
            1-sphere S^1(c) of radius 1 smoothly embedded in the Euclidean plane
             E^2 centered at the Point c

        By default, the shortcut variables belong to the stereographic
        projection from the north pole::

            sage: S1.coordinate_charts('stereographic')
            [Chart (S^1(c)-{NP}, (a,)), Chart (S^1(c)-{SP}, (ap,))]
            sage: S1.embedding().display()
            iota: S^1(c) → E^2
            on S^1(c)-{NP}: a ↦ (x, y) = (2*a/(a^2 + 1) + 1, (a^2 - 1)/(a^2 + 1))
            on S^1(c)-{SP}: ap ↦ (x, y) = (2*ap/(ap^2 + 1) + 1, -(ap^2 - 1)/(ap^2 + 1))

        Initialize a 2-sphere from scratch::

            sage: S2 = manifolds.Sphere(2)
            sage: S2.atlas()
            [Chart (A, (theta, phi))]

        In the previous block, the stereographic coordinates have not been
        initialized. This happens subsequently with the invocation of
        ``stereographic_coordinates``::

            sage: stereoS.<u,v> = S2.stereographic_coordinates(pole='south')
            sage: S2.coordinate_charts('stereographic')
            [Chart (S^2-{NP}, (up, vp)), Chart (S^2-{SP}, (u, v))]

        If not specified by the user, the default coordinate names are given by
        `(y_1, \ldots, y_n)` and `(y'_1, \ldots, y'_n)` respectively::

            sage: S3 = manifolds.Sphere(3, coordinates='stereographic')
            sage: S3.stereographic_coordinates(pole='north')
            Chart (S^3-{NP}, (y1, y2, y3))
            sage: S3.stereographic_coordinates(pole='south')
            Chart (S^3-{SP}, (yp1, yp2, yp3))

        """
        coordinates = 'stereographic'
        if coordinates not in self._coordinates:
            self._init_coordinates[coordinates](names, default_pole=pole)
        if pole == 'north':
            return self._coordinates[coordinates][0]
        elif pole == 'south':
            return self._coordinates[coordinates][1]
        else:
            raise ValueError("pole must be 'north' or 'south'")

    def spherical_coordinates(self, names=None):
        r"""
        Return the spherical coordinates of ``self``.

        INPUT:

        - ``names`` -- (default: ``None``) must be a tuple containing
          the coordinate symbols (this guarantees the usage of the shortcut
          operator ``<,>``)

        OUTPUT:

        - the chart of spherical coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        Let `\mathbb{S}^n_r(c)` be an `n`-sphere of radius `r` smoothly
        embedded in the Euclidean space `E^{n+1}` centered at `c \in E^{n+1}`.
        We say that `(\varphi_1, \ldots, \varphi_n)` define *spherical
        coordinates* on the open subset `A \subset \mathbb{S}^n_r(c)` for the
        Cartesian coordinates `(x_1, \ldots, x_{n+1})` on `E^{n+1}` (not
        necessarily centered at `c`) if

        .. MATH::

            \begin{aligned}
                \left. x_1 \circ \iota \right|_{A} &= r \cos(\varphi_n)\sin(
                \varphi_{n-1}) \cdots
                \sin(\varphi_1)
                + x_1(c), \\
                \left. x_1 \circ \iota \right|_{A} &= r \sin(\varphi_n)\sin(
                \varphi_{n-1}) \cdots
                \sin(\varphi_1)
                + x_1(c), \\
                \left. x_2 \circ \iota \right|_{A} &= r \cos(\varphi_{
                n-1})\sin(\varphi_{n-2}) \cdots
                \sin(\varphi_1)
                + x_2(c), \\
                \left. x_3 \circ \iota \right|_{A} &= r \cos(\varphi_{
                n-2})\sin(\varphi_{n-3}) \cdots
                \sin(\varphi_1)
                + x_3(c), \\
                \vdots & \\
                \left. x_{n+1} \circ \iota \right|_{A} &= r \cos(\varphi_1) +
                x_{n+1}(c),
            \end{aligned}

        where `\varphi_i` has range `(0, \pi)` for `i=1, \ldots, n-1` and
        `\varphi_n` lies in `(-\pi, \pi)`. Notice that the above expressions
        together with the ranges of the `\varphi_i` fully determine the open
        set `A`.

        .. NOTE::

            Notice that our convention slightly differs from the one given on
            the :wikipedia:`N-sphere#Spherical_coordinates`. The definition
            above ensures that the conventions for the most common cases
            `n=1` and `n=2` are maintained.

        EXAMPLES:

        The spherical coordinates on a 2-sphere follow the common conventions::

            sage: S2 = manifolds.Sphere(2)
            sage: spher = S2.spherical_coordinates(); spher
            Chart (A, (theta, phi))

        The coordinate range of spherical coordinates::

            sage: spher.coord_range()
            theta: (0, pi); phi: [-pi, pi] (periodic)

        Spherical coordinates do not cover the 2-sphere entirely::

            sage: A = spher.domain(); A
            Open subset A of the 2-sphere S^2 of radius 1 smoothly embedded in
             the Euclidean space E^3

        The embedding of a 2-sphere in Euclidean space via spherical
        coordinates::

            sage: S2.embedding().display()
            iota: S^2 → E^3
             on A: (theta, phi) ↦ (x, y, z) =
                                     (cos(phi)*sin(theta),
                                      sin(phi)*sin(theta),
                                      cos(theta))

        Now, consider spherical coordinates on a 3-sphere::

            sage: S3 = manifolds.Sphere(3)
            sage: spher = S3.spherical_coordinates(); spher
            Chart (A, (chi, theta, phi))
            sage: S3.embedding().display()
            iota: S^3 → E^4
            on A: (chi, theta, phi) ↦ (x1, x2, x3, x4) =
                                         (cos(phi)*sin(chi)*sin(theta),
                                          sin(chi)*sin(phi)*sin(theta),
                                          cos(theta)*sin(chi),
                                          cos(chi))

        By convention, the last coordinate is periodic::

            sage: spher.coord_range()
            chi: (0, pi); theta: (0, pi); phi: [-pi, pi] (periodic)

        """
        coordinates = 'spherical'
        if coordinates not in self._coordinates:
            self._init_coordinates[coordinates](names)
        return self._coordinates[coordinates][0]

    def _init_stereographic(self, names, default_pole='north'):
        r"""
        Construct the charts of stereographic coordinates.

        TESTS:

        Stereographic coordinates on the 2-sphere::

            sage: S2.<x,y> = manifolds.Sphere(2, coordinates='stereographic')
            sage: S2.atlas()
            [Chart (S^2-{NP}, (x, y)),
             Chart (S^2-{SP}, (xp, yp)),
             Chart (S^2-{NP,SP}, (x, y)),
             Chart (S^2-{NP,SP}, (xp, yp))]

        Stereographic coordinates on the 1-sphere::

            sage: S1 = manifolds.Sphere(1)
            sage: stereoS.<x> = S1.stereographic_coordinates(pole='south')
            sage: S1.atlas()
            [Chart (A, (phi,)),
             Chart (S^1-{NP}, (xp,)),
             Chart (S^1-{SP}, (x,)),
             Chart (S^1-{NP,SP}, (xp,)),
             Chart (S^1-{NP,SP}, (x,)),
             Chart (A, (xp,)),
             Chart (A, (x,))]

        """
        # speed-up via simplification method...
        self.set_simplify_function(lambda expr: expr.simplify_rational())
        # TODO: More speed-up?

        # get domains...
        U = self._stereoN_dom
        V = self._stereoS_dom

        # initialize coordinates...
        symbols_N = ''
        symbols_S = ''
        if names:
            for x in names:
                if default_pole == 'north':
                    symbols_N += x + ' '
                    symbols_S += "{}p".format(x) + ":{}' ".format(x)
                elif default_pole == 'south':
                    symbols_S += x + ' '
                    symbols_N += "{}p".format(x) + ":{}' ".format(x)
        else:
            for i in self.irange():
                symbols_N += "y{}".format(i) + r":y_{" + str(i) + r"} "
                symbols_S += "yp{}".format(i) + r":y'_{" + str(i) + r"} "
        symbols_N = symbols_N[:-1]
        symbols_S = symbols_S[:-1]
        stereoN = U.chart(coordinates=symbols_N)
        stereoS = V.chart(coordinates=symbols_S)
        coordN = stereoN[:]
        coordS = stereoS[:]

        # predefine variables...
        r2_N = sum(y ** 2 for y in coordN)
        r2_S = sum(yp ** 2 for yp in coordS)
        R = self._radius
        R2 = R**2

        # define transition map...
        coordN_to_S = tuple(R*y/r2_N for y in coordN)
        coordS_to_N = tuple(R*yp/r2_S for yp in coordS)
        stereoN_to_S = stereoN.transition_map(stereoS, coordN_to_S,
                                              restrictions1=r2_N != 0,
                                              restrictions2=r2_S != 0)
        stereoN_to_S.set_inverse(*coordS_to_N, check=False)

        # manage embedding...
        coordfuncN = [2*y*R2 / (R2+r2_N) for y in coordN]
        coordfuncN += [(R*r2_N-R*R2)/(R2+r2_N)]
        coordfuncS = [2*yp*R2 / (R2+r2_S) for yp in coordS]
        coordfuncS += [(R*R2-R*r2_S)/(R2+r2_S)]
        cart = self._ambient.cartesian_coordinates()
        # shift coordinates to barycenter:
        coordfuncN = self._shift_coords(coordfuncN, s='+')
        coordfuncS = self._shift_coords(coordfuncS, s='+')
        # add expressions to embedding:
        self._immersion.add_expr(stereoN, cart, coordfuncN)
        self._immersion.add_expr(stereoS, cart, coordfuncS)
        self.clear_cache()  # clear cache of extrinsic information

        # define orientation...
        eN = stereoN.frame()
        eS = stereoS.frame()  # oriented w.r.t. embedding
        frame_comp = list(eN[:])
        frame_comp[0] = -frame_comp[0]  # reverse orientation
        f = U.vector_frame('f', frame_comp)
        self.set_orientation([eS, f])

        # finish process...
        self._coordinates['stereographic'] = [stereoN, stereoS]

        # adapt other coordinates...
        if 'spherical' in self._coordinates:
            self._transition_spher_stereo()

        # reset simplification method...
        self.set_simplify_function('default')

    def _transition_spher_stereo(self):
        r"""
        Initialize the transition map between spherical and stereographic
        coordinates.

        TESTS::

            sage: S1 = manifolds.Sphere(1)
            sage: spher = S1.spherical_coordinates(); spher
            Chart (A, (phi,))
            sage: A = spher.domain()
            sage: stereoN = S1.stereographic_coordinates(pole='north'); stereoN
            Chart (S^1-{NP}, (y1,))
            sage: S1.coord_change(spher, stereoN.restrict(A))
            Change of coordinates from Chart (A, (phi,)) to Chart (A, (y1,))

        """
        # speed-up via simplification method...
        self.set_simplify_function(lambda expr: expr.simplify())

        # configure preexisting charts...
        W = self._stereoN_dom.intersection(self._stereoS_dom)
        A = self._spher_dom
        stereoN, stereoS = self._coordinates['stereographic'][:]
        coordN = stereoN[:]
        coordS = stereoS[:]
        rstN = (coordN[0] != 0,)
        rstS = (coordS[0] != 0,)
        if self._dim > 1:
            rstN += (coordN[0] > 0,)
            rstS += (coordS[0] > 0,)
        stereoN_A = stereoN.restrict(A, rstN)
        stereoS_A = stereoS.restrict(A, rstS)
        self._coord_changes[(stereoN.restrict(W),
                             stereoS.restrict(W))].restrict(A)
        self._coord_changes[(stereoS.restrict(W),
                             stereoN.restrict(W))].restrict(A)
        spher = self._coordinates['spherical'][0]

        R = self._radius
        n = self._dim

        # transition: spher to stereoN...
        imm = self.embedding()
        cart = self._ambient.cartesian_coordinates()
        # get ambient coordinates and shift to coordinate origin:
        x = self._shift_coords(imm.expr(spher, cart), s='-')
        coordfunc = [(R*x[i])/(R-x[-1]) for i in range(n)]
        # define transition map:
        spher_to_stereoN = spher.transition_map(stereoN_A, coordfunc)

        # transition: stereoN to spher...
        from sage.functions.trig import acos, atan2
        from sage.misc.functional import sqrt
        # get ambient coordinates and shift to coordinate origin:
        x = self._shift_coords(imm.expr(stereoN, cart), s='-')
        coordfunc = [atan2(x[1],x[0])]
        for k in range(2, n+1):
            c = acos(x[k]/sqrt(sum(x[i]**2 for i in range(k+1))))
            coordfunc.append(c)
        coordfunc = reversed(coordfunc)
        spher_to_stereoN.set_inverse(*coordfunc, check=False)

        # transition spher <-> stereoS...
        stereoN_to_S_A = self.coord_change(stereoN_A, stereoS_A)
        stereoN_to_S_A * spher_to_stereoN  # generates spher_to_stereoS
        stereoS_to_N_A = self.coord_change(stereoS_A, stereoN_A)
        spher_to_stereoN.inverse() * stereoS_to_N_A  # generates stereoS_to_spher

    def dist(self, p, q):
        r"""
        Return the great circle distance between the points ``p`` and ``q`` on
        ``self``.

        INPUT:

        - ``p`` -- an element of ``self``
        - ``q`` -- an element of ``self``

        OUTPUT:

        - the great circle distance `d(p, q)` on ``self``

        The great circle distance `d(p, q)` of the points
        `p, q \in \mathbb{S}^n_r(c)` is the length of the shortest great circle
        segment on `\mathbb{S}^n_r(c)` that joins `p` and `q`. If we choose
        Cartesian coordinates `(x_1, \ldots, x_{n+1})` of the ambient Euclidean
        space such that the center lies in the coordinate origin, i.e.
        `x(c)=0`, the great circle distance can be expressed in terms of the
        following formula:

        .. MATH::

            d(p,q) = r \, \arccos\left(\frac{x(\iota(p)) \cdot
                x(\iota(q))}{r^2}\right).

        EXAMPLES:

        Define a 2-sphere with unspecified radius::

            sage: r = var('r')
            sage: S2_r = manifolds.Sphere(2, radius=r); S2_r
            2-sphere S^2_r of radius r smoothly embedded in the Euclidean space E^3

        Given two antipodal points in spherical coordinates::

            sage: p = S2_r.point((pi/2, pi/2), name='p'); p
            Point p on the 2-sphere S^2_r of radius r smoothly embedded in the
             Euclidean space E^3
            sage: q = S2_r.point((pi/2, -pi/2), name='q'); q
            Point q on the 2-sphere S^2_r of radius r smoothly embedded in the
             Euclidean space E^3

        The distance is determined as the length of the half great circle::

            sage: S2_r.dist(p, q)
            pi*r

        """
        from sage.functions.trig import acos
        # get Euclidean points:
        x = self._immersion(p)
        y = self._immersion(q)
        cart = self._ambient.cartesian_coordinates()
        # get ambient coordinates and shift to coordinate origin:
        x_coord = self._shift_coords(x.coord(chart=cart), s='-')
        y_coord = self._shift_coords(y.coord(chart=cart), s='-')

        n = self._dim + 1
        r = self._radius
        inv_angle = sum(x_coord[i]*y_coord[i] for i in range(n)) / r**2
        return (r * acos(inv_angle)).simplify()

    def radius(self):
        r"""
        Return the radius of ``self``.

        EXAMPLES:

        3-sphere with radius 3::

            sage: S3_2 = manifolds.Sphere(3, radius=2); S3_2
            3-sphere S^3_2 of radius 2 smoothly embedded in the 4-dimensional
             Euclidean space E^4
            sage: S3_2.radius()
            2

        2-sphere with unspecified radius::

            sage: r = var('r')
            sage: S2_r = manifolds.Sphere(3, radius=r); S2_r
            3-sphere S^3_r of radius r smoothly embedded in the 4-dimensional
             Euclidean space E^4
            sage: S2_r.radius()
            r

        """
        return self._radius

    def minimal_triangulation(self):
        r"""
        Return the minimal triangulation of ``self`` as a simplicial complex.

        EXAMPLES:

        Minimal triangulation of the 2-sphere::

            sage: S2 = manifolds.Sphere(2)
            sage: S = S2.minimal_triangulation(); S
            Minimal triangulation of the 2-sphere

        The Euler characteristic of a 2-sphere::

            sage: S.euler_characteristic()
            2

        """
        from sage.topology.simplicial_complex_examples import Sphere as SymplicialSphere
        return SymplicialSphere(self._dim)

    def center(self):
        r"""
        Return the barycenter of ``self`` in the ambient Euclidean space.

        EXAMPLES:

        2-sphere embedded in Euclidean space centered at `(1,2,3)` in
        Cartesian coordinates::

            sage: E3 = EuclideanSpace(3)
            sage: c = E3.point((1,2,3), name='c')
            sage: S2c = manifolds.Sphere(2, ambient_space=E3, center=c); S2c
            2-sphere S^2(c) of radius 1 smoothly embedded in the Euclidean space
             E^3 centered at the Point c
            sage: S2c.center()
            Point c on the Euclidean space E^3

        We can see that the embedding is shifted accordingly::

            sage: S2c.embedding().display()
            iota: S^2(c) → E^3
            on A: (theta, phi) ↦ (x, y, z) = (cos(phi)*sin(theta) + 1,
                                                 sin(phi)*sin(theta) + 2,
                                                 cos(theta) + 3)

        """
        return self._center
