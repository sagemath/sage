r"""
Degenerate submanifolds

An *embedded (resp. immersed) degenerate submanifold of a proper
pseudo-Riemannian manifold* `(M,g)` is an embedded (resp. immersed)
submanifold `H` of `M` as a differentiable manifold such that pull
back of the metric tensor `g` via the embedding (resp. immersion)
endows `H` with the structure of a degenerate manifold.

Degenerate submanifolds are study in many fields of mathematics and physics,
for instance in Differential Geometry (especially in geometry of
lightlike submanifold) and in General Relativity. In geometry of lightlike
submanifolds, according to the dimension `r` of the radical distribution
(see below for definition of radical distribution), degenerate submanifolds
have been classified into 4 subgroups: `r`-lightlike submanifolds, Coisotropic
submanifolds, Isotropic submanifolds and Totally lightlike submanifolds.
(See the book of Krishan L. Duggal and Aurel Bejancu [DS2010]_.)

In the present module, you can define any of the 4 types but most of the
methods are implemented only for degenerate hypersurfaces who belong to
`r`-lightlike submanifolds. However, they might be generalized to
`1`-lightlike submanifolds. In the literature there is a new approach
(the rigging technique) for studying `1`-lightlike submanifolds but
here we use the method of Krishan L. Duggal and Aurel Bejancu based on
the screen distribution.

Let `H` be a lightlike hypersurface of a pseudo-Riemannian manifold
`(M,g)`. Then the normal bundle `TH^\perp` intersect the tangent
bundle `TH`. The radical distribution is defined as
`Rad(TH)=TH\cap TH^\perp`. In case of hypersurfaces, and more
generally `1`-lightlike submanifolds, this is a rank 1 distribution.
A screen distribution `S(TH)` is a complementary of `Rad(TH)` in `TH`.

Giving a screen distribution `S(TH)` and a null vector field `\xi`
locally defined and spanning `Rad(TH)`, there exists a unique
transversal null vector field 'N' locally defined and such that
`g(N,\xi)=1`. From a transverse vector 'v', the null rigging 'N'
is giving by the formula

.. MATH::

    N = \frac{1}{g(\xi, v)}\left(v-\frac{g(v,v)}{2g(\xi, v)}\xi\right)

Tensors on the ambient manifold `M` are projected on `H` along `N`
to obtain induced objects. For instance, induced connection is the
linear connection defined on H through the Levi-Civitta connection of
`g` along `N`.

To work on a degenerate submanifold, after defining `H` as an instance
of :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`,
with the keyword ``structure='degenerate_metric'``, you have to set a
transvervector `v` and screen distribution together with the radical
distribution.

An example of degenerate submanifold from General Relativity is the
horizon of the Schwarzschild black hole. Allow us to recall that
Schwarzschild black hole is the first non-trivial solution of Einstein's
equations. It describes the metric inside a star of radius `R = 2m`,
being `m` the inertial mass of the star. It can be seen as an open
ball in a Lorentzian manifold structure on `\RR^4`::

    sage: M = Manifold(4, 'M', structure='Lorentzian')
    sage: X_M.<t, r, th, ph> = \
    ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi")
    sage: var('m'); assume(m>0)
    m
    sage: g = M.metric()
    sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
    ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2

Let us define the horizon as a degenerate hypersurface::

    sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
    sage: H
    degenerate hypersurface H embedded in 4-dimensional differentiable
    manifold M

A `2`-dimensional degenerate submanifold of a Lorentzian manifold::

    sage: M = Manifold(4, 'M', structure='Lorentzian')
    sage: X.<t,x,y,z> = M.chart()
    sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
    sage: S
    2-dimensional degenerate submanifold S embedded in 4-dimensional
    differentiable manifold M
    sage: X_S.<u,v> = S.chart()
    sage: Phi = S.diff_map(M, {(X_S, X): [u, u, u, v]},
    ....:         name='Phi', latex_name=r'\Phi');
    sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y]}, name='Phi_inv',
    ....:           latex_name=r'\Phi^{-1}');
    sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
    sage: g = M.metric()
    sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
    sage: S.set_transverse(rigging=[x,y])
    sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
    sage: V = M.vector_field(); V[3] = 1
    sage: Sc = S.screen('Sc', V, xi)

    sage: S.default_screen()
    screen distribution Sc along the 2-dimensional degenerate submanifold
    S embedded in 4-dimensional differentiable manifold M mapped into
    the 4-dimensional Lorentzian manifold M

    sage: S.ambient_metric()
    Lorentzian metric g on the 4-dimensional Lorentzian manifold M

    sage: S.induced_metric()
    degenerate metric gamma on the 2-dimensional degenerate submanifold S
    embedded in 4-dimensional differentiable manifold M

    sage: S.first_fundamental_form()
    Field of symmetric bilinear forms g|S along the 2-dimensional
    degenerate submanifold S embedded in 4-dimensional differentiable manifold M
    with values on the 4-dimensional Lorentzian manifold M

    sage: S.adapted_frame()
    Vector frame (S, (vv_0,vv_1,vv_2,vv_3)) with values on the 4-dimensional Lorentzian manifold M

    sage: S.projection(V)
    Tensor field of type (1,0) along the 2-dimensional degenerate submanifold S
    embedded in 4-dimensional differentiable manifold M
    with values on the 4-dimensional Lorentzian manifold M

    sage: S.weingarten_map()  # long time
    Tensor field nabla_g(xi)|X(S) of type (1,1) along the 2-dimensional
    degenerate submanifold S embedded in 4-dimensional differentiable manifold M
    with values on the 4-dimensional Lorentzian manifold M

    sage: SO = S.shape_operator()  # long time
    sage: SO.display()             # long time
    A^* = 0

    sage: S.is_tangent(xi.along(Phi))
    True
    sage: v = M.vector_field(); v[1] = 1
    sage: S.is_tangent(v.along(Phi))
    False

AUTHORS:

- Hans Fotsing Tetsing (2019) : initial version

REFERENCES:

- [DB1996]_
- [DS2010]_
- [FNO2019]_
"""
# *****************************************************************************
#  Copyright (C) 2019 Hans Fotsing Tetsing <hans.fotsing@aims-cameroon.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING

from sage.manifolds.differentiable.degenerate import DegenerateManifold, TangentTensor
from sage.manifolds.differentiable.differentiable_submanifold import (
    DifferentiableSubmanifold,
)
from sage.manifolds.differentiable.pseudo_riemannian import PseudoRiemannianManifold
from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
from sage.matrix.constructor import matrix
from sage.rings.infinity import infinity
from sage.symbolic.expression import Expression

if TYPE_CHECKING:
    from sage.manifolds.differentiable.metric import DegenerateMetric


class DegenerateSubmanifold(DegenerateManifold, DifferentiableSubmanifold):
    r"""
    Degenerate submanifolds.

    An *embedded (resp. immersed) degenerate submanifold of a proper
    pseudo-Riemannian manifold* `(M,g)` is an embedded (resp. immersed)
    submanifold `H` of `M` as a differentiable manifold such that pull
    back of the metric tensor `g` via the embedding (resp. immersion)
    endows `H` with the structure of a degenerate manifold.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``ambient`` -- (default: ``None``) pseudo-Riemannian manifold `M` in
      which the submanifold is embedded (or immersed). If ``None``, it is set
      to ``self``
    - ``metric_name`` -- (default: ``None``) string; name (symbol) given to the
      metric; if ``None``, ``'g'`` is used
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      tuple: `S = (n_+, n_-, n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the
      number of positive terms (resp. negative terms, resp. zero tems) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``base_manifold``
    - ``diff_degree`` -- (default: ``infinity``) degree of differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the metric; if none is provided, it is set to ``metric_name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g., coordinates
      in a chart
      - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`
      would return the previously constructed object corresponding to these
      arguments)

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold` and
        :mod:`~sage.manifolds.differentiable.differentiable_submanifold`
    """
    def __init__(self, n, name, ambient=None, metric_name=None, signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a degenerate submanifold.

        EXAMPLES:

        A `2`-dimensional degenerate submanifold of a Lorentzian manifold::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: S
            2-dimensional degenerate submanifold S embedded in 4-dimensional
            differentiable manifold M
        """
        DegenerateManifold.__init__(self, n, name=name,
                                          metric_name=metric_name,
                                          signature=signature,
                                          base_manifold=base_manifold,
                                          diff_degree=diff_degree,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        DifferentiableSubmanifold.__init__(self, n, name, self._field,
                                           self._structure, ambient=ambient,
                                           base_manifold=base_manifold,
                                           latex_name=latex_name,
                                           start_index=start_index,
                                           category=category)
        self._normal = None
        self._first_fundamental_form = None
        self._induced_metric = None
        self._second_fundamental_form = {}
        self._ambient_metric = None
        self._gauss_curvature = {}
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = {}
        self._screens = {}
        self._default_screen = None
        self._adapted_frame = {}
        self._shape_operator = {}
        self._rotation_one_form = {}
        signature = self._ambient.metric().signature()
        ndim = self._ambient._dim
        try:
            if signature[0] == ndim or signature[1] == ndim:
                raise ValueError("ambient must be a proper pseudo-Riemannian"
                                 " or a degenerate manifold")
        except TypeError:
            if signature == ndim or signature == -ndim:
                raise ValueError("ambient must be a proper pseudo-Riemannian"
                                 " or a degenerate manifold")
        self._transverse = {}

    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered
        as a manifold.

        TESTS:

        A `2`-dimensional degenerate submanifold of a Lorentzian manifold::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: S.__repr__()
            '2-dimensional degenerate submanifold S embedded in 4-dimensional
            differentiable manifold M'
        """
        if self._ambient is None:
            return super(DegenerateManifold, self).__repr__()
        if self._ambient._dim - self._dim == 1:
            return "degenerate hypersurface {} embedded " \
                "in {}-dimensional differentiable " \
                "manifold {}".format(self._name, self._ambient._dim,
                                     self._ambient._name)
        return "{}-dimensional degenerate submanifold {} embedded " \
               "in {}-dimensional differentiable " \
               "manifold {}".format(self._dim, self._name, self._ambient._dim,
                                    self._ambient._name)

    def ambient_metric(self):
        r"""
        Return the metric of the ambient manifold. The submanifold has to be
        embedded

        OUTPUT: the metric of the ambient manifold

        EXAMPLES:

        The lightcone of the 3D Minkowski space::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: X.<t,x,y> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v]},
            ....:           name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:               latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: S.ambient_metric()
            Lorentzian metric g on the 3-dimensional Lorentzian manifold M
        """
        if self._ambient_metric is None:
            if not self._embedded or not isinstance(self._ambient,
                            (PseudoRiemannianManifold, DegenerateManifold)):
                raise ValueError("degenerate submanifold must be embedded in a "
                                 "pseudo-Riemannian or degenerate manifold")
            self._ambient_metric = self._ambient.metric()
        return self._ambient_metric

    def default_screen(self):
        r"""
        Return the default screen distribution.

        OUTPUT:

        - an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi)  # long time
            sage: S.default_screen()              # long time
            screen distribution Sc along the degenerate hypersurface S embedded
            in 4-dimensional differentiable manifold M mapped into the 4-dimensional
            Lorentzian manifold M
        """
        return self._default_screen

    def list_of_screens(self):
        r"""
        Return the default screen distribution.

        OUTPUT:

        - an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi)  # long time
            sage: S.list_of_screens()             # long time
            {'Sc': screen distribution Sc along the degenerate hypersurface S
            embedded in 4-dimensional differentiable manifold M mapped into the
            4-dimensional Lorentzian manifold M}
        """
        return self._screens

    def set_transverse(self, rigging=None, normal=None):
        r"""
        For setting a transversal distribution of the degenerate submanifold.

        According to the type of the submanifold among the 4 possible types,
        one must enter a list of normal transversal vector fields and/or a
        list of transversal and not normal vector fields spanning a transverse
        distribution.

        INPUT:

        - ``rigging`` -- list or tuple (default: ``None``); list of vector fields
          of the ambient manifold or chart function; of the ambient manifold in
          the latter case, the corresponding gradient vector field with respect to
          the ambient metric is calculated; the vectors must be linearly independent,
          transversal to the submanifold but not normal
        - ``normal`` -- list or tuple (default: ``None``); list of vector fields
          of the ambient manifold or chart function; of the ambient manifold in
          the latter case, the corresponding gradient vector field with respect to
          the ambient metric is calculated; the vectors must be linearly independent,
          transversal and normal to the submanifold

        EXAMPLES:

        The lightcone of the 3-dimensional Minkowski space `\RR^3_1`::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: X.<t,x,y> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v]},
            ....:                   name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                           latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2] = -1,1,1
            sage: S.set_transverse(rigging=t)
        """
        if isinstance(rigging, (list, tuple)):
            rigging = list(rigging)
        else:
            if rigging is not None:
                rigging = [rigging]
        if isinstance(normal, (list, tuple)):
            normal = list(normal)
        else:
            if normal is not None:
                normal = [normal]
        frame = self.default_frame()
        im = self.immersion()
        g = self.ambient_metric().along(im)
        nor = []
        rig = []
        l1, l2 = 0, 0
        if normal is not None:
            for u in normal:
                if isinstance(u, Expression):
                    u = self._ambient.scalar_field(u).gradient()
                for v in frame:
                    v = im.pushforward(v)
                    if not g(u.along(im), v).is_zero():
                        raise ValueError("{} is not tangent to {}".format(u.display(), self._name))
                nor.append(u)
                l1 += 1
        if rigging is not None:
            for u in rigging:
                if isinstance(u, Expression):
                    u = self._ambient.scalar_field(u).gradient()
                rigg = False
                for v in frame:
                    v = im.pushforward(v)
                    if not g(u.along(im), v).is_zero():
                        rigg = True
                if not rigg:
                    raise ValueError("{} is normal to {}".format(u.display(), self._name))
                rig.append(u)
                l2 += 1
        if l1+l2 != self._codim:
            raise ValueError("length of the transverse must be {}".format(self._codim))
        self._transverse['normal'] = tuple(nor)
        self._transverse['rigging'] = tuple(rig)

    def screen(self, name, screen, rad, latex_name=None):
        r"""
        For setting a screen distribution and vector fields of the radical distribution
        that will be used for computations

        INPUT:

        - ``name`` -- string (default: ``None``); name given to the screen
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the screen; if ``None``, the LaTeX symbol is set to ``name``
        - ``screen`` -- list or tuple  of vector fields
          of the ambient manifold or chart function; of the ambient manifold in
          the latter case, the corresponding gradient vector field with respect to
          the ambient metric is calculated; the vectors must be linearly independent,
          tangent to the submanifold but not normal
        - ``rad`` -- -- list or tuple  of vector fields
          of the ambient manifold or chart function; of the ambient manifold in
          the latter case, the corresponding gradient vector field with respect to
          the ambient metric is calculated; the vectors must be linearly independent,
          tangent and normal to the submanifold

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi); Sc  # long time
            screen distribution Sc along the degenerate hypersurface S embedded
            in 4-dimensional differentiable manifold M mapped into the 4-dimensional
            Lorentzian manifold M
        """
        if isinstance(screen, (list, tuple)):
            screen = list(screen)
        else:
            screen = [screen]
        if isinstance(rad, (list, tuple)):
            rad = list(rad)
        else:
            rad = [rad]
        if name in self._screens:
            if list(screen) == self._screens[name]._screen and list(rad) == self._screens[name]._rad:
                return self._screens[name]
            else:
                raise ValueError("a different screen distribution with the "
                                 "same name had already been set")
        if len(screen)+len(rad) != self._dim:
            raise ValueError("total length screen+rad must be {}".format(self._dim))
        frame = self.default_frame()
        im = self.immersion()
        g = self.ambient_metric().along(im)
        for (i, u) in enumerate(screen):
            if isinstance(u, Expression):
                u = self._ambient.scalar_field(u).gradient()
            screen[i] = u
            sc = False
            for v in frame:
                v = im.pushforward(v)
                if not g(u.along(im), v).is_zero():
                    sc = True
            if not sc:
                raise ValueError("{} cannot belong to a screen distribution".format(u.display()))
        for (i, u) in enumerate(rad):
            if isinstance(u, Expression):
                u = self._ambient.scalar_field(u).gradient()
            rad[i] = u
            for v in frame:
                v = im.pushforward(v)
                if not g(u.along(im), v).is_zero():
                    raise ValueError("{} is not orthogonal to {}".format(u.display(), self._name))
        for u in self._transverse['rigging']:
            for v in screen:
                if not g(u.along(im), v.along(im)).is_zero():
                    raise ValueError("{} is not orthogonal to the rigging {}".format(v.display(), u.display()))
        self._screens[name] = Screen(self, name, tuple(screen), tuple(rad), latex_name=latex_name)
        self._default_screen = self._screens[name]
        return self._screens[name]

    def induced_metric(self) -> DegenerateMetric:
        r"""
        Return the pullback of the ambient metric.

        OUTPUT:

        - induced metric, as an instance of
          :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane
        passing through the origin::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]},
            ....:             name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: h = S.induced_metric(); h  # long time
            degenerate metric gamma on the 2-dimensional degenerate
            submanifold S embedded in 4-dimensional differentiable manifold M
        """
        if self._induced_metric is None or self._induced_metric._components == {}:
            self._induced_metric = self.metric()
            self._induced_metric.set(
                               self.immersion().pullback(self.ambient_metric()))
            self._induced_metric.set_name("gamma", r"\gamma")
        return self._induced_metric

    def first_fundamental_form(self):
        r"""
        Return the restriction of the ambient metric on vector field
        along the submanifold and tangent to it. It is difference from
        induced metric who gives the pullback of the ambient metric
        on the submanifold.

        OUTPUT:

        - the first fundamental form, as an instance of
          :class:`~sage.manifolds.differentiable.degenerate.TangentTensor`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: h = S.first_fundamental_form()   # long time
        """
        if self._first_fundamental_form is None:
            g = self.ambient_metric()
            h = g.copy()
            h.set_name(g._name+"|"+self._name, g._latex_name+r"|_"+self._latex_name)
            h = TangentTensor(h, self.immersion())
            self._first_fundamental_form = h
        return self._first_fundamental_form

    def _ambient_decomposition(self, screen=None):
        r"""
        Return a list ``[screen, rad, normal, rig]`` where ``screen``
        is a list a vector fields on the ambient manifold `M` spanning
        the giving screen distribution, ``rad`` a list of vector fields
        spanning radical distribution, ``normal`` list of normal transversal
        vector fields, and ``rig`` list of rigging vector fields
        corresponding to the giving screen.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`;
          if ``None``, the default screen is used.

        OUTPUT:

        - a list of 4 lists, this 1st one being independent vector fields
          spanning the screen distribution, the 2nd one independent vector fields
          spanning the radical distribution, the 3rd one independent vector fields
          spanning the transversal normal distribution, the 4th one being a list
          of independent riggings in `Rig(T\Sigma)` according to the decomposition

        .. MATH::

            TM_{|\Sigma}=S(T\Sigma)\oplus_{orth}((Rad(T\Sigma)\oplus_{orth}(
            T\sigma^\perp\cap tr(TM))\oplus Rig(T\Sigma))

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);   # long time
            sage: S._ambient_decomposition()        # long time
            [[Vector field on the 4-dimensional Lorentzian manifold M,
              Vector field on the 4-dimensional Lorentzian manifold M],
             [Vector field on the 4-dimensional Lorentzian manifold M],
             (),
             [Vector field N on the 4-dimensional Lorentzian manifold M]]
        """
        try:
            normal = self._transverse['normal']
        except KeyError:
            raise ValueError("Set first transverse by using set_transverse")
        if screen is None:
            screen = self.default_screen()
        if isinstance(screen, Screen):
            rad = screen._rad
            screen = screen._screen
            rig = self._transverse['rigging']
        else:
            raise ValueError("set first a screen distribution")
        if self._codim == 1:
            xi = rad[0]
            v = rig[0]
            g = self.ambient_metric()
            N = (1/g(xi, v))*(v-(g(v,v)/(2*g(xi, v)))*xi)
            if not self._adapted_frame:
                N.set_name(name='N')
            else:
                n = len(self._adapted_frame)
                N.set_name(name='N'+str(n))
            rig = [N]
        return [screen, rad, normal, rig]

    def _adapted_frame_(self, screen=None):
        r"""
        Return a frame
        `(e_1,\ldots,e_p, \xi_1,\ldots, \xi_r, v_1,\ldots, v_q, N_1, \ldots, N_n)`
        of the ambient manifold, being `e_i` vector fields
        spanning the giving screen distribution, `\xi_i` vector fields spanning
        radical distribution, `v_i` normal transversal vector fields, `N_i`
        rigging vector fields corresponding to the giving screen.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`;
          if ``None`` default screen is used.

        OUTPUT: a frame on the ambient manifold

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: T = S._adapted_frame_();         # long time
        """

        if screen is None:
            screen = self.default_screen()
        if screen is None:
            raise ValueError("set first a screen distribution")
        if screen._name in self._adapted_frame:
            return self._adapted_frame[screen._name]
        decomposition = self._ambient_decomposition(screen)
        sc, rad = decomposition[0], decomposition[1]
        normal, rigging = decomposition[2], decomposition[3]
        A = self._ambient.automorphism_field()
        i = self._ambient._sindex
        for u in sc:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in rad:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in normal:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in rigging:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        f = self._ambient.default_frame()
        GLHPhi = f.along(self.immersion())[0].parent().general_linear_group()
        if not self._adapted_frame:
            e = f.new_frame(A, 'vv')
        else:
            n = len(self._adapted_frame)
            e = f.new_frame(A, 'vv'+str(n))
        self.set_change_of_frame(f.along(self.immersion()), e.along(
                  self.immersion()), GLHPhi(A.along(self.immersion())))
        b = e.dual_basis()
        if self._codim == 1:
            if not self._adapted_frame:
                e[self._dim-self._sindex].set_name('N')
            else:
                n = len(self._adapted_frame)
                e[self._dim-self._sindex].set_name('N'+str(n))
            e[self._dim-self._sindex-1].set_name('xi', latex_name=r'\xi')
            if not self._adapted_frame:
                b[self._dim-self._sindex].set_name('N^b', latex_name=r'N^\flat')
            else:
                b[self._dim-self._sindex].set_name('N'+str(n)+'^b', latex_name=r'N'+str(n)+r'^\flat')
            b[self._dim-self._sindex-1].set_name('xi^b', latex_name=r'\xi^\flat')
        self._adapted_frame[screen._name] = e
        return e

    def adapted_frame(self, screen=None):
        r"""
        Return a frame
        `(e_1,\ldots,e_p, \xi_1,\ldots, \xi_r, v_1,\ldots, v_q, N_1, \ldots, N_n)`
        of the ambient manifold along the submanifold, being `e_i` vector fields
        spanning the giving screen distribution, `\xi_i` vector fields spanning
        radical distribution, `v_i` normal transversal vector fields, `N_i`
        rigging vector fields corresponding to the giving screen.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          if ``None`` default screen is used.

        OUTPUT: a frame on the ambient manifold along the submanifold

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: T = S.adapted_frame(); T         # long time
            Vector frame (S, (vv_0,vv_1,vv_2,vv_3)) with values on the 4-dimensional
            Lorentzian manifold M
        """
        e = self._adapted_frame_(screen).along(self.immersion())
        b = e.dual_basis()
        if self._codim == 1:
            if not self._adapted_frame:
                e[self._dim-self._sindex].set_name('N')
            else:
                n = len(self._adapted_frame)
                e[self._dim-self._sindex].set_name('N'+str(n))
            e[self._dim-self._sindex-1].set_name('xi', latex_name=r'\xi')
            if not self._adapted_frame:
                b[self._dim-self._sindex].set_name('N^b', latex_name=r'N^\flat')
            else:
                b[self._dim-self._sindex].set_name('N'+str(n)+'^b', latex_name=r'N'+str(n)+r'^\flat')
            b[self._dim-self._sindex-1].set_name('xi^b', latex_name=r'\xi^\flat')
        return self._adapted_frame_(screen).along(self.immersion())

    def second_fundamental_form(self, screen=None):
        r"""

        This method is implemented only for null hypersurfaces. The method
        returns a tensor `B` of type `(0,2)` instance of
        :class:`~sage.manifolds.differentiable.degenerate.TangentTensor`
        such that for two vector fields `U, V` on the ambient manifold along
        the null hypersurface, one has:

        .. MATH::

            \nabla_UV=D(U,V)+B(U,V)N

        being `\nabla` the ambient connection, `D` the induced connection
        and `N` the chosen rigging.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None``, the default screen is used

        OUTPUT:

        - an instance of
          :class:`~sage.manifolds.differentiable.degenerate.TangentTensor`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);   # long time
            sage: B = S.second_fundamental_form();  # long time
            sage: B.display()                       # long time
            B = 0
        """
        if self._ambient._dim-self._dim != 1:
            raise ValueError("'second_fundamental_form' is defined" +
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name not in self._second_fundamental_form:
            resu = self._ambient.vector_field_module() \
                .tensor((0, 2), name='B', latex_name='B', sym=[(0, 1)], antisym=[]) \
                .along(self.immersion())
            f = self.adapted_frame(screen)
            rad = self._ambient_decomposition(screen)[1]
            nab = self.ambient_metric().connection()
            xi = rad[0]
            A = self.screen_projection(nab(xi))
            resu[f, :] = A[f, :]
            im = self.immersion()
            self._second_fundamental_form[screen._name] = TangentTensor(resu, im)
        return self._second_fundamental_form[screen._name]

    extrinsic_curvature = second_fundamental_form

    def projection(self, tensor, screen=None):
        r"""

        For a given tensor `T` of type `(r, 1)` on the ambient manifold, this
        method returns the tensor `T'` of type `(r,1)` such that for `r`
        vector fields `v_1,\ldots,v_r`, `T'(v_1,\ldots,v_r)` is the projection
        of  `T(v_1,\ldots,v_r)` on ``self`` along the bundle spanned by the
        transversal vector fields provided by :meth:`set_transverse`.

        INPUT:

        - ``tensor`` -- a tensor of type `(r,1)` on the ambient manifold

        OUTPUT:

        - a tensor of type `(r,1)` on the ambient manifold along ``self``

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: U1 = S.projection(U)             # long time
        """
        if tensor.tensor_type()[0] != 1:
            raise NotImplementedError("``projection`` is implemented only for "
                                      "tensors with 1 as contravariant order")
        return TangentTensor(tensor, self.immersion(), screen)

    def screen_projection(self, tensor, screen=None):
        r"""
        For a given tensor `T` of type `(r, 1)` on the ambient manifold, this
        method returns the tensor `T'` of type `(r,1)` such that for `r`
        vector fields `v_1,\ldots,v_r`, `T'(v_1,\ldots,v_r)` is the projection
        of  `T(v_1,\ldots,v_r)` on the bundle spanned by ``screen`` along the
        bundle spanned by the transversal plus the radical vector fields provided.

        INPUT:

        - ``tensor`` -- a tensor of type `(r,1)` on the ambient manifold

        OUTPUT:

        - a tensor of type `(r,1)` on the ambient manifold

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: U1 = S.screen_projection(U);     # long time
        """
        if tensor.tensor_type()[0] != 1:
            raise NotImplementedError("``projection`` is implemented only for " +
                                      "tensors with 1 as contravariant order")
        frame = self.adapted_frame(screen)
        T = tensor.copy()
        try:
            T = T.along(self.immersion())
        except ValueError:
            pass
        T.display(frame)
        for i in self._ambient.index_generator(T.tensor_rank()):
            if i[0] in range(self._dim-self._sindex-1,self._ambient._dim-self._sindex):
                T[frame, i] = 0
        if tensor._latex_name is None:
            T.set_name(tensor._name)
        else:
            T.set_name("P"+tensor._name, latex_name=r'P'+tensor._latex_name)
        return TangentTensor(T, self.immersion(), screen)

    def weingarten_map(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces.
        *Weigarten map* is the `1`-form `W` defined for a vector field
        `U` tangent to ``self`` by

        .. MATH::

            W(U)=\nabla_U\xi

        being `\nabla` the Levi-Civita connection of the ambient manifold
        and `\xi` the chosen vector field spanning the radical distribution.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None`` the default screen is used.

        OUTPUT:

        - tensor of type `(1,1)` instance of
          :class:`~sage.manifolds.differentiable.degenerate.TangentTensor`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: W = S.weingarten_map();          # long time
            sage: W.display()                      # long time
            nabla_g(xi)|X(S) = 0
        """

        im = self.immersion()
        rad = self._ambient_decomposition(screen)[1]
        nab = self.ambient_metric().connection()
        xi = rad[0]
        T = self.projection(nab(xi)).extension()
        try:
            T = T.along(im)
        except ValueError:
            pass
        T.set_name("nabla_g(xi)|X("+self._name+")",
                   latex_name=r'\nabla_g(\xi)|_{\mathfrak{X}('+self._latex_name+r')}')
        return TangentTensor(T, im, screen)

    def shape_operator(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces.
        *shape operator* is the projection of the Weingarten map
        on the screen distribution along the radical distribution.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None`` the default screen is used.

        OUTPUT:

        - tensor of type `(1,1)` instance of
          :class:`~sage.manifolds.differentiable.degenerate.TangentTensor`

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: SO = S.shape_operator();         # long time
            sage: SO.display()                     # long time
            A^* = 0
        """
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._shape_operator:
            return self._shape_operator[screen._name]
        im = self.immersion()
        rad = self._ambient_decomposition(screen)[1]
        nab = self.ambient_metric().connection()
        xi = rad[0]
        T = self.screen_projection(nab(-xi), screen=screen).extension()
        try:
            T = T.along(im)
        except ValueError:
            pass
        T.set_name("A^*", latex_name=r'A^\ast')
        A = TangentTensor(T, im)
        self._shape_operator[screen._name] = A
        return A

    def gauss_curvature(self, screen=None):
        r"""
        Gauss curvature is the product of all eigenfunctions of the shape operator.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None`` the default screen is used.

        OUTPUT: a scalar function on ``self``

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: K = S.gauss_curvature();         # long time
            sage: K.display()                      # long time
            S → ℝ
            (u, v, w) ↦ 0
        """
        if self._ambient._dim-self._dim != 1:
            raise ValueError("'gauss_curvature' is defined"
                             " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name not in self._gauss_curvature:
            f = self.adapted_frame()
            A = self.shape_operator()
            self._gauss_curvature[screen._name] = self.scalar_field(
              {chart: A[f,:,chart].determinant()
                for chart in self.top_charts()})
        return self._gauss_curvature[screen._name]

    def principal_directions(self, screen=None):
        r"""

        Principal directions are eigenvectors of the shape operator. This
        method is implemented only for hypersurfaces.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None`` default screen is used.

        OUTPUT:

        - list of pairs (vector field, scalar field) representing the
          principal directions and the associated principal curvatures

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi); T = S.adapted_frame()  # long time
            sage: PD = S.principal_directions()                          # long time
            sage: PD[2][0].display(T)                                    # long time
            e_2 = xi
        """
        if self._codim != 1:
            raise ValueError("'principal directions' is defined" +
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._principal_directions:
            return self._principal_directions[screen._name]
        a = self.shape_operator(screen)
        frame = self.adapted_frame(screen)
        pr_d = matrix(
            [[a[frame, :][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvectors_right()
        res = []
        counter = self.irange()
        for eigen_space in pr_d:
            for eigen_vector in eigen_space[1]:
                v = self._ambient.vector_field(name="e_{}".format(next(counter))
                                                            ).along(self.immersion())
                v[frame, :] = list(eigen_vector) + [0]
                res.append((TangentTensor(v, self.immersion()), self.scalar_field(
                  {chart: eigen_space[0] for chart in self.top_charts()})))
                #res[-1][0].set_name("e_{}".format(next(counter)))
        self._principal_directions[screen._name] = res
        return res

    def mean_curvature(self, screen=None):
        r"""

        Mean curvature is the sum of principal curvatures. This
        method is implemented only for hypersurfaces.

        INPUT:

        - ``screen`` -- (default: ``None``) an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`.
          If ``None`` the default screen is used.

        OUTPUT: the mean curvature, as a scalar field on the submanifold

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: S.set_transverse(rigging=x)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: m = S.mean_curvature(); m        # long time
            Scalar field on the degenerate hypersurface S embedded in 4-dimensional
            differentiable manifold M
            sage: m.display()                      # long time
            S → ℝ
            (u, v, w) ↦ 0
        """
        if self._codim != 1:
            raise ValueError("'mean_curvature' is defined" +
                    " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._mean_curvature:
            return self._mean_curvature[screen._name]
        pc = [elt[-1] for elt in self.principal_directions(screen)]
        self._mean_curvature[screen._name] = self.scalar_field({chart: sum(
            pc).expr(chart)/self._dim for chart in self.top_charts()})
        return self._mean_curvature[screen._name]

    def is_tangent(self, v):
        r"""
        Determine whether a vector field on the ambient manifold along ``self``
        is tangent to ``self`` or not.

        INPUT:

        - ``v`` -- field on the ambient manifold along ``self``

        OUTPUT:

        - ``True`` if ``v`` is everywhere tangent to ``self`` or ``False`` if
          not

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);  # long time
            sage: S.is_tangent(xi.along(Phi))      # long time
            True
            sage: S.is_tangent(v.along(Phi))       # long time
            False
        """
        g = self.ambient_metric()
        im = self.immersion()
        decomposition = self._ambient_decomposition()
        rad, normal = decomposition[1], decomposition[2]
        for u in rad:
            if not g.along(im)(u.along(im),v).is_zero():
                return False
        for u in normal:
            if not g.along(im)(u.along(im),v).is_zero():
                return False
        return True


#**************************************************************************************

class Screen(VectorFieldModule):
    r"""
    Let `H` be a lightlike submanifold embedded in a pseudo-Riemannian
    manifold `(M,g)` with `\Phi` the embedding map. A screen distribution
    is a complementary `S(TH)` of the radical distribution `Rad(TM)=TH\cap
    TH^\perp` in `TH`. One then has

    .. MATH::

        TH=S(TH)\oplus_{orth}Rad(TH)

    INPUT:

    - ``submanifold`` -- a lightlike submanifold, as an instance of
      :class:`DegenerateSubmanifold`
    - ``name`` -- name given to the screen distribution
    - ``screen`` -- vector fields of the ambient manifold which
      span the screen distribution
    - ``rad`` -- vector fields of the ambient manifold which
      span the radical distribution
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      screen distribution; if ``None``, it is formed from ``name``

    EXAMPLES:

    The horizon of the Schwarzschild black hole::

        sage: M = Manifold(4, 'M', structure='Lorentzian')
        sage: X_M.<t, r, th, ph> = \
        ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi")
        sage: var('m'); assume(m>0)
        m
        sage: g = M.metric()
        sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
        ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
        sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
        sage: X_H.<ht,hth,hph> = \
        ....: H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta hph:(0,2*pi):\phi")
        sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, \
        ....: name='Phi', latex_name=r'\Phi')
        sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
        ....: name='Phi_inv', latex_name=r'\Phi^{-1}')
        sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
        sage: xi = M.vector_field(-1, 0, 0, 0)
        sage: v = M.vector_field(r, -r, 0, 0)
        sage: e1 = M.vector_field(0, 0, 1, 0)
        sage: e2 = M.vector_field(0, 0, 0, 1)

    A screen distribution for the Schwarzschild black hole horizon::

        sage: H.set_transverse(rigging=v)
        sage: S = H.screen('S', [e1, e2], (xi)); S  # long time
        screen distribution S along the degenerate hypersurface H embedded
        in 4-dimensional differentiable manifold M mapped into the
        4-dimensional Lorentzian manifold M

    The corresponding normal tangent null vector field and null
    transversal vector field::

        sage: xi = S.normal_tangent_vector(); xi.display()  # long time
        xi = -∂/∂t
        sage: N = S.rigging(); N.display()  # long time
        N = ∂/∂t - ∂/∂r

    Those vector fields are normalized by `g(\xi,N)=1`::

        sage: g.along(Phi)(xi, N).display()  # long time
        g(xi,N): H → ℝ
        (ht, hth, hph) ↦ 1

    """

    def __init__(self, submanifold, name, screen, rad, latex_name=None):
        r"""

        TESTS::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: X.<t,x,y> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2] = -1,1,1
            sage: S.set_transverse(rigging=t)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = -y; U[2] = x
            sage: Sc = S.screen('Sc', U, xi);
        """
        if not isinstance(submanifold, DegenerateSubmanifold):
            raise TypeError("the first argument must be a null submanifold")
        VectorFieldModule.__init__(self, submanifold, dest_map=submanifold.immersion())
        try:
            self._rigging = submanifold._transverse['rigging']
        except KeyError:
            raise ValueError("set first transverse by using `set_transverse` method")
        self._screen = list(screen)
        self._rad = list(rad)
        self._name = name
        if not latex_name:
            latex_name = name
        self._latex_name = latex_name

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: X.<t,x,y> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2] = -1,1,1
            sage: S.set_transverse(rigging=t)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = -y; U[2] = x
            sage: Sc = S.screen('Sc', U, xi); Sc._repr_()
            'screen distribution Sc along the degenerate hypersurface S embedded in
            3-dimensional differentiable manifold M mapped into the 3-dimensional
            Lorentzian manifold M'
        """
        description = "screen distribution "+self._name
        if self._dest_map is self._domain.identity_map():
            description += " on the {}".format(self._domain)
        else:
            description += (" along the {}".format(self._domain)
                            + " mapped into the {}".format(self._ambient_domain))
        return description

    def __getitem__(self, i):
        r"""
        Access vector fields spanning the screen distribution.

        INPUT:

        - ``i`` -- index of the coordinate; if the slice ``[:]``, then all
          the coordinates are returned

        OUTPUT:

        - a vector field on the ambient manifold that belong to
          the screen distribution

        TESTS::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: X.<t,x,y> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2] = -1,1,1
            sage: S.set_transverse(rigging=t)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = -y; U[2] = x
            sage: Sc = S.screen('Sc', U, xi);
            sage: Sc.__getitem__(0)
            Vector field along the degenerate hypersurface S embedded in
            3-dimensional differentiable manifold M with values on the 3-dimensional
            Lorentzian manifold M
        """
        sc = [elt.along(self._domain.immersion()) for elt in self._screen]
        return sc[i-self._domain._sindex]

    def normal_tangent_vector(self):
        r"""
        Return either a list ``Rad`` of vector fields spanning the radical
        distribution or (in case of a hypersurface) a normal tangent null
        vector field spanning the radical distribution.

        OUTPUT:

        - either a list of vector fields or a single vector field in
          case of a hypersurface

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:                  name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);                  # long time
            sage: Rad = Sc.normal_tangent_vector(); Rad.display()  # long time
            xi = ∂/∂t + ∂/∂x
        """
        rad = [elt.along(self._domain.immersion()) for elt in self._rad]
        if self._domain._codim == 1:
            xi = rad[0]
            xi.set_name(name='xi', latex_name=r'\xi')
            return xi
        else:
            return rad

    def rigging(self):
        r"""
        Return either a list ``Rad`` of vector fields spanning the
        complementary of the normal distribution `TH^\perp` in the
        transverse bundle or (when `H` is a null hypersurface) the
        null transversal vector field defined in [DB1996]_.

        OUTPUT:

        - either a list made by vector fields or a vector field in
          case of hypersurface

        EXAMPLES:

        A degenerate hyperplane the 4-dimensional Minkowski space `\RR^4_1`::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);    # long time
            sage: rig = Sc.rigging(); rig.display()  # long time
            N = -1/2 ∂/∂t + 1/2 ∂/∂x
        """
        im = self._domain.immersion()
        rig = [elt.along(im) for elt in self._domain._transverse['rigging']]
        if self._domain._codim != 1:
            return rig
        xi = self.normal_tangent_vector()
        v = rig[0]
        g = self._domain.ambient_metric().along(im)
        N = (1/g(xi, v))*(v-(g(v,v)/(2*g(xi, v)))*xi)
        N.set_name(name='N')
        return N
