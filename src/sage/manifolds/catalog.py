r"""
Manifolds Catalog

A catalog of manifolds to rapidly create various simple manifolds.

The current entries to the catalog are obtained by typing
``manifolds.<tab>``, where ``<tab>`` indicates pressing the tab key.
They are:

- :class:`~sage.manifolds.differentiable.examples.euclidean.EuclideanSpace`: Euclidean space
- :class:`~sage.manifolds.differentiable.examples.real_line.RealLine`: real line
- :class:`~sage.manifolds.differentiable.examples.real_line.OpenInterval`: open interval on the real line
- :class:`~sage.manifolds.differentiable.examples.sphere.Sphere`: sphere embedded in Euclidean space
- :func:`Torus`: torus embedded in Euclidean space
- :func:`Minkowski`: 4-dimensional Minkowski space
- :func:`Kerr`: Kerr spacetime

AUTHORS:

- Florentin Jaffredo (2018) : initial version
"""

# *****************************************************************************
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Lazy import from examples folders:
from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import('sage.manifolds.differentiable.examples.real_line', 'OpenInterval')
_lazy_import('sage.manifolds.differentiable.examples.real_line', 'RealLine')
_lazy_import('sage.manifolds.differentiable.examples.euclidean', 'EuclideanSpace')
_lazy_import('sage.manifolds.differentiable.examples.sphere', 'Sphere')

def Minkowski(positive_spacelike=True, names=None):
    """
    Generate a Minkowski space of dimension 4.

    By default the signature is set to `(- + + +)`, but can be changed to
    `(+ - - -)` by setting the optional argument ``positive_spacelike`` to
    ``False``. The shortcut operator ``.<,>`` can be used to
    specify the coordinates.

    INPUT:

    - ``positive_spacelike`` -- (default: ``True``) if ``False``, then
      the spacelike vectors yield a negative sign (i.e., the signature
      is `(+ - - - )`)
    - ``names`` -- (default: ``None``) name of the coordinates,
      automatically set by the shortcut operator

    OUTPUT:

    - Lorentzian manifold of dimension 4 with (flat) Minkowskian metric

    EXAMPLES::

        sage: M.<t, x, y, z> = manifolds.Minkowski()
        sage: M.metric()[:]
        [-1  0  0  0]
        [ 0  1  0  0]
        [ 0  0  1  0]
        [ 0  0  0  1]

        sage: M.<t, x, y, z> = manifolds.Minkowski(False)
        sage: M.metric()[:]
        [ 1  0  0  0]
        [ 0 -1  0  0]
        [ 0  0 -1  0]
        [ 0  0  0 -1]
    """
    from sage.manifolds.manifold import Manifold
    M = Manifold(4, 'M', structure='Lorentzian')
    if names is None:
        names = ("t", "x", "y", "z")
    C = M.chart(names=names)
    M._first_ngens = C._first_ngens

    g = M.metric('g')
    sgn = 1 if positive_spacelike else -1
    g[0,0] = -sgn
    g[1,1], g[2,2], g[3,3] = sgn, sgn, sgn
    return M

def Kerr(m=1, a=0, coordinates="BL", names=None):
    """
    Generate a Kerr spacetime.

    A Kerr spacetime is a 4 dimensional manifold describing a rotating black
    hole. Two coordinate systems are implemented: Boyer-Lindquist and Kerr
    (3+1 version).

    The shortcut operator ``.<,>`` can be used to specify the coordinates.

    INPUT:

    - ``m`` -- (default: ``1``) mass of the black hole in natural units
      (`c=1`, `G=1`)
    - ``a`` -- (default: ``0``) angular momentum in natural units; if set to
      ``0``, the resulting spacetime corresponds to a Schwarzschild black hole
    - ``coordinates`` -- (default: ``"BL"``) either ``"BL"`` for
      Boyer-Lindquist coordinates or ``"Kerr"`` for Kerr coordinates (3+1
      version)
    - ``names`` -- (default: ``None``) name of the coordinates,
      automatically set by the shortcut operator

    OUTPUT:

    - Lorentzian manifold

    EXAMPLES::

        sage: m, a = var('m, a')
        sage: K = manifolds.Kerr(m, a)
        sage: K
        4-dimensional Lorentzian manifold M
        sage: K.atlas()
        [Chart (M, (t, r, th, ph))]
        sage: K.metric().display()
        g = (2*m*r/(a^2*cos(th)^2 + r^2) - 1) dt⊗dt
         + 2*a*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2) dt⊗dph
         + (a^2*cos(th)^2 + r^2)/(a^2 - 2*m*r + r^2) dr⊗dr
         + (a^2*cos(th)^2 + r^2) dth⊗dth
         + 2*a*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2) dph⊗dt
         + (2*a^2*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2) + a^2 + r^2)*sin(th)^2 dph⊗dph

        sage: K.<t, r, th, ph> = manifolds.Kerr()
        sage: K
        4-dimensional Lorentzian manifold M
        sage: K.metric().display()
        g = (2/r - 1) dt⊗dt + r^2/(r^2 - 2*r) dr⊗dr
         + r^2 dth⊗dth + r^2*sin(th)^2 dph⊗dph
        sage: K.default_chart().coord_range()
        t: (-oo, +oo); r: (0, +oo); th: (0, pi); ph: [-pi, pi] (periodic)

        sage: m, a = var('m, a')
        sage: K.<t, r, th, ph> = manifolds.Kerr(m, a, coordinates="Kerr")
        sage: K
        4-dimensional Lorentzian manifold M
        sage: K.atlas()
        [Chart (M, (t, r, th, ph))]
        sage: K.metric().display()
        g = (2*m*r/(a^2*cos(th)^2 + r^2) - 1) dt⊗dt
         + 2*m*r/(a^2*cos(th)^2 + r^2) dt⊗dr
         - 2*a*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2) dt⊗dph
         + 2*m*r/(a^2*cos(th)^2 + r^2) dr⊗dt
         + (2*m*r/(a^2*cos(th)^2 + r^2) + 1) dr⊗dr
         - a*(2*m*r/(a^2*cos(th)^2 + r^2) + 1)*sin(th)^2 dr⊗dph
         + (a^2*cos(th)^2 + r^2) dth⊗dth
         - 2*a*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2) dph⊗dt
         - a*(2*m*r/(a^2*cos(th)^2 + r^2) + 1)*sin(th)^2 dph⊗dr
         + (2*a^2*m*r*sin(th)^2/(a^2*cos(th)^2 + r^2)
         + a^2 + r^2)*sin(th)^2 dph⊗dph
        sage: K.default_chart().coord_range()
        t: (-oo, +oo); r: (0, +oo); th: (0, pi); ph: [-pi, pi] (periodic)
    """
    from sage.misc.functional import sqrt
    from sage.functions.trig import cos, sin
    from sage.manifolds.manifold import Manifold
    M = Manifold(4, 'M', structure="Lorentzian")
    if coordinates == "Kerr":
        if names is None:
            names = (r't:(-oo,+oo)', r'r:(0,+oo)', r'th:(0,pi):\theta',
                     r'ph:(-pi,pi):periodic:\phi')
        else:
            names = (names[0]+r':(-oo,+oo)', names[1]+r':(0,+oo)',
                     names[2]+r':(0,pi):\theta',
                     names[3]+r':(-pi,pi):periodic:\phi')
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        t, r, th, ph = C[:]
        rho = sqrt(r**2+a**2*cos(th)**2)
        g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -(1-2*m*r/rho**2), 1+2*m*r/rho**2,\
                rho**2, (r**2+a**2+2*a**2*m*r*sin(th)**2/rho**2)*sin(th)**2
        g[0, 1] = 2*m*r/rho**2
        g[0, 3] = -2*a*m*r/rho**2*sin(th)**2
        g[1, 3] = -a*sin(th)**2*(1+2*m*r/rho**2)
        return M

    if coordinates == "BL":
        if names is None:
            names = (r't:(-oo,+oo)', r'r:(0,+oo)', r'th:(0,pi):\theta',
                     r'ph:(-pi,pi):periodic:\phi')
        else:
            names = (names[0]+r':(-oo,+oo)', names[1]+r':(0,+oo)',
                     names[2]+r':(0,pi):\theta',
                     names[3]+r':(-pi,pi):periodic:\phi')
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        t, r, th, ph = C[:]
        rho = sqrt(r**2+a**2*cos(th)**2)
        g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -(1-2*m*r/rho**2), \
            rho**2/(r**2-2*m*r+a**2), rho**2, \
            (r**2+a**2+2*m*r*a**2/rho**2*sin(th)**2)*sin(th)**2
        g[0, 3] = 2*m*r*a*sin(th)**2/rho**2
        return M

    raise NotImplementedError("coordinates system not implemented, see help"
                              " for details")

def Torus(R=2, r=1, names=None):
    """
    Generate a 2-dimensional torus embedded in Euclidean space.

    The shortcut operator ``.<,>`` can be used to specify the coordinates.

    INPUT:

    - ``R`` -- (default: ``2``) distance form the center to the
      center of the tube
    - ``r`` -- (default: ``1``) radius of the tube
    - ``names`` -- (default: ``None``) name of the coordinates,
      automatically set by the shortcut operator

    OUTPUT:

    - Riemannian manifold

    EXAMPLES::

        sage: T.<theta, phi> = manifolds.Torus(3, 1)
        sage: T
        2-dimensional Riemannian submanifold T embedded in the Euclidean
         space E^3
        sage: T.atlas()
        [Chart (T, (theta, phi))]
        sage: T.embedding().display()
        T → E^3
           (theta, phi) ↦ (X, Y, Z) = ((cos(theta) + 3)*cos(phi),
                                          (cos(theta) + 3)*sin(phi),
                                          sin(theta))
        sage: T.metric().display()
        gamma = dtheta⊗dtheta + (cos(theta)^2 + 6*cos(theta) + 9) dphi⊗dphi
    """
    from sage.functions.trig import cos, sin
    from sage.manifolds.manifold import Manifold
    from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
    E = EuclideanSpace(3, symbols='X Y Z')
    M = Manifold(2, 'T', ambient=E, structure="Riemannian")
    if names is None:
        names = ("th", "ph")
    names = tuple([names[i] + ":(-pi,pi):periodic" for i in range(2)])
    C = M.chart(names=names)
    M._first_ngens = C._first_ngens
    th, ph = C[:]
    coordfunc = [(R+r*cos(th))*cos(ph), (R+r*cos(th))*sin(ph), r*sin(th)]
    imm = M.diff_map(E, coordfunc)
    M.set_embedding(imm)
    M.induced_metric()
    return M

