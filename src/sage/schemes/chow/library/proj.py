# -*- coding: utf-8 -*-
r"""
The Proj ChowScheme

AUTHORS:

- Manfred Lehn (2013)

- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.chow.library.grass import Grass, GrassBundle


def Proj(n, hyperplane_class='h', names=None, name=None, latex_name=None):
    r"""
    Return the projective space `\mathbb{P}^n` of dimension n.
    INPUT:
    - ``n`` -- An integer, the dimension of the projective space.
    - ``hyperplane_class`` - An (optional) name for the hyperplane class
    - ``name`` -- An optional string, the name of the ChowScheme
    - ``latex_name``-- An optional string, the latex representation of the ChowScheme
    OUTPUT:
    - The ChowScheme corresponding to the projective space in the sense of
    Grothendieck, i.e. the rank 1 quotients of `mathbb{C}^{n+1}`.

    EXAMPLES::

        sage: X = Proj(3)  # P3 of rank 1 quotients of a 4 dim. vector space
        sage: X.sheaves["universal_sub"]
        Bundle(Proj(3, 'h'), 3, [1, -h, h^2, -h^3])
    """
    if name is None:
        name = "Proj(%s, '%s')" % (str(n), hyperplane_class)
    if latex_name is None:
        latex_name = r'\mathbb{P}^{%s}' % str(n)
    return Grass(n + 1, 1, chern_class=hyperplane_class, names=names,
                 name=name, latex_name=latex_name)


def ProjBundle(E, hyperplane_class='h',
               names=None, name=None, latex_name=None):
    r"""
    Return the *Proj* of a bundle E.
    INPUT:
    - ``E`` -- A sheaf on a ChowScheme
    - ``hyperplane_class`` - An (optional) name for the hyperplane class
    - ``name`` -- An optional string, the name of the ProjBundle
    - ``latex_name``-- An optional string, the latex representation of the ProjBundle
    OUTPUT:
    - The ChowScheme corresponding to `\mathbb{P}(E)` in the sense of
    Grothendieck, i.e. the rank 1 quotient modules of E.

    EXAMPLES::

        sage: P3 = Proj(3, name='P3')
        sage: S = P3.sheaves['universal_sub']
        sage: PS = ProjBundle(S); str(PS)
        'Proj(Bundle(P3, 3, [1, -h, h^2, -h^3]))'
    """
    if name is None:
        name = "Proj(%s)" % str(E)
    if latex_name is None:
        latex_name = r'\mathbb{P}^{%s}' % str(E)

    return GrassBundle(E, 1, chern_class=hyperplane_class, names=names,
                       name=name, latex_name=latex_name)
