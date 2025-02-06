# -*- coding: utf-8 -*-
r"""
The Grassmannian ChowScheme
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

from sage.schemes.chow.bundle import Bundle, BundleDiffRelations
from sage.schemes.chow.scheme import PointChowScheme
from sage.schemes.chow.sheaf import SHom, Sheaf
from sage.schemes.chow.scheme import ChowScheme
from sage.rings.integer import Integer


def GrassBundle(A, B, chern_class='c', names=None, name=None, latex_name=None):
    r"""
    Return either
    - the Grassmannian of quotients of rank B of A if A is a sheaf or
    - the Grassmannian of subbundles of B of rank A if B is a bundle.
    """
    if isinstance(A, Sheaf):
        E, n = A, A.rank()
        r = int(B) if isinstance(B, Integer) else B
        if not isinstance(r, int):
            err = "Expect an integer as second argument if first is a Sheaf."
            raise TypeError(err)
        if (r <= 0) or (r >= n):
            err = "Expect 0 < %d < %d as second argument." % (r, n)
            raise ValueError(err)
        k = n - r
    elif isinstance(B, Bundle):
        E, n = B, B.rank()
        k = int(A) if isinstance(A, Integer) else A
        if not isinstance(A, int):
            err = "Expect an integer as first argument if second is a Sheaf."
            raise TypeError(err)
        if (k <= 0) or (k >= n):
            err = "Expect 0 < %d < %d as first argument." % (k, n)
            raise ValueError(err)
        r = n - k
    else:
        raise TypeError("Expect first or second argument to be a Sheaf.")
    if not isinstance(chern_class, str):
        raise ValueError("Need a string as Chern class.")

    # Now E is a sheaf of rank n and we are looking for quotients of rank r.
    # Keep the number of variables small.

    X = E.chowscheme()
    d = r * k + X.dimension()  # Dimension of the Grassmannian
    s = min(r, k)  # s new variables: c_1, ..., c_s

    # Find the new variables
    new_vars = [chern_class]
    if s > 1:
        new_vars = [chern_class + '%d' % i for i in range(1, s + 1)]
    while [c for c in new_vars if c in X.variable_names()]:
        new_vars = [chern_class + '%d' % i for i in range(1, s + 1)]
        chern_class += chern_class

    # Find the new degrees
    new_degs = [i for i in range(1, s + 1)]

    # Find the new relations

    Z = ChowScheme(d, generators=list(X.variable_names()) + new_vars,
                   degrees=list(X.degs()) + new_degs,
                   relations=[str(rel) for rel in X.rels()])
    #
    #  The universal exact sequence 0->S->p^{*}E->Q->0 gives the relations...
    #

    EZ = Bundle(Z, n, [str(c) for c in E.chern_classes()])  # p^{*}E
    QZ = Bundle(Z, r, [str(1)] + new_vars)  # if s==k : c_{s+1} =...= c_{r} = 0
    SZ = Bundle(Z, k, (EZ - QZ).chern_classes()[0: k + 1])
    QZ = Bundle(Z, r, (EZ - SZ).chern_classes()[0: r + 1]) if s == k else QZ

    new_rels = [str(rel) for rel in BundleDiffRelations(EZ, SZ)]

    # Construct GrE
    GrE = ChowScheme(d, generators=new_vars + list(X.variable_names()),
                     degrees=new_degs + list(X.degs()),
                     relations=new_rels + [str(rel) for rel in X.rels()],
                     names=names, name=name, latex_name=latex_name)

    if not X.is_point():
        f = GrE.hom(X.gens(), X)
        GrE = GrE.base_change(f)

    S = Bundle(GrE, k, SZ.chern_classes())
    Q = Bundle(GrE, r, QZ.chern_classes())

    point_class = GrE.chowring()(X.point_class())
    GrE.set_point_class(Q.chern_classes()[r] ** (n - r) * point_class)
    GrE.sheaves['universal_sub'] = S
    GrE.sheaves['universal_quotient'] = Q
    TGrE_X = SHom(S, Q)
    TGrE = TGrE_X + GrE.base_morphism().upperstar(X.tangent_bundle())
    GrE.sheaves['tangent'] = TGrE
    GrE.sheaves['o1'] = Q.determinant()
    return GrE


def Grass(n, r, chern_class='c', names=None, name=None, latex_name=None):
    """
    Return either depending respectively whether `n > r` or `n < r`:
    -  The Grassmannian of quotients of rank r of an n dimensional vector space;
    -  The Grassmannian of subspaces of rank n of an r dimensional vector space.
    EXAMPLES::
        sage: G = Grass(6, 4, chern_class='w')
        sage: G.dimension()
        8
        sage: G.rels()
        [w2^5, w1*w2^4, w1^2*w2^3 - 4/3*w2^4, w1^6 - 3*w1^2*w2^2 + w2^3, w1^3*w2 - 3/2*w1*w2^2]
        sage: G.sheaves["universal_sub"]
        Bundle(Grass(6, 4), 2, [1, -w1, w1^2 - w2])
        sage: G.sheaves["universal_quotient"]
        Bundle(Grass(6, 4), 4, [1, w1, w2, -w1^3 + 2*w1*w2, -w1^4 + w1^2*w2 + w2^2])
        sage: H = Grass(6, 2, chern_class='v')
        sage: H.dimension()
        8
        sage: H.rels()
        [v2^5, v1*v2^4, v1^2*v2^3 - v2^4, v1^3*v2^2 - 2*v1*v2^3, v1^4*v2 - 3*v1^2*v2^2 + v2^3, v1^5 - 4*v1^3*v2 + 3*v1*v2^2]
        sage: H.sheaves["universal_sub"]
        Bundle(Grass(6, 2), 4, [1, -v1, v1^2 - v2, -v1^3 + 2*v1*v2, v1^4 - 3*v1^2*v2 + v2^2])
        sage: H.sheaves["universal_quotient"]
        Bundle(Grass(6, 2), 2, [1, v1, v2])
    """
    n = int(n) if isinstance(n, Integer) else n
    r = int(r) if isinstance(r, Integer) else r
    if not (isinstance(n, int) and isinstance(r, int)):
        raise TypeError("Need integers, got %s and %s:" % (n, r))
    if r == n:
        return PointChowScheme
    if name is None:
        name = "Grass(%d, %d)" % (n, r)
    if latex_name is None:
        latex_name = r"Grass(%d, %d)" % (n, r)
    if r < n:
        # We are talking about quotients
        E = Bundle(PointChowScheme, n, [1])
        return GrassBundle(E, r, chern_class,
                           names=names, name=name, latex_name=latex_name)
    # We are talking about subspaces
    E = Bundle(PointChowScheme, r, [1])
    return GrassBundle(n, E, chern_class,
                       names=names, name=name, latex_name=latex_name)
