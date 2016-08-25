"""
Homomorphisms Between Matrix Groups

AUTHORS:

- David Joyner and William Stein (2006-03): initial version

- David Joyner (2006-05): examples

- Simon King (2011-01): cleaning and improving code

- Volker Braun (2013-1) port to new Parent, libGAP.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.categories.morphism import Morphism
from sage.misc.latex import latex


def to_libgap(x):
    """
    Helper to convert ``x`` to a LibGAP matrix or matrix group
    element.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.morphism import to_libgap
        sage: to_libgap(GL(2,3).gen(0))
        [ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        sage: to_libgap(matrix(QQ, [[1,2],[3,4]]))
        [ [ 1, 2 ], [ 3, 4 ] ]
    """
    try:
        return x.gap()
    except AttributeError:
        from sage.libs.gap.libgap import libgap
        return libgap(x)


class MatrixGroupMap(Morphism):

    def __init__(self, parent):
        """
        Set-theoretic map between matrix groups.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))   # mathematical nonsense
            MatrixGroup endomorphism of Integer Ring
        """
        Morphism.__init__(self, parent)

    def _repr_type(self):
        """
        Part of the implementation of :meth:`_repr_`

        EXAMPLES::

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))._repr_type()
            'MatrixGroup'
        """
        return "MatrixGroup"


class MatrixGroupMorphism(MatrixGroupMap):
    pass


class MatrixGroupMorphism_im_gens(MatrixGroupMorphism):

    def __init__(self, homset, imgsH, check=True):
        """
        Group morphism specified by images of generators.

        Some Python code for wrapping GAP's GroupHomomorphismByImages
        function but only for matrix groups. Can be expensive if G is
        large.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi
            Homomorphism : Matrix group over Finite Field of size 5 with 1 generators (
            [1 1]
            [0 1]
            ) --> Matrix group over Finite Field of size 5 with 1 generators (
            [1 0]
            [1 1]
            )
            sage: phi(MS([1,1,0,1]))
            [1 0]
            [1 1]
            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])

        TESTS:

        Check that :trac:`19406` is fixed::

            sage: G = GL(2, GF(3))
            sage: H = GL(3, GF(2))
            sage: mat1 = H([[-1,0,0],[0,0,-1],[0,-1,0]])
            sage: mat2 = H([[1,1,1],[0,0,-1],[-1,0,0]])
            sage: phi = G.hom([mat1, mat2])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a group homomorphism
        """
        MatrixGroupMorphism.__init__(self, homset)   # sets the parent
        from sage.libs.gap.libgap import libgap
        G = homset.domain()
        H = homset.codomain()
        gens = [x.gap() for x in G.gens()]
        imgs = [to_libgap(x) for x in imgsH]
        self._phi = phi = libgap.GroupHomomorphismByImages(G.gap(), H.gap(), gens, imgs)
        if not phi.IsGroupHomomorphism():
            raise ValueError('the map {}-->{} is not a homomorphism'.format(G.gens(), imgsH))

    def gap(self):
        """
        Return the underlying LibGAP group homomorphism

        OUTPUT:

        A LibGAP element.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi.gap()
            CompositionMapping( [ (6,7,8,10,9)(11,13,14,12,15)(16,19,20,18,17)(21,25,22,24,23) ]
            -> [ [ [ Z(5)^0, 0*Z(5) ], [ Z(5)^0, Z(5)^0 ] ] ], <action isomorphism> )
            sage: type(_)
            <type 'sage.libs.gap.element.GapElement'>
        """
        return self._phi

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: H = MatrixGroup([MS([1,0,1,1])])
            sage: phi = G.hom(H.gens())
            sage: phi
            Homomorphism : Matrix group over Finite Field of size 5 with 1 generators (
            [1 1]
            [0 1]
            ) --> Matrix group over Finite Field of size 5 with 1 generators (
            [1 0]
            [1 1]
            )
            sage: phi(MS([1,1,0,1]))
            [1 0]
            [1 1]
        """
        return "Homomorphism : %s --> %s"%(self.domain(),self.codomain())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS([1,1,0,1])])
            sage: phi = G.hom(G.gens())
            sage: print(latex(phi))
            \left\langle \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle \rightarrow{} \left\langle \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle
        """
        return "%s \\rightarrow{} %s"%(latex(self.domain()), latex(self.codomain()))

    def kernel(self):
        """
        Return the kernel of ``self``, i.e., a matrix group.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.kernel()
            Matrix group over Finite Field of size 7 with 1 generators (
            [6 0]
            [0 1]
            )
        """
        gap_ker = self.gap().Kernel()
        F = self.domain().base_ring()
        from sage.groups.matrix_gps.all import MatrixGroup
        return MatrixGroup([x.matrix(F) for x in gap_ker.GeneratorsOfGroup()])

    def pushforward(self, J, *args,**kwds):
        """
        The image of an element or a subgroup.

        INPUT:

        ``J`` -- a subgroup or an element of the domain of ``self``

        OUTPUT:

        The image of ``J`` under ``self``.

        .. NOTE::

            ``pushforward`` is the method that is used when a map is called
            on anything that is not an element of its domain. For historical
            reasons, we keep the alias ``image()`` for this method.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: F.multiplicative_generator()
            3
            sage: G = MatrixGroup([MS([3,0,0,1])])
            sage: a = G.gens()[0]^2
            sage: phi = G.hom([a])
            sage: phi.image(G.gens()[0]) # indirect doctest
            [2 0]
            [0 1]
            sage: H = MatrixGroup([MS(a.list())])
            sage: H
            Matrix group over Finite Field of size 7 with 1 generators (
            [2 0]
            [0 1]
            )

        The following tests against :trac:`10659`::

            sage: phi(H)   # indirect doctest
            Matrix group over Finite Field of size 7 with 1 generators (
            [4 0]
            [0 1]
            )
        """
        phi = self.gap()
        F = self.codomain().base_ring()
        gapJ = to_libgap(J)
        if gapJ.IsGroup():
            from sage.groups.matrix_gps.all import MatrixGroup
            img_gens = [x.matrix(F) for x in phi.Image(gapJ).GeneratorsOfGroup()]
            return MatrixGroup(img_gens)
        C = self.codomain()
        return C(phi.Image(gapJ).matrix(F))

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function for a
        matrix group ``G``. Returns an error if ``g`` is not in ``G``.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: g = MS([1,1,0,1])
            sage: G = MatrixGroup([g])
            sage: phi = G.hom(G.gens())
            sage: phi(G.0)
            [1 1]
            [0 1]
            sage: phi(G(g^2))
            [1 2]
            [0 1]

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([1,2,  -1,1]),MS([1,1,  0,1])]
            sage: G = MatrixGroup(gens)
            sage: phi = G.hom(G.gens())
            sage: phi(G.0)
            [1 2]
            [4 1]
            sage: phi(G.1)
            [1 1]
            [0 1]

        TESTS:

        The following tests that the call method was successfully
        improved in :trac:`10659`::

            sage: O = WeylGroup(['D',6])
            sage: r = prod(O.gens())
            sage: r_ = r^-1
            sage: f = O.hom([r*x*r_ for x in O.gens()])  # long time (19s on sage.math, 2011)
            sage: [f(x) for x in O.gens()]  # long time
            [
            [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
            [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
            <BLANKLINE>
            [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
            [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
            [1 0 0 0 0 0], [-1  0  0  0  0  0]
            ]
            sage: f(O)  # long time
            Matrix group over Rational Field with 6 generators
            sage: f(O).gens()   # long time
            (
            [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
            [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
            <BLANKLINE>
            [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
            [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
            [1 0 0 0 0 0], [-1  0  0  0  0  0]
            )

        We check that :trac:`19780` is fixed::

            sage: G = groups.matrix.SO(3, 3)
            sage: H = groups.matrix.GL(3, 3)
            sage: phi = G.hom([H(x) for x in G.gens()])
            sage: phi(G.one()).parent()
            General Linear Group of degree 3 over Finite Field of size 3
        """
        phi = self.gap()
        G = self.domain()
        C = self.codomain()
        F = C.base_ring()
        h = g.gap()
        return C(phi.Image(h).matrix(F))

