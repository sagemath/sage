"""
Homomorphisms of abelian groups

TODO:

- there must be a homspace first

- there should be hom and Hom methods in abelian group

AUTHORS:

- David Joyner (2006-03-03): initial version
"""

#*****************************************************************************
#   Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#
#                    http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.groups.perm_gps.permgroup import *
from sage.interfaces.gap import *
from sage.categories.morphism import *
from sage.categories.homset import *

from sage.misc.all import prod

def is_AbelianGroupMorphism(f):
    return isinstance(f, AbelianGroupMorphism);

class AbelianGroupMap(Morphism):
    """
    A set-theoretic map between AbelianGroups.
    """
    def __init__(self, parent):
        """
        The Python constructor.
        """
        Morphism.__init__(self, parent)

    def _repr_type(self):
        return "AbelianGroup"

class AbelianGroupMorphism_id(AbelianGroupMap):
    """
    Return the identity homomorphism from X to itself.

    EXAMPLES:
    """
    def __init__(self, X):
        AbelianGroupMorphism.__init__(self, X.Hom(X))

    def _repr_defn(self):
        return "Identity map of "+str(X)

class AbelianGroupMorphism(Morphism):
    """
    Some python code for wrapping GAP's GroupHomomorphismByImages
    function for abelian groups. Returns "fail" if gens does not
    generate self or if the map does not extend to a group
    homomorphism, self - other.

    EXAMPLES::

        sage: G = AbelianGroup(3,[2,3,4],names="abc"); G
        Multiplicative Abelian group isomorphic to C2 x C3 x C4
        sage: a,b,c = G.gens()
        sage: H = AbelianGroup(2,[2,3],names="xy"); H
        Multiplicative Abelian group isomorphic to C2 x C3
        sage: x,y = H.gens()

    ::

        sage: from sage.groups.abelian_gps.abelian_group_morphism import AbelianGroupMorphism
        sage: phi = AbelianGroupMorphism(H,G,[x,y],[a,b])

    AUTHORS:

    - David Joyner (2006-02)
    """

#    There is a homomorphism from H to G but not from G to H:
#
#        sage: phi = AbelianGroupMorphism_im_gens(G,H,[a*b,a*c],[x,y])
#------------------------------------------------------------
#Traceback (most recent call last):
#  File "<ipython console>", line 1, in ?
#  File ".abeliangp_hom.sage.py", line 737, in __init__
#    raise TypeError, "Sorry, the orders of the corresponding elements in %s, %s must be equal."%(genss,imgss)
#TypeError: Sorry, the orders of the corresponding elements in [a*b, a*c], [x, y] must be equal.
#
#        sage: phi = AbelianGroupMorphism_im_gens(G,H,[a*b,(a*c)^2],[x*y,y])
#------------------------------------------------------------
#Traceback (most recent call last):
#  File "<ipython console>", line 1, in ?
#  File ".abeliangp_hom.sage.py", line 730, in __init__
#    raise TypeError, "Sorry, the list %s must generate G."%genss
#TypeError: Sorry, the list [a*b, c^2] must generate G.

    def __init__(self, G, H, genss, imgss ):
        from sage.categories.homset import Hom
        Morphism.__init__(self, Hom(G, H))
        if len(genss) != len(imgss):
            raise TypeError("Sorry, the lengths of %s, %s must be equal."%(genss,imgss))
        self._domain = G
        self._codomain = H
        if not(G.is_abelian()):
            raise TypeError("Sorry, the groups must be abelian groups.")
        if not(H.is_abelian()):
            raise TypeError("Sorry, the groups must be abelian groups.")
        G_domain = G.subgroup(genss)
        if G_domain.order() != G.order():
            raise TypeError("Sorry, the list %s must generate G."%genss)
        # self.domain_invs = G.gens_orders()
        # self.codomaininvs = H.gens_orders()
        self.domaingens = genss
        self.codomaingens = imgss
        for i in range(len(self.domaingens)):
            if (self.domaingens[i]).order() != (self.codomaingens[i]).order():
                raise TypeError("Sorry, the orders of the corresponding elements in %s, %s must be equal."%(genss,imgss))

    def _gap_init_(self):
        """
        Only works for finite groups.

        EXAMPLES::

            sage: G = AbelianGroup(3,[2,3,4],names="abc"); G
            Multiplicative Abelian group isomorphic to C2 x C3 x C4
            sage: a,b,c = G.gens()
            sage: H = AbelianGroup(2,[2,3],names="xy"); H
            Multiplicative Abelian group isomorphic to C2 x C3
            sage: x,y = H.gens()
            sage: phi = AbelianGroupMorphism(H,G,[x,y],[a,b])
            sage: phi._gap_init_()
            'phi := GroupHomomorphismByImages(G,H,[x, y],[a, b])'
        """
        G  = (self.domain())._gap_init_()
        H  = (self.codomain())._gap_init_()
        s3 = 'G:=%s; H:=%s'%(G,H)
        gap.eval(s3)
        gensG = self.domain().variable_names()                    ## the Sage group generators
        gensH = self.codomain().variable_names()
        s1 = "gensG := GeneratorsOfGroup(G)"          ## the GAP group generators
        gap.eval(s1)
        s2 = "gensH := GeneratorsOfGroup(H)"
        gap.eval(s2)
        for i in range(len(gensG)):                      ## making the Sage group gens
           cmd = ("%s := gensG["+str(i+1)+"]")%gensG[i]  ## correspond to the Sage group gens
           gap.eval(cmd)
        for i in range(len(gensH)):
           cmd = ("%s := gensH["+str(i+1)+"]")%gensH[i]
           gap.eval(cmd)
        args = str(self.domaingens)+","+ str(self.codomaingens)
        gap.eval("phi := GroupHomomorphismByImages(G,H,"+args+")")
        self.gap_hom_string = "phi := GroupHomomorphismByImages(G,H,"+args+")"
        return self.gap_hom_string

    def _repr_type(self):
        return "AbelianGroup"

    def kernel(self):
        """
        Only works for finite groups.

        TODO: not done yet; returns a gap object but should return a Sage
        group.

        EXAMPLES::

            sage: H = AbelianGroup(3,[2,3,4],names="abc"); H
            Multiplicative Abelian group isomorphic to C2 x C3 x C4
            sage: a,b,c = H.gens()
            sage: G = AbelianGroup(2,[2,3],names="xy"); G
            Multiplicative Abelian group isomorphic to C2 x C3
            sage: x,y = G.gens()
            sage: phi = AbelianGroupMorphism(G,H,[x,y],[a,b])
            sage: phi.kernel()
            'Group([  ])'
        """
        cmd = self._gap_init_()
        gap.eval(cmd)
        return gap.eval("Kernel(phi)")

    def image(self, J):
        """
        Only works for finite groups.

        J must be a subgroup of G. Computes the subgroup of H which is the
        image of J.

        EXAMPLES::

            sage: G = AbelianGroup(2,[2,3],names="xy")
            sage: x,y = G.gens()
            sage: H = AbelianGroup(3,[2,3,4],names="abc")
            sage: a,b,c = H.gens()
            sage: phi = AbelianGroupMorphism(G,H,[x,y],[a,b])
        """
        G = self.domain()
        gensJ = J.gens()
        for g in gensJ:
            print(g)
            print(self(g), "\n")
        L = [self(g) for g in gensJ]
        return G.subgroup(L)

    def _call_( self, g ):
        """
        Some python code for wrapping GAP's Images function but only for
        permutation groups. Returns an error if g is not in G.

        EXAMPLES::

            sage: H = AbelianGroup(3, [2,3,4], names="abc")
            sage: a,b,c = H.gens()
            sage: G = AbelianGroup(2, [2,3], names="xy")
            sage: x,y = G.gens()
            sage: phi = AbelianGroupMorphism(G,H,[x,y],[a,b])
            sage: phi(y*x)
            a*b
            sage: phi(y^2)
            b^2
        """
        G = g.parent()
        w = g.word_problem(self.domaingens)
        n = len(w)
        # g.word_problem is faster in general than word_problem(g)
        gens = self.codomaingens
        h = prod([gens[(self.domaingens).index(w[i][0])]**(w[i][1]) for i in range(n)])
        return h
