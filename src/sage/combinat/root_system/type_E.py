"""
Root system data for type E
"""
# ****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import ambient_space
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family


class AmbientSpace(ambient_space.AmbientSpace):
    """
    The lattice behind E6, E7, or E8.  The computations are based on Bourbaki,
    Groupes et Algèbres de Lie, Ch. 4,5,6 (planche V-VII).
    """

    def __init__(self, root_system, baseRing):
        """
        Create the ambient space for the root system for E6, E7, E8.
        Specify the Base, i.e., the simple roots w.r. to the canonical
        basis for R^8.

        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: e == loads(dumps(e))
            True
            sage: [e.weyl_dimension(v) for v in e.fundamental_weights()]
            [27, 78, 351, 2925, 351, 27]
            sage: e = RootSystem(['E',7]).ambient_space()
            sage: [e.weyl_dimension(v) for v in e.fundamental_weights()]
            [133, 912, 8645, 365750, 27664, 1539, 56]
            sage: e = RootSystem(['E',8]).ambient_space()
            sage: [e.weyl_dimension(v) for v in e.fundamental_weights()]
            [3875, 147250, 6696000, 6899079264, 146325270, 2450240, 30380, 248]
        """
        v = ZZ(1)/ZZ(2)
        self.rank = root_system.cartan_type().rank()
        ambient_space.AmbientSpace.__init__(self, root_system, baseRing)
        if self.rank == 6:
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1)]
        elif self.rank == 7:
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1),
                         self.root(4,5,p1=1)]
        elif self.rank == 8:
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1),
                         self.root(4,5,p1=1),
                         self.root(5,6,p1=1)]
        else:
            raise NotImplementedError("Type \'E\' root systems only come in flavors 6, 7, 8.  Please make another choice")

    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: e.dimension()
            8
        """
        return 8

    def root(self, i1, i2=None, i3=None, i4=None, i5=None, i6=None, i7=None, i8=None, p1=0, p2=0, p3=0, p4=0, p5=0, p6=0, p7=0, p8=0):
        """
        Compute an element of the underlying lattice, using the specified elements of
        the standard basis, with signs dictated by the corresponding 'pi' arguments.
        We rely on the caller to provide the correct arguments.
        This is typically used to generate roots, although the generated elements
        need not be roots themselves.
        We assume that if one of the indices is not given, the rest are not as well.
        This should work for E6, E7, E8.

        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: [ e.root(i, j, p3=1) for i in range(e.n) for j in range(i+1, e.n) ]
            [(1, 1, 0, 0, 0, 0, 0, 0),
             (1, 0, 1, 0, 0, 0, 0, 0),
             (1, 0, 0, 1, 0, 0, 0, 0),
             (1, 0, 0, 0, 1, 0, 0, 0),
             (1, 0, 0, 0, 0, 1, 0, 0),
             (1, 0, 0, 0, 0, 0, 1, 0),
             (1, 0, 0, 0, 0, 0, 0, 1),
             (0, 1, 1, 0, 0, 0, 0, 0),
             (0, 1, 0, 1, 0, 0, 0, 0),
             (0, 1, 0, 0, 1, 0, 0, 0),
             (0, 1, 0, 0, 0, 1, 0, 0),
             (0, 1, 0, 0, 0, 0, 1, 0),
             (0, 1, 0, 0, 0, 0, 0, 1),
             (0, 0, 1, 1, 0, 0, 0, 0),
             (0, 0, 1, 0, 1, 0, 0, 0),
             (0, 0, 1, 0, 0, 1, 0, 0),
             (0, 0, 1, 0, 0, 0, 1, 0),
             (0, 0, 1, 0, 0, 0, 0, 1),
             (0, 0, 0, 1, 1, 0, 0, 0),
             (0, 0, 0, 1, 0, 1, 0, 0),
             (0, 0, 0, 1, 0, 0, 1, 0),
             (0, 0, 0, 1, 0, 0, 0, 1),
             (0, 0, 0, 0, 1, 1, 0, 0),
             (0, 0, 0, 0, 1, 0, 1, 0),
             (0, 0, 0, 0, 1, 0, 0, 1),
             (0, 0, 0, 0, 0, 1, 1, 0),
             (0, 0, 0, 0, 0, 1, 0, 1),
             (0, 0, 0, 0, 0, 0, 1, 1)]
        """
        if i1 == i2 or i2 is None:
            return (-1)**p1*self.monomial(i1)
        if i3 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)
        if i4 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)
        if i5 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)+(-1)**p4*self.monomial(i4)
        if i6 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)+(-1)**p4*self.monomial(i4)+(-1)**p5*self.monomial(i5)
        if i7 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)+(-1)**p4*self.monomial(i4)+(-1)**p5*self.monomial(i5)+(-1)**p6*self.monomial(i6)
        if i8 is None:
            return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)+(-1)**p4*self.monomial(i4)+(-1)**p5*self.monomial(i5)+(-1)**p6*self.monomial(i6)+(-1)**p7*self.monomial(i7)
        return (-1)**p1*self.monomial(i1) + (-1)**p2*self.monomial(i2)+(-1)**p3*self.monomial(i3)+(-1)**p4*self.monomial(i4)+(-1)**p5*self.monomial(i5)+(-1)**p6*self.monomial(i6)+(-1)**p7*self.monomial(i7)+(-1)**p8*self.monomial(i8)

    def simple_root(self, i):
        """
        There are computed as what Bourbaki calls the Base:
            a1 = e2-e3, a2 = e3-e4, a3 = e4, a4 = 1/2*(e1-e2-e3-e4)

        EXAMPLES::

            sage: LE6 = RootSystem(['E',6]).ambient_space()
            sage: LE6.simple_roots()
            Finite family {1: (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2), 2: (1, 1, 0, 0, 0, 0, 0, 0), 3: (-1, 1, 0, 0, 0, 0, 0, 0), 4: (0, -1, 1, 0, 0, 0, 0, 0), 5: (0, 0, -1, 1, 0, 0, 0, 0), 6: (0, 0, 0, -1, 1, 0, 0, 0)}
        """
        if i not in self.index_set():
            raise ValueError("{} is not in the index set".format(i))
        return self.Base[i-1]

    def negative_roots(self):
        """
        The negative roots.

        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: e.negative_roots()
            [(-1, -1, 0, 0, 0, 0, 0, 0),
             (-1, 0, -1, 0, 0, 0, 0, 0),
             (-1, 0, 0, -1, 0, 0, 0, 0),
             (-1, 0, 0, 0, -1, 0, 0, 0),
             (0, -1, -1, 0, 0, 0, 0, 0),
             (0, -1, 0, -1, 0, 0, 0, 0),
             (0, -1, 0, 0, -1, 0, 0, 0),
             (0, 0, -1, -1, 0, 0, 0, 0),
             (0, 0, -1, 0, -1, 0, 0, 0),
             (0, 0, 0, -1, -1, 0, 0, 0),
             (1, -1, 0, 0, 0, 0, 0, 0),
             (1, 0, -1, 0, 0, 0, 0, 0),
             (1, 0, 0, -1, 0, 0, 0, 0),
             (1, 0, 0, 0, -1, 0, 0, 0),
             (0, 1, -1, 0, 0, 0, 0, 0),
             (0, 1, 0, -1, 0, 0, 0, 0),
             (0, 1, 0, 0, -1, 0, 0, 0),
             (0, 0, 1, -1, 0, 0, 0, 0),
             (0, 0, 1, 0, -1, 0, 0, 0),
             (0, 0, 0, 1, -1, 0, 0, 0),
             (-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2)]
        """
        return [ -a for a in self.positive_roots()]

    def positive_roots(self):
        """
        These are the roots positive w.r. to lexicographic ordering of the
        basis elements (e1<...<e4).

        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: e.positive_roots()
            [(1, 1, 0, 0, 0, 0, 0, 0),
             (1, 0, 1, 0, 0, 0, 0, 0),
             (1, 0, 0, 1, 0, 0, 0, 0),
             (1, 0, 0, 0, 1, 0, 0, 0),
             (0, 1, 1, 0, 0, 0, 0, 0),
             (0, 1, 0, 1, 0, 0, 0, 0),
             (0, 1, 0, 0, 1, 0, 0, 0),
             (0, 0, 1, 1, 0, 0, 0, 0),
             (0, 0, 1, 0, 1, 0, 0, 0),
             (0, 0, 0, 1, 1, 0, 0, 0),
             (-1, 1, 0, 0, 0, 0, 0, 0),
             (-1, 0, 1, 0, 0, 0, 0, 0),
             (-1, 0, 0, 1, 0, 0, 0, 0),
             (-1, 0, 0, 0, 1, 0, 0, 0),
             (0, -1, 1, 0, 0, 0, 0, 0),
             (0, -1, 0, 1, 0, 0, 0, 0),
             (0, -1, 0, 0, 1, 0, 0, 0),
             (0, 0, -1, 1, 0, 0, 0, 0),
             (0, 0, -1, 0, 1, 0, 0, 0),
             (0, 0, 0, -1, 1, 0, 0, 0),
             (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2),
             (1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2),
             (1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2),
             (1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2),
             (1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2),
             (1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2),
             (1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2),
             (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2),
             (-1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2),
             (-1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2),
             (-1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2),
             (-1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2),
             (-1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2),
             (-1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2),
             (-1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2),
             (-1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2)]
            sage: e.rho()
            (0, 1, 2, 3, 4, -4, -4, 4)
            sage: E8 = RootSystem(['E',8])
            sage: e = E8.ambient_space()
            sage: e.negative_roots()
            [(-1, -1, 0, 0, 0, 0, 0, 0),
             (-1, 0, -1, 0, 0, 0, 0, 0),
             (-1, 0, 0, -1, 0, 0, 0, 0),
             (-1, 0, 0, 0, -1, 0, 0, 0),
             (-1, 0, 0, 0, 0, -1, 0, 0),
             (-1, 0, 0, 0, 0, 0, -1, 0),
             (-1, 0, 0, 0, 0, 0, 0, -1),
             (0, -1, -1, 0, 0, 0, 0, 0),
             (0, -1, 0, -1, 0, 0, 0, 0),
             (0, -1, 0, 0, -1, 0, 0, 0),
             (0, -1, 0, 0, 0, -1, 0, 0),
             (0, -1, 0, 0, 0, 0, -1, 0),
             (0, -1, 0, 0, 0, 0, 0, -1),
             (0, 0, -1, -1, 0, 0, 0, 0),
             (0, 0, -1, 0, -1, 0, 0, 0),
             (0, 0, -1, 0, 0, -1, 0, 0),
             (0, 0, -1, 0, 0, 0, -1, 0),
             (0, 0, -1, 0, 0, 0, 0, -1),
             (0, 0, 0, -1, -1, 0, 0, 0),
             (0, 0, 0, -1, 0, -1, 0, 0),
             (0, 0, 0, -1, 0, 0, -1, 0),
             (0, 0, 0, -1, 0, 0, 0, -1),
             (0, 0, 0, 0, -1, -1, 0, 0),
             (0, 0, 0, 0, -1, 0, -1, 0),
             (0, 0, 0, 0, -1, 0, 0, -1),
             (0, 0, 0, 0, 0, -1, -1, 0),
             (0, 0, 0, 0, 0, -1, 0, -1),
             (0, 0, 0, 0, 0, 0, -1, -1),
             (1, -1, 0, 0, 0, 0, 0, 0),
             (1, 0, -1, 0, 0, 0, 0, 0),
             (1, 0, 0, -1, 0, 0, 0, 0),
             (1, 0, 0, 0, -1, 0, 0, 0),
             (1, 0, 0, 0, 0, -1, 0, 0),
             (1, 0, 0, 0, 0, 0, -1, 0),
             (1, 0, 0, 0, 0, 0, 0, -1),
             (0, 1, -1, 0, 0, 0, 0, 0),
             (0, 1, 0, -1, 0, 0, 0, 0),
             (0, 1, 0, 0, -1, 0, 0, 0),
             (0, 1, 0, 0, 0, -1, 0, 0),
             (0, 1, 0, 0, 0, 0, -1, 0),
             (0, 1, 0, 0, 0, 0, 0, -1),
             (0, 0, 1, -1, 0, 0, 0, 0),
             (0, 0, 1, 0, -1, 0, 0, 0),
             (0, 0, 1, 0, 0, -1, 0, 0),
             (0, 0, 1, 0, 0, 0, -1, 0),
             (0, 0, 1, 0, 0, 0, 0, -1),
             (0, 0, 0, 1, -1, 0, 0, 0),
             (0, 0, 0, 1, 0, -1, 0, 0),
             (0, 0, 0, 1, 0, 0, -1, 0),
             (0, 0, 0, 1, 0, 0, 0, -1),
             (0, 0, 0, 0, 1, -1, 0, 0),
             (0, 0, 0, 0, 1, 0, -1, 0),
             (0, 0, 0, 0, 1, 0, 0, -1),
             (0, 0, 0, 0, 0, 1, -1, 0),
             (0, 0, 0, 0, 0, 1, 0, -1),
             (0, 0, 0, 0, 0, 0, 1, -1),
             (-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2),
             (-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2),
             (-1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2),
             (-1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2),
             (-1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2),
             (-1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2),
             (-1/2, -1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2),
             (-1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2),
             (-1/2, -1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2),
             (-1/2, -1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2),
             (-1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2),
             (-1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2),
             (-1/2, 1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2),
             (-1/2, 1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2),
             (-1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2),
             (-1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2),
             (-1/2, 1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (-1/2, 1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2),
             (-1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2),
             (-1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2),
             (-1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2),
             (-1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2),
             (-1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2),
             (1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2),
             (1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2),
             (1/2, -1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2),
             (1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2),
             (1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2),
             (1/2, -1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2),
             (1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2),
             (1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2),
             (1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2),
             (1/2, 1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2),
             (1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2),
             (1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2),
             (1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2),
             (1/2, 1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2),
             (1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2),
             (1/2, 1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2),
             (1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2),
             (1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2)]
            sage: e.rho()
            (0, 1, 2, 3, 4, 5, 6, 23)
        """
        v = ZZ(1)/ZZ(2)
        # Note that
        if not hasattr(self, 'PosRoots'):
            if self.rank == 6:
                self.PosRoots = ( [ self.root(i,j) for i in range(self.rank-1) for j in range(i+1,self.rank-1) ] +
                                  [ self.root(i,j,p1=1) for i in range(self.rank-1) for j in range(i+1,self.rank-1) ] +
                                  [ v*(self.root(7)-self.root(6)-self.root(5)+self.root(0,1,2,3,4,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] if (p1+p2+p3+p4+p5) % 2 == 0 ])
            elif self.rank == 7:
                self.PosRoots = ( [ self.root(i,j) for i in range(self.rank-1) for j in range(i+1,self.rank-1) ] +
                                  [ self.root(i,j,p1=1) for i in range(self.rank-1) for j in range(i+1,self.rank-1) ] +
                                  [ self.root(6,7,p1=1) ] +
                                  [ v*(self.root(7)-self.root(6)+self.root(0,1,2,3,4,5,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] for p6 in [0,1] if (p1+p2+p3+p4+p5+p6) % 2 == 1 ])
            elif self.rank == 8:
                self.PosRoots = ( [ self.root(i,j) for i in range(self.rank) for j in range(i+1,self.rank) ] +
                                  [ self.root(i,j,p1=1) for i in range(self.rank) for j in range(i+1,self.rank) ] +
                                  [ v*(self.root(7)+self.root(0,1,2,3,4,5,6,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] for p6 in [0,1] for p7 in [0,1] if (p1+p2+p3+p4+p5+p6+p7) % 2 == 0 ])

        return self.PosRoots

    def fundamental_weights(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['E',6]).ambient_space()
            sage: e.fundamental_weights()
            Finite family {1: (0, 0, 0, 0, 0, -2/3, -2/3, 2/3), 2: (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2), 3: (-1/2, 1/2, 1/2, 1/2, 1/2, -5/6, -5/6, 5/6), 4: (0, 0, 1, 1, 1, -1, -1, 1), 5: (0, 0, 0, 1, 1, -2/3, -2/3, 2/3), 6: (0, 0, 0, 0, 1, -1/3, -1/3, 1/3)}
        """
        v2 = ZZ(1)/ZZ(2)
        v3 = ZZ(1)/ZZ(3)
        if self.rank == 6:
            return Family({ 1: 2*v3*self.root(7,6,5,p2=1,p3=1),
                            2: v2*self.root(0,1,2,3,4,5,6,7,p6=1,p7=1),
                            3: 5*v2*v3*self.root(7,6,5,p2=1,p3=1)+v2*self.root(0,1,2,3,4,p1=1),
                            4: self.root(2,3,4,5,6,7,p4=1,p5=1),
                            5: 2*v3*self.root(7,6,5,p2=1,p3=1)+self.root(3,4),
                            6: v3*self.root(7,6,5,p2=1,p3=1)+self.root(4)})
        elif self.rank == 7:
            return Family({ 1: self.root(7,6,p2=1),
                            2: v2*self.root(0,1,2,3,4,5)+self.root(6,7,p1=1),
                            3: v2*(self.root(0,1,2,3,4,5,p1=1)+3*self.root(6,7,p1=1)),
                            4: self.root(2,3,4,5)+2*self.root(6,7,p1=1),
                            5: 3*v2*self.root(6,7,p1=1)+self.root(3,4,5),
                            6: self.root(4,5,6,7,p3=1),
                            7: self.root(5)+v2*self.root(6,7,p1=1)})
        elif self.rank == 8:
            return Family({ 1: 2*self.root(7),
                            2: v2*(self.root(0,1,2,3,4,5,6)+5*self.root(7)),
                            3: v2*(self.root(0,1,2,3,4,5,6,p1=1)+7*self.root(7)),
                            4: self.root(2,3,4,5,6)+5*self.root(7),
                            5: self.root(3,4,5,6)+4*self.root(7),
                            6: self.root(4,5,6)+3*self.root(7),
                            7: self.root(5,6)+2*self.root(7),
                            8: self.root(6,7)})


from .cartan_type import CartanType_standard_finite, CartanType_simple, CartanType_simply_laced


class CartanType(CartanType_standard_finite, CartanType_simple, CartanType_simply_laced):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['E',6])
            sage: ct
            ['E', 6]
            sage: ct._repr_(compact = True)
            'E6'
            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['E', 6, 1]
            sage: ct.dual()
            ['E', 6]

        TESTS::

            sage: TestSuite(ct).run()
        """
        if n < 6 or n > 8:
            raise ValueError("Invalid Cartan Type for Type E")
        CartanType_standard_finite.__init__(self, "E", n)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['E',7]))
            E_7
        """
        return "E_%s" % self.n

    AmbientSpace = AmbientSpace

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['E',6]).coxeter_number()
            12
            sage: CartanType(['E',7]).coxeter_number()
            18
            sage: CartanType(['E',8]).coxeter_number()
            30
        """
        if self.n == 6:
            return 12
        if self.n == 7:
            return 18
        # n == 8
        return 30

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['E',6]).dual_coxeter_number()
            12
            sage: CartanType(['E',7]).dual_coxeter_number()
            18
            sage: CartanType(['E',8]).dual_coxeter_number()
            30
        """
        if self.n == 6:
            return 12
        if self.n == 7:
            return 18
        # n == 8
        return 30

    def dynkin_diagram(self):
        """
        Return a Dynkin diagram for type E.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: e = CartanType(['E',6]).dynkin_diagram(); e
                    O 2
                    |
                    |
            O---O---O---O---O
            1   3   4   5   6
            E6
            sage: e.edges(sort=True)
            [(1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1), (4, 2, 1),
             (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1), (6, 5, 1)]
            sage: e = CartanType(['E',7]).dynkin_diagram(); e
                    O 2
                    |
                    |
            O---O---O---O---O---O
            1   3   4   5   6   7
            E7
            sage: e.edges(sort=True)
            [(1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1), (4, 2, 1),
             (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1), (6, 5, 1),
             (6, 7, 1), (7, 6, 1)]
            sage: e = CartanType(['E',8]).dynkin_diagram(); e
                    O 2
                    |
                    |
            O---O---O---O---O---O---O
            1   3   4   5   6   7   8
            E8
            sage: e.edges(sort=True)
            [(1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1), (4, 2, 1),
             (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1), (6, 5, 1),
             (6, 7, 1), (7, 6, 1), (7, 8, 1), (8, 7, 1)]
        """
        from .dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        g.add_edge(1,3)
        g.add_edge(2,4)
        for i in range(3, self.n):
            g.add_edge(i, i+1)
        return g

    def _latex_dynkin_diagram(self, label=None, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['E',7])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (10 cm,0);
            \draw (4 cm, 0 cm) -- +(0,2 cm);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$5$};
            \draw[fill=white] (8 cm, 0 cm) circle (.25cm) node[below=4pt]{$6$};
            \draw[fill=white] (10 cm, 0 cm) circle (.25cm) node[below=4pt]{$7$};
            \draw[fill=white] (4 cm, 2 cm) circle (.25cm) node[right=3pt]{$2$};
            <BLANKLINE>
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._latex_draw_node
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n" % ((self.n-2)*node_dist)
        ret += "\\draw (%s cm, 0 cm) -- +(0,%s cm);\n" % (2*node_dist, node_dist)
        ret += node(0, 0, label(1))
        for i in range(1, self.n-1):
            ret += node(i*node_dist, 0, label(i+2))
        ret += node(2*node_dist, node_dist, label(2), 'right=3pt')
        return ret

    def ascii_art(self, label=None, node=None):
        """
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['E',6]).ascii_art(label = lambda x: x+2))
                    O 4
                    |
                    |
            O---O---O---O---O
            3   5   6   7   8
            sage: print(CartanType(['E',7]).ascii_art(label = lambda x: x+2))
                    O 4
                    |
                    |
            O---O---O---O---O---O
            3   5   6   7   8   9
            sage: print(CartanType(['E',8]).ascii_art(label = lambda x: x+1))
                    O 3
                    |
                    |
            O---O---O---O---O---O---O
            2   4   5   6   7   8   9
        """
        if label is None:
            label = lambda i: i
        if node is None:
            node = self._ascii_art_node
        labels = [label(i) for i in [1,3,4,5,6] + list(range(7, self.n+1))] # We exclude 2 because of the special case
        ret = "        {} {}\n        |\n        |\n".format(node(label(2)), label(2))
        return ret + '---'.join(node(i) for i in labels) + '\n' + "".join("{!s:4}".format(i) for i in labels)


# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_E',
                           'ambient_space',
                           AmbientSpace)
