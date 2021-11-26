r"""
Links

A knot is defined as embedding of the circle `\mathbb{S}^1` in the
3-dimensional sphere `\mathbb{S}^3`, considered up to ambient isotopy.
They represent the physical idea of a knotted rope, but with the
particularity that the rope is closed. That is, the ends of the
rope are joined.

A link is an embedding of one or more copies of `\mathbb{S}^1` in
`\mathbb{S}^3`, considered up to ambient isotopy. That is, a link
represents the idea of one or more tied ropes. Every knot is a link,
but not every link is a knot.

Generically, the projection of a link on `\RR^2` is a curve with
crossings. The crossings are represented to show which strand goes
over the other. This curve is called a planar diagram of the link.
If we remove the crossings, the resulting connected components are
segments. These segments are called the edges of the diagram.

REFERENCES:

- :wikipedia:`Knot_(mathematics)`
- [Col2013]_
- [KnotAtlas]_

.. SEEALSO::

    There are also tables of link and knot invariants at web-pages
    `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
    `LinkInfo <https://linkinfo.sitehost.iu.edu>`__. These can be
    used inside Sage after installing the optional package
    ``database_knotinfo`` (type ``sage -i database_knotinfo`` in a command shell,
    see :mod:`~sage.knots.knotinfo`).

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
- Sebastian Oehms (October 2020, add :meth:`get_knotinfo` and meth:`is_isotopic`)
"""

# ****************************************************************************
#       Copyright (C) 2014  Miguel Angel Marco Buzunariz
#                           Amit Jamadagni
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.symbolic.ring import SR
from sage.rings.integer import Integer
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.functions.generalized import sign
from sage.homology.chain_complex import ChainComplex
from sage.misc.flatten import flatten
from sage.misc.cachefunc import cached_method
from copy import deepcopy, copy
from itertools import combinations
from sage.structure.sage_object import SageObject


class Link(SageObject):
    r"""
    A link.

    A link is an embedding of one or more copies of `\mathbb{S}^1` in
    `\mathbb{S}^3`, considered up to ambient isotopy. That is, a link
    represents the idea of one or more tied ropes. Every knot is a link,
    but not every link is a knot.

    A link can be created by using one of the conventions mentioned below:

    Braid:

    - The closure of a braid is a link::

        sage: B = BraidGroup(8)
        sage: L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 3]))
        sage: L
        Link with 1 component represented by 9 crossings
        sage: L = Link(B([1, 2, 1, -2, -1]))
        sage: L
        Link with 2 components represented by 5 crossings

      .. NOTE::

          The strands of the braid that have no crossings at all
          are removed.

    - Oriented Gauss Code:

      Label the crossings from `1` to `n` (where `n` is the number of
      crossings) and start moving along the link. Trace every component of
      the link, by starting at a particular point on one component of the
      link and writing down each of the crossings that you encounter until
      returning to the starting point. The crossings are written with sign
      depending on whether we cross them as over or undercrossing. Each
      component is then represented as a list whose elements are the
      crossing numbers. A second list of `+1` and `-1`'s keeps track of
      the orientation of each crossing::

        sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
        ....:           [-1, -1, -1, -1, 1, 1, -1, 1]])
        sage: L
        Link with 1 component represented by 8 crossings

      For links there may be more than one component and the input is
      as follows::

        sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
        sage: L
        Link with 3 components represented by 4 crossings

    - Planar Diagram (PD) Code:

      The diagram of the link is formed by segments that are adjacent to
      the crossings. Label each one of this segments with a positive number,
      and for each crossing, write down the four incident segments. The
      order of these segments is clockwise, starting with the incoming
      undercrossing.

      There is no particular distinction between knots and links for
      this input.

    EXAMPLES:

    One of the representations of the trefoil knot::

        sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sage: L
        Link with 1 component represented by 3 crossings

    .. PLOT::
        :width: 300 px

        L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sphinx_plot(L.plot())

    One of the representations of the Hopf link::

        sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sage: L
        Link with 2 components represented by 2 crossings

    .. PLOT::
        :width: 300 px

        L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sphinx_plot(L.plot())

    We can construct links from the braid group::

        sage: B = BraidGroup(4)
        sage: L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
        sage: L
        Link with 2 components represented by 8 crossings

    .. PLOT::
        :width: 300 px

        B = BraidGroup(4)
        L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
        sphinx_plot(L.plot())

    ::

        sage: L = Link(B([1, 2, 1, 3]))
        sage: L
        Link with 2 components represented by 4 crossings

    .. PLOT::
        :width: 300 px

        B = BraidGroup(4)
        L = Link(B([1, 2, 1, 3]))
        sphinx_plot(L.plot())

    We construct the "monster" unknot using a planar code, and
    then construct the oriented Gauss code and braid representation::

        sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
        ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
        ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
        sage: L.oriented_gauss_code()
        [[[1, -4, 3, -1, 10, -9, 6, -7, 8, 5, 4, -3, 2, -6, 7, -8, 9, -10, -5, -2]],
         [1, -1, 1, 1, 1, -1, -1, -1, -1, -1]]
        sage: L.braid()
        s0*s1^-3*s2^-1*s1*s3*s2^2*s1^-1*s0^-1*s2*s1^-1*s3^-1*s2*s1^-1

    .. PLOT::
        :width: 300 px

        L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
                  [17,19,8,18],[9,10,11,14],[10,12,13,11],
                  [12,19,15,13],[20,16,14,15],[16,20,17,2]])
        sphinx_plot(L.plot())

    We construct the Ochiai unknot by using an oriented Gauss code::

        sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
        ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
        ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
        sage: L.pd_code()
        [[10, 2, 11, 1], [2, 12, 3, 11], [3, 20, 4, 21], [12, 19, 13, 20],
         [21, 32, 22, 1], [31, 22, 32, 23], [9, 25, 10, 24], [4, 29, 5, 30],
         [23, 30, 24, 31], [28, 14, 29, 13], [17, 14, 18, 15], [5, 17, 6, 16],
         [15, 7, 16, 6], [7, 27, 8, 26], [25, 9, 26, 8], [18, 28, 19, 27]]

    .. PLOT::
        :width: 300 px

        L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
                    -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
                  [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
        sphinx_plot(L.plot())

    We construct the knot `7_1` and compute some invariants::

        sage: B = BraidGroup(2)
        sage: L = Link(B([1]*7))

    .. PLOT::
        :width: 300 px

        B = BraidGroup(2)
        L = Link(B([1]*7))
        sphinx_plot(L.plot())

    ::

        sage: L.alexander_polynomial()
        t^-3 - t^-2 + t^-1 - 1 + t - t^2 + t^3
        sage: L.jones_polynomial()
        -t^10 + t^9 - t^8 + t^7 - t^6 + t^5 + t^3
        sage: L.determinant()
        7
        sage: L.signature()
        -6

    The links here have removed components in which no strand is used::

        sage: B = BraidGroup(8)
        sage: b = B([1])
        sage: L = Link(b)
        sage: b.components_in_closure()
        7
        sage: L.number_of_components()
        1
        sage: L.braid().components_in_closure()
        1
        sage: L.braid().parent()
        Braid group on 2 strands

    .. WARNING::

        Equality of knots is done by comparing the corresponding braids,
        which may give false negatives.

    .. NOTE::

        The behavior of removing unused strands from an element of a
        braid group may change without notice in the future. Do not
        rely on this feature.

    .. TODO::

        Implement methods to creating new links from previously created links.
    """

    def __init__(self, data):
        r"""
        Initialize ``self``.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, -1, -1, -2,1, -2, 3, -2]))
            sage: TestSuite(L).run()
            sage: L = Link(B([1, 2, 1]))
            sage: TestSuite(L).run()
            sage: L = Link([[1, 1, 2, 2]])
            sage: TestSuite(L).run()

            sage: L = Link(B.one())
            sage: L = Link([])
            sage: L = Link([[], []])

            sage: Link([[[-1, 2, -1, 2]],  [1, 1, 1, 1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid input: data is not a valid oriented Gauss code

            sage: Link([[[-1, 2, 3, 4]]])
            Traceback (most recent call last):
            ...
            ValueError: invalid PD code: crossings must be represented by four segments

            sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 3]])
            Traceback (most recent call last):
            ...
            ValueError: invalid PD code: each segment must appear twice

        Segments in PD code must be labelled by positive integers::

            sage: code = [(2, 5, 3, 0), (4, 1, 5, 2), (0, 3, 1, 4)]
            sage: Knot(code)
            Traceback (most recent call last):
            ...
            ValueError: invalid PD code: segment label 0 not allowed

            sage: L = Link(5)
            Traceback (most recent call last):
            ...
            ValueError: invalid input: data must be either a list or a braid

        Verify that :trac:`29692` is fixed::

            sage: B = BraidGroup(5)
            sage: L = Link(B([3,4,3,-4]))
            sage: L
            Link with 1 component represented by 4 crossings
            sage: L.braid()
            s0*s1*s0*s1^-1

        PD code can be a list of 4-tuples::

            sage: code = [(2, 5, 3, 6), (4, 1, 5, 2), (6, 3, 1, 4)]
            sage: K = Knot(code); K.alexander_polynomial()
            t^-1 - 1 + t
        """
        if isinstance(data, list):
            # either oriented Gauss or PD code
            if len(data) != 2 or not all(isinstance(i, list) for i in data[0]):
                # PD code
                if any(len(i) != 4 for i in data):
                    raise ValueError("invalid PD code: crossings must be represented by four segments")
                flat = flatten(data)
                if 0 in flat:
                    raise ValueError("invalid PD code: segment label 0 not allowed")
                if any(flat.count(i) != 2 for i in set(flat)):
                    raise ValueError("invalid PD code: each segment must appear twice")
                self._pd_code = [list(vertex) for vertex in data]
                self._oriented_gauss_code = None
                self._braid = None
            else:
                # oriented Gauss code
                flat = flatten(data[0])
                if flat:
                    a, b = max(flat), min(flat)
                    if 2 * len(data[1]) != len(flat) or set(range(b, a + 1)) - set([0]) != set(flat):
                        raise ValueError("invalid input: data is not a valid oriented Gauss code")
                self._oriented_gauss_code = data
                self._pd_code = None
                self._braid = None

        else:
            from sage.groups.braid import Braid, BraidGroup
            if isinstance(data, Braid):
                # Remove all unused strands
                support = sorted(set().union(*((abs(x), abs(x) + 1) for x in data.Tietze())))
                d = {}
                for i, s in enumerate(support):
                    d[s] = i + 1
                    d[-s] = -i - 1
                if not support:
                    B = BraidGroup(2)
                else:
                    B = BraidGroup(len(support))
                self._braid = B([d[x] for x in data.Tietze()])
                self._oriented_gauss_code = None
                self._pd_code = None

            else:
                raise ValueError("invalid input: data must be either a list or a braid")

        self._mirror  = None  # set on invocation of :meth:`mirror_image`
        self._reverse = None  # set on invocation of :meth:`reverse`

    def arcs(self, presentation='pd'):
        r"""
        Return the arcs of ``self``.

        Arcs are the connected components of the planar diagram.

        INPUT:

        - ``presentation`` -- one of the following:

          * ``'pd'`` - the arcs are returned as lists of parts in the PD code
          * ``'gauss_code'`` - the arcs are returned as pieces of the Gauss
            code that start with a negative number, and end with the
            following negative one; of there exist a closed arc,
            it is returned as a list of positive numbers only

        OUTPUT:

        A list of lists representing the arcs based upon ``presentation``.

        EXAMPLES::

            sage: K = Knot([[[1,-2,3,-1,2,-3]],[1,1,1]])
            sage: K.arcs()
            [[1, 2], [3, 4], [5, 6]]
            sage: K.arcs(presentation='gauss_code')
            [[-3, 1, -2], [-2, 3, -1], [-1, 2, -3]]

        ::

            sage: L = Link([[1, 2, 3, 4], [3, 2, 1, 4]])
            sage: L.arcs()
            [[2, 4], [1], [3]]
            sage: L.arcs(presentation='gauss_code')
            [[-2, -1], [-1, -2], [2, 1]]
            sage: L.gauss_code()
            [[-1, -2], [2, 1]]
        """
        if presentation == 'pd':
            G = DiGraph()
            for e in set(flatten(self.pd_code())):
                G.add_vertex(e)
            for cr in zip(self.pd_code(), self.orientation()):
                if cr[1] == 1:
                    G.add_edge(cr[0][1], cr[0][3])
                else:
                    G.add_edge(cr[0][3], cr[0][1])
            res = []
            for S in G.connected_components_subgraphs():
                check = S.is_directed_acyclic(certificate=True)
                if check[0]:
                    source = S.sources()[0]
                    sink = S.sinks()[0]
                    res.append(S.shortest_path(source, sink))
                else:
                    res.append(check[1])
            return res
        elif presentation == 'gauss_code':
            res = []
            for comp in self.gauss_code():
                if not any(i < 0 for i in comp):
                    res.append(comp)
                else:
                    rescom = []
                    par = []
                    for i in comp:
                        par.append(i)
                        if i < 0:
                            rescom.append(copy(par))
                            par = [i]
                    rescom[0] = par + rescom[0]
                    res = res + rescom
            return res

    def fundamental_group(self, presentation='wirtinger'):
        r"""
        Return the fundamental group of the complement of ``self``.

        INPUT:

        - ``presentation`` -- string; one of the following:

          * ``'wirtinger'`` - (default) the Wirtinger presentation
            (see :wikipedia:`Link_group`)
          * ``'braid'`` - the presentation is given by the braid action
            on the free group (see chapter 2 of [Bir1975]_)

        OUTPUT:

        - a finitely presented group

        EXAMPLES::

            sage: L = Link([[1, 2, 3, 4], [3, 2, 1, 4]])
            sage: L.fundamental_group()
            Finitely presented group < x0, x1, x2 | x1*x0^-1*x2^-1*x0, x2*x0*x1^-1*x0^-1 >
            sage: L.fundamental_group('braid')
            Finitely presented group < x0, x1 | 1, 1 >

        We can see, for instance, that the  two presentations of the group
        of the figure eight knot correspond to isomorphic groups::

            sage: K8 = Knot([[[1, -2, 4, -3, 2, -1, 3, -4]], [1, 1, -1, -1]])
            sage: GA = K8.fundamental_group()
            sage: GA
            Finitely presented group < x0, x1, x2, x3 |
             x2*x0*x3^-1*x0^-1, x0*x2*x1^-1*x2^-1,
             x1*x3^-1*x2^-1*x3, x3*x1^-1*x0^-1*x1 >
            sage: GB = K8.fundamental_group(presentation='braid')
            sage: GB
            Finitely presented group < x0, x1, x2 | x1*x2^-1*x1^-1*x0*x1*x2*x1*x2^-1*x1^-1*x0^-1*x1*x2*x1^-1*x0^-1, x1*x2^-1*x1^-1*x0*x1*x2*x1^-1*x2^-1*x1^-1*x0^-1*x1*x2*x1^-1*x0*x1*x2*x1*x2^-1*x1^-1*x0^-1*x1*x2*x1^-2, x1*x2^-1*x1^-1*x0*x1*x2*x1^-1*x2^-1 >
            sage: GA.simplified()
            Finitely presented group < x0, x1 |
             x1^-1*x0*x1*x0^-1*x1*x0*x1^-1*x0^-1*x1*x0^-1 >
            sage: GB.simplified()
            Finitely presented group < x0, x2 |
             x2^-1*x0*x2^-1*x0^-1*x2*x0*x2^-1*x0*x2*x0^-1 >
        """
        from sage.groups.free_group import FreeGroup
        if presentation == 'braid':
            b = self.braid()
            F = FreeGroup(b.strands())
            rels = []
            for x in F.gens():
                rels.append(x * b / x)
            return F.quotient(rels)
        elif presentation == 'wirtinger':
            arcs = self.arcs(presentation='pd')
            F = FreeGroup(len(arcs))
            rels = []
            for crossing, orientation in zip(self.pd_code(), self.orientation()):
                a = arcs.index([i for i in arcs if crossing[0] in i][0])
                b = arcs.index([i for i in arcs if crossing[1] in i][0])
                c = arcs.index([i for i in arcs if crossing[2] in i][0])
                ela = F.gen(a)
                elb = F.gen(b)
                if orientation < 0:
                    elb = elb.inverse()
                elc = F.gen(c)
                rels.append(ela * elb / elc / elb)
            return F.quotient(rels)

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT: string representation

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L
            Link with 1 component represented by 4 crossings
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L
            Link with 3 components represented by 4 crossings
        """
        number_of_components = self.number_of_components()
        if number_of_components > 1:
            plural = 's'
        else:
            plural = ''
        pd_len = len(self.pd_code())
        return 'Link with {} component{} represented by {} crossings'.format(number_of_components, plural, pd_len)

    def __eq__(self, other):
        r"""
        Check equality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 == L2
            True
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 == L3
            False
        """
        if not isinstance(other, self.__class__):
            return False
        if self._pd_code is not None:
            if self.pd_code() == other.pd_code():
                return True
        if self._oriented_gauss_code is not None:
            if self.oriented_gauss_code() == other.oriented_gauss_code():
                return True
        return self.braid() == other.braid()

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: H = hash(L1)
        """
        return hash(self.braid())

    def __ne__(self, other):
        r"""
        Check inequality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 != L2
            False
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 != L3
            True
        """
        return not self.__eq__(other)

    def braid(self):
        r"""
        Return a braid representation of ``self``.

        OUTPUT: an element in the braid group

        .. WARNING::

            For the unknot with no crossings, this returns the identity
            of the braid group with 2 strands because this disregards
            strands with no crossings.

        EXAMPLES::

            sage: L = Link([[2, 3, 1, 4], [4, 1, 3, 2]])
            sage: L.braid()
            s^2
            sage: L = Link([[[-1, 2, -3, 1, -2, 3]], [-1, -1, -1]])
            sage: L.braid()
            s^-3
            sage: L = Link([[1,8,2,7], [8,4,9,5], [3,9,4,10], [10,1,7,6], [5,3,6,2]])
            sage: L.braid()
            (s0*s1^-1)^2*s1^-1

        TESTS::

            sage: L = Link([])
            sage: L.braid()
            1
            sage: L = Link([[], []])
            sage: L.braid()
            1

        Check that :trac:`25050` is solved::

            sage: A = Link([[[1, 2, -2, -1, -3, -4, 4, 3]], [1, 1, 1, 1]])
            sage: A.braid()
            s0*s1*s2*s3
        """
        if self._braid is not None:
            return self._braid

        from sage.groups.braid import BraidGroup
        comp = self._isolated_components()
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            b1 = L1.braid()
            b2 = L2.braid()
            n1 = b1.parent().strands()
            n2 = b2.parent().strands()
            t1 = list(b1.Tietze())
            t2 = [sign(x) * (abs(x) + n1) for x in b2.Tietze()]
            B = BraidGroup(n1 + n2)
            self._braid = B(t1 + t2)
            return self._braid

        # look for possible Vogel moves, perform them and call recursively to the modified link
        pd_code = self.pd_code()
        if not pd_code:
            B = BraidGroup(2)
            self._braid = B.one()
            return self._braid
        seifert_circles = self.seifert_circles()
        newedge = max(flatten(pd_code)) + 1
        for region in self.regions():
            n = len(region)
            for i in range(n - 1):
                a = region[i]
                seifcirca = [x for x in seifert_circles if abs(a) in x]
                for j in range(i + 1, n):
                    b = region[j]
                    seifcircb = [x for x in seifert_circles if abs(b) in x]
                    if seifcirca != seifcircb and sign(a) == sign(b):
                        tails, heads = self._directions_of_edges()

                        newPD = [list(vertex) for vertex in pd_code]
                        if sign(a) == 1:
                            C1 = newPD[newPD.index(heads[a])]
                            C1[C1.index(a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[b])]
                            C2[C2.index(b)] = newedge + 2
                            newPD.append([newedge + 3, a, b, newedge])
                            newPD.append([newedge + 2, newedge + 1, newedge + 3, newedge])
                            self._braid = Link(newPD).braid()
                            return self._braid
                        else:
                            C1 = newPD[newPD.index(heads[-a])]
                            C1[C1.index(-a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[-b])]
                            C2[C2.index(-b)] = newedge + 2
                            newPD.append([newedge + 2, newedge, newedge + 3, newedge + 1])
                            newPD.append([newedge + 3, newedge, -b, -a])
                            self._braid = Link(newPD).braid()
                            return self._braid

        # We are in the case where no Vogel moves are necessary.
        G = DiGraph()
        G.add_vertices([tuple(c) for c in seifert_circles])
        for i,c in enumerate(pd_code):
            if self.orientation()[i] == 1:
                a = [x for x in seifert_circles if c[1] in x][0]
                b = [x for x in seifert_circles if c[0] in x][0]
            else:
                a = [x for x in seifert_circles if c[0] in x][0]
                b = [x for x in seifert_circles if c[3] in x][0]
            G.add_edge(tuple(a), tuple(b))

        # Get a simple path from a source to a sink in the digraph
        it = G.all_paths_iterator(starting_vertices=G.sources(), ending_vertices=G.sinks(), simple=True)
        ordered_cycles = next(it)

        B = BraidGroup(len(ordered_cycles))
        available_crossings = copy(pd_code)
        oc_set = set(ordered_cycles[0])
        for i,x in enumerate(pd_code):
            if any(elt in oc_set for elt in x):
                crossing = x
                crossing_index = i
                break
        available_crossings.remove(crossing)
        status = [None for i in ordered_cycles]
        orientation = self.orientation()
        if orientation[crossing_index] == 1:
            b = B([1])
            status[0] = crossing[2]
            status[1] = crossing[3]
        else:
            b = B([-1])
            status[0] = crossing[1]
            status[1] = crossing[2]
        counter = 0
        while available_crossings:
            possibles = [x for x in available_crossings if status[counter] in x]
            if len(status) < counter + 2 or status[counter + 1] is not None:
                possibles = [x for x in possibles if status[counter + 1] in x]
            if possibles:
                added = possibles[0]
                if orientation[pd_code.index(added)] == 1:
                    b *= B([counter + 1])
                    status[counter] = added[2]
                    status[counter + 1] = added[3]
                else:
                    b *= B([-counter - 1])
                    status[counter] = added[1]
                    status[counter + 1] = added[2]
                if counter > 0:
                    counter -= 1
                available_crossings.remove(added)
            else:
                counter += 1
        self._braid = b
        return b

    def _directions_of_edges(self):
        r"""
        Return the directions of the edges given by the PD code of ``self``.

        OUTPUT:

        A tuple of two dictionaries. The first one assigns
        each edge of the PD code to the crossing where it starts.
        The second dictionary assigns it to where it ends.

        EXAMPLES::

            sage: L = Link([[1, 3, 2, 4], [2, 3, 1, 4]])
            sage: tails, heads = L._directions_of_edges()
            sage: tails
            {1: [2, 3, 1, 4], 2: [1, 3, 2, 4], 3: [1, 3, 2, 4], 4: [2, 3, 1, 4]}
            sage: heads
            {1: [1, 3, 2, 4], 2: [2, 3, 1, 4], 3: [2, 3, 1, 4], 4: [1, 3, 2, 4]}

        ::

            sage: L = Link([[1,5,2,4], [5,3,6,2], [3,1,4,6]])
            sage: tails, heads = L._directions_of_edges()
            sage: tails
            {1: [3, 1, 4, 6],
             2: [1, 5, 2, 4],
             3: [5, 3, 6, 2],
             4: [3, 1, 4, 6],
             5: [1, 5, 2, 4],
             6: [5, 3, 6, 2]}
            sage: heads
            {1: [1, 5, 2, 4],
             2: [5, 3, 6, 2],
             3: [3, 1, 4, 6],
             4: [1, 5, 2, 4],
             5: [5, 3, 6, 2],
             6: [3, 1, 4, 6]}

        ::

            sage: L = Link([[1,2,3,3], [2,4,5,5], [4,1,7,7]])
            sage: tails, heads = L._directions_of_edges()
            sage: tails
            {1: [4, 1, 7, 7],
             2: [1, 2, 3, 3],
             3: [1, 2, 3, 3],
             4: [2, 4, 5, 5],
             5: [2, 4, 5, 5],
             7: [4, 1, 7, 7]}
            sage: heads
            {1: [1, 2, 3, 3],
             2: [2, 4, 5, 5],
             3: [1, 2, 3, 3],
             4: [4, 1, 7, 7],
             5: [2, 4, 5, 5],
             7: [4, 1, 7, 7]}
        """
        tails = {}
        heads = {}
        pd_code = self.pd_code()
        for C in pd_code:
            tails[C[2]] = C
            a = C[2]
            D = C
            while a not in heads:
                next_crossing = [x for x in pd_code if a in x and x != D]
                if not next_crossing:
                    heads[a] = D
                    tails[a] = D
                    if D[0] == a:
                        a = D[2]
                    elif D[1] == a:
                        a = D[3]
                    else:
                        a = D[1]
                else:
                    heads[a] = next_crossing[0]
                    tails[a] = D
                    D = next_crossing[0]
                    a = D[(D.index(a)+2) % 4]

        unassigned = set(flatten(pd_code)).difference(set(tails))
        while unassigned:
            a = unassigned.pop()
            for x in pd_code:
                if a in x:
                    D = x
                    break
            while a not in heads:
                tails[a] = D
                for x in pd_code:
                    if a in x and x != D:
                        next_crossing = x
                        break
                heads[a] = next_crossing
                D = next_crossing
                a = D[(D.index(a)+2) % 4]
                if a in unassigned:
                    unassigned.remove(a)
        return tails, heads

    @cached_method
    def _enhanced_states(self):
        r"""
        Return the enhanced states of the diagram.

        Each enhanced state is represented as a tuple containing:

        - A tuple with the type of smoothing made at each crossing (0 represents
          a A-type smoothing, and 1 represents B-type).

        - A tuple with the circles marked as negative. Each circle is
          represented by the smoothings it goes through. Each smoothing
          is represented by the indices of the two strands, and the
          index of the chord, counted clockwise.

        - A tuple with the circles marked as negative.

        - The i-index (degree) corresponding to the state.

        - the j-index (height) corresponding to the state.

        EXAMPLES::

            sage: K = Link([[[1,-2,3,-1,2,-3]],[-1,-1,-1]])
            sage: K.pd_code()
            [[4, 2, 5, 1], [2, 6, 3, 5], [6, 4, 1, 3]]
            sage: K._enhanced_states()
            (((0, 0, 0),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (5, 2, 8)), ((3, 6, 9), (6, 3, 8))),
              (),
              -3,
              -9),
             ((0, 0, 0),
              (((2, 5, 7), (5, 2, 8)), ((3, 6, 9), (6, 3, 8))),
              (((1, 4, 7), (4, 1, 9)),),
              -3,
              -7),
             ((0, 0, 0),
              (((1, 4, 7), (4, 1, 9)), ((3, 6, 9), (6, 3, 8))),
              (((2, 5, 7), (5, 2, 8)),),
              -3,
              -7),
             ((0, 0, 0),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (5, 2, 8))),
              (((3, 6, 9), (6, 3, 8)),),
              -3,
              -7),
             ((0, 0, 0),
              (((3, 6, 9), (6, 3, 8)),),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (5, 2, 8))),
              -3,
              -5),
             ((0, 0, 0),
              (((2, 5, 7), (5, 2, 8)),),
              (((1, 4, 7), (4, 1, 9)), ((3, 6, 9), (6, 3, 8))),
              -3,
              -5),
             ((0, 0, 0),
              (((1, 4, 7), (4, 1, 9)),),
              (((2, 5, 7), (5, 2, 8)), ((3, 6, 9), (6, 3, 8))),
              -3,
              -5),
             ((0, 0, 0),
              (),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (5, 2, 8)), ((3, 6, 9), (6, 3, 8))),
              -3,
              -3),
             ((1, 0, 0),
              (((3, 6, 9), (6, 3, 8)), ((4, 1, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8))),
              (),
              -2,
              -7),
             ((1, 0, 0),
              (((4, 1, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8)),),
              (((3, 6, 9), (6, 3, 8)),),
              -2,
              -5),
             ((1, 0, 0),
              (((3, 6, 9), (6, 3, 8)),),
              (((4, 1, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8)),),
              -2,
              -5),
             ((1, 0, 0),
              (),
              (((3, 6, 9), (6, 3, 8)), ((4, 1, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8))),
              -2,
              -3),
             ((0, 1, 0),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (2, 6, 8), (3, 5, 8), (3, 6, 9))),
              (),
              -2,
              -7),
             ((0, 1, 0),
              (((2, 5, 7), (2, 6, 8), (3, 5, 8), (3, 6, 9)),),
              (((1, 4, 7), (4, 1, 9)),),
              -2,
              -5),
             ((0, 1, 0),
              (((1, 4, 7), (4, 1, 9)),),
              (((2, 5, 7), (2, 6, 8), (3, 5, 8), (3, 6, 9)),),
              -2,
              -5),
             ((0, 1, 0),
              (),
              (((1, 4, 7), (4, 1, 9)), ((2, 5, 7), (2, 6, 8), (3, 5, 8), (3, 6, 9))),
              -2,
              -3),
             ((1, 1, 0),
              (((2, 6, 8), (3, 5, 8), (3, 6, 9), (4, 1, 9), (4, 2, 7), (5, 1, 7)),),
              (),
              -1,
              -5),
             ((1, 1, 0),
              (),
              (((2, 6, 8), (3, 5, 8), (3, 6, 9), (4, 1, 9), (4, 2, 7), (5, 1, 7)),),
              -1,
              -3),
             ((0, 0, 1),
              (((1, 3, 9), (1, 4, 7), (6, 3, 8), (6, 4, 9)), ((2, 5, 7), (5, 2, 8))),
              (),
              -2,
              -7),
             ((0, 0, 1),
              (((2, 5, 7), (5, 2, 8)),),
              (((1, 3, 9), (1, 4, 7), (6, 3, 8), (6, 4, 9)),),
              -2,
              -5),
             ((0, 0, 1),
              (((1, 3, 9), (1, 4, 7), (6, 3, 8), (6, 4, 9)),),
              (((2, 5, 7), (5, 2, 8)),),
              -2,
              -5),
             ((0, 0, 1),
              (),
              (((1, 3, 9), (1, 4, 7), (6, 3, 8), (6, 4, 9)), ((2, 5, 7), (5, 2, 8))),
              -2,
              -3),
             ((1, 0, 1),
              (((1, 3, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8), (6, 3, 8), (6, 4, 9)),),
              (),
              -1,
              -5),
             ((1, 0, 1),
              (),
              (((1, 3, 9), (4, 2, 7), (5, 1, 7), (5, 2, 8), (6, 3, 8), (6, 4, 9)),),
              -1,
              -3),
             ((0, 1, 1),
              (((1, 3, 9), (1, 4, 7), (2, 5, 7), (2, 6, 8), (3, 5, 8), (6, 4, 9)),),
              (),
              -1,
              -5),
             ((0, 1, 1),
              (),
              (((1, 3, 9), (1, 4, 7), (2, 5, 7), (2, 6, 8), (3, 5, 8), (6, 4, 9)),),
              -1,
              -3),
             ((1, 1, 1),
              (((1, 3, 9), (3, 5, 8), (5, 1, 7)), ((2, 6, 8), (4, 2, 7), (6, 4, 9))),
              (),
              0,
              -5),
             ((1, 1, 1),
              (((2, 6, 8), (4, 2, 7), (6, 4, 9)),),
              (((1, 3, 9), (3, 5, 8), (5, 1, 7)),),
              0,
              -3),
             ((1, 1, 1),
              (((1, 3, 9), (3, 5, 8), (5, 1, 7)),),
              (((2, 6, 8), (4, 2, 7), (6, 4, 9)),),
              0,
              -3),
             ((1, 1, 1),
              (),
              (((1, 3, 9), (3, 5, 8), (5, 1, 7)), ((2, 6, 8), (4, 2, 7), (6, 4, 9))),
              0,
              -1))
        """
        writhe = self.writhe()
        crossings = self.pd_code()
        ncross = len(crossings)
        smoothings = []
        nmax = max(flatten(crossings)) + 1
        for i in range(2 ** ncross):
            v = Integer(i).bits()
            v = v + (ncross - len(v))*[0]
            G = Graph()
            for j, cr in enumerate(crossings):
                n = nmax + j
                if not v[j]:
                    # For negative crossings, we go from undercrossings to the left
                    G.add_edge((cr[3], cr[0], n), cr[0])
                    G.add_edge((cr[3], cr[0], n), cr[3])
                    G.add_edge((cr[1], cr[2], n), cr[2])
                    G.add_edge((cr[1], cr[2], n), cr[1])
                else:
                    # positive crossings, from undercrossing to the right
                    G.add_edge((cr[0], cr[1], n), cr[0])
                    G.add_edge((cr[0], cr[1], n), cr[1])
                    G.add_edge((cr[2], cr[3], n), cr[2])
                    G.add_edge((cr[2], cr[3], n), cr[3])
            sm = set(tuple(sorted(x for x in b if isinstance(x, tuple)))
                     for b in G.connected_components(sort=False))
            iindex = (writhe - ncross + 2 * sum(v)) // 2
            jmin = writhe + iindex - len(sm)
            jmax = writhe + iindex + len(sm)
            smoothings.append((tuple(v), sm, iindex, jmin, jmax))
        states = []  # we got all the smoothings, now find all the states
        for sm in smoothings:
            for k in range(len(sm[1])+1):
                for circpos in combinations(sorted(sm[1]), k):  # Add each state
                    circneg = sm[1].difference(circpos)
                    j = writhe + sm[2] + len(circpos) - len(circneg)
                    states.append((sm[0], tuple(sorted(circneg)), tuple(circpos), sm[2], j))
        return tuple(states)

    @cached_method
    def _khovanov_homology_cached(self, height, ring=ZZ):
        r"""
        Return the Khovanov homology of the link.

        INPUT:

        - ``height`` -- the height of the homology to compute
        - ``ring`` -- (default: ``ZZ``) the coefficient ring

        OUTPUT:

        The Khovanov homology of the Link in the given height. It is given
        as a tuple of key-value pairs, whose keys are the degrees.

        .. NOTE::

            This method is intended only as the cache for
            :meth:`khovanov_homology`.

        EXAMPLES::

            sage: K = Link([[[1, -2, 3, -1, 2, -3]],[-1, -1, -1]])
            sage: K._khovanov_homology_cached(-5)
            ((-3, 0), (-2, Z), (-1, 0), (0, 0))

        The figure eight knot::

            sage: L = Link([[1, 6, 2, 7], [5, 2, 6, 3], [3, 1, 4, 8], [7, 5, 8, 4]])
            sage: L._khovanov_homology_cached(-1)
            ((-2, 0), (-1, Z), (0, Z), (1, 0), (2, 0))
        """
        crossings = self.pd_code()
        ncross = len(crossings)
        states = [(_0, set(_1), set(_2), _3, _4)
                  for (_0, _1, _2, _3, _4) in self._enhanced_states()]
        bases = {}  # arrange them by (i,j)
        for st in states:
            i, j = st[3], st[4]
            if j == height:
                if (i,j) in bases:
                    bases[i,j].append(st)
                else:
                    bases[i,j] = [st]
        complexes = {}
        for (i, j) in bases:
            if (i+1, j) in bases:
                m = matrix(ring, len(bases[(i,j)]), len(bases[(i+1,j)]))
                for ii in range(m.nrows()):
                    V1 = bases[(i,j)][ii]
                    for jj in range(m.ncols()):
                        V2 = bases[(i+1, j)][jj]
                        V20 = V2[0]
                        difs = [index for index,value in enumerate(V1[0]) if value != V20[index]]
                        if len(difs) == 1 and not (V2[2].intersection(V1[1]) or V2[1].intersection(V1[2])):
                            m[ii,jj] = (-1)**sum(V2[0][x] for x in range(difs[0]+1, ncross))
                            # Here we have the matrix constructed, now we have to put it in the dictionary of complexes
            else:
                m = matrix(ring, len(bases[(i,j)]), 0)
            complexes[i] = m.transpose()
            if not (i-1, j) in bases:
                complexes[i-1] = matrix(ring, len(bases[(i,j)]), 0)
        homologies = ChainComplex(complexes).homology()
        return tuple(sorted(homologies.items()))

    def khovanov_homology(self, ring=ZZ, height=None, degree=None):
        r"""
        Return the Khovanov homology of the link.

        INPUT:

        - ``ring`` -- (default: ``ZZ``) the coefficient ring

        - ``height`` -- the height of the homology to compute,
          if not specified, all the heights are computed

        - ``degree`` -- the degree of the homology to compute,
          if not specified, all the degrees are computed

        OUTPUT:

        The Khovanov homology of the Link. It is given as a dictionary
        whose keys are the different heights. For each height, the
        homology is given as another dictionary whose keys are the degrees.

        EXAMPLES::

            sage: K = Link([[[1, -2, 3, -1, 2, -3]],[-1, -1, -1]])
            sage: K.khovanov_homology()
            {-9: {-3: Z},
             -7: {-3: 0, -2: C2},
             -5: {-3: 0, -2: Z, -1: 0, 0: 0},
             -3: {-3: 0, -2: 0, -1: 0, 0: Z},
             -1: {0: Z}}

        The figure eight knot::

            sage: L = Link([[1, 6, 2, 7], [5, 2, 6, 3], [3, 1, 4, 8], [7, 5, 8, 4]])
            sage: L.khovanov_homology(height=-1)
            {-1: {-2: 0, -1: Z, 0: Z, 1: 0, 2: 0}}

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([1, 1])
            sage: K = Link(b)
            sage: K.khovanov_homology(degree = 2)
            {2: {2: 0}, 4: {2: Z}, 6: {2: Z}}

        TESTS:

        Check that :trac:`31001` is fixed::

            sage: L = Link([])
            sage: L.khovanov_homology()
            {-1: {0: Z}, 1: {0: Z}}
            sage: L.khovanov_homology(height=-1)
            {-1: {0: Z}}
            sage: L.khovanov_homology(height=0)
            {}
            sage: L.khovanov_homology(QQ, height=1)
            {1: {0: Vector space of dimension 1 over Rational Field}}
            sage: L.khovanov_homology(GF(2), degree=0)
            {-1: {0: Vector space of dimension 1 over Finite Field of size 2},
             1: {0: Vector space of dimension 1 over Finite Field of size 2}}
            sage: L.khovanov_homology(degree=1)
            {}
            sage: L.khovanov_homology(degree=0, height=1)
            {1: {0: Z}}
            sage: L.khovanov_homology(degree=1, height=1)
            {}
        """
        if not self.pd_code():  # special case for the unknot with no crossings
            from sage.homology.homology_group import HomologyGroup
            homs = {-1: {0: HomologyGroup(1, ring, [0])},
                    1: {0: HomologyGroup(1, ring, [0])}}
            if height is not None:
                if height not in homs:
                    return {}
                homs = {height: homs[height]}
            if degree is not None:
                homs = {ht: {degree: homs[ht][degree]} for ht in homs if degree in homs[ht]}
            return homs

        if height is not None:
            heights = [height]
        else:
            heights = sorted(set(state[-1] for state in self._enhanced_states()))
        if degree is not None:
            homs = {j: dict(self._khovanov_homology_cached(j, ring)) for j in heights}
            homologies = {j: {degree: homs[j][degree]} for j in homs if degree in homs[j]}
        else:
            homologies = {j: dict(self._khovanov_homology_cached(j, ring)) for j in heights}
        return homologies

    def oriented_gauss_code(self):
        r"""
        Return the oriented Gauss code of ``self``.

        The oriented Gauss code has two parts:

        a. the Gauss code

        b. the orientation of each crossing

        The following orientation was taken into consideration for
        construction of knots:

        From the outgoing of the overcrossing if we move in the clockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `-1`.

        From the outgoing of the overcrossing if we move in the anticlockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `+1`.

        One more consideration we take in while constructing the orientation
        is the order of the orientation is same as the ordering of the
        crossings in the Gauss code.

        .. NOTE::

            Convention: under is denoted by `-1`, and over by `+1` in the
            crossing info.

        EXAMPLES::

            sage: L = Link([[1, 11, 2, 10], [6, 2, 7, 3], [3, 12, 4, 9], [9, 5, 10, 6], [8, 1, 5, 4], [11, 8, 12, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]]
            sage: L = Link([[1, 4, 2, 3], [6, 1, 3, 2], [7, 4, 8, 5], [5, 8, 6, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]]
            sage: B = BraidGroup(8)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.oriented_gauss_code()
            [[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]]

        TESTS::

            sage: L = Link([])
            sage: L.oriented_gauss_code()
            [[], []]
            sage: L = Link(BraidGroup(2).one())
            sage: L.oriented_gauss_code()
            [[], []]
        """
        if self._oriented_gauss_code is not None:
            return self._oriented_gauss_code

        pd = self.pd_code()
        orient = self.orientation()
        crossing_info = {}
        for i, j in enumerate(pd):
            if orient[i] == -1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[3], 1, i + 1)] = j[1]
            elif orient[i] == 1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[1], 1, i + 1)] = j[3]
        edges = {}
        cross_number = {}
        for i, j in crossing_info.items():
            edges[i[0]] = [j]
            if i[1] == 1:
                cross_number[i[0]] = i[2]
            elif i[1] == -1:
                cross_number[i[0]] = -i[2]
        edges_graph = DiGraph(edges)
        d = edges_graph.all_simple_cycles()
        code = []
        for i in d:
            l = []
            for j in i:
                l.append(cross_number[j])
            del l[-1]
            code.append(l)
        oriented_code = [code, orient]
        self._oriented_gauss_code = oriented_code
        return self._oriented_gauss_code

    def pd_code(self):
        r"""
        Return the planar diagram code of ``self``.

        The planar diagram is returned in the following format.

        We construct the crossing by starting with the entering component
        of the undercrossing, move in the clockwise direction and then
        generate the list. If the crossing is given by `[a, b, c, d]`,
        then we interpret this information as:

        1. `a` is the entering component of the undercrossing;
        2. `b, d` are the components of the overcrossing;
        3. `c` is the leaving component of the undercrossing.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]], [1, 1, -1, -1]])
            sage: L.pd_code()
            [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]]
            sage: B = BraidGroup(2)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.pd_code()
            [[2, 1, 3, 4], [4, 3, 5, 6], [6, 5, 7, 8], [8, 7, 9, 10], [10, 9, 1, 2]]
            sage: L = Link([[[2, -1], [1, -2]], [1, 1]])
            sage: L.pd_code()
            [[2, 3, 1, 4], [4, 1, 3, 2]]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.pd_code()
            [[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]]

        TESTS::

            sage: L = Link([[], []])
            sage: L.pd_code()
            []
            sage: L = Link(BraidGroup(2).one())
            sage: L.pd_code()
            []
        """
        if self._pd_code is not None:
            return self._pd_code

        if self._oriented_gauss_code is not None:
            oriented_gauss_code = self._oriented_gauss_code
            d_dic = {}
            if len(oriented_gauss_code[0]) > 1:
                d = flatten(oriented_gauss_code[0])
                for i, j in enumerate(d):
                    d_dic[j] = [i + 1, i + 2]
                # here we collect the final component in each Gauss code
                last_component = [i[-1] for i in oriented_gauss_code[0]]
                first_component = [i[0] for i in oriented_gauss_code[0]]
                # here we correct the last_component
                for i, j in zip(last_component, first_component):
                    d_dic[i][1] = d_dic[j][0]
                crossing_dic = {}
                for i,x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]
            elif len(oriented_gauss_code[0]) == 1:
                for i, j in enumerate(oriented_gauss_code[0][0]):
                    d_dic[j] = [i + 1, i + 2]
                d_dic[oriented_gauss_code[0][0][-1]][1] = 1
                crossing_dic = {}
                for i, x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]
            else:
                crossing_dic = {}

            pd = list(crossing_dic.values())
            self._pd_code = pd
            return self._pd_code

        if self._braid is not None:
            strings = list(range(1, self._braid.strands() + 1))
            b = list(self._braid.Tietze())
            pd = []
            strings_max = strings[-1]
            for i in b:
                if i > 0:
                    pd.append(
                        [strings[i], strings[i - 1], strings_max + 1, strings_max + 2])
                else:
                    pd.append(
                        [strings[abs(i) - 1], strings_max + 1, strings_max + 2, strings[abs(i)]])
                strings[abs(i) - 1] = strings_max + 1
                strings[abs(i)] = strings_max + 2
                strings_max = strings_max + 2
            for i in pd:
                for j in range(4):
                    if i[j] in strings:
                        i[j] = strings.index(i[j]) + 1
            self._pd_code = pd
            return pd

        raise AssertionError("invalid state")

    def gauss_code(self):
        r"""
        Return the Gauss code of ``self``.

        The Gauss code is generated by the following procedure:

        a. Number the crossings from `1` to `n`.
        b. Select a point on the knot and start moving along the component.
        c. At each crossing, take the number of the crossing, along with
           sign, which is `-` if it is an undercrossing and `+` if it is a
           overcrossing.

        EXAMPLES::

            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sage: L.gauss_code()
            [[-1, 2], [1, -2]]
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, -2, 1, -2, -2]))
            sage: L.gauss_code()
            [[-1, 3, -4, 5], [1, -2, 4, -5, 2, -3]]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L.gauss_code()
            [[-1, 2], [-3, 4], [1, 3, -4, -2]]
        """
        return self.oriented_gauss_code()[0]

    def dowker_notation(self):
        r"""
        Return the Dowker notation of ``self``.

        Similar to the PD code we number the components, so every crossing
        is represented by four numbers. We focus on the incoming entities
        of the under and the overcrossing. It is the pair of incoming
        undercrossing and the incoming overcrossing. This information at
        every crossing gives the Dowker notation.

        OUTPUT:

        A list containing the pair of incoming under cross and the incoming
        over cross.

        EXAMPLES::

            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6,-5]], [-1, -1, -1, -1, 1, -1, 1]])
            sage: L.dowker_notation()
            [(1, 6), (7, 2), (3, 10), (11, 4), (14, 5), (13, 8), (12, 9)]
            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.dowker_notation()
            [(2, 1), (3, 5), (6, 4), (7, 9)]
            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sage: L.dowker_notation()
            [(1, 3), (4, 2)]
        """
        pd = self.pd_code()
        orient = self.orientation()
        dn = [(i[0], i[3]) if orient[j] == -1 else (i[0], i[1])
              for j, i in enumerate(pd)]
        return dn

    def _braid_word_components(self):
        r"""
        Return the disjoint braid components, if any, else return the braid
        of ``self``.

        For example consider the braid ``[-1, 3, 1, 3]`` this can be viewed
        as a braid with components as ``[-1, 1]`` and ``[3, 3]``. There is no
        common crossing to these two (in sense there is a crossing between
        strand `1` and `2`, crossing between `3` and `4` but no crossing
        between strand `2` and `3`, so these can be viewed as independent
        components in the braid).

        OUTPUT: list containing the components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braid_word_components()
            ([-1, 1], [3, 3])
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braid_word_components()
            ([-1, 1, 1, 1], [3], [5, 7, 6])
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braid_word_components()
            ([-2, 1, 1], [4, 4], [6])
        """
        ml = list(self.braid().Tietze())
        if not ml:
            return tuple()

        l = set(abs(k) for k in ml)
        missing1 = set(range(min(l), max(l) + 1)) - l
        if not missing1:
            return (ml,)

        missing = sorted(missing1)
        x = [[] for i in range(len(missing) + 1)]
        for i,a in enumerate(missing):
            for j, mlj in enumerate(ml):
                if mlj != 0 and abs(mlj) < a:
                    x[i].append(mlj)
                    ml[j] = 0
                elif mlj != 0 and abs(mlj) > missing[-1]:
                    x[-1].append(mlj)
                    ml[j] = 0
        return tuple([a for a in x if a])

    def _braid_word_components_vector(self):
        r"""
        The list from the :meth:`_braid_word_components` is flattened to
        give out the vector form.

        OUTPUT: list containing braid word components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braid_word_components_vector()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braid_word_components_vector()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braid_word_components_vector()
            [-2, 1, 1, 4, 4, 6]
        """
        return flatten(self._braid_word_components())

    def _homology_generators(self):
        r"""
        The set of generators for the first homology group of the connected
        Seifert surface of the given link.

        This method uses the :meth:`_braid_word_components_vector` to generate
        the homology generators. The position of the repeated element w.r.t.
        the braid word component vector list is compiled into a list.

        This is based on Lemma 3.1 in [Col2013]_.

        OUTPUT:

        A list of integers `i \in \{1, 2, \ldots, n-1\}` corresponding
        to the simple generators `s_i` that gives a homology generator or
        `0` if the position does not represent a generator.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._homology_generators()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._homology_generators()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._homology_generators()
            [0, 2, 0, 4, 0]
        """
        x = self._braid_word_components_vector()
        hom_gen = []
        for j in range(len(x) - 1):
            a = abs(x[j])
            for i in range(j + 1, len(x)):
                if a == abs(x[i]):
                    hom_gen.append(i)
                    break
            else:
                hom_gen.append(0)
        return hom_gen

    @cached_method
    def seifert_matrix(self):
        r"""
        Return the Seifert matrix associated with ``self``.

        ALGORITHM:

        This is the algorithm presented in Section 3.3 of [Col2013]_.

        OUTPUT:

        The intersection matrix of a (not necessarily minimal) Seifert surface.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.seifert_matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.seifert_matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.seifert_matrix()
            [-1  0]
            [ 0 -1]
        """
        x = self._braid_word_components_vector()
        h = self._homology_generators()
        indices = [i for i, hi in enumerate(h) if hi]
        N = len(indices)
        A = matrix(ZZ, N, N, 0)
        for ni, i in enumerate(indices):
            hi = h[i]
            A[ni, ni] = -(x[i] + x[hi]).sign()
            for nj in range(ni + 1, N):
                j = indices[nj]
                if hi > h[j] or hi < j:
                    continue
                if hi == j:
                    if x[j] > 0:
                        A[nj, ni] = 1
                    else:
                        A[ni, nj] = -1
                elif abs(x[i]) - abs(x[j]) == 1:
                    A[nj, ni] = -1
                elif abs(x[j]) - abs(x[i]) == 1:
                    A[ni, nj] = 1
        A.set_immutable()
        return A

    @cached_method
    def number_of_components(self):
        r"""
        Return the number of connected components of ``self``.

        OUTPUT: number of connected components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.number_of_components()
            4
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.number_of_components()
            5
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.number_of_components()
            1
            sage: L = Link(B.one())
            sage: L.number_of_components()
            1
        """
        G = Graph()
        pd = self.pd_code()
        if not pd:
            return ZZ.one()
        G.add_vertices(set(flatten(pd)))
        for c in pd:
            G.add_edge(c[0], c[2])
            G.add_edge(c[1], c[3])
        return G.connected_components_number()

    def is_knot(self):
        r"""
        Return ``True`` if ``self`` is a knot.

        Every knot is a link but the converse is not true.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 3, 1, -3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        return self.number_of_components() == 1

    def genus(self):
        r"""
        Return the genus of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: L = Link(B([1,3]))
            sage: L.genus()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.genus()
            0
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.genus()
            1
        """
        b = self.braid().Tietze()
        if not b:
            return ZZ.zero()

        B = self.braid().parent()
        x = self._braid_word_components()
        q = []
        s_tmp = []
        for xi in x:
            tmp = []
            b1 = min(abs(k) for k in xi)
            for xij in xi:
                if xij > 0:
                    xij = xij - b1 + 1
                else:
                    xij = xij + b1 - 1
                tmp.append(xij)
            s_tmp.append(B(tmp))
        s = []
        for i in s_tmp:
            b = i.Tietze()
            s.append(list(b))
        t = [Link(B(si)).number_of_components() for si in s]
        for i, j in enumerate(s):
            if not j:
                s[i].append(-2)
        for i in s:
            q2 = max(abs(k) + 1 for k in i)
            q.append(q2)
        g = [((2 - t[i]) + len(x[i]) - q[i]) / 2 for i in range(len(x))]
        return sum(g, ZZ.zero())

    def signature(self):
        r"""
        Return the signature of ``self``.

        This is defined as the signature of the symmetric matrix

        .. MATH::

             V + V^{t},

        where `V` is the :meth:`Seifert matrix <seifert_matrix>`.

        .. SEEALSO:: :meth:`omega_signature`, :meth:`seifert_matrix`

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.signature()
            -1
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.signature()
            -2
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.signature()
            -2
        """
        V = self.seifert_matrix()
        m = V + V.transpose()
        return ZZ.sum(j.real().sign() for j in m.eigenvalues())

    def omega_signature(self, omega):
        r"""
        Compute the `\omega`-signature of ``self``.

        INPUT:

        - `\omega` -- a complex number of modulus 1. This is assumed to be
          coercible to ``QQbar``.

        This is defined as the signature of the Hermitian matrix

        .. MATH::

             (1 - \omega) V + (1 - \omega^{-1}) V^{t},

        where `V` is the :meth:`Seifert matrix <seifert_matrix>`,
        as explained on page 122 of [Liv1993]_.

        According to [Con2018]_, this is also known as the
        Levine-Tristram signature, the equivariant signature or the
        Tristram-Levine signature.

        .. SEEALSO:: :meth:`signature`, :meth:`seifert_matrix`

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))
            sage: omega = QQbar.zeta(3)
            sage: K.omega_signature(omega)
            -2
        """
        from sage.rings.qqbar import QQbar
        omega = QQbar(omega)
        V = self.seifert_matrix()
        m = (1 - omega) * V + (1 - omega.conjugate()) * V.transpose()
        return ZZ.sum(j.real().sign() for j in m.eigenvalues())

    def alexander_polynomial(self, var='t'):
        r"""
        Return the Alexander polynomial of ``self``.

        INPUT:

        - ``var`` -- (default: ``'t'``) the variable in the polynomial

        EXAMPLES:

        We begin by computing the Alexander polynomial for the
        figure-eight knot::

            sage: B = BraidGroup(3)
            sage: L = Link(B([1, -2, 1, -2]))
            sage: L.alexander_polynomial()
            -t^-1 + 3 - t

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.alexander_polynomial()
            1

        Some additional examples::

            sage: B = BraidGroup(2)
            sage: L = Link(B([1]))
            sage: L.alexander_polynomial()
            1
            sage: L = Link(B.one())
            sage: L.alexander_polynomial()
            1
            sage: B = BraidGroup(3)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^-1 - 1 + t

        When the Seifert surface is disconnected, the Alexander
        polynomial is defined to be `0`::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1,3]))
            sage: L.alexander_polynomial()
            0

        TESTS::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: L = Link(B([1,3,1,1,3,3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            0
        """
        R = LaurentPolynomialRing(ZZ, var)
        # The Alexander polynomial of disjoint links are defined to be 0
        if len(self._braid_word_components()) > 1:
            return R.zero()
        t = R.gen()
        seifert_matrix = self.seifert_matrix()
        f = (seifert_matrix - t * seifert_matrix.transpose()).determinant()
        # could we use a charpoly here ? or faster determinant ?
        if f != 0:
            exp = f.exponents()
            return t ** ((-max(exp) - min(exp)) // 2) * f
        return f

    def determinant(self):
        r"""
        Return the determinant of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.determinant()
            3
            sage: L = Link(B([1]*16 + [2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.determinant()
            65
            sage: B = BraidGroup(3)
            sage: Link(B([1, 2, 1, 1, 2])).determinant()
            4

        TESTS::

            sage: B = BraidGroup(3)
            sage: Link(B([1, 2, 1, -2, -1])).determinant()
            0

        REFERENCES:

        - Definition 6.6.3 in [Cro2004]_
        """
        V = self.seifert_matrix()
        m = V + V.transpose()
        return Integer(abs(m.det()))

    def is_alternating(self):
        r"""
        Return whether the given knot diagram is alternating.

        Alternating diagram implies every overcross is followed by an
        undercross or the vice-versa.

        We look at the Gauss code if the sign is alternating, ``True``
        is returned else the knot is not alternating ``False`` is returned.

        .. WARNING::

            This does not check if a knot admits an alternating diagram
            or not. Thus, this term is used differently than in some of
            the literature, such as in Hoste-Thistlethwaite table.

        .. NOTE::

            Links with more than one component are considered to not
            be alternating (knots) even when such a diagram exists.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, -1, -1, -1]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1, -2, -1, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1, 3, 1, 3, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1]*16 + [2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1,2,-1,2]))
            sage: L.is_alternating()
            True

        We give the `5_2` knot with an alternating diagram and a
        non-alternating diagram::

            sage: K5_2 = Link([[1, 4, 2, 5], [3, 8, 4, 9], [5, 10, 6, 1],
            ....:              [7, 2, 8, 3], [9, 6, 10, 7]])
            sage: K5_2.is_alternating()
            True

            sage: K5_2b = Link(K5_2.braid())
            sage: K5_2b.is_alternating()
            False

        TESTS:

        Check that :trac:`31001` is fixed::

            sage: L = Knot([])
            sage: L.is_alternating()
            True
        """
        if not self.is_knot():
            return False
        x = self.gauss_code()
        if not x:
            return True
        s = [Integer(i).sign() for i in x[0]]
        return (s == [(-1) ** (i + 1) for i in range(len(x[0]))]
                or s == [(-1) ** i for i in range(len(x[0]))])

    def orientation(self):
        r"""
        Return the orientation of the crossings of the link diagram
        of ``self``.

        EXAMPLES::

            sage: L = Link([[1, 4, 5, 2], [3, 5, 6, 7], [4, 8, 9, 6], [7, 9, 10, 11], [8, 1, 13, 10], [11, 13, 2, 3]])
            sage: L.orientation()
            [-1, 1, -1, 1, -1, 1]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.orientation()
            [-1, -1, -1, -1, 1, -1, 1]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.orientation()
            [-1, -1, -1]
        """
        directions = self._directions_of_edges()[0]
        orientation = []
        for C in self.pd_code():
            if C[0] == C[1] or C[2] == C[3]:
                orientation.append(-1)
            elif C[1] == C[2] or C[0] == C[3]:
                orientation.append(1)
            elif directions[C[1]] == C:
                orientation.append(-1)
            else:
                orientation.append(1)
        return orientation

    def seifert_circles(self):
        r"""
        Return the Seifert circles from the link diagram of ``self``.

        Seifert circles are the circles obtained by smoothing all crossings
        respecting the orientation of the segments.

        Each Seifert circle is represented as a list of the segments
        that form it.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]], [1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 7, 5, 3], [2, 6], [4, 8]]
            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]], [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[1, 13, 9, 3, 15, 5, 11, 7], [2, 10, 6, 12], [4, 16, 8, 14]]
            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6, -5]], [-1, -1, -1, -1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 11, 8], [2, 7, 12, 4, 5, 10], [3, 9, 6]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([1, 1, 1]))
            sage: L.seifert_circles()
            [[1, 3, 5], [2, 4, 6]]

        TESTS:

        Check that :trac:`25050` is solved::

            sage: A = Link([[[1, 2, -2, -1, -3, -4, 4, 3]], [1, 1, 1, 1]])
            sage: A.seifert_circles()
            [[3], [7], [1, 5], [2, 4], [6, 8]]
        """
        available_segments = set(flatten(self.pd_code()))
        result = []
        # detect looped segments. They must be their own seifert circles
        for a in available_segments:
            if any(C.count(a) > 1 for C in self.pd_code()):
                result.append([a])
        # remove the looped segments from the available
        for a in result:
            available_segments.remove(a[0])
        tails, heads = self._directions_of_edges()
        while available_segments:
            a = available_segments.pop()
            if heads[a] == tails[a]:
                result.append([a])
            else:
                C = heads[a]
                par = []
                while a not in par:
                    par.append(a)
                    posnext = C[(C.index(a) + 1) % 4]
                    if tails[posnext] == C and not [posnext] in result:
                        a = posnext
                    else:
                        a = C[(C.index(a) - 1) % 4]
                    if a in available_segments:
                        available_segments.remove(a)
                    C = heads[a]
                result.append(par)
        return result

    def regions(self):
        r"""
        Return the regions from the link diagram of ``self``.

        Regions are obtained always turning left at each crossing.

        Then the regions are represented as a list with the segments that form
        its boundary, with a sign depending on the orientation of the segment
        as part of the boundary.

        EXAMPLES::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.regions()
            [[14, -5, 12, -9], [13, 9], [11, 5, 1, 7, 3], [10, -3, 8, -13], [6, -1], [4, -11], [2, -7], [-2, -6, -14, -8], [-4, -10, -12]]
            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.regions()
            [[8, 4], [7, -4, 1], [6, -1, -3], [5, 3, -8], [2, -5, -7], [-2, -6]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.regions()
            [[16, 8, 14, 4], [15, -4], [13, -8, 1], [12, -1, -7], [11, 7, -16, 5], [10, -5, -15, -3], [9, 3, -14], [6, -11], [2, -9, -13], [-2, -12, -6, -10]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([-1, -1, -1]))
            sage: L.regions()
            [[6, -5], [5, 1, 3], [4, -3], [2, -1], [-2, -6, -4]]
            sage: L = Link([[[1, -2, 3, -4], [-1, 5, -3, 2, -5, 4]], [-1, 1, 1, -1, -1]])
            sage: L.regions()
            [[10, -4, -7], [9, 7, -3], [8, 3], [6, -9, -2], [5, 2, -8, 4], [1, -5], [-1, -10, -6]]
            sage: L = Link([[1, 2, 3, 3], [2, 5, 4, 4], [5, 7, 6, 6], [7, 1, 8, 8]])
            sage: L.regions()
            [[-3], [-4], [-6], [-8], [7, 1, 2, 5], [-1, 8, -7, 6, -5, 4, -2, 3]]

        .. NOTE::

            The link diagram is assumed to have only one completely isolated
            component. This is because otherwise some regions would have
            disconnected boundary.

        TESTS::

            sage: B = BraidGroup(6)
            sage: L = Link(B([1, 3, 5]))
            sage: L.regions()
            Traceback (most recent call last):
            ...
            NotImplementedError: can only have one isolated component
        """
        if len(self._isolated_components()) != 1:
            raise NotImplementedError("can only have one isolated component")
        pd = self.pd_code()
        if len(pd) == 1:
            if pd[0][0] == pd[0][1]:
                return [[-pd[0][2]], [pd[0][0]], [pd[0][2], -pd[0][0]]]
            else:
                return [[pd[0][2]], [-pd[0][0]], [-pd[0][2], pd[0][0]]]

        tails, heads = self._directions_of_edges()
        available_edges = set(flatten(pd))

        loops = [i for i in available_edges if heads[i] == tails[i]]
        available_edges = available_edges.union({-i for i in available_edges})
        regions = []

        for edge in loops:
            cros = heads[edge]
            if cros[1] == edge:
                regions.append([edge])
            else:
                regions.append([-edge])
            available_edges.remove(edge)
            available_edges.remove(-edge)
        available_edges = sorted(available_edges)

        while available_edges:
            edge = available_edges.pop()
            region = []
            while edge not in region:
                region.append(edge)
                if edge > 0:
                    cros = heads[edge]
                    ind = cros.index(edge)
                else:
                    cros = tails[-edge]
                    ind = cros.index(-edge)
                next_edge = cros[(ind + 1) % 4]
                if [next_edge] in regions:
                    region.append(-next_edge)
                    next_edge = cros[(ind - 1) % 4]
                elif [-next_edge] in regions:
                    region.append(next_edge)
                    next_edge = cros[(ind - 1) % 4]
                if tails[next_edge] == cros:
                    edge = next_edge
                else:
                    edge = -next_edge
                if edge in available_edges:
                    available_edges.remove(edge)
            regions.append(region)
        return regions

    @cached_method
    def mirror_image(self):
        r"""
        Return the mirror image of ``self``.

        EXAMPLES::

            sage: g = BraidGroup(2).gen(0)
            sage: K = Link(g^3)
            sage: K2 = K.mirror_image(); K2
            Link with 1 component represented by 3 crossings
            sage: K2.braid()
            s^-3

        .. PLOT::
            :width: 300 px

            g = BraidGroup(2).gen(0)
            K = Link(g**3)
            sphinx_plot(K.plot())

        .. PLOT::
            :width: 300 px

            g = BraidGroup(2).gen(0)
            K = Link(g**3)
            sphinx_plot(K.mirror_image().plot())

        ::

            sage: K = Knot([[[1, -2, 3, -1, 2, -3]], [1, 1, 1]])
            sage: K2 = K.mirror_image(); K2
            Knot represented by 3 crossings
            sage: K.pd_code()
            [[4, 1, 5, 2], [2, 5, 3, 6], [6, 3, 1, 4]]
            sage: K2.pd_code()
            [[4, 2, 5, 1], [2, 6, 3, 5], [6, 4, 1, 3]]

        .. PLOT::
            :width: 300 px

            K = Link([[[1,-2,3,-1,2,-3]],[1,1,1]])
            sphinx_plot(K.plot())

        .. PLOT::
            :width: 300 px

            K = Link([[[1,-2,3,-1,2,-3]],[1,1,1]])
            K2 = K.mirror_image()
            sphinx_plot(K2.plot())

        TESTS:

        check that :trac:`30997` is fixed::

            sage: L = Link([[6, 2, 7, 1], [5, 13, 6, 12], [8, 3, 9, 4],
            ....:           [2, 13, 3, 14], [14, 8, 15, 7], [11, 17, 12, 16],
            ....:           [9, 18, 10, 11], [17, 10, 18, 5], [4, 16, 1, 15]]) # L9n25{0}{0} from KnotInfo
            sage: Lmm = L.mirror_image().mirror_image()
            sage: L == Lmm
            True
        """
        # Use the braid information if it is the shortest version
        #   of what we have already computed
        if self._mirror:
            return self._mirror

        if self._braid:
            lb = len(self._braid.Tietze())

            if self._pd_code:
                lpd = len(self.pd_code())
            else:
                lpd = float('inf')

            if self._oriented_gauss_code:
                logc = len(self.oriented_gauss_code()[-1])
            else:
                logc = float('inf')

            if lb <= logc and lb <= lpd:
                self._mirror = type(self)(self._braid.mirror_image())
                self._mirror._mirror = self
                return self._mirror

        # Otherwise we fallback to the PD code
        pd = [[a[0], a[3], a[2], a[1]] for a in self.pd_code()]
        self._mirror = type(self)(pd)
        self._mirror._mirror = self
        return self._mirror

    def reverse(self):
        r"""
        Return the reverse of ``self``. This is the link obtained from ``self``
        by reverting the orientation on all components.

        EXAMPLES::

           sage: K3 = Knot([[5, 2, 4, 1], [3, 6, 2, 5], [1, 4, 6, 3]])
           sage: K3r = K3.reverse(); K3r.pd_code()
           [[4, 1, 5, 2], [2, 5, 3, 6], [6, 3, 1, 4]]
           sage: K3 == K3r
           True

        a non reversable knot::

           sage: K8_17 = Knot([[6, 2, 7, 1], [14, 8, 15, 7], [8, 3, 9, 4],
           ....:               [2, 13, 3, 14], [12, 5, 13, 6], [4, 9, 5, 10],
           ....:               [16, 12, 1, 11], [10, 16, 11, 15]])
           sage: K8_17r = K8_17.reverse()
           sage: b = K8_17.braid(); b
           s0^2*s1^-1*(s1^-1*s0)^2*s1^-1
           sage: br = K8_17r.braid(); br
           s0^-1*s1*s0^-2*s1^2*s0^-1*s1
           sage: b.is_conjugated(br)
           False
           sage: b == br.reverse()
           False
           sage: b.is_conjugated(br.reverse())
           True
           sage: K8_17b = Link(b)
           sage: K8_17br = K8_17b.reverse()
           sage: bbr = K8_17br.braid(); bbr
           (s1^-1*s0)^2*s1^-2*s0^2
           sage: br == bbr
           False
           sage: br.is_conjugated(bbr)
           True
        """
        if self._reverse:
            return self._reverse

        if self._braid:
            self._reverse = type(self)(self._braid.reverse())
            self._reverse._reverse = self
            return self._reverse

        # Otherwise we fallback to the PD code
        pd = [[a[2], a[3], a[0], a[1]] for a in self.pd_code()]
        self._reverse = type(self)(pd)
        self._reverse._reverse = self
        return self._reverse

    def writhe(self):
        r"""
        Return the writhe of ``self``.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.writhe()
            0
            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6,-5]],
            ....:            [-1, -1, -1, -1, 1, -1, 1]])
            sage: L.writhe()
            -3
            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
            ....:            [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.writhe()
            -2
        """
        x = self.oriented_gauss_code()
        pos = x[1].count(1)
        neg = (-1) * x[1].count(-1)
        return pos + neg

    def jones_polynomial(self, variab=None, skein_normalization=False, algorithm='jonesrep'):
        r"""
        Return the Jones polynomial of ``self``.

        The normalization is so that the unknot has Jones polynomial `1`.
        If ``skein_normalization`` is ``True``, the variable of the result
        is replaced by a itself to the power of `4`, so that the result
        agrees with the conventions of [Lic1997]_ (which in particular differs
        slightly from the conventions used otherwise in this class), had
        one used the conventional Kauffman bracket variable notation directly.

        If ``variab`` is ``None`` return a polynomial in the variable `A`
        or `t`, depending on the value ``skein_normalization``. In
        particular, if ``skein_normalization`` is ``False``, return the
        result in terms of the variable `t`, also used in [Lic1997]_.

        ALGORITHM:

        The calculation goes through one of two possible algorithms,
        depending on the value of ``algorithm``. Possible values are
        ``'jonesrep'`` which uses the Jones representation of a braid
        representation of ``self`` to compute the polynomial of the
        trace closure of the braid, and ``statesum`` which recursively
        computes the Kauffman bracket of ``self``. Depending on how the
        link is given, there might be significant time gains in using
        one over the other. When the trace closure of the braid is
        ``self``, the algorithms give the same result.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if unspecified, use either a default variable
          in `\ZZ[A,A^{-1}]` or the variable `t` in the symbolic ring

        - ``skein_normalization`` -- boolean (default: ``False``); determines
          the variable of the resulting polynomial

        - ``algorithm`` -- string (default: ``'jonesrep'``); algorithm to use
          and can be one of the following:

          * ``'jonesrep'`` - use the Jones representation of the braid
            representation

          * ``'statesum'`` - recursively computes the Kauffman bracket

        OUTPUT:

        If ``skein_normalization`` if ``False``, this returns an element
        in the symbolic ring as the Jones polynomial of the link might
        have fractional powers when the link is not a knot. Otherwise the
        result is a Laurent polynomial in ``variab``.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: Link(b).jones_polynomial()
            1

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.jones_polynomial()
            1

        The Ochiai unknot::

            sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
            ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
            ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sage: L.jones_polynomial()  # long time
            1

        Two unlinked unknots::

            sage: B = BraidGroup(4)
            sage: b = B([1, 3])
            sage: Link(b).jones_polynomial()
            -sqrt(t) - 1/sqrt(t)

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: Link(b).jones_polynomial()
            -1/sqrt(t) - 1/t^(5/2)

        Different representations of the trefoil and one of its mirror::

            sage: B = BraidGroup(2)
            sage: b = B([-1, -1, -1])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: Link(b).jones_polynomial()
            1/t + 1/t^3 - 1/t^4
            sage: B = BraidGroup(3)
            sage: b = B([-1, -2, -1, -2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: Link(b).jones_polynomial(skein_normalization=True, variab=x)
            x^-16 + x^-12 + x^-4
            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1, 2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            A^4 + A^12 - A^16

        `K11n42` (the mirror of the "Kinoshita-Terasaka" knot) and `K11n34`
        (the mirror of the "Conway" knot) in [KnotAtlas]_::

            sage: B = BraidGroup(4)
            sage: K11n42 = Link(B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2]))
            sage: K11n34 = Link(B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3]))
            sage: bool(K11n42.jones_polynomial() == K11n34.jones_polynomial())
            True

        The two algorithms for computation give the same result when the
        trace closure of the braid representation is the link itself::

            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6, -5]],
            ....:           [-1, -1, -1, -1, 1, -1, 1]])
            sage: jonesrep = L.jones_polynomial(algorithm='jonesrep')
            sage: statesum = L.jones_polynomial(algorithm='statesum')
            sage: bool(jonesrep == statesum)
            True

        When we have thrown away unknots so that the trace closure of the
        braid is not necessarily the link itself, this is only true up to a
        power of the Jones polynomial of the unknot::

            sage: B = BraidGroup(3)
            sage: b = B([1])
            sage: L = Link(b)
            sage: b.components_in_closure()
            2
            sage: L.number_of_components()
            1
            sage: b.jones_polynomial()
            -sqrt(t) - 1/sqrt(t)
            sage: L.jones_polynomial()
            1
            sage: L.jones_polynomial(algorithm='statesum')
            1

        TESTS::

            sage: L = Link([])
            sage: L.jones_polynomial(algorithm='statesum')
            1

            sage: L.jones_polynomial(algorithm='other')
            Traceback (most recent call last):
            ...
            ValueError: bad value of algorithm

        Check that :trac:`31001` is fixed::

            sage: L.jones_polynomial()
            1
        """
        if algorithm == 'statesum':
            poly = self._bracket()
            t = poly.parent().gens()[0]
            writhe = self.writhe()
            jones = poly * (-t)**(-3 * writhe)
            # Switch to the variable A to have the result agree with the output
            # of the jonesrep algorithm
            A = LaurentPolynomialRing(ZZ, 'A').gen()
            jones = jones(A**-1)

            if skein_normalization:
                if variab is None:
                    return jones
                else:
                    return jones(variab)
            else:
                if variab is None:
                    variab = 't'
                # We force the result to be in the symbolic ring because of the expand
                return jones(SR(variab)**(ZZ(1)/ZZ(4))).expand()
        elif algorithm == 'jonesrep':
            braid = self.braid()
            # Special case for the trivial knot with no crossings
            if not braid.Tietze():
                if skein_normalization:
                    return LaurentPolynomialRing(ZZ, 'A').one()
                else:
                    return SR.one()
            return braid.jones_polynomial(variab, skein_normalization)

        raise ValueError("bad value of algorithm")

    @cached_method
    def _bracket(self):
        r"""
        Return the Kaufmann bracket polynomial of the diagram of ``self``.

        Note that this is not an invariant of the link, but of the diagram.
        In particular, it is not invariant under Reidemeister I moves.

        EXAMPLES::

            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
            ....:           [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._bracket()
            -t^-10 + 2*t^-6 - t^-2 + 2*t^2 - t^6 + t^10 - t^14
            sage: L = Link([[2, 1, 3, 4], [4, 3, 1, 2]])
            sage: L._bracket()
            -t^-4 - t^4
        """
        t = LaurentPolynomialRing(ZZ, 't').gen()
        pd_code = self.pd_code()
        if not pd_code:
            return t.parent().one()
        if len(pd_code) == 1:
            if pd_code[0][0] == pd_code[0][1]:
                return -t**(-3)
            else:
                return -t**3

        cross = pd_code[0]
        rest = [list(vertex) for vertex in pd_code[1:]]
        [a, b, c, d] = cross
        if a == b and c == d and len(rest) > 0:
            return (~t + t**(-5)) * Link(rest)._bracket()
        elif a == d and c == b and len(rest) > 0:
            return (t + t**5) * Link(rest)._bracket()
        elif a == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = c
            return -t**(-3) * Link(rest)._bracket()
        elif a == d:
            for cross in rest:
                if c in cross:
                    cross[cross.index(c)] = b
            return -t**3 * Link(rest)._bracket()
        elif c == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
            return -t**3 * Link(rest)._bracket()
        elif c == d:
            for cross in rest:
                if b in cross:
                    cross[cross.index(b)] = a
            return -t**(-3) * Link(rest)._bracket()
        else:
            rest_2 = [list(vertex) for vertex in rest]
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
                if c in cross:
                    cross[cross.index(c)] = b
            for cross in rest_2:
                if d in cross:
                    cross[cross.index(d)] = c
                if b in cross:
                    cross[cross.index(b)] = a
            return t * Link(rest)._bracket() + ~t * Link(rest_2)._bracket()

    @cached_method
    def _isolated_components(self):
        r"""
        Return the PD codes of the isolated components of ``self``.

        Isolated components are links corresponding to subdiagrams that
        do not have any common crossing.

        EXAMPLES::

            sage: L = Link([[1, 1, 2, 2], [3, 3, 4, 4]])
            sage: L._isolated_components()
            [[[1, 1, 2, 2]], [[3, 3, 4, 4]]]
        """
        G = Graph()
        for c in self.pd_code():
            G.add_vertex(tuple(c))
        V = G.vertices()
        setV = [set(c) for c in V]
        for i in range(len(V) - 1):
            for j in range(i + 1, len(V)):
                if setV[i].intersection(setV[j]):
                    G.add_edge(V[i], V[j])
        return [[list(i) for i in j]
                for j in G.connected_components(sort=False)]

    @cached_method
    def homfly_polynomial(self, var1=None, var2=None, normalization='lm'):
        r"""
        Return the HOMFLY polynomial of ``self``.

        The HOMFLY polynomial `P(K)` of a link `K` is a Laurent polynomial
        in two variables defined using skein relations and for the unknot
        `U`, we have `P(U) = 1`.

        INPUT:

        - ``var1`` -- (default: ``'L'``) the first variable. If ``normalization``
          is set to ``az`` resp. ``vz`` the default is ``a`` resp. ``v``
        - ``var2`` -- (default: ``'M'``) the second variable. If ``normalization``
          is set to ``az`` resp. ``vz`` the default is ``z``
        - ``normalization`` -- (default: ``lm``) the system of coordinates
          and can be one of the following:

          * ``'lm'`` -- corresponding to the Skein relation
            `L\cdot P(K _+) + L^{-1}\cdot P(K _-) + M\cdot P(K _0) = 0`

          * ``'az'`` -- corresponding to the Skein relation
            `a\cdot P(K _+) - a^{-1}\cdot P(K _-) = z  \cdot P(K _0)`

          * ``'vz'`` -- corresponding to the Skein relation
            `v^{-1}\cdot P(K _+) - v\cdot P(K _-) = z  \cdot P(K _0)`

          where `P(K _+)`, `P(K _-)` and `P(K _0)` represent the HOMFLY
          polynomials of three links that vary only in one crossing;
          that is the positive, negative, or smoothed links respectively

        OUTPUT:

        A Laurent polynomial over the integers.

        .. NOTE::

            Use the ``'az'`` normalization to agree with the data
            in [KnotAtlas]_

            Use the ``'vz'`` normalization to agree with the data
            `KnotInfo <http://www.indiana.edu/~knotinfo/>`__.

        EXAMPLES:

        We give some examples::

            sage: g = BraidGroup(2).gen(0)
            sage: K = Knot(g^5)
            sage: K.homfly_polynomial()
            L^-4*M^4 - 4*L^-4*M^2 + 3*L^-4 - L^-6*M^2 + 2*L^-6

        The Hopf link::

            sage: L = Link([[1,3,2,4],[4,2,3,1]])
            sage: L.homfly_polynomial('x', 'y')
            -x^-1*y + x^-1*y^-1 + x^-3*y^-1

        Another version of the Hopf link where the orientation
        has been changed. Therefore we substitute `x \mapsto L^{-1}`
        and `y \mapsto M`::

            sage: L = Link([[1,4,2,3], [4,1,3,2]])
            sage: L.homfly_polynomial()
            L^3*M^-1 - L*M + L*M^-1
            sage: L = Link([[1,4,2,3], [4,1,3,2]])
            sage: L.homfly_polynomial(normalization='az')
            a^3*z^-1 - a*z - a*z^-1

        The figure-eight knot::

            sage: L = Link([[2,1,4,5], [5,6,7,3], [6,4,1,9], [9,2,3,7]])
            sage: L.homfly_polynomial()
            -L^2 + M^2 - 1 - L^-2
            sage: L.homfly_polynomial('a', 'z', 'az')
            a^2 - z^2 - 1 + a^-2

        The "monster" unknot::

            sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
            ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
            ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
            sage: L.homfly_polynomial()
            1

        Comparison with KnotInfo::

            sage: KI, m = K.get_knotinfo(); KI, m
             (<KnotInfo.K5_1: '5_1'>, False)
            sage: K.homfly_polynomial(normalization='vz') == KI.homfly_polynomial()
            True

        The knot `9_6`::

            sage: B = BraidGroup(3)
            sage: K = Knot(B([-1,-1,-1,-1,-1,-1,-2,1,-2,-2]))
            sage: K.homfly_polynomial()
            L^10*M^4 - L^8*M^6 - 3*L^10*M^2 + 4*L^8*M^4 + L^6*M^6 + L^10
             - 3*L^8*M^2 - 5*L^6*M^4 - L^8 + 7*L^6*M^2 - 3*L^6
            sage: K.homfly_polynomial('a', 'z', normalization='az')
            -a^10*z^4 + a^8*z^6 - 3*a^10*z^2 + 4*a^8*z^4 + a^6*z^6 - a^10
             + 3*a^8*z^2 + 5*a^6*z^4 - a^8 + 7*a^6*z^2 + 3*a^6

        TESTS:

        This works with isolated components::

            sage: L = Link([[[1, -1], [2, -2]], [1, 1]])
            sage: L2 = Link([[1, 3, 2, 4], [2, 3, 1, 4]])
            sage: L2.homfly_polynomial()
            -L*M^-1 - L^-1*M^-1
            sage: L.homfly_polynomial()
            -L*M^-1 - L^-1*M^-1
            sage: L.homfly_polynomial(normalization='az')
            a*z^-1 - a^-1*z^-1
            sage: L2.homfly_polynomial('α', 'ζ', 'az')
            α*ζ^-1 - α^-1*ζ^-1
            sage: L.homfly_polynomial(normalization='vz')
            -v*z^-1 + v^-1*z^-1
            sage: L2.homfly_polynomial('ν', 'ζ', 'vz')
            -ν*ζ^-1 + ν^-1*ζ^-1

        Check that :trac:`30346` is fixed::

            sage: L = Link([])
            sage: L.homfly_polynomial()
            1

        REFERENCES:

        - :wikipedia:`HOMFLY_polynomial`
        - http://mathworld.wolfram.com/HOMFLYPolynomial.html
        """
        if not var1:
            if   normalization == 'az':
                var1 = 'a'
            elif normalization == 'vz':
                var1 = 'v'
            else:
                var1 = 'L'
        if not var2:
            if  normalization == 'lm':
                var2 = 'M'
            else:
                var2 = 'z'

        L = LaurentPolynomialRing(ZZ, [var1, var2])
        if len(self._isolated_components()) > 1:
            if normalization == 'lm':
                fact = L({(1, -1):-1, (-1, -1):-1})
            elif normalization == 'az':
                fact = L({(1, -1):1, (-1, -1):-1})
            elif normalization == 'vz':
                fact = L({(1, -1):-1, (-1, -1):1})
            else:
                raise ValueError('normalization must be either `lm`, `az` or `vz`')
            fact = fact ** (len(self._isolated_components())-1)
            for i in self._isolated_components():
                fact = fact * Link(i).homfly_polynomial(var1, var2, normalization)
            return fact
        s = '{}'.format(self.number_of_components())
        ogc = self.oriented_gauss_code()
        if not ogc[0]:
            return L.one()
        for comp in ogc[0]:
            s += ' {}'.format(len(comp))
            for cr in comp:
                s += ' {} {}'.format(abs(cr)-1, sign(cr))
        for i, cr in enumerate(ogc[1]):
            s += ' {} {}'.format(i, cr)
        from sage.libs.homfly import homfly_polynomial_dict
        dic = homfly_polynomial_dict(s)
        if normalization == 'lm':
            return L(dic)
        elif normalization == 'az':
            auxdic = {}
            for a in dic:
                if (a[0] + a[1]) % 4 == 0:
                    auxdic[a] = dic[a]
                else:
                    auxdic[a] = -dic[a]
            if self.number_of_components() % 2:
                return L(auxdic)
            else:
                return -L(auxdic)
        elif normalization == 'vz':
            h_az = self.homfly_polynomial(var1=var1, var2=var2, normalization='az')
            a, z = h_az.parent().gens()
            v = ~a
            return h_az.subs({a:v})
        else:
            raise ValueError('normalization must be either `lm`, `az` or `vz`')

    def _coloring_matrix(self, n):
        r"""
        Return the coloring matrix of ``self``.

        The coloring matrix is a matrix over a prime field
        whose right kernel gives the colorings of the diagram.

        INPUT:

        - ``n`` -- the number of colors to consider

        If ``n`` is not a prime number, it is replaced by the smallest
        prime number that is larger than ``n``.

        OUTPUT:

        a matrix over the smallest prime field with cardinality
        larger than or equal to ``n``.

        EXAMPLES::

            sage: K = Link([[[1, -2, 3, -1, 2, -3]], [1, 1, 1]])
            sage: K._coloring_matrix(3)
            [2 2 2]
            [2 2 2]
            [2 2 2]
            sage: K8 = Knot([[[1, -2, 4, -3, 2, -1, 3, -4]], [1, 1, -1, -1]])
            sage: K8._coloring_matrix(4)
            [2 0 4 4]
            [4 4 2 0]
            [0 4 4 2]
            [4 2 0 4]

        REFERENCES:

        - :wikipedia:`Fox_n-coloring`
        """
        from sage.arith.misc import next_prime
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        p = next_prime(n - 1)
        F = FiniteField(p)
        arcs = self.arcs(presentation='pd')
        di = len(arcs)
        M = matrix(F, di, di)
        crossings = self.pd_code()
        for i in range(di):
            crossing = crossings[i]
            for j in range(di):
                arc = arcs[j]
                if crossing[1] in arc:
                    M[i, j] += 2
                if crossing[0] in arc:
                    M[i, j] -= 1
                if crossing[2] in arc:
                    M[i, j] -= 1
        return M

    def is_colorable(self, n):
        r"""
        Return whether the link is ``n``-colorable.

        A link is ``n``-colorable if its arcs can be painted with
        ``n`` colours, labeled from ``0`` to ``n - 1``, in such a way
        that at any crossing, the average of the indices of the
        undercrossings equals twice the index of the overcrossing.

        INPUT:

        - ``n`` -- the number of colors to consider

        If ``n`` is not a prime number, it is replaced by the smallest
        prime number that is larger than ``n``.

        EXAMPLES:

        We show that the trefoil knot is 3-colorable::

            sage: K = Link([[[1, -2, 3, -1, 2, -3]], [1, 1, 1]])
            sage: K.is_colorable(3)
            True

        But the figure eight knot is not::

            sage: K8 = Link([[[1, -2, 4, -3, 2, -1, 3, -4]], [1, 1, -1, -1]])
            sage: K8.is_colorable(3)
            False

        REFERENCES:

        - :wikipedia:`Fox_n-coloring`

        - Chapter 3 of [Liv1993]_

        .. SEEALSO:: :meth:`colorings`
        """
        return self._coloring_matrix(n).nullity() > 1

    def colorings(self, n):
        r"""
        Return the ``n``-colorings of ``self``.

        INPUT:

        - ``n`` -- the number of colors to consider

        If ``n`` is not a prime number, it is replaced by the smallest
        prime number that is larger than ``n``.

        OUTPUT:

        a list with the colorings. Each coloring is represented as
        a dictionary that maps a tuple of the edges forming each arc
        (as in the PD code) to the index of the corresponding color.

        EXAMPLES::

            sage: K = Link([[[1, -2, 3, -1, 2, -3]], [1, 1, 1]])
            sage: K.colorings(3)
            [{(1, 2): 0, (3, 4): 1, (5, 6): 2},
             {(1, 2): 0, (3, 4): 2, (5, 6): 1},
             {(1, 2): 1, (3, 4): 0, (5, 6): 2},
             {(1, 2): 1, (3, 4): 2, (5, 6): 0},
             {(1, 2): 2, (3, 4): 0, (5, 6): 1},
             {(1, 2): 2, (3, 4): 1, (5, 6): 0}]
            sage: K.pd_code()
            [[4, 1, 5, 2], [2, 5, 3, 6], [6, 3, 1, 4]]
            sage: K.arcs('pd')
            [[1, 2], [3, 4], [5, 6]]

        REFERENCES:

        - :wikipedia:`Fox_n-coloring`

        - Chapter 3 of [Liv1993]_

        .. SEEALSO:: :meth:`is_colorable`
        """
        from sage.arith.misc import next_prime
        p = next_prime(n-1)
        M = self._coloring_matrix(n)
        K = M.right_kernel()
        res = set([])
        arcs = self.arcs('pd')
        for coloring in K:
            colors = sorted(set(coloring))
            if len(colors) == p:
                colors = {b: a for a, b in enumerate(colors)}
                res.add(tuple(colors[c] for c in coloring))
        return [{tuple(arc): col for arc, col in zip(arcs, c)}
                for c in sorted(res)]

    def plot(self, gap=0.1, component_gap=0.5, solver=None,
             color='blue', **kwargs):
        r"""
        Plot ``self``.

        INPUT:

        - ``gap`` -- (default: 0.1) the size of the blank gap left for
          the crossings

        - ``component_gap`` -- (default: 0.5) the gap between isolated
          components

        - ``solver`` -- the linear solver to use, see
          :class:`~sage.numerical.mip.MixedIntegerLinearProgram`.

        - ``color`` -- (default: 'blue') a color or a coloring (as returned
          by :meth:`colorings`.

        The usual keywords for plots can be used here too.

        EXAMPLES:

        We construct the simplest version of the unknot::

            sage: L = Link([[2, 1, 1, 2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            B = BraidGroup(2)
            L = Link([[2, 1, 1, 2]])
            sphinx_plot(L.plot())

        We construct a more interesting example of the unknot::

            sage: L = Link([[2, 1, 4, 5], [3, 5, 6, 7], [4, 1, 9, 6], [9, 2, 3, 7]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[2,1,4,5], [3,5,6,7], [4,1,9,6], [9,2,3,7]])
            sphinx_plot(L.plot())

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
                      [17,19,8,18],[9,10,11,14],[10,12,13,11],
                      [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sphinx_plot(L.plot())

        The Ochiai unknot::

            sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
            ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
            ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
                        -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
                      [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sphinx_plot(L.plot())

        One of the representations of the trefoil knot::

            sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
            sage: L.plot()
            Graphics object consisting of 14 graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
            sphinx_plot(L.plot())

        The figure-eight knot::

            sage: L = Link([[2, 1, 4, 5], [5, 6, 7, 3], [6, 4, 1, 9], [9, 2, 3, 7]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[2,1,4,5], [5,6,7,3], [6,4,1,9], [9,2,3,7]])
            sphinx_plot(L.plot())

        The knot `K11n121` in [KnotAtlas]_::

            sage: L = Link([[4,2,5,1], [10,3,11,4], [5,16,6,17], [7,12,8,13],
            ....:           [18,9,19,10], [2,11,3,12], [13,20,14,21], [15,6,16,7],
            ....:           [22,18,1,17], [8,19,9,20], [21,14,22,15]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[4,2,5,1], [10,3,11,4], [5,16,6,17], [7,12,8,13],
                      [18,9,19,10], [2,11,3,12], [13,20,14,21], [15,6,16,7],
                      [22,18,1,17], [8,19,9,20], [21,14,22,15]])
            sphinx_plot(L.plot())

        One of the representations of the Hopf link::

            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sphinx_plot(L.plot())

        Plotting links with multiple isolated components::

            sage: L = Link([[[-1, 2, -3, 1, -2, 3], [4, -5, 6, -4, 5, -6]], [1, 1, 1, 1, 1, 1]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[[-1,2,-3,1,-2,3], [4,-5,6,-4,5,-6]], [1,1,1,1,1,1]])
            sphinx_plot(L.plot())

        If a coloring is passed, the different arcs are plotted with
        the corresponding colors::

            sage: B = BraidGroup(4)
            sage: b = B([1,2,3,1,2,-1,-3,2,3])
            sage: L = Link(b)
            sage: L.plot(color=L.colorings(3)[0])
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            B = BraidGroup(4)
            b = B([1, 2, 3, 1, 2, -1, -3, 2, 3])
            L = Link(b)
            sphinx_plot(L.plot(color=L.colorings(3)[0]))

        TESTS:

        Check that :trac:`20315` is fixed::

            sage: L = Link([[2,1,4,5], [5,6,7,3], [6,4,1,9], [9,2,3,7]])
            sage: L.plot(solver='GLPK')
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='Coin')    # optional - sage_numerical_backends_coin
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='CPLEX')   # optional - CPLEX
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='Gurobi')  # optional - Gurobi
            Graphics object consisting of ... graphics primitives
        """
        if type(color) is not dict:
            coloring = {int(i): color for i in set(flatten(self.pd_code()))}
        else:
            from sage.plot.colors import rainbow
            ncolors = len(set(color.values()))
            arcs = self.arcs()
            if len(color) != len(arcs):
                raise ValueError("Number of entries in the color vector must match the number of arcs")
            rainb = rainbow(ncolors)
            coloring = {int(i): rainb[color[tuple(j)]] for j in arcs for i in j}
        comp = self._isolated_components()
        # Handle isolated components individually
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            P1 = L1.plot(gap, **kwargs)
            P2 = L2.plot(gap, **kwargs)
            xtra = P1.get_minmax_data()['xmax'] + component_gap - P2.get_minmax_data()['xmin']
            for P in P2:
                if hasattr(P, 'path'):
                    for p in P.path[0]:
                        p[0] += xtra
                    for p in P.vertices:
                        p[0] += xtra
                else:
                    P.xdata = [p + xtra for p in P.xdata]
            return P1 + P2

        if 'axes' not in kwargs:
            kwargs['axes'] = False
        if 'aspect_ratio' not in kwargs:
            kwargs['aspect_ratio'] = 1

        from sage.plot.line import line
        from sage.plot.bezier_path import bezier_path
        from sage.plot.circle import circle

        # Special case for the unknot
        if not self.pd_code():
            return circle((0,0), ZZ(1)/ZZ(2), color=color, **kwargs)

        # The idea is the same followed in spherogram, but using MLP instead of
        # network flows.
        # We start by computing a way to bend the edges left or right
        # such that the resulting regions are in fact closed regions
        # with straight angles, and using the minimal number of bends.
        regions = sorted(self.regions(), key=len)
        regions = regions[:-1]
        edges = list(set(flatten(self.pd_code())))
        edges.sort()
        MLP = MixedIntegerLinearProgram(maximization=False, solver=solver)
        # v will be the list of variables in the MLP problem. There will be
        # two variables for each edge: number of right bendings and number of
        # left bendings (at the end, since we are minimizing the total, only one
        # of each will be nonzero
        v = MLP.new_variable(nonnegative=True, integer=True)
        # one condition for each region
        for i in range(len(regions)):
            cond = 0
            r = regions[i]
            for e in r:
                if e > 0:
                    cond = cond + v[2*edges.index(e)] - v[2*edges.index(e) + 1]
                else:
                    cond = cond - v[2*edges.index(-e)] + v[2*edges.index(-e) + 1]
            MLP.add_constraint(cond == 4 - len(r))
        MLP.set_objective(MLP.sum(v.values()))
        MLP.solve()
        # we store the result in a vector s packing right bends as negative left ones
        values = MLP.get_values(v, convert=ZZ, tolerance=1e-3)
        s = [values[2*i] - values[2*i + 1]
             for i in range(len(edges))]
        # segments represents the different parts of the previous edges after bending
        segments = {e: [(e,i) for i in range(abs(s[edges.index(e)])+1)] for e in edges}
        pieces = {tuple(i): [i] for j in segments.values() for i in j}
        nregions = []
        for r in regions:
            nregion = []
            for e in r:
                if e > 0:
                    rev = segments[e][:-1]
                    sig = sign(s[edges.index(e)])
                    nregion += [[a, sig] for a in rev]
                    nregion.append([segments[e][-1], 1])
                else:
                    rev = segments[-e][1:]
                    rev.reverse()
                    sig = sign(s[edges.index(-e)])
                    nregion += [[a, -sig] for a in rev]
                    nregion.append([segments[-e][0], 1])
            nregions.append(nregion)
        N = max(segments) + 1
        segments = [i for j in segments.values() for i in j]
        badregions = [nr for nr in nregions if any(-1 == x[1] for x in nr)]
        while badregions:
            badregion = badregions[0]
            a = 0
            while badregion[a][1] != -1:
                a += 1
            c = -1
            b = a
            while c != 2:
                if b == len(badregion)-1:
                    b = 0
                else:
                    b += 1
                c += badregion[b][1]
            otherregion = [nr for nr in nregions
                           if any(badregion[b][0] == x[0] for x in nr)]
            if len(otherregion) == 1:
                otherregion = None
            elif otherregion[0] == badregion:
                otherregion = otherregion[1]
            else:
                otherregion = otherregion[0]
            N1 = N
            N = N + 2
            N2 = N1 + 1
            segments.append(N1)
            segments.append(N2)
            if type(badregion[b][0]) in (int, Integer):
                segmenttoadd = [x for x in pieces
                                if badregion[b][0] in pieces[x]]
                if len(segmenttoadd) > 0:
                    pieces[segmenttoadd[0]].append(N2)
            else:
                pieces[tuple(badregion[b][0])].append(N2)

            if a < b:
                r1 = badregion[:a] + [[badregion[a][0],0], [N1,1]] + badregion[b:]
                r2 = badregion[a+1:b] + [[N2,1],[N1,1]]
            else:
                r1 = badregion[b:a] + [[badregion[a][0],0], [N1,1]]
                r2 = badregion[:b] + [[N2,1],[N1,1]] + badregion[a+1:]

            if otherregion:
                c = [x for x in otherregion if badregion[b][0] == x[0]]
                c = otherregion.index(c[0])
                otherregion.insert(c+1, [N2,otherregion[c][1]])
                otherregion[c][1] = 0
            nregions.remove(badregion)
            nregions.append(r1)
            nregions.append(r2)
            badregions = [nr for nr in nregions if any(x[1] == -1 for x in nr)]
        MLP = MixedIntegerLinearProgram(maximization=False, solver=solver)
        v = MLP.new_variable(nonnegative=True, integer=True)
        for e in segments:
            MLP.set_min(v[e], 1)
        for r in nregions:
            horp = []
            horm = []
            verp = []
            verm = []
            direction = 0
            for se in r:
                if direction % 4 == 0:
                    horp.append(v[se[0]])
                elif direction == 1:
                    verp.append(v[se[0]])
                elif direction == 2:
                    horm.append(v[se[0]])
                elif direction == 3:
                    verm.append(v[se[0]])
                if se[1] == 1:
                    direction += 1
            MLP.add_constraint(MLP.sum(horp) - MLP.sum(horm) == 0)
            MLP.add_constraint(MLP.sum(verp) - MLP.sum(verm) == 0)
        MLP.set_objective(MLP.sum(v.values()))
        MLP.solve()
        v = MLP.get_values(v)
        lengths = {piece: sum(v[a] for a in pieces[piece]) for piece in pieces}
        image = line([], **kwargs)
        crossings = {tuple(self.pd_code()[0]): (0, 0, 0)}
        availables = self.pd_code()[1:]
        used_edges = []
        ims = line([], **kwargs)
        while len(used_edges) < len(edges):
            cross_keys = list(crossings.keys())
            i = 0
            j = 0
            while cross_keys[i][j] in used_edges:
                if j < 3:
                    j += 1
                else:
                    j = 0
                    i += 1
            c = cross_keys[i]
            e = c[j]
            kwargs['color'] = coloring[e]
            used_edges.append(e)
            direction = (crossings[c][2] - c.index(e)) % 4
            orien = self.orientation()[self.pd_code().index(list(c))]
            if s[edges.index(e)] < 0:
                turn = -1
            else:
                turn = 1
            lengthse = [lengths[(e,k)] for k in range(abs(s[edges.index(e)])+1)]
            if c.index(e) == 0 or (c.index(e) == 1 and orien == 1) or (c.index(e) == 3 and orien == -1):
                turn = -turn
                lengthse.reverse()
            tailshort = (c.index(e) % 2 == 0)
            x0 = crossings[c][0]
            y0 = crossings[c][1]
            im = []
            for l in lengthse:
                if direction == 0:
                    x1 = x0 + l
                    y1 = y0
                elif direction == 1:
                    x1 = x0
                    y1 = y0 + l
                elif direction == 2:
                    x1 = x0 - l
                    y1 = y0
                elif direction == 3:
                    x1 = x0
                    y1 = y0 - l
                im.append(([[x0, y0], [x1, y1]], l, direction))
                direction = (direction + turn) % 4
                x0 = x1
                y0 = y1
            direction = (direction - turn) % 4
            c2 = [ee for ee in availables if e in ee]
            if len(c2) == 1:
                availables.remove(c2[0])
                crossings[tuple(c2[0])] = (x1, y1, (direction + c2[0].index(e) + 2) % 4)
            c2 = [ee for ee in self.pd_code() if e in ee and ee != list(c)]
            if not c2:
                headshort = not tailshort
            else:
                headshort = (c2[0].index(e) % 2 == 0)
            a = deepcopy(im[0][0])
            b = deepcopy(im[-1][0])

            def delta(u, v):
                if u < v:
                    return -gap
                if u > v:
                    return gap
                return 0

            if tailshort:
                im[0][0][0][0] += delta(a[1][0], im[0][0][0][0])
                im[0][0][0][1] += delta(a[1][1], im[0][0][0][1])
            if headshort:
                im[-1][0][1][0] -= delta(b[1][0], im[-1][0][0][0])
                im[-1][0][1][1] -= delta(b[1][1], im[-1][0][0][1])
            l = line([], **kwargs)
            c = 0
            p = im[0][0][0]
            if len(im) == 4 and max(x[1] for x in im) == 1:
                l = bezier_path([[im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]], **kwargs)
                p = im[-1][0][1]
            else:
                while c < len(im)-1:
                    if im[c][1] > 1:
                        (a, b) = im[c][0]
                        if b[0] > a[0]:
                            e = [b[0] - 1, b[1]]
                        elif b[0] < a[0]:
                            e = [b[0] + 1, b[1]]
                        elif b[1] > a[1]:
                            e = [b[0], b[1] - 1]
                        elif b[1] < a[1]:
                            e = [b[0], b[1] + 1]
                        l += line((p, e), **kwargs)
                        p = e
                    if im[c+1][1] == 1 and c < len(im) - 2:
                        xr = round(im[c+2][0][1][0])
                        yr = round(im[c+2][0][1][1])
                        xp = xr - im[c+2][0][1][0]
                        yp = yr - im[c+2][0][1][1]
                        q = [p[0] + im[c+1][0][1][0] - im[c+1][0][0][0] - xp,
                             p[1] + im[c+1][0][1][1] - im[c+1][0][0][1] - yp]
                        l += bezier_path([[p, im[c+1][0][0], im[c+1][0][1], q]], **kwargs)
                        c += 2
                        p = q
                    else:
                        if im[c+1][1] == 1:
                            q = im[c+1][0][1]
                        else:
                            q = [im[c+1][0][0][0] + sign(im[c+1][0][1][0] - im[c+1][0][0][0]),
                                 im[c+1][0][0][1] + sign(im[c+1][0][1][1] - im[c+1][0][0][1])]
                        l += bezier_path([[p, im[c+1][0][0], q]], **kwargs)
                        p = q
                        c += 1
            l += line([p, im[-1][0][1]], **kwargs)
            image += l
            ims += sum(line(a[0], **kwargs) for a in im)
        return image


    def _markov_move_cmp(self, braid):
        r"""
        Return whether ``self`` can be transformed to the closure of ``braid``
        by a sequence of Markov moves.

        More precisely it is checked whether the braid of ``self`` is conjugated
        to the given braid in the following sense. If both braids have the same
        number of strands it is checked if they are conjugated to each other in
        their common braid group (Markov move I).  If the number of strands differs,
        the braid having less strands is extended by Markov moves II (appendening
        the largest generator or its inverse recursively) until a common braid
        group can be achieved, where conjugation is tested.

        Be aware, that a negative result does not ensure that ``self`` is not
        isotopic to the closure of ``braid``.

        EXAMPLES::

            sage: b = BraidGroup(4)((1, 2, -3, 2, 2, 2, 2, 2, 2, -1, 2, 3, 2))
            sage: L = Link([[2, 1, 4, 5], [5, 4, 6, 7], [7, 6, 8, 9], [9, 8, 10, 11],
            ....:           [11, 10, 12, 13], [13, 12, 14, 15], [15, 14, 16, 17],
            ....:           [3, 17, 18, 19], [16, 1, 21, 18], [19, 21, 2, 3]])
            sage: L._markov_move_cmp(b)  # both are isotopic to ``9_3``
            True
            sage: bL = L.braid(); bL
            s0^7*s1*s0^-1*s1
            sage: Lb = Link(b); Lb
            Link with 1 component represented by 13 crossings
            sage: Lb._markov_move_cmp(bL)
            True
            sage: L == Lb
            False
            sage: b.strands() > bL.strands()
            True

        REFERENCES:

        - :wikipedia:`Markov_theorem`
        """
        sb      = self.braid()
        sb_ind  = sb.strands()

        ob      = braid
        ob_ind  = ob.strands()

        if sb_ind == ob_ind:
            return sb.is_conjugated(ob)

        if sb_ind > ob_ind:
            # if the braid of self has more strands we have to perfom
            # Markov II moves
            B = sb.parent()
            g = B.gen(ob_ind-1)
            ob = B(ob)
            if sb_ind > ob_ind+1:
                # proceed by recursion
                res = self._markov_move_cmp(ob*g)
                if not res:
                    res = self._markov_move_cmp(ob*~g)
            else:
                res = sb.is_conjugated(ob*g)
                if not res:
                    res = sb.is_conjugated(ob*~g)
            return res
        else:
            L = Link(ob)
            return L._markov_move_cmp(sb)

    @cached_method
    def _knotinfo_matching_list(self):
        r"""
        Return a list of links from the KnotInfo and LinkInfo databases which match
        the properties of ``self`` as much as possible.

        OUTPUT:

        A tuple ``(l, proved)`` where ``l`` is the matching list and ``proved`` a boolean
        telling if the entries of ``l`` are checked to be isotopic to self or not.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: KnotInfo.L5a1_0.inject()
            Defining L5a1_0
            sage: L5a1_0.link()._knotinfo_matching_list()
            ([<KnotInfo.L5a1_0: 'L5a1{0}'>], True)
            sage: Link(L5a1_0.braid())._knotinfo_matching_list()
            ([<KnotInfo.L5a1_0: 'L5a1{0}'>, <KnotInfo.L5a1_1: 'L5a1{1}'>], True)

        Care is needed for links having non irreducible HOMFLY-PT polynomials::

            sage: k4_1 = KnotInfo.K4_1.link()
            sage: k5_2 = KnotInfo.K5_2.link()
            sage: k = k4_1.connected_sum(k5_2)
            sage: k._knotinfo_matching_list()   # optional - database_knotinfo
            ([<KnotInfo.K9_12: '9_12'>], False)

        """
        from sage.knots.knotinfo import KnotInfoSeries
        cr  = len(self.pd_code())
        co  = self.number_of_components()

        # set the limits for the KnotInfoSeries
        if cr > 11 and co > 1:
            cr = 11
        if cr > 12:
            cr = 12

        Hp  = self.homfly_polynomial(normalization='vz')

        det = None
        if cr > 6:
            # for larger crossing numbers the KnotInfoSeries become very
            # large, as well. For better performance we restrict the cached
            # lists by the determinant and number of components.
            #
            # Since :meth:`determinant` is not implemented for proper links
            # we have to go back to the roots.
            ap = self.alexander_polynomial()
            det = Integer(abs(ap(-1)))

        is_knot = self.is_knot()
        if is_knot and cr < 11:
            S = KnotInfoSeries(cr, True, None)
            l = S.lower_list(oriented=True, comp=co, det=det, homfly=Hp)
        else:
            # the result of :meth:`is_alternating` depends on the specific
            # diagram of the link. For example ``K11a_2`` is an alternating
            # knot but ``Link(KnotInfo.K11a_2.braid()).is_alternating()``
            # gives ``False``. Therefore, we have to take both series
            # into consideration.
            Sa = KnotInfoSeries(cr, is_knot, True)
            Sn = KnotInfoSeries(cr, is_knot, False)
            la = Sa.lower_list(oriented=True, comp=co, det=det, homfly=Hp)
            ln = Sn.lower_list(oriented=True, comp=co, det=det, homfly=Hp)
            l = sorted(list(set(la + ln)))

        pdm = [[a[0], a[3], a[2], a[1]] for a in self.pd_code() ]
        br  = self.braid()
        br_ind = br.strands()

        res = []
        for L in l:
            if L.pd_notation() == pdm:
                # note that KnotInfo pd_notation works counter clockwise. Therefore,
                # to compensate this we compare with the mirrored pd_code. See also,
                # docstring of :meth:`link` of :class:`~sage.knots.knotinfo.KnotInfoBase`.
                return[L], True  # pd_notation is unique in the KnotInfo database

            if L.braid_index() <= br_ind:
                if self._markov_move_cmp(L.braid()):
                    res.append(L)

        if res:
            if len(res) > 1 or res[0].is_unique():
                return res, True
        return l, False

    def get_knotinfo(self, mirror_version=True, unique=True):
        r"""
        Identify this link as an item of the KnotInfo database (if possible).

        INPUT:

        - ``mirror_version`` -- boolean (default is ``True``). If set to ``False``
          the result of the method will be just the instance of :class:`~sage.knots.knotinfo.KnotInfoBase`
          (by default the result is a tuple of the instance and a boolean, see
          explanation of the output below)

        - ``unique`` -- boolean (default is ``True``). This only affects the case
          where a unique identification is not possible. If set to ``False`` you
          can obtain a matching list (see explanation of the output below)

        OUTPUT:

        A tuple ``(K, m)`` where ``K`` is an instance of :class:`~sage.knots.knotinfo.KnotInfoBase`
        and ``m`` a boolean (for chiral links) telling if ``self`` corresponds
        to the mirrored version of ``K`` or not. The value of ``m`` is ``None``
        for amphicheiral links and ``?`` if it cannot be determined uniquely
        and the keyword option ``unique=False`` is given.

        For proper links, if the orientation mutant cannot be uniquely determined,
        K will be a series of links gathering all links having the same unoriented
        name, that is an instance of :class:`~sage.knots.knotinfo.KnotInfoSeries`.

        If ``mirror_version`` is set to ``False`` then the result is just ``K``
        (that is: ``m`` is suppressed).

        If it is not possible to determine a unique result a ``NotImplementedError``
        will be raised. To avoid this you can set ``unique`` to ``False``. You
        will get a list of matching candidates instead.

        .. NOTE::

            The identification of proper links may fail to be unique due to the
            following fact: In opposite to the database for knots, there are pairs
            of oriented mutants of an unoriented link which are isotopic to each
            other. For example ``L5a1_0`` and ``L5a1_1`` is such a pair.

            This is because all combinatorial possible oriented mutants are
            listed with individual names regardless whether they are pairwise
            non isotopic or not. In such a case the identification is not
            unique and therefore a series of the links will be returned which
            gathers all having the same unoriented name.

            To obtain the individual oriented links being isotopic to ``self``
            use the keyword ``unique`` (see the examples for ``L2a1_1`` and
            ``L5a1_0`` below).

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = Link([[4,2,5,1], [10,3,11,4], [5,16,6,17], [7,12,8,13],
            ....:           [18,9,19,10], [2,11,3,12], [13,20,14,21], [15,6,16,7],
            ....:           [22,18,1,17], [8,19,9,20], [21,14,22,15]])
            sage: L.get_knotinfo()           # optional - database_knotinfo
            (<KnotInfo.K11n_121: '11n_121'>, True)

            sage: K = KnotInfo.K10_25        # optional - database_knotinfo
            sage: l = K.link()               # optional - database_knotinfo
            sage: l.get_knotinfo()           # optional - database_knotinfo
            (<KnotInfo.K10_25: '10_25'>, False)

        Knots with more than 12 and proper links having more than 11 crossings
        cannot be identified. In addition non prime links or even links whose
        HOMFLY-PT polynomial is not irreducible cannot be identified::

            sage: b, = BraidGroup(2).gens()
            sage: Link(b**13).get_knotinfo()
            Traceback (most recent call last):
            ...
            NotImplementedError: this knot having more than 12 crossings cannot be determined

            sage: Link([[1, 5, 2, 4], [3, 1, 4, 8], [5, 3, 6, 2], [6, 9, 7, 10], [10, 7, 9, 8]])
            Link with 2 components represented by 5 crossings
            sage: _.get_knotinfo()
            Traceback (most recent call last):
            ...
            NotImplementedError: this (possibly non prime) link cannot be determined

        Lets identify the monster unknot::

            sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
            ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
            ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
            sage: L.get_knotinfo()
            (<KnotInfo.K0_1: '0_1'>, None)

        Usage of option ``mirror_version``::

            sage: L.get_knotinfo(mirror_version=False) == KnotInfo.K0_1
            True

        Usage of option ``unique``::

            sage: l = K.link(K.items.gauss_notation)  # optional - database_knotinfo
            sage: l.get_knotinfo()                    # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: this link cannot be uniquely determined
            use keyword argument `unique` to obtain more details

            sage: l.get_knotinfo(unique=False)        # optional - database_knotinfo
            [(<KnotInfo.K10_25: '10_25'>, False), (<KnotInfo.K10_56: '10_56'>, False)]

            sage: k11  = KnotInfo.K11n_82.link()      # optional - database_knotinfo
            sage: k11m = k11.mirror_image()           # optional - database_knotinfo
            sage: k11mr = k11m.reverse()              # optional - database_knotinfo
            sage: k11mr.get_knotinfo()                # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: mirror type of this link cannot be uniquely determined
            use keyword argument `unique` to obtain more details

            sage: k11mr.get_knotinfo(unique=False)    # optional - database_knotinfo
            [(<KnotInfo.K11n_82: '11n_82'>, '?')]

            sage: t = (1, -2, 1, 1, -2, 1, -2, -2)
            sage: l8 = Link(BraidGroup(3)(t))
            sage: l8.get_knotinfo()                   # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: this link cannot be uniquely determined
            use keyword argument `unique` to obtain more details

            sage: l8.get_knotinfo(unique=False)       # optional - database_knotinfo
            [(<KnotInfo.L8a19_0_0: 'L8a19{0,0}'>, None),
             (<KnotInfo.L8a19_1_1: 'L8a19{1,1}'>, None)]

            sage: t = (2, -3, -3, -2, 3, 3, -2, 3, 1, -2, -2, 1)
            sage: l12 = Link(BraidGroup(5)(t))
            sage: l12.get_knotinfo()                  # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: this link having more than 11 crossings cannot be uniquely determined
            use keyword argument `unique` to obtain more details

            sage: l12.get_knotinfo(unique=False)      # optional - database_knotinfo
            [(<KnotInfo.L10n36_0: 'L10n36{0}'>, '?'),
             (<KnotInfo.L10n36_1: 'L10n36{1}'>, None),
             (<KnotInfo.L10n59_0: 'L10n59{0}'>, None),
             (<KnotInfo.L10n59_1: 'L10n59{1}'>, None)]

        Furthermore, if the result is a complete  series of oriented links having
        the same unoriented name (according to the note above) the option can be
        used to achieve more detailed information::

            sage: L2a1 = Link(b**2)
            sage: L2a1.get_knotinfo()
            (Series of links L2a1, None)
            sage: L2a1.get_knotinfo(unique=False)
            [(<KnotInfo.L2a1_0: 'L2a1{0}'>, True), (<KnotInfo.L2a1_1: 'L2a1{1}'>, False)]

            sage: KnotInfo.L5a1_0.inject()
            Defining L5a1_0
            sage: l5 = Link(L5a1_0.braid())
            sage: l5.get_knotinfo()
            (Series of links L5a1, False)
            sage: _[0].inject()
            Defining L5a1
            sage: list(L5a1)
            [<KnotInfo.L5a1_0: 'L5a1{0}'>, <KnotInfo.L5a1_1: 'L5a1{1}'>]
            sage: l5.get_knotinfo(unique=False)
            [(<KnotInfo.L5a1_0: 'L5a1{0}'>, False), (<KnotInfo.L5a1_1: 'L5a1{1}'>, False)]

        Clarifying the series around the Perko pair (:wikipedia:`Perko_pair`)::

            sage: for i in range(160, 166):           # optional - database_knotinfo
            ....:     K = Knots().from_table(10, i)
            ....:     print('%s_%s' %(10, i), '--->', K.get_knotinfo())
            10_160 ---> (<KnotInfo.K10_160: '10_160'>, False)
            10_161 ---> (<KnotInfo.K10_161: '10_161'>, True)
            10_162 ---> (<KnotInfo.K10_162: '10_162'>, False)
            10_163 ---> (<KnotInfo.K10_163: '10_163'>, False)
            10_164 ---> (<KnotInfo.K10_164: '10_164'>, False)
            10_165 ---> (<KnotInfo.K10_165: '10_165'>, True)

        Clarifying ther Perko series against `SnapPy
        <https://snappy.math.uic.edu/index.html>`__::

            sage: import snappy                    # optional - snappy
            Plink failed to import tkinter.

            sage: from sage.knots.knotinfo import KnotInfoSeries
            sage: KnotInfoSeries(10, True, True)   # optional - database_knotinfo
            Series of knots K10
            sage: _.inject()                       # optional - database_knotinfo
            Defining K10
            sage: for i in range(160, 166):        # optional - database_knotinfo snappy
            ....:     K = K10(i)
            ....:     k = K.link(K.items.name, snappy=True)
            ....:     print(k, '--->', k.sage_link().get_knotinfo())
            <Link 10_160: 1 comp; 10 cross> ---> (<KnotInfo.K10_160: '10_160'>, False)
            <Link 10_161: 1 comp; 10 cross> ---> (<KnotInfo.K10_161: '10_161'>, True)
            <Link 10_162: 1 comp; 10 cross> ---> (<KnotInfo.K10_161: '10_161'>, False)
            <Link 10_163: 1 comp; 10 cross> ---> (<KnotInfo.K10_162: '10_162'>, False)
            <Link 10_164: 1 comp; 10 cross> ---> (<KnotInfo.K10_163: '10_163'>, False)
            <Link 10_165: 1 comp; 10 cross> ---> (<KnotInfo.K10_164: '10_164'>, False)

            sage: snappy.Link('10_166')            # optional - snappy
            <Link 10_166: 1 comp; 10 cross>
            sage: _.sage_link().get_knotinfo()     # optional - database_knotinfo snappy
            (<KnotInfo.K10_165: '10_165'>, True)

        Another pair of confusion (see the corresponding `Warning
        <http://katlas.math.toronto.edu/wiki/10_86>`__)::

           sage: Ks10_86 = snappy.Link('10_86')     # optional - snappy
           sage: Ks10_83 = snappy.Link('10_83')     # optional - snappy
           sage: Ks10_86.sage_link().get_knotinfo() # optional - snappy
           (<KnotInfo.K10_83: '10_83'>, True)
           sage: Ks10_83.sage_link().get_knotinfo() # optional - snappy
           (<KnotInfo.K10_86: '10_86'>, False)

        TESTS:

            sage: L = KnotInfo.L10a171_1_1_0         # optional - database_knotinfo
            sage: l = L.link(L.items.braid_notation) # optional - database_knotinfo
            sage: l.get_knotinfo(unique=False)       # optional - database_knotinfo
            [(<KnotInfo.L10a171_0_1_0: 'L10a171{0,1,0}'>, True),
             (<KnotInfo.L10a171_1_0_1: 'L10a171{1,0,1}'>, True),
             (<KnotInfo.L10a171_1_1_0: 'L10a171{1,1,0}'>, False),
             (<KnotInfo.L10a171_1_1_1: 'L10a171{1,1,1}'>, False)]
        """
        # ToDo: extension to non prime links in which case an element of the monoid
        # over :class:`KnotInfo` should be returned

        non_unique_hint = '\nuse keyword argument `unique` to obtain more details'
        def answer(L):
            r"""
            Return a single item of the KnotInfo database according to the keyword
            arguments ``mirror_version``.
            """
            if not mirror_version:
                return L

            chiral = True
            ach = L.is_amphicheiral()
            achp = L.is_amphicheiral(positive=True)
            if ach is None and achp is None:
                if unique:
                    raise NotImplementedError('this link cannot be uniquely determined (unknown chirality)%s' %non_unique_hint)
            elif  L.is_amphicheiral() or L.is_amphicheiral(positive=True):
                chiral = False

            if not chiral:
                mirrored = None
            elif proved_m and not proved_s and L in lm:
                mirrored = True
            elif proved_s and not proved_m and L in l:
                mirrored = False
            else:
                # nothing proved
                if L in ls and L in lm:
                    # In case of a chiral link this means that the HOMFLY-PT
                    # polynomial does not distinguish mirror images (see the above
                    # example ``k11m``).
                    if unique:
                        raise NotImplementedError('mirror type of this link cannot be uniquely determined%s' %non_unique_hint)
                    mirrored = '?'
                elif L in lm:
                    mirrored = True
                else:
                    mirrored = False

            return L, mirrored

        def answer_unori(S):
            r"""
            Return a series of oriented links having the same unoriented name
            according to the keyword ``mirror_version``.
            """
            if not mirror_version:
                return S

            mirrored = [answer(L)[1] for L in S]
            if all(mirrored):
                # all matching links are mirrored to self
                return S, True
            if any(i == '?' for i in mirrored):
                # unknown chirality for a matching link
                return S, '?'
            if any(i is None for i in mirrored):
                # an amphicheiral link matches
                return S, None
            if not any(mirrored):
                # no matching link is mirrored to self
                return S, False
            # finally both mirror types match
            return S, None

        def answer_list(l):
            r"""
            Return a list of items of the KnotInfo database according to the keyword
            argument ``unique``.
            """
            if not unique:
                return sorted([answer(L) for L in l])

            if len(l) == 1:
                return answer(l[0])

            if not l[0].is_knot():
                S = l[0].series(oriented=True)
                if set(list(S)) == set(l):
                    return answer_unori(S)

            raise NotImplementedError('this link cannot be uniquely determined%s' %non_unique_hint)


        self_m = self.mirror_image()
        ls, proved_s = self._knotinfo_matching_list()
        lm, proved_m = self_m._knotinfo_matching_list()
        l = list(set(ls + lm))

        if l and not unique:
            return answer_list(l)

        if proved_s and proved_m:
            return answer_list(l)

        if proved_s:
            return answer_list(ls)

        if proved_m:
            return answer_list(lm)

        # here we come if we cannot be sure about the found result

        uniq_txt = ('', '')
        if l:
            uniq_txt = (' uniquely', non_unique_hint)

        cr = len(self.pd_code())
        if self.is_knot() and cr > 12:
            # we cannot not be sure if this link is recorded in the KnotInfo database
            raise NotImplementedError('this knot having more than 12 crossings cannot be%s determined%s' %uniq_txt)

        if not self.is_knot() and cr > 11:
            # we cannot not be sure if this link is recorded in the KnotInfo database
            raise NotImplementedError('this link having more than 11 crossings cannot be%s determined%s' %uniq_txt)

        H = self.homfly_polynomial(normalization='vz')

        if sum(exp for f, exp in H.factor()) > 1:
            # we cannot be sure if this is a prime link (see the example for the connected
            # sum of K4_1 and K5_2 in the doctest of :meth:`_knotinfo_matching_list`)
            raise NotImplementedError('this (possibly non prime) link cannot be%s determined%s' %uniq_txt)

        if not l:
            from sage.features.databases import DatabaseKnotInfo
            DatabaseKnotInfo().require()
            return l

        return answer_list(l)


    def is_isotopic(self, other):
        r"""
        Check whether ``self`` is isotopic to ``other``.

        INPUT:

        - ``other`` -- another instance of :class:`Link`

        EXAMPLES::

            sage: l1 = Link([[2, 9, 3, 10], [4, 13, 5, 14], [6, 11, 7, 12],
            ....:            [8, 1, 9, 2], [10, 7, 11, 8], [12, 5, 13, 6],
            ....:            [14, 3, 1, 4]])
            sage: l2 = Link([[1, 8, 2, 9], [9, 2, 10, 3], [3, 14, 4, 1],
            ....:            [13, 4, 14, 5], [5, 12, 6, 13], [11, 6, 12, 7],
            ....:            [7, 10, 8, 11]])
            sage: l1.is_isotopic(l2)
            True

            sage: l3 = l2.mirror_image()
            sage: l1.is_isotopic(l3)
            False

            sage: from sage.knots.knotinfo import KnotInfo
            sage: L = KnotInfo.L7a7_0_0             # optional - database_knotinfo
            sage: L.series(oriented=True).inject()  # optional - database_knotinfo
            Defining L7a7
            sage: L == L7a7(0)                      # optional - database_knotinfo
            True
            sage: l = L.link()                      # optional - database_knotinfo
            sage: l.is_isotopic(L7a7(1).link())     # optional - database_knotinfo
            Traceback (most recent call last):
            ...
            NotImplementedError: comparison not possible!
            sage: l.is_isotopic(L7a7(2).link())     # optional - database_knotinfo
            True
            sage: l.is_isotopic(L7a7(3).link())     # optional - database_knotinfo
            False
        """
        from sage.misc.verbose import verbose
        if not isinstance(other, Link):
            verbose('other is not a link')
            return False

        if self == other:
            # surely isotopic
            verbose('identified by representation')
            return True

        if self.homfly_polynomial() != other.homfly_polynomial():
            # surely non isotopic
            verbose('different Homfly-PT polynomials')
            return False

        if self._markov_move_cmp(other.braid()):
            # surely isotopic
            verbose('identified via Markov moves')
            return True

        try:
            ki, m = self.get_knotinfo()
            verbose('KnotInfo self: %s mirrored %s' %(ki, m))
            try:
                if ki.is_unique():
                    try:
                        kio = other.get_knotinfo()
                        verbose('KnotInfo other: %s mirrored %s' %kio)
                        return  (ki, m) == kio
                    except NotImplementedError:
                        pass
            except AttributeError:
                # ki is a series
                pass
        except NotImplementedError:
            pass

        raise NotImplementedError('comparison not possible!')
