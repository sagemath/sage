r"""
Diagram and Partition Algebras

AUTHORS:

- Mike Hansen (2007): Initial version
- Stephen Doty, Aaron Lauve, George H. Seelinger (2012): Implementation of
  partition, Brauer, Temperley--Lieb, and ideal partition algebras
"""

#*****************************************************************************
#  Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                2012 Stephen Doty <doty@math.luc.edu>,
#                     Aaron Lauve <lauve@math.luc.edu>,
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.all import FiniteDimensionalAlgebrasWithBasis
from sage.structure.element import generic_power
from sage.combinat.free_module import (CombinatorialFreeModule,
    CombinatorialFreeModuleElement)
from sage.combinat.set_partition import SetPartitions, SetPartition
from sage.sets.set import Set
from sage.graphs.graph import Graph
from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ
import math

def partition_diagrams(k):
    r"""
    Return a list of all partition diagrams of order ``k``.

    A partition diagram of order `k \in \ZZ` to is a set partition of
    `\{1, \dots, k, -1, \ldots, -k\}`. If we have `k - 1/2 \in ZZ`, then
    a partition diagram of order `k \in 1/2 \ZZ` is a set partition of
    `\{1, \ldots, k+1/2, -1, \ldots, -(k+1/2)\}` with `k+1/2` and `-(k+1/2)`
    in the same block. See [HR2005]_.

    INPUT:

    - ``k`` -- the order of the partition diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.partition_diagrams(2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}},
         {{-2}, {-1, 1, 2}}, {{-2, 1, 2}, {-1}}, {{-2, 1}, {-1, 2}},
         {{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}},
         {{-2}, {-1, 2}, {1}}, {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}},
         {{-2, 1}, {-1}, {2}}, {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: da.partition_diagrams(3/2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 2}, {-1, 1}},
         {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    if k in ZZ:
        return SetPartitions( range(1, k+1) + [-j for j in range(1, k+1)] ).list()
    # Else k in 1/2 ZZ
    L = []
    k += ZZ(1) / ZZ(2)
    for sp in SetPartitions( range(1, k+1) + [-j for j in range(1, k)] ):
        sp = list(sp)
        for i in range(len(sp)):
            if k in sp[i]:
                sp[i] += Set([-k])
                break
        L.append(SetPartition(sp))
    return L

def brauer_diagrams(k):
    r"""
    Return a list of all Brauer diagrams of order ``k``.

    A Brauer diagram of order `k` is a partition diagram of order `k`
    with block size 2.

    INPUT:

     - ``k`` -- the order of the Brauer diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.brauer_diagrams(2)
        [{{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}]
        sage: da.brauer_diagrams(5/2)
        [{{-3, 3}, {-2, 1}, {-1, 2}}, {{-3, 3}, {-2, 2}, {-1, 1}}, {{-3, 3}, {-2, -1}, {1, 2}}]
    """
    if k in ZZ:
        return [SetPartition(list(x)) for x in
                SetPartitions( range(1,k+1) + [-j for j in range(1,k+1)],
                               [2 for j in range(1,k+1)] )]
    # Else k in 1/2 ZZ
    L = []
    k += ZZ(1) / ZZ(2)
    for i in SetPartitions( range(1, k) + [-j for j in range(1, k)],
                            [2 for j in range(1, k)] ):
        L.append(SetPartition(list(i) + [Set([k, -k])]))
    return L

def temperley_lieb_diagrams(k):
    r"""
    Return a list of all Temperley--Lieb diagrams of order ``k``.

    A Temperley--Lieb diagram of order `k` is a partition diagram of order `k`
    with block size  2 and is planar.

    INPUT:

    - ``k`` -- the order of the Temperley--Lieb diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.temperley_lieb_diagrams(2)
        [{{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}]
        sage: da.temperley_lieb_diagrams(5/2)
        [{{-3, 3}, {-2, 2}, {-1, 1}}, {{-3, 3}, {-2, -1}, {1, 2}}]
    """
    B = brauer_diagrams(k)
    T = []
    for i in B:
        if is_planar(i) == True:
            T.append(i)
    return T

def planar_diagrams(k):
    r"""
    Return a list of all planar diagrams of order ``k``.

    A planar diagram of order `k` is a partition diagram of order `k`
    that has no crossings.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.planar_diagrams(2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}},
         {{-2}, {-1, 1, 2}}, {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1, 1}},
         {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}}, {{-2}, {-1, 2}, {1}},
         {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}}, {{-2, 1}, {-1}, {2}},
         {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: da.planar_diagrams(3/2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 2}, {-1, 1}},
         {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    A = partition_diagrams(k)
    P = []
    for i in A:
        if is_planar(i) == True:
            P.append(i)
    return P

def ideal_diagrams(k):
    r"""
    Return a list of all "ideal" diagrams of order ``k``.

    An ideal diagram of order `k` is a partition diagram of order `k` with
    propagating number less than `k`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.ideal_diagrams(2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}}, {{-2}, {-1, 1, 2}},
         {{-2, 1, 2}, {-1}}, {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}},
         {{-2}, {-1, 2}, {1}}, {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}}, {{-2, 1},
         {-1}, {2}}, {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: da.ideal_diagrams(3/2)
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    A = partition_diagrams(k)
    I = []
    for i in A:
        if propagating_number(i) < k:
            I.append(i)
    return I

class DiagramAlgebra(CombinatorialFreeModule):
    r"""
    Abstract class for diagram algebras and is not designed to be used
    directly. If used directly, the class could create an "algebra"
    that is not actually an algebra.

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: R.<x> = QQ[]
        sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
        sage: sorted(D.basis())
        [P[{{-2}, {-1}, {1}, {2}}],
         P[{{-2}, {-1}, {1, 2}}],
         P[{{-2}, {-1, 1}, {2}}],
         P[{{-2}, {-1, 1, 2}}],
         P[{{-2}, {-1, 2}, {1}}],
         P[{{-2, -1}, {1}, {2}}],
         P[{{-2, -1}, {1, 2}}],
         P[{{-2, -1, 1}, {2}}],
         P[{{-2, -1, 1, 2}}],
         P[{{-2, -1, 2}, {1}}],
         P[{{-2, 1}, {-1}, {2}}],
         P[{{-2, 1}, {-1, 2}}],
         P[{{-2, 1, 2}, {-1}}],
         P[{{-2, 2}, {-1}, {1}}],
         P[{{-2, 2}, {-1, 1}}]]
    """
    def __init__(self, k, q, base_ring, prefix, diagrams, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``k`` -- the rank
        - ``q`` -- the deformation parameter
        - ``base_ring`` -- the base ring
        - ``prefix`` -- the prefix of our monomials
        - ``diagrams`` -- the *function* which will generate all diagrams
          (i.e. indices for the basis elements)

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: TestSuite(D).run()
        """
        self._prefix = prefix
        self._q = base_ring(q)
        self._k = k
        if category is None:
            category = FiniteDimensionalAlgebrasWithBasis(base_ring)
        CombinatorialFreeModule.__init__(self, base_ring, diagrams(k),
                    category=category, prefix=prefix)

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: sp = da.to_set_partition( [[1,2], [-1,-2]] )
            sage: b_elt = D(sp); b_elt
            P[{{-2, -1}, {1, 2}}]
            sage: b_elt in D
            True
            sage: D([[1,2],[-1,-2]]) == b_elt
            True
            sage: D([{1,2},{-1,-2}]) == b_elt
            True
        """
        if set_partition in self.basis().keys():
            return CombinatorialFreeModule._element_constructor_(self, set_partition)

        sp = SetPartition(set_partition) # attempt conversion
        if sp in self.basis().keys():
            return self.basis()[sp]

        raise ValueError("invalid input of {0}".format(set_partition))

    def __getitem__(self, i):
        """
        Get the basis item of ``self`` indexed by ``i``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: sp = da.to_set_partition( [[1,2], [-1,-2]] )
            sage: D[sp]
            P[{{-2, -1}, {1, 2}}]
        """
        i = to_set_partition(i)
        if i in self.basis().keys():
            return self.basis()[i]
        raise ValueError("{0} is not an index of a basis element".format(i))

    def order(self):
        r"""
        Return the order of ``self``.

        The order of a partition algebra is defined as half of the number
        of nodes in the diagrams.

        EXAMPLES::

            sage: q = var('q')
            sage: PA = PartitionAlgebra(2, q)
            sage: PA.order()
            2
        """
        return self._k

    def set_partitions(self):
        r"""
        Return the collection of underlying set partitions indexing the
        basis elements of a given diagram algebra.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: list(D.set_partitions()) == da.partition_diagrams(2)
            True
        """
        return self.basis().keys()

    def product_on_basis(self, d1, d2):
        r"""
        Returns the product `D_{d_1} D_{d_2}` by two basis diagrams.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: sp = SetPartition([[1,2],[-1,-2]])
            sage: D.product_on_basis(sp, sp)
            x*P[{{-2, -1}, {1, 2}}]
        """
        (composite_diagram, loops_removed) = set_partition_composition(d1, d2)
        return self.term(composite_diagram, self._q**loops_removed)

    @cached_method
    def one_basis(self):
        r"""
        The following constructs the identity element of the diagram algebra.

        It is not called directly; instead one should use ``DA.one()`` if
        ``DA`` is a defined diagram algebra.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.partition_diagrams)
            sage: D.one_basis()
            {{-2, 2}, {-1, 1}}
        """
        return identity_set_partition(self._k)

    # The following subclass provides a few additional methods for
    # partition algebra elements.
    class Element(CombinatorialFreeModuleElement):
        r"""
        This subclass provides a few additional methods for
        partition algebra elements. Most element methods are
        already implemented elsewhere.
        """
        def diagram(self):
            r"""
            Return the underlying diagram of ``self`` if ``self`` is a basis
            element. Raises an error if ``self`` is not a basis element.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: elt = 3*P([[1,2],[-2,-1]])
                sage: elt.diagram()
                {{-2, -1}, {1, 2}}
            """
            if len(self) != 1:
                raise ValueError("this is only defined for basis elements")
            PA = self.parent()
            ans = self.support_of_term()
            if ans not in partition_diagrams(PA.order()):
                raise ValueError("element should be keyed by a diagram")
            return ans

        def diagrams(self):
            r"""
            Return the diagrams in the support of ``self``.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: elt = 3*P([[1,2],[-2,-1]]) + P([[1,2],[-2], [-1]])
                sage: elt.diagrams()
                [{{-2}, {-1}, {1, 2}}, {{-2, -1}, {1, 2}}]
            """
            return self.support()

        def _latex_(self):
            r"""
            Return `\LaTeX` representation of ``self`` to draw single
            diagrams in latex using tikz.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: latex(P([[1,2],[-2,-1]]))
                \begin{tikzpicture}[scale = 0.9,thick]
                \tikzstyle{vertex} = [shape = circle, minimum size = 9pt, inner sep = 1pt]
                \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
                \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
                \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
                \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
                \draw (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
                \draw (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
                \end{tikzpicture}
            """
            # these allow the view command to work (maybe move them somewhere more appropriate?)
            from sage.misc.latex import latex
            latex.add_to_mathjax_avoid_list('tikzpicture')
            latex.add_package_to_preamble_if_available('tikz')

            # Define the sign function
            def sgn(x):
                if x > 0:
                    return 1
                if x < 0:
                    return -1
                return 0
            diagram = self.diagram()
            l1 = [] #list of blocks
            l2 = [] #lsit of nodes
            for i in list(diagram):
                l1.append(list(i))
                for j in list(i):
                    l2.append(j)
            #setup beginning of picture
            output = "\\begin{tikzpicture}[scale = 0.9,thick] \n\\tikzstyle{vertex} = [shape = circle, minimum size = 9pt, inner sep = 1pt] \n"
            for i in l2: #add nodes
                output = output + "\\node[vertex] (G-%s) at (%s, %s) [shape = circle, draw] {}; \n" % (i, (abs(i)-1)*1.5, sgn(i))
            for i in l1: #add edges
                if (len(i) > 1):
                    l4 = list(i)
                    posList = []
                    negList = []
                    for i in l4: #sort list so rows are grouped together
                        if i > 0:
                            posList.append(i)
                        elif i < 0:
                            negList.append(i)
                    posList.sort()
                    negList.sort()
                    l4 = posList + negList
                    l5 = l4[:] #deep copy
                    for j in range(len(l5)):
                        l5[j-1] = l4[j] #create a permuted list
                    if (len(l4) == 2):
                        l4.pop()
                        l5.pop() #pops to prevent duplicating edges
                    for j in zip(l4, l5):
                        xdiff = abs(j[1])-abs(j[0])
                        y1 = sgn(j[0])
                        y2 = sgn(j[1])
                        if ((y2-y1) == 0 and abs(xdiff) < 5): #if nodes are close to each other on same row
                            diffCo = (0.5+0.1*(abs(xdiff)-1)) #gets bigger as nodes are farther apart; max value of 1; min value of 0.5.
                            outVec = (sgn(xdiff)*diffCo, -1*diffCo*y1)
                            inVec = (-1*diffCo*sgn(xdiff), -1*diffCo*y2)
                        elif ((y2-y1) != 0 and abs(xdiff) == 1): #if nodes are close enough curviness looks bad.
                            outVec = (sgn(xdiff)*0.75, -1*y1)
                            inVec = (-1*sgn(xdiff)*0.75, -1*y2)
                        else:
                            outVec = (sgn(xdiff)*1, -1*y1)
                            inVec = (-1*sgn(xdiff), -1*y2)
                        output = output + "\\draw (G-%s) .. controls +%s and +%s .. (G-%s); \n" % (j[0], outVec, inVec,j[1])
            output = output + "\\end{tikzpicture}" #end picture
            return output

class PartitionAlgebra(DiagramAlgebra):
    r"""
    A partition algebra.

    A partition algebra of rank `k` over a given ground ring `R` is an
    algebra with (`R`-module) basis indexed by the collection of set
    partitions of `\{1, \ldots, k, -1, \ldots, -k\}`. Each such set
    partition can be represented by a graph on nodes `\{1, \ldots, k, -1,
    \ldots, -k\}` arranged in two rows, with nodes `1, \dots, k` in the
    top row from left to right and with nodes `-1, \ldots, -k` in the
    bottom row from left to right, and edges drawn such that the connected
    components of the graph are precisely the parts of the set partition.
    (This choice of edges is often not unique, and so there are often many
    graphs representing one and the same set partition; the representation
    nevertheless is useful and vivid. We often speak of "diagrams" to mean
    graphs up to such equivalence of choices of edges; of course, we could
    just as well speak of set partitions.)

    There is not just one partition algebra of given rank over a given
    ground ring, but rather a whole family of them, indexed by the
    elements of `R`. More precisely, for every `q \in R`, the partition
    algebra of rank `k` over `R` with parameter `q` is defined to be the
    `R`-algebra with basis the collection of all set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}`, where the product of two basis
    elements is given by the rule

    .. MATH::

        a \cdot b = q^N (a \circ b),

    where `a \circ b` is the composite set partition obtained by placing
    the diagram (i.e., graph) of `a` above the diagram of `b`, identifying
    the bottom row nodes of `a` with the top row nodes of `b`, and
    omitting any closed "loops" in the middle. The number `N` is the
    number of connected components formed by the omitted loops.

    The parameter `q` is a deformation parameter. Taking `q = 1` produces
    the semigroup algebra (over the base ring) of the partition monoid,
    in which the product of two set partitions is simply given by their
    composition.

    The Iwahori--Hecke algebra of type `A` (with a single parameter) is
    naturally a subalgebra of the partition algebra.

    The partition algebra is regarded as an example of a "diagram algebra"
    due to the fact that its natural basis is given by certain graphs
    often called diagrams.

    An excellent reference for partition algebras and their various
    subalgebras (Brauer algebra, Temperley--Lieb algebra, etc) is the
    paper [HR2005]_.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if
      ``None``, then Sage automatically chooses the parent of ``q``

    - ``prefix`` -- (default ``"P"``) a label for the basis elements

    EXAMPLES:

    The following shorthand simultaneously defines the univariate polynomial
    ring over the rationals as well as the variable ``x``::

        sage: R.<x> = PolynomialRing(QQ)
        sage: R
        Univariate Polynomial Ring in x over Rational Field
        sage: x
        x
        sage: x.parent() is R
        True

    We now define the partition algebra of rank `2` with parameter ``x``
    over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: P = PartitionAlgebra(2, x, R)
        sage: P
        Partition Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: P.basis().list()
        [P[{{-2, -1, 1, 2}}], P[{{-2, -1, 2}, {1}}],
         P[{{-2, -1, 1}, {2}}], P[{{-2}, {-1, 1, 2}}],
         P[{{-2, 1, 2}, {-1}}], P[{{-2, 1}, {-1, 2}}],
         P[{{-2, 2}, {-1, 1}}], P[{{-2, -1}, {1, 2}}],
         P[{{-2, -1}, {1}, {2}}], P[{{-2}, {-1, 2}, {1}}],
         P[{{-2, 2}, {-1}, {1}}], P[{{-2}, {-1, 1}, {2}}],
         P[{{-2, 1}, {-1}, {2}}], P[{{-2}, {-1}, {1, 2}}],
         P[{{-2}, {-1}, {1}, {2}}]]
        sage: E = P([[1,2],[-2,-1]]); E
        P[{{-2, -1}, {1, 2}}]
        sage: E in P.basis()
        True
        sage: E^2
        x*P[{{-2, -1}, {1, 2}}]
        sage: E^5
        x^4*P[{{-2, -1}, {1, 2}}]
        sage: (P([[2,-2],[-1,1]]) - 2*P([[1,2],[-1,-2]]))^2
        (4*x-4)*P[{{-2, -1}, {1, 2}}] + P[{{-2, 2}, {-1, 1}}]

    One can work with partition algebras using a symbol for the parameter,
    leaving the base ring unspecified. This implies that the underlying
    base ring is Sage's symbolic ring.

    ::

        sage: q = var('q')
        sage: PA = PartitionAlgebra(2, q); PA
        Partition Algebra of rank 2 with parameter q over Symbolic Ring
        sage: PA([[1,2],[-2,-1]])^2 == q*PA([[1,2],[-2,-1]])
        True
        sage: (PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2 == (4*q-4)*PA([[1, 2], [-2, -1]]) + PA([[2, -2], [1, -1]])
        True

    The identity element of the partition algebra is the set
    partition `\{\{1,-1\}, \{2,-2\}, \ldots, \{k,-k\}\}`::

        sage: P = PA.basis().list()
        sage: PA.one()
        P[{{-2, 2}, {-1, 1}}]
        sage: PA.one()*P[7] == P[7]
        True
        sage: P[7]*PA.one() == P[7]
        True

    We now give some further examples of the use of the other arguments.
    One may wish to "specialize" the parameter to a chosen element of
    the base ring::

        sage: R.<q> = RR[]
        sage: PA = PartitionAlgebra(2, q, R, prefix='B')
        sage: PA
        Partition Algebra of rank 2 with parameter q over
         Univariate Polynomial Ring in q over Real Field with 53 bits of precision
        sage: PA([[1,2],[-1,-2]])
        1.00000000000000*B[{{-2, -1}, {1, 2}}]
        sage: PA = PartitionAlgebra(2, 5, base_ring=ZZ, prefix='B')
        sage: PA
        Partition Algebra of rank 2 with parameter 5 over Integer Ring
        sage: (PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2 == 16*PA([[-2, -1], [1, 2]]) + PA([[2, -2], [1, -1]])
        True

    TESTS:

    A computation that returned an incorrect result until :trac:`15958`::

        sage: A = PartitionAlgebra(1,17)
        sage: g = SetPartitionsAk(1).list()
        sage: a = A[g[1]]
        sage: a
        P[{{-1}, {1}}]
        sage: a*a
        17*P[{{-1}, {1}}]

    REFERENCES:

    .. [HR2005] Tom Halverson and Arun Ram, *Partition algebras*. European
       Journal of Combinatorics **26** (2005), 869--921.
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="P"):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PA1 = PartitionAlgebra(2, q)
            sage: PA2 = PartitionAlgebra(2, q, R, 'P')
            sage: PA1 is PA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super(PartitionAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    # The following is the basic constructor method for the class.
    # The purpose of the "prefix" is to label the basis elements
    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PA = PartitionAlgebra(2, q, R)
            sage: TestSuite(PA).run()
        """
        self._k = k
        self._prefix = prefix
        self._q = base_ring(q)
        DiagramAlgebra.__init__(self, k, q, base_ring, prefix, partition_diagrams)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: PartitionAlgebra(2, q, R)
            Partition Algebra of rank 2 with parameter q over Univariate Polynomial Ring in q over Rational Field
        """
        return "Partition Algebra of rank %s with parameter %s over %s"%(self._k,
                self._q, self.base_ring())

class SubPartitionAlgebra(DiagramAlgebra):
    """
    A subalgebra of the partition algebra indexed by a subset of the diagrams.
    """
    def __init__(self, k, q, base_ring, prefix, diagrams, category=None):
        """
        Initialize ``self`` by adding a coercion to the ambient space.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: BA.ambient().has_coerce_map_from(BA)
            True
        """
        DiagramAlgebra.__init__(self, k, q, base_ring, prefix, diagrams, category)
        amb = self.ambient()
        self.module_morphism(self.lift, codomain=amb,
                             category=self.category()).register_as_coercion()

    #These methods allow for a sub algebra to be correctly identified in a partition algebra
    def ambient(self):
        r"""
        Return the partition algebra ``self`` is a sub-algebra of.
        Generally, this method is not called directly.

        EXAMPLES::

            sage: x = var('x')
            sage: BA = BrauerAlgebra(2, x)
            sage: BA.ambient()
            Partition Algebra of rank 2 with parameter x over Symbolic Ring
        """
        return PartitionAlgebra(self._k, self._q, self.base_ring(), prefix=self._prefix)

    def lift(self, x):
        r"""
        Lift a diagram subalgebra element to the corresponding element
        in the ambient space. This method is not intended to be called
        directly.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: E = BA([[1,2],[-1,-2]])
            sage: lifted = BA.lift(E); lifted
            B[{{-2, -1}, {1, 2}}]
            sage: lifted.parent() is BA.ambient()
            True
        """
        if x not in self:
            raise ValueError("{0} is not in {1}".format(x, self))
        monomial_indices = x.support()
        monomial_coefficients = x.coefficients()
        result = 0
        for i in xrange(len(monomial_coefficients)):
            result += monomial_coefficients[i]*self.ambient().monomial(monomial_indices[i])
        return result

    def retract(self, x):
        r"""
        Retract an appropriate partition algebra element to the
        corresponding element in the partition subalgebra. This method
        is not intended to be called directly.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: PA = BA.ambient()
            sage: E = PA([[1,2], [-1,-2]])
            sage: BA.retract(E) in BA
            True
        """
        if x not in self.ambient() or not set(x.support()).issubset(set(self.basis().keys())):
            raise ValueError("{0} cannot retract to {1}".format(x, self))
        monomial_indices = x.support()
        monomial_coefficients = x.coefficients()
        result = self.zero()
        for i in xrange(len(monomial_coefficients)):
            result += monomial_coefficients[i]*self.monomial(monomial_indices[i])
        return result

class BrauerAlgebra(SubPartitionAlgebra):
    r"""
    A Brauer algebra.

    The Brauer algebra of rank `k` is an algebra with basis indexed by the
    collection of set partitions of `\{1, \ldots, k, -1, \ldots, -k\}`
    with block size 2.

    This algebra is a subalgebra of the partition algebra.
    For more information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"B"``) a label for the basis elements

    EXAMPLES:

    We now define the Brauer algebra of rank `2` with parameter ``x`` over
    `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: B = BrauerAlgebra(2, x, R)
        sage: B
        Brauer Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: B.basis()
        Finite family {{{-2, -1}, {1, 2}}: B[{{-2, -1}, {1, 2}}], {{-2, 1}, {-1, 2}}: B[{{-2, 1}, {-1, 2}}], {{-2, 2}, {-1, 1}}: B[{{-2, 2}, {-1, 1}}]}
        sage: b = B.basis().list()
        sage: b
        [B[{{-2, 1}, {-1, 2}}], B[{{-2, 2}, {-1, 1}}], B[{{-2, -1}, {1, 2}}]]
        sage: b[2]
        B[{{-2, -1}, {1, 2}}]
        sage: b[2]^2
        x*B[{{-2, -1}, {1, 2}}]
        sage: b[2]^5
        x^4*B[{{-2, -1}, {1, 2}}]
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="B"):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: BA1 = BrauerAlgebra(2, q)
            sage: BA2 = BrauerAlgebra(2, q, R, 'B')
            sage: BA1 is BA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super(BrauerAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q, R)
            sage: TestSuite(BA).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, brauer_diagrams)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BrauerAlgebra(2, q, R)
            Brauer Algebra of rank 2 with parameter q over Univariate Polynomial Ring in q over Rational Field
        """
        return "Brauer Algebra of rank %s with parameter %s over %s"%(self._k, self._q, self.base_ring())

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q, R)
            sage: sp = SetPartition([[1,2], [-1,-2]])
            sage: b_elt = BA(sp); b_elt
            B[{{-2, -1}, {1, 2}}]
            sage: b_elt in BA
            True
            sage: BA([[1,2],[-1,-2]]) == b_elt
            True
            sage: BA([{1,2},{-1,-2}]) == b_elt
            True
        """
        set_partition = to_Brauer_partition(set_partition, k = self.order())
        return DiagramAlgebra._element_constructor_(self, set_partition)

class TemperleyLiebAlgebra(SubPartitionAlgebra):
    r"""
    A Temperley--Lieb algebra.

    The Temperley--Lieb algebra of rank `k` is an algebra with basis
    indexed by the collection of planar set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}` with block size 2.

    This algebra is thus a subalgebra of the partition algebra.
    For more information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"T"``) a label for the basis elements

    EXAMPLES:

    We define the Temperley--Lieb algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: T = TemperleyLiebAlgebra(2, x, R); T
        Temperley-Lieb Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: T.basis()
        Finite family {{{-2, 2}, {-1, 1}}: T[{{-2, 2}, {-1, 1}}], {{-2, -1}, {1, 2}}: T[{{-2, -1}, {1, 2}}]}
        sage: b = T.basis().list()
        sage: b
        [T[{{-2, 2}, {-1, 1}}], T[{{-2, -1}, {1, 2}}]]
        sage: b[1]
        T[{{-2, -1}, {1, 2}}]
        sage: b[1]^2 == x*b[1]
        True
        sage: b[1]^5 == x^4*b[1]
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="T"):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: T1 = TemperleyLiebAlgebra(2, q)
            sage: T2 = TemperleyLiebAlgebra(2, q, R, 'T')
            sage: T1 is T2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super(TemperleyLiebAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``

        TESTS::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(2, q, R)
            sage: TestSuite(TL).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, temperley_lieb_diagrams)

    def _repr_(self):
        """
        Return a string represetation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TemperleyLiebAlgebra(2, q, R)
            Temperley-Lieb Algebra of rank 2 with parameter q over Univariate Polynomial Ring in q over Rational Field
        """
        return "Temperley-Lieb Algebra of rank %s with parameter %s over %s"%(self._k,
                self._q, self.base_ring())

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(2, q, R)
            sage: sp = SetPartition([[1,2], [-1,-2]])
            sage: b_elt = TL(sp); b_elt
            T[{{-2, -1}, {1, 2}}]
            sage: b_elt in TL
            True
            sage: TL([[1,2],[-1,-2]]) == b_elt
            True
            sage: TL([{1,2},{-1,-2}]) == b_elt
            True
        """
        set_partition = to_Brauer_partition(set_partition, k = self.order())
        return SubPartitionAlgebra._element_constructor_(self, set_partition)

class PlanarAlgebra(SubPartitionAlgebra):
    """
    A planar algebra.

    The planar algebra of rank `k` is an algebra with basis indexed by the
    collection of all planar set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}`.

    This algebra is thus a subalgebra of the partition algebra. For more
    information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"Pl"``) a label for the basis elements

    EXAMPLES:

    We define the planar algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: Pl = PlanarAlgebra(2, x, R); Pl
        Planar Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: Pl.basis().list()
        [Pl[{{-2, -1, 1, 2}}], Pl[{{-2, -1, 2}, {1}}],
         Pl[{{-2, -1, 1}, {2}}], Pl[{{-2}, {-1, 1, 2}}],
         Pl[{{-2, 1, 2}, {-1}}], Pl[{{-2, 2}, {-1, 1}}],
         Pl[{{-2, -1}, {1, 2}}], Pl[{{-2, -1}, {1}, {2}}],
         Pl[{{-2}, {-1, 2}, {1}}], Pl[{{-2, 2}, {-1}, {1}}],
         Pl[{{-2}, {-1, 1}, {2}}], Pl[{{-2, 1}, {-1}, {2}}],
         Pl[{{-2}, {-1}, {1, 2}}], Pl[{{-2}, {-1}, {1}, {2}}]]
        sage: E = Pl([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="Pl"):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: Pl1 = PlanarAlgebra(2, q)
            sage: Pl2 = PlanarAlgebra(2, q, R, 'Pl')
            sage: Pl1 is Pl2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super(PlanarAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PlA = PlanarAlgebra(2, q, R)
            sage: TestSuite(PlA).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, planar_diagrams)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: Pl = PlanarAlgebra(2, x, R); Pl
            Planar Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        """
        return "Planar Algebra of rank %s with parameter %s over %s"%(self._k,
                self._q, self.base_ring())

class PropagatingIdeal(SubPartitionAlgebra):
    r"""
    A propagating ideal.

    The propagating ideal of rank `k` is a non-unital algebra with basis
    indexed by the collection of ideal set partitions of `\{1, \ldots, k, -1,
    \ldots, -k\}`. We say a set partition is *ideal* if its propagating
    number is less than `k`.

    This algebra is a non-unital subalgebra and an ideal of the partition
    algebra.
    For more information, see :class:`PartitionAlgebra`.

    EXAMPLES:

    We now define the propagating ideal of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = QQ[]
        sage: I = PropagatingIdeal(2, x, R); I
        Propagating Ideal of rank 2 with parameter x over Univariate Polynomial Ring in x over Rational Field
        sage: I.basis().list()
        [I[{{-2, -1, 1, 2}}], I[{{-2, -1, 2}, {1}}],
         I[{{-2, -1, 1}, {2}}], I[{{-2}, {-1, 1, 2}}],
         I[{{-2, 1, 2}, {-1}}], I[{{-2, -1}, {1, 2}}],
         I[{{-2, -1}, {1}, {2}}], I[{{-2}, {-1, 2}, {1}}],
         I[{{-2, 2}, {-1}, {1}}], I[{{-2}, {-1, 1}, {2}}],
         I[{{-2, 1}, {-1}, {2}}], I[{{-2}, {-1}, {1, 2}}],
         I[{{-2}, {-1}, {1}, {2}}]]
        sage: E = I([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="I"):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: IA1 = PropagatingIdeal(2, q)
            sage: IA2 = PropagatingIdeal(2, q, R, 'I')
            sage: IA1 is IA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super(PropagatingIdeal, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: I = PropagatingIdeal(2, q, R)
            sage: TestSuite(I).run() # Not tested -- needs non-unital algebras category
        """
        # This should be the category of non-unital fin-dim algebras with basis
        category = FiniteDimensionalAlgebrasWithBasis(base_ring)
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, ideal_diagrams, category)

    @cached_method
    def one_basis(self):
        r"""
        The propagating ideal is a non-unital algebra, i.e. it does not have a
        multiplicative identity.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: I = PropagatingIdeal(2, q, R)
            sage: I.one_basis()
            Traceback (most recent call last):
            ...
            ValueError: The ideal partition algebra is not unital
            sage: I.one()
            Traceback (most recent call last):
            ...
            ValueError: The ideal partition algebra is not unital
        """
        raise ValueError("The ideal partition algebra is not unital")
        #return identity_set_partition(self._k)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: PropagatingIdeal(2, x, R)
            Propagating Ideal of rank 2 with parameter x over Univariate Polynomial Ring in x over Rational Field
        """
        return "Propagating Ideal of rank %s with parameter %s over %s"%(self._k,
                self._q, self.base_ring())

    class Element(SubPartitionAlgebra.Element):
        """
        Need to take care of exponents since we are not unital.
        """
        def __pow__(self, n):
            """
            Return ``self`` to the `n`-th power.

            INPUT::

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: I = PropagatingIdeal(2, x, R)
                sage: E = I([[1,2],[-1,-2]])
                sage: E^2
                x*I[{{-2, -1}, {1, 2}}]
                sage: E^0
                Traceback (most recent call last):
                ...
                ValueError: can only take positive integer powers
            """
            if n <= 0:
                raise ValueError("can only take positive integer powers")
            return generic_power(self, n)

#########################################################################
# START BORROWED CODE
#########################################################################
# Borrowed from Mike Hansen's original code -- global methods for dealing
# with partition diagrams, compositions of partition diagrams, and so on.
# --> CHANGED 'identity' to 'identity_set_partition' for enhanced clarity.
#########################################################################

def is_planar(sp):
    r"""
    Return ``True`` if the diagram corresponding to the set partition ``sp``
    is planar; otherwise, return ``False``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.is_planar( da.to_set_partition([[1,-2],[2,-1]]))
        False
        sage: da.is_planar( da.to_set_partition([[1,-1],[2,-2]]))
        True
    """
    to_consider = map(list, sp)

    #Singletons don't affect planarity
    to_consider = [x for x in to_consider if len(x) > 1]
    n = len(to_consider)

    for i in range(n):
        #Get the positive and negative entries of this part
        ap = [x for x in to_consider[i] if x>0]
        an = [abs(x) for x in to_consider[i] if x<0]

        #Check if a includes numbers in both the top and bottom rows
        if len(ap) > 0 and len(an) > 0:
            for j in range(n):
                if i == j:
                    continue
                #Get the positive and negative entries of this part
                bp = [x for x in to_consider[j] if x>0]
                bn = [abs(x) for x in to_consider[j] if x<0]

                #Skip the ones that don't involve numbers in both
                #the bottom and top rows
                if len(bn) == 0 or len(bp) == 0:
                    continue

                #Make sure that if min(bp) > max(ap)
                #then min(bn) >  max(an)
                if max(bp) > max(ap):
                    if min(bn) < min(an):
                        return False

        #Go through the bottom and top rows
        for row in [ap, an]:
            if len(row) > 1:
                row.sort()
                for s in range(len(row)-1):
                    if row[s] + 1 == row[s+1]:
                        #No gap, continue on
                        continue

                    rng = range(row[s] + 1, row[s+1])

                    #Go through and make sure any parts that
                    #contain numbers in this range are completely
                    #contained in this range
                    for j in range(n):
                        if i == j:
                            continue

                        #Make sure we make the numbers negative again
                        #if we are in the bottom row
                        if row is ap:
                            sr = Set(rng)
                        else:
                            sr = Set(map(lambda x: -1*x, rng))

                        sj = Set(to_consider[j])
                        intersection = sr.intersection(sj)
                        if intersection:
                            if sj != intersection:
                                return False

    return True


def to_graph(sp):
    r"""
    Return a graph representing the set partition ``sp``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: g = da.to_graph( da.to_set_partition([[1,-2],[2,-1]])); g
        Graph on 4 vertices

        sage: g.vertices()
        [-2, -1, 1, 2]
        sage: g.edges()
        [(-2, 1, None), (-1, 2, None)]
    """
    g = Graph()
    for part in sp:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex(part_list[0])
        for i in range(1, len(part_list)):
            g.add_vertex(part_list[i])
            g.add_edge(part_list[i-1], part_list[i])
    return g

def pair_to_graph(sp1, sp2):
    r"""
    Return a graph consisting of the disjoint union of the graphs of set
    partitions ``sp1`` and ``sp2`` along with edges joining the bottom
    row (negative numbers) of ``sp1`` to the top row (positive numbers)
    of ``sp2``.

    The vertices of the graph ``sp1`` appear in the result as pairs
    ``(k, 1)``, whereas the vertices of the graph ``sp2`` appear as
    pairs ``(k, 2)``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,-2],[2,-1]])
        sage: g = da.pair_to_graph( sp1, sp2 ); g
        Graph on 8 vertices

        sage: g.vertices()
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges()
        [((-2, 1), (1, 1), None), ((-2, 1), (2, 2), None),
         ((-2, 2), (1, 2), None), ((-1, 1), (1, 2), None),
         ((-1, 1), (2, 1), None), ((-1, 2), (2, 2), None)]

    Another example which used to be wrong until :trac:`15958`::

        sage: sp3 = da.to_set_partition([[1, -1], [2], [-2]])
        sage: sp4 = da.to_set_partition([[1], [-1], [2], [-2]])
        sage: g = da.pair_to_graph( sp3, sp4 ); g
        Graph on 8 vertices

        sage: g.vertices()
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges()
        [((-2, 1), (2, 2), None), ((-1, 1), (1, 1), None),
         ((-1, 1), (1, 2), None)]
    """
    g = Graph()

    #Add the first set partition to the graph
    for part in sp1:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex( (part_list[0],1) )

            #Add the edge to the second part of the graph
            if part_list[0] < 0:
                g.add_edge( (part_list[0], 1), (abs(part_list[0]), 2)  )

        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 1) )

            #Add the edge to the second part of the graph
            if part_list[i] < 0:
                g.add_edge( (part_list[i], 1), (abs(part_list[i]), 2) )

            #Add the edge between adjacent elements of a part
            g.add_edge( (part_list[i-1], 1), (part_list[i], 1) )

    #Add the second set partition to the graph
    for part in sp2:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex( (part_list[0], 2) )
        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 2) )
            g.add_edge( (part_list[i-1], 2), (part_list[i], 2) )

    return g

def propagating_number(sp):
    r"""
    Return the propagating number of the set partition ``sp``.

    The propagating number is the number of blocks with both a positive and
    negative number.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,2],[-2,-1]])
        sage: da.propagating_number(sp1)
        2
        sage: da.propagating_number(sp2)
        0
    """
    pn = 0
    for part in sp:
        if min(part) < 0  and max(part) > 0:
            pn += 1
    return pn

def to_set_partition(l, k=None):
    r"""
    Convert a list of a list of numbers to a set partitions. Each list
    of numbers in the outer list specifies the numbers contained in one
    of the blocks in the set partition.

    If `k` is specified, then the set partition will be a set partition
    of `\{1, \ldots, k, -1, \ldots, -k\}`. Otherwise, `k` will default to
    the minimum number needed to contain all of the specified numbers.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.to_set_partition([[1,-1],[2,-2]]) == da.identity_set_partition(2)
        True
    """
    if k is None:
        if l == []:
            return SetPartition([])
        else:
            k = max( map( lambda x: max( map(abs, x) ), l) )

    to_be_added = Set( range(1, k+1) + map(lambda x: -1*x, range(1, k+1) ) )

    sp = []
    for part in l:
        spart = Set(part)
        to_be_added -= spart
        sp.append(spart)

    for singleton in to_be_added:
        sp.append(Set([singleton]))

    return SetPartition(sp)

def to_Brauer_partition(l, k=None):
    r"""
    Same as :func:`to_set_partition` but assumes omitted elements are
    connected straight through.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.to_Brauer_partition([[1,2],[-1,-2]]) == SetPartition([[1,2],[-1,-2]])
        True
        sage: da.to_Brauer_partition([[1,3],[-1,-3]]) == SetPartition([[1,3],[-3,-1],[2,-2]])
        True
        sage: da.to_Brauer_partition([[1,2],[-1,-2]], k=4) == SetPartition([[1,2],[-1,-2],[3,-3],[4,-4]])
        True
        sage: da.to_Brauer_partition([[1,-4],[-3,-1],[3,4]]) == SetPartition([[-3,-1],[2,-2],[1,-4],[3,4]])
        True
    """
    L = to_set_partition(l, k=k)
    L2 = []
    paired = []
    not_paired = []
    for i in L:
        L2.append(list(i))
    for i in L2:
        if len(i) >= 3:
            raise ValueError("blocks must have size at most 2, but {0} has {1}".format(i, len(i)))
        if (len(i) == 2):
            paired.append(i)
        if (len(i) == 1):
            not_paired.append(i)
    if any(i[0] in j or -1*i[0] in j for i in not_paired for j in paired):
        raise ValueError("unable to convert {0} to a Brauer partition due to the invalid block {1}".format(l, i))
    for i in not_paired:
        if [-1*i[0]] in not_paired:
            not_paired.remove([-1*i[0]])
        paired.append([i[0], -1*i[0]])
    return to_set_partition(paired)

def identity_set_partition(k):
    """
    Return the identity set partition `\{\{1, -1\}, \ldots, \{k, -k\}\}`

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.identity_set_partition(2)
        {{-2, 2}, {-1, 1}}
    """
    if k in ZZ:
        return SetPartition( [[i,-i] for i in range(1, k + 1)] )
    # Else k in 1/2 ZZ
    return SetPartition( [[i, -i] for i in range(1, k + ZZ(3)/ZZ(2))] )

def set_partition_composition(sp1, sp2):
    r"""
    Return a tuple consisting of the composition of the set partitions
    ``sp1`` and ``sp2`` and the number of components removed from the middle
    rows of the graph.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,-2],[2,-1]])
        sage: da.set_partition_composition(sp1, sp2) == (da.identity_set_partition(2), 0)
        True
    """
    g = pair_to_graph(sp1, sp2)
    connected_components = g.connected_components()

    res = []
    total_removed = 0
    for cc in connected_components:
        #Remove the vertices that live in the middle two rows
        new_cc = [x for x in cc if not ( (x[0]<0 and x[1] == 1) or (x[0]>0 and x[1]==2))]

        if new_cc == []:
            if len(cc) > 1:
                total_removed += 1
        else:
            res.append( Set(map(lambda x: x[0], new_cc)) )

    return (SetPartition(Set(res)), total_removed)

##########################################################################
# END BORROWED CODE
##########################################################################

