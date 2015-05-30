r"""
Finite semilattices and lattices
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.elements import (LatticePosetElement,
                                           MeetSemilatticeElement,
                                           JoinSemilatticeElement)

####################################################################################

def MeetSemilattice(data, *args, **options):
    r"""
    Construct a meet semi-lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a meet semilattice.

    .. seealso:: :func:`Poset`, :func:`JoinSemilattice`, :func:`LatticePoset`

    EXAMPLES:

    Using data that defines a poset::

          sage: MeetSemilattice([[1,2],[3],[3]])
          Finite meet-semilattice containing 4 elements

          sage: MeetSemilattice([[1,2],[3],[3]], cover_relations = True)
          Finite meet-semilattice containing 4 elements

    Using a previously constructed poset::

          sage: P = Poset([[1,2],[3],[3]])
          sage: L = MeetSemilattice(P); L
          Finite meet-semilattice containing 4 elements
          sage: type(L)
          <class 'sage.combinat.posets.lattices.FiniteMeetSemilattice_with_category'>

    If the data is not a lattice, then an error is raised::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: MeetSemilattice((elms, rels))
          Traceback (most recent call last):
          ...
          ValueError: Not a meet semilattice.
    """
    if isinstance(data,FiniteMeetSemilattice) and len(args) == 0 and len(options) == 0:
        return data
    P = Poset(data, *args, **options)
    if not P.is_meet_semilattice():
        raise ValueError("Not a meet semilattice.")
    return FiniteMeetSemilattice(P)

class FiniteMeetSemilattice(FinitePoset):
    """
    .. note::
        We assume that the argument passed to MeetSemilattice is the poset
        of a meet-semilattice (i.e. a poset with greatest lower bound for
        each pair of elements).

    TESTS::

        sage: M = MeetSemilattice([[1,2],[3],[3]])
        sage: TestSuite(M).run()

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: M = MeetSemilattice(P)
        sage: TestSuite(M).run()

    """
    Element = MeetSemilatticeElement

    def _repr_(self):
        r"""
        TESTS::

            sage: M = MeetSemilattice([[1,2],[3],[3]])
            sage: M._repr_()
            'Finite meet-semilattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: M = MeetSemilattice(P)
            sage: M._repr_()
            'Finite meet-semilattice containing 4 elements'
        """
        s  = "Finite meet-semilattice containing %s elements" %self._hasse_diagram.order()
        if self._with_linear_extension:
            s += " with distinguished linear extension"
        return s

    def meet_matrix(self):
        """
        Return a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the meet (greatest lower bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = LatticePoset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: M = P.meet_matrix(); M
            [0 0 0 0 0 0 0 0]
            [0 1 0 1 0 0 0 1]
            [0 0 2 2 2 0 2 2]
            [0 1 2 3 2 0 2 3]
            [0 0 2 2 4 0 2 4]
            [0 0 0 0 0 5 5 5]
            [0 0 2 2 2 5 6 6]
            [0 1 2 3 4 5 6 7]
            sage: M[P(4).vertex,P(3).vertex] == P(0).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(2).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(5).vertex
            False
        """
        return self._hasse_diagram.meet_matrix()

    def meet(self, x, y=None):
        r"""
        Return the meet of given elements in the lattice.

        EXAMPLES::

            sage: D = Posets.DiamondPoset(5)
            sage: D.meet(1, 2)
            0
            sage: D.meet(1, 1)
            1
            sage: D.meet(1, 0)
            0
            sage: D.meet(1, 4)
            1

        Using list of elements as an argument. Meet of empty list is
        the bottom element::

            sage: B4=Posets.BooleanLattice(4)
            sage: B4.meet([3,5,6])
            0
            sage: B4.meet([])
            15

        Test that this method also works for facade lattices::

            sage: L = LatticePoset([[1,2],[3],[3]], facade = True)
            sage: L.meet(2, 3)
            2
            sage: L.meet(1, 2)
            0

        """
        if y is not None: # Handle basic case fast
            i, j = map(self._element_to_vertex, (x,y))
            return self._vertex_to_element(self._hasse_diagram._meet[i,j])
        L = map(self._element_to_vertex, x)
        m = self.cardinality()-1 # m = top element
        for i in L:
            m = self._hasse_diagram._meet[i, m]
        return self._vertex_to_element(m)

####################################################################################

def JoinSemilattice(data, *args, **options):
    r"""
    Construct a join semi-lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a join semilattice.

    .. seealso:: :func:`Poset`, :func:`MeetSemilattice`, :func:`LatticePoset`

    EXAMPLES:

    Using data that defines a poset::

          sage: JoinSemilattice([[1,2],[3],[3]])
          Finite join-semilattice containing 4 elements

          sage: JoinSemilattice([[1,2],[3],[3]], cover_relations = True)
          Finite join-semilattice containing 4 elements

    Using a previously constructed poset::

          sage: P = Poset([[1,2],[3],[3]])
          sage: J = JoinSemilattice(P); J
          Finite join-semilattice containing 4 elements
          sage: type(J)
          <class 'sage.combinat.posets.lattices.FiniteJoinSemilattice_with_category'>

    If the data is not a lattice, then an error is raised::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: JoinSemilattice((elms, rels))
          Traceback (most recent call last):
          ...
          ValueError: Not a join semilattice.
    """
    if isinstance(data,FiniteJoinSemilattice) and len(args) == 0 and len(options) == 0:
        return data
    P = Poset(data, *args, **options)
    if not P.is_join_semilattice():
        raise ValueError("Not a join semilattice.")
    return FiniteJoinSemilattice(P)

class FiniteJoinSemilattice(FinitePoset):
    """
    We assume that the argument passed to FiniteJoinSemilattice is the
    poset of a join-semilattice (i.e. a poset with least upper bound
    for each pair of elements).

    TESTS::

        sage: J = JoinSemilattice([[1,2],[3],[3]])
        sage: TestSuite(J).run()

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: J = JoinSemilattice(P)
        sage: TestSuite(J).run()

    """
    Element = JoinSemilatticeElement

    def _repr_(self):
        r"""
        TESTS::

            sage: J = JoinSemilattice([[1,2],[3],[3]])
            sage: J._repr_()
            'Finite join-semilattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: J = JoinSemilattice(P)
            sage: J._repr_()
            'Finite join-semilattice containing 4 elements'
        """
        s = "Finite join-semilattice containing %s elements"%self._hasse_diagram.order()
        if self._with_linear_extension:
            s += " with distinguished linear extension"
        return s

    def join_matrix(self):
        """
        Return a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the join (least upper bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = LatticePoset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: J = P.join_matrix(); J
            [0 1 2 3 4 5 6 7]
            [1 1 3 3 7 7 7 7]
            [2 3 2 3 4 6 6 7]
            [3 3 3 3 7 7 7 7]
            [4 7 4 7 4 7 7 7]
            [5 7 6 7 7 5 6 7]
            [6 7 6 7 7 6 6 7]
            [7 7 7 7 7 7 7 7]
            sage: J[P(4).vertex,P(3).vertex] == P(7).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(5).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(2).vertex
            False
        """
        return self._hasse_diagram.join_matrix()

    def join(self, x, y=None):
        r"""
        Return the join of given elements in the lattice.

        INPUT:

        -  ``x, y`` - two elements of the (semi)lattice OR

        -  ``x`` - a list or tuple of elements

        EXAMPLES::

            sage: D = Posets.DiamondPoset(5)
            sage: D.join(1, 2)
            4
            sage: D.join(1, 1)
            1
            sage: D.join(1, 4)
            4
            sage: D.join(1, 0)
            1

        Using list of elements as an argument. Join of empty list is
        the bottom element::

            sage: B4=Posets.BooleanLattice(4)
            sage: B4.join([2,4,8])
            14
            sage: B4.join([])
            0

        Test that this method also works for facade lattices::

            sage: L = LatticePoset([[1,2],[3],[3]], facade = True)
            sage: L.join(1, 0)
            1
            sage: L.join(1, 2)
            3

        """
        if y is not None: # Handle basic case fast
            i, j = map(self._element_to_vertex, (x,y))
            return self._vertex_to_element(self._hasse_diagram._join[i,j])
        L = map(self._element_to_vertex, x)
        j = 0 # j = bottom element
        for i in L:
            j = self._hasse_diagram._join[i, j]
        return self._vertex_to_element(j)
        

####################################################################################

def LatticePoset(data, *args, **options):
    r"""
    Construct a lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a lattice.

    OUTPUT:

        FiniteLatticePoset -- an instance of :class:`FiniteLatticePoset`

    .. seealso:: :class:`Posets`, :class:`FiniteLatticePosets`, :func:`JoinSemiLattice`, :func:`MeetSemiLattice`

    EXAMPLES:

    Using data that defines a poset::

        sage: LatticePoset([[1,2],[3],[3]])
        Finite lattice containing 4 elements

        sage: LatticePoset([[1,2],[3],[3]], cover_relations = True)
        Finite lattice containing 4 elements

    Using a previously constructed poset::

        sage: P = Poset([[1,2],[3],[3]])
        sage: L = LatticePoset(P); L
        Finite lattice containing 4 elements
        sage: type(L)
        <class 'sage.combinat.posets.lattices.FiniteLatticePoset_with_category'>

    If the data is not a lattice, then an error is raised::

        sage: elms = [1,2,3,4,5,6,7]
        sage: rels = [[1,2],[3,4],[4,5],[2,5]]
        sage: LatticePoset((elms, rels))
        Traceback (most recent call last):
        ...
        ValueError: Not a lattice.

    Creating a facade lattice::

        sage: L = LatticePoset([[1,2],[3],[3]], facade = True)
        sage: L.category()
        Join of Category of finite lattice posets and Category of finite enumerated sets and Category of facade sets
        sage: parent(L[0])
        Integer Ring
        sage: TestSuite(L).run(skip = ['_test_an_element']) # is_parent_of is not yet implemented

    """
    if isinstance(data,FiniteLatticePoset) and len(args) == 0 and len(options) == 0:
        return data
    P = Poset(data, *args, **options)
    if not P.is_lattice():
        raise ValueError("Not a lattice.")

    return FiniteLatticePoset(P, category = FiniteLatticePosets(), facade = P._is_facade)

class FiniteLatticePoset(FiniteMeetSemilattice, FiniteJoinSemilattice):
    """
    We assume that the argument passed to FiniteLatticePoset is the
    poset of a lattice (i.e. a poset with greatest lower bound and
    least upper bound for each pair of elements).

    TESTS::

        sage: L = LatticePoset([[1,2],[3],[3]])
        sage: TestSuite(L).run()

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: L = LatticePoset(P)
        sage: TestSuite(L).run()

    """
    Element = LatticePosetElement

    def _repr_(self):
        r"""
        TESTS::

            sage: L = LatticePoset([[1,2],[3],[3]])
            sage: L._repr_()
            'Finite lattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: L = LatticePoset(P)
            sage: L._repr_()
            'Finite lattice containing 4 elements'
        """
        s = "Finite lattice containing %s elements"%self._hasse_diagram.order()
        if self._with_linear_extension:
            s += " with distinguished linear extension"
        return s

    def is_distributive(self):
        r"""
        Return ``True`` if the lattice is distributive, and ``False``
        otherwise.

        A lattice `(L, \vee, \wedge)` is distributive if meet
        distributes over join: `x \wedge (y \vee z) = (x \wedge y)
        \vee (x \wedge z)` for every `x,y,z \in L` just like `x \cdot
        (y+z)=x \cdot y + x \cdot z` in normal arithmetic. For duality
        in lattices it follows that then also join distributes over
        meet.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3]})
            sage: L.is_distributive()
            True
            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.is_distributive()
            False
        """
        if self.cardinality() == 0: return True
        return (self.is_graded() and 
         self.rank() == len(self.join_irreducibles()) ==
         len(self.meet_irreducibles()))

    def is_complemented(self):
        r"""
        Returns ``True`` if ``self`` is a complemented lattice, and
        ``False`` otherwise.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.is_complemented()
            True

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.is_complemented()
            False
        """
        return self._hasse_diagram.is_complemented_lattice()

    def complements(self, element=None):
        r"""
        Return the list of complements of an element in the lattice,
        or the dictionary of complements for all elements.

        Elements `x` and `y` are complements if their meet and join
        are respectively the bottom and the top element of the lattice.

        INPUT:

        - ``element`` - an element of the poset whose complement is
          returned. If ``None`` (default) then dictionary of
          complements for all elements having at least one
          complement is returned.

        EXAMPLES::

            sage: L=LatticePoset({0:['a','b','c'], 'a':[1], 'b':[1], 'c':[1]})
            sage: C = L.complements()

        Let us check that `'a'` and `'b'` are complements of each other::

            sage: 'a' in C['b']
            True
            sage: 'b' in C['a']
            True

        Full list of complements::

            sage: L.complements() # random
            {0: [1], 1: [0], 'a': ['b', 'c'], 'b': ['c', 'a'], 'c': ['b', 'a']}

            sage: L=LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.complements() # random
            {0: [4], 4: [0]}
            sage: L.complements(1)
            []

        TESTS::

            sage: L=LatticePoset({0:['a','b','c'], 'a':[1], 'b':[1], 'c':[1]})
            sage: for v,v_complements in L.complements().items():
            ....:     for v_c in v_complements:
            ....:         assert L.meet(v,v_c) == L.bottom()
            ....:         assert L.join(v,v_c) == L.top()
        """
        if element is None:
            jn = self.join_matrix()
            mt = self.meet_matrix()
            n = self.cardinality()
            zero = 0
            one = n-1
            c = [[] for x in range(n)]
            for x in range(n):
                for y in range(x,n):
                    if jn[x][y]==one and mt[x][y]==zero:
                        c[x].append(y)
                        c[y].append(x)

            comps={}
            for i in range(n):
                if len(c[i]) > 0:
                    comps[self._vertex_to_element(i)] = (
                        [self._vertex_to_element(x) for x in c[i]] )
            return comps

        # Looking for complements of one element.
        if not element in self:
            raise ValueError("element (=%s) not in poset"%element)
        return [x for x in self if
         self.meet(x, element)==self.bottom() and
         self.join(x, element)==self.top()]

    def is_atomic(self):
        r"""
        Returns ``True`` if ``self`` is an atomic lattice and ``False`` otherwise.

        A lattice is atomic if every element can be written as a join of atoms.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.is_atomic()
            True

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.is_atomic()
            False

        NOTES:

        See [Sta97]_, Section 3.3 for a discussion of atomic lattices.

        REFERENCES:

        .. [Sta97] Stanley, Richard.
           Enumerative Combinatorics, Vol. 1.
           Cambridge University Press, 1997

        """
        bottom_element = self.bottom()
        for x in self:
            if x == bottom_element:
                continue
            lcovers = self.lower_covers(x)
            if bottom_element in lcovers:
                continue
            if len(lcovers)<=1:
                return False
        return True

    def is_modular(self, L=None):
        r"""
        Return ``True`` if the lattice is modular and ``False`` otherwise.

        Using the parameter ``L``, this can also be used to check that
        some subset of elements are all modular.

        INPUT:

        - ``L`` -- (default: ``None``) a list of elements to check being
          modular, if ``L`` is ``None``, then this checks the entire lattice

        An element `x` in a lattice `L` is *modular* if `x \leq b` implies

        .. MATH::

            x \vee (a \wedge b) = (x \vee a) \wedge b

        for every `a, b \in L`. We say `L` is modular if `x` is modular
        for all `x \in L`. There are other equivalent definitions,
        see :wikipedia:`Modular_lattice`.

        .. SEEALSO::

            :meth:`is_upper_semimodular`, :meth:`is_lower_semimodular`
            and :meth:`is_modular_element`

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_modular()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_modular()
            False

            sage: L = posets.ChainPoset(6)
            sage: L.is_modular()
            True

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_modular()
            False
            sage: [L.is_modular([x]) for x in L]
            [True, True, False, True, True, False, True]

        ALGORITHM:

        Based on pp. 286-287 of Enumerative Combinatorics, Vol 1 [EnumComb1]_.
        """
        if not self.is_ranked():
            return False
        H = self._hasse_diagram
        n = H.order()
        if L is None:
            return all(H._rank[a] + H._rank[b] ==
                       H._rank[H._meet[a, b]] + H._rank[H._join[a, b]]
                       for a in range(n) for b in range(a + 1, n))

        L = [self._element_to_vertex_dict[x] for x in L]
        return all(H._rank[a] + H._rank[b] ==
                   H._rank[H._meet[a, b]] + H._rank[H._join[a, b]]
                   for a in L for b in range(n))

    def is_modular_element(self, x):
        r"""
        Return ``True`` if ``x`` is a modular element and ``False`` otherwise.

        INPUT:

        - ``x`` -- an element of the lattice

        An element `x` in a lattice `L` is *modular* if `x \leq b` implies

        .. MATH::

            x \vee (a \wedge b) = (x \vee a) \wedge b

        for every `a, b \in L`.

        .. SEEALSO::

            :meth:`is_modular` to check modularity for the full lattice or
            some set of elements

        EXAMPLES::

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_modular()
            False
            sage: [L.is_modular_element(x) for x in L]
            [True, True, False, True, True, False, True]
        """
        return self.is_modular([x])

    def is_upper_semimodular(self):
        r"""
        Return ``True`` if the lattice is upper semimodular and
        ``False`` otherwise.

        A lattice is upper semimodular if for any `x` in the poset that is
        covered by `y` and `z`, both `y` and `z` are covered by their join.

        See also :meth:`is_modular` and :meth:`is_lower_semimodular`.

        See :wikipedia:`Semimodular_lattice`

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_upper_semimodular()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_upper_semimodular()
            False

            sage: L = posets.ChainPoset(6)
            sage: L.is_upper_semimodular()
            True

            sage: L = LatticePoset(posets.IntegerPartitions(4))
            sage: L.is_upper_semimodular()
            True

        ALGORITHM:

        Based on pp. 286-287 of Enumerative Combinatorics, Vol 1 [EnumComb1]_.
        """
        if not self.is_ranked():
            return False
        H = self._hasse_diagram
        n = H.order()
        return all(H._rank[a] + H._rank[b] >=
                   H._rank[H._meet[a, b]] + H._rank[H._join[a, b]]
                   for a in range(n) for b in range(a + 1, n))

    def is_lower_semimodular(self):
        r"""
        Return ``True`` if the lattice is lower semimodular and
        ``False`` otherwise.

        A lattice is lower semimodular if for any `x` in the poset that covers
        `y` and `z`, both `y` and `z` cover their meet.

        See also :meth:`is_modular` and :meth:`is_upper_semimodular`.

        See :wikipedia:`Semimodular_lattice`

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_lower_semimodular()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_lower_semimodular()
            False

            sage: L = posets.ChainPoset(6)
            sage: L.is_lower_semimodular()
            True

        ALGORITHM:

        Based on pp. 286-287 of Enumerative Combinatorics, Vol 1 [EnumComb1]_.
        """
        if not self.is_ranked():
            return False
        H = self._hasse_diagram
        n = H.order()
        return all(H._rank[a] + H._rank[b] <=
                   H._rank[H._meet[a,b]] + H._rank[H._join[a,b]]
                   for a in range(n) for b in range(a+1, n))

    def is_supersolvable(self):
        """
        Return ``True`` if ``self`` is a supersolvable lattice and
        ``False`` otherwise.

        A lattice `L` is *supersolvable* if there exists a maximal chain `C`
        such that every `x \in C` is a modular element in `L`.

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_supersolvable()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_supersolvable()
            False

            sage: L = posets.ChainPoset(6)
            sage: L.is_supersolvable()
            True

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_supersolvable()
            True
            sage: L.is_modular()
            False

            sage: L = LatticePoset({0: [1, 2, 3, 4], 1: [5, 6, 7],
            ....:                   2: [5, 8, 9], 3: [6, 8, 10], 4: [7, 9, 10],
            ....:                   5: [11], 6: [11], 7: [11], 8: [11],
            ....:                   9: [11], 10: [11]})
            sage: L.is_supersolvable()
            False
        """
        from sage.misc.cachefunc import cached_function

        if not self.is_ranked():
            return False

        H = self._hasse_diagram
        height = self.height()
        n = H.order()
        cur = H.maximal_elements()[0]
        next = [H.neighbor_in_iterator(cur)]

        @cached_function
        def is_modular_elt(a):
            return all(H._rank[a] + H._rank[b] ==
                       H._rank[H._meet[a, b]] + H._rank[H._join[a, b]]
                       for b in range(n))

        if not is_modular_elt(cur):
            return False
        while len(next) < height:
            try:
                cur = next[-1].next()
            except StopIteration:
                next.pop()
                if not next:
                    return False
                continue
            if is_modular_elt(cur):
                next.append(H.neighbor_in_iterator(cur))
        return True

############################################################################

FiniteMeetSemilattice._dual_class = FiniteJoinSemilattice
FiniteJoinSemilattice._dual_class = FiniteMeetSemilattice
FiniteLatticePoset   ._dual_class = FiniteLatticePoset
