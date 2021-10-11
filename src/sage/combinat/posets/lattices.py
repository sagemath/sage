# -*- coding: utf-8 -*-
r"""
Finite lattices and semilattices

This module implements finite (semi)lattices. It defines:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`LatticePoset` | Construct a lattice.
    :meth:`MeetSemilattice` | Construct a meet semi-lattice.
    :meth:`JoinSemilattice` | Construct a join semi-lattice.
    :class:`FiniteLatticePoset` | A class for finite lattices.
    :class:`FiniteMeetSemilattice` | A class for finite meet semilattices.
    :class:`FiniteJoinSemilattice` | A class for finite join semilattices.

List of (semi)lattice methods
-----------------------------

**Meet and join**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteMeetSemilattice.meet` | Return the meet of given elements.
    :meth:`~FiniteJoinSemilattice.join` | Return the join of given elements.
    :meth:`~FiniteMeetSemilattice.meet_matrix` | Return the matrix of meets of all elements of the meet semi-lattice.
    :meth:`~FiniteJoinSemilattice.join_matrix` | Return the matrix of joins of all elements of the join semi-lattice.

**Properties of the lattice**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteLatticePoset.is_distributive` | Return ``True`` if the lattice is distributive.
    :meth:`~FiniteLatticePoset.is_modular` | Return ``True`` if the lattice is modular.
    :meth:`~FiniteLatticePoset.is_lower_semimodular` | Return ``True`` if all elements with common upper cover have a common lower cover.
    :meth:`~FiniteLatticePoset.is_upper_semimodular` | Return ``True`` if all elements with common lower cover have a common upper cover.
    :meth:`~FiniteLatticePoset.is_semidistributive` | Return ``True`` if the lattice is both join- and meet-semidistributive.
    :meth:`~FiniteLatticePoset.is_join_semidistributive` | Return ``True`` if the lattice is join-semidistributive.
    :meth:`~FiniteLatticePoset.is_meet_semidistributive` | Return ``True`` if the lattice is meet-semidistributive.
    :meth:`~FiniteLatticePoset.is_join_distributive` | Return ``True`` if the lattice is join-distributive.
    :meth:`~FiniteLatticePoset.is_meet_distributive` | Return ``True`` if the lattice is meet-distributive.
    :meth:`~FiniteLatticePoset.is_atomic` | Return ``True`` if every element of the lattice can be written as a join of atoms.
    :meth:`~FiniteLatticePoset.is_coatomic` | Return ``True`` if every element of the lattice can be written as a meet of coatoms.
    :meth:`~FiniteLatticePoset.is_geometric` | Return ``True`` if the lattice is atomic and upper semimodular.
    :meth:`~FiniteLatticePoset.is_extremal` | Return ``True`` if the lattice is extremal.
    :meth:`~FiniteLatticePoset.is_complemented` | Return ``True`` if every element of the lattice has at least one complement.
    :meth:`~FiniteLatticePoset.is_sectionally_complemented` | Return ``True`` if every interval from the bottom is complemented.
    :meth:`~FiniteLatticePoset.is_cosectionally_complemented` | Return ``True`` if every interval to the top is complemented.
    :meth:`~FiniteLatticePoset.is_relatively_complemented` | Return ``True`` if every interval of the lattice is complemented.
    :meth:`~FiniteLatticePoset.is_pseudocomplemented` | Return ``True`` if every element of the lattice has a (meet-)pseudocomplement.
    :meth:`~FiniteLatticePoset.is_join_pseudocomplemented` | Return ``True`` if every element of the lattice has a join-pseudocomplement.
    :meth:`~FiniteLatticePoset.is_orthocomplemented` | Return ``True`` if the lattice has an orthocomplementation.
    :meth:`~FiniteLatticePoset.is_supersolvable` | Return ``True`` if the lattice is supersolvable.
    :meth:`~FiniteLatticePoset.is_planar` | Return ``True`` if the lattice has an upward planar drawing.
    :meth:`~FiniteLatticePoset.is_dismantlable` | Return ``True`` if the lattice is dismantlable.
    :meth:`~FiniteLatticePoset.is_interval_dismantlable` | Return ``True`` if the lattice is interval dismantlable.
    :meth:`~FiniteLatticePoset.is_sublattice_dismantlable` | Return ``True`` if the lattice is sublattice dismantlable.
    :meth:`~FiniteLatticePoset.is_stone` | Return ``True`` if the lattice is a Stone lattice.
    :meth:`~FiniteLatticePoset.is_trim` | Return ``True`` if the lattice is a trim lattice.
    :meth:`~FiniteLatticePoset.is_vertically_decomposable` | Return ``True`` if the lattice is vertically decomposable.
    :meth:`~FiniteLatticePoset.is_simple` | Return ``True`` if the lattice has no nontrivial congruences.
    :meth:`~FiniteLatticePoset.is_isoform` | Return ``True`` if all congruences of the lattice consists of isoform blocks.
    :meth:`~FiniteLatticePoset.is_uniform` | Return ``True`` if all congruences of the lattice consists of equal-sized blocks.
    :meth:`~FiniteLatticePoset.is_regular` | Return ``True`` if all congruences of lattice are determined by any of the congruence blocks.
    :meth:`~FiniteLatticePoset.is_subdirectly_reducible` | Return ``True`` if the lattice is a sublattice of the product of smaller lattices.
    :meth:`~FiniteLatticePoset.is_constructible_by_doublings` | Return ``True`` if the lattice is constructible by doublings from the one-element lattice.
    :meth:`~FiniteLatticePoset.breadth` | Return the breadth of the lattice.

**Specific elements**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteMeetSemilattice.atoms()` | Return elements covering the bottom element.
    :meth:`~FiniteJoinSemilattice.coatoms()` | Return elements covered by the top element.
    :meth:`~FiniteLatticePoset.double_irreducibles` | Return double irreducible elements.
    :meth:`~FiniteLatticePoset.join_primes` | Return the join prime elements.
    :meth:`~FiniteLatticePoset.meet_primes` | Return the meet prime elements.
    :meth:`~FiniteLatticePoset.complements` | Return the list of complements of an element, or the dictionary of complements for all elements.
    :meth:`~FiniteMeetSemilattice.pseudocomplement` | Return the pseudocomplement of an element.
    :meth:`~FiniteLatticePoset.is_modular_element` | Return ``True`` if given element is modular in the lattice.
    :meth:`~FiniteLatticePoset.is_left_modular_element` | Return ``True`` if given element is left modular in the lattice.
    :meth:`~FiniteLatticePoset.neutral_elements` | Return neutral elements of the lattice.
    :meth:`~FiniteLatticePoset.canonical_joinands` | Return the canonical joinands of an element.
    :meth:`~FiniteLatticePoset.canonical_meetands` | Return the canonical meetands of an element.

**Sublattices**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteLatticePoset.sublattice` | Return sublattice generated by list of elements.
    :meth:`~FiniteMeetSemilattice.submeetsemilattice` | Return meet-subsemilattice generated by list of elements.
    :meth:`~FiniteMeetSemilattice.subjoinsemilattice` | Return join-subsemilattice generated by list of elements.
    :meth:`~FiniteLatticePoset.is_sublattice` | Return ``True`` if the lattice is a sublattice of given lattice.
    :meth:`~FiniteLatticePoset.sublattices` | Return all sublattices of the lattice.
    :meth:`~FiniteLatticePoset.sublattices_lattice` | Return the lattice of sublattices.
    :meth:`~FiniteLatticePoset.isomorphic_sublattices_iterator` | Return an iterator over the sublattices isomorphic to given lattice.
    :meth:`~FiniteLatticePoset.maximal_sublattices` | Return maximal sublattices of the lattice.
    :meth:`~FiniteLatticePoset.frattini_sublattice` | Return the intersection of maximal sublattices of the lattice.
    :meth:`~FiniteLatticePoset.skeleton` | Return the skeleton of the lattice.
    :meth:`~FiniteLatticePoset.center` | Return the sublattice of complemented neutral elements.
    :meth:`~FiniteLatticePoset.vertical_decomposition` | Return the vertical decomposition of the lattice.

**Miscellaneous**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteLatticePoset.moebius_algebra` | Return the Möbius algebra of the lattice.
    :meth:`~FiniteLatticePoset.quantum_moebius_algebra` | Return the quantum Möbius algebra of the lattice.
    :meth:`~FiniteLatticePoset.vertical_composition` | Return ordinal sum of lattices with top/bottom element unified.
    :meth:`~FiniteLatticePoset.day_doubling` | Return the lattice with Alan Day's doubling construction of a subset.
    :meth:`~FiniteLatticePoset.adjunct` | Return the adjunct with other lattice.
    :meth:`~FiniteLatticePoset.subdirect_decomposition` | Return the subdirect decomposition of the lattice.
    :meth:`~FiniteLatticePoset.congruence` | Return the congruence generated by lists of elements.
    :meth:`~FiniteLatticePoset.quotient` | Return the quotient lattice by a congruence.
    :meth:`~FiniteLatticePoset.congruences_lattice` | Return the lattice of congruences.
"""
# *****************************************************************************
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.elements import (LatticePosetElement,
                                           MeetSemilatticeElement,
                                           JoinSemilatticeElement)
from sage.combinat.posets.hasse_diagram import LatticeError


####################################################################################

def MeetSemilattice(data=None, *args, **options):
    r"""
    Construct a meet semi-lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a meet semilattice.

    .. SEEALSO:: :func:`Poset`, :func:`JoinSemilattice`, :func:`LatticePoset`

    EXAMPLES:

    Using data that defines a poset::

        sage: MeetSemilattice([[1,2],[3],[3]])
        Finite meet-semilattice containing 3 elements

        sage: MeetSemilattice([[1,2],[3],[3]], cover_relations = True)
        Finite meet-semilattice containing 3 elements

    Using a previously constructed poset::

        sage: P = Poset([[1,2],[3],[3]])
        sage: L = MeetSemilattice(P); L
        Finite meet-semilattice containing 3 elements
        sage: type(L)
        <class 'sage.combinat.posets.lattices.FiniteMeetSemilattice_with_category'>

    If the data is not a lattice, then an error is raised::

        sage: MeetSemilattice({'a': ['b', 'c'], 'b': ['d', 'e'],
        ....:                  'c': ['d', 'e'], 'd': ['f'], 'e': ['f']})
        Traceback (most recent call last):
        ...
        LatticeError: no meet for e and d
    """
    if isinstance(data, FiniteMeetSemilattice) and not args and not options:
        return data
    if "check" in options:
        check = options.pop('check')
    else:
        check = True
    P = Poset(data, *args, **options)
    if check:
        try:
            P._hasse_diagram.meet_matrix()
        except LatticeError as error:
            error.x = P._vertex_to_element(error.x)
            error.y = P._vertex_to_element(error.y)
            raise
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
    _desc = 'Finite meet-semilattice'

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

        INPUT:

        -  ``x, y`` -- two elements of the (semi)lattice OR
        -  ``x`` -- a list or tuple of elements

        EXAMPLES::

            sage: D = posets.DiamondPoset(5)
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

            sage: B4=posets.BooleanLattice(4)
            sage: B4.meet([3,5,6])
            0
            sage: B4.meet([])
            15

        For non-facade lattices operator ``*`` works for meet::

            sage: L = posets.PentagonPoset(facade=False)
            sage: L(1)*L(2)
            0

        .. SEEALSO::

            - Dual function: :meth:`~sage.combinat.posets.lattices.FiniteJoinSemilattice.join`
        """
        mt = self._hasse_diagram.meet_matrix()
        if y is not None:  # Handle basic case fast
            i, j = map(self._element_to_vertex, (x, y))
            return self._vertex_to_element(mt[i, j])
        m = self.cardinality() - 1  # m = top element
        for i in (self._element_to_vertex(_) for _ in x):
            m = mt[i, m]
        return self._vertex_to_element(m)

    def atoms(self):
        """
        Return the list atoms of this (semi)lattice.

        An *atom* of a lattice is an element covering the bottom element.

        EXAMPLES::

            sage: L = posets.DivisorLattice(60)
            sage: sorted(L.atoms())
            [2, 3, 5]

        .. SEEALSO::

            - Dual function: :meth:`~FiniteJoinSemilattice.coatoms`

        TESTS::

            sage: LatticePoset().atoms()
            []
            sage: LatticePoset({0: []}).atoms()
            []
        """
        if self.cardinality() == 0:
            return []
        return self.upper_covers(self.bottom())

    def submeetsemilattice(self, elms):
        r"""
        Return the smallest meet-subsemilattice containing elements on the given list.

        INPUT:

        - ``elms`` -- a list of elements of the lattice.

        EXAMPLES::

            sage: L = posets.DivisorLattice(1000)
            sage: L_ = L.submeetsemilattice([200, 250, 125]); L_
            Finite meet-semilattice containing 5 elements
            sage: L_.list()
            [25, 50, 200, 125, 250]

        .. SEEALSO::

            - Dual function: :meth:`subjoinsemilattice`

        TESTS::

            sage: L.submeetsemilattice([])
            Finite meet-semilattice containing 0 elements
        """
        gens_remaining = set(elms)
        current_set = set()

        # We add elements one by one in 'current_set'.
        #
        # When adding a point g to 'current_set', we add to 'gens_remaining' all
        # meet obtained from g and another point of 'current_set'.
        while gens_remaining:
            g = gens_remaining.pop()
            if g in current_set:
                continue
            for x in current_set:
                gens_remaining.add(self.meet(x, g))
            current_set.add(g)

        return MeetSemilattice(self.subposet(current_set))

    def subjoinsemilattice(self, elms):
        r"""
        Return the smallest join-subsemilattice containing elements on the given list.

        INPUT:

        - ``elms`` -- a list of elements of the lattice.

        EXAMPLES::

            sage: L = posets.DivisorLattice(1000)
            sage: L_ = L.subjoinsemilattice([2, 25, 125]); L_
            Finite join-semilattice containing 5 elements
            sage: sorted(L_.list())
            [2, 25, 50, 125, 250]

        .. SEEALSO::

            - Dual function: :meth:`submeetsemilattice`

        TESTS::

            sage: L.subjoinsemilattice([])
            Finite join-semilattice containing 0 elements
        """
        gens_remaining = set(elms)
        current_set = set()

        # We add elements one by one in 'current_set'.
        #
        # When adding a point g to 'current_set', we add to 'gens_remaining' all
        # join obtained from g and another point of 'current_set'.
        while gens_remaining:
            g = gens_remaining.pop()
            if g in current_set:
                continue
            for x in current_set:
                gens_remaining.add(self.join(x, g))
            current_set.add(g)

        return JoinSemilattice(self.subposet(current_set))

    def pseudocomplement(self, element):
        r"""
        Return the pseudocomplement of ``element``, if it exists.

        The (meet-)pseudocomplement is the greatest element whose
        meet with given element is the bottom element. I.e.
        in a meet-semilattice with bottom element `\hat{0}`
        the pseudocomplement of an element `e` is the element
        `e^\star` such that `e \wedge e^\star = \hat{0}` and
        `e' \le e^\star` if `e \wedge e' = \hat{0}`.

        See :wikipedia:`Pseudocomplement`.

        INPUT:

        - ``element`` -- an element of the lattice.

        OUTPUT:

        An element of the lattice or ``None`` if the pseudocomplement does
        not exist.

        EXAMPLES:

        The pseudocomplement's pseudocomplement is not always the original
        element::

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [5], 4: [6], 5: [6]})
            sage: L.pseudocomplement(2)
            5
            sage: L.pseudocomplement(5)
            4

        An element can have complements but no pseudocomplement, or vice
        versa::

            sage: L = LatticePoset({0: [1, 2], 1: [3, 4, 5], 2: [5], 3: [6],
            ....:                   4: [6], 5: [6]})
            sage: L.complements(1), L.pseudocomplement(1)
            ([], 2)
            sage: L.complements(2), L.pseudocomplement(2)
            ([3, 4], None)

        .. SEEALSO:: :meth:`~sage.combinat.posets.lattices.FiniteLatticePoset.is_pseudocomplemented`

        TESTS::

            sage: L = LatticePoset({'a': []})
            sage: L.pseudocomplement('a')
            'a'
            sage: L = LatticePoset({'a': ['b'], 'b': ['c']})
            sage: [L.pseudocomplement(e) for e in ['a', 'b', 'c']]
            ['c', 'a', 'a']
        """
        v = self._element_to_vertex(element)
        e = self._hasse_diagram.pseudocomplement(v)
        if e is None:
            return None
        return self._vertex_to_element(e)

####################################################################################


def JoinSemilattice(data=None, *args, **options):
    r"""
    Construct a join semi-lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a join semilattice

    .. SEEALSO:: :func:`Poset`, :func:`MeetSemilattice`, :func:`LatticePoset`

    EXAMPLES:

    Using data that defines a poset::

        sage: JoinSemilattice([[1,2],[3],[3]])
        Finite join-semilattice containing 3 elements

        sage: JoinSemilattice([[1,2],[3],[3]], cover_relations = True)
        Finite join-semilattice containing 3 elements

    Using a previously constructed poset::

        sage: P = Poset([[1,2],[3],[3]])
        sage: J = JoinSemilattice(P); J
        Finite join-semilattice containing 3 elements
        sage: type(J)
        <class 'sage.combinat.posets.lattices.FiniteJoinSemilattice_with_category'>

    If the data is not a lattice, then an error is raised::

        sage: JoinSemilattice({'a': ['b', 'c'], 'b': ['d', 'e'],
        ....:                  'c': ['d', 'e'], 'd': ['f'], 'e': ['f']})
        Traceback (most recent call last):
        ...
        LatticeError: no join for b and c
    """
    if isinstance(data, FiniteJoinSemilattice) and not args and not options:
        return data
    if "check" in options:
        check = options.pop('check')
    else:
        check = True
    P = Poset(data, *args, **options)
    if check:
        try:
            P._hasse_diagram.join_matrix()
        except LatticeError as error:
            error.x = P._vertex_to_element(error.x)
            error.y = P._vertex_to_element(error.y)
            raise
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
    _desc = 'Finite join-semilattice'

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

        -  ``x, y`` -- two elements of the (semi)lattice OR
        -  ``x`` -- a list or tuple of elements

        EXAMPLES::

            sage: D = posets.DiamondPoset(5)
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

            sage: B4=posets.BooleanLattice(4)
            sage: B4.join([2,4,8])
            14
            sage: B4.join([])
            0

        For non-facade lattices operator ``+`` works for join::

            sage: L = posets.PentagonPoset(facade=False)
            sage: L(1)+L(2)
            4

        .. SEEALSO::

            - Dual function: :meth:`~sage.combinat.posets.lattices.FiniteMeetSemilattice.meet`
        """
        jn = self._hasse_diagram.join_matrix()
        if y is not None:  # Handle basic case fast
            i, j = map(self._element_to_vertex, (x, y))
            return self._vertex_to_element(jn[i, j])
        j = 0  # j = bottom element
        for i in (self._element_to_vertex(_) for _ in x):
            j = jn[i, j]
        return self._vertex_to_element(j)

    def coatoms(self):
        """
        Return the list of co-atoms of this (semi)lattice.

        A *co-atom* of a lattice is an element covered by the top element.

        EXAMPLES::

            sage: L = posets.DivisorLattice(60)
            sage: sorted(L.coatoms())
            [12, 20, 30]

        .. SEEALSO::

            - Dual function: :meth:`~FiniteMeetSemilattice.atoms`

        TESTS::

            sage: LatticePoset().coatoms()
            []
            sage: LatticePoset({0: []}).coatoms()
            []
        """
        if self.cardinality() == 0:
            return []
        return self.lower_covers(self.top())

###############################################################################


def LatticePoset(data=None, *args, **options):
    r"""
    Construct a lattice from various forms of input data.

    INPUT:

    - ``data``, ``*args``, ``**options`` -- data and options that will
      be passed down to :func:`Poset` to construct a poset that is
      also a lattice.

    OUTPUT:

    An instance of :class:`FiniteLatticePoset`.

    .. SEEALSO::

        :class:`Posets`, :class:`FiniteLatticePosets`,
        :func:`JoinSemiLattice`, :func:`MeetSemiLattice`

    EXAMPLES:

    Using data that defines a poset::

        sage: LatticePoset([[1,2],[3],[3]])
        Finite lattice containing 3 elements

        sage: LatticePoset([[1,2],[3],[3]], cover_relations = True)
        Finite lattice containing 3 elements

    Using a previously constructed poset::

        sage: P = Poset([[1,2],[3],[3]])
        sage: L = LatticePoset(P); L
        Finite lattice containing 3 elements
        sage: type(L)
        <class 'sage.combinat.posets.lattices.FiniteLatticePoset_with_category'>

    If the data is not a lattice, then an error is raised::

        sage: elms = [1,2,3,4,5,6,7]
        sage: rels = [[1,2],[3,4],[4,5],[2,5]]
        sage: LatticePoset((elms, rels))
        Traceback (most recent call last):
        ...
        ValueError: not a meet-semilattice: no bottom element

    Creating a facade lattice::

        sage: L = LatticePoset([[1,2],[3],[3]], facade = True)
        sage: L.category()
        Category of facade finite enumerated lattice posets
        sage: parent(L[0])
        Integer Ring
        sage: TestSuite(L).run(skip = ['_test_an_element']) # is_parent_of is not yet implemented
    """
    if isinstance(data, FiniteLatticePoset) and not args and not options:
        return data
    if "check" in options:
        check = options.pop('check')
    else:
        check = True
    P = Poset(data, *args, **options)
    if P.cardinality() != 0:
        if not P.has_bottom():
            raise ValueError("not a meet-semilattice: no bottom element")
        if check:
            try:
                P._hasse_diagram.join_matrix()
            except LatticeError as error:
                error.x = P._vertex_to_element(error.x)
                error.y = P._vertex_to_element(error.y)
                raise
    return FiniteLatticePoset(P, category=FiniteLatticePosets(), facade=P._is_facade)


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
            'Finite lattice containing 3 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: L = LatticePoset(P)
            sage: L._repr_()
            'Finite lattice containing 3 elements'
        """
        s = "Finite lattice containing %s elements" % self._hasse_diagram.order()
        if self._with_linear_extension:
            s += " with distinguished linear extension"
        return s

    def double_irreducibles(self):
        """
        Return the list of double irreducible elements of this lattice.

        A *double irreducible* element of a lattice is an element
        covering and covered by exactly one element. In other words
        it is neither a meet nor a join of any elements.

        EXAMPLES::

            sage: L = posets.DivisorLattice(12)
            sage: sorted(L.double_irreducibles())
            [3, 4]

            sage: L = posets.BooleanLattice(3)
            sage: L.double_irreducibles()
            []

        .. SEEALSO::

            :meth:`~sage.categories.finite_lattice_posets.FiniteLatticePosets.ParentMethods.meet_irreducibles`,
            :meth:`~sage.categories.finite_lattice_posets.FiniteLatticePosets.ParentMethods.join_irreducibles`

        TESTS::

            sage: LatticePoset().double_irreducibles()
            []
            sage: posets.ChainPoset(2).double_irreducibles()
            []
        """
        H = self._hasse_diagram
        return [self._vertex_to_element(e) for e in H
                if H.in_degree(e) == 1 and H.out_degree(e) == 1]

    def join_primes(self):
        r"""
        Return the join-prime elements of the lattice.

        An element `x` of a lattice `L` is *join-prime* if `x \le a \vee b`
        implies `x \le a` or `x \le b` for every `a, b \in L`.

        These are also called *coprime* in some books. Every join-prime
        is join-irreducible; converse holds if and only if the lattice
        is distributive.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5, 6], 3: [5],
            ....:                   4: [6], 5: [7], 6: [7]})
            sage: L.join_primes()
            [3, 4]

            sage: D12 = posets.DivisorLattice(12)  # Distributive lattice
            sage: D12.join_irreducibles() == D12.join_primes()
            True

        .. SEEALSO::

            - Dual function: :meth:`meet_primes`
            - Other: :meth:`~sage.categories.finite_lattice_posets.FiniteLatticePosets.ParentMethods.join_irreducibles`

        TESTS::

            sage: LatticePoset().join_primes()
            []
            sage: posets.DiamondPoset(5).join_primes()
            []
        """
        return [self._vertex_to_element(v) for
                v in self._hasse_diagram.prime_elements()[0]]

    def meet_primes(self):
        r"""
        Return the meet-prime elements of the lattice.

        An element `x` of a lattice `L` is *meet-prime* if `x \ge a \wedge b`
        implies `x \ge a` or `x \ge b` for every `a, b \in L`.

        These are also called just *prime* in some books. Every meet-prime
        is meet-irreducible; converse holds if and only if the lattice
        is distributive.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5, 6], 3: [5],
            ....:                   4: [6], 5: [7], 6: [7]})
            sage: L.meet_primes()
            [6, 5]

            sage: D12 = posets.DivisorLattice(12)
            sage: sorted(D12.meet_primes())
            [3, 4, 6]

        .. SEEALSO::

            - Dual function: :meth:`join_primes`
            - Other: :meth:`~sage.categories.finite_lattice_posets.FiniteLatticePosets.ParentMethods.meet_irreducibles`

        TESTS::

            sage: LatticePoset().meet_primes()
            []
            sage: posets.DiamondPoset(5).meet_primes()
            []
        """
        return [self._vertex_to_element(v) for
                v in self._hasse_diagram.prime_elements()[1]]

    def neutral_elements(self):
        r"""
        Return the list of neutral elements of the lattice.

        An element `e` of the lattice `L` is *neutral* if the sublattice
        generated by `e`, `x` and `y` is distributive for all `x, y \in L`.
        It can also be characterized as an element of intersection of
        maximal distributive sublattices.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3], 2: [6], 3: [4, 5, 6], 4: [8],
            ....:                   5: [7], 6: [7], 7: [8, 9], 8: [10], 9: [10]})
            sage: L.neutral_elements()
            [1, 3, 8, 10]

        TESTS::

            sage: all(posets.ChainPoset(i).neutral_elements() == list(range(i))
            ....:     for i in range(4))
            True

            sage: posets.BooleanLattice(3).neutral_elements()
            [0, 1, 2, 3, 4, 5, 6, 7]

            sage: L = LatticePoset(DiGraph('QQG?LA??__?OG@C??p???O??A?E??@??@g??Q??S??@??E??@??@???'))
            sage: L.neutral_elements()
            [0, 1, 4, 5, 15, 17]
        """
        t = sorted(self._hasse_diagram.neutral_elements())
        return [self._vertex_to_element(v) for v in t]

    def is_join_distributive(self, certificate=False):
        """
        Return ``True`` if the lattice is join-distributive and ``False``
        otherwise.

        A lattice is *join-distributive* if every interval from an element
        to the join of the element's upper covers is a distributive lattice.
        Actually this distributive sublattice is then a Boolean lattice.

        They are also called as *Dilworth's lattices* and *upper locally
        distributive lattices*. They can be characterized in many other
        ways, see [Dil1940]_.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where `e` is an element such that the interval
          from `e` to the meet of upper covers of `e` is not distributive.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5, 6], 3: [5, 7],
            ....:                   4: [6, 7], 5: [8, 9], 6: [9], 7: [9, 10],
            ....:                   8: [11], 9: [11], 10: [11]})
            sage: L.is_join_distributive()
            True

            sage: L = LatticePoset({1: [2], 2: [3, 4], 3: [5], 4: [6],
            ....:                   5: [7], 6: [7]})
            sage: L.is_join_distributive()
            False
            sage: L.is_join_distributive(certificate=True)
            (False, 2)

        .. SEEALSO::

            - Dual property: :meth:`is_meet_distributive`
            - Weaker properties: :meth:`is_meet_semidistributive`,
              :meth:`is_upper_semimodular`
            - Stronger properties: :meth:`is_distributive`

        TESTS::

            sage: E = LatticePoset()
            sage: E.is_join_distributive()
            True
            sage: E.is_join_distributive(certificate=True)
            (True, None)

            sage: L = LatticePoset({1: []})
            sage: L.is_join_distributive()
            True
            sage: L.is_join_distributive(certificate=True)
            (True, None)

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3: [5, 6],
            ....:                   4: [6], 5: [7], 6:[7]})
            sage: L.is_join_distributive()
            False
            sage: L.is_join_distributive(certificate=True)
            (False, 1)

            sage: L = LatticePoset({1: [2], 2: [3, 4, 5], 3: [6], 4: [6], 5: [6]})
            sage: L.is_join_distributive(certificate=True)
            (False, 2)
        """
        if ((self.is_ranked() and len(self.meet_irreducibles()) == self.rank()) or
            self.cardinality() == 0):
            return (True, None) if certificate else True
        if not certificate:
            return False

        # A lattice that is not join-distributive is either not upper
        # semimodular or contains a diamond as a covering sublattice.
        result = self.is_upper_semimodular(certificate=True)
        if not result[0]:
            return (False, self.meet(result[1]))

        from sage.graphs.digraph import DiGraph
        M3 = DiGraph({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
        diamond = next(self._hasse_diagram.subgraph_search_iterator(M3))
        return (False, self[diamond[0]])

    def is_meet_distributive(self, certificate=False):
        """
        Return ``True`` if the lattice is meet-distributive and ``False``
        otherwise.

        A lattice is *meet-distributive* if every interval to an element
        from the meet of the element's lower covers is a distributive lattice.
        Actually this distributive sublattice is then a Boolean lattice.

        They are also called as *lower locally distributive lattices*.
        They can be characterized in many other ways, see [Dil1940]_.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where `e` is an element such that the interval
          to `e` from the meet of lower covers of `e` is not distributive.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3: [5, 6, 7],
            ....:                   4: [7], 5: [9, 8], 6: [10, 8], 7:
            ....:                   [9, 10], 8: [11], 9: [11], 10: [11]})
            sage: L.is_meet_distributive()
            True

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [5], 4: [6],
            ....:                   5: [6], 6: [7]})
            sage: L.is_meet_distributive()
            False
            sage: L.is_meet_distributive(certificate=True)
            (False, 6)

        .. SEEALSO::

            - Dual property: :meth:`is_join_distributive`
            - Weaker properties: :meth:`is_join_semidistributive`,
              :meth:`is_lower_semimodular`
            - Stronger properties: :meth:`is_distributive`

        TESTS::

            sage: E = LatticePoset()
            sage: E.is_meet_distributive()
            True
            sage: E.is_meet_distributive(certificate=True)
            (True, None)

            sage: L = LatticePoset({1: []})
            sage: L.is_meet_distributive()
            True
            sage: L.is_meet_distributive(certificate=True)
            (True, None)

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [5, 6], 4: [7],
            ....:                   5: [7], 6: [7]})
            sage: L.is_meet_distributive()
            False
            sage: L.is_meet_distributive(certificate=True)
            (False, 7)

            sage: L = LatticePoset({1: [2], 2: [3, 4, 5], 3: [6], 4: [6], 5: [6]})
            sage: L.is_meet_distributive(certificate=True)
            (False, 6)
        """
        if ((self.is_ranked() and len(self.join_irreducibles()) == self.rank()) or
            self.cardinality() == 0):
            return (True, None) if certificate else True
        if not certificate:
            return False

        # A lattice that is not meet-distributive is either not lower
        # semimodular or contains a diamond as a covering sublattice.
        result = self.is_lower_semimodular(certificate=True)
        if not result[0]:
            return (False, self.join(result[1]))

        from sage.graphs.digraph import DiGraph
        M3 = DiGraph({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
        diamond = next(self._hasse_diagram.subgraph_search_iterator(M3))
        return (False, self[diamond[4]])

    def is_stone(self, certificate=False):
        r"""
        Return ``True`` if the lattice is a Stone lattice, and ``False``
        otherwise.

        The lattice is expected to be distributive (and hence
        pseudocomplemented).

        A pseudocomplemented lattice is a Stone lattice if

        .. MATH::

            e^* \vee e^{**} = \top

        for every element `e` of the lattice, where `^*` is the
        pseudocomplement and `\top` is the top element of the lattice.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)`` such that `e^* \vee e^{**} \neq \top`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES:

        Divisor lattices are canonical example::

            sage: D72 = posets.DivisorLattice(72)
            sage: D72.is_stone()
            True

        A non-example::

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [4], 4: [5]})
            sage: L.is_stone()
            False

        .. SEEALSO::

            - Weaker properties: :meth:`is_distributive`

        TESTS::

            sage: LatticePoset().is_stone()  # Empty lattice
            True

            sage: L = LatticePoset(DiGraph('GW?_W@?W@?O?'))
            sage: L.is_stone()  # Pass the fast check, but not a Stone lattice
            False
        """
        # TODO: For now we can factor only undirected graphs. When that
        # is extended to directed, use that; see comment below.

        if not self.is_distributive():
            raise ValueError("the lattice is not distributive")

        from sage.arith.misc import factor
        ok = (True, None) if certificate else True

        # Needed for the empty lattice that has no bottom element.
        if self.cardinality() < 5:
            return ok

        # Quick check:
        # A Stone lattice is direct product of distributive lattices with
        # one atom. Return False if for example the lattice has two atoms
        # and odd number of elements.
        atoms_n = self._hasse_diagram.out_degree(0)
        if atoms_n == 1:
            return ok
        if not certificate:
            if sum([x[1] for x in factor(self.cardinality())]) < atoms_n:
                return False
            if self._hasse_diagram.in_degree(self.cardinality() - 1) < atoms_n:
                return False

        # Quick check failed
        one = self.top()
        tested = set()
        for e in self:
            e_ = self.pseudocomplement(e)
            if e_ not in tested:
                if self.join(e_, self.pseudocomplement(e_)) != one:
                    if certificate:
                        return (False, e)
                    return False
                tested.add(e_)

        return ok

    def is_distributive(self, certificate=False):
        r"""
        Return ``True`` if the lattice is distributive, and ``False``
        otherwise.

        A lattice `(L, \vee, \wedge)` is distributive if meet
        distributes over join: `x \wedge (y \vee z) = (x \wedge y)
        \vee (x \wedge z)` for every `x,y,z \in L` just like `x \cdot
        (y+z)=x \cdot y + x \cdot z` in normal arithmetic. For duality
        in lattices it follows that then also join distributes over
        meet.

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (x, y, z))``, where `x`, `y` and `z` are elements
          of the lattice such that
          `x \wedge (y \vee z) \neq (x \wedge y) \vee (x \wedge z)`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [4], 4: [5]})
            sage: L.is_distributive()
            True
            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3: [6], 4: [6], 5: [6]})
            sage: L.is_distributive()
            False
            sage: L.is_distributive(certificate=True)
            (False, (5, 3, 2))

        .. SEEALSO::

            - Weaker properties: :meth:`is_modular`,
              :meth:`is_semidistributive`, :meth:`is_join_distributive`,
              :meth:`is_meet_distributive`, :meth:`is_subdirectly_reducible`,
              :meth:`is_trim`,
              :meth:`is_constructible_by_doublings` (by interval doubling),
              :meth:`is_extremal`

            - Stronger properties: :meth:`is_stone`

        TESTS::

            sage: [posets.ChainPoset(i).is_distributive() for i in range(3)]
            [True, True, True]
        """
        from sage.graphs.digraph import DiGraph

        ok = (True, None) if certificate else True

        if self.cardinality() == 0:
            return ok

        if (self.is_graded() and
         self.rank() == len(self.join_irreducibles()) ==
         len(self.meet_irreducibles())):
            return ok

        if not certificate:
            return False

        result, cert = self.is_modular(certificate=True)
        if not result:
            return (False, (cert[2], cert[1], cert[0]))
        M3 = DiGraph({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
        diamond = next(self._hasse_diagram.subgraph_search_iterator(M3))
        return (False, (self._vertex_to_element(diamond[1]),
                        self._vertex_to_element(diamond[2]),
                        self._vertex_to_element(diamond[3])))

    def is_semidistributive(self):
        """
        Return ``True`` if the lattice is both join- and meet-semidistributive,
        and ``False`` otherwise.

        EXAMPLES:

        Tamari lattices are typical examples of semidistributive but not
        distributive (and hence not modular) lattices::

            sage: T4 = posets.TamariLattice(4)
            sage: T4.is_semidistributive(), T4.is_distributive()
            (True, False)

        Smallest non-selfdual example::

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [5], 4: [6], 5: [7], 6: [7]})
            sage: L.is_semidistributive()
            True

        The diamond is not semidistributive::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_semidistributive()
            False

        .. SEEALSO::

            - Weaker properties: :meth:`is_join_semidistributive`,
              :meth:`is_meet_semidistributive`
            - Stronger properties: :meth:`is_distributive`

        TESTS::

            sage: LatticePoset().is_semidistributive()
            True
            sage: LatticePoset({1: []}).is_semidistributive()
            True
        """
        H = self._hasse_diagram
        # See trac #21528 for explanation.
        return ( (H.in_degree_sequence().count(1) ==
                 H.out_degree_sequence().count(1)) and
                 self.is_meet_semidistributive() )

    def is_meet_semidistributive(self, certificate=False):
        r"""
        Return ``True`` if the lattice is meet-semidistributive, and ``False``
        otherwise.

        A lattice is meet-semidistributive if for all elements
        `e, x, y` in the lattice we have

        .. MATH::

            e \wedge x = e \wedge y \implies
            e \wedge x = e \wedge (x \vee y)

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (e, x, y))`` such that `e \wedge x = e \wedge y`
          but `e \wedge x \neq e \wedge (x \vee y)`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1:[2, 3, 4], 2:[4, 5], 3:[5, 6],
            ....:                   4:[7], 5:[7], 6:[7]})
            sage: L.is_meet_semidistributive()
            True
            sage: L_ = L.dual()
            sage: L_.is_meet_semidistributive()
            False
            sage: L_.is_meet_semidistributive(certificate=True)
            (False, (5, 4, 6))

        .. SEEALSO::

            - Dual property: :meth:`is_join_semidistributive`
            - Weaker properties: :meth:`is_pseudocomplemented`,
              :meth:`is_interval_dismantlable`
            - Stronger properties: :meth:`is_semidistributive`,
              :meth:`is_join_distributive`,
              :meth:`is_constructible_by_doublings` (by upper pseudo-intervals)

        TESTS::

            sage: LatticePoset().is_meet_semidistributive()
            True

        Smallest lattice that fails the quick check::

            sage: L = LatticePoset(DiGraph('IY_T@A?CC_@?W?O@??'))
            sage: L.is_meet_semidistributive()
            False

        Confirm that :trac:`21340` is fixed::

            sage: posets.BooleanLattice(4).is_meet_semidistributive()
            True
        """
        # See http://www.math.hawaii.edu/~ralph/Preprints/algorithms-survey.pdf
        # for explanation of this
        n = self.cardinality()
        if n == 0:
            if certificate:
                return (True, None)
            return True
        H = self._hasse_diagram
        if not certificate and H.size() * 2 > n * _log_2(n):
            return False

        for v in H:
            if H.in_degree(v) == 1 and H.kappa(v) is None:
                if not certificate:
                    return False
                v_ = next(H.neighbor_in_iterator(v))
                t1 = set(H.depth_first_search(v_))
                t2 = set(H.depth_first_search(v))
                tmp = sorted(t1.difference(t2), reverse=True)
                x = tmp[0]
                for y in tmp:
                    if H.are_incomparable(x, y):
                        return (False,
                                (self._vertex_to_element(v),
                                 self._vertex_to_element(x),
                                 self._vertex_to_element(y)))
        if certificate:
            return (True, None)
        return True

    def is_join_semidistributive(self, certificate=False):
        r"""
        Return ``True`` if the lattice is join-semidistributive, and ``False``
        otherwise.

        A lattice is join-semidistributive if for all elements `e, x, y` in
        the lattice we have

        .. MATH::

            e \vee x = e \vee y \implies
            e \vee x = e \vee (x \wedge y)

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (e, x, y))`` such that `e \vee x = e \vee y`
          but `e \vee x \neq e \vee (x \wedge y)`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: T4 = posets.TamariLattice(4)
            sage: T4.is_join_semidistributive()
            True
            sage: L = LatticePoset({1:[2, 3], 2:[4, 5], 3:[5, 6],
            ....:                   4:[7], 5:[7], 6:[7]})
            sage: L.is_join_semidistributive()
            False
            sage: L.is_join_semidistributive(certificate=True)
            (False, (5, 4, 6))

        .. SEEALSO::

            - Dual property: :meth:`is_meet_semidistributive`
            - Weaker properties: :meth:`is_join_pseudocomplemented`,
              :meth:`is_interval_dismantlable`
            - Stronger properties: :meth:`is_semidistributive`,
              :meth:`is_meet_distributive`,
              :meth:`is_constructible_by_doublings` (by lower pseudo-intervals)

        TESTS::

            sage: LatticePoset().is_join_semidistributive()
            True

        Smallest lattice that fails the quick check::

            sage: L = LatticePoset(DiGraph('IY_T@A?CC_@?W?O@??'))
            sage: L.is_join_semidistributive()
            False

        Confirm that :trac:`21340` is fixed::

            sage: posets.BooleanLattice(3).is_join_semidistributive()
            True
        """
        # See http://www.math.hawaii.edu/~ralph/Preprints/algorithms-survey.pdf
        # for explanation of this
        n = self.cardinality()
        if n == 0:
            if certificate:
                return (True, None)
            return True
        H = self._hasse_diagram
        if not certificate and H.size() * 2 > n * _log_2(n):
            return False

        for v in H:
            if H.out_degree(v) == 1 and H.kappa_dual(v) is None:
                if not certificate:
                    return False
                v_ = next(H.neighbor_out_iterator(v))
                it = H.neighbor_in_iterator
                t1 = set(H.depth_first_search(v_, neighbors=it))
                t2 = set(H.depth_first_search(v, neighbors=it))
                tmp = sorted(t1.difference(t2))
                x = tmp[0]
                for y in tmp:
                    if H.are_incomparable(x, y):
                        return (False,
                                (self._vertex_to_element(v),
                                 self._vertex_to_element(x),
                                 self._vertex_to_element(y)))
        if certificate:
            return (True, None)
        return True

        return all(H.kappa_dual(v) is not None
                   for v in H if H.out_degree(v) == 1)

    def is_extremal(self):
        """
        Return ``True`` if the lattice is extremal, and ``False``
        otherwise.

        A lattice is *extremal* if the number of join-irreducibles is equal
        to the number of meet-irreducibles and to the number of
        cover relations in the longest chains.

        EXAMPLES::

            sage: posets.PentagonPoset().is_extremal()
            True

            sage: P = LatticePoset(posets.SymmetricGroupWeakOrderPoset(3))
            sage: P.is_extremal()
            False

        .. SEEALSO::

            - Stronger properties: :meth:`is_distributive`, :meth:`is_trim`

        REFERENCES:

        - [Mark1992]_
        """
        ji = len(self.join_irreducibles())
        mi = len(self.meet_irreducibles())
        return ji == mi == self.height() - 1

    def is_trim(self, certificate=False):
        """
        Return whether a lattice is trim.

        A lattice is trim if it is extremal and left modular.

        This notion is defined in [Thom2006]_.

        INPUT:

        - certificate -- boolean (default ``False``) whether to return
          instead a maximum chain of left modular elements

        EXAMPLES::

            sage: P = posets.PentagonPoset()
            sage: P.is_trim()
            True

            sage: Q = LatticePoset(posets.SymmetricGroupWeakOrderPoset(3))
            sage: Q.is_trim()
            False

        TESTS::

            sage: LatticePoset({1:[]}).is_trim(True)
            (True, [1])

        .. SEEALSO::

            - Weaker properties: :meth:`is_extremal`
            - Stronger properties: :meth:`is_distributive`

        REFERENCES:

        .. [Thom2006] Hugh Thomas, *An analogue of distributivity for
           ungraded lattices*. Order 23 (2006), no. 2-3, 249-269.
        """
        ji = len(self.join_irreducibles())
        mi = len(self.meet_irreducibles())
        h, chain = self.height(certificate=True)
        if not (ji == mi == h - 1):
            return (False, None) if certificate else False

        if all(self.is_left_modular_element(e) for e in chain):
            return (True, chain) if certificate else True
        else:
            return (False, None) if certificate else False

    def is_complemented(self, certificate=False):
        r"""
        Return ``True`` if the lattice is complemented, and
        ``False`` otherwise.

        A lattice is complemented if every element has at least one
        complement.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where ``e`` is an element without a complement.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: L.is_complemented()
            True

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5, 6], 3: [5], 4: [6],
            ....:                   5: [7], 6: [7]})
            sage: L.is_complemented()
            False
            sage: L.is_complemented(certificate=True)
            (False, 2)

        .. SEEALSO::

            - Stronger properties: :meth:`is_sectionally_complemented`,
              :meth:`is_cosectionally_complemented`,
              :meth:`is_orthocomplemented`
            - Other: :meth:`complements`

        TESTS::

            sage: [posets.ChainPoset(i).is_complemented() for i in range(5)]
            [True, True, True, False, False]
        """
        e = self._hasse_diagram.is_complemented()
        if not certificate:
            return e is None
        if e is None:
            return (True, None)
        return (False, self._vertex_to_element(e))

    def is_cosectionally_complemented(self, certificate=False):
        """
        Return ``True`` if the lattice is cosectionally complemented, and
        ``False`` otherwise.

        A lattice is *cosectionally complemented* if all intervals to
        the top element interpreted as sublattices are complemented
        lattices.

        INPUT:

        - ``certificate`` -- (default: ``False``) Whether to return
          a certificate if the lattice is not cosectionally complemented.

        OUTPUT:

        - If ``certificate=False`` return ``True`` or ``False``.
          If ``certificate=True`` return either ``(True, None)``
          or ``(False, (b, e))``, where `b` is an element so that in the
          sublattice from `b` to the top element has no complement
          for element `e`.

        EXAMPLES:

        The smallest sectionally but not cosectionally complemented lattice::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3: [5], 4: [6], 5: [6]})
            sage: L.is_sectionally_complemented(), L.is_cosectionally_complemented()
            (True, False)

        A sectionally and cosectionally but not relatively complemented
        lattice::

            sage: L = LatticePoset(DiGraph('MYi@O?P??D?OG?@?O_?C?Q??O?W?@??O??'))
            sage: L.is_sectionally_complemented() and L.is_cosectionally_complemented()
            True
            sage: L.is_relatively_complemented()
            False

        Getting a certificate::

            sage: L = LatticePoset(DiGraph('HW?@D?Q?GE?G@??'))
            sage: L.is_cosectionally_complemented(certificate=True)
            (False, (2, 7))

        .. SEEALSO::

            - Dual property: :meth:`is_sectionally_complemented`
            - Weaker properties: :meth:`is_complemented`, :meth:`is_coatomic`,
              :meth:`is_regular`
            - Stronger properties: :meth:`is_relatively_complemented`

        TESTS::

            sage: [posets.ChainPoset(i).is_cosectionally_complemented() for i in range(5)]
            [True, True, True, False, False]
        """
        # Quick check: every cosectionally complemented lattice is coatomic.
        if not certificate and not self.is_coatomic():
            return False

        H = self._hasse_diagram
        jn = H.join_matrix()
        n = H.order()
        for e in range(n-2, -1, -1):
            t = 0
            for uc in H.neighbors_out(e):
                t = jn[t, uc]
                if t == n-1:
                    break
            else:
                if certificate:
                    return (False, (self[e], self[t]))
                return False
        return (True, None) if certificate else True

    def is_relatively_complemented(self, certificate=False):
        """
        Return ``True`` if the lattice is relatively complemented, and
        ``False`` otherwise.

        A lattice is relatively complemented if every interval of it
        is a complemented lattice.

        INPUT:

        - ``certificate`` -- (default: ``False``) Whether to return
          a certificate if the lattice is not relatively complemented.

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (a, b, c))``, where `b` is the only element that
          covers `a` and is covered by `c`. If ``certificate=False``
          return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4, 8], 2: [5, 6], 3: [5, 7],
            ....:                   4: [6, 7], 5: [9], 6: [9], 7: [9], 8: [9]})
            sage: L.is_relatively_complemented()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_relatively_complemented()
            False

        Relatively complemented lattice must be both atomic and coatomic.
        Implication to other direction does not hold::

            sage: L = LatticePoset({0: [1, 2, 3, 4, 5], 1: [6, 7], 2: [6, 8],
            ....:                   3: [7, 8, 9], 4: [9, 11], 5: [9, 10],
            ....:                   6: [10, 11], 7: [12], 8: [12], 9: [12],
            ....:                   10: [12], 11: [12]})
            sage: L.is_atomic() and L.is_coatomic()
            True
            sage: L.is_relatively_complemented()
            False

        We can also get a non-complemented 3-element interval::

            sage: L.is_relatively_complemented(certificate=True)
            (False, (1, 6, 11))

        .. SEEALSO::

            - Weaker properties: :meth:`is_sectionally_complemented`,
              :meth:`is_cosectionally_complemented`, :meth:`is_isoform`
            - Stronger properties: :meth:`is_geometric`

        TESTS::

            sage: [posets.ChainPoset(i).is_relatively_complemented() for
            ....:  i in range(5)]
            [True, True, True, False, False]

        Usually a lattice that is not relatively complemented contains elements
        `l`, `m`, and `u` such that `r(l) + 1 = r(m) = r(u) - 1`, where `r` is
        the rank function and `m` is the only element in the interval `[l, u]`.
        We construct an example where this does not hold::

            sage: B3 = posets.BooleanLattice(3)
            sage: B5 = posets.BooleanLattice(5)
            sage: B3 = B3.subposet([e for e in B3 if e not in [0, 7]])
            sage: B5 = B5.subposet([e for e in B5 if e not in [0, 31]])
            sage: B3 = B3.hasse_diagram()
            sage: B5 = B5.relabel(lambda x: x+10).hasse_diagram()
            sage: G = B3.union(B5)
            sage: G.add_edge(B3.sources()[0], B5.neighbors_in(B5.sinks()[0])[0])
            sage: L = LatticePoset(Poset(G).with_bounds())
            sage: L.is_relatively_complemented()
            False

        Confirm that :trac:`22292` is fixed::

            sage: L = LatticePoset(DiGraph('IYOS`G?CE?@?C?_@??'))
            sage: L.is_relatively_complemented(certificate=True)
            (False, (7, 8, 9))
        """
        from sage.misc.flatten import flatten
        from collections import Counter

        # Work directly with Hasse diagram
        H = self._hasse_diagram
        n = H.order()
        if n < 3:
            return (True, None) if certificate else True

        # Quick check: the lattice must be atomic and coatomic.
        if not certificate:
            if H.out_degree(0) != H.in_degree().count(1):
                return False
            if H.in_degree(n - 1) != H.out_degree().count(1):
                return False

        for e1 in range(n - 1):
            C = Counter(flatten([H.neighbors_out(e2) for e2 in H.neighbors_out(e1)]))
            for e3, c in C.items():
                if c == 1 and len(H.closed_interval(e1, e3)) == 3:
                    if not certificate:
                        return False
                    for e2 in H.neighbors_in(e3):
                        if e2 in H.neighbors_out(e1):
                            break
                    return (False, (self._vertex_to_element(e1),
                                    self._vertex_to_element(e2),
                                    self._vertex_to_element(e3)))
        return (True, None) if certificate else True

    def is_sectionally_complemented(self, certificate=False):
        """
        Return ``True`` if the lattice is sectionally complemented, and
        ``False`` otherwise.

        A lattice is sectionally complemented if all intervals from
        the bottom element interpreted as sublattices are complemented
        lattices.

        INPUT:

        - ``certificate`` -- (default: ``False``) Whether to return
          a certificate if the lattice is not sectionally complemented.

        OUTPUT:

        - If ``certificate=False`` return ``True`` or ``False``.
          If ``certificate=True`` return either ``(True, None)``
          or ``(False, (t, e))``, where `t` is an element so that in the
          sublattice from the bottom element to `t` has no complement
          for element `e`.

        EXAMPLES:

        Smallest examples of a complemented but not sectionally complemented
        lattice and a sectionally complemented but not relatively complemented
        lattice::

            sage: L = posets.PentagonPoset()
            sage: L.is_complemented()
            True
            sage: L.is_sectionally_complemented()
            False

            sage: L = LatticePoset({0: [1, 2, 3], 1: [4], 2: [4], 3: [5], 4: [5]})
            sage: L.is_sectionally_complemented()
            True
            sage: L.is_relatively_complemented()
            False

        Getting a certificate::

            sage: L = LatticePoset(DiGraph('HYOgC?C@?C?G@??'))
            sage: L.is_sectionally_complemented(certificate=True)
            (False, (6, 1))

        .. SEEALSO::

            - Dual property: :meth:`is_cosectionally_complemented`
            - Weaker properties: :meth:`is_complemented`, :meth:`is_atomic`,
              :meth:`is_regular`
            - Stronger properties: :meth:`is_relatively_complemented`

        TESTS::

            sage: [posets.ChainPoset(i).is_sectionally_complemented() for i in range(5)]
            [True, True, True, False, False]
        """
        # Quick check: every sectionally complemented lattice is atomic.
        if not certificate and not self.is_atomic():
            return False

        H = self._hasse_diagram
        mt = H.meet_matrix()
        n = H.order()-1
        for e in range(2, n+1):
            t = n
            for lc in H.neighbors_in(e):
                t = mt[t, lc]
                if t == 0:
                    break
            else:
                if certificate:
                    return (False, (self[e], self[t]))
                return False
        return (True, None) if certificate else True

    def breadth(self, certificate=False):
        r"""
        Return the breadth of the lattice.

        The breadth of a lattice is the largest integer `n` such that
        any join of elements `x_1, x_2, \ldots, x_{n+1}` is join of a
        proper subset of `x_i`.

        This can be also characterized by sublattices: a lattice
        of breadth at least `n` contains a sublattice isomorphic to the
        Boolean lattice of `2^n` elements.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return the pair `(b, a)` where `b` is
          the breadth and `a` is an antichain such that the join of `a`
          differs from the join of any proper subset of `a`.
          If ``certificate=False`` return just the breadth.

        EXAMPLES::

            sage: D10 = posets.DiamondPoset(10)
            sage: D10.breadth()
            2

            sage: B3 = posets.BooleanLattice(3)
            sage: B3.breadth()
            3
            sage: B3.breadth(certificate=True)
            (3, [1, 2, 4])

        ALGORITHM:

        For a lattice to have breadth at least `n`, it must have an
        `n`-element antichain `A` with join `j`. Element `j` must
        cover at least `n` elements. There must also be `n-2` levels
        of elements between `A` and `j`.  So we start by searching
        elements that could be our `j` and then just check possible
        antichains `A`.

        .. NOTE::

            Prior to version 8.1 this function returned just an
            antichain with ``certificate=True``.

        TESTS::

            sage: posets.ChainPoset(0).breadth()
            0
            sage: posets.ChainPoset(1).breadth()
            1
        """
        # A place for optimization: Adding a doubly irreducible element to
        # a lattice does not change the breadth, except from 1 to 2.
        # Hence we could start by removing double irreducibles.

        from sage.combinat.subsets_pairwise import PairwiseCompatibleSubsets

        # First check if breadth is zero (empty lattice) or one (a chain).
        n = self.cardinality()
        if n == 0:
            return (0, []) if certificate else 0
        if self.is_chain():
            return (1, [self.bottom()]) if certificate else 1
        # Breadth is at least two.

        # Work directly with the Hasse diagram
        H = self._hasse_diagram

        # Helper function: Join of elements in the list L.
        jn = H.join_matrix()

        def join(L):
            j = 0
            for i in L:
                j = jn[i, j]
            return j

        indegs = [H.in_degree(i) for i in range(n)]
        max_breadth = max(indegs)

        for B in range(max_breadth, 1, -1):
            for j in H:
                if indegs[j] < B:
                    continue

                # Get elements more than B levels below it.
                too_close = set(H.breadth_first_search(j,
                                                      neighbors=H.neighbors_in,
                                                      distance=B - 2))
                elems = [e for e in H.order_ideal([j]) if e not in too_close]

                achains = PairwiseCompatibleSubsets(elems,
                                                    H.are_incomparable)
                achains_n = achains.elements_of_depth_iterator(B)

                for A in achains_n:
                    if join(A) == j:
                        if all(join(A[:i] + A[i + 1:]) != j for i in range(B)):
                            if certificate:
                                return (B, [self._vertex_to_element(e) for e in A])
                            else:
                                return B
        assert False, "BUG: breadth() in lattices.py have an error."

    def complements(self, element=None):
        r"""
        Return the list of complements of an element in the lattice,
        or the dictionary of complements for all elements.

        Elements `x` and `y` are complements if their meet and join
        are respectively the bottom and the top element of the lattice.

        INPUT:

        - ``element`` -- an element of the lattice whose complement is
          returned. If ``None`` (default) then dictionary of
          complements for all elements having at least one
          complement is returned.

        EXAMPLES::

            sage: L=LatticePoset({0:['a','b','c'], 'a':[1], 'b':[1], 'c':[1]})
            sage: C = L.complements()

        Let us check that 'a' and 'b' are complements of each other::

            sage: 'a' in C['b']
            True
            sage: 'b' in C['a']
            True

        Full list of complements::

            sage: L.complements() # random order
            {0: [1], 1: [0], 'a': ['b', 'c'], 'b': ['c', 'a'], 'c': ['b', 'a']}

            sage: L=LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.complements() # random order
            {0: [4], 4: [0]}
            sage: L.complements(1)
            []

        .. SEEALSO:: :meth:`is_complemented`

        TESTS::

            sage: L=LatticePoset({0:['a','b','c'], 'a':[1], 'b':[1], 'c':[1]})
            sage: for v,v_complements in L.complements().items():
            ....:     for v_c in v_complements:
            ....:         assert L.meet(v,v_c) == L.bottom()
            ....:         assert L.join(v,v_c) == L.top()

            sage: posets.ChainPoset(0).complements()
            {}
            sage: posets.ChainPoset(1).complements()
            {0: [0]}
            sage: posets.ChainPoset(2).complements()
            {0: [1], 1: [0]}
        """
        if element is None:
            n = self.cardinality()
            if n == 1:
                return {self[0]: [self[0]]}
            jn = self.join_matrix()
            mt = self.meet_matrix()
            zero = 0
            one = n - 1
            c = [[] for x in range(n)]
            for x in range(n):
                for y in range(x, n):
                    if jn[x][y] == one and mt[x][y] == zero:
                        c[x].append(y)
                        c[y].append(x)

            comps = {}
            for i in range(n):
                if c[i]:
                    comps[self._vertex_to_element(i)] = (
                        [self._vertex_to_element(x) for x in c[i]])
            return comps

        # Looking for complements of one element.
        if element not in self:
            raise ValueError("element (=%s) not in poset" % element)
        return [x for x in self
                if self.meet(x, element) == self.bottom() and
                self.join(x, element) == self.top()]

    def is_pseudocomplemented(self, certificate=False):
        r"""
        Return ``True`` if the lattice is pseudocomplemented, and ``False``
        otherwise.

        A lattice is (meet-)pseudocomplemented if every element `e` has a
        pseudocomplement `e^\star`, i.e. the greatest element such that
        the meet of `e` and `e^\star` is the bottom element.

        See :wikipedia:`Pseudocomplement`.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where ``e`` is an element without a
          pseudocomplement. If ``certificate=False`` return ``True``
          or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 5], 2: [3, 6], 3: [4], 4: [7],
            ....:                   5: [6], 6: [7]})
            sage: L.is_pseudocomplemented()
            True

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5, 6], 3: [6], 4: [7],
            ....:                   5: [7], 6: [7]})
            sage: L.is_pseudocomplemented()
            False
            sage: L.is_pseudocomplemented(certificate=True)
            (False, 3)

        .. SEEALSO::

            - Dual property: :meth:`is_join_pseudocomplemented`
            - Stronger properties: :meth:`is_meet_semidistributive`
            - Other: :meth:`~sage.combinat.posets.lattices.FiniteMeetSemilattice.pseudocomplement()`.

        ALGORITHM:

        According to [Cha92]_ a lattice is pseudocomplemented if and
        only if every atom has a pseudocomplement. So we only check those.

        TESTS::

            sage: LatticePoset({}).is_pseudocomplemented()
            True
        """
        H = self._hasse_diagram
        if H.order() == 0:
            if certificate:
                return (True, None)
            return True
        for e in H.neighbor_out_iterator(0):
            if H.kappa(e) is None:
                if certificate:
                    return (False, self._vertex_to_element(e))
                return False
        if certificate:
            return (True, None)
        return True

    def is_join_pseudocomplemented(self, certificate=False):
        """
        Return ``True`` if the lattice is join-pseudocomplemented, and
        ``False`` otherwise.

        A lattice is join-pseudocomplemented if every element `e` has a
        join-pseudocomplement `e'`, i.e. the least element such that
        the join of `e` and `e'` is the top element.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where ``e`` is an element without a
          join-pseudocomplement. If ``certificate=False`` return ``True``
          or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 5], 2: [3, 6], 3: [4], 4: [7],
            ....:                   5: [6], 6: [7]})
            sage: L.is_join_pseudocomplemented()
            True

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5, 6], 3: [6], 4: [7],
            ....:                   5: [7], 6: [7]})
            sage: L.is_join_pseudocomplemented()
            False
            sage: L.is_join_pseudocomplemented(certificate=True)
            (False, 4)

        .. SEEALSO::

            - Dual property: :meth:`is_pseudocomplemented`
            - Stronger properties: :meth:`is_join_semidistributive`

        TESTS::

            sage: LatticePoset({}).is_pseudocomplemented()
            True
        """
        H = self._hasse_diagram
        if H.order() == 0:
            if certificate:
                return (True, None)
            return True
        for e in H.neighbor_in_iterator(H.order() - 1):
            if H.kappa_dual(e) is None:
                if certificate:
                    return (False, self._vertex_to_element(e))
                return False
        if certificate:
            return (True, None)
        return True

    def skeleton(self):
        r"""
        Return the skeleton of the lattice.

        The lattice is expected to be pseudocomplemented.

        The *skeleton* of a pseudocomplemented lattice `L`, where `^*` is
        the pseudocomplementation operation, is the subposet induced by
        `\{e^* \mid e \in L\}`. Actually this poset is a Boolean lattice.

        EXAMPLES::

            sage: D12 = posets.DivisorLattice(12)
            sage: S = D12.skeleton(); S
            Finite lattice containing 4 elements
            sage: S.cover_relations()
            [[1, 3], [1, 4], [3, 12], [4, 12]]

            sage: T4 = posets.TamariLattice(4)
            sage: T4.skeleton().is_isomorphic(posets.BooleanLattice(3))
            True

        .. SEEALSO:: :meth:`sage.combinat.posets.lattices.FiniteMeetSemilattice.pseudocomplement`.

        TESTS::

            sage: posets.ChainPoset(0).skeleton()
            Finite lattice containing 0 elements
            sage: posets.ChainPoset(1).skeleton()
            Finite lattice containing 1 elements
            sage: posets.ChainPoset(2).skeleton()
            Finite lattice containing 2 elements
            sage: posets.ChainPoset(3).skeleton()
            Finite lattice containing 2 elements

            sage: L = posets.BooleanLattice(3)
            sage: L == L.skeleton()
            True

            sage: posets.DiamondPoset(5).skeleton()
            Traceback (most recent call last):
            ...
            ValueError: lattice is not pseudocomplemented
        """
        # TODO: What about non-facade lattices and lattices with
        # given linear extension?
        if self.cardinality() < 3:
            return self
        elms = [self._vertex_to_element(v) for v in
                self._hasse_diagram.skeleton()]
        return LatticePoset(self.subposet(elms))

    def is_orthocomplemented(self, unique=False):
        """
        Return ``True`` if the lattice admits an orthocomplementation, and
        ``False`` otherwise.

        An orthocomplementation of a lattice is a function defined for
        every element `e` and marked as `e^{\\bot}` such that
        1) they are complements, i.e. `e \\vee e^{\\bot}` is the top element
        and `e \\wedge e^{\\bot}` is the bottom element, 2) it is involution,
        i.e. `{(e^{\\bot})}^{\\bot} = e`, and 3) it is order-reversing, i.e.
        if `a < b` then `b^{\\bot} < a^{\\bot}`.

        INPUT:

        - ``unique``, a Boolean -- If ``True``, return ``True`` only
          if the lattice has exactly one orthocomplementation. If
          ``False`` (the default), return ``True`` when the lattice
          has at least one orthocomplementation.

        EXAMPLES::

            sage: D5 = posets.DiamondPoset(5)
            sage: D5.is_orthocomplemented()
            False

            sage: D6 = posets.DiamondPoset(6)
            sage: D6.is_orthocomplemented()
            True
            sage: D6.is_orthocomplemented(unique=True)
            False

            sage: hexagon = LatticePoset({0:[1, 2], 1:[3], 2:[4], 3:[5], 4:[5]})
            sage: hexagon.is_orthocomplemented(unique=True)
            True

        .. SEEALSO::

            - Weaker properties: :meth:`is_complemented`,
              :meth:`~sage.categories.finite_posets.FinitePosets.ParentMethods.is_self_dual`

        TESTS::

            sage: [posets.ChainPoset(i).is_orthocomplemented() for i in range(4)]
            [True, True, True, False]
        """
        it = self._hasse_diagram.orthocomplementations_iterator()
        try:
            next(it)
            if not unique:
                return True
        except StopIteration:
            return False
        try:
            next(it)
            return False
        except StopIteration:
            return True
        raise AssertionError("bug in is_orthocomplemented()")

    def is_atomic(self, certificate=False):
        r"""
        Return ``True`` if the lattice is atomic, and ``False`` otherwise.

        A lattice is atomic if every element can be written as a join of atoms.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where `e` is a join-irreducible element
          that is not an atom. If ``certificate=False`` return
          ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3:[5], 4:[6], 5:[6]})
            sage: L.is_atomic()
            True

            sage: L = LatticePoset({0: [1, 2], 1: [3], 2: [3], 3:[4]})
            sage: L.is_atomic()
            False
            sage: L.is_atomic(certificate=True)
            (False, 4)

        TESTS::

            sage: LatticePoset({}).is_atomic()
            True

        .. NOTE::

            See [EnumComb1]_, Section 3.3 for a discussion of atomic lattices.

        .. SEEALSO::

            - Dual property: :meth:`~FiniteLatticePoset.is_coatomic`
            - Stronger properties: :meth:`is_sectionally_complemented`
            - Mutually exclusive properties: :meth:`is_vertically_decomposable`
        """
        if not certificate:
            return (self.cardinality() == 0 or
                    self._hasse_diagram.out_degree(0) ==
                    self._hasse_diagram.in_degree().count(1))
        if self.cardinality() < 3:
            return (True, None)
        H = self._hasse_diagram
        atoms = set(H.neighbors_out(0))
        for v in H:
            if H.in_degree(v) == 1 and v not in atoms:
                return (False, self._vertex_to_element(v))
        return (True, None)

    def is_coatomic(self, certificate=False):
        r"""
        Return ``True`` if the lattice is coatomic, and ``False`` otherwise.

        A lattice is coatomic if every element can be written as a meet
        of coatoms; i.e. if the dual of the lattice is atomic.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, e)``, where `e` is a meet-irreducible element
          that is not a coatom. If ``certificate=False`` return
          ``True`` or ``False``.

        EXAMPLES::

            sage: L = posets.BooleanLattice(3)
            sage: L.is_coatomic()
            True

            sage: L = LatticePoset({1: [2], 2: [3, 4], 3: [5], 4:[5]})
            sage: L.is_coatomic()
            False
            sage: L.is_coatomic(certificate=True)
            (False, 1)

        TESTS::

            sage: LatticePoset({}).is_coatomic()
            True

        .. SEEALSO::

            - Dual property: :meth:`~FiniteLatticePoset.is_atomic`
            - Stronger properties: :meth:`is_cosectionally_complemented`
            - Mutually exclusive properties: :meth:`is_vertically_decomposable`
        """
        n = self.cardinality()
        if not certificate:
            if n == 0:
                return True
            return (self._hasse_diagram.in_degree(n - 1) ==
                    self._hasse_diagram.out_degree().count(1))

        if self.cardinality() < 3:
            return (True, None)
        H = self._hasse_diagram
        coatoms = set(H.neighbors_in(n - 1))
        for v in H:
            if H.out_degree(v) == 1 and v not in coatoms:
                return (False, self._vertex_to_element(v))
        return (True, None)

    def is_geometric(self):
        """
        Return ``True`` if the lattice is geometric, and ``False`` otherwise.

        A lattice is geometric if it is both atomic and upper semimodular.

        EXAMPLES:

        Canonical example is the lattice of partitions of finite set
        ordered by refinement::

            sage: L = posets.SetPartitions(4)
            sage: L.is_geometric()
            True

        Smallest example of geometric lattice that is not modular::

            sage: L = LatticePoset(DiGraph('K]?@g@S?q?M?@?@?@?@?@?@??'))
            sage: L.is_geometric()
            True
            sage: L.is_modular()
            False

        Two non-examples::

            sage: L = LatticePoset({1:[2, 3, 4], 2:[5, 6], 3:[5], 4:[6], 5:[7], 6:[7]})
            sage: L.is_geometric()  # Graded, but not upper semimodular
            False
            sage: L = posets.ChainPoset(3)
            sage: L.is_geometric()  # Modular, but not atomic
            False

        .. SEEALSO::

            - Weaker properties: :meth:`is_upper_semimodular`, :meth:`is_relatively_complemented`

        TESTS::

            sage: LatticePoset({}).is_geometric()
            True
            sage: LatticePoset({1:[]}).is_geometric()
            True
        """
        return self.is_atomic() and self.is_upper_semimodular()

    def is_planar(self):
        r"""
        Return ``True`` if the lattice is *upward* planar, and ``False``
        otherwise.

        A lattice is upward planar if its Hasse diagram has a planar drawing in
        the `\mathbb{R}^2` plane, in such a way that `x` is strictly below `y`
        (on the vertical axis) whenever `x<y` in the lattice.

        Note that the scientific literature on posets often omits "upward" and
        shortens it to "planar lattice" (e.g. [GW2014]_), which can cause
        confusion with the notion of graph planarity in graph theory.

        .. NOTE::

            Not all lattices which are planar -- in the sense of graph planarity
            -- admit such a planar drawing (see example below).

        ALGORITHM:

        Using the result from [Platt1976]_, this method returns its result by
        testing that the Hasse diagram of the lattice is planar (in the sense of
        graph theory) when an edge is added between the top and bottom elements.

        EXAMPLES:

        The Boolean lattice of `2^3` elements is not upward planar, even if
        its covering relations graph is planar::

            sage: B3 = posets.BooleanLattice(3)
            sage: B3.is_planar()
            False
            sage: G = B3.cover_relations_graph()
            sage: G.is_planar()
            True

        Ordinal product of planar lattices is obviously planar. Same does
        not apply to Cartesian products::

            sage: P = posets.PentagonPoset()
            sage: Pc = P.product(P)
            sage: Po = P.ordinal_product(P)
            sage: Pc.is_planar()
            False
            sage: Po.is_planar()
            True

        .. SEEALSO::

            - Weaker properties: :meth:`is_dismantlable`

        TESTS::

            sage: posets.ChainPoset(0).is_planar()
            True
            sage: posets.ChainPoset(1).is_planar()
            True
        """
        # The 8-element Boolean lattice is the smallest non-planar lattice.
        if self.cardinality() < 8:
            return True
        g = self._hasse_diagram.copy(immutable=False)
        g.add_edge(0, self.cardinality() - 1)
        return g.is_planar()

    def is_modular(self, L=None, certificate=False):
        r"""
        Return ``True`` if the lattice is modular and ``False`` otherwise.

        An element `b` of a lattice is *modular* if

        .. MATH::

            x \vee (a \wedge b) = (x \vee a) \wedge b

        for every element `x \leq b` and `a`. A lattice is modular if every
        element is modular. There are other equivalent definitions, see
        :wikipedia:`Modular_lattice`.

        With the parameter ``L`` this can be used to check that
        some subset of elements are all modular.

        INPUT:

        - ``L`` -- (default: ``None``) a list of elements to check being
          modular, if ``L`` is ``None``, then this checks the entire lattice

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (x, a, b))``, where `a`, `b` and `x` are elements
          of the lattice such that `x < b` but
          `x \vee (a \wedge b) \neq (x \vee a) \wedge b`. If also
          `L` is given then `b` in the certificate will be an element
          of `L`. If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_modular()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_modular()
            False

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_modular(certificate=True)
            (False, (2, 6, 4))
            sage: [L.is_modular([x]) for x in L]
            [True, True, False, True, True, False, True]

        .. SEEALSO::

            - Weaker properties: :meth:`is_upper_semimodular`,
              :meth:`is_lower_semimodular`, :meth:`is_supersolvable`
            - Stronger properties: :meth:`is_distributive`
            - Other: :meth:`is_modular_element`

        TESTS::

            sage: all(posets.ChainPoset(i).is_modular() for i in range(4))
            True

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_modular(L=[1, 4, 2], certificate=True)
            (False, (2, 6, 4))
            sage: L.is_modular(L=[1, 6, 2], certificate=True)
            (False, (3, 4, 6))
        """
        if not certificate and L is None:
            return self.is_upper_semimodular() and self.is_lower_semimodular()

        if certificate and L is None:
            tmp = self.is_lower_semimodular(certificate=True)
            if not tmp[0]:
                a, b = tmp[1]
                t = self.meet(a, b)
                for x in self.upper_covers(t):
                    if self.is_less_than(x, b):
                        return (False, (x, a, b))
                    if self.is_less_than(x, a):
                        return (False, (x, b, a))
            tmp = self.is_upper_semimodular(certificate=True)
            if not tmp[0]:
                x, a = tmp[1]
                t = self.join(x, a)
                for b in self.lower_covers(t):
                    if self.is_greater_than(b, x):
                        return (False, (x, a, b))
                    if self.is_greater_than(b, a):
                        return (False, (a, x, b))
            return (True, None)

        # L is not None
        for b in L:
            for x in self.principal_lower_set(b):
                for a in self:
                    if (self.join(x, self.meet(a, b)) !=
                        self.meet(self.join(x, a), b)):
                        if certificate:
                            return (False, (x, a, b))
                        return False
        if certificate:
            return (True, None)
        return True

    def is_modular_element(self, x):
        r"""
        Return ``True`` if ``x`` is a modular element and ``False`` otherwise.

        INPUT:

        - ``x`` -- an element of the lattice

        An element `x` in a lattice `L` is *modular* if `x \leq b` implies

        .. MATH::

            x \vee (a \wedge b) = (x \vee a) \wedge b

        for every `a, b \in L`.

        EXAMPLES::

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_modular()
            False
            sage: [L.is_modular_element(x) for x in L]
            [True, True, False, True, True, False, True]

        .. SEEALSO::

            - Weaker properties: :meth:`is_left_modular_element`
            - Other: :meth:`is_modular` to check modularity for the full
              lattice or some set of elements
        """
        return self.is_modular([x])

    def is_left_modular_element(self, x):
        r"""
        Return ``True`` if ``x`` is a left modular element
        and ``False`` otherwise.

        INPUT:

        - ``x`` -- an element of the lattice

        An element `x` in a lattice `L` is *left modular* if

        .. MATH::

            (y \vee x) \wedge z = y \vee (x \wedge z)

        for every `y \leq z \in L`.

        It is enough to check this condition on all cover relations `y < z`.

        EXAMPLES::

            sage: P = posets.PentagonPoset()
            sage: [i for i in P if P.is_left_modular_element(i)]
            [0, 2, 3, 4]

        .. SEEALSO::

            - Stronger properties: :meth:`is_modular_element`
        """
        return all(self.meet(self.join(y, x), z) ==
                   self.join(y, self.meet(x, z))
                   for y, z in self.cover_relations_iterator())

    def is_upper_semimodular(self, certificate=False):
        r"""
        Return ``True`` if the lattice is upper semimodular and
        ``False`` otherwise.

        A lattice is upper semimodular if any pair of elements with
        a common lower cover have also a common upper cover.

        INPUT:

        - ``certificate`` -- (default: ``False``) Whether to return
          a certificate if the lattice is not upper semimodular.

        OUTPUT:

        - If ``certificate=False`` return ``True`` or ``False``.
          If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (a, b))``, where `a` and `b` covers their meet but
          are not covered by their join.

        See :wikipedia:`Semimodular_lattice`

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_upper_semimodular()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_upper_semimodular()
            False

            sage: L = LatticePoset(posets.IntegerPartitions(4))
            sage: L.is_upper_semimodular()
            True

            sage: L = LatticePoset({1:[2, 3, 4], 2: [5], 3:[5, 6], 4:[6], 5:[7], 6:[7]})
            sage: L.is_upper_semimodular(certificate=True)
            (False, (4, 2))

        .. SEEALSO::

            - Dual property: :meth:`is_lower_semimodular`
            - Weaker properties: :meth:`~sage.combinat.posets.posets.FinitePoset.is_graded`
            - Stronger properties: :meth:`is_modular`,
              :meth:`is_join_distributive`, :meth:`is_geometric`

        TESTS::

            sage: all(posets.ChainPoset(i).is_upper_semimodular() for i in range(5))
            True
        """
        nonmodular = self._hasse_diagram.find_nonsemimodular_pair(upper=True)
        if nonmodular is None:
            return (True, None) if certificate else True
        if certificate:
            return (False, (self._vertex_to_element(nonmodular[0]),
                            self._vertex_to_element(nonmodular[1])))
        return False

    def is_lower_semimodular(self, certificate=False):
        r"""
        Return ``True`` if the lattice is lower semimodular and
        ``False`` otherwise.

        A lattice is lower semimodular if any pair of elements with
        a common upper cover have also a common lower cover.

        INPUT:

        - ``certificate`` -- (default: ``False``) Whether to return
          a certificate if the lattice is not lower semimodular.

        OUTPUT:

        - If ``certificate=False`` return ``True`` or ``False``.
          If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (a, b))``, where `a` and `b` are covered by their
          join but do no cover their meet.

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

            sage: L = LatticePoset(DiGraph('IS?`?AAOE_@?C?_@??'))
            sage: L.is_lower_semimodular(certificate=True)
            (False, (4, 2))

        .. SEEALSO::

            - Dual property: :meth:`is_upper_semimodular`
            - Weaker properties: :meth:`~sage.combinat.posets.posets.FinitePoset.is_graded`
            - Stronger properties: :meth:`is_modular`,
              :meth:`is_meet_distributive`
        """
        nonmodular = self._hasse_diagram.find_nonsemimodular_pair(upper=False)
        if nonmodular is None:
            return (True, None) if certificate else True
        if certificate:
            return (False, (self._vertex_to_element(nonmodular[0]),
                            self._vertex_to_element(nonmodular[1])))
        return False

    def is_supersolvable(self, certificate=False):
        r"""
        Return ``True`` if the lattice is supersolvable, and
        ``False`` otherwise.

        A lattice `L` is *supersolvable* if there exists a maximal chain `C`
        such that every `x \in C` is a modular element in `L`. Equivalent
        definition is that the sublattice generated by `C` and any other chain
        is distributive.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(False, None)`` or
          ``(True, C)``, where ``C`` is a maximal chain of modular elements.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = posets.DiamondPoset(5)
            sage: L.is_supersolvable()
            True

            sage: L = posets.PentagonPoset()
            sage: L.is_supersolvable()
            False

            sage: L = LatticePoset({1:[2,3],2:[4,5],3:[5,6],4:[7],5:[7],6:[7]})
            sage: L.is_supersolvable()
            True
            sage: L.is_supersolvable(certificate=True)
            (True, [1, 2, 5, 7])
            sage: L.is_modular()
            False

            sage: L = LatticePoset({0: [1, 2, 3, 4], 1: [5, 6, 7],
            ....:                   2: [5, 8, 9], 3: [6, 8, 10], 4: [7, 9, 10],
            ....:                   5: [11], 6: [11], 7: [11], 8: [11],
            ....:                   9: [11], 10: [11]})
            sage: L.is_supersolvable()
            False

        .. SEEALSO::

            - Weaker properties: :meth:`~sage.combinat.posets.posets.FinitePoset.is_graded`
            - Stronger properties: :meth:`is_modular`

        TESTS::

            sage: LatticePoset().is_supersolvable()
            True
        """
        from sage.misc.cachefunc import cached_function

        not_ok = (False, None) if certificate else False

        if not self.is_ranked():
            return not_ok

        if self.cardinality() == 0:
            if certificate:
                return (True, [])
            return True

        H = self._hasse_diagram
        mt = H.meet_matrix()
        jn = H.join_matrix()
        height = self.height()
        n = H.order()
        cur = H.maximal_elements()[0]
        cert = [cur]
        next_ = [H.neighbor_in_iterator(cur)]

        @cached_function
        def is_modular_elt(a):
            return all(H._rank[a] + H._rank[b] ==
                       H._rank[mt[a, b]] + H._rank[jn[a, b]]
                       for b in range(n))

        if not is_modular_elt(cur):
            return not_ok
        while len(next_) < height:
            try:
                cur = next(next_[-1])
            except StopIteration:
                next_.pop()
                cert.pop()
                if not next_:
                    return not_ok
                continue
            if is_modular_elt(cur):
                next_.append(H.neighbor_in_iterator(cur))
                cert.append(cur)
        if certificate:
            return (True, [self._vertex_to_element(e) for e in reversed(cert)])
        return True

    def vertical_composition(self, other, labels='pairs'):
        r"""
        Return the vertical composition of the lattice with ``other``.

        Let `L` and `K` be lattices and `b_K` the bottom element
        of `K`. The vertical composition of `L` and `K` is the ordinal
        sum of `L` and `K \setminus \{b_K\}`. Informally said this is
        lattices "glued" together with a common element.

        Mathematically, it is only defined when `L` and `K` have no
        common element; here we force that by giving them different
        names in the resulting poset.

        INPUT:

        - ``other`` -- a lattice

        - ``labels`` -- a string (default ``'pairs'``); can be one of
          the following:

          * ``'pairs'`` - each element ``v`` in this poset will be
            named ``(0, v)`` and each element ``u`` in ``other`` will
            be named ``(1, u)`` in the result
          * ``'integers'`` - the elements of the result will be
            relabeled with consecutive integers

        EXAMPLES::

            sage: L = LatticePoset({'a': ['b', 'c'], 'b': ['d'], 'c': ['d']})
            sage: K = LatticePoset({'e': ['f', 'g'], 'f': ['h'], 'g': ['h']})
            sage: M = L.vertical_composition(K)
            sage: M.list()
            [(0, 'a'), (0, 'b'), (0, 'c'), (0, 'd'), (1, 'f'), (1, 'g'), (1, 'h')]
            sage: M.upper_covers((0, 'd'))
            [(1, 'f'), (1, 'g')]

            sage: C2 = posets.ChainPoset(2)
            sage: M3 = posets.DiamondPoset(5)
            sage: L = C2.vertical_composition(M3, labels='integers')
            sage: L.cover_relations()
            [[0, 1], [1, 2], [1, 3], [1, 4], [2, 5], [3, 5], [4, 5]]

        .. SEEALSO::

            :meth:`vertical_decomposition`,
            :meth:`sage.combinat.posets.posets.FinitePoset.ordinal_sum`

        TESTS::

            sage: C0 = LatticePoset()
            sage: C1 = LatticePoset({'a': []})
            sage: C2 = LatticePoset({'b': ['c']})
            sage: C2.vertical_composition(C2)
            Finite lattice containing 3 elements
            sage: C0.vertical_composition(C0)
            Finite lattice containing 0 elements
            sage: C0.vertical_composition(C1).list()
            [(1, 'a')]
            sage: C1.vertical_composition(C0).list()
            [(0, 'a')]
            sage: C1.vertical_composition(C1).list()
            [(0, 'a')]
            sage: C1.vertical_composition(C2).list()
            [(0, 'a'), (1, 'c')]
            sage: C2.vertical_composition(C1).list()
            [(0, 'b'), (0, 'c')]
        """
        from copy import copy

        # Todo: This and ordinal_sum() of posets could keep
        # distinguished linear extension, if it is defined
        # for both posets/lattices. That can be done after
        # trac ticket #21607.

        if labels not in ['integers', 'pairs']:
            raise ValueError("labels must be either 'pairs' or 'integers'")
        if not isinstance(self, FiniteLatticePoset):
            raise ValueError("the input is not a finite lattice")
        if self._is_facade != other._is_facade:
            raise ValueError("mixing facade and non-facade lattices is not defined")

        if labels == 'integers':
            g_self = copy(self._hasse_diagram)
            g_other = other._hasse_diagram.copy(immutable=False)
            n = max(g_self.order(), 1)  # max() takes care of empty 'self'.
            g_other.relabel(lambda v: v + n - 1)
            g_result = g_self.union(g_other)
            return FiniteLatticePoset(g_result, elements=range(g_result.order()),
                                      facade=self._is_facade, category=FiniteLatticePosets())

        if self.cardinality() == 0:
            return other.relabel(lambda e: (1, e))
        S = other.subposet([e for e in other if e != other.bottom()])
        return LatticePoset(self.ordinal_sum(S), facade=self._is_facade)

    def vertical_decomposition(self, elements_only=False):
        r"""
        Return sublattices from the vertical decomposition of the lattice.

        Let `d_1, \ldots, d_n` be elements (excluding the top and bottom
        elements) comparable to every element of the lattice. Let `b`
        be the bottom element and `t` be the top element. This function
        returns either a list `d_1, \ldots, d_n`, or the list of
        intervals `[b, d_1], [d_1, d_2], \ldots, [d_{n-1}, d_n], [d_n,
        t]` as lattices.

        Informally said, this returns the lattice split into parts at
        every single-element "cutting point".

        INPUT:

        - ``elements_only`` - if ``True``, return the list of decomposing
          elements as defined above; if ``False`` (the default),
          return the list of sublattices so that the lattice is a
          vertical composition of them.

        EXAMPLES:

        Number 6 is divided by 1, 2, and 3, and it divides 12, 18 and 36::

            sage: L = LatticePoset( ([1, 2, 3, 6, 12, 18, 36],
            ....:     attrcall("divides")) )
            sage: parts = L.vertical_decomposition()
            sage: [lat.list() for lat in parts]
            [[1, 2, 3, 6], [6, 12, 18, 36]]
            sage: L.vertical_decomposition(elements_only=True)
            [6]

        .. SEEALSO::

            :meth:`vertical_composition`,
            :meth:`is_vertically_decomposable`

        TESTS::

            sage: [posets.ChainPoset(i).vertical_decomposition(elements_only=True)
            ....:     for i in range(5)]
            [[], [], [], [1], [1, 2]]
        """
        if self.cardinality() <= 2:
            if not elements_only:
                return [self]
            else:
                return []
        if elements_only:
            return [self[e] for e in
                    self._hasse_diagram.vertical_decomposition(return_list=True)]
        elms = ( [0] +
                 self._hasse_diagram.vertical_decomposition(return_list=True) +
                 [self.cardinality() - 1] )
        n = len(elms)
        result = []
        for i in range(n - 1):
            result.append(LatticePoset(
                self.subposet([self[e] for e in range(elms[i], elms[i + 1] + 1)])))
        return result

    def is_vertically_decomposable(self, certificate=False):
        r"""
        Return ``True`` if the lattice is vertically decomposable, and
        ``False`` otherwise.

        A lattice is vertically decomposable if it has an element that
        is comparable to all elements and is neither the bottom nor
        the top element.

        Informally said, a lattice is vertically decomposable if it
        can be seen as two lattices "glued" by unifying the top
        element of first lattice to the bottom element of second one.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - If ``certificate=True`` return either ``(False, None)`` or
          ``(True, e)``, where `e` is an element that is comparable to all
          other elements and is neither the bottom nor the top element.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: posets.TamariLattice(4).is_vertically_decomposable()
            False
            sage: L = LatticePoset( ([1, 2, 3, 6, 12, 18, 36],
            ....:     attrcall("divides")) )
            sage: L.is_vertically_decomposable()
            True
            sage: L.is_vertically_decomposable(certificate=True)
            (True, 6)

        .. SEEALSO::

            - Weaker properties: :meth:`is_subdirectly_reducible`
            - Mutually exclusive properties: :meth:`is_atomic`, :meth:`is_coatomic`,
              :meth:`is_regular`
            - Other: :meth:`vertical_decomposition`

        TESTS::

            sage: [posets.ChainPoset(i).is_vertically_decomposable() for i in
            ....:  range(5)]
            [False, False, False, True, True]
        """
        e = self._hasse_diagram.vertical_decomposition()
        if e is None:
            if certificate:
                return (False, None)
            return False
        if certificate:
            return (True, self._vertex_to_element(e))
        return True

    def sublattice(self, elms):
        r"""
        Return the smallest sublattice containing elements on the given list.

        INPUT:

        - ``elms`` -- a list of elements of the lattice.

        EXAMPLES::

            sage: L=LatticePoset(( [], [[1,2],[1,17],[1,8],[2,3],[2,22],[2,5],[2,7],[17,22],[17,13],[8,7],[8,13],[3,16],[3,9],[22,16],[22,18],[22,10],[5,18],[5,14],[7,9],[7,14],[7,10],[13,10],[16,6],[16,19],[9,19],[18,6],[18,33],[14,33],[10,19],[10,33],[6,4],[19,4],[33,4]] ))
            sage: L.sublattice([14, 13, 22]).list()
            [1, 2, 8, 7, 14, 17, 13, 22, 10, 33]

            sage: L = posets.BooleanLattice(3)
            sage: L.sublattice([3,5,6,7])
            Finite lattice containing 8 elements
        """
        gens_remaining = set(elms)
        current_set = set()

        # We add elements one by one in 'current_set'.
        #
        # When adding a point g to 'current_set', we add to 'gens_remaining' all
        # meet/join obtained from g and another point of 'current_set'
        while gens_remaining:
            g = gens_remaining.pop()
            if g in current_set:
                continue
            for x in current_set:
                gens_remaining.add(self.join(x, g))
                gens_remaining.add(self.meet(x, g))
            current_set.add(g)

        return LatticePoset(self.subposet(current_set))

    def is_sublattice(self, other):
        """
        Return ``True`` if the lattice is a sublattice of ``other``,
        and ``False`` otherwise.

        Lattice `K` is a sublattice of `L` if `K` is an (induced) subposet
        of `L` and closed under meet and join of `L`.

        .. NOTE::

            This method does not check whether the lattice is a
            *isomorphic* (i.e., up to relabeling) sublattice of ``other``,
            but only if ``other`` directly contains the lattice as an
            sublattice.

        EXAMPLES:

        A pentagon sublattice in a non-modular lattice::

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [5, 6], 4: [7], 5: [7], 6: [7]})
            sage: N5 = LatticePoset({1: [2, 6], 2: [4], 4: [7], 6: [7]})
            sage: N5.is_sublattice(L)
            True

        This pentagon is a subposet but not closed under join, hence not a sublattice::

            sage: N5_ = LatticePoset({1: [2, 3], 2: [4], 3: [7], 4: [7]})
            sage: N5_.is_induced_subposet(L)
            True
            sage: N5_.is_sublattice(L)
            False

        .. SEEALSO::

            :meth:`isomorphic_sublattices_iterator`

        TESTS::

            sage: E = LatticePoset({})
            sage: P = posets.PentagonPoset()
            sage: E.is_sublattice(P)
            True

            sage: P1 = LatticePoset({'a':['b']})
            sage: P2 = P1.dual()
            sage: P1.is_sublattice(P2)
            False

            sage: P = MeetSemilattice({0: [1]})
            sage: E.is_sublattice(P)
            True
            sage: P = JoinSemilattice({0: [1]})
            sage: E.is_sublattice(P)
            True
            sage: P = Poset({0: [1]})
            sage: E.is_sublattice(P)
            True
        """
        try:
            o_meet = other.meet
            o_join = other.join
        except (AttributeError):
            raise TypeError('other is not a lattice')
        if not self.is_induced_subposet(other):
            return False

        n = self.cardinality()
        for i in range(n):
            for j in range(i):
                if (o_meet(self[i], self[j]) not in self or
                    o_join(self[i], self[j]) not in self):
                    return False
        return True

    def sublattices(self):
        """
        Return all sublattices of the lattice.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2:[5], 3:[5, 6], 4:[6],
            ....:                   5:[7], 6:[7]})
            sage: sublats = L.sublattices(); len(sublats)
            54
            sage: sublats[3]
            Finite lattice containing 4 elements
            sage: sublats[3].list()
            [1, 2, 3, 5]

        TESTS:

        A subposet that is a lattice but not a sublattice::

            sage: L = LatticePoset({1: [2, 3], 2:[4], 3:[4], 4:[5]})
            sage: sl = L.sublattices()
            sage: LatticePoset({1: [2, 3], 2:[5], 3:[5]}) in sl
            False

        `n`-element chain has `2^n` sublattices (also tests empty lattice)::

            sage: [len(posets.ChainPoset(n).sublattices()) for n in range(4)]
            [1, 2, 4, 8]
        """
        return [LatticePoset(self.subposet(map(self._vertex_to_element, elms)))
                for elms in self._hasse_diagram.sublattices_iterator(set(), 0)]

    def sublattices_lattice(self, labels='lattice'):
        """
        Return the lattice of sublattices.

        Every element of the returned lattice is a sublattice and
        they are ordered by containment; that is, atoms are one-element
        lattices, coatoms are maximal sublattices of the original
        lattice and so on.

        INPUT:

        - ``labels`` -- string; can be one of the following:

          * ``'lattice'`` (default) elements of the lattice will be
            lattices that correspond to sublattices of the original lattice

          * ``'tuple'`` - elements are tuples of elements of the sublattices
            of the original lattice

          * ``'integer'`` - elements are plain integers

        EXAMPLES::

            sage: D4 = posets.DiamondPoset(4)
            sage: sll = D4.sublattices_lattice(labels='tuple')
            sage: sll.coatoms()  # = maximal sublattices of the original lattice
            [(0, 1, 3), (0, 2, 3)]

            sage: L = posets.DivisorLattice(12)
            sage: sll = L.sublattices_lattice()
            sage: L.is_dismantlable() == (len(sll.atoms()) == sll.rank())
            True

        TESTS::

            sage: E = posets.ChainPoset(0)
            sage: E.sublattices_lattice()
            Finite lattice containing 1 elements

            sage: C3 = posets.ChainPoset(3)
            sage: sll = C3.sublattices_lattice(labels='integer')
            sage: sll.is_isomorphic(posets.BooleanLattice(3))
            True
        """
        if labels not in ['lattice', 'tuple', 'integer']:
            raise ValueError("labels must be one of 'lattice', 'tuple' or 'integer'")
        sublats = [frozenset(x) for x in self._hasse_diagram.sublattices_iterator(set(), 0)]
        L = LatticePoset([sublats, lambda a, b: a != b and a.issubset(b)])
        if labels == 'integer':
            return L.canonical_label()
        L = L.relabel(lambda x: tuple(self._vertex_to_element(y) for y in x))
        if labels == 'lattice':
            return L.relabel(self.sublattice)
        return L

    def isomorphic_sublattices_iterator(self, other):
        """
        Return an iterator over the sublattices of the lattice isomorphic to ``other``.

        INPUT:

        - other --  a finite lattice

        EXAMPLES:

        A non-modular lattice contains a pentagon sublattice::

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [5, 6], 4: [7], 5: [7], 6: [7]})
            sage: L.is_modular()
            False
            sage: N5 = posets.PentagonPoset()
            sage: N5_in_L = next(L.isomorphic_sublattices_iterator(N5)); N5_in_L
            Finite lattice containing 5 elements
            sage: N5_in_L.list()
            [1, 3, 6, 4, 7]

        A divisor lattice is modular, hence does not contain the
        pentagon as sublattice, even if it has the pentagon
        subposet::

            sage: D12 = posets.DivisorLattice(12)
            sage: D12.has_isomorphic_subposet(N5)
            True
            sage: list(D12.isomorphic_sublattices_iterator(N5))
            []

        .. SEEALSO::

            :meth:`sage.combinat.posets.posets.FinitePoset.isomorphic_subposets_iterator`

        .. WARNING::

            This function will return same sublattice as many times as
            there are automorphism on it. This is due to
            :meth:`~sage.graphs.generic_graph.GenericGraph.subgraph_search_iterator`
            returning labelled subgraphs.

        TESTS::

            sage: E = LatticePoset()
            sage: P = LatticePoset({1: []})
            sage: list(N5.isomorphic_sublattices_iterator(E))
            [Finite lattice containing 0 elements]
            sage: len(list(N5.isomorphic_sublattices_iterator(P)))
            5
        """
        from itertools import combinations
        if not isinstance(other, FiniteLatticePoset):
            raise TypeError('the input is not a finite lattice')
        H = self._hasse_diagram
        mt = H.meet_matrix()
        jn = H.join_matrix()
        self_closure = H.transitive_closure()
        other_closure = other._hasse_diagram.transitive_closure()
        for g in self_closure.subgraph_search_iterator(other_closure, induced=True):
            if all(mt[a, b] in g and jn[a, b] in g for a, b in combinations(g, 2)):
                yield self.sublattice([self._vertex_to_element(v) for v in g])

    def maximal_sublattices(self):
        r"""
        Return maximal (proper) sublattices of the lattice.

        EXAMPLES::

            sage: L = LatticePoset(( [], [[1,2],[1,17],[1,8],[2,3],[2,22],
            ....:                         [2,5],[2,7],[17,22],[17,13],[8,7],
            ....:                         [8,13],[3,16],[3,9],[22,16],[22,18],
            ....:                         [22,10],[5,18],[5,14],[7,9],[7,14],
            ....:                         [7,10],[13,10],[16,6],[16,19],[9,19],
            ....:                         [18,6],[18,33],[14,33],[10,19],
            ....:                         [10,33],[6,4],[19,4],[33,4]] ))
            sage: maxs = L.maximal_sublattices()
            sage: len(maxs)
            7
            sage: sorted(maxs[0].list())
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 16, 18, 19, 22, 33]
        """
        n = self.cardinality()
        if n < 2:
            return []
        if n == 2:
            return [self.sublattice([self.bottom()]), self.sublattice([self.top()])]
        return [self.sublattice([self[x] for x in d]) for d in self._hasse_diagram.maximal_sublattices()]

    def frattini_sublattice(self):
        r"""
        Return the Frattini sublattice of the lattice.

        The Frattini sublattice `\Phi(L)` is the intersection of all
        proper maximal sublattices of `L`. It is also the set of
        "non-generators" - if the sublattice generated by set `S` of
        elements is whole lattice, then also `S \setminus \Phi(L)`
        generates whole lattice.

        EXAMPLES::

            sage: L = LatticePoset(( [], [[1,2],[1,17],[1,8],[2,3],[2,22],
            ....:                         [2,5],[2,7],[17,22],[17,13],[8,7],
            ....:                         [8,13],[3,16],[3,9],[22,16],[22,18],
            ....:                         [22,10],[5,18],[5,14],[7,9],[7,14],
            ....:                         [7,10],[13,10],[16,6],[16,19],[9,19],
            ....:                         [18,6],[18,33],[14,33],[10,19],
            ....:                         [10,33],[6,4],[19,4],[33,4]] ))
            sage: sorted(L.frattini_sublattice().list())
            [1, 2, 4, 10, 19, 22, 33]
        """
        return LatticePoset(self.subposet([self[x] for x in
                self._hasse_diagram.frattini_sublattice()]))

    def moebius_algebra(self, R):
        """
        Return the Möbius algebra of ``self`` over ``R``.

        OUTPUT:

        An instance of :class:`sage.combinat.posets.moebius_algebra.MoebiusAlgebra`.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: L.moebius_algebra(QQ)
            Moebius algebra of Finite lattice containing 16 elements over Rational Field
        """
        from sage.combinat.posets.moebius_algebra import MoebiusAlgebra
        return MoebiusAlgebra(R, self)

    def quantum_moebius_algebra(self, q=None):
        """
        Return the quantum Möbius algebra of ``self`` with parameter ``q``.

        INPUT:

        - ``q`` -- (optional) the deformation parameter `q`

        OUTPUT:

        An instance of :class:`sage.combinat.posets.moebius_algebra.QuantumMoebiusAlgebra`.

        EXAMPLES::

            sage: L = posets.BooleanLattice(4)
            sage: L.quantum_moebius_algebra()
            Quantum Moebius algebra of Finite lattice containing 16 elements
             with q=q over Univariate Laurent Polynomial Ring in q over Integer Ring
        """
        from sage.combinat.posets.moebius_algebra import QuantumMoebiusAlgebra
        return QuantumMoebiusAlgebra(self, q)

    def day_doubling(self, S):
        r"""
        Return the lattice with Alan Day's doubling construction of subset `S`.

        The subset `S` is assumed to be convex (i.e. if
        `a, c \in S` and `a < b < c` in the lattice, then `b \in S`)
        and connected (i.e. if `a, b \in S` then there is a chain
        `a=e_1, e_2, \ldots, e_n=b` such that `e_i` either covers or
        is covered by `e_{i+1}`).

        .. image:: ../../../media/day-doubling.png

        Alan Day's doubling construction is a specific extension of
        the lattice. Here we formulate it in a format more suitable
        for computation.

        Let `L` be a lattice and `S` a convex subset of it. The resulting
        lattice `L[S]` has elements `(e, 0)` for each `e \in L` and
        `(e, 1)` for each `e \in S`. If `x \le y` in `L`, then in the
        new lattice we have

        * `(x, 0), (x, 1) \le (y, 0), (y, 1)`
        * `(x, 0) \le (x, 1)`

        INPUT:

        - ``S`` -- a subset of the lattice

        EXAMPLES::

            sage: L = LatticePoset({1: ['a', 'b', 2], 'a': ['c'], 'b': ['c', 'd'],
            ....:                   2: [3], 'c': [4], 'd': [4], 3: [4]})
            sage: L2 = L.day_doubling(['a', 'b', 'c', 'd']); L2
            Finite lattice containing 12 elements
            sage: set(L2.upper_covers((1, 0))) == set([(2, 0), ('a', 0), ('b', 0)])
            True
            sage: set(L2.upper_covers(('b', 0))) == set([('d', 0), ('b', 1), ('c', 0)])
            True

        .. SEEALSO::

            :meth:`is_constructible_by_doublings`

        TESTS::

            sage: L2._hasse_diagram.is_isomorphic(DiGraph('KSCH??_BO?g?_?@?G?@?A?@??'))
            True

            sage: L = LatticePoset({'a': ['b']})
            sage: set(L.day_doubling([]).list()) == set([('a', 0), ('b', 0)])
            True
            sage: set(L.day_doubling(['a', 'b']).list()) == set([('a', 0), ('a', 1), ('b', 0), ('b', 1)])
            True
        """
        # Rationale for naming of elements: a lattice can have
        # elements 1, (1, 1), (1, (1, 1)) and so on. We can't just
        # make a copy of S with elements (s, 1).

        # The construction could be defined for any convex
        # subset S, but we assume that the user made an error
        # if S is not also connected.

        S = sorted(set(S))
        S_ = [self._element_to_vertex(e) for e in S]
        if not self._hasse_diagram.is_convex_subset(S_):
            raise ValueError("subset S is not convex")
        if not self._hasse_diagram.subgraph(S_).is_connected():
            raise ValueError("subset S is not connected")

        g = self.hasse_diagram()
        g.relabel(lambda e: (e, 0))

        for e in S:
            g.add_edge((e, 0), (e, 1))
            for e_up in self.upper_covers(e):
                if e_up in S:
                    g.add_edge((e, 1), (e_up, 1))
                else:
                    g.delete_edge((e, 0), (e_up, 0))
                    g.add_edge((e, 1), (e_up, 0))

        return LatticePoset(g)

    def adjunct(self, other, a, b):
        r"""
        Return the adjunct of the lattice by ``other`` on the pair `(a, b)`.

        It is assumed that `a < b` but `b` does not cover `a`.

        The adjunct of a lattice `K` to `L` with respect to pair
        `(a, b)` of `L` is defined such that `x < y` if

        - `x, y \in K` and `x < y` in `K`,
        - `x, y \in L` and `x < y` in `L`,
        - `x \in L`, `y \in K` and `x \le a` in `L`, or
        - `x \in K`, `y \in L` and `b \le y` in `L`.

        Informally this can be seen as attaching the lattice `K` to `L`
        as a new block between `a` and `b`. Dismantlable lattices are exactly
        those that can be created from chains with this function.

        Mathematically, it is only defined when `L` and `K` have no
        common element; here we force that by giving them different
        names in the resulting lattice.

        EXAMPLES::

            sage: Pnum = posets.PentagonPoset()
            sage: Palp = Pnum.relabel(lambda x: chr(ord('a')+x))
            sage: PP = Pnum.adjunct(Palp, 0, 3)
            sage: PP.atoms()
            [(0, 1), (0, 2), (1, 'a')]
            sage: PP.coatoms()
            [(0, 3), (0, 1)]

        TESTS::

            sage: P = posets.PentagonPoset()
            sage: E = LatticePoset()
            sage: PE = P.adjunct(E, 0, 3); PE.is_isomorphic(P)
            True
            sage: PE.bottom()
            (0, 0)
            sage: C4 = posets.ChainPoset(4)
            sage: C1 = posets.ChainPoset(1)
            sage: C4.adjunct(C1, 0, 3).is_isomorphic(P)
            True
        """
        if not isinstance(other, FiniteLatticePoset):
            raise ValueError("other is not a finite lattice")
        if not self.is_greater_than(b, a):
            raise ValueError("element %s is not greater than %s in the lattice" % (b, a))
        if self.covers(a, b):
            raise ValueError("element %s covers element %s in the lattice" % (b, a))

        if other.cardinality() == 0:
            return self.relabel(lambda e: (0, e))

        g_self = self.hasse_diagram()
        g_other = other.hasse_diagram()
        g = g_self.disjoint_union(g_other, labels='pairs')
        g.add_edge((0, a), (1, other.bottom()))
        g.add_edge((1, other.top()), (0, b))
        return LatticePoset(g)

    def center(self):
        """
        Return the center of the lattice.

        An element of a lattice is *central* if it is neutral and has a
        complement. The subposet induced by central elements is a *center* of
        the lattice. Actually it is a Boolean lattice.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [6, 7], 3: [8, 9, 7],
            ....:                   4: [5, 6], 5: [8, 10], 6: [10], 7: [13, 11],
            ....:                   8: [13, 12], 9: [11, 12], 10: [13],
            ....:                   11: [14], 12: [14], 13: [14]})
            sage: C = L.center(); C
            Finite lattice containing 4 elements
            sage: C.cover_relations()
            [[1, 2], [1, 12], [2, 14], [12, 14]]

            sage: L = posets.DivisorLattice(60)
            sage: sorted(L.center().list())
            [1, 3, 4, 5, 12, 15, 20, 60]

        .. SEEALSO::

            :meth:`neutral_elements`, :meth:`complements`

        TESTS::

            sage: LatticePoset().center()
            Finite lattice containing 0 elements

            sage: posets.ChainPoset(1).center()
            Finite lattice containing 1 elements

            sage: L = posets.BooleanLattice(3)
            sage: L.center() == L
            True
        """
        neutrals = self.neutral_elements()
        comps = self.complements()
        return self.sublattice([e for e in neutrals if e in comps])

    def is_dismantlable(self, certificate=False):
        r"""
        Return ``True`` if the lattice is dismantlable, and ``False``
        otherwise.

        An `n`-element lattice `L_n` is dismantlable if there is a sublattice
        chain `L_{n-1} \supset L_{n-2}, \supset \cdots, \supset L_0` so that
        every `L_i` is a sublattice of `L_{i+1}` with one element less, and
        `L_0` is the empty lattice. In other words, a dismantlable lattice
        can be reduced to empty lattice removing doubly irreducible
        element one by one.

        INPUT:

        - ``certificate`` (boolean) -- Whether to return a certificate.

          * If ``certificate = False`` (default), returns ``True`` or
            ``False`` accordingly.

          * If ``certificate = True``, returns:

            * ``(True, elms)`` when the lattice is dismantlable, where
              ``elms`` is elements listed in a possible removing order.

            * ``(False, crown)`` when the lattice is not dismantlable,
              where ``crown`` is a subposet of `2k` elements
              `a_1, \ldots, a_k, b_1, \ldots, b_k` with covering
              relations `a_i \lessdot b_i` and `a_i \lessdot b_{i+1}`
              for `i \in [1, \ldots, k-1]`, and `a_k \lessdot b_1`.

        EXAMPLES::

            sage: DL12 = LatticePoset((divisors(12), attrcall("divides")))
            sage: DL12.is_dismantlable()
            True
            sage: DL12.is_dismantlable(certificate=True)
            (True, [4, 2, 1, 3, 6, 12])

            sage: B3 = posets.BooleanLattice(3)
            sage: B3.is_dismantlable()
            False
            sage: B3.is_dismantlable(certificate=True)
            (False, Finite poset containing 6 elements)

        Every planar lattice is dismantlable. Converse is not true::

            sage: L = LatticePoset( ([], [[0, 1], [0, 2], [0, 3], [0, 4],
            ....:                         [1, 7], [2, 6], [3, 5], [4, 5],
            ....:                         [4, 6], [4, 7], [5, 8], [6, 8],
            ....:                         [7, 8]]) )
            sage: L.is_dismantlable()
            True
            sage: L.is_planar()
            False

        .. SEEALSO::

            - Stronger properties: :meth:`is_planar`
            - Weaker properties: :meth:`is_sublattice_dismantlable`

        TESTS::

            sage: posets.ChainPoset(0).is_dismantlable()
            True
            sage: posets.ChainPoset(1).is_dismantlable()
            True

            sage: L = LatticePoset(DiGraph('L@_?W?E?@CCO?A?@??_?O?Jw????C?'))
            sage: L.is_dismantlable()
            False
            sage: c = L.is_dismantlable(certificate=True)[1]
            sage: (3 in c, 12 in c, 9 in c)
            (True, False, True)
        """
        from sage.graphs.digraph import DiGraph
        from copy import copy

        H = copy(self._hasse_diagram)
        cert = []
        # Smallest lattice that is not dismantlable is the
        # Boolean lattice with 2^3=8 elements. Hence the limit 7.
        limit = 0 if certificate else 7

        while H.order() > limit:
            for e in H:
                i = H.in_degree(e)
                o = H.out_degree(e)
                if i < 2 and o < 2:
                    if certificate:
                        cert.append(e)
                    if i == 1 and o == 1:  # Remove inside the lattice
                        lower = H.neighbors_in(e)[0]
                        upper = H.neighbors_out(e)[0]
                        H.delete_vertex(e)
                        if upper not in H.depth_first_search(lower):
                            H.add_edge(lower, upper)
                    else:  # Remove the top or bottom element
                        H.delete_vertex(e)
                    break
            else:
                if not certificate:
                    return False
                k = 3
                while True:
                    crown = DiGraph({i: [k + i, k + (i + 1) % k]
                                     for i in range(k)})
                    sg = H.transitive_closure().subgraph_search(crown, True)
                    if sg:
                        elms = [self[e] for e in sg]
                        return (False, self.subposet(elms))
                    k += 1
        if not certificate:
            return True
        return (True, [self[e] for e in cert])

    def is_interval_dismantlable(self, certificate=False):
        """
        Return ``True`` if the lattice is interval dismantlable, and
        ``False`` otherwise.

        An interval dismantling is a subdivision of a lattice to a principal
        upper set and a principal lower set. A lattice is *interval
        dismantlable* if it can be decomposed into 1-element lattices by
        consecutive interval distmantlings.

        A lattice is *minimally interval non-dismantlable* if it is not
        interval dismantlable, but all of its sublattices are interval
        dismantlable.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - if ``certificate=False``, return only ``True`` or ``False``
        - if ``certificate=True``, return either

          * ``(True, list)`` where ``list`` is a nested list showing the
            decomposition; for example ``list[1][0]`` is a lower part of
            upper part of the lattice when decomposed twice.
          * ``(False, M)`` where `M` is a minimally interval non-dismantlable
            sublattice of the lattice.

        EXAMPLES::

            sage: L1 = LatticePoset({1: [2, 3], 3: [4, 5], 2: [6], 4: [6], 5: [6]})
            sage: L1.is_interval_dismantlable()
            True

            sage: L2 = LatticePoset({1: [2, 3, 4, 5], 2: [6], 3: [6], 4: [6],
            ....:                    5: [6, 7], 6: [8], 7: [9, 10], 8:[10], 9:[10]})
            sage: L2.is_interval_dismantlable()
            False

        To get certificates::

            sage: L1.is_interval_dismantlable(certificate=True)
            (True, [[[1], [2]], [[[3], [5]], [[4], [6]]]])
            sage: L2.is_interval_dismantlable(certificate=True)
            (False, Finite lattice containing 5 elements)

        .. SEEALSO::

            - Stronger properties: :meth:`is_join_semidistributive`,
              :meth:`is_meet_semidistributive`
            - Weaker properties: :meth:`is_sublattice_dismantlable`

        TESTS::

            sage: LatticePoset().is_interval_dismantlable()
            True
            sage: LatticePoset().is_interval_dismantlable(certificate=True)
            (True, [])
        """
        def minimal_non_int_dismant(L):
            """
            Return a minimally interval non-dismantlable sublattice.

            Assumes that L is interval non-dismantlable.
            """
            while True:
                H = L._hasse_diagram
                for sl in H.sublattices_iterator(set(), 0):
                    if len(sl) < 5 or len(sl) == H.order():
                        continue
                    L_ = L.sublattice([L._vertex_to_element(v) for v in sl])
                    if not L_.is_interval_dismantlable():
                        L = L_
                        break
                else:
                    return L

        def recursive_is_interval_dismantlable(self):
            if self.cardinality() == 1:
                return (True, [self[0]])

            mp = self.meet_primes()
            if not mp:
                return (False, self)

            e = mp[0]
            e_lower = self.principal_lower_set(e)

            # First check the lower part...
            L_down = LatticePoset(self.subposet(e_lower))
            result_down = recursive_is_interval_dismantlable(L_down)
            if not result_down[0]:
                return result_down

            # ...and then check the upper.
            L_up = LatticePoset(self.subposet([x for x in self if x not in e_lower]))
            result_up = recursive_is_interval_dismantlable(L_up)
            if not result_up[0]:
                return result_up

            # OK, both parts are dismantlable, so then is the whole lattice.
            return (True, [result_down[1], result_up[1]])

        if self.cardinality() == 0:
            return (True, []) if certificate else True

        result = recursive_is_interval_dismantlable(self)
        if result[0]:
            return result if certificate else True
        return (False, minimal_non_int_dismant(self)) if certificate else False

    def is_sublattice_dismantlable(self):
        """
        Return ``True`` if the lattice is sublattice dismantlable, and
        ``False`` otherwise.

        A sublattice dismantling is a subdivision of a lattice into two
        non-empty sublattices. A lattice is  *sublattice dismantlable*
        if it can be decomposed into 1-element lattices by consecutive
        sublattice dismantlings.

        EXAMPLES:

        The smallest non-example is this (and the dual)::

            sage: P = Poset({1: [11, 12, 13], 2: [11, 14, 15],
            ....:            3: [12, 14, 16], 4: [13, 15, 16]})
            sage: L = LatticePoset(P.with_bounds())
            sage: L.is_sublattice_dismantlable()
            False

        Here we adjoin a (double-irreducible-)dismantlable lattice
        as a part to an interval-dismantlable lattice::

            sage: B3 = posets.BooleanLattice(3)
            sage: N5 = posets.PentagonPoset()
            sage: L = B3.adjunct(N5, 1, 7)
            sage: L.is_dismantlable(), L.is_interval_dismantlable()
            (False, False)
            sage: L.is_sublattice_dismantlable()
            True

        .. SEEALSO::

            - Stronger properties: :meth:`is_dismantlable`, :meth:`is_interval_dismantlable`

        TESTS::

            sage: LatticePoset().is_sublattice_dismantlable()
            True

        .. TODO::

            Add a certificate-option.
        """
        # Todo: This can be made much faster, if we don't regenerate meet- and
        # join-matrices every time, but instead remove some rows and columns
        # from them.

        from sage.combinat.subset import Subsets

        # All lattices up to 11 elements have been tested to be sublattice
        # dismantlable, hence the limit.
        if self.cardinality() <= 11:
            return True

        for low in self:
            if len(self.lower_covers(low)) <= 1:
                for up in self.principal_upper_set(low):
                    if len(self.upper_covers(up)) <= 1:
                        if low == self.bottom() and up == self.top():
                            continue
                        S1 = self.interval(low, up)
                        S2 = [e for e in self if e not in S1]
                        if all(self.meet(a, b) in S2 and
                               self.join(a, b) in S2
                               for a, b in Subsets(S2, 2)):
                            sub1 = self.sublattice(S1)
                            sub2 = self.sublattice(S2)
                            return (sub1.is_sublattice_dismantlable() and
                                    sub2.is_sublattice_dismantlable())

        return False

    def is_subdirectly_reducible(self, certificate=False):
        r"""
        Return ``True`` if the lattice is subdirectly reducible.

        A lattice `M` is a *subdirect product* of `K` and `L` if it
        is a sublattice of `K \times L`. Lattice `M` is *subdirectly
        reducible* if there exists such lattices `K` and `L` so that
        `M` is not a sublattice of either.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate

        OUTPUT:

        - if ``certificate=False``, return only ``True`` or ``False``
        - if ``certificate=True``, return either

          * ``(True, (K, L))`` such that the lattice is isomorphic to
            a sublattice of `K \times L`.
          * ``(False, (a, b))``, where `a` and `b` are elements that are
            in the same congruence class for every nontrivial congruence
            of the lattice. Special case: If the lattice has zero or one element,
            return ``(False, None)``.

        EXAMPLES::

            sage: N5 = posets.PentagonPoset()
            sage: N5.is_subdirectly_reducible()
            False

            sage: hex = LatticePoset({1: [2, 3], 2: [4], 3: [5], 4: [6], 5: [6]})
            sage: hex.is_subdirectly_reducible()
            True

            sage: hex.is_subdirectly_reducible(certificate=True)
            (True,
             (Finite lattice containing 5 elements, Finite lattice containing 5 elements))

            sage: N5.is_subdirectly_reducible(certificate=True)
            (False, (2, 3))
            sage: res, cert = hex.is_subdirectly_reducible(certificate=True)
            sage: cert[0].is_isomorphic(N5)
            True

        .. SEEALSO::

            - Stronger properties: :meth:`is_distributive`,
              :meth:`is_vertically_decomposable`
            - Other: :meth:`subdirect_decomposition`

        TESTS::

            sage: [posets.ChainPoset(i).is_subdirectly_reducible() for i in range(5)]
            [False, False, False, True, True]
        """
        H = self._hasse_diagram
        A = H.atoms_of_congruence_lattice()

        if not certificate:
            return len(A) > 1

        # Kind of special cases. How should we define this for empty,
        # one-element and two-element lattices?
        if self.cardinality() < 2:
            return (False, None)

        if len(A) == 1:
            for a in A[0]:
                if len(a) > 1:
                    x, y = min(a), max(a)
                    return (False, (self._vertex_to_element(x),
                                    self._vertex_to_element(y)))

        H_closure = H.transitive_closure()
        a0 = [min(v) for v in A[0]]
        a1 = [min(v) for v in A[1]]
        K0 = LatticePoset(H_closure.subgraph(a0).transitive_reduction())
        K1 = LatticePoset(H_closure.subgraph(a1).transitive_reduction())
        return (True, (K0, K1))

    def canonical_meetands(self, e):
        r"""
        Return the canonical meetands of `e`.

        The canonical meetands of an element `e` in the lattice `L` is the
        subset `S \subseteq L` such that 1) the meet of `S` is `e`, and
        2) if the meet of some other subset `S'` of is also `e`, then for
        every element `s \in S` there is an element `s' \in S'` such that
        `s \ge s'`.

        Informally said this is the set of greatest possible elements
        with given meet. It exists for every element if and only if
        the lattice is meet-semidistributive. Canonical meetands are
        always meet-irreducibles.

        INPUT:

        - ``e`` -- an element of the lattice

        OUTPUT:

        - canonical meetands as a list, if it exists; if not, ``None``

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [5, 6], 4: [6],
            ....:                   5: [7], 6: [7]})
            sage: L.canonical_meetands(1)
            [5, 4]

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [6], 4: [6],
            ....: 5: [6]})
            sage: L.canonical_meetands(1) is None
            True

        .. SEEALSO::

            :meth:`canonical_joinands`

        TESTS::

            LatticePoset({1: []}).canonical_meetands(1)
            [1]
        """
        # Algorithm: Make interval from e to the top element.
        # Now compute kappa function for every atom of that lattice, i.e.
        # kind of "restricted" kappa for elements covering e.
        # This is done implicitly here.
        H = self._hasse_diagram
        e = self._element_to_vertex(e)
        meetands = []
        for a in H.neighbors_out(e):
            above_a = list(H.depth_first_search(a))
            go_up = lambda v: [v_ for v_ in H.neighbors_out(v) if v_ not in above_a]
            result = None
            for v in H.depth_first_search(e, neighbors=go_up):
                if H.out_degree(v) == 1 and next(H.neighbor_out_iterator(v)) in above_a:
                    if result is not None:
                        return None
                    result = v
            meetands.append(result)
        return [self._vertex_to_element(v) for v in meetands]

    def canonical_joinands(self, e):
        r"""
        Return the canonical joinands of `e`.

        The canonical joinands of an element `e` in the lattice `L` is the
        subset `S \subseteq L` such that 1) the join of `S` is `e`, and
        2) if the join of some other subset `S'` of is also `e`, then for
        every element `s \in S` there is an element `s' \in S'` such that
        `s \le s'`.

        Informally said this is the set of lowest possible elements
        with given join. It exists for every element if and only if
        the lattice is join-semidistributive. Canonical joinands are
        always join-irreducibles.

        INPUT:

        - ``e`` -- an element of the lattice

        OUTPUT:

        - canonical joinands as a list, if it exists; if not, ``None``

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [5], 4: [6],
            ....:                   5: [7], 6: [7]})
            sage: L.canonical_joinands(7)
            [3, 4]

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [6], 4: [6],
            ....: 5: [6]})
            sage: L.canonical_joinands(6) is None
            True

        .. SEEALSO::

            :meth:`canonical_meetands`

        TESTS::

            LatticePoset({1: []}).canonical_joinands(1)
            [1]
        """
        # Algorithm: Make dual of interval from the bottom element to e.
        # Now compute kappa function for every atom of that lattice, i.e.
        # kind of "restricted" dual kappa for elements covered by e.
        # This is done implicitly here.
        H = self._hasse_diagram
        e = self._element_to_vertex(e)
        joinands = []
        for a in H.neighbors_in(e):
            below_a = list(H.depth_first_search(a, neighbors=H.neighbors_in))
            go_down = lambda v: [v_ for v_ in H.neighbors_in(v) if v_ not in below_a]
            result = None
            for v in H.depth_first_search(e, neighbors=go_down):
                if H.in_degree(v) == 1 and next(H.neighbor_in_iterator(v)) in below_a:
                    if result is not None:
                        return None
                    result = v
            joinands.append(result)
        return [self._vertex_to_element(v) for v in joinands]

    def is_constructible_by_doublings(self, type):
        r"""
        Return ``True`` if the lattice is constructible by doublings, and
        ``False`` otherwise.

        We call a lattice doubling constructible if it can be constructed
        from the one element lattice by a sequence of Alan Day's doubling
        constructions.

        Lattices constructible by interval doubling are also called
        *bounded*. Lattices constructible by lower and upper pseudo-interval
        are called *lower bounded* and *upper bounded*. Lattices
        constructible by any convex set doubling are called *congruence
        normal*.

        INPUT:

        - ``type`` -- a string; can be one of the following:

          * ``'interval'`` - allow only doublings of an interval
          * ``'lower'`` - allow doublings of lower pseudo-interval; that is, a
            subset of the lattice with a unique minimal element
          * ``'upper'`` - allow doublings of upper pseudo-interval; that is, a
            subset of the lattice with a unique maximal element
          * ``'convex'`` - allow doubling of any convex set
          * ``'any'`` - allow doubling of any set

        EXAMPLES:

        The pentagon can be constructed by doubling intervals; the 5-element
        diamond can not be constructed by any doublings::

            sage: posets.PentagonPoset().is_constructible_by_doublings('interval')
            True

            sage: posets.DiamondPoset(5).is_constructible_by_doublings('any')
            False

        After doubling both upper and lower pseudo-interval a lattice is
        constructible by convex subset doubling::

            sage: L = posets.BooleanLattice(2)
            sage: L = L.day_doubling([0, 1, 2])  # A lower pseudo-interval
            sage: L.is_constructible_by_doublings('interval')
            False
            sage: L.is_constructible_by_doublings('lower')
            True
            sage: L = L.day_doubling([(3,0), (1,1), (2,1)])  # An upper pseudo-interval
            sage: L.is_constructible_by_doublings('upper')
            False
            sage: L.is_constructible_by_doublings('convex')
            True

        An example of a lattice that can be constructed by doublings
        of a non-convex subsets::

            sage: L = LatticePoset(DiGraph('OQC?a?@CO?G_C@?GA?O??_??@?BO?A_?G??C??_?@???'))
            sage: L.is_constructible_by_doublings('convex')
            False
            sage: L.is_constructible_by_doublings('any')
            True

        .. SEEALSO::

            - Stronger properties: :meth:`is_distributive` (doubling by interval),
              :meth:`is_join_semidistributive` (doubling by lower pseudo-intervals),
              :meth:`is_meet_semidistributive` (doubling by upper pseudo-intervals)
            - Mutually exclusive properties: :meth:`is_simple` (doubling by any set)
            - Other: :meth:`day_doubling`

        TESTS::

            sage: LatticePoset().is_constructible_by_doublings('interval')
            True

        The congruence lattice of this lattice has maximal chains satisfying the needed
        property, but also maximal chains not satisfying that; this shows that the code
        can't be optimized to test just some maximal chain::

            sage: L = LatticePoset(DiGraph('QSO?I?_?_GBG??_??a???@?K??A??B???C??s??G??I??@??A??@???'))
            sage: L.is_constructible_by_doublings('convex')
            False
            sage: L.is_constructible_by_doublings('any')
            True

        ALGORITHM:

        According to [HOLM2016]_ a lattice `L` is lower bounded if and only if
        `|\mathrm{Ji}(L)| = |\mathrm{Ji}(\mathrm{Con}\ L)|`, and so dually
        `|\mathrm{Mi}(L)| = |\mathrm{Mi}(\mathrm{Con}\ L)|` in upper bounded
        lattices. The same reference gives a test for being constructible by
        convex or by any subset.
        """
        if type not in ['interval', 'lower', 'upper', 'convex', 'any']:
            raise ValueError("type must be one of 'interval', 'lower', 'upper', 'convex' or 'any'")

        if self.cardinality() < 5:
            return True

        if (type == 'interval' and len(self.join_irreducibles()) !=
            len(self.meet_irreducibles())):
            return False

        if type == 'upper' or type == 'interval':
            H = self._hasse_diagram
            found = set()
            for v in H:
                if H.out_degree(v) == 1:
                    S = frozenset(map(frozenset, H.congruence([[v, next(H.neighbor_out_iterator(v))]])))
                    if S in found:
                        return False
                    found.add(S)
            return True

        if type == 'lower':
            H = self._hasse_diagram
            found = set()
            for v in H:
                if H.in_degree(v) == 1:
                    S = frozenset(map(frozenset, H.congruence([[v, next(H.neighbor_in_iterator(v))]])))
                    if S in found:
                        return False
                    found.add(S)
            return True

        if type == 'convex':
            return self._hasse_diagram.is_congruence_normal()
        # type == 'any'

        def splitting_depth_2(a, b):
            """
            Return ``True`` if every block of `b` is made from
            combining at most two blocks of `a`.
            """
            return all(len([x for x in a if x.issubset(y)]) <= 2 for y in b)

        conL = self.congruences_lattice()
        todo = [conL[0]]
        reachable = []

        while todo:
            e = todo.pop()
            for e_up in conL.upper_covers(e):
                if e_up not in reachable and splitting_depth_2(e, e_up):
                    if len(e_up) == 1:  # = the top of the cong. lattice
                        return True
                    reachable.append(e_up)
                    todo.append(e_up)
        return False

    def is_isoform(self, certificate=False):
        """
        Return ``True`` if the lattice is isoform and ``False`` otherwise.

        A congruence is *isoform* (or *isotype*) if all blocks are isomorphic
        sublattices. A lattice is isoform if it has only isoform
        congruences.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate if the lattice is not isoform

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, C)``, where `C` is a non-isoform congruence as a
          :class:`sage.combinat.set_partition.SetPartition`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1:[2, 3, 4], 2: [5, 6], 3: [6, 7], 4: [7], 5: [8], 6: [8], 7: [8]})
            sage: L.is_isoform()
            True

        Every isoform lattice is (trivially) uniform, but the converse is
        not true::

            sage: L = LatticePoset({1: [2, 3, 6], 2: [4, 5], 3: [5], 4: [9, 8], 5: [7, 8], 6: [9], 7: [10], 8: [10], 9: [10]})
            sage: L.is_isoform(), L.is_uniform()
            (False, True)

            sage: L.is_isoform(certificate=True)
            (False, {{1, 2, 4, 6, 9}, {3, 5, 7, 8, 10}})

        .. SEEALSO::

            - Weaker properties: :meth:`is_uniform`
            - Stronger properties: :meth:`is_simple`,
              :meth:`is_relatively_complemented`
            - Other: :meth:`congruence`

        TESTS::

            sage: [posets.ChainPoset(i).is_isoform() for i in range(5)]
            [True, True, True, False, False]

            sage: posets.DiamondPoset(5).is_isoform()  # Simple, so trivially isoform
            True
        """
        ok = (True, None) if certificate else True

        H = self._hasse_diagram
        if H.order() == 0:
            return ok
        for c in H.congruences_iterator():
            cong = list(c)
            if len(cong) in [1, H.order()]:
                continue
            if not certificate:
                if any(len(x) != len(cong[0]) for x in cong):
                    return False
            d = H.subgraph(cong[0])
            for part in cong:
                if not H.subgraph(part).is_isomorphic(d):
                    if certificate:
                        from sage.combinat.set_partition import SetPartition
                        return (False,
                                SetPartition([[self._vertex_to_element(v) for v in p] for p in cong]))
                    return False
        return ok

    def is_uniform(self, certificate=False):
        """
        Return ``True`` if the lattice is uniform and ``False`` otherwise.

        A congruence is *uniform* if all blocks have equal number
        of elements. A lattice is uniform if it has only uniform
        congruences.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate if the lattice is not uniform

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, C)``, where `C` is a non-uniform congruence as a
          :class:`sage.combinat.set_partition.SetPartition`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [6, 7], 3: [5], 4: [5], 5: [9, 8], 6: [9], 7: [10], 8: [10], 9: [10]})
            sage: L.is_uniform()
            True

        Every uniform lattice is regular, but the converse is not true::

            sage: N6 = LatticePoset({1: [2, 3, 5], 2: [4], 3: [4], 5: [6], 4: [6]})
            sage: N6.is_uniform(), N6.is_regular()
            (False, True)

            sage: N6.is_uniform(certificate=True)
            (False, {{1, 2, 3, 4}, {5, 6}})

        .. SEEALSO::

            - Weaker properties: :meth:`is_regular`
            - Stronger properties: :meth:`is_isoform`
            - Other: :meth:`congruence`

        TESTS::

            sage: [posets.ChainPoset(i).is_uniform() for i in range(5)]
            [True, True, True, False, False]

            sage: posets.DiamondPoset(5).is_uniform()  # Simple, so trivially uniform
            True
        """
        ok = (True, None) if certificate else True

        H = self._hasse_diagram
        if H.order() < 3:
            return ok

        # Check for trivial cases
        x = H._trivial_nonregular_congruence()
        if x is not None:
            if certificate:
                return (False, self.congruence([[self._vertex_to_element(x[0]), self._vertex_to_element(x[1])]]))
            return False
        x = self.is_vertically_decomposable(certificate=True)
        if x[0]:
            if certificate:
                return (False, self.congruence([self.vertical_decomposition()[0]]))
            return False

        for c in H.congruences_iterator():
            cong = list(c)
            n = len(cong[0])
            for part in cong:
                if len(part) != n:
                    if certificate:
                        from sage.combinat.set_partition import SetPartition
                        return (False,
                                SetPartition([[self._vertex_to_element(v) for v in p] for p in c]))
                    return False
        return ok

    def is_regular(self, certificate=False):
        """
        Return ``True`` if the lattice is regular and ``False`` otherwise.

        A congruence of a lattice is *regular* if it is generated
        by any of its parts. A lattice is regular if it has only
        regular congruences.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate if the lattice is not regular

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, (C, p))``, where `C` is a non-regular congruence as a
          :class:`sage.combinat.set_partition.SetPartition` and `p` is a
          congruence class of `C` such that the congruence generated by `p`
          is not `C`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5, 6], 3: [8, 7], 4: [6, 7], 5: [8], 6: [9], 7: [9], 8: [9]})
            sage: L.is_regular()
            True

            sage: N5 = posets.PentagonPoset()
            sage: N5.is_regular()
            False
            sage: N5.is_regular(certificate=True)
            (False, ({{0}, {1}, {2, 3}, {4}}, [0]))

        .. SEEALSO::

            - Stronger properties: :meth:`is_uniform`,
              :meth:`is_sectionally_complemented`,
              :meth:`is_cosectionally_complemented`
            - Mutually exclusive properties: :meth:`is_vertically_decomposable`
            - Other: :meth:`congruence`

        TESTS::

            sage: [posets.ChainPoset(i).is_regular() for i in range(5)]
            [True, True, True, False, False]
        """
        ok = (True, None) if certificate else True

        H = self._hasse_diagram
        if H.order() < 3:
            return ok
        for cong in H.congruences_iterator():
            x = iter(cong.root_to_elements_dict().values())
            ell = len(next(x))
            if all(len(p) == ell for p in x):
                continue
            for part in cong:
                if H.congruence([part]) != cong:
                    if certificate:
                        from sage.combinat.set_partition import SetPartition
                        return (False,
                                (SetPartition([[self._vertex_to_element(v) for v in p] for p in cong]),
                                 [self._vertex_to_element(v) for v in part]))
                    return False
        return ok

    def is_simple(self, certificate=False):
        """
        Return ``True`` if the lattice is simple and ``False`` otherwise.

        A lattice is *simple* if it has no nontrivial congruences; in
        other words, for every two distinct elements `a` and `b` the
        principal congruence generated by `(a, b)` has only one
        component, i.e. the whole lattice.

        INPUT:

        - ``certificate`` -- (default: ``False``) whether to return
          a certificate if the lattice is not simple

        OUTPUT:

        - If ``certificate=True`` return either ``(True, None)`` or
          ``(False, c)``, where `c` is a nontrivial congruence as a
          :class:`sage.combinat.set_partition.SetPartition`.
          If ``certificate=False`` return ``True`` or ``False``.

        EXAMPLES::

            sage: posets.DiamondPoset(5).is_simple()  # Smallest nontrivial example
            True
            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [6], 4: [6], 5: [6]})
            sage: L.is_simple()
            False
            sage: L.is_simple(certificate=True)
            (False, {{1, 3}, {2, 4, 5, 6}})

        Two more examples. First is a non-simple lattice without any
        2-element congruences::

            sage: L = LatticePoset({1: [2, 3, 4], 2: [5], 3: [5], 4: [6, 7],
            ....:                   5: [8], 6: [8], 7: [8]})
            sage: L.is_simple()
            False
            sage: L = LatticePoset({1: [2, 3], 2: [4, 5], 3: [6, 7], 4: [8],
            ....:                   5: [8], 6: [8], 7: [8]})
            sage: L.is_simple()
            True

        .. SEEALSO::

            - Weaker properties: :meth:`is_isoform`
            - Mutually exclusive properties: :meth:`is_constructible_by_doublings`
              (by any set)
            - Other: :meth:`congruence`

        TESTS::

            sage: [posets.ChainPoset(i).is_simple() for i in range(5)]
            [True, True, True, False, False]
        """
        from sage.combinat.set_partition import SetPartition
        cong = self._hasse_diagram.find_nontrivial_congruence()
        if cong is None:
            return (True, None) if certificate else True
        if not certificate:
            return False
        return (False, SetPartition([[self._vertex_to_element(v) for v in s]
                                     for s in cong]))

    def subdirect_decomposition(self):
        r"""
        Return the subdirect decomposition of the lattice.

        The subdirect decomposition of a lattice `L` is the list
        of smaller lattices `L_1, \ldots, L_n` such that `L` is
        a sublattice of `L_1 \times \ldots \times L_n`, none
        of `L_i` can be decomposed further and `L` is not a sublattice
        of  any `L_i`. (Except when the list has only one element, i.e.
        when the lattice is subdirectly irreducible.)

        EXAMPLES::

            sage: posets.ChainPoset(3).subdirect_decomposition()
            [Finite lattice containing 2 elements, Finite lattice containing 2 elements]

            sage: L = LatticePoset({1: [2, 4], 2: [3], 3: [6, 7], 4: [5, 7],
            ....:                   5: [9, 8], 6: [9], 7: [9], 8: [10], 9: [10]})
            sage: Ldecomp = L.subdirect_decomposition()
            sage: [fac.cardinality() for fac in Ldecomp]
            [2, 5, 7]
            sage: Ldecomp[1].is_isomorphic(posets.PentagonPoset())
            True

        TESTS::

            sage: posets.ChainPoset(0).subdirect_decomposition()
            [Finite lattice containing 0 elements]
            sage: posets.ChainPoset(1).subdirect_decomposition()
            [Finite lattice containing 1 elements]
            sage: posets.ChainPoset(2).subdirect_decomposition()
            [Finite lattice containing 2 elements]

        The pentagon is subdirectly irreducible, i.e. the decomposition
        has only one element::

            sage: N5 = posets.PentagonPoset()
            sage: N5.subdirect_decomposition()
            [Finite lattice containing 5 elements]
        """
        H = self._hasse_diagram
        cong_ji, congs = H.principal_congruences_poset()

        if self.cardinality() <= 2 or cong_ji.has_bottom():
            return [self.relabel(self._element_to_vertex_dict)]

        L_ = cong_ji.order_ideals_lattice()
        c = L_.canonical_meetands(L_.bottom())
        L = L_.subposet(L_.order_ideal(c))

        C = {}
        for e in L:
            low = L.lower_covers(e)
            if len(low) == 1:  # a join-irreducible element
                C[e] = congs[max(e, key=cong_ji._element_to_vertex)]
            elif low:  # "extending" congruence to avoid re-computation
                low_0 = min(low, key=lambda x: C[x].number_of_subsets())
                for new_pair in e:
                    if new_pair not in low_0:
                        break
                C[e] = self._hasse_diagram.congruence([new_pair], start=C[low_0])

        decomposing_congruences = [C[m] for m in L.maximal_elements()]
        decomposing_congruences.sort(key=lambda x: x.number_of_subsets())

        result = []
        for congruence in decomposing_congruences:
            part_bottoms = [min(part) for part in congruence]
            F = H.transitive_closure().subgraph(part_bottoms)
            result.append(LatticePoset(F))

        return result

    def congruence(self, S):
        r"""
        Return the congruence generated by set of sets `S`.

        A congruence of a lattice is an equivalence relation `\cong` that is
        compatible with meet and join; i.e. if `a_1 \cong a_2` and
        `b_1 \cong b_2`, then `(a_1 \\vee b_1) \cong (a_2 \\vee b_2)` and
        `(a_1 \wedge b_1) \cong (a_2 \wedge b_2)`.

        By the congruence generated by set of sets `\{S_1, \ldots, S_n\}` we
        mean the least congruence `\cong` such that for every `x, y \in S_i`
        for some `i` we have `x \cong y`.

        INPUT:

        - ``S`` -- a list of lists; list of element blocks that the congruence
          will contain

        OUTPUT:

        Congruence of the lattice as a
        :class:`sage.combinat.set_partition.SetPartition`.

        EXAMPLES::

            sage: L = posets.DivisorLattice(12)
            sage: cong = L.congruence([[1, 3]])
            sage: sorted(sorted(c) for c in cong)
            [[1, 3], [2, 6], [4, 12]]
            sage: L.congruence([[1, 2], [6, 12]])
            {{1, 2, 4}, {3, 6, 12}}

            sage: L = LatticePoset({1: [2, 3], 2: [4], 3: [4], 4: [5]})
            sage: L.congruence([[1, 2]])
            {{1, 2}, {3, 4}, {5}}

            sage: L = LatticePoset({1: [2, 3], 2: [4, 5, 6], 4: [5], 5: [7, 8],
            ....:                   6: [8], 3: [9], 7: [10], 8: [10], 9:[10]})
            sage: cong = L.congruence([[1, 2]])
            sage: cong[0]
            frozenset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10})

        .. SEEALSO:: :meth:`quotient`

        TESTS::

            sage: P = posets.PentagonPoset()
            sage: P.congruence([])
            {{0}, {1}, {2}, {3}, {4}}
            sage: P.congruence([[2]])
            {{0}, {1}, {2}, {3}, {4}}
            sage: P.congruence([[2, 2]])
            {{0}, {1}, {2}, {3}, {4}}
            sage: P.congruence([[0, 4]])
            {{0, 1, 2, 3, 4}}
            sage: LatticePoset().congruence([])
            {}

        "Double zigzag" to ensure that up-down propagation
        works::

            sage: L = LatticePoset(DiGraph('P^??@_?@??B_?@??B??@_?@??B_?@??B??@??A??C??G??O???'))
            sage: sorted(sorted(p) for p in L.congruence([[1,6]]))
            [[0], [1, 6], [2], [3, 8], [4], [5, 10], [7, 12], [9, 14], [11], [13], [15], [16]]

        Simple lattice, i.e. a lattice without any nontrivial congruence::

            sage: L = LatticePoset(DiGraph('GPb_@?OC@?O?'))
            sage: L.congruence([[1,2]])
            {{0, 1, 2, 3, 4, 5, 6, 7}}
        """
        from sage.combinat.set_partition import SetPartition
        S = [[self._element_to_vertex(e) for e in s] for s in S]
        cong = self._hasse_diagram.congruence(S)
        return SetPartition([[self._vertex_to_element(v) for v in s]
                             for s in cong])

    def quotient(self, congruence, labels='tuple'):
        r"""
        Return the quotient lattice by ``congruence``.

        Let `L` be a lattice and `\Theta` be a congruence of `L` with
        congruence classes `\Theta_1, \Theta_2, \ldots`. The quotient
        lattice `L/\Theta` is the lattice with elements
        `\{\Theta_1, \Theta_2, \ldots\}` and meet and join given by the
        original lattice. Explicitly, if `e_1 \in \Theta_1` and
        `e_2 \in \Theta_2`, such that `e_1 \vee e_2 \in \Theta_3` then
        `\Theta_1 \vee \Theta_2 = \Theta_3` in `L/\Theta` and similarly
        for meets.

        INPUT:

        - ``congruence`` -- list of lists; a congruence

        - ``labels`` -- string; the elements of the resulting
          lattice and can be one of the following:

          * ``'tuple'`` - elements are tuples of elements of the original
            lattice
          * ``'lattice'`` - elements are sublattices of the original lattice
          * ``'integer'`` - elements are labeled by integers

        .. WARNING::

            ``congruence`` is expected to be a valid congruence of the
            lattice. This is *not* checked.

        EXAMPLES::

            sage: L = posets.PentagonPoset()
            sage: c = L.congruence([[0, 1]])
            sage: I = L.quotient(c); I
            Finite lattice containing 2 elements
            sage: I.top()
            (2, 3, 4)
            sage: I = L.quotient(c, labels='lattice')
            sage: I.top()
            Finite lattice containing 3 elements

            sage: B3 = posets.BooleanLattice(3)
            sage: c = B3.congruence([[0,1]])
            sage: B2 = B3.quotient(c, labels='integer')
            sage: B2.is_isomorphic(posets.BooleanLattice(2))
            True

        .. SEEALSO:: :meth:`congruence`

        TESTS::

            sage: E = LatticePoset()
            sage: E.quotient([])
            Finite lattice containing 0 elements

            sage: L = posets.PentagonPoset()
            sage: L.quotient(L.congruence([[1]])).is_isomorphic(L)
            True
        """
        if labels not in ['lattice', 'tuple', 'integer']:
            raise ValueError("labels must be one of 'lattice', 'tuple' or 'integer'")

        parts_H = [sorted([self._element_to_vertex(e) for e in part]) for
                   part in congruence]
        minimal_vertices = [part[0] for part in parts_H]
        H = self._hasse_diagram.transitive_closure().subgraph(minimal_vertices).transitive_reduction()
        if labels == 'integer':
            H.relabel(list(range(len(minimal_vertices))))
            return LatticePoset(H)
        part_dict = {m[0]: [self._vertex_to_element(x) for x in m] for m
                     in parts_H}
        if labels == 'tuple':
            H.relabel(lambda m: tuple(part_dict[m]))
            return LatticePoset(H)
        H.relabel(lambda m: self.sublattice(part_dict[m]))
        return LatticePoset(H)

    def congruences_lattice(self, labels='congruence'):
        r"""
        Return the lattice of congruences.

        A congruence of a lattice is a partition of elements to classes
        compatible with both meet- and join-operation; see :meth:`congruence`.
        Elements of the *congruence lattice* are congruences ordered by
        refinement; i.e. if every class of a congruence `\Theta` is
        contained in some class of `\Phi`, then `\Theta \le \Phi`
        in the congruence lattice.

        INPUT:

        - ``labels`` -- a string; the type of elements in the resulting lattice

        OUTPUT:

        A distributive lattice.

        - If ``labels='congruence'``, then elements of the
          result will be congruences given as
          :class:`sage.combinat.set_partition.SetPartition`.
        - If ``labels='integers'``, result is a lattice on
          integers isomorphic to the congruence lattice.

        EXAMPLES::

            sage: N5 = posets.PentagonPoset()
            sage: CL = N5.congruences_lattice(); CL
            Finite lattice containing 5 elements
            sage: CL.atoms()
            [{{0}, {1}, {2, 3}, {4}}]
            sage: CL.coatoms()
            [{{0, 1}, {2, 3, 4}}, {{0, 2, 3}, {1, 4}}]

            sage: C4 = posets.ChainPoset(4)
            sage: CL = C4.congruences_lattice(labels='integer')
            sage: CL.is_isomorphic(posets.BooleanLattice(3))
            True

        TESTS::

            sage: posets.ChainPoset(0).congruences_lattice()
            Finite lattice containing 1 elements
            sage: posets.ChainPoset(1).congruences_lattice()
            Finite lattice containing 1 elements
            sage: posets.ChainPoset(2).congruences_lattice()
            Finite lattice containing 2 elements
            sage: posets.ChainPoset(3).congruences_lattice()
            Finite lattice containing 4 elements
        """
        from sage.sets.set import Set
        from sage.sets.disjoint_set import DisjointSet
        from sage.combinat.set_partition import SetPartition
        if labels not in ['integer', 'congruence']:
            raise ValueError("'labels' must be 'integer' or 'congruence'")

        cong_ji, congs = self._hasse_diagram.principal_congruences_poset()

        # Form of the lattice of congruences can be computed much faster than
        # all congruences.
        if labels == 'integer':
            tmp = Poset(cong_ji).order_ideals_lattice(as_ideals=False)
            return tmp.relabel(tmp._element_to_vertex_dict)

        # To compute full lattice of congruences we "extend" already computed
        # parts of a congruence.
        L = cong_ji.order_ideals_lattice()
        C = {}
        C[Set()] = DisjointSet(self.cardinality())  # the bottom element
        for e in L:
            low = L.lower_covers(e)
            if len(low) == 1:  # a join-irreducible element
                C[e] = congs[max(e, key=cong_ji._element_to_vertex)]
            if len(low) > 1:  # "extending" congruence to avoid re-computation
                low_0 = min(low, key=lambda x: C[x].number_of_subsets())
                for new_pair in e:
                    if new_pair not in low_0:
                        break
                C[e] = self._hasse_diagram.congruence([new_pair], start=C[low_0])

        return L.relabel(lambda e: SetPartition([[self._vertex_to_element(v)
                                                  for v in p] for p in C[e]]))


def _log_2(n):
    """
    Return the 2-based logarithm of `n` rounded up.

    `n` is assumed to be a positive integer.

    EXAMPLES::

        sage: sage.combinat.posets.lattices._log_2(10)
        4

    TESTS::

        sage: sage.combinat.posets.lattices._log_2(15)
        4
        sage: sage.combinat.posets.lattices._log_2(16)
        4
        sage: sage.combinat.posets.lattices._log_2(17)
        5
    """
    bits = -1
    i = n
    while i:
        i = i >> 1
        bits += 1
    if 1 << bits == n:
        return bits
    return bits + 1

############################################################################

FiniteMeetSemilattice._dual_class = FiniteJoinSemilattice
FiniteJoinSemilattice._dual_class = FiniteMeetSemilattice
FiniteLatticePoset   ._dual_class = FiniteLatticePoset
