# sage.doctest: needs sage.combinat sage.modules
r"""
PBW data

This contains helper classes and functions which encode PBW data
in finite type.

AUTHORS:

- Dinakar Muthiah (2015-05): initial version
- Travis Scrimshaw (2016-06): simplified code and converted to Cython
"""

# ****************************************************************************
#       Copyright (C) 2015 Dinakar Muthiah <muthiah at ualberta.ca>
#                          Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator

cimport cython


class PBWDatum():
    """
    Helper class which represents a PBW datum.
    """
    def __init__(self, parent, long_word, lusztig_datum):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L = PBWDatum(P, (1,2,1), (1,4,7))
            sage: TestSuite(L).run(skip='_test_pickling')
        """
        self.parent = parent
        self.long_word = tuple(long_word)
        self.lusztig_datum = tuple(lusztig_datum)

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: PBWDatum(P, (1,2,1), (1,4,7))
            PBW Datum element of type ['A', 2] with long word (1, 2, 1)
             and Lusztig datum (1, 4, 7)
        """
        return_str = "PBW Datum element of type {cartan_type} with ".format(
                     cartan_type=self.parent.cartan_type)
        return_str += "long word {long_word} and Lusztig datum {lusztig_datum}".format(
                      long_word=self.long_word,
                      lusztig_datum=self.lusztig_datum)
        return return_str

    def __eq__(self, other_PBWDatum):
        """
        Check equality.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,4,7))
            sage: L2 = PBWDatum(P, (1,2,1), (1,4,7))
            sage: L1 == L2
            True
        """
        return (self.parent == other_PBWDatum.parent and
                self.long_word == other_PBWDatum.long_word and
                self.lusztig_datum == other_PBWDatum.lusztig_datum)

    def is_equivalent_to(self, other_pbw_datum):
        r"""
        Return whether ``self`` is equivalent to ``other_pbw_datum``.
        modulo the tropical Plücker relations.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,0,1))
            sage: L2 = PBWDatum(P, (2,1,2), (0,1,0))
            sage: L1.is_equivalent_to(L2)
            True
            sage: L1 == L2
            False
        """
        other_long_word = other_pbw_datum.long_word
        other_lusztig_datum = other_pbw_datum.lusztig_datum
        equiv_pbw_datum = self.convert_to_new_long_word(other_long_word)
        return equiv_pbw_datum.lusztig_datum == other_lusztig_datum

    def convert_to_long_word_with_first_letter(self, i):
        r"""
        Return a new PBWDatum equivalent to ``self``
        whose long word begins with ``i``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A3")
            sage: datum = PBWDatum(P, (1,2,1,3,2,1), (1,0,1,4,2,3))
            sage: datum.convert_to_long_word_with_first_letter(1)
            PBW Datum element of type ['A', 3] with long word (1, 2, 3, 1, 2, 1)
             and Lusztig datum (1, 0, 4, 1, 2, 3)
            sage: datum.convert_to_long_word_with_first_letter(2)
            PBW Datum element of type ['A', 3] with long word (2, 1, 2, 3, 2, 1)
             and Lusztig datum (0, 1, 0, 4, 2, 3)
            sage: datum.convert_to_long_word_with_first_letter(3)
            PBW Datum element of type ['A', 3] with long word (3, 1, 2, 3, 1, 2)
             and Lusztig datum (8, 1, 0, 4, 1, 2)
        """
        return self.convert_to_new_long_word(self.parent._long_word_begin_with(i))

    def convert_to_new_long_word(self, new_long_word):
        r"""
        Return a new PBWDatum equivalent to ``self``
        whose long word is ``new_long_word``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: datum = PBWDatum(P, (1,2,1), (1,0,1))
            sage: new_datum = datum.convert_to_new_long_word((2,1,2))
            sage: new_datum.long_word
            (2, 1, 2)
            sage: new_datum.lusztig_datum
            (0, 1, 0)
        """
        return self.parent.convert_to_new_long_word(self, new_long_word)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L = PBWDatum(P, (1,2,1), (1,1,1))
            sage: L.weight()
            -2*alpha[1] - 2*alpha[2]
        """
        root_list = self.parent._root_list_from(tuple(self.long_word))
        R = self.parent.root_lattice
        return R.linear_combination((root_list[i], -coeff)
                                    for i, coeff in enumerate(self.lusztig_datum))

    def star(self):
        """
        Return the starred version of ``self``, i.e.,
        with reversed ``long_word`` and ``lusztig_datum``

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: L1 = PBWDatum(P, (1,2,1), (1,2,3))
            sage: L1.star() == PBWDatum(P, (2,1,2), (3,2,1))
            True
        """
        aut = self.parent.cartan_type.opposition_automorphism()
        reversed_long_word = [aut[i] for i in reversed(self.long_word)]
        reversed_lusztig_datum = reversed(self.lusztig_datum)
        return PBWDatum(self.parent, reversed_long_word, reversed_lusztig_datum)


class PBWData(): # UniqueRepresentation?
    """
    Helper class for the set of PBW data.
    """
    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData
            sage: P = PBWData(["A",2])
            sage: TestSuite(P).run(skip='_test_pickling')
        """
        self.cartan_type = CartanType(cartan_type)
        self.root_system = RootSystem(self.cartan_type)
        self.root_lattice = self.root_system.root_lattice()
        self.weyl_group = self.root_lattice.weyl_group()
        self._braid_move_calc = BraidMoveCalculator(self.weyl_group)

    def convert_to_new_long_word(self, pbw_datum, new_long_word):
        """
        Convert the PBW datum ``pbw_datum`` from its long word to
        ``new_long_word``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum
            sage: P = PBWData("A2")
            sage: datum = PBWDatum(P, (1,2,1), (1,0,1))
            sage: new_datum = P.convert_to_new_long_word(datum,(2,1,2))
            sage: new_datum
            PBW Datum element of type ['A', 2] with long word (2, 1, 2)
             and Lusztig datum (0, 1, 0)
            sage: new_datum.long_word
            (2, 1, 2)
            sage: new_datum.lusztig_datum
            (0, 1, 0)
        """
        assert pbw_datum.parent is self
        chain = self._braid_move_calc.chain_of_reduced_words(pbw_datum.long_word,
                                                             new_long_word)
        cdef list enhanced_braid_chain = enhance_braid_move_chain(chain, self.cartan_type)
        new_lusztig_datum = compute_new_lusztig_datum(enhanced_braid_chain,
                                                      pbw_datum.lusztig_datum)
        return PBWDatum(self, new_long_word, new_lusztig_datum)

    @cached_method
    def _root_list_from(self, reduced_word):
        """
        Return the list of positive roots in the order determined by
        ``reduced_word``.

        .. WARNING::

            No error checking is done to verify that ``reduced_word``
            is reduced.

        INPUT:

        - ``reduced_word`` -- tuple corresponding to a reduced word

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData
            sage: P = PBWData(["A",2])
            sage: P._root_list_from((1,2,1))
            [alpha[1], alpha[1] + alpha[2], alpha[2]]
        """
        al = self.root_lattice.simple_roots()
        cur = []
        for i in reversed(reduced_word):
            cur = [al[i]] + [x.simple_reflection(i) for x in cur]
        return cur

    @cached_method
    def _long_word_begin_with(self, i):
        """
        Return a reduced expression of the long word which begins with ``i``.

        EXAMPLES::

            sage: from sage.combinat.crystals.pbw_datum import PBWData
            sage: P = PBWData(["C",3])
            sage: P._long_word_begin_with(1)
            (1, 3, 2, 3, 1, 2, 3, 1, 2)
            sage: P._long_word_begin_with(2)
            (2, 3, 2, 3, 1, 2, 3, 2, 1)
            sage: P._long_word_begin_with(3)
            (3, 2, 3, 1, 2, 3, 1, 2, 1)
        """
        si = self.weyl_group.simple_reflection(i)
        w0 = self.weyl_group.long_element()
        return tuple([i] + (si * w0).reduced_word())


# enhanced_braid_chain is an ugly data structure.
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple compute_new_lusztig_datum(list enhanced_braid_chain, initial_lusztig_datum):
    """
    Return the Lusztig datum obtained by applying tropical Plücker
    relations along ``enhanced_braid_chain`` starting with
    ``initial_lusztig_datum``.

    EXAMPLES::

        sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: from sage.combinat.crystals.pbw_datum import compute_new_lusztig_datum
        sage: ct = CartanType(['A', 2])
        sage: W = CoxeterGroup(ct)
        sage: B = BraidMoveCalculator(W)
        sage: chain = B.chain_of_reduced_words((1,2,1),(2,1,2))
        sage: enhanced_braid_chain = enhance_braid_move_chain(chain, ct)
        sage: compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1))
        (0, 1, 0)

    TESTS::

        sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: from sage.combinat.crystals.pbw_datum import compute_new_lusztig_datum
        sage: ct = CartanType(['A', 2])
        sage: W = CoxeterGroup(ct)
        sage: B = BraidMoveCalculator(W)
        sage: chain = B.chain_of_reduced_words((1,2,1), (2,1,2))
        sage: enhanced_braid_chain = enhance_braid_move_chain(chain, ct)
        sage: compute_new_lusztig_datum(enhanced_braid_chain,(1,0,1)) == (0,1,0)
        True
    """
    cdef tuple interval_of_change
    # Does not currently check that len(initial_lusztig_datum) is appropriate
    cdef list new_lusztig_datum = list(initial_lusztig_datum)  # shallow copy
    cdef int i
    for i in range(1, len(enhanced_braid_chain)):
        interval_of_change, type_data = enhanced_braid_chain[i]
        a, b = interval_of_change
        old_interval_datum = new_lusztig_datum[a:b]
        new_interval_datum = tropical_plucker_relation(type_data, old_interval_datum)
        new_lusztig_datum[a:b] = new_interval_datum
    return tuple(new_lusztig_datum)


# The tropical Plücker relations
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple tropical_plucker_relation(tuple a, lusztig_datum):
    r"""
    Apply the tropical Plücker relation of type ``a`` to ``lusztig_datum``.

    The relations are obtained by tropicalizing the relations in
    Proposition 7.1 of [BZ01]_.

    INPUT:

    - ``a`` -- a pair ``(x, y)`` of the off-diagonal entries of a
      `2 \times 2` Cartan matrix

    EXAMPLES::

        sage: from sage.combinat.crystals.pbw_datum import tropical_plucker_relation
        sage: tropical_plucker_relation((0,0), (2,3))
        (3, 2)
        sage: tropical_plucker_relation((-1,-1), (1,2,3))
        (4, 1, 2)
        sage: tropical_plucker_relation((-1,-2), (1,2,3,4))
        (8, 1, 2, 3)
        sage: tropical_plucker_relation((-2,-1), (1,2,3,4))
        (6, 1, 2, 3)
    """
    if a == (0, 0): # A1xA1
        t1, t2 = lusztig_datum
        return (t2, t1)
    elif a == (-1, -1): # A2
        t1,t2,t3 = lusztig_datum
        return (t2+t3-min(t1,t3),
                min(t1,t3),
                t1+t2-min(t1,t3))
    elif a == (-1, -2): # B2
        t1,t2,t3,t4 = lusztig_datum
        pi1 = min(t1+t2,min(t1,t3)+t4)
        pi2 = min(2*t1+t2,2*min(t1,t3)+t4)
        return (t2+2*t3+t4-pi2,
                pi2-pi1,
                2*pi1-pi2,
                t1+t2+t3-pi1)
    elif a == (-1, -3): # G2
        t1,t2,t3,t4,t5,t6 = lusztig_datum
        pi1 = min(t1+t2+2*t3+t4,
                  t1+t2+2*min(t3,t5)+t6,
                  min(t1,t3)+t4+2*t5+t6)
        pi2 = min(2*t1+2*t2+3*t3+t4,
                  2*t1+2*t2+3*min(t3,t5)+t6,
                  2*min(t1,t3)+2*t4+3*t5+t6,
                  t1+t2+t4+2*t5+t6+min(t1+t3,2*t3,t3+t5,t1+t5))
        pi3 = min(3*t1+2*t2+3*t3+t4,
                  3*t1+2*t2+3*min(t3,t5)+t6,
                  3*min(t1,t3)+2*t4+3*t5+t6,
                  2*t1+t2+t4+2*t5+t6+min(t1+t3,2*t3,t3+t5,t1+t5))
        pi4 = min(2*t1+2*t2+3*t3+t4+min(t1+t2+3*t3+t4,
                                        t1+t2+3*min(t3,t5)+t6,
                                        min(t1+t3,2*t3,t3+t5,t1+t5)+t4+2*t5+t6),
                  2*t6+3*min(t1+t2+2*min(t3,t5),min(t1,t3)+t4+2*t5))
        return (t2+3*t3+2*t4+3*t5+t6-pi3,
                pi3-pi2,
                3*pi2-pi3-pi4,
                pi4-pi1-pi2,
                3*pi1-pi4,
                t1+t2+2*t3+t4+t5-pi1)
    else: # (-1,-2) and (-1,-3)
        reversed_lusztig_datum = tuple(reversed(lusztig_datum))
        return tuple(reversed(tropical_plucker_relation((a[1], a[0]),
                                                        reversed_lusztig_datum)))


# Maybe we need to be more specific, and pass not the Cartan type, but the root lattice?
# TODO: Move to PBW_data?
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list enhance_braid_move_chain(braid_move_chain, cartan_type):
    r"""
    Return a list of tuples that records the data of the long words in
    ``braid_move_chain`` plus the data of the intervals where the braid moves
    occur and the data of the off-diagonal entries of the `2 \times 2` Cartan
    submatrices of each braid move.

    INPUT:

    - ``braid_move_chain`` -- a chain of reduced words in the Weyl group
      of ``cartan_type``
    - ``cartan_type`` -- a finite Cartan type

    OUTPUT:

    A list of 2-tuples
    ``(interval_of_change, cartan_sub_matrix)`` where

    - ``interval_of_change`` is the (half-open) interval of indices where
      the braid move occurs; this is ``None`` for the first tuple
    - ``cartan_sub_matrix`` is the off-diagonal entries of the `2 \times 2`
      submatrix of the Cartan matrix corresponding to the braid move;
      this is ``None`` for the first tuple

    For a matrix::

        [2 a]
        [b 2]

    the ``cartan_sub_matrix`` is the pair ``(a, b)``.

    TESTS::

        sage: from sage.combinat.crystals.pbw_datum import enhance_braid_move_chain
        sage: braid_chain = [(1, 2, 1, 3, 2, 1),
        ....:                (1, 2, 3, 1, 2, 1),
        ....:                (1, 2, 3, 2, 1, 2),
        ....:                (1, 3, 2, 3, 1, 2),
        ....:                (3, 1, 2, 3, 1, 2),
        ....:                (3, 1, 2, 1, 3, 2),
        ....:                (3, 2, 1, 2, 3, 2),
        ....:                (3, 2, 1, 3, 2, 3)]
        sage: enhanced_chain = enhance_braid_move_chain(braid_chain, CartanType(["A",5]))
        sage: enhanced_chain[0]
        (None, None)
        sage: enhanced_chain[7]
        ((3, 6), (-1, -1))
    """
    cdef int i, j
    cdef int k, pos, first, last
    cdef tuple cartan_sub_matrix
    cdef list output_list = []
    output_list.append( (None, None) )
    cdef tuple previous_word = <tuple> (braid_move_chain[0])
    cdef tuple current_word
    cartan_matrix = cartan_type.cartan_matrix()
    cdef int ell = len(previous_word)
    # TODO - Optimize this by avoiding calls to here?
    # This likely could be done when performing chain_of_reduced_words
    # Things in here get called the most (about 50x more than enhance_braid_move_chain)
    for pos in range(1, len(braid_move_chain)):
        # This gets the smallest contiguous half-open interval [a, b)
        # that contains the indices where current_word and previous_word differ.
        current_word = <tuple> (braid_move_chain[pos])
        for k in range(ell):
            i = previous_word[k]
            j = current_word[k]
            if i != j:
                i -= 1  # -1 for indexing
                j -= 1  # -1 for indexing
                first = k
                break
        for k in range(ell-1, k-1, -1):
            if previous_word[k] != current_word[k]:
                last = k + 1
                break

        cartan_sub_matrix = (cartan_matrix[i, j], cartan_matrix[j, i])
        output_list.append(((first, last), cartan_sub_matrix))
        previous_word = current_word
    return output_list
