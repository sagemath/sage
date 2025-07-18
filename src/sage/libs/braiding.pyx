# distutils: libraries = braiding
# distutils: language = c++
r"""
Cython wrapper for the libbraiding library.

The libbraiding library is a modification of the braiding program
by Juan Gonzalez-Meneses (https://github.com/jeanluct/cbraid)
to expose the functions as a C++ library instead of an interactive
program.

Braids are returned in left normal form as a list of lists. The
first list contains only an integer, representing the power of
`\Delta`. The subsequent lists are the Tietze lists of the elementary
permutation braids.
"""

# ****************************************************************************
#       Copyright (C) 2016 Miguel Marco  <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at youroption) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from cysignals.signals cimport sig_on, sig_off

from libcpp.list cimport list


cdef extern from "braiding.h" namespace "Braiding":
    list[list[int]] ConjugatingBraid (int n, list[int] word, list[int] word2)
    list[list[int]] LeftNormalForm (int n, list[int] word)
    list[list[int]] RightNormalForm (int n, list[int] word)
    list[list[int]] GreatestCommonDivisor(int n, list[int] word1, list[int] word2)
    list[list[int]] LeastCommonMultiple(int n, list[int] word1, list[int] word2)
    list[list[list[int]]] CentralizerGenerators(int n, list[int] word)
    list[list[list[int]]] SuperSummitSet(int n, list[int] word)
    list[list[list[list[int]]]] UltraSummitSet(int n, list[int] word)
    int thurstontype(int n, list[int] word)
    int Rigidity_ext(int n, list[int] word)
    list[list[list[list[int]]]] SlidingCircuits(int n, list[int] word)
    list[list[list[int]]] SendToSSS(int n, list[int] word)
    list[list[list[int]]] SendToUSS(int n, list[int] word)
    list[list[list[int]]] SendToSC(int n, list[int] word)
    list[list[list[int]]] Trajectory(int n, list[int] word)
    list[list[list[list[int]]]] CyclicSlidings(int n, list[int] word)


def conjugatingbraid(braid1, braid2):
    r"""
    Return a braid that conjugates ``braid1`` to ``braid2`` if
    such a braid exists.

    INPUT:

    - ``braid1`` -- the braid to be conjugated
    - ``braid2`` -- the braid to conjugate to

    OUTPUT:

    The list of lists that represent a conjugating braid. If the input braids
    are not conjugate, an empty list is returned.

    EXAMPLES::

        sage: from sage.libs.braiding import conjugatingbraid
        sage: B = BraidGroup(3)
        sage: b = B([1,2,1,-2])
        sage: c = B([1,2])
        sage: conjugatingbraid(b,c)
        [[0], [2]]
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = ConjugatingBraid(nstrands, l1, l2)
    sig_off()
    return rop


def leftnormalform(braid):
    r"""
    Return the left normal form of a braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of lists with the left normal form. The first list contains
    the power of `\Delta`. The subsequent lists are the elementary
    permutation braids.

    EXAMPLES::

        sage: from sage.libs.braiding import leftnormalform
        sage: B = BraidGroup(3)
        sage: b = B([1,2,1,-2])
        sage: leftnormalform(b)
        [[0], [2, 1]]
    """
    nstrands = braid.parent().strands()
    l1 = braid.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeftNormalForm(nstrands, l1)
    sig_off()
    return rop


def rightnormalform(braid):
    r"""
    Return the right normal form of a braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of lists with the right normal form. The first list contains
    the power of `\Delta`. The subsequent lists are the elementary
    permutation braids.

    EXAMPLES::

        sage: from sage.libs.braiding import rightnormalform
        sage: B = BraidGroup(3)
        sage: b = B([1,2,1,-2])
        sage: rightnormalform(b)
        [[2, 1], [0]]
    """
    nstrands = braid.parent().strands()
    l1 = braid.Tietze()
    sig_on()
    cdef list[list[int]] rop = RightNormalForm(nstrands, l1)
    sig_off()
    return rop


def greatestcommondivisor(braid1, braid2):
    r"""
    Return the greatest common divisor of two braids.

    INPUT:

    - ``braid1`` -- a braid
    - ``braid2`` -- a braid

    OUTPUT: list of lists representing the gcd of ``braid1`` and ``braid2``

    EXAMPLES::

        sage: from sage.libs.braiding import greatestcommondivisor
        sage: B = BraidGroup(3)
        sage: b1 = B([1, 2, -1])
        sage: b2 = B([2, 2, 2])
        sage: greatestcommondivisor(b1, b2)
        [[-1], [2, 1]]
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = GreatestCommonDivisor(nstrands, l1, l2)
    sig_off()
    return rop


def leastcommonmultiple(braid1, braid2):
    r"""
    Return the least common multiple of two braids.

    INPUT:

    - ``braid1`` -- a braid
    - ``braid2`` -- a braid

    OUTPUT: list of lists representing the lcm of ``braid1`` and ``braid2``

    EXAMPLES::

        sage: from sage.libs.braiding import leastcommonmultiple
        sage: B = BraidGroup(3)
        sage: b1 = B([1, 2, -1])
        sage: b2 = B([2, 2, 2])
        sage: leastcommonmultiple(b1, b2)
        [[1], [1], [1]]
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeastCommonMultiple(nstrands, l1, l2)
    sig_off()
    return rop


def centralizer(braid):
    r"""
    Return a list of generators of the centralizer of a braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of lists representing the generators of the centralizer
    of ``braid``.

    EXAMPLES::

        sage: from sage.libs.braiding import centralizer
        sage: B = BraidGroup(3)
        sage: b = B([1,2,-1])
        sage: centralizer(b)
        [[[-1], [2, 1], [1, 2]], [[0], [1], [1, 2], [2]]]
    """
    nstrands = braid.parent().strands()
    lnf = leftnormalform(braid)
    if len(lnf) == 1: # (lib)braiding crashes when the input is a power of Delta.
        if lnf[0][0] % 2 == 0:
            return [[[0], [i+1]] for i in range(nstrands)]
        elif nstrands % 2:
            return [[[0], [i+1, nstrands - i -1]] for i in range(nstrands//2)]
        else:
            return [[[0], [i+1, nstrands - i -1]] for i in range(nstrands//2-1)] + [[[0], [nstrands//2]]]
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = CentralizerGenerators(nstrands, l)
    sig_off()
    return rop


def supersummitset(braid):
    r"""
    Return a list with the super-summit-set of a braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT: list of lists representing the super summit set of ``braid``

    EXAMPLES::

        sage: from sage.libs.braiding import supersummitset
        sage: B = BraidGroup(3)
        sage: b = B([1,2,-1])
        sage: supersummitset(b)
        [[[0], [2]], [[0], [1]]]
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SuperSummitSet(nstrands, l)
    sig_off()
    return rop


def ultrasummitset(braid):
    r"""
    Return a list with the orbits forming the ultra-summit-set of the braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of lists of lists representing the orbits of the ultra summit
    set of ``braid``.

    EXAMPLES::

        sage: from sage.libs.braiding import ultrasummitset
        sage: B = BraidGroup(3)
        sage: b = B([1,2,-1])
        sage: ultrasummitset(b)
        [[[[0], [2]]], [[[0], [1]]]]
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[list[int]]]] rop = UltraSummitSet(nstrands, l)
    sig_off()
    return rop


def thurston_type(braid):
    r"""
    Return the Thurston type of the braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    One of ``'periodic'``, ``'reducible'`` or ``'pseudo-anosov'``.

    EXAMPLES::

        sage: from sage.libs.braiding import thurston_type
        sage: B = BraidGroup(3)
        sage: b = B([1,2,-1])
        sage: thurston_type(b)
        'reducible'
        sage: c = B([1,2,1])
        sage: thurston_type(c)
        'periodic'
        sage: d = B([1,1,1,2,2])
        sage: thurston_type(d)
        'pseudo-anosov'
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef int i = thurstontype(nstrands, l)
    sig_off()
    if i == 1:
        return 'periodic'
    if i == 2:
        return 'reducible'
    if i == 3:
        return 'pseudo-anosov'


def rigidity(braid):
    r"""
    Return the rigidity of the braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT: the rigidity of the braid

    EXAMPLES::

        sage: from sage.libs.braiding import rigidity
        sage: B = BraidGroup(3)
        sage: c = B([1,1,1,2,2])
        sage: rigidity(c)
        3
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef int i = Rigidity_ext(nstrands, l)
    sig_off()
    return i


def sliding_circuits(braid):
    r"""
    Return the set of sliding circuits of the braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list with the sliding circuits of ``braid``. Each sliding circuit
    is a list of braids.

    EXAMPLES::

        sage: from sage.libs.braiding import sliding_circuits
        sage: B = BraidGroup(3)
        sage: c = B([1,1,1,2,2])
        sage: sliding_circuits(c)
        [[[[0], [1], [1, 2], [2, 1]]],
        [[[0], [2], [2, 1], [1, 2]]],
        [[[0], [1, 2], [2, 1], [1]]],
        [[[0], [2, 1], [1, 2], [2]]],
        [[[0], [1, 2], [2], [2, 1]]],
        [[[0], [2, 1], [1], [1, 2]]]]
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[list[int]]]] rop = SlidingCircuits(nstrands, l)
    sig_off()
    return rop


def send_to_sss(braid):
    r"""
    Return an element of the braid's super summit set and the conjugating braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list with two braids, the first one is an element of ``braid`` super summit
    set, the second one is the corresponding conjugating braid.

    EXAMPLES::

        sage: from sage.libs.braiding import send_to_sss
        sage: B = BraidGroup(4)
        sage: d = B([1, 2, 1, 2, 3, -1, 2,- 3])
        sage: send_to_sss(d)
        [[[0], [1, 2, 1, 3]], [[-1], [1, 2, 3, 2]]]

    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SendToSSS(nstrands, l)
    sig_off()
    return rop


def send_to_uss(braid):
    r"""
    Return an element of the braid's ultra summit set and the conjugating braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list with two braids, the first one is an element of ``braid`` ultra summit
    set, the second one is the corresponding conjugating braid.

    EXAMPLES::

        sage: from sage.libs.braiding import send_to_uss
        sage: B = BraidGroup(4)
        sage: d = B([1, 2, 1, 2, 3, -1, 2,- 1])
        sage: send_to_uss(d)
        [[[0], [1, 2, 3, 2]], [[0], [1]]]

    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SendToUSS(nstrands, l)
    sig_off()
    return rop


def send_to_sc(braid):
    r"""
    Return an element of the braid's sliding circuits and the conjugating braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list with two braids, the first one is an element of ``braid`` sliding
    circuits, the second one is the corresponding conjugating braid.

    EXAMPLES::

        sage: from sage.libs.braiding import send_to_sc
        sage: B = BraidGroup(4)
        sage: d = B([1, 2, 1, 2, 3, -1, 2, 2])
        sage: send_to_sc(d)
        [[[0], [1, 2, 3, 2], [2, 1]], [[0], [1]]]


    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SendToSC(nstrands, l)
    sig_off()
    return rop


def trajectory(braid):
    r"""
    Return the braid's trajectory.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of braids, formed by the ``braid```trajectory.

    EXAMPLES::

        sage: from sage.libs.braiding import trajectory
        sage: B = BraidGroup(4)
        sage: d = B([1, 2, 1, 2, 3, -1, 2, 2])
        sage: trajectory(d)
        [[[0], [1, 3], [1, 2, 3, 2]],
         [[0], [1, 2, 3, 2, 1], [3]],
         [[0], [1, 2, 3, 2], [2, 1]],
         [[0], [2, 1, 3], [1, 2, 3]]]


    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = Trajectory(nstrands, l)
    sig_off()
    return rop


def cyclic_slidings(braid):
    r"""
    Return the cyclic slidings of the braid.

    INPUT:

    - ``braid`` -- a braid

    OUTPUT:

    A list of lists of braids, given by the input's cyclic slidings.

    EXAMPLES::

        sage: from sage.libs.braiding import cyclic_slidings
        sage: B = BraidGroup(4)
        sage: d = B([1, 2, 1, 2, 3, -1, 2, 2])
        sage: cyclic_slidings(d)
        [[[[0], [1, 2, 3, 2], [2, 1]],
          [[0], [1, 2, 3, 2, 1], [3]],
          [[0], [2, 1, 3], [1, 2, 3]]],
         [[[0], [1, 3, 2, 1], [2, 3]],
          [[0], [1, 2, 3, 2, 1], [1]],
          [[0], [2, 1, 3], [3, 2, 1]]]]

    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[list[int]]]] rop = CyclicSlidings(nstrands, l)
    sig_off()
    return rop
