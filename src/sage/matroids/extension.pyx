"""
Iterators for linear subclasses

The classes below are iterators returned by the functions
:func:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`
and :func:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`.
See the documentation of these methods for more detail.
For direct access to these classes, run::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version

Methods
=======
"""
# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.data_structures.bitset_base cimport *
from sage.matroids.basis_matroid cimport BasisMatroid


cdef class CutNode:
    """
    An internal class used for creating linear subclasses of a matroids in a
    depth-first manner.

    A linear subclass is a set of hyperplanes `mc` with the property that
    certain sets of hyperplanes must either be fully contained in `mc` or
    intersect `mc` in at most 1 element. The way we generate them is by a
    depth-first search. This class represents a node in the search tree.

    It contains the set of hyperplanes selected so far, as well as a
    collection of hyperplanes whose insertion has been explored elsewhere in
    the search tree.

    The class has methods for selecting a hyperplane to insert, for inserting
    hyperplanes and closing the set to become a linear subclass again, and for
    adding a hyperplane to the set of *forbidden* hyperplanes, and similarly
    closing that set.
    """
    def __cinit__(self, MC, N=None):
        """
        Internal data structure init.

        EXAMPLES::

            sage: len(list(matroids.catalog.Fano().linear_subclasses()))  # indirect doctest
            16
        """
        cdef CutNode node
        self._MC = MC
        bitset_init(self._p_free, self._MC._hyperplanes_count + 1)
        bitset_init(self._p_in, self._MC._hyperplanes_count + 1)
        bitset_init(self._l0, self._MC._hyperlines_count + 1)
        bitset_init(self._l1, self._MC._hyperlines_count + 1)
        if N is None:
            bitset_set_first_n(self._p_free, self._MC._hyperplanes_count)
            bitset_clear(self._p_in)
            bitset_set_first_n(self._l0, self._MC._hyperlines_count)
            bitset_clear(self._l1)
        else:
            node = N
            bitset_copy(self._p_free, node._p_free)
            bitset_copy(self._p_in, node._p_in)
            bitset_copy(self._l0, node._l0)
            bitset_copy(self._l1, node._l1)
            self._ml = node._ml

    def __dealloc__(self):
        bitset_free(self._p_free)
        bitset_free(self._p_in)
        bitset_free(self._l0)
        bitset_free(self._l1)

    cdef CutNode copy(self):
        return CutNode(self._MC, self)

    cdef bint insert_plane(self, long p0) noexcept:
        """
        Add a hyperplane to the linear subclass.
        """
        cdef long l, p
        cdef list p_stack, l_stack
        if bitset_in(self._p_in, p0):
            return True
        if not bitset_in(self._p_free, p0):
            return False
        bitset_discard(self._p_free, p0)
        bitset_add(self._p_in, p0)
        p_stack = [p0]
        l_stack = []
        while len(p_stack) > 0:
            while len(p_stack) > 0:
                p = p_stack.pop()
                for l in self._MC._lines_on_plane[p]:
                    if bitset_in(self._l0, l):
                        bitset_discard(self._l0, l)
                        bitset_add(self._l1, l)
                    else:
                        if bitset_in(self._l1, l):
                            bitset_discard(self._l1, l)
                            l_stack.append(l)
            while len(l_stack) > 0:
                l = l_stack.pop()
                for p in self._MC._planes_on_line[l]:
                    if bitset_in(self._p_in, p):
                        continue
                    if bitset_in(self._p_free, p):
                        bitset_discard(self._p_free, p)
                        bitset_add(self._p_in, p)
                        p_stack.append(p)
                    else:
                        return False
        return True

    cdef bint remove_plane(self, long p0) noexcept:
        """
        Remove a hyperplane from the linear subclass.
        """
        cdef long p, l
        cdef list p_stack, l_stack
        bitset_discard(self._p_free, p0)
        p_stack = [p0]
        l_stack = []
        while len(p_stack) > 0:
            while len(p_stack) > 0:
                p = p_stack.pop()
                for l in self._MC._lines_on_plane[p]:
                    if bitset_in(self._l1, l):
                        l_stack.append(l)
            while len(l_stack) > 0:
                l = l_stack.pop()
                for p in self._MC._planes_on_line[l]:
                    if bitset_in(self._p_free, p):
                        bitset_discard(self._p_free, p)
                        p_stack.append(p)
                    else:
                        return False
        return True

    cdef select_plane(self):
        """
        Choose a hyperplane from the linear subclass.
        """
        cdef long l, p
        while self._ml >= 0:
            l = self._MC._mandatory_lines[self._ml]
            if bitset_in(self._l0, l):
                for p in self._MC._planes_on_line[l]:
                    if bitset_in(self._p_free, p):
                        return p
                return None
            self._ml -= 1

        return bitset_first(self._p_free)

    cdef list planes(self):
        """
        Return all hyperplanes from the linear subclass.
        """
        return bitset_list(self._p_in)

cdef class LinearSubclassesIter:
    """
    An iterator for a set of linear subclass.
    """
    def __init__(self, MC):
        """
        Create a linear subclass iterator.

        Auxiliary to class LinearSubclasses.

        INPUT:

        - ``MC`` -- a member of class LinearSubclasses

        EXAMPLES::

            sage: from sage.matroids.extension import LinearSubclasses
            sage: M = matroids.Uniform(3, 6)
            sage: type(LinearSubclasses(M).__iter__())
            <... 'sage.matroids.extension.LinearSubclassesIter'>
        """
        cdef CutNode first_cut = CutNode(MC)
        self._MC = MC
        self._nodes = []

        for i in self._MC._forbidden_planes:
            if not first_cut.remove_plane(i):
                return

        for i in self._MC._mandatory_planes:
            if not first_cut.insert_plane(i):
                return

        first_cut._ml = len(self._MC._mandatory_lines) - 1

        self._nodes = [first_cut]

    def __next__(self):
        """
        Return the next linear subclass.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: from sage.matroids.extension import LinearSubclasses
            sage: M = BasisMatroid(matroids.Uniform(3, 6))
            sage: I = LinearSubclasses(M).__iter__()
            sage: M.extension('x', I.__next__())
            Matroid of rank 3 on 7 elements with 35 bases
            sage: M.extension('x', I.__next__())
            Matroid of rank 3 on 7 elements with 34 bases
        """
        cdef CutNode node, node2
        while True:
            if len(self._nodes) == 0:
                raise StopIteration
            node = self._nodes.pop()
            p0 = node.select_plane()
            while p0 is not None and p0 >= 0:
                node2 = node.copy()
                if node2.insert_plane(p0):
                    self._nodes.append(node2)
                node.remove_plane(p0)
                p0 = node.select_plane()

            if p0 is not None:
                res = self._MC[node]
                if res is not None:
                    return res


cdef class LinearSubclasses:
    r"""
    An iterable set of linear subclasses of a matroid.

    Enumerate linear subclasses of a given matroid. A *linear subclass* is a
    collection of hyperplanes (flats of rank `r - 1` where `r` is the rank of
    the matroid) with the property that no modular triple of hyperplanes has
    exactly two members in the linear subclass. A triple of hyperplanes in a
    matroid of rank `r` is *modular* if its intersection has rank `r - 2`.

    INPUT:

    - ``M`` -- matroid
    - ``line_length`` -- integer (default: ``None``)
    - ``subsets`` -- (default: ``None``) a set of subsets of the groundset of
      ``M``
    - ``splice`` -- (default: ``None``) a matroid `N` such that for some
      `e \in E(N)` and some `f \in E(M)`, we have
      `N\setminus e= M\setminus f`

    OUTPUT: an enumerator for the linear subclasses of M

    If ``line_length`` is not ``None``, the enumeration is restricted to
    linear subclasses ``mc`` so containing at least one of each set of
    ``line_length`` hyperplanes which have a common intersection of
    rank `r - 2`.

    If ``subsets`` is not ``None``, the enumeration is restricted to linear
    subclasses ``mc`` containing all hyperplanes which fully contain some set
    from ``subsets``.

    If ``splice`` is not ``None``, then the enumeration is restricted to
    linear subclasses `mc` such that if `M'` is the extension of `M` by `e`
    that arises from `mc`, then `M'\setminus f = N`.

    EXAMPLES::

        sage: from sage.matroids.extension import LinearSubclasses
        sage: M = matroids.Uniform(3, 6)
        sage: len([mc for mc in LinearSubclasses(M)])
        83
        sage: len([mc for mc in LinearSubclasses(M, line_length=5)])
        22
        sage: for mc in LinearSubclasses(M, subsets=[[0, 1], [2, 3], [4, 5]]):
        ....:     print(len(mc))
        3
        15

    Note that this class is intended for runtime, internal use, so no
    loads/dumps mechanism was implemented.
    """
    def __init__(self, M, line_length=None, subsets=None, splice=None):
        """
        See class docstring for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *  # LinearSubclasses, BasisMatroid
            sage: M = matroids.Uniform(3, 6)
            sage: len([mc for mc in LinearSubclasses(M)])
            83
            sage: len([mc for mc in LinearSubclasses(M, line_length=5)])
            22
            sage: for mc in LinearSubclasses(M,
            ....:                           subsets=[[0, 1], [2, 3], [4, 5]]):
            ....:     print(len(mc))
            3
            15
            sage: M = BasisMatroid(matroids.catalog.BetsyRoss()); M
            Matroid of rank 3 on 11 elements with 140 bases
            sage: e = 'k'; f = 'h'; Me = M.delete(e); Mf=M.delete(f)
            sage: for mc in LinearSubclasses(Mf, splice=Me):
            ....:     print(Mf.extension(f, mc))
            Matroid of rank 3 on 11 elements with 141 bases
            Matroid of rank 3 on 11 elements with 140 bases
            sage: for mc in LinearSubclasses(Me, splice=Mf):
            ....:     print(Me.extension(e, mc))
            Matroid of rank 3 on 11 elements with 141 bases
            Matroid of rank 3 on 11 elements with 140 bases
        """
        # set up hyperplanes/ hyperlines
        E = M.groundset()
        R = M.full_rank()

        self._hyperlines = M.flats(R - 2)
        self._hyperlines_count = len(self._hyperlines)

        self._hyperplanes = M.flats(R - 1)
        self._hyperplanes_count = len(self._hyperplanes)

        self._planes_on_line = [[] for l in range(self._hyperlines_count)]
        self._lines_on_plane = [[] for l in range(self._hyperplanes_count)]
        for l in range(self._hyperlines_count):
            for h in range(self._hyperplanes_count):
                if self._hyperplanes[h] >= self._hyperlines[l]:
                    self._lines_on_plane[h].append(l)
                    self._planes_on_line[l].append(h)

        self._mandatory_planes = []
        self._forbidden_planes = []
        self._mandatory_lines = []
        self._line_length = -1

        if line_length is not None:
            self._line_length = line_length
            self._mandatory_lines = [l for l in range(self._hyperlines_count) if len(self._planes_on_line[l]) >= line_length]

        if subsets is not None:
            for p in range(self._hyperplanes_count):
                H = self._hyperplanes[p]
                for S in subsets:
                    if frozenset(S).issubset(H):
                        self._mandatory_planes.append(p)
                        break

        if splice is not None:
            E2 = splice.groundset()
            F = frozenset(E - E2)
            F2 = frozenset(E2 - E)
            if len(E) != len(E2) or len(F) != 1:
                raise ValueError("LinearSubclasses: the ground set of the splice matroid is not of the form E + e-f")

            for p in range(self._hyperplanes_count):
                X = self._hyperplanes[p] - F
                if splice._rank(X) == splice.full_rank() - 1:
                    if splice._rank(X | F2) == splice.full_rank() - 1:
                        self._mandatory_planes.append(p)
                    else:
                        self._forbidden_planes.append(p)

    def __iter__(self):
        """
        Return an iterator for the linear subclasses.

        EXAMPLES::

            sage: from sage.matroids.extension import LinearSubclasses
            sage: M = matroids.Uniform(3, 6)
            sage: for mc in LinearSubclasses(M, subsets=[[0, 1], [2, 3], [4, 5]]):
            ....:     print(len(mc))
            3
            15
        """
        return LinearSubclassesIter(self)

    def __getitem__(self, CutNode node):
        """
        Return a linear subclass stored in a given CutNode.

        Internal function.

        EXAMPLES::

            sage: from sage.matroids.extension import LinearSubclasses
            sage: M = matroids.Uniform(3, 6)
            sage: len([mc for mc in LinearSubclasses(M)])
            83
        """
        cdef long p
        return [self._hyperplanes[p] for p in node.planes()]


cdef class MatroidExtensions(LinearSubclasses):
    r"""
    An iterable set of single-element extensions of a given matroid.

    INPUT:

    - ``M`` -- matroid
    - ``e`` -- an element
    - ``line_length`` -- integer (default: ``None``)
    - ``subsets`` -- (default: ``None``) a set of subsets of the groundset of
      ``M``
    - ``splice`` -- a matroid `N` such that for some `f \in E(M)`, we have
      `N\setminus e= M\setminus f`

    OUTPUT:

    An enumerator for the extensions of ``M`` to a matroid ``N`` so that
    `N\setminus e = M`. If ``line_length`` is not ``None``, the enumeration
    is restricted to extensions `N` without `U(2, k)`-minors, where
    ``k > line_length``.

    If ``subsets`` is not ``None``, the enumeration is restricted to
    extensions `N` of `M` by element `e` so that all hyperplanes of `M`
    which fully contain some set from ``subsets``, will also span `e`.

    If ``splice`` is not ``None``, then the enumeration is restricted to
    extensions `M'` such that `M'\setminus f = N`, where
    `E(M)\setminus E(N)=\{f\}`.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.Uniform(3, 6)
        sage: len([N for N in MatroidExtensions(M, 'x')])
        83
        sage: len([N for N in MatroidExtensions(M, 'x', line_length=5)])
        22
        sage: for N in MatroidExtensions(M, 'x', subsets=[[0, 1], [2, 3],
        ....:                                             [4, 5]]): print(N)
        Matroid of rank 3 on 7 elements with 32 bases
        Matroid of rank 3 on 7 elements with 20 bases
        sage: M = BasisMatroid(matroids.catalog.BetsyRoss()); M
        Matroid of rank 3 on 11 elements with 140 bases
        sage: e = 'k'; f = 'h'; Me = M.delete(e); Mf=M.delete(f)
        sage: for N in MatroidExtensions(Mf, f, splice=Me): print(N)
        Matroid of rank 3 on 11 elements with 141 bases
        Matroid of rank 3 on 11 elements with 140 bases
        sage: for N in MatroidExtensions(Me, e, splice=Mf): print(N)
        Matroid of rank 3 on 11 elements with 141 bases
        Matroid of rank 3 on 11 elements with 140 bases

    Note that this class is intended for runtime, internal use, so no
    loads/dumps mechanism was implemented.
    """
    def __init__(self, M, e, line_length=None, subsets=None, splice=None, orderly=False):
        """
        See class docstring for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.Uniform(3, 6)
            sage: len([N for N in MatroidExtensions(M, 'x')])
            83
            sage: len([N for N in MatroidExtensions(M, 'x', line_length=5)])
            22
            sage: for N in MatroidExtensions(M, 'x', subsets=[[0, 1], [2, 3],
            ....:                                            [4, 5]]): print(N)
            Matroid of rank 3 on 7 elements with 32 bases
            Matroid of rank 3 on 7 elements with 20 bases
        """
        if M.full_rank() == 0:
            pass
        if isinstance(M, BasisMatroid):
            BM = M
        else:
            BM = BasisMatroid(M)
        LinearSubclasses.__init__(self, BM, line_length=line_length, subsets=subsets, splice=splice)
        self._BX = BM._extension(e, [])
        self._BH = [BM._extension(e, [self._hyperplanes[i]]) for i in range(len(self._hyperplanes))]
        if orderly:
            self._orderly = True

    def __getitem__(self, CutNode node):
        """
        Return a single-element extension determined by a given CutNode.

        Internal function.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.Uniform(3, 6)
            sage: len([N for N in MatroidExtensions(M, 'x')])
            83
        """
        cdef long i
        X = BasisMatroid(self._BX)
        for i in node.planes():
            bitset_intersection(X._bb, X._bb, (<BasisMatroid>self._BH[i])._bb)
        X._reset_invariants()
        if self._orderly and not X.is_distinguished(len(X) - 1):
            return None
        if self._line_length > 0 and X.has_line_minor(self._line_length + 1):
            return None
        return X
