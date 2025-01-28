"""
Quiver Paths
"""
# ****************************************************************************
#  Copyright (C) 2012    Jim Stark <jstarx@gmail.com>
#                2013/14 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.data_structures.bounded_integer_sequences cimport *
from cpython.slice cimport PySlice_GetIndicesEx
from sage.structure.richcmp cimport rich_to_bool
from sage.data_structures.bitset_base cimport *


cdef class QuiverPath(MonoidElement):
    r"""
    Class for paths in a quiver.

    A path is given by two vertices, ``start`` and ``end``, and a finite
    (possibly empty) list of edges `e_1, e_2, \ldots, e_n` such that the
    initial vertex of `e_1` is ``start``, the final vertex of `e_i` is
    the initial vertex of `e_{i+1}`, and the final vertex of `e_n` is
    ``end``.  In the case where no edges are specified, we must have
    ``start = end`` and the path is called the trivial path at the given
    vertex.

    .. NOTE::

        Do *not* use this constructor directly! Instead, pass the input to the
        path semigroup that shall be the parent of this path.

    EXAMPLES:

    Specify a path by giving a list of edges::

        sage: Q = DiGraph({1:{2:['a','d'], 3:['e']}, 2:{3:['b']}, 3:{1:['f'], 4:['c']}})
        sage: F = Q.path_semigroup()
        sage: p = F([(1, 2, 'a'), (2, 3, 'b')])
        sage: p
        a*b

    Paths are not *unique*, but different representations of "the same" path
    yield *equal* paths::

        sage: q = F([(1, 1)]) * F([(1, 2, 'a'), (2, 3, 'b')]) * F([(3, 3)])
        sage: p is q
        False
        sage: p == q
        True

    The ``*`` operator is concatenation of paths. If the two paths do not
    compose, its result is ``None``::

        sage: print(p*q)
        None
        sage: p*F([(3, 4, 'c')])
        a*b*c
        sage: F([(2,3,'b'), (3,1,'f')])*p
        b*f*a*b

    The length of a path is the number of edges in that path.  Trivial paths
    are therefore length-`0`::

        sage: len(p)
        2
        sage: triv = F([(1, 1)])
        sage: len(triv)
        0

    List index and slice notation can be used to access the edges in a path.
    QuiverPaths can also be iterated over.  Trivial paths have no elements::

        sage: for x in p: print(x)
        (1, 2, 'a')
        (2, 3, 'b')
        sage: list(triv)
        []

    There are methods giving the initial and terminal vertex of a path::

        sage: p.initial_vertex()
        1
        sage: p.terminal_vertex()
        3
    """
    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1)]) * Q([(1, 1)])
            sage: del p    # indirect doctest
        """
        biseq_dealloc(self._path)

    cdef QuiverPath _new_(self, int start, int end):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])    # indirect doctest
        """
        cdef QuiverPath out = QuiverPath.__new__(self._parent.element_class)
        out._parent = self._parent
        out._start = start
        out._end = end
        return out

    def __init__(self, parent, start, end, path):
        """
        Create a path object.  Type ``QuiverPath?`` for more information.

        INPUT:

        - ``parent`` -- a path semigroup
        - ``start`` -- integer; the label of the initial vertex
        - ``end`` -- integer; the label of the terminal vertex
        - ``path`` -- list of integers, providing the list of arrows
          occurring in the path, labelled according to the position in
          the list of all arrows (resp. the list of outgoing arrows at
          each vertex)

        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1)]) * Q([(1, 1)])
            sage: Q([(1,3,'x')])
            Traceback (most recent call last):
            ...
            ValueError: (1, 3, 'x') is not an edge

        Note that QuiverPath should not be called directly, because
        the elements of the path semigroup associated with a quiver
        may use a sub-class of QuiverPath. Nonetheless, just for test, we
        show that it *is* possible to create a path in a deprecated way::

            sage: p == QuiverPath(Q, 1, 1, [])
            True
            sage: list(Q([(1, 1)])*Q([(1, 2, 'a')])*Q([(2, 2)])*Q([(2, 3, 'b')])*Q([(3, 3)]))
            [(1, 2, 'a'), (2, 3, 'b')]
        """
        MonoidElement.__init__(self, parent=parent)
        self._start = start
        self._end = end
        biseq_init_list(self._path, path, parent._nb_arrows)

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])
            sage: loads(dumps(p)) == p   # indirect doctest
            True
            sage: loads(dumps(p)) is p
            False
        """
        return NewQuiverPath, (self._parent, self._start, self._end,
                               biseq_pickle(self._path))

    def __hash__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])
            sage: q = Q([(1, 1)])
            sage: {p:1, q:2}[Q(['a','b'])]    # indirect doctest
            1
        """
        if self._path.length == 0:
            return hash(self._start)
        cdef Py_hash_t h = self._start*(<Py_hash_t>1073807360) + biseq_hash(self._path)
        if h == -1:
            return -2
        return h
        # bitset_hash is not a good hash either
        # We should consider using FNV-1a hash, see http://www.isthe.com/chongo/tech/comp/fnv/,
        # Or the hash defined in http://burtleburtle.net/bob/hash/doobs.html
        # Or http://www.azillionmonkeys.com/qed/hash.html

    def _repr_(self):
        r"""
        Default representation of a path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]) # indirect doctest
            a*b
            sage: Q([(1, 1)]) # indirect doctest
            e_1
        """
        cdef mp_size_t i
        if not self._path.length:
            return 'e_{0}'.format(self._start)
        L = self._parent._labels
        return '*'.join(L[biseq_getitem(self._path, i)]
                        for i in range(self._path.length))

    def __len__(self):
        """
        Return the length of the path.

        ``length()`` and ``degree()`` are aliases

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: len(Q([(1, 2, 'a'), (2, 3, 'b')]))
            2
            sage: Q([(1, 1)]).degree()
            0
            sage: Q([(1, 2, 'a')]).length()
            1
        """
        return self._path.length

    degree = __len__
    length = __len__

    def __bool__(self):
        """
        Implement boolean values for paths.

        .. NOTE::

            The boolean value is ``True`` if and only if this path is of
            positive length.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(2, 3, 'b')])
            sage: bool(a*b)
            True
            sage: bool(Q.idempotents()[0])
            False
        """
        return self._path.length != 0

    cpdef _richcmp_(left, right, int op):
        """
        Comparison for :class:`QuiverPaths`.

        The following data (listed in order of preference) is used for
        comparison:

        - **Negative** length of the paths
        - initial and terminal vertices of the paths
        - Edge sequence of the paths, by reverse lexicographical ordering.

        .. NOTE::

            This code is used by :class:`CombinatorialFreeModule` to order
            the monomials when printing elements of path algebras.

        EXAMPLES:

        A path semigroup of a quiver with a single vertex has a multiplicative
        unit::

            sage: D = DiGraph({0:{0:['x','y','z']}}).path_semigroup()
            sage: D(1)
            e_0
            sage: D in Monoids()
            True

        However, there is no coercion from the ring of integers and hence the
        generic comparison code finds that the multiplicative units in `D` and
        in the ring of integers evaluate unequal::

            sage: D.has_coerce_map_from(ZZ)
            False
            sage: D(1) == 1
            False

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c'], 4:['d']}}).path_semigroup()
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(1, 2, 'b')])
            sage: c = Q([(2, 3, 'c')])
            sage: d = Q([(2, 4, 'd')])
            sage: e = Q.idempotents()[3]
            sage: e < a  # e is shorter than a
            False
            sage: a < e
            True
            sage: d < a*c
            False
            sage: a*c < d
            True
            sage: a < b
            True
            sage: b < a
            False
            sage: a*c < a*d
            True
            sage: a*d < a*c
            False
            sage: a < a
            False
        """
        # Since QuiverPath inherits from Element, it is guaranteed that
        # both arguments are elements of the same path semigroup
        cdef QuiverPath cself, other
        cself = left
        other = right
        # we want *negative* degree reverse lexicographical order
        if other._path.length < cself._path.length:
            return rich_to_bool(op, -1)
        if other._path.length > cself._path.length:
            return rich_to_bool(op, 1)
        if cself._start < other._start:
            return rich_to_bool(op, -1)
        if cself._start > other._start:
            return rich_to_bool(op, 1)
        if cself._end < other._end:
            return rich_to_bool(op, -1)
        if cself._end > other._end:
            return rich_to_bool(op, 1)
        if cself._path.length == 0:
            return rich_to_bool(op, 0)
        return biseq_richcmp(cself._path, other._path, op)

    def __getitem__(self, index):
        """
        Implement index notation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c'], 1:['d']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 1, 'd'), (1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*d*a*b*c

        A single index returns the arrow that appears in the path at this index::

            sage: p[0]
            a
            sage: p[-1]
            c

        A slice with step 1 returns a sub-path of this path::

            sage: p[1:5]
            b*d*a*b

        A slice with step -1 return a sub-path of the reversed path::

            sage: p[4:1:-1]
            b*a*d

        If the start index is greater than the terminal index and the step
        -1 is not explicitly given, then a path of length zero is returned,
        which is compatible with Python lists::

            sage: list(range(6))[4:1]
            []

        The following was fixed in :issue:`22278`. A path slice of length
        zero of course has a specific start- and endpoint. It is always
        the startpoint of the arrow corresponding to the first item of
        the range::

            sage: p[4:1]
            e_2
            sage: p[4:1].initial_vertex() == p[4].initial_vertex()
            True

        If the slice boundaries are out of bound, then no error is raised,
        which is compatible with Python lists::

            sage: list(range(6))[20:40]
            []

        In that case, the startpoint of the slice of length zero is the
        endpoint of the path::

            sage: p[20:40]
            e_4
            sage: p[20:40].initial_vertex() == p.terminal_vertex()
            True
        """
        cdef tuple E
        cdef Py_ssize_t start, stop, step, slicelength
        cdef int init, end
        cdef QuiverPath OUT
        if isinstance(index, slice):
            PySlice_GetIndicesEx(index, self._path.length,
                                 &start, &stop, &step,
                                 &slicelength)
            if step != 1 and step != -1:
                raise ValueError("slicing only possible for step +/-1")
            if step == -1:
                return self.reversal()[self._path.length-1-start:self._path.length-1-stop]
            if start == 0 and stop == self._path.length:
                return self
            if start > stop:
                stop=start
            E = self._parent._sorted_edges
            if start < self._path.length:
                init = E[biseq_getitem(self._path, start)][0]
            else:
                init = self._end
            if start < stop:
                end = E[biseq_getitem(self._path, stop-1)][1]
            else:  # the result will be a path of length 0
                end = init
            OUT = self._new_(init, end)
            biseq_init_slice(OUT._path, self._path, start, stop, step)
            return OUT
        if index < 0:
            index = self._path.length+index
        if index < 0 or index >= self._path.length:
            raise IndexError("list index out of range")
        E = self._parent._sorted_edges
        init = E[biseq_getitem(self._path, index)][0]
        end = E[biseq_getitem(self._path, index)][1]
        OUT = self._new_(init, end)
        biseq_init_slice(OUT._path, self._path, index, index+1, 1)
        return OUT

    def __iter__(self):
        """
        Iteration over the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: for e in p: print(e)
            (1, 2, 'a')
            (2, 3, 'b')
            (3, 4, 'c')
        """
        # Return an iterator over an empty tuple for trivial paths, otherwise
        # return an iterator for _path as a list
        cdef mp_size_t i
        cdef tuple E = self._parent._sorted_edges
        for i in range(self._path.length):
            yield E[biseq_getitem(self._path, i)]

    cpdef _mul_(self, other):
        """
        Compose two paths.

        .. NOTE::

            ``None`` is returned if the terminal vertex of the first path
            does not coincide with the initial vertex of the second path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}, 4:{5:['d']}}).path_semigroup()
            sage: x = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = Q([(3, 4, 'c'), (4, 5, 'd')])
            sage: print(y*x)
            None
            sage: x*y
            a*b*c*d
            sage: x*Q([(3, 4, 'c')])
            a*b*c
            sage: x*Q([(3, 4, 'c'), (4, 5, 'd')])
            a*b*c*d
            sage: x*6
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
             'Partial semigroup formed by the directed paths of Multi-digraph on 5 vertices'
             and 'Integer Ring'
        """
        # By Sage's coercion model, both paths belong to the same quiver
        # In particular, both are QuiverPath
        cdef QuiverPath right = other
        if self._end != right._start:
            return None
        cdef QuiverPath OUT = self._new_(self._start, right._end)
        biseq_init_concat(OUT._path, self._path, right._path)
        return OUT

    cpdef _mod_(self, other):
        """
        Return what remains of this path after removing the initial segment ``other``.

        If ``other`` is not an initial segment of this path then ``None`` is
        returned.  Deleting the trivial path at vertex `v` from a path that
        begins at `v` does not change the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(2, 3, 'b')])
            sage: e1 = Q([(1, 1)])
            sage: e2 = Q([(2, 2)])
            sage: p % a
            b
            sage: print(p % b)
            None
            sage: p % e1
            a*b
            sage: print(p % e2)
            None
        """
        cdef QuiverPath right = <QuiverPath>other
        # Handle trivial case
        if self._start != right._start:
            return None
        if right._path.length == 0:
            return self

        # If other is the beginning, return the rest
        cdef QuiverPath OUT
        if (self._start == right._start) and biseq_startswith(self._path, right._path):
            OUT = self._new_(right._end, self._end)
            biseq_init_slice(OUT._path, self._path, right._path.length, self._path.length, 1)
            return OUT
        else:
            return None

    def gcd(self, QuiverPath P):
        """
        Greatest common divisor of two quiver paths, with co-factors.

        For paths, by "greatest common divisor", we mean the largest terminal
        segment of the first path that is an initial segment of the second
        path.

        INPUT:

        - ``P`` -- a :class:`QuiverPath`

        OUTPUT:

        - :class:`QuiverPath`s ``(C1,G,C2)`` such that ``self = C1*G`` and ``P = G*C2``, or
        - ``(None, None, None)``, if the paths do not overlap (or belong to different quivers).

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{1:['b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: p1 = Q(['c','d','a','b','a','c','d'])
            sage: p1
            c*d*a*b*a*c*d
            sage: p2 = Q(['a','b','a','c','d','a','c','d','a','b'])
            sage: p2
            a*b*a*c*d*a*c*d*a*b
            sage: S1, G, S2 = p1.gcd(p2)
            sage: S1, G, S2
            (c*d, a*b*a*c*d, a*c*d*a*b)
            sage: S1*G == p1
            True
            sage: G*S2 == p2
            True
            sage: p2.gcd(p1)
            (a*b*a*c*d*a, c*d*a*b, a*c*d)

        We test that a full overlap is detected::

            sage: p2.gcd(p2)
            (e_1, a*b*a*c*d*a*c*d*a*b, e_1)

        The absence of an overlap is detected::

            sage: p2[2:-1]
            a*c*d*a*c*d*a
            sage: p2[1:]
            b*a*c*d*a*c*d*a*b
            sage: print(p2[2:-1].gcd(p2[1:]))
            (None, None, None)
        """
        if self._parent is not P._parent:
            return (None, None, None)
        cdef size_t i
        sig_on()
        i = biseq_startswith_tail(P._path, self._path, 0)
        sig_off()
        if i == <size_t>-1:
            return (None, None, None)
        return (self[:i], self[i:], P[self._path.length-i:])

    cpdef tuple complement(self, QuiverPath subpath):
        """
        Return a pair ``(a,b)`` of paths s.t. ``self = a*subpath*b``,
        or ``(None, None)`` if ``subpath`` is not a subpath of this path.

        .. NOTE::

            ``a`` is chosen of minimal length.

        EXAMPLES::

            sage: S = DiGraph({1:{1:['a','b','c','d']}}).path_semigroup()
            sage: S.inject_variables()
            Defining e_1, a, b, c, d
            sage: (b*c*a*d*b*a*d*d).complement(a*d)
            (b*c, b*a*d*d)
            sage: (b*c*a*d*b).complement(a*c)
            (None, None)
        """
        cdef mp_size_t i = biseq_contains(self._path, subpath._path, 0)
        if i == -1:
            return (None, None)
        return self[:i], self[i+len(subpath):]

    cpdef bint has_subpath(self, QuiverPath subpath) except -1:
        """
        Tells whether this path contains a given sub-path.

        INPUT:

        - ``subpath`` -- a path of positive length in the same path semigroup
          as this path

        EXAMPLES::

            sage: S = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup()
            sage: S.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: (c*b*e*a).has_subpath(b*e)
            1
            sage: (c*b*e*a).has_subpath(b*f)
            0
            sage: (c*b*e*a).has_subpath(e_1)
            Traceback (most recent call last):
            ...
            ValueError: we only consider sub-paths of positive length
            sage: (c*b*e*a).has_subpath(None)
            Traceback (most recent call last):
            ...
            ValueError: the given sub-path is empty
        """
        if subpath is None:
            raise ValueError("the given sub-path is empty")
        if subpath._parent is not self._parent:
            raise ValueError("the two paths belong to different quivers")
        if subpath._path.length == 0:
            raise ValueError("we only consider sub-paths of positive length")
        if self._path.length < subpath._path.length:
            return 0
        if biseq_contains(self._path, subpath._path, 0) == -1:
            return 0
        return 1

    cpdef bint has_prefix(self, QuiverPath subpath) except -1:
        """
        Tells whether this path starts with a given sub-path.

        INPUT:

        - ``subpath`` -- a path in the same path semigroup as this path

        OUTPUT: ``0`` or ``1``, which stands for ``False`` resp. ``True``

        EXAMPLES::

            sage: S = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup()
            sage: S.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: (c*b*e*a).has_prefix(b*e)
            0
            sage: (c*b*e*a).has_prefix(c*b)
            1
            sage: (c*b*e*a).has_prefix(e_1)
            1
            sage: (c*b*e*a).has_prefix(e_2)
            0
        """
        if subpath._parent is not self._parent:
            raise ValueError("the two paths belong to different quivers")
        if self._start != subpath._start:
            return 0
        if subpath._path.length == 0:
            return 1
        if biseq_startswith(self._path, subpath._path):
            return 1
        return 0

    def initial_vertex(self):
        """
        Return the initial vertex of the path.

        OUTPUT: integer; the label of the initial vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.initial_vertex()
            1
        """
        return self._start

    def terminal_vertex(self):
        """
        Return the terminal vertex of the path.

        OUTPUT: integer; the label of the terminal vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.terminal_vertex()
            3
        """
        return self._end

    def reversal(self):
        """
        Return the path along the same edges in reverse order in the
        opposite quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c'], 1:['d']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 1, 'd'), (1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*d*a*b*c
            sage: p.reversal()
            c*b*a*d*b*a
            sage: e = Q.idempotents()[0]
            sage: e
            e_1
            sage: e.reversal()
            e_1
        """
        Q = self._parent.reverse()
        # Handle trivial paths
        if self._path.length == 0:
            return Q.element_class(Q, self._end, self._start, [])

        # Reverse all the edges in the path, then reverse the path
        cdef mp_size_t i
        cdef QuiverPath out = QuiverPath.__new__(Q.element_class)
        out._parent = Q
        out._start = self._end
        out._end = self._start
        sig_check()
        biseq_init(out._path, self._path.length, self._path.itembitsize)
        cdef mp_size_t l = self._path.length - 1
        for i in range(self._path.length):
            sig_check()
            biseq_inititem(out._path, i, biseq_getitem(self._path, l-i))
        return out


@cython.binding(True)
def NewQuiverPath(Q, start, end, biseq_data):
    """
    Return a new quiver path for given defining data.

    INPUT:

    - ``Q`` -- the path semigroup of a quiver
    - ``start`` -- integer; the label of the startpoint
    - ``end`` -- integer; the label of the endpoint
    - ``biseq_data`` -- tuple formed by

      - a string, encoding a bitmap representing the path as integer
        at base `32`
      - the number of bits used to store the path
      - the number of bits used to store a single item
      - the number of items in the path

    TESTS::

        sage: from sage.quivers.paths import QuiverPath
        sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: p = Q(['a']) * Q(['b'])
        sage: loads(dumps(p)) == p   # indirect doctest
        True
        sage: p.__reduce__()
        (<cyfunction NewQuiverPath at ...>,
         (Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices,
          1,
          3,
          ((0, 4, 1, ..., (4,)), 2, 2)))
    """
    cdef QuiverPath out = QuiverPath.__new__(Q.element_class)
    out._parent = Q
    out._start = start
    out._end = end
    biseq_unpickle(out._path, biseq_data[0], biseq_data[1], biseq_data[2])
    return out
