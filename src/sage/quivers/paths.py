"""
Quiver Paths
"""

#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import MonoidElement
from sage.rings.integer_ring import ZZ

class QuiverPath(MonoidElement):
    r"""
    Class for paths in a quiver.

    A path is given by two vertices, ``start`` and ``end``, and a finite
    (possibly empty) list of edges `e_1, e_2, \ldots, e_n` such that the
    initial vertex of `e_1` is ``start``, the final vertex of `e_i` is
    the initial vertex of `e_{i+1}`, and the final vertex of `e_n` is
    ``end``.  In the case where no edges are specified, we must have
    ``start = end`` and the path is called the trivial path at the given
    vertex.

    INPUT:

    - ``parent`` -- the path semigroup associated with a quiver; this is
      where the path will live
    - ``path`` -- tuple or iterable. If ``path`` is a tuple then it is
      assumed to be of the form ``(vertex, vertex)`` or
      ``(start, end, label)``.  In the first case the trivial path at the
      given vertex is created.  In the second case a path consisting of
      just the given edge is created.  If ``path`` is not a tuple then it
      is assumed to be an iterable variable giving the edges of a path,
      where each edge is in one of the two forms above.
    - ``check`` -- boolean (default: ``True``); if it is ``False``, no
      sanity check will be performed on the given iterable

    OUTPUT:

    - :class:`QuiverPath`

    .. NOTE::

        Do *not* use this constructor directly! Instead, pass the input to the
        path semigroup that shall be the parent of this path.

    EXAMPLES:

    Specify a path by giving a list of edges::

        sage: Q = DiGraph({1:{2:['a','d'], 3:['e']}, 2:{3:['b']}, 3:{1:['f'], 4:['c']}})
        sage: F = Q.path_semigroup()
        sage: p = F([(1, 2, 'a'), (2, 2), (2, 3, 'b')])
        sage: p
        a*b

    Paths are not *unique*, but different representations of "the same" path
    yield *equal* paths::

        sage: q = F([(1, 1), (1, 2, 'a'), (2, 3, 'b'), (3, 3)])
        sage: p is q
        False
        sage: p == q
        True

    The ``*`` operator is concatenation of paths. If the two paths do not
    compose, its result is ``None``::

        sage: print(p*q)
        None
        sage: p*F((3, 4, 'c'))
        a*b*c
        sage: F([(2,3,'b'), (3,1,'f')])*p
        b*f*a*b

    The length of a path is the number of edges in that path.  Trivial paths
    are therefore length-`0`::

        sage: len(p)
        2
        sage: triv = F((1, 1))
        sage: len(triv)
        0

    List index and slice notation can be used to access the edges in a path.
    QuiverPaths can also be iterated over.  Trivial paths have no elements::

        sage: for x in p: print x
        (1, 2, 'a')
        (2, 3, 'b')
        sage: triv[:]
        []

    There are methods giving the initial and terminal vertex of a path::

        sage: p.initial_vertex()
        1
        sage: p.terminal_vertex()
        3
    """
    def __init__(self, parent, path, check=True):
        """
        Creates a path object.  Type ``QuiverPath?`` for more information.

        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1), (1, 1)])
            sage: Q([(1,3,'x')])
            Traceback (most recent call last):
            ...
            ValueError: Cannot interpret [(1, 3, 'x')] as element of
            Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices

        Note that QuiverPath should not be called directly, because
        the elements of the path semigroup associated with a quiver
        use a sub-class of QuiverPath. Nonetheless, just for test, we
        show that it *is* possible to create a path in a deprecated way::

            sage: p == QuiverPath(Q, (1, 1))
            True
            sage: Q([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'b'), (3, 3)])._path
            ((1, 2, 'a'), (2, 3, 'b'))
        """
        MonoidElement.__init__(self, parent=parent)

        # Normalise the given path, unless it is asserted that the input is
        # fine
        if not check:
            self._path = tuple(path)
            return
        if path == 1:
            # We do not guarantee that there is only one vertex.
            # However, this element certainly exists.
            v = parent.quiver().vertices()[0]
            self._path = ((v,v),)
            return
        E = parent.quiver().edges()
        if isinstance(path, QuiverPath):
            if path.parent() is parent:
                self._path = path._path
                return
            new_path = list(path._path)
        # A tuple is assumed to be an edge, anything else is assumed to be a
        # list of edges
        elif isinstance(path, tuple):
            new_path = [path]
        else:
            new_path = list(path)

        # Check that each edge in the path is valid
        good = True
        for x in new_path:
            if (len(x) < 2 or x[0] not in ZZ or x[1] not in ZZ
                           or len(x) == 2 and x[0] != x[1]
                           or len(x) == 3 and x not in E
                           or len(x) > 3):
                good = False
                break
        if not good:
            raise ValueError("Cannot interpret %s as element of %s"%(path,parent))
        # Delete trivial edges, and clear the path if not valid
        i = 0
        while i + 1 < len(new_path):
            if new_path[i][1] != new_path[i + 1][0]:
                raise ValueError("Cannot interpret %s as element of %s"%(path,parent))
            elif len(new_path[i])!=3:
                del new_path[i]
            else:
                i += 1
        if len(new_path) > 1 and len(new_path[-1])!=3:
            del new_path[-1]
        self._path = tuple(new_path)

    def _repr_(self):
        r"""
        Default representation of a path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]) # indirect doctest
            a*b
            sage: Q((1, 1)) # indirect doctest
            e_1
        """
        if len(self._path[0])!=3:
            return 'e_{0}'.format(self._path[0][0])
        else:
            return '*'.join([e[2] for e in self._path])

    def __len__(self):
        """
        Return the length of the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: len(Q([(1, 2, 'a'), (2, 3, 'b')]))
            2
            sage: len(Q((1, 1)))
            0
            sage: len(Q((1, 2, 'a')))
            1
        """
        if (not self._path) or (len(self._path[0])==2):
            return 0
        else:
            return len(self._path)

    def deg(self):
        """
        Return the degree of the path, which is the same as its length.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]).deg()
            2
            sage: Q((1, 1)).deg()
            0
            sage: Q((1, 2, 'a')).deg()
            1
        """
        return len(self)

    def __nonzero__(self):
        """
        Implement boolean values for the object.

        .. NOTE::

            The boolean value is always ``True``, since the partial semigroup
            formed by the paths of a quiver does not contain zero.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((2, 3, 'b'))
            sage: bool(a*b)
            True
            sage: bool(Q.idempotents()[0])
            True
        """
        return True

    def __cmp__(self, other):
        """
        Comparison for :class:`QuiverPaths`.

        As usual in Sage, the ``__cmp__`` method of a Python sub-class of
        :class:`sage.structure.element.Element` can assume that both arguments
        belong to the same parent.

        If the QuiverPaths are unequal then one of the following data (listed
        in order of preferance) is unequal and used for comparison:

        - Length of the paths
        - Edge sequence of the paths.

        .. NOTE::

            This code is used by :class:`CombinatorialFreeModule` to order
            the monomials when printing elements of path algebras.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c'], 4:['d']}}).path_semigroup()
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((1, 2, 'b'))
            sage: c = Q((2, 3, 'c'))
            sage: d = Q((2, 4, 'd'))
            sage: e = Q.idempotents()[3]
            sage: e < a
            True
            sage: a < e
            False
            sage: d < a*c
            True
            sage: a*c < d
            False
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
        s_p = self._path
        o_p = other._path

        # Compare lengths if different
        c = cmp(len(s_p),len(o_p))
        if c:
            return c

        # Trivial paths should be smaller than non-trivial paths
        c = cmp(len(s_p[0]),len(o_p[0]))
        if c:
            return c

        # Now we have two non-trivial paths. Compare internal tuple
        return cmp(s_p, o_p)

    def __getitem__(self, *args):
        """
        Implement index notation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*c
            sage: p[0]
            (1, 2, 'a')
            sage: p[-1]
            (3, 4, 'c')
            sage: p[1:]
            ((2, 3, 'b'), (3, 4, 'c'))
        """
        if self._path and self._path[0][0] == self._path[0][1]:
            return list().__getitem__(*args)
        else:
            return self._path.__getitem__(*args)

    def __iter__(self):
        """
        Iteration over the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: for e in p: print e
            (1, 2, 'a')
            (2, 3, 'b')
            (3, 4, 'c')
        """
        # Return an iterator over an empty tuple for trivial paths, otherwise
        # return an iterator for _path as a list
        if not len(self):
            return list().__iter__()
        else:
            return list(self._path).__iter__()

    def _mul_(self, other):
        """
        Compose two paths.

        .. NOTE::

            ``None`` is returned if the terminal vertex of the first path
            does not coincide with the initial vertex of the second path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}, 4:{5:['d']}}).path_semigroup()
            sage: x = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = Q([(3, 4, 'c'), (4, 5, 'd')])
            sage: print y*x
            None
            sage: x*y
            a*b*c*d
            sage: x*Q((3, 4, 'c'))
            a*b*c
            sage: x*Q([(3, 4, 'c'), (4, 5, 'd')])
            a*b*c*d
            sage: x*6
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*':
             'Partial semigroup formed by the directed paths of Multi-digraph on 5 vertices'
             and 'Integer Ring'
        """
        # By Sage's coercion model, both paths belong to the same quiver
        # In particular, both are QuiverPath
        Q = self.parent()
        if self.terminal_vertex() != other.initial_vertex():
            return None
        if len(self._path[0]) < 3:
            return other
        elif len(other._path[0]) < 3:
            return self
        return Q(self._path+other._path, check=False)

    def __mod__(self, other):
        """
        Return ``self`` with ``other`` deleted from the beginning.

        If ``other`` is not the beginning of ``self`` then ``None`` is
        returned.  Deleting the trivial path at vertex `v` from a path that
        begins at `v` does not change the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((2, 3, 'b'))
            sage: e1 = Q((1, 1))
            sage: e2 = Q((2, 2))
            sage: p % a
            b
            sage: print p % b
            None
            sage: p % e1
            a*b
            sage: print p % e2
            None
        """
        Q = self.parent()
        # Convert other to a QuiverPath
        oth = Q(other)
        if self._path == oth._path:
            v = self._path[-1][1]
            return Q(((v, v),), check=False)

        # Handle trivial paths
        if oth._path[0][0] == oth._path[0][1]:
            if self._path[0][0] == oth._path[0][0]:
                return self
            else:
                return None

        # If other is the beginning, return the rest
        if self._path[:len(oth._path)] == oth._path:
            return Q(list(self._path[len(oth._path):]), check=False)
        else:
            return None

    def initial_vertex(self):
        """
        Return the initial vertex of the path.

        OUTPUT:

        - integer, the label of the initial vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.initial_vertex()
            1
        """
        return self._path[0][0]

    def terminal_vertex(self):
        """
        Return the terminal vertex of the path.

        OUTPUT:

        - integer, the label of the terminal vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.terminal_vertex()
            3
        """
        return self._path[-1][1]

    def reverse(self):
        """
        Return the path along the same edges in reverse order in the
        opposite quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: p.reverse()
            b*a
            sage: p.reverse().parent() is Q.reverse()
            True
            sage: e = Q.idempotents()[0]
            sage: e
            e_1
            sage: e.reverse()
            e_1
        """
        Q = self.parent().reverse()
        # Handle trivial paths
        if self._path[0][0] == self._path[0][1]:
            return Q(((self._path[0][0],self._path[0][0]),), check=False)

        # Reverse all the edges in the path, then reverse the path
        new_path = [(e[1], e[0], e[2]) for e in self._path]
        return Q(reversed(new_path), check=False)

