"""
Data structures for maps between finite sets

This module implements several fast Cython data structures for maps
between two finite set. Those classes are not intended to be used
directly. Instead, such a map should be constructed via its parent,
using the class :class:`~sage.sets.finite_set_maps.FiniteSetMaps`.

EXAMPLES:

To create a map between two sets, one first creates the set of such maps::

    sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5])

The map can then be constructed either from a function::

    sage: f1 = M(lambda c: ord(c)-94); f1
    map: a -> 3, b -> 4

or from a dictionary::

    sage: f2 = M.from_dict({'a':3, 'b':4}); f2
    map: a -> 3, b -> 4

The two created maps are equal::

    sage: f1 == f2
    True

Internally, maps are represented as the list of the ranks of the
images ``f(x)`` in the co-domain, in the order of the domain::

    sage: list(f2)
    [0, 1]

A third fast way to create a map it to use such a list. it should be
kept for internal use::

    sage: f3 = M._from_list_([0, 1]); f3
    map: a -> 3, b -> 4
    sage: f1 == f3
    True

AUTHORS:

- Florent Hivert
"""
# ****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.list_clone cimport ClonableIntArray
from sage.structure.parent cimport Parent
from sage.arith.power cimport generic_power
from sage.sets.set import Set_object_enumerated


cpdef fibers(f, domain):
    r"""
    Return the fibers of the function ``f`` on the finite set ``domain``.

    INPUT:

    - ``f`` -- a function or callable
    - ``domain`` -- a finite iterable

    OUTPUT:

    - a dictionary ``d`` such that ``d[y]`` is the set of all ``x`` in
      ``domain`` such that ``f(x) = y``

    EXAMPLES::

        sage: from sage.sets.finite_set_map_cy import fibers, fibers_args
        sage: fibers(lambda x: 1, [])
        {}
        sage: fibers(lambda x: x^2, [-1, 2, -3, 1, 3, 4])
        {1: {1, -1}, 4: {2}, 9: {3, -3}, 16: {4}}
        sage: fibers(lambda x: 1,   [-1, 2, -3, 1, 3, 4])
        {1: {1, 2, 3, 4, -3, -1}}
        sage: fibers(lambda x: 1, [1,1,1])
        {1: {1}}

    .. SEEALSO:: :func:`fibers_args` if one needs to pass extra
       arguments to ``f``.
    """
    result = {}
    for x in domain:
        y = f(x)
        result.setdefault(y, set()).add(x)
    for x, v in result.iteritems():
        result[x] = Set_object_enumerated(v)
    return result


def fibers_args(f, domain, *args, **opts):
    r"""
    Return the fibers of the function ``f`` on the finite set ``domain``.

    It is the same as :func:`fibers` except that one can pass extra
    argument for ``f`` (with a small overhead)

    EXAMPLES::

        sage: from sage.sets.finite_set_map_cy import fibers_args
        sage: fibers_args(operator.pow, [-1, 2, -3, 1, 3, 4], 2)
        {1: {1, -1}, 4: {2}, 9: {3, -3}, 16: {4}}
    """
    return fibers(lambda x: f(x, *args, **opts), domain)


cdef class FiniteSetMap_MN(ClonableIntArray):
    """
    Data structure for maps from ``range(m)`` to ``range(n)``.

    We assume that the parent given as argument is such that:

    - ``m`` -- stored in ``self.parent()._m``
    - ``n`` -- stored in ``self.parent()._n``
    - the domain is in ``self.parent().domain()``
    - the codomain is in ``self.parent().codomain()``
    """

    def __call__(self, int i):
        """
        Return the image of ``i`` under the map ``self``.

        INPUT:

        - ``i`` -- integer

        OUTPUT: integer

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: fs(0), fs(1), fs(2), fs(3)
            (1, 0, 2, 1)
        """
        return self._getitem(i)

    # Needed by generic power which refuses to compute 0^0
    def __bool__(self):
        """
        Return whether ``self`` is nonzero; this is always ``True``.

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: bool(fs)
            True
        """
        return True

    cpdef domain(self):
        """
        Return the domain of ``self``.

        EXAMPLES::

            sage: FiniteSetMaps(4, 3)([1, 0, 2, 1]).domain()
            {0, 1, 2, 3}
        """
        return self._parent.domain()

    cpdef codomain(self):
        """
        Return the codomain of ``self``.

        EXAMPLES::

            sage: FiniteSetMaps(4, 3)([1, 0, 2, 1]).codomain()
            {0, 1, 2}
        """
        return self._parent.codomain()

    cpdef _setimage(self, int i, int j):
        """
        Set the image of ``i`` as ``j`` in ``self``.

        This is a fast internal version where ``i`` and ``j`` are both int's.

        .. warning:: ``self`` must be mutable; otherwise an exception is raised.

        INPUT:

        - ``i``, ``j`` -- integers

        OUTPUT: none

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: fs2 = copy(fs)
            sage: fs2._setimage(2, 1)
            sage: fs2
            [1, 0, 1, 1]
            sage: with fs.clone() as fs3:
            ....:     fs3._setimage(0, 2)
            ....:     fs3._setimage(1, 2)
            sage: fs3
            [2, 2, 2, 1]

        TESTS::

            sage: with fs.clone() as fs3:
            ....:     fs3._setimage(6, 2)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: with fs.clone() as fs3:
            ....:     fs3._setimage(1, 4)
            Traceback (most recent call last):
            ...
            AssertionError: Wrong value self(1) = 4
        """
        self._setitem(i, j)

    cpdef _getimage(self, int i):
        """
        Return the image of ``i`` by ``self``.

        This is a fast internal version where ``i`` is integer.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: fs._getimage(0), fs._getimage(1), fs._getimage(2), fs._getimage(3)
            (1, 0, 2, 1)
        """
        return self._getitem(i)

    cpdef setimage(self, i, j):
        """
        Set the image of ``i`` as ``j`` in ``self``.

        .. warning:: ``self`` must be mutable; otherwise an exception is raised.

        INPUT:

        - ``i``, ``j`` -- two ``object``'s

        OUTPUT: none

        .. NOTE:: if you need speed, please use instead :meth:`_setimage`

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: fs2 = copy(fs)
            sage: fs2.setimage(2, 1)
            sage: fs2
            [1, 0, 1, 1]
            sage: with fs.clone() as fs3:
            ....:     fs3.setimage(0, 2)
            ....:     fs3.setimage(1, 2)
            sage: fs3
            [2, 2, 2, 1]
        """
        self._setitem(int(i), int(j))

    cpdef getimage(self, i):
        """
        Return the image of ``i`` by ``self``.

        INPUT:

        - ``i`` -- any object

        .. NOTE:: if you need speed, please use instead :meth:`_getimage`

        EXAMPLES::

            sage: fs = FiniteSetMaps(4, 3)([1, 0, 2, 1])
            sage: fs.getimage(0), fs.getimage(1), fs.getimage(2), fs.getimage(3)
            (1, 0, 2, 1)
        """
        return self._getitem(int(i))

    cpdef image_set(self):
        """
        Return the image set of ``self``.

        EXAMPLES::

            sage: FiniteSetMaps(4, 3)([1, 0, 2, 1]).image_set()
            {0, 1, 2}
            sage: FiniteSetMaps(4, 3)([1, 0, 0, 1]).image_set()
            {0, 1}
        """
        return Set_object_enumerated(self)

    cpdef fibers(self):
        """
        Return the fibers of ``self``.

        OUTPUT:

            a dictionary ``d`` such that ``d[y]`` is the set of all ``x`` in
            ``domain`` such that ``f(x) = y``

        EXAMPLES::

            sage: FiniteSetMaps(4, 3)([1, 0, 2, 1]).fibers()
            {0: {1}, 1: {0, 3}, 2: {2}}
            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: F.from_dict({"a": "b", "b": "a", "c": "b"}).fibers() == {'a': {'b'}, 'b': {'a', 'c'}}
            True
        """
        return fibers(self, self.domain())

    cpdef items(self):
        """
        The items of ``self``.

        Return the list of the ordered pairs ``(x, self(x))``

        EXAMPLES::

            sage: FiniteSetMaps(4, 3)([1, 0, 2, 1]).items()
            [(0, 1), (1, 0), (2, 2), (3, 1)]
        """
        return [(i, self._getimage(i)) for i in self.domain()]

    cpdef check(self):
        """
        Perform checks on ``self``.

        Check that ``self`` is a proper function and then calls
        ``parent.check_element(self)`` where ``parent`` is the parent
        of ``self``.

        TESTS::

            sage: fs = FiniteSetMaps(3, 2)
            sage: for el in fs: el.check()
            sage: fs([1,1])
            Traceback (most recent call last):
            ...
            AssertionError: Wrong number of values
            sage: fs([0,0,2])
            Traceback (most recent call last):
            ...
            AssertionError: Wrong value self(2) = 2
        """
        cdef int i, m, n
        m = self._parent._m
        n = self._parent._n
        assert self._len == m, "Wrong number of values"
        for i in range(m):
            assert 0 <= self._list[i] < n, "Wrong value self(%i) = %i" % (i, self._list[i])
        if hasattr(self._parent, 'check_element'):
            self._parent.check_element(self)

    cpdef FiniteSetMap_MN _compose_internal_(self, FiniteSetMap_MN other,
                                             Parent resParent):
        """
        TESTS::

            sage: FSM = FiniteSetMaps(3)
            sage: fs1 = FSM([1, 0, 2])
            sage: fs2 = FSM([0, 1, 2])
            sage: el = fs1*fs2; el
            [1, 0, 2]
            sage: el.check()
        """
        cdef type t = type(self)
        cdef FiniteSetMap_MN res = t.__new__(t)
        res._parent = resParent
        res._alloc_(self._len)
        for i in range(self._len):
            res._list[i] = other._list[self._list[i]]
        res.set_immutable()
        return res


cpdef FiniteSetMap_Set FiniteSetMap_Set_from_list(t, parent, lst):
    """
    Create a ``FiniteSetMap`` from a list.

    .. warning:: no check is performed !

    TESTS::

        sage: from sage.sets.finite_set_map_cy import FiniteSetMap_Set_from_list as from_list
        sage: F = FiniteSetMaps(["a", "b"], [3, 4, 5])
        sage: f = from_list(F.element_class, F, [0,2]); f.check(); f
        map: a -> 3, b -> 5
        sage: f.parent() is F
        True
        sage: f.is_immutable()
        True
    """
    cdef FiniteSetMap_MN res
    cdef type cls = <type>t
    res = cls.__new__(cls)
    super(FiniteSetMap_MN, res).__init__(parent, lst)
    return res

cpdef FiniteSetMap_Set FiniteSetMap_Set_from_dict(t, parent, d):
    """
    Create a ``FiniteSetMap`` from a dictionary.

    .. warning:: no check is performed !

    TESTS::

        sage: from sage.sets.finite_set_map_cy import FiniteSetMap_Set_from_dict as from_dict
        sage: F = FiniteSetMaps(["a", "b"], [3, 4, 5])
        sage: f = from_dict(F.element_class, F, {"a": 3, "b": 5}); f.check(); f
        map: a -> 3, b -> 5
        sage: f.parent() is F
        True
        sage: f.is_immutable()
        True
    """
    cdef FiniteSetMap_Set res
    cdef type cls = <type>t
    res = cls.__new__(cls)
    res.__init__(parent, d.__getitem__)
    return res


cdef class FiniteSetMap_Set(FiniteSetMap_MN):
    """
    Data structure for maps.

    We assume that the parent given as argument is such that:

    - the domain is in ``parent.domain()``
    - the codomain is in ``parent.codomain()``
    - ``parent._m`` contains the cardinality of the domain
    - ``parent._n`` contains the cardinality of the codomain
    - ``parent._unrank_domain`` and ``parent._rank_domain`` is a pair of
      reciprocal rank and unrank functions between the domain and
      ``range(parent._m)``.
    - ``parent._unrank_codomain`` and ``parent._rank_codomain`` is a pair of
      reciprocal rank and unrank functions between the codomain and
      ``range(parent._n)``.
    """

    def __init__(self, parent, fun, check=True):
        """
        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c", "d"], ["u", "v", "w"])
            sage: F(lambda x: "v")
            map: a -> v, b -> v, c -> v, d -> v
        """
        # For speed we initialize self by hand.
        # super(FiniteSetMap_Set, self).__init__(parent, [], check=False)
        self._parent = parent
        self._alloc_(int(parent._m))
        for i, el in enumerate(parent.domain().list()):
            self._setitem(i, parent._rank_codomain(fun(el)))
        self.set_immutable()
        if check:
            self.check()

    from_list = classmethod(FiniteSetMap_Set_from_list)
    from_dict = classmethod(FiniteSetMap_Set_from_dict)

    def __call__(self, i):
        """
        Return the image of ``i`` under the map ``self``.

        INPUT:

        - ``i`` -- integer

        OUTPUT: integer

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b"], [3, 4, 5])
            sage: fs = F.from_dict({"a": 3, "b": 5})
            sage: fs('a'), fs('b')
            (3, 5)
        """
        parent = self._parent
        return parent._unrank_codomain(self._getitem(parent._rank_domain(i)))

    cpdef image_set(self):
        """
        Return the image set of ``self``.

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: sorted(F.from_dict({"a": "b", "b": "a", "c": "b"}).image_set())
            ['a', 'b']
            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: F(lambda x: "c").image_set()
            {'c'}
        """
        image_i = self._parent._unrank_codomain
        return Set_object_enumerated([image_i(i) for i in self])

    cpdef setimage(self, i, j):
        """
        Set the image of ``i`` as ``j`` in ``self``.

        .. warning:: ``self`` must be mutable otherwise an exception is raised.

        INPUT:

        - ``i``, ``j`` -- two ``object``'s

        OUTPUT: none

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c", "d"], ["u", "v", "w"])
            sage: fs = F(lambda x: "v")
            sage: fs2 = copy(fs)
            sage: fs2.setimage("a", "w")
            sage: fs2
            map: a -> w, b -> v, c -> v, d -> v
            sage: with fs.clone() as fs3:
            ....:     fs3.setimage("a", "u")
            ....:     fs3.setimage("c", "w")
            sage: fs3
            map: a -> u, b -> v, c -> w, d -> v

        TESTS::

            sage: with fs.clone() as fs3:
            ....:     fs3.setimage("z", 2)
            Traceback (most recent call last):
            ...
            ValueError: 'z' is not in dict

            sage: with fs.clone() as fs3:
            ....:     fs3.setimage(1, 4)
            Traceback (most recent call last):
            ...
            ValueError: 1 is not in dict
        """
        parent = self._parent
        return self._setitem(parent._rank_domain(i), parent._rank_codomain(j))

    cpdef getimage(self, i):
        """
        Return the image of ``i`` by ``self``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c", "d"], ["u", "v", "w"])
            sage: fs = F._from_list_([1, 0, 2, 1])
            sage: list(map(fs.getimage, ["a", "b", "c", "d"]))
            ['v', 'u', 'w', 'v']
        """
        parent = self._parent
        return parent._unrank_codomain(self._getitem(parent._rank_domain(i)))

    cpdef items(self):
        """
        The items of ``self``.

        Return the list of the couple ``(x, self(x))``

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: F.from_dict({"a": "b", "b": "a", "c": "b"}).items()
            [('a', 'b'), ('b', 'a'), ('c', 'b')]

        TESTS::

            sage: all(F.from_dict(dict(f.items())) == f for f in F)
            True
        """
        parent = self._parent
        return [(parent._unrank_domain(i),
                 parent._unrank_codomain(self._getitem(i)))
                for i in range(parent._m)]

    def _repr_(self):
        """
        TESTS::

            sage: F = FiniteSetMaps(["a", "b"], [3, 4, 5])
            sage: F._from_list_([0, 2])
            map: a -> 3, b -> 5
        """
        return "map: " + ", ".join("%s -> %s" % (i, self(i))
                                   for i in self.domain())


cdef class FiniteSetEndoMap_N(FiniteSetMap_MN):
    """
    Map from ``range(n)`` to itself.

    .. SEEALSO:: :class:`FiniteSetMap_MN` for assumptions on the parent

    TESTS::

        sage: fs = FiniteSetMaps(3)([1, 0, 2])
        sage: TestSuite(fs).run()
    """

    def __mul__(FiniteSetEndoMap_N self, FiniteSetEndoMap_N other):
        """
        TESTS::

            sage: F = FiniteSetMaps(3)
            sage: F([1, 0, 2]) * F([2, 1, 0])
            [2, 0, 1]
            sage: F = FiniteSetMaps(3, action='right')
            sage: F([1, 0, 2]) * F([2, 1, 0])
            [1, 2, 0]
        """
        assert (self._parent is other._parent), "Parent mismatch"
        if self._parent._action == "right":
            return self._compose_internal_(other, self._parent)
        return other._compose_internal_(self, self._parent)

    def __pow__(self, n, dummy):
        """
        Return the ``n``-th power of ``self``.

        INPUT:

        - ``n`` -- positive integer
        - ``dummy`` -- not used; must be set to ``None`` (for compatibility only)

        EXAMPLES::

            sage: fs = FiniteSetMaps(3)([1,0,2])
            sage: fs^2
            [0, 1, 2]
            sage: fs^0
            [0, 1, 2]
            sage: fs.__pow__(2)
            [0, 1, 2]
        """
        if dummy is not None:
            raise RuntimeError("__pow__ dummy argument not used")
        return generic_power(self, n)


cdef class FiniteSetEndoMap_Set(FiniteSetMap_Set):
    """
    Map from a set to itself.

    .. SEEALSO:: :class:`FiniteSetMap_Set` for assumptions on the parent

    TESTS::

        sage: F = FiniteSetMaps(["a", "b", "c"])
        sage: f = F.from_dict({"a": "b", "b": "a", "c": "b"}); f
        map: a -> b, b -> a, c -> b
        sage: TestSuite(f).run()
    """

    def __mul__(FiniteSetEndoMap_Set self, FiniteSetEndoMap_Set other):
        """
        TESTS::

            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: f = F.from_dict({"a": "b", "b": "a", "c": "b"}); f
            map: a -> b, b -> a, c -> b
            sage: g = F.from_dict({"a": "c", "b": "c", "c": "a"}); g
            map: a -> c, b -> c, c -> a
            sage: f * g
            map: a -> b, b -> b, c -> b
            sage: g * f
            map: a -> c, b -> c, c -> c
        """
        assert (self._parent is other._parent), "Parent mismatch"
        if self._parent._action == "right":
            return self._compose_internal_(other, self._parent)
        return other._compose_internal_(self, self._parent)

    def __pow__(self, n, dummy):
        """
        Return the ``n``-th power of ``self``.

        INPUT:

        - ``n`` -- positive integer
        - ``dummy`` -- not used; must be set to None (for compatibility only)

        EXAMPLES::

            sage: F = FiniteSetMaps(["a", "b", "c"])
            sage: f = F.from_dict({"a": "b", "b": "a", "c": "b"}); f
            map: a -> b, b -> a, c -> b
            sage: f^2
            map: a -> a, b -> b, c -> a
            sage: f^3 == f
            True
        """
        if dummy is not None:
            raise RuntimeError("__pow__ dummy argument not used")
        return generic_power(self, n)
