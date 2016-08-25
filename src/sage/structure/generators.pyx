"""
(Algebraic) generators of a Sage object
"""

from sage.structure.sage_object cimport SageObject

class GenIter(SageObject):
    """
    An iterator over a set of generators.
    """
    def __init__(self, Generators gens):
        if gens._index_set is None:
            self._index_iter = iter(xrange(gens.count()))
        else:
            self._index_iter = iter(gens._index_set)
        self._gens = gens

    def next(self):
        return self._gens.get_from_index(next(self._index_iter))

    def __iter__(self):
        return self

cdef class Generators(SageObject):
    """
    This class implements generators that can be attached to CategoryObjects.
    """
    def __init__(self, obj, index_set, category):
        self._obj = obj
        self._index_set = index_set
        self._category = category

    cpdef get_from_index(self, i):
        return self._obj._gen_(i)

    def __contains__(self, x):
        for a in self:
            if x == a:
                return True
        return False

    def __call__(self, i):
        return self.get_from_index(i)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return self.list()[i]
        return self.get_from_index(i)

    def __iter__(self):
        return GenIter(self)

    def __len__(self):
        return self.count()

    cpdef index_set(self):
        if self._index_set is None:
            return range(self.count())
        else:
            return self._index_set

    cpdef category(self):
        return self._category

    cpdef obj(self):
        return self._obj

    cpdef count(self):
        return len(self._index_set)

    cpdef list(self):
        try:
            if not self._index_set.is_finite():
                raise ValueError("index set must be finite to compute list")
        except AttributeError:
            pass
        return [self.get_from_index(i) for i in self.index_set()]

    def _repr_(self):
        if not isinstance(self._obj, Generators):
            return "Set of generators of %s"%(self.obj())
        else:
            raise RuntimeError("Set of generators of a generators object")

cdef class Generators_finite(Generators):
    def __init__(self, obj, n, index_set, category):
        """
        Creates a generator object representing finitely many generators.
        """
        self._n = n
        Generators.__init__(self, obj, index_set, category)

    cpdef count(self):
        return self._n

    def __cmp__(self, other_unty):
        cdef Generators_finite other = other_unty
        return cmp((self._obj, self._n, self._index_set, self._category),
                   (other._obj, other._n, other._index_set, other._category))


    def _repr_(self):
        return "(" + ", ".join([repr(self[i]) for i in self.index_set()]) +  ")"

    def __reduce__(self):
        return make_finite_gens, (self._obj, self._n, self._index_set, self._category)

def make_finite_gens(obj, n, index_set, category):
    return Generators_finite(obj, n, index_set, category)


cdef class Generators_list(Generators_finite):
    """
    This class represents a set of generators as a tuple of elements, indexed
    by the integers 0 to len(self)-1.

    It is the easiest to use of all the generator classes, and gets constructed
    implicitly when a list is passed into \code{_populate_generators_}.
    """
    def __init__(self, obj, L, category):
        """
        EXAMPLES:
            sage: from sage.structure.generators import Generators_list
            sage: gens = Generators_list(ZZ, [2,3], Rings)
            sage: gens.count()
            2
            sage: gens(1)
            3
            sage: list(gens)
            [2, 3]
        """
        self._List = tuple(L)
        Generators_finite.__init__(self, obj, len(self._List), None, category)

    cpdef get_from_index(self, i):
        """
        EXAMPLES:
            sage: from sage.structure.generators import Generators_list
            sage: gens = Generators_list(ZZ, [5,9], Rings)
            sage: gens.get_from_index(0)
            5
        """
        try:
            return self._List[i]
        except IndexError:
            raise IndexError("No such generator: %s must be less than %s" % (i, len(self._List)))

    def __iter__(self):
        """
        EXAMPLES:
            sage: from sage.structure.generators import Generators_list
            sage: gens = Generators_list(ZZ, [3,4,5], Rings)
            sage: [x^2 for x in gens]
            [9, 16, 25]
        """
        return iter(self._List)

    cpdef list(self):
        """
        Actually returns a tuple, which is immutable.

        EXAMPLES:
            sage: from sage.structure.generators import Generators_list
            sage: gens = Generators_list(ZZ, (1,2,3), Rings)
            sage: gens.list()
            (1, 2, 3)
        """
        return self._List

    def __reduce__(self):
        return make_list_gens, (self._obj, self._List, self._category)

    def __cmp__(left, _right):
        cdef Generators_list right
        try:
            right = _right
            return cmp(left._n, right._n) or \
                   cmp(left._List, right._List) or \
                   cmp(left._category, right._category)
        except TypeError:
            return cmp(type(left), type(right))

    def _repr_(self):
        return repr(self._List)

def make_list_gens(*args):
    """
    TEST:
        sage: from sage.structure.generators import Generators_list
        sage: gens = Generators_list(QQ, [-3,17], Fields())
        sage: loads(dumps(gens)) == gens
        True
    """
    return Generators_list(*args)

cdef class Generators_naturals(Generators):
    def __init__(self, obj, category):
        import sage.rings.natural_parents
        Generators.__init__(self, obj, sage.rings.natural_parents.NaturalSemiring, category)

    cpdef count(self):
        import sage.rings.infinity
        return sage.rings.infinity.infinity

    def __getitem__(self, i):
        if isinstance(i, slice):
            raise NotImplementedError
        return self.get_from_index(i)

cdef class Generators_none(Generators):
    def __init__(self, obj, category = None):
        Generators.__init__(self, obj, None, category)

    cpdef get_from_index(self, i):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")

    def __contains__(self, x):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")

    def __iter__(self):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")

    cpdef index_set(self):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")

    cpdef count(self):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")

    cpdef list(self):
        if self._category is None:
            raise ValueError("this object does not have generators")
        else:
            raise ValueError("this object does not have generators in that category")



cdef class Generators_lazy_all(Generators):
    """
    Use this generators class if there is a finite list of generators that is expensive to compute but which is computed all at once.
    """
    def __init__(self, G, index_set, category, computing_function):
        self._f = computing_function
        Generators.__init__(self, G, index_set, category)

    cdef int _compute_gens(self) except -1:
        self._gens = self._f()

    cpdef get_from_index(self, i):
        try:
            return self._gens[i]
        except AttributeError:
            self._compute_gens()
        return self._gens[i]

    cpdef count(self):
        try:
            return len(self._gens)
        except AttributeError:
            self._compute_gens()
        return len(self._gens)

    cpdef list(self):
        try:
            return list(self._gens)
        except AttributeError:
            self._compute_gens()
        return list(self._gens)

    def _repr_(self):
        try:
            return repr(self._gens)
        except AttributeError:
            self._compute_gens()
        return repr(self._gens)

