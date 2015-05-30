"""
Low level part of the interface to Fokko Ducloux's Coxeter 3 library

.. TODO::

    - Write a more efficient method for converting polynomials in
      Coxeter to Sage polynomials.
"""
#*****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"
include "decl.pxi"

initConstants()

from sage.rings.all import Integer, ZZ

cdef class String:
    def __cinit__(self, s=""):
        """
        Construct a Coxeter string from a Python string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String("hello"); s                              # optional - coxeter3
            hello
        """
        String_construct_str(&self.x, s)

    def __dealloc__(self):
        """
        Deallocate the memory for this string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String("hello")                                 # optional - coxeter3
            sage: del s                                               # optional - coxeter3
        """
        String_destruct(&self.x)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String('Hi')                                    # optional - coxeter3
            sage: s                                                   # optional - coxeter3
            Hi
        """
        return self.x.ptr()

    def __hash__(self):
        """
        Return the hash of this String

        This is the hash of the tuple consisting of the class name and
        the name of this type.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String        # optional - coxeter3
            sage: s = String('hello')                                  # optional - coxeter3
            sage: hash(s) == hash('hello')                             # optional - coxeter3
            True
        """
        return hash(repr(self))

    def __richcmp__(String self, other, int op):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String        # optional - coxeter3
            sage: ta1 = String('A')                                    # optional - coxeter3
            sage: ta2 = String('A')                                    # optional - coxeter3
            sage: tb = String('b')                                     # optional - coxeter3
            sage: ta1 == ta2                                           # optional - coxeter3
            True
            sage: tb != ta1                                            # optional - coxeter3
            True
            sage: all([ta1 < tb, ta1 <= tb, ta1 <= ta1])               # optional - coxeter3
            True
            sage: all([tb > ta1, tb >= ta1, tb >= tb])                 # optional - coxeter3
            True
        """
        if type(other) != type(self):
            return False

        s = repr(self)
        o = repr(other)

        if op == 2: # ==
            return s == o
        elif op == 3: # !=
            return s != o
        elif op == 0: # <
            return s < o
        elif op == 1: # <=
            return s <= o
        elif op == 4: # >
            return s > o
        elif op == 5: # >=
            return s >= o

    def __len__(self):
        """
        Return the length of this string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String('Hi')                                    # optional - coxeter3
            sage: len(s)                                              # optional - coxeter3
            2
        """
        return self.x.length()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String('Hi')                                    # optional - coxeter3
            sage: TestSuite(s).run()                                  # optional - coxeter3
        """
        return (String, (repr(self),) )

cdef class Type:
    def __cinit__(self, s):
        """
        Construct a Coxeter Type from a Python string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A'); t                                    # optional - coxeter3
            A
        """
        Type_construct_str(&self.x, s)

    def __dealloc__(self):
        """
        Deallocate the memory for this Type.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A')                                       # optional - coxeter3
            sage: del t                                               # optional - coxeter3
        """
        Type_destruct(&self.x)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A'); t                                    # optional - coxeter3
            A
        """
        return self.x.name().ptr()

    def name(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A')                                       # optional - coxeter3
            sage: t.name()                                            # optional - coxeter3
            A
        """
        return String(self.x.name().ptr())

    def __hash__(self):
        """
        Return the hash of this Type.

        This is the hash of the tuple consisting of the class name and
        the name of this type.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type          # optional - coxeter3
            sage: a = Type('A')                                        # optional - coxeter3
            sage: b = Type('B')                                        # optional - coxeter3
            sage: hash(a) == hash(b)                                   # optional - coxeter3
            False
            sage: d = {a: 1, b: 2}                                     # optional - coxeter3
        """
        return hash(('Type', self.name()))

    def __richcmp__(Type self, other, int op):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type          # optional - coxeter3
            sage: ta1 = Type('A')                                      # optional - coxeter3
            sage: ta2 = Type('A')                                      # optional - coxeter3
            sage: tb = Type('b')                                       # optional - coxeter3
            sage: ta1 == ta2                                           # optional - coxeter3
            True
            sage: tb != ta1                                            # optional - coxeter3
            True
            sage: all([ta1 < tb, ta1 <= tb, ta1 <= ta1])               # optional - coxeter3
            True
            sage: all([tb > ta1, tb >= ta1, tb >= tb])                 # optional - coxeter3
            True
        """
        if type(other) != type(self):
            return False

        s = repr(self)
        o = repr(other)

        if op == 2: # ==
            return s == o
        elif op == 3: # !=
            return s != o
        elif op == 0: # <
            return s < o
        elif op == 1: # <=
            return s <= o
        elif op == 4: # >
            return s > o
        elif op == 5: # >=
            return s >= o

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type           # optional - coxeter3
            sage: t = Type('A')                                         # optional - coxeter3
            sage: TestSuite(t).run()                                    # optional - coxeter3
        """
        return (Type, (repr(self), ))

cdef class CoxGroup(SageObject):
    def __cinit__(self, cartan_type):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup   # optional - coxeter3
            sage: W = CoxGroup(['A', 5]); W                                         # optional - coxeter3
            Coxeter group of type A and rank 5

        Coxeter 3 segfault's on the trivial Coxeter group; so we catch
        this and raise a not implemented error::

            sage: W = CoxGroup(['A', 0]); W                                         # optional - coxeter3
            Traceback (most recent call last):
            ...
            NotImplementedError: Coxeter group of type ['A',0] using Coxeter 3 not yet implemented
        """
        from sage.combinat.root_system.all import CartanType, coxeter_matrix
        self.cartan_type = CartanType(cartan_type)
        ordering = self._ordering_from_cartan_type(self.cartan_type)

        if len(cartan_type) == 2:
            type, rank = cartan_type
        else:
            type, rank, affine = cartan_type
            if affine != 1:
                raise NotImplementedError

            type = type.lower()
            rank = rank + 1

        type = 'B' if type == 'C' else type

        if rank == 0:
            raise NotImplementedError("Coxeter group of type ['A',0] using Coxeter 3 not yet implemented")
        cdef Type t = Type(type)
        cdef c_CoxGroup* c_W = coxeterGroup(t.x, rank)
        self.x = c_W
        self.out_ordering = dict(zip(range(1, rank+1), ordering))
        self.in_ordering = dict([(b,a) for a,b in self.out_ordering.items()])

        # Check that the Coxeter matrices match up.
        if self.coxeter_matrix() != coxeter_matrix(self.cartan_type):
            print "Warning, differing Coxeter matrices"

    @classmethod
    def _ordering_from_cartan_type(cls, cartan_type):
        """
        Return an ordering of the index set associated to the Cartan type.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W._ordering_from_cartan_type(CartanType(['A',5]))                     # optional - coxeter3
            [1, 2, 3, 4, 5]
        """
        from sage.misc.all import srange
        from sage.rings.all import Integer
        t = cartan_type.type()
        r = cartan_type.rank()
        is_affine = cartan_type.is_affine()

        if t in ['B', 'C', 'D', 'F', 'H']:
            return srange(r-1 if is_affine else r,
                          -1 if is_affine else 0, -1)
        elif t in ['A', 'I']:
            return srange(0 if is_affine else 1, r+1)
        elif t in ['G']:
            if is_affine:
                raise NotImplementedError
            else:
                return map(Integer, [1, 2])
        elif t in ['E']:
            if is_affine:
                return srange(1, r) + [Integer(0)]
            else:
                return srange(1, r+1)
        else:
            raise NotImplementedError

    def __hash__(self):
        """
        Return the hash of this CoxGroup.

        This is the hash of the tuple of the class's name, the type,
        and the rank.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: A4 = CoxGroup(['A', 4])                                               # optional - coxeter3
            sage: d = {A4: True}                                                        # optional - coxeter3
        """
        return hash((self.__class__.__name__, self.type(), self.rank()))

    def __richcmp__(CoxGroup self, other, int op):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: A4 = CoxGroup(['A', 4])                                               # optional - coxeter3
            sage: A5 = CoxGroup(['A', 5])                                               # optional - coxeter3
            sage: B4 = CoxGroup(['B', 4])                                               # optional - coxeter3
            sage: A4 == A4                                                              # optional - coxeter3
            True
            sage: A4 != B4                                                              # optional - coxeter3
            True
            sage: A4 < B4                                                               # optional - coxeter3
            True
            sage: A5 > A4                                                               # optional - coxeter3
            True
            sage: A4 >= A4                                                              # optional - coxeter3
            True
            sage: B4 >= A5                                                              # optional - coxeter3
            True
        """
        if type(other) != type(self):
            return False

        s_t = self.type()
        o_t = other.type()
        s_r = self.rank()
        o_r = other.rank()

        if op == 2: # ==
            return s_t == o_t and s_r == o_r
        elif op == 3: # !=
            return s_t != o_t or s_r != o_r
        elif op == 0: # <
            return s_t < o_t or (s_t == o_t and s_r < o_r)
        elif op == 1: # <=
            return s_t < o_t or (s_t == o_t and s_r <= o_r)
        elif op == 4: # >
            return s_t > o_t or (s_t == o_t and s_r > o_r)
        elif op == 5: # >=
            return s_t > o_t or (s_t == o_t and s_r >= o_r)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                               # optional - coxeter3
            sage: TestSuite((W)).run()                                                 # optional - coxeter3
        """
        return (CoxGroup, (self.cartan_type,))

    def __dealloc__(self):
        """
        Deallocate the memory for this CoxGroup.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                               # optional - coxeter3
            sage: del W                                                                # optional - coxeter3
        """
        CoxGroup_delete(self.x)

    def __repr__(self):
        """
        Return a string representation of this Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 5]); W                                            # optional - coxeter3
            Coxeter group of type A and rank 5
        """
        return "Coxeter group of type %s and rank %s"%(self.type(), self.rank())

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                               # optional - coxeter3
            sage: list(iter(W))                                                        # optional - coxeter3
            [[], [1], [2], [1, 2], [2, 1], [1, 2, 1]]
        """
        return CoxGroupIterator(self)

    def bruhat_interval(self, w, v):
        """
        Return the list of the elements in the Bruhat interval between `w` and `v`.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                               # optional - coxeter3
            sage: W.bruhat_interval([], [1,2])                                         # optional - coxeter3
            [[], [1], [2], [1, 2]]
        """
        cdef CoxGroupElement ww = CoxGroupElement(self, w)
        cdef CoxGroupElement vv = CoxGroupElement(self, v)
        cdef c_List_CoxWord l = c_List_CoxWord_factory(0)
        interval(l, self.x[0], ww.word, vv.word)
        bruhat_interval = []
        cdef int j = 0
        cdef CoxGroupElement u
        cdef CoxGroupElement gg = CoxGroupElement(self, [])
        for j from 0 <= j < l.size():
            u = gg._new()
            u.word = l.get_index(j)
            bruhat_interval.append(u)

        # This destruction most likely does not be belong there, and
        # it causes a segfault. See discussion on #12912.
        # List_CoxWord_destruct(&l)

        return bruhat_interval

    def orderings(self):
        """
        Return two dictionaries specifying the mapping of the labels
        of the Dynkin diagram between Sage and Coxeter3 and the its
        inverse.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.orderings()                                                         # optional - coxeter3
            ({1: 1, 2: 2, 3: 3, 4: 4, 5: 5}, {1: 1, 2: 2, 3: 3, 4: 4, 5: 5})
        """
        return self.in_ordering, self.out_ordering

    def type(self):
        """
        Return the type of this Coxeter group.

        Note that the type does not include the rank.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.type()                                                              # optional - coxeter3
            A
        """
        return Type(self.x.type().name().ptr())

    def rank(self):
        """
        Return the rank of this Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.rank()                                                              # optional - coxeter3
            5
        """
        return self.x.rank()

    def order(self):
        """
        Return the order of this Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.order()                                                             # optional - coxeter3
            720
            sage: W = CoxGroup(['A', 3, 1])                                             # optional - coxeter3
            sage: W.order()                                                             # optional - coxeter3
            +Infinity
        """
        if self.is_finite():
            return Integer(self.x.order())
        else:
            from sage.all import infinity
            return infinity

    def is_finite(self):
        """
        Return whether this Coxeter group is finite.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.is_finite()                                                         # optional - coxeter3
            True
            sage: W = CoxGroup(['A', 3, 1])                                             # optional - coxeter3
            sage: W.is_finite()                                                         # optional - coxeter3
            False
        """
        return isFiniteType(self.x)

    cpdef full_context(self):
        """
        Make all of the elements of a finite Coxeter group available.

        Raises an error if W is not finite

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                # optional - coxeter3
            sage: W.full_context()                                                      # optional - coxeter3
            sage: W = CoxGroup(['A', 2,1])                                              # optional - coxeter3
            sage: W.full_context()                                                      # optional - coxeter3
            Traceback (most recent call last):
            ...
            TypeError: Group needs to be finite.
        """
        if not self.is_finite():
            raise TypeError, "Group needs to be finite."
        cdef c_FiniteCoxGroup* fcoxgroup = <c_FiniteCoxGroup*>(self.x)
        if not fcoxgroup.isFullContext():
            fcoxgroup.fullContext()


    def long_element(self):
        """
        Return the longest word in a finite Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                # optional - coxeter3
            sage: W.long_element()                                                      # optional - coxeter3
            [1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1]

            sage: W = CoxGroup(['A', 3, 1])                                             # optional - coxeter3
            sage: W.long_element()                                                      # optional - coxeter3
            Traceback (most recent call last):
            ...
            TypeError: Group needs to be finite.
        """
        self.full_context()
        cdef c_FiniteCoxGroup* fcoxgroup = <c_FiniteCoxGroup*>(self.x)
        cdef CoxGroupElement w0 = CoxGroupElement(self, [])
        w0.word = fcoxgroup.longest_coxword()
        return w0

    def __call__(self, w):
        """
        Return a reduced expression for `w`.

        INPUT:

        - ``w`` -- a word for an element of ``self``, not necessarily reduced

        OUTPUT:

        - a reduced expression for ``w``

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                               # optional - coxeter3
            sage: w = [1,1,3,5,4,5,4]                                                  # optional - coxeter3
            sage: W.__call__(w)                                                        # optional - coxeter3
            [3, 4, 5]
        """
        return CoxGroupElement(self, w).reduced()

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix for this Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup      # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                               # optional - coxeter3
            sage: W.coxeter_matrix()                                                   # optional - coxeter3
            [1 3 2 2 2]
            [3 1 3 2 2]
            [2 3 1 3 2]
            [2 2 3 1 3]
            [2 2 2 3 1]

        """
        from sage.all import matrix, ZZ
        rank = self.rank()
        m = matrix(ZZ, rank, rank)
        for i, ii in enumerate(self.cartan_type.index_set()):
            ii = self.in_ordering[ii]-1
            for j, jj in enumerate(self.cartan_type.index_set()):
                jj = self.in_ordering[jj]-1
                m[i,j] = self.x.M(ii, jj)
        return m

    def coxeter_graph(self):
        """
        Return the Coxeter graph for this Coxeter group.

        OUTPUT:: a Sage graph

        .. NOTE::

           This uses the labels native to Coxeter3. This is useful
           when trying to obtain the mapping between the labels of
           Sage and Coxeter3.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup# optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                         # optional - coxeter3
            sage: W.coxeter_graph()                                              # optional - coxeter3
            Graph on 5 vertices
            sage: sorted(W.coxeter_graph().edges())                              # optional - coxeter3
            [(1, 2, None), (2, 3, None), (3, 4, None), (4, 5, None)]
        """
        from sage.all import Graph
        g = Graph()
        m = self.coxeter_matrix()
        rank = self.rank()
        for i, row in enumerate(m.rows()):
            for j in range(i+1,rank):
                if row[j] == 3:
                    g.add_edge(i+1, j+1)
                elif row[j] > 4:
                    g.add_edge(i+1, j+1, row[j])
        return g



cdef class CoxGroupElement:
    def __init__(self, CoxGroup group, w, normal_form=True):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [2,1,2,1,1], normal_form=False); w                         # optional - coxeter3
            [2, 1, 2, 1, 1]
            sage: w = CoxGroupElement(W, [1,1,4,5,4], normal_form=False); w                         # optional - coxeter3
            [1, 1, 4, 5, 4]
            sage: w = CoxGroupElement(W, [1,1,4,5,4]); w                                            # optional - coxeter3
            [4, 5, 4]
        """
        self.group = (<CoxGroup>group).x
        self._parent = group
        self.word.reset()
        for i in w:
            self.word.append_letter(self._parent.in_ordering[i])

        if normal_form:
            self.group.normalForm(self.word)


    def __cinit__(self):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: CoxGroupElement(W, [1,2,3,2,3])                                                   # optional - coxeter3
            [1, 3, 2]
        """
        CoxWord_construct(&self.word)

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [1,2,3,2,3])                                               # optional - coxeter3
            sage: del w                                                                             # optional - coxeter3
        """
        CoxWord_destruct(&self.word)

    def _coxnumber(self):
        """
        Return the internal integer used by Coxeter3 to represent this element.

        TESTS::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [1,2,3,2,3])                                               # optional - coxeter3
            sage: w._coxnumber()                                                                    # optional - coxeter3
            7L
        """
        return self.group.extendContext(self.word)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: TestSuite(w).run()                                                                # optional - coxeter3
        """
        return (CoxGroupElement, (self._parent, list(self)))

    def __invert__(self):
        """
        Return the inverse of this element.

        EXAMPLES::


            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: ~w                                                                                # optional - coxeter3
            [3, 2, 1]
        """
        return CoxGroupElement(self._parent, reversed(self))

    inverse = __invert__

    cpdef CoxGroup parent(self):
        """
        Return the parent Coxeter group for this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: w.parent()                                                                        # optional - coxeter3
            Coxeter group of type A and rank 5

        """
        return self._parent

    def __getitem__(self, i):
        """
        Return the `i^{th}` entry of this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: w[0]                                                                              # optional - coxeter3
            1
            sage: w[2]                                                                              # optional - coxeter3
            3
            sage: w[:-2]                                                                            # optional - coxeter3
            [1]
            sage: w[-2:]                                                                            # optional - coxeter3
            [2, 3]
            sage: w[3:0:-1]                                                                         # optional - coxeter3
            [3, 2]
            sage: w[4]                                                                              # optional - coxeter3
            Traceback (most recent call last):
            ...
            IndexError: The index (4) is out of range.
        """
        if isinstance(i, slice):
            #Get the start, stop, and step from the slice
            return [self[ii] for ii in xrange(*i.indices(len(self)))]
        if i < 0:
            i += len(self)
        if i >= len(self):
            raise IndexError, "The index (%d) is out of range."%i

        return self._parent.out_ordering[self.word.get_index(i)]

    def __repr__(self):
        """
        Return a string representation of this CoxGroupElement as a list of generators.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                             # optional - coxeter3
            sage: w = W([1,2,3]); w                                 # optional - coxeter3
            [1, 2, 3]

        """
        return repr(list(self))

    def __hash__(self):
        """
        Return the hash of this element.

        This is a hash of the tuple of the class name, the parent, and
        a tuple of the reduced word.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *         # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                           # optional - coxeter3
            sage: w = W([1,2,3])                                   # optional - coxeter3
            sage: v = W([2,3,4])                                   # optional - coxeter3
            sage: hash(w) == hash(v)                               # optional - coxeter3
            False
        """
        return hash((self.__class__.__name__, self.parent(), tuple(self)))

    def __richcmp__(CoxGroupElement self, other, int op):
        """
        EXAMPLES:

            sage: from sage.libs.coxeter3.coxeter import *        # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                          # optional - coxeter3
            sage: V = CoxGroup(['A', 6])                          # optional - coxeter3
            sage: w1 = W([1,2,3])                                 # optional - coxeter3
            sage: w2 = W([2,3,4])                                 # optional - coxeter3
            sage: v1 = V([1,2,3])                                 # optional - coxeter3
            sage: w1 == w1                                        # optional - coxeter3
            True
            sage: w1 != w2                                        # optional - coxeter3
            True
            sage: all([w1 < w2, w1 <= w2, w1 <= w1])              # optional - coxeter3
            True
            sage: all([w2 > w1, w2 >= w1, w2 >= w2])              # optional - coxeter3
            True
            sage: w1 == v1                                        # optional - coxeter3
            False
            sage: w1 != v1                                        # optional - coxeter3
            True
        """
        if type(other) != type(self):
            return False

        s_p = self.parent()
        o_p = other.parent()
        s_l = list(self)
        o_l = list(other)

        if op == 2: # ==
            return s_p == o_p and s_l == o_l
        elif op == 3: # !=
            return s_p != o_p or s_l != o_l
        elif op == 0: # <
            return s_p < o_p or (s_p == o_p and s_l < o_l)
        elif op == 1: # <=
            return s_p < o_p or (s_p == o_p and s_l <= o_l)
        elif op == 4: # >
            return s_p > o_p or (s_p == o_p and s_l > o_l)
        elif op == 5: # >=
            return s_p > o_p or (s_p == o_p and s_l >= o_l)


    def __iter__(self):
        """
        Return an iterator for the letters in the reduced word for this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *      # optional - coxeter3
            sage: W = CoxGroup(['A',5])                         # optional - coxeter3
            sage: w = W([1,2,3])                                # optional - coxeter3
            sage: [a for a in w]                                # optional - coxeter3
            [1, 2, 3]
        """
        return (self[i] for i in xrange(len(self)))

    def __len__(self):
        """
        Return the length of this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *       # optional - coxeter3
            sage: W = CoxGroup(['A',5])                          # optional - coxeter3
            sage: w = W([1,2,3])                                 # optional - coxeter3
            sage: len(w)                                         # optional - coxeter3
            3
        """
        return self.word.length()

    def left_descents(self):
        """
        Return the left descent set of this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *      # optional - coxeter3
            sage: W = CoxGroup(['A',5])                         # optional - coxeter3
            sage: w = W([1,2,1])                                # optional - coxeter3
            sage: w.left_descents()                             # optional - coxeter3
            [1, 2]
        """
        return LFlags_to_list(self._parent, self.group.ldescent(self.word))

    def right_descents(self):
        """
        Return the right descent set of this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *      # optional - coxeter3
            sage: W = CoxGroup(['A',5])                         # optional - coxeter3
            sage: w = W([1,2,1])                                # optional - coxeter3
            sage: w.right_descents()                            # optional - coxeter3
            [1, 2]
        """
        return LFlags_to_list(self._parent, self.group.rdescent(self.word))

    def bruhat_le(self, w):
        """
        Return whether u = (self) is less than w in Bruhat order.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *       # optional - coxeter3
            sage: W = CoxGroup(['A',5])                          # optional - coxeter3
            sage: w = W([1,2,3,4,5,4])                           # optional - coxeter3
            sage: v = W([1,2,4,5,4])                             # optional - coxeter3
            sage: v.bruhat_le(w)                                 # optional - coxeter3
            True
            sage: w.bruhat_le(w)                                 # optional - coxeter3
            True
            sage: w.bruhat_le(v)                                 # optional - coxeter3
            False
        """
        cdef CoxGroupElement ww = CoxGroupElement(self._parent, w)
        return self.group.inOrder_word(self.word, ww.word)

    def is_two_sided_descent(self, s):
        """
        Return whether ``s`` is a two sided descent of ``self``.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *       # optional - coxeter3
            sage: W = CoxGroup(['A',2])                          # optional - coxeter3
            sage: x = W([1,2,1])                                 # optional - coxeter3
            sage: x.is_two_sided_descent(1)                      # optional - coxeter3
            True
        """
        cdef Generator ss = self._parent.in_ordering[s]
        return self.group.isDescent(self.word, s)

    cdef CoxGroupElement _new(self):
        """
        Return a new copy of this element.
        """
        cdef CoxGroupElement res = CoxGroupElement(self.parent(), [])
        res.word.set(self.word)
        return res

    def coatoms(self):
        """
        Return the coatoms of this element in Bruhat order.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *          # optional - coxeter3
            sage: W = CoxGroup(['A',2])                             # optional - coxeter3
            sage: W([1,2,1]).coatoms()                              # optional - coxeter3
            [[2, 1], [1, 2]]
            sage: W([]).coatoms()                                   # optional - coxeter3
            []
        """
        cdef c_List_CoxWord list = c_List_CoxWord_factory(0)
        self.group.coatoms(list, self.word)

        coatoms = []

        cdef Length i = 0
        cdef CoxGroupElement res
        for i from 0 <= i < list.size():
            res = self._new()
            res.word = list.get_index(i)
            coatoms.append(res)
        return coatoms

    def normal_form(self):
        """
        Return ``self`` in normal form.

        This is the lexicographically minimal reduced word for
        ``self``.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [2,1,2], normal_form=False); w                             # optional - coxeter3
            [2, 1, 2]
            sage: w.normal_form()                                                                   # optional - coxeter3
            [1, 2, 1]

        """
        cdef CoxGroupElement res = self._new()
        self.group.normalForm(res.word)
        return res

    def reduced(self):
        """
        Return a reduced word for this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [2,1,2,1,1], normal_form=False); w                         # optional - coxeter3
            [2, 1, 2, 1, 1]
            sage: w.reduced()                                                                       # optional - coxeter3
            [1, 2, 1]
        """
        cdef CoxGroupElement res = self._new()
        self.group.reduced(res.word, self.word)
        return res

    def __mul__(CoxGroupElement self, CoxGroupElement y):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup                    # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                             # optional - coxeter3
            sage: W([1]) * W([1])                                                                    # optional - coxeter3
            []
            sage: W([1,2]) * W([1])                                                                  # optional - coxeter3
            [1, 2, 1]
        """
        cdef CoxGroupElement res = self._new()
        self.group.prod(res.word, y.word)
        return res

    def poincare_polynomial(self):
        """
        Return the Poincare polynomial associated with the Bruhat
        interval between the identity element and this one.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup                    # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                             # optional - coxeter3
            sage: W([]).poincare_polynomial()                                                        # optional - coxeter3
            1
            sage: W([1,2,1]).poincare_polynomial()                                                   # optional - coxeter3
            t^3 + 2*t^2 + 2*t + 1
        """
        cdef CoxGroup W = self.parent()
        cdef c_List_CoxWord result = c_List_CoxWord_factory(0)
        cdef CoxGroupElement id = CoxGroupElement(W, [])
        cdef CoxGroupElement ww = CoxGroupElement(W, self)
        interval(result, W.x[0], id.word, ww.word)

        cdef int j = 0
        cdef list coefficients = [0]*(len(ww)+1)
        for j from 0 <= j < result.size():
            coefficients[result.get_index(j).length()] += 1
        return ZZ['t'](coefficients)


    def kazhdan_lusztig_polynomial(self, v):
        """
        Return the Kazhdan-Lusztig polynomial `P_{u,v}` where `u` is this element.

        Currently this is a bit inefficient as it constructs the
        polynomial from its string representation.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup                    # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                             # optional - coxeter3
            sage: W([]).kazhdan_lusztig_polynomial([1,2,1])                                          # optional - coxeter3
            1
            sage: W([1,2,1]).kazhdan_lusztig_polynomial([])                                          # optional - coxeter3
            0
        """
        from sage.all import ZZ
        cdef CoxGroupElement vv
        if not isinstance(v, CoxGroupElement):
            vv = CoxGroupElement(self._parent, v)
        else:
            vv = v

        ZZq = ZZ['q']
        if not self.group.inOrder_word(self.word, vv.word):
            return ZZq(0)

        cdef CoxNbr x = self.group.extendContext(self.word)
        cdef CoxNbr y = self.group.extendContext(vv.word)
        cdef c_KLPol kl_poly = self.group.klPol(x, y)

        cdef String s = String()
        klpoly_append(s.x, kl_poly, "q")
        return ZZq(str(s))

    def mu_coefficient(self, v):
        r"""
        Return the mu coefficient `\mu(u,v)` where `u` is this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                             # optional - coxeter3
            sage: w = W([1,2,3,4,5,4])                              # optional - coxeter3
            sage: v = W([1,2,4,5,4])                                # optional - coxeter3
            sage: w.mu_coefficient(v)                               # optional - coxeter3
            0
            sage: w.mu_coefficient(w)                               # optional - coxeter3
            0
            sage: v.mu_coefficient(w)                               # optional - coxeter3
            1
        """
        from sage.all import ZZ
        cdef CoxGroupElement vv = CoxGroupElement(self._parent, v)
        cdef CoxNbr x = self.group.extendContext(self.word)
        cdef CoxNbr y = self.group.extendContext(vv.word)
        return ZZ(self.group.mu(x,y))

cdef LFlags_to_list(CoxGroup parent, LFlags f):
    """
    Return the right descent set of this element.

    EXAMPLES::

        sage: from sage.libs.coxeter3.coxeter import *         # optional - coxeter3
        sage: W = CoxGroup(['A',5])                            # optional - coxeter3
        sage: w = W([1,2,1])                                   # optional - coxeter3
        sage: w.right_descents()                               # optional - coxeter3
        [1, 2]
    """
    cdef Generator s
    cdef LFlags f1 = f
    l = []
    while f1:
        s = firstBit(f1)
        l.append(parent.out_ordering[s+1])
        f1 = f1 & (f1-1)
    return l

class CoxGroupIterator(object):
    def __init__(self, group):
        """
        A class used to iterate over all of the elements of a Coxeter group.

        .. note::

           This will construct all of the elements of the group within
           Coxeter3.  For some groups, this may be too large to fit
           into memory.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                            # optional - coxeter3
            sage: it = CoxGroupIterator(W)                                                          # optional - coxeter3
            sage: [it.next() for i in range(W.order())]                                             # optional - coxeter3
            [[], [1], [2], [1, 2], [2, 1], [1, 2, 1]]
        """
        self.group = group
        self.order = group.order()
        self.n = 0
        self.group.full_context()

    def __iter__(self):
        """
        Return self, as per the iterator protocol.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                            # optional - coxeter3
            sage: it = iter(W)                                                                      # optional - coxeter3
            sage: it is iter(it)                                                                    # optional - coxeter3
            True
        """
        return self

    def next(self):
        """
        Return the next element in the associated Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                            # optional - coxeter3
            sage: it = CoxGroupIterator(W)                                                          # optional - coxeter3
            sage: it.next()                                                                         # optional - coxeter3
            []
        """
        if self.n >= self.order:
            raise StopIteration
        cdef CoxGroupElement w = self.group([])

        (<CoxGroup>self.group).x.prod_nbr(w.word, self.n)
        self.n += 1
        return w

CoxGroup_cache = {}
def get_CoxGroup(cartan_type):
    """
    TESTS::

        sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator  # optional - coxeter3
        sage: W = CoxGroup(['A', 2])                                                             # optional - coxeter3
    """
    from sage.all import CartanType
    cartan_type = CartanType(cartan_type)
    if cartan_type not in CoxGroup_cache:
        CoxGroup_cache[cartan_type] = CoxGroup(cartan_type)
    return CoxGroup_cache[cartan_type]
