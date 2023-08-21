# distutils: language = c++
# distutils: libraries = coxeter3
# sage_setup: distribution = sagemath-coxeter3

"""
Low level part of the interface to Fokko Ducloux's Coxeter 3 library

.. TODO::

    - Write a more efficient method for converting polynomials in
      Coxeter to Sage polynomials.
"""
# ****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .decl cimport *
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from sage.cpython.string cimport str_to_bytes, bytes_to_str

initConstants()

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


cdef class String:
    def __init__(self, s=""):
        """
        Construct a Coxeter string from a Python string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String("hello"); s                              # optional - coxeter3
            hello
            sage: del s                                               # optional - coxeter3
        """
        self.x = c_String(str_to_bytes(s))

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import String       # optional - coxeter3
            sage: s = String('Hi')                                    # optional - coxeter3
            sage: s                                                   # optional - coxeter3
            Hi
        """
        return bytes_to_str(self.x.ptr())

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
        if type(other) is not type(self):
            if op in (Py_LT, Py_LE, Py_GT, Py_GE):
                return NotImplemented
            return op == Py_NE

        s = repr(self)
        o = repr(other)

        if op == Py_EQ:
            return s == o
        elif op == Py_NE:
            return s != o
        elif op == Py_LT:
            return s < o
        elif op == Py_LE:
            return s <= o
        elif op == Py_GT:
            return s > o
        elif op == Py_GE:
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
    def __init__(self, s):
        """
        Construct a Coxeter Type from a Python string.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A'); t                                    # optional - coxeter3
            A
            sage: del t                                               # optional - coxeter3
        """
        self.x = c_Type(str_to_bytes(s))

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A'); t                                    # optional - coxeter3
            A
        """
        return bytes_to_str(self.x.name().ptr())

    def name(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import Type         # optional - coxeter3
            sage: t = Type('A')                                       # optional - coxeter3
            sage: t.name()                                            # optional - coxeter3
            A
        """
        return String(bytes_to_str(self.x.name().ptr()))

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
        if type(other) is not type(self):
            if op in (Py_LT, Py_LE, Py_GT, Py_GE):
                return NotImplemented
            return op == Py_NE

        s = repr(self)
        o = repr(other)

        if op == Py_EQ:
            return s == o
        elif op == Py_NE:
            return s != o
        elif op == Py_LT:
            return s < o
        elif op == Py_LE:
            return s <= o
        elif op == Py_GT:
            return s > o
        elif op == Py_GE:
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

        Successfully initializes from a relabeled Cartan type::

            sage: ctype = CartanType(['B', 3]).relabel({1: 3, 2: 2, 3: 1})
            sage: W = CoxGroup(ctype)                                               # optional - coxeter3
            sage: CoxeterMatrix(W.coxeter_matrix(), ctype.index_set()) == CoxeterMatrix(ctype) # optional - coxeter3
            True
        """
        from sage.combinat.root_system.cartan_type import CartanType
        from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
        self.cartan_type = CartanType(cartan_type)
        ordering = self._ordering_from_cartan_type(self.cartan_type)

        type, rank = self.cartan_type.type(), self.cartan_type.rank()
        if self.cartan_type.is_affine():
            # Only untwisted affine groups are supported
            try:
                if not self.cartan_type.is_untwisted_affine():
                    raise NotImplementedError('twisted affine groups are not supported in coxeter3')
            except AttributeError:
                pass
            type = type.lower()

        type = 'B' if type == 'C' else type

        if rank == 0:
            raise NotImplementedError("Coxeter group of type ['A',0] using Coxeter 3 not yet implemented")
        cdef Type t = Type(type)
        cdef c_CoxGroup* c_W = coxeterGroup(t.x, rank)
        self.x = c_W
        self.out_ordering = {i+1: o for i,o in enumerate(ordering)}
        self.in_ordering = {self.out_ordering[a]: a for a in self.out_ordering}

        # If the Cartan type supplied is relabeled, compose these orderings
        # with the relabelling on the appropriate sides:
        if hasattr(self.cartan_type, '_relabelling'):
            r = self.cartan_type._relabelling
            r_inv = {v: k for (k, v) in r.items()}
            # Pre-compose in_ordering with r
            self.in_ordering = {i: self.in_ordering[r[i]] for i in self.in_ordering}
            # Post-compose out_ordering with r inverse
            self.out_ordering = {i: r_inv[self.out_ordering[i]] for i in self.out_ordering}

        # Check that the Coxeter matrices match up.
        cox_mat = CoxeterMatrix(self.coxeter_matrix(), self.cartan_type.index_set())
        if cox_mat != CoxeterMatrix(self.cartan_type):
            print("Warning, differing Coxeter matrices")

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
        from sage.arith.srange import srange
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
                return [Integer(1), Integer(2)]
        elif t in ['E']:
            if is_affine:
                return srange(1, r) + [ZZ.zero()]
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
        if type(other) is not type(self):
            if op in (Py_LT, Py_LE, Py_GT, Py_GE):
                return NotImplemented
            return op == Py_NE

        s_t = self.type()
        o_t = other.type()
        s_r = self.rank()
        o_r = other.rank()

        if op == Py_EQ:
            return s_t == o_t and s_r == o_r
        elif op == Py_NE:
            return s_t != o_t or s_r != o_r
        elif op == Py_LT:
            return s_t < o_t or (s_t == o_t and s_r < o_r)
        elif op == Py_LE:
            return s_t < o_t or (s_t == o_t and s_r <= o_r)
        elif op == Py_GT:
            return s_t > o_t or (s_t == o_t and s_r > o_r)
        elif op == Py_GE:
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
        del self.x

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
        cdef c_List_CoxWord l = c_List_CoxWord(0)
        interval(l, self.x[0], ww.word, vv.word)
        bruhat_interval = []
        cdef CoxGroupElement u
        cdef CoxGroupElement gg = CoxGroupElement(self, [])
        cdef size_t j
        for j in range(l.size()):
            u = gg._new()
            u.word = l[j]
            bruhat_interval.append(u)

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
        return Type(bytes_to_str(self.x.type().name().ptr()))

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
            from sage.rings.infinity import infinity
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

        Raises an error if ``self`` is not finite.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup       # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                # optional - coxeter3
            sage: W.full_context()                                                      # optional - coxeter3
            sage: W = CoxGroup(['A', 2,1])                                              # optional - coxeter3
            sage: W.full_context()                                                      # optional - coxeter3
            Traceback (most recent call last):
            ...
            TypeError: group needs to be finite
        """
        if not self.is_finite():
            raise TypeError("group needs to be finite")
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
            TypeError: group needs to be finite
        """
        self.full_context()
        cdef c_FiniteCoxGroup* fcoxgroup = <c_FiniteCoxGroup*>(self.x)
        cdef CoxGroupElement w0 = CoxGroupElement(self, [])
        w0.word = fcoxgroup.longest_coxword()
        return w0

    def __call__(self, w):
        """
        Return a reduced expression for ``w``.

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
        from sage.matrix.constructor import matrix
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

        OUTPUT: a Sage graph

        .. NOTE::

           This uses the labels native to Coxeter3. This is useful
           when trying to obtain the mapping between the labels of
           Sage and Coxeter3.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup# optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                         # optional - coxeter3
            sage: W.coxeter_graph()                                              # optional - coxeter3
            Graph on 5 vertices
            sage: W.coxeter_graph().edges(sort=True)                             # optional - coxeter3
            [(1, 2, None), (2, 3, None), (3, 4, None), (4, 5, None)]
        """
        from sage.graphs.graph import Graph
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
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: CoxGroupElement(W, [1,2,3,2,3])                                                   # optional - coxeter3
            [1, 3, 2]
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [1,2,3,2,3])                                               # optional - coxeter3
            sage: del w                                                                             # optional - coxeter3
        """
        self.group = (<CoxGroup>group).x
        self._parent_group = group
        self.word.reset()
        for i in w:
            self.word.append(self._parent_group.in_ordering[i])

        if normal_form:
            self.group.normalForm(self.word)

    def _coxnumber(self):
        """
        Return the internal integer used by Coxeter3 to represent this element.

        TESTS::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupElement  # optional - coxeter3
            sage: W = CoxGroup(['A', 4])                                                            # optional - coxeter3
            sage: w = CoxGroupElement(W, [1,2,3,2,3])                                               # optional - coxeter3
            sage: w._coxnumber()                                                                    # optional - coxeter3
            7
        """
        return int(self.group.extendContext(self.word))

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: TestSuite(w).run()                                                                # optional - coxeter3
        """
        return (CoxGroupElement, (self._parent_group, list(self)))

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
        return CoxGroupElement(self._parent_group, reversed(self))

    inverse = __invert__

    cpdef CoxGroup parent_group(self):
        """
        Return the parent Coxeter group for this element.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import *                                          # optional - coxeter3
            sage: W = CoxGroup(['A',5])                                                             # optional - coxeter3
            sage: w = W([1,2,3])                                                                    # optional - coxeter3
            sage: w.parent_group()                                                                        # optional - coxeter3
            Coxeter group of type A and rank 5

        """
        return self._parent_group

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
            return [self[ii] for ii in range(*i.indices(len(self)))]
        if i < 0:
            i += len(self)
        if i >= len(self):
            raise IndexError("The index (%d) is out of range." % i)

        return self._parent_group.out_ordering[self.word[i]]

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
        return hash((self.__class__.__name__, self.parent_group(), tuple(self)))

    def __richcmp__(CoxGroupElement self, other, int op):
        """
        EXAMPLES::

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
        if type(other) is not type(self):
            if op in (Py_LT, Py_LE, Py_GT, Py_GE):
                return NotImplemented
            return op == Py_NE

        s_p = self.parent_group()
        o_p = other.parent_group()
        s_l = list(self)
        o_l = list(other)

        if op == Py_EQ:
            return s_p == o_p and s_l == o_l
        elif op == Py_NE:
            return s_p != o_p or s_l != o_l
        elif op == Py_LT:
            return s_p < o_p or (s_p == o_p and s_l < o_l)
        elif op == Py_LE:
            return s_p < o_p or (s_p == o_p and s_l <= o_l)
        elif op == Py_GT:
            return s_p > o_p or (s_p == o_p and s_l > o_l)
        elif op == Py_GE:
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
        return (self[i] for i in range(len(self)))

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
        return LFlags_to_list(self._parent_group, self.group.ldescent(self.word))

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
        return LFlags_to_list(self._parent_group, self.group.rdescent(self.word))

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
        cdef CoxGroupElement ww = CoxGroupElement(self._parent_group, w)
        return self.group.inOrder(self.word, ww.word)

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
        cdef Generator ss = self._parent_group.in_ordering[s]
        return self.group.isDescent(self.word, s)

    cdef CoxGroupElement _new(self):
        """
        Return a new copy of this element.
        """
        cdef CoxGroupElement res = CoxGroupElement(self.parent_group(), [])
        res.word = self.word
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
        cdef c_List_CoxWord list = c_List_CoxWord(0)
        self.group.coatoms(list, self.word)

        coatoms = []

        cdef Length i = 0
        cdef CoxGroupElement res
        for i in range(list.size()):
            res = self._new()
            res.word = list[i]
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
        Return the Poincaré polynomial associated with the Bruhat
        interval between the identity element and this one.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup                    # optional - coxeter3
            sage: W = CoxGroup(['A', 5])                                                             # optional - coxeter3
            sage: W([]).poincare_polynomial()                                                        # optional - coxeter3
            1
            sage: W([1,2,1]).poincare_polynomial()                                                   # optional - coxeter3
            t^3 + 2*t^2 + 2*t + 1
        """
        cdef CoxGroup W = self.parent_group()
        cdef c_List_CoxWord result = c_List_CoxWord(0)
        cdef CoxGroupElement id = CoxGroupElement(W, [])
        cdef CoxGroupElement ww = CoxGroupElement(W, self)
        interval(result, W.x[0], id.word, ww.word)

        cdef list coefficients = [0]*(len(ww)+1)
        cdef size_t j
        for j in range(result.size()):
            coefficients[result[j].length()] += 1
        return ZZ['t'](coefficients)


    def kazhdan_lusztig_polynomial(self, v):
        """
        Return the Kazhdan-Lusztig polynomial `P_{u,v}` where `u` is ``self``.

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
        cdef CoxGroupElement vv
        if not isinstance(v, CoxGroupElement):
            vv = CoxGroupElement(self._parent_group, v)
        else:
            vv = v

        ZZq = PolynomialRing(ZZ, 'q')
        if not self.group.inOrder(self.word, vv.word):
            return ZZq.zero()

        cdef CoxNbr x = self.group.extendContext(self.word)
        cdef CoxNbr y = self.group.extendContext(vv.word)
        cdef c_KLPol kl_poly = self.group.klPol(x, y)
        if kl_poly.isZero():
            return ZZq.zero()
        cdef size_t i
        l = [kl_poly[i] for i in range(kl_poly.deg()+1)]
        return ZZq(l)

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
        cdef CoxGroupElement vv = CoxGroupElement(self._parent_group, v)
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

class CoxGroupIterator():
    def __init__(self, group):
        """
        A class used to iterate over all of the elements of a Coxeter group.

        .. NOTE::

           This will construct all of the elements of the group within
           Coxeter3.  For some groups, this may be too large to fit
           into memory.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                            # optional - coxeter3
            sage: it = CoxGroupIterator(W)                                                          # optional - coxeter3
            sage: [next(it) for i in range(W.order())]                                              # optional - coxeter3
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

    def __next__(self):
        """
        Return the next element in the associated Coxeter group.

        EXAMPLES::

            sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator # optional - coxeter3
            sage: W = CoxGroup(['A', 2])                                                            # optional - coxeter3
            sage: it = CoxGroupIterator(W)                                                          # optional - coxeter3
            sage: next(it)                                                                          # optional - coxeter3
            []
        """
        if self.n >= self.order:
            raise StopIteration
        cdef CoxGroupElement w = self.group([])

        (<CoxGroup>self.group).x.prod_nbr(w.word, self.n)
        self.n += 1
        return w

    next = __next__


CoxGroup_cache = {}
def get_CoxGroup(cartan_type):
    """
    TESTS::

        sage: from sage.libs.coxeter3.coxeter import get_CoxGroup as CoxGroup, CoxGroupIterator  # optional - coxeter3
        sage: W = CoxGroup(['A', 2])                                                             # optional - coxeter3
    """
    from sage.combinat.root_system.cartan_type import CartanType
    cartan_type = CartanType(cartan_type)
    if cartan_type not in CoxGroup_cache:
        CoxGroup_cache[cartan_type] = CoxGroup(cartan_type)
    return CoxGroup_cache[cartan_type]
