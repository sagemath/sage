"""
Permutation group elements

AUTHORS:

- David Joyner (2006-02)

- David Joyner (2006-03): word problem method and reorganization

- Robert Bradshaw (2007-11): convert to Cython

EXAMPLES: The Rubik's cube group::

    sage: f= [(17,19,24,22),(18,21,23,20),(6,25,43,16),(7,28,42,13),(8,30,41,11)]
    sage: b=[(33,35,40,38),(34,37,39,36),( 3, 9,46,32),( 2,12,47,29),( 1,14,48,27)]
    sage: l=[( 9,11,16,14),(10,13,15,12),( 1,17,41,40),( 4,20,44,37),( 6,22,46,35)]
    sage: r=[(25,27,32,30),(26,29,31,28),( 3,38,43,19),( 5,36,45,21),( 8,33,48,24)]
    sage: u=[( 1, 3, 8, 6),( 2, 5, 7, 4),( 9,33,25,17),(10,34,26,18),(11,35,27,19)]
    sage: d=[(41,43,48,46),(42,45,47,44),(14,22,30,38),(15,23,31,39),(16,24,32,40)]
    sage: cube = PermutationGroup([f,b,l,r,u,d])
    sage: F=cube.gens()[0]
    sage: B=cube.gens()[1]
    sage: L=cube.gens()[2]
    sage: R=cube.gens()[3]
    sage: U=cube.gens()[4]
    sage: D=cube.gens()[5]
    sage: cube.order()
    43252003274489856000
    sage: F.order()
    4

The interested user may wish to explore the following commands:
move = cube.random_element() and time word_problem([F,B,L,R,U,D],
move, False). This typically takes about 5 minutes (on a 2 Ghz
machine) and outputs a word ('solving' the cube in the position
move) with about 60 terms or so.

OTHER EXAMPLES: We create element of a permutation group of large
degree.

::

    sage: G = SymmetricGroup(30)
    sage: s = G(srange(30,0,-1)); s
    (1,30)(2,29)(3,28)(4,27)(5,26)(6,25)(7,24)(8,23)(9,22)(10,21)(11,20)(12,19)(13,18)(14,17)(15,16)
"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


import random

import sage.groups.old as group

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
from cpython.list cimport *

from sage.rings.all      import ZZ, Integer, is_MPolynomial, is_Polynomial
from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
import sage.structure.coerce as coerce
from sage.misc.superseded import deprecated_function_alias

import operator

from sage.rings.fast_arith cimport arith_llong
cdef arith_llong arith = arith_llong()
cdef extern from *:
    long long LONG_LONG_MAX

#import permgroup_named

def make_permgroup_element(G, x):
    """
    Returns a PermutationGroupElement given the permutation group
    ``G`` and the permutation ``x`` in list notation.

    This is function is used when unpickling old (pre-domain) versions
    of permutation groups and their elements.  This now does a bit of
    processing and calls :func:`make_permgroup_element_v2` which is
    used in unpickling the current PermutationGroupElements.

    EXAMPLES::

        sage: from sage.groups.perm_gps.permgroup_element import make_permgroup_element
        sage: S = SymmetricGroup(3)
        sage: make_permgroup_element(S, [1,3,2])
        (2,3)
    """
    domain = FiniteEnumeratedSet(range(1, len(x)+1))
    return make_permgroup_element_v2(G, x, domain)

def make_permgroup_element_v2(G, x, domain):
    """
    Returns a PermutationGroupElement given the permutation group
    ``G``, the permutation ``x`` in list notation, and the domain
    ``domain`` of the permutation group.

    This is function is used when unpickling permutation groups and
    their elements.

    EXAMPLES::

        sage: from sage.groups.perm_gps.permgroup_element import make_permgroup_element_v2
        sage: S = SymmetricGroup(3)
        sage: make_permgroup_element_v2(S, [1,3,2], S.domain())
        (2,3)
    """
    # Note that it has to be in-sync with the __init__ method of
    # PermutationGroup_generic since the elements have to be created
    # before the PermutationGroup_generic is initialized.  The
    # constructor for PermutationGroupElement requires that
    # G._domain_to_gap be set.
    G._domain = domain
    G._deg = len(domain)
    G._domain_to_gap = dict([(key, i+1) for i, key in enumerate(domain)])
    G._domain_from_gap = dict([(i+1, key) for i, key in enumerate(domain)])
    return G(x, check=False)


def is_PermutationGroupElement(x):
    """
    Returns True if ``x`` is a PermutationGroupElement.

    EXAMPLES::

        sage: p = PermutationGroupElement([(1,2),(3,4,5)])
        sage: from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
        sage: is_PermutationGroupElement(p)
        True
    """
    return isinstance(x, PermutationGroupElement)

def string_to_tuples(g):
    """
    EXAMPLES::

        sage: from sage.groups.perm_gps.permgroup_element import string_to_tuples
        sage: string_to_tuples('(1,2,3)')
        [(1, 2, 3)]
        sage: string_to_tuples('(1,2,3)(4,5)')
        [(1, 2, 3), (4, 5)]
        sage: string_to_tuples(' (1,2, 3) (4,5)')
        [(1, 2, 3), (4, 5)]
        sage: string_to_tuples('(1,2)(3)')
        [(1, 2), (3,)]
    """
    from sage.misc.all import sage_eval

    if not isinstance(g, str):
        raise ValueError, "g (= %s) must be a string"%g
    elif g == '()':
        return []
    g = g.replace('\n','').replace(' ', '').replace(')(', '),(').replace(')', ',)')
    g = '[' + g + ']'
    return sage_eval(g, preparse=False)

def standardize_generator(g, convert_dict=None):
    """
    Standardizes the input for permutation group elements to a list of
    tuples.  This was factored out of the
    PermutationGroupElement.__init__ since
    PermutationGroup_generic.__init__ needs to do the same computation
    in order to compute the domain of a group when it's not explicitly
    specified.

    INPUT:

    - ``g`` - a list, tuple, string, GapElement,
      PermutationGroupElement, Permutation

    - ``convert_dict`` - (optional) a dictionary used to convert the
      points to a number compatible with GAP.

    OUTPUT:

    The permutation in as a list of cycles.

    EXAMPLES::

        sage: from sage.groups.perm_gps.permgroup_element import standardize_generator
        sage: standardize_generator('(1,2)')
        [(1, 2)]

        sage: p = PermutationGroupElement([(1,2)])
        sage: standardize_generator(p)
        [(1, 2)]
        sage: standardize_generator(p._gap_())
        [(1, 2)]
        sage: standardize_generator((1,2))
        [(1, 2)]
        sage: standardize_generator([(1,2)])
        [(1, 2)]
        sage: standardize_generator(Permutation([2,1,3]))
        [(1, 2), (3,)]

    ::

        sage: d = {'a': 1, 'b': 2}
        sage: p = SymmetricGroup(['a', 'b']).gen(0); p
        ('a','b')
        sage: standardize_generator(p, convert_dict=d)
        [(1, 2)]
        sage: standardize_generator(p._gap_(), convert_dict=d)
        [(1, 2)]
        sage: standardize_generator(('a','b'), convert_dict=d)
        [(1, 2)]
        sage: standardize_generator([('a','b')], convert_dict=d)
        [(1, 2)]
    """
    from sage.interfaces.gap import GapElement
    from sage.combinat.permutation import Permutation
    from sage.libs.pari.all import pari_gen

    if isinstance(g, pari_gen):
        g = list(g)

    needs_conversion = True
    if isinstance(g, GapElement):
        g = str(g)
        needs_conversion = False
    if isinstance(g, Permutation):
        return g.cycle_tuples()
    if isinstance(g, PermutationGroupElement):
        g = g.cycle_tuples()
    if isinstance(g, str):
        g = string_to_tuples(g)
    if isinstance(g, tuple) and (len(g) == 0 or not isinstance(g[0], tuple)):
        g = [g]

    #Get the permutation in list notation
    if PyList_CheckExact(g) and (len(g) == 0 or not PY_TYPE_CHECK(g[0], tuple)):
        if convert_dict is not None and needs_conversion:
            g = [convert_dict[x] for x in g]


        #FIXME: Currently, we need to verify that g defines an actual
        #permutation since the constructor Permutation does not
        #perform this check.  When it does, we should remove this code.
        #See #8392
        if set(g) != set(range(1, len(g)+1)):
            raise ValueError, "Invalid permutation vector: %s"%g
        return Permutation(g).cycle_tuples()
    else:
        if convert_dict is not None and needs_conversion:
            g = [tuple([convert_dict[x] for x in cycle])for cycle in g]

    return g

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    """
    An element of a permutation group.

    EXAMPLES::

        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5)]
        sage: g = G.random_element()
        sage: g in G
        True
        sage: g = G.gen(0); g
        (1,2,3)(4,5)
        sage: print g
        (1,2,3)(4,5)
        sage: g*g
        (1,3,2)
        sage: g**(-1)
        (1,3,2)(4,5)
        sage: g**2
        (1,3,2)
        sage: G = PermutationGroup([(1,2,3)])
        sage: g = G.gen(0); g
        (1,2,3)
        sage: g.order()
        3

    This example illustrates how permutations act on multivariate
    polynomials.

    ::

        sage: R = PolynomialRing(RationalField(), 5, ["x","y","z","u","v"])
        sage: x, y, z, u, v = R.gens()
        sage: f = x**2 - y**2 + 3*z**2
        sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
        sage: sigma = G.gen(0)
        sage: f * sigma
        3*x^2 + y^2 - z^2
    """
    def __init__(self, g, parent = None, check = True):
        r"""
        Create element of a permutation group.

        There are several ways to define a permutation group element:


        -  Define a permutation group `G`, then use
           ``G.gens()`` and multiplication \* to construct
           elements.

        -  Define a permutation group `G`, then use e.g.,
           ``G([(1,2),(3,4,5)])`` to construct an element of the
           group. You could also use ``G('(1,2)(3,4,5)')``

        -  Use e.g.,
           ``PermutationGroupElement([(1,2),(3,4,5)])`` or
           ``PermutationGroupElement('(1,2)(3,4,5)')`` to make a
           permutation group element with parent `S_5`.


        INPUT:


        -  ``g`` - defines element

        -  ``parent (optional)`` - defines parent group (g must
           be in parent if specified, or a TypeError is raised).

        -  ``check`` - bool (default: True), if False assumes g
           is a gap element in parent (if specified).


        EXAMPLES: We illustrate construction of permutation using several
        different methods.

        First we construct elements by multiplying together generators for
        a group.

        ::

            sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'], canonicalize=False)
            sage: s = G.gens()
            sage: s[0]
            (1,2)(3,4)
            sage: s[1]
            (3,4,5,6)
            sage: s[0]*s[1]
            (1,2)(3,5,6)
            sage: (s[0]*s[1]).parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]

        Next we illustrate creation of a permutation using coercion into an
        already-created group.

        ::

            sage: g = G([(1,2),(3,5,6)])
            sage: g
            (1,2)(3,5,6)
            sage: g.parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]
            sage: g == s[0]*s[1]
            True

        We can also use a string instead of a list to specify the
        permutation.

        ::

            sage: h = G('(1,2)(3,5,6)')
            sage: g == h
            True

        We can also make a permutation group element directly using the
        ``PermutationGroupElement`` command. Note that the
        parent is then the full symmetric group `S_n`, where
        `n` is the largest integer that is moved by the
        permutation.

        ::

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k
            (1,2)(3,5,6)
            sage: k.parent()
            Symmetric group of order 6! as a permutation group

        Note the comparison of permutations doesn't require that the parent
        groups are the same.

        ::

            sage: k == g
            True

        Arithmetic with permutations having different parents is also
        defined::

            sage: k*g
            (3,6,5)
            sage: (k*g).parent()
            Symmetric group of order 6! as a permutation group

        ::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: loads(dumps(G.0)) == G.0
            True

        EXAMPLES::

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k._gap_()
            (1,2)(3,5,6)
            sage: k._gap_().parent()
            Gap

        List notation::

            sage: PermutationGroupElement([1,2,4,3,5])
            (3,4)

        TESTS::

            sage: PermutationGroupElement(())
            ()
            sage: PermutationGroupElement([()])
            ()
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.groups.perm_gps.permgroup import PermutationGroup_generic
        from sage.combinat.permutation import from_cycles

        convert_dict = parent._domain_to_gap if parent is not None else None
        try:
            v = standardize_generator(g, convert_dict)
        except KeyError:
            raise ValueError, "Invalid permutation vector: %s" % g


        degree = max([1] + [max(cycle+(1,)) for cycle in v])
        v = from_cycles(degree, v)

        self.__gap = 'PermList(%s)'%v

        if parent is None:
            parent = SymmetricGroup(len(v))

        if check and parent.__class__ != SymmetricGroup:
            if not (parent is None or isinstance(parent, PermutationGroup_generic)):
                raise TypeError, 'parent must be a permutation group'
            if parent is not None:
                P = parent._gap_()
                if not P.parent()(self.__gap) in P:
                    raise TypeError, 'permutation %s not in %s'%(g, parent)

        Element.__init__(self, parent)

        self.n = max(parent.degree(), 1)

        if self.perm is NULL or self.perm is self.perm_buf:
            self.perm = <int *>sage_malloc(sizeof(int) * self.n)
        else:
            self.perm = <int *>sage_realloc(self.perm, sizeof(int) * self.n)


        cdef int i, vn = len(v)
        assert(vn <= self.n)
        for i from 0 <= i < vn:
            self.perm[i] = v[i] - 1
        for i from vn <= i < self.n:
            self.perm[i] = i

        # We do this check even if check=False because it's fast
        # (relative to other things in this function) and the
        # rest of the code is assumes that self.perm specifies
        # a valid permutation (else segfaults, infinite loops may occur).
        if not is_valid_permutation(self.perm, vn):
            raise ValueError, "Invalid permutation vector: %s" % v

    def __cinit__(self, g = None, parent = None, check = True):
        self.perm = NULL

    def __dealloc__(self):
        if self.perm is not NULL and self.perm is not self.perm_buf:
            sage_free(self.perm)

    def __reduce__(self):
        """
        Returns a function and its arguments needed to create this
        permutation group element.  This is used in pickling.

        EXAMPLES::

           sage: g = PermutationGroupElement([(1,2,3),(4,5)]); g
           (1,2,3)(4,5)
           sage: func, args = g.__reduce__()
           sage: func(*args)
           (1,2,3)(4,5)
        """
        return make_permgroup_element_v2, (self._parent, self.domain(), self._parent.domain())

    cdef PermutationGroupElement _new_c(self):
        cdef PermutationGroupElement other = PY_NEW_SAME_TYPE(self)
        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__
        other._parent = self._parent
        other.n = self.n
        if other.n <= sizeof(other.perm_buf) / sizeof(int):
            other.perm = other.perm_buf
        else:
            other.perm = <int *>sage_malloc(sizeof(int) * other.n)
        return other

    def _gap_(self, gap=None):
        """
        Returns

        EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)]); g
            (1,2,3)(4,5)
            sage: a = g._gap_(); a
            (1,2,3)(4,5)
            sage: g._gap_() is g._gap_()
            True

        Note that only one GapElement is cached:

            sage: gap2 = Gap()
            sage: b = g._gap_(gap2)
            sage: c = g._gap_()
            sage: a is c
            False
        """
        if (self._gap_element is None or
            (gap is not None and self._gap_element._parent is not gap)):
            if gap is None:
                from sage.interfaces.gap import gap
            self._gap_element = gap(self._gap_init_())
        return self._gap_element

    def _gap_init_(self):
        """
        Returns a GAP string representation for this
        PermutationGroupElement.

        EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._gap_init_()
            'PermList([2, 3, 1, 5, 4])'
        """
        return 'PermList(%s)'%self._gap_list()

    def _repr_(self):
        """
        Return string representation of this permutation.

        EXAMPLES:

        We create the permutation `(1,2,3)(4,5)` and
        print it. ::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._repr_()
            '(1,2,3)(4,5)'

        Permutation group elements support renaming them so they print
        however you want, as illustrate below::

            sage: g.rename('sigma')
            sage: g
            sigma
            sage: g.rename()
            sage: g
            (1,2,3)(4,5)
        """
        return self.cycle_string()

    def _latex_(self):
        r"""
        Returns a latex representation of this permutation.

        EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: latex(g)
            (1,2,3)(4,5)

            sage: S = SymmetricGroup(['a', 'b'])
            sage: latex(S.gens())
            \left[(\text{\texttt{a}},\text{\texttt{b}})\right]
        """
        from sage.misc.latex import latex
        return "".join(["(" + ",".join([latex(x) for x in cycle])+")" for cycle in self.cycle_tuples()])

    def __getitem__(self, i):
        """
        Return the ith permutation cycle in the disjoint cycle
        representation of self.

        INPUT:


        -  ``i`` - integer


        OUTPUT: a permutation group element

        EXAMPLE::

            sage: G = PermutationGroup([[(1,2,3),(4,5)]],5)
            sage: g = G.gen(0)
            sage: g[0]
            (1,2,3)
            sage: g[1]
            (4,5)
        """
        return self.cycles()[i]

    def __cmp__(PermutationGroupElement self, PermutationGroupElement right):
        """
        Compare group elements self and right.

        EXAMPLES::

            sage: G = PermutationGroup([[(3,4)], [(1,2,3),(4,5)]])
            sage: G.gen(0) != G.gen(1)
            True
            sage: G.gen(0) == G.gen(1)
            False

        Permutations are ordered left lexicographically on their
        associated 'lists'; thus in the symmetric group S(5), the
        permutation (1,2)(3,4), which corresponds to the list
        [2,1,4,3,5], is larger than (1,2), which corresponds to the
        list [2,1,3,4,5].

        ::

            sage: S = SymmetricGroup(5)
            sage: S("(1,2)(3,4)") < S("(1,2)")
            False
            sage: S("(1,2)(3,4)") > S("(1,2)")
            True

        TESTS:

        Verify that we fixed bug #5537::

            sage: h = PermutationGroupElement('(1,3,2)')
            sage: k = PermutationGroupElement('(1,2,3),(4,5)')
            sage: k^2 == h, h == k^2
            (True, True)
            sage: k^6 == PermutationGroupElement('()')
            True
        """
        cdef int i
        cdef int swap = 1
        if right.n < self.n:
            self, right = right, self
            swap = -1
        for i from 0 <= i < self.n:
            if self.perm[i] < right.perm[i]:
                return -swap
            elif self.perm[i] > right.perm[i]:
                return swap
        for i from self.n <= i < right.n:
            if i < right.perm[i]:
                return -swap
            elif i > right.perm[i]:
                return swap
        return 0

    def __call__(self, i):
        """
        Returns the image of the integer i under this permutation.
        Alternately, if i is a list, tuple or string, returns the result of
        self acting on i.

        EXAMPLE::

            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G
            Permutation Group with generators [(1,2,3)(4,5)]
            sage: g = G.gen(0)
            sage: g(5)
            4
            sage: g('abcde')
            'bcaed'
            sage: g([0,1,2,3,4])
            [1, 2, 0, 4, 3]
            sage: g(('who','what','when','where','why'))
            ('what', 'when', 'who', 'why', 'where')

        ::

            sage: g(x)
            Traceback (most recent call last):
            ...
            ValueError: Must be in the domain or a list, tuple or string.
            sage: g(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Must be in the domain or a list, tuple or string.
        """
        to_gap = self._parent._domain_to_gap
        from_gap = self._parent._domain_from_gap
        cdef int j

        try:
            i = to_gap[i]
        except (KeyError, TypeError):
            # We currently have to include this to maintain the
            # current behavior where if you pass in an integer which
            # is not in the domain of the permutation group, then that
            # integer itself will be returned.
            if isinstance(i, (long, int, Integer)):
                return i


            if not isinstance(i,(list,tuple,str)):
                raise ValueError, "Must be in the domain or a list, tuple or string."

            permuted = [i[self.perm[j]] for j from 0 <= j < self.n]
            if PY_TYPE_CHECK(i, tuple):
                permuted = tuple(permuted)
            elif PY_TYPE_CHECK(i, str):
                permuted = ''.join(permuted)
            permuted += i[self.n:]
            return permuted
        else:
            j = i
            if 1 <= j <= self.n:
                return from_gap[self.perm[j-1]+1]
            else:
                return from_gap[i]

    cpdef list _act_on_list_on_position(self, list x):
        """
        Returns the right action of ``self`` on the list ``x``. This is the
        action on positions.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3,4,5,6)]])
            sage: p = G.gen()^2; p
            (1,3,5)(2,4,6)
            sage: p._act_on_list_on_position([1,2,3,4,5,6])
            [3, 4, 5, 6, 1, 2]
            sage: p._act_on_list_on_position(['a','b','c','d','e','f'])
            ['c', 'd', 'e', 'f', 'a', 'b']
            sage: p._act_on_list_on_position(['c','d','e','f','a','b'])
            ['e', 'f', 'a', 'b', 'c', 'd']
            sage: p._act_on_list_on_position([])
            Traceback (most recent call last):
            ...
            AssertionError: (1,3,5)(2,4,6) and [] should have the same length
            sage: p._act_on_list_on_position([1,2,3,4,5,6,7])
            Traceback (most recent call last):
            ...
            AssertionError: (1,3,5)(2,4,6) and [1, 2, 3, 4, 5, 6, 7] should have the same length
        """
        assert len(x) == self.n, '%s and %s should have the same length'%(self, x)
        return [ x[self.perm[i]] for i in range(self.n) ]

    cpdef ClonableIntArray _act_on_array_on_position(self, ClonableIntArray x):
        """
        Returns the right action of ``self`` on the ClonableIntArray
        ``x``. This is the action on positions.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: v = IncreasingIntArrays()([1,2,3,4])
            sage: G = PermutationGroup([[(1,2,3,4)]])
            sage: id = G.identity()
            sage: id._act_on_array_on_position(v)
            [1, 2, 3, 4]
        """
        cdef int i
        cdef ClonableIntArray y
        cdef int l = self.n
        assert x._len == l, '%s and %s should have the same length'%(self, x)
        y = x.clone()
        for i in range(l):
            y._list[i] = x._list[self.perm[i]]
        y.set_immutable()
        return y

    cpdef _act_on_(self, x, bint self_on_left):
        """
        Return the right action of self on left.

        For example, if f=left is a polynomial, then this function returns
        f(sigma\*x), which is image of f under the right action of sigma on
        the indeterminates. This is a right action since the image of
        f(sigma\*x) under tau is f(sigma\*tau\*x).

        INPUT:


        -  ``left`` - element of space on which permutations
           act from the right


        EXAMPLES::

            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: R.<x,y,z,u,v> = PolynomialRing(QQ,5)
            sage: f = x^2 + y^2 - z^2 + 2*u^2
            sage: sigma, tau = G.gens()
            sage: f*sigma
            -x^2 + y^2 + z^2 + 2*v^2
            sage: f*tau
            y^2 + z^2 - u^2 + 2*v^2
            sage: f*(sigma*tau)
            2*x^2 - y^2 + z^2 + u^2
            sage: (f*sigma)*tau
            2*x^2 - y^2 + z^2 + u^2
        """
        if not self_on_left:
            left = x
            if is_Polynomial(left):
                if self != 1:
                    raise ValueError, "%s does not act on %s"%(self, left.parent())
                return left
            elif is_MPolynomial(left):
                R = left.parent()
                vars = R.gens()
                try:
                    sigma_x  = [vars[self(i+1)-1] for i in range(R.ngens())]
                except IndexError:
                    raise TypeError, "%s does not act on %s"%(self, left.parent())
                return left(tuple(sigma_x))

    cpdef MonoidElement _mul_(left, MonoidElement _right):
        """
        EXAMPLES::

            sage: S = SymmetricGroup(['a', 'b'])
            sage: s = S([('a', 'b')]); s
            ('a','b')
            sage: s*s
            ()
        """
        cdef PermutationGroupElement prod = left._new_c()
        cdef PermutationGroupElement right = <PermutationGroupElement>_right
        cdef int i
        for i from 0 <= i < left.n:
            prod.perm[i] = right.perm[left.perm[i]]
        return prod

    def __invert__(self):
        """
        Return the inverse of this permutation.

        EXAMPLES::

            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: ~g
            (1,3,2)(4,5)
            sage: (~g) * g
            ()
        """
        cdef PermutationGroupElement inv = self._new_c()
        cdef int i
        for i from 0 <= i < self.n:
            inv.perm[self.perm[i]] = i
        return inv

    cpdef _gap_list(self):
        """
        Returns this permutation in list notation compatible with the
        GAP numbering.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: s = S.gen(0); s
            (1,2,3)
            sage: s._gap_list()
            [2, 3, 1]

        ::

            sage: S = SymmetricGroup(['a', 'b', 'c'])
            sage: s = S.gen(0); s
            ('a','b','c')
            sage: s._gap_list()
            [2, 3, 1]
        """
        cdef int i
        return [self.perm[i]+1 for i from 0 <= i < self.n]

    def _gap_cycle_string(self):
        """
        Returns a cycle string for this permutation compatible with
        the GAP numbering.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: s = S.gen(0); s
            (1,2,3)
            sage: s._gap_cycle_string()
            '(1,2,3)'

        ::

            sage: S = SymmetricGroup(['a', 'b', 'c'])
            sage: s = S.gen(0); s
            ('a','b','c')
            sage: s._gap_cycle_string()
            '(1,2,3)'
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self._gap_list()).cycle_string()

    list = deprecated_function_alias(14319, domain)

    cpdef domain(self):
        """
        Returns the domain of self.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: x = G([2,1,4,3]); x
            (1,2)(3,4)
            sage: v = x.domain(); v
            [2, 1, 4, 3]
            sage: type(v[0])
            <type 'int'>
            sage: x = G([2,1]); x
            (1,2)
            sage: x.domain()
            [2, 1, 3, 4]

        TESTS::

            sage: S = SymmetricGroup(0)
            sage: x = S.one(); x
            ()
            sage: x.domain()
            []
        """
        cdef int i

        #We need to do this to handle the case of SymmetricGroup(0)
        #where the domain is (), but the permutation group element has
        #an underlying representation of [1].  The 1 doesn't
        #correspond to anything in the domain
        if len(self._parent._domain) == 0:
            return []
        else:
            from_gap = self._parent._domain_from_gap
            return [from_gap[self.perm[i]+1] for i from 0 <= i < self.n]

    def __hash__(self):
        """
        Return hash of this permutation.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: s = G([2,1,5,3,4])
            sage: s.tuple()
            (2, 1, 5, 3, 4)
            sage: hash(s)
            1592966088          # 32-bit
            2865702456085625800 # 64-bit
            sage: hash(s.tuple())
            1592966088          # 32-bit
            2865702456085625800 # 64-bit
        """
        return hash(tuple(self._gap_list()))

    def tuple(self):
        """
        Return tuple of images of the domain under self.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: s = G([2,1,5,3,4])
            sage: s.tuple()
            (2, 1, 5, 3, 4)

            sage: S = SymmetricGroup(['a', 'b'])
            sage: S.gen().tuple()
            ('b', 'a')
        """
        if self.__tuple is None:
            self.__tuple = tuple(self.domain())
        return self.__tuple

    def dict(self):
        """
        Returns a dictionary associating each element of the domain with its
        image.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: g = G((1,2,3,4)); g
            (1,2,3,4)
            sage: v = g.dict(); v
            {1: 2, 2: 3, 3: 4, 4: 1}
            sage: type(v[1])
            <type 'int'>
            sage: x = G([2,1]); x
            (1,2)
            sage: x.dict()
            {1: 2, 2: 1, 3: 3, 4: 4}
        """
        from_gap = self._parent._domain_from_gap
        to_gap = self._parent._domain_to_gap
        cdef int i
        return {e:from_gap[self.perm[i-1]+1] for e,i in to_gap.iteritems()}

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer `n` for which `g^n = 1`.

        EXAMPLES::

            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.order()
            6

        TESTS::

            sage: prod(primes(150))
            1492182350939279320058875736615841068547583863326864530410
            sage: L = [tuple(range(sum(primes(p))+1, sum(primes(p))+1+p)) for p in primes(150)]
            sage: PermutationGroupElement(L).order()
            1492182350939279320058875736615841068547583863326864530410
        """
        order = None
        cdef long long order_c = 1
        cdef int cycle_len
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            k = self.perm[i]
            cycle_len = 1
            while k != i:
                seen[k] = 1
                k = self.perm[k]
                cycle_len += 1
            if order is not None:
                order = order.lcm(cycle_len)
            else:
                order_c = (order_c * cycle_len) / arith.c_gcd_longlong(order_c, cycle_len)
                if order_c > LONG_LONG_MAX / (self.n - i):
                    order = Integer(order_c)
        sage_free(seen)
        return int(order_c) if order is None else order

    def inverse(self):
        r"""
        Returns the inverse permutation.

        OUTPUT:

        For an element of a permutation group, this method returns the inverse
        element, which is both the inverse function and the inverse as an
        element of a group.

        EXAMPLES::

            sage: s = PermutationGroupElement("(1,2,3)(4,5)")
            sage: s.inverse()
            (1,3,2)(4,5)

            sage: A = AlternatingGroup(4)
            sage: t = A("(1,2,3)")
            sage: t.inverse()
            (1,3,2)

        There are several ways (syntactically) to get an inverse
        of a permutation group element.  ::

            sage: s = PermutationGroupElement("(1,2,3,4)(6,7,8)")
            sage: s.inverse() == s^-1
            True
            sage: s.inverse() == ~s
            True
        """
        return ~self

    def sign(self):
        """
        Returns the sign of self, which is `(-1)^{s}`, where
        `s` is the number of swaps.

        EXAMPLES::

            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.sign()
            -1

        ALGORITHM: Only even cycles contribute to the sign, thus

        .. math::

            sign(sigma) = (-1)^{\sum_c len(c)-1}


        where the sum is over cycles in self.
        """
        cdef int cycle_len_sum = 0
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            k = self.perm[i]
            while k != i:
                seen[k] = 1
                k = self.perm[k]
                cycle_len_sum += 1
        sage_free(seen)
        return 1 - 2*(cycle_len_sum % 2) # == (-1)^cycle_len


    def orbit(self, n, bint sorted=True):
        """
        Returns the orbit of the integer `n` under this group
        element, as a sorted list.

        EXAMPLES::

            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.orbit(4)
            [4, 5]
            sage: g.orbit(3)
            [1, 2, 3]
            sage: g.orbit(10)
            [10]

        ::

            sage: s = SymmetricGroup(['a', 'b']).gen(0); s
            ('a','b')
            sage: s.orbit('a')
            ['a', 'b']
        """
        to_gap = self._parent._domain_to_gap
        from_gap = self._parent._domain_from_gap
        try:
            n = to_gap[n]
        except KeyError:
            return [n]

        cdef int i = n
        cdef int start = i
        if 1 <= i <= self.n:
            L = [from_gap[i]]
            i = self.perm[i-1]+1
            while i != start:
                PyList_Append(L,from_gap[i])
                i = self.perm[i-1]+1
            if sorted:
                L.sort()
            return L
        else:
            return from_gap[n]

    def cycles(self):
        """
        Return self as a list of disjoint cycles.

        EXAMPLES::

            sage: G = PermutationGroup(['(1,2,3)(4,5,6,7)'])
            sage: g = G.0
            sage: g.cycles()
            [(1,2,3), (4,5,6,7)]
            sage: a, b = g.cycles()
            sage: a(1), b(1)
            (2, 1)
        """
        L = []
        cdef PermutationGroupElement cycle
        cdef int i, j, k, next_k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            cycle = self._new_c()
            for j from 0 <= j < self.n: cycle.perm[j] = j
            k = cycle.perm[i] = self.perm[i]
            while k != i:
                seen[k] = 1
                next_k = cycle.perm[k] = self.perm[k]
                k = next_k
            PyList_Append(L, cycle)
        sage_free(seen)
        return L

    def cycle_tuples(self, singletons=False):
        """
        Return self as a list of disjoint cycles, represented as tuples
        rather than permutation group elements.

        INPUT:

        - ``singletons`` - boolean (default: False) whether or not consider the
          cycle that correspond to fixed point

        EXAMPLES::

            sage: p = PermutationGroupElement('(2,6)(4,5,1)')
            sage: p.cycle_tuples()
            [(1, 4, 5), (2, 6)]
            sage: p.cycle_tuples(singletons=True)
            [(1, 4, 5), (2, 6), (3,)]

        EXAMPLES::

            sage: S = SymmetricGroup(4)
            sage: S.gen(0).cycle_tuples()
            [(1, 2, 3, 4)]

        ::

            sage: S = SymmetricGroup(['a','b','c','d'])
            sage: S.gen(0).cycle_tuples()
            [('a', 'b', 'c', 'd')]
            sage: S([('a', 'b'), ('c', 'd')]).cycle_tuples()
            [('a', 'b'), ('c', 'd')]
        """
        from_gap = self._parent._domain_from_gap
        L = []
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i]:
                continue
            if self.perm[i] == i:
                if singletons:
                    PyList_Append(L, (from_gap[i+1],))
                    # it is not necessary to put seen[i] to 1 as we will never
                    # see i again
                else:
                    continue
            else:
                cycle = [from_gap[i+1]]
                k = self.perm[i]
                while k != i:
                    PyList_Append(cycle, from_gap[k+1])
                    seen[k] = 1
                    k = self.perm[k]
                PyList_Append(L, tuple(cycle))
        sage_free(seen)
        return L

    def cycle_string(self, singletons=False):
        """
        Return string representation of this permutation.

       EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g.cycle_string()
            '(1,2,3)(4,5)'

            sage: g = PermutationGroupElement([3,2,1])
            sage: g.cycle_string(singletons=True)
            '(1,3)(2)'
        """
        cycles = self.cycle_tuples(singletons)
        if len(cycles) == 0:
            return '()'
        return ''.join([repr(c) for c in cycles]).replace(', ',',').replace(',)',')')

    def has_descent(self, i, side = "right", positive = False):
        """
        INPUT:

         - ``i``: an element of the index set
         - ``side``: "left" or "right" (default: "right")
         - ``positive``: a boolean (default: False)

        Returns whether ``self`` has a left (resp. right) descent at
        position ``i``. If ``positive`` is True, then test for a non
        descent instead.

        Beware that, since permutations are acting on the right, the
        meaning of descents is the reverse of the usual
        convention. Hence, ``self`` has a left descent at position
        ``i`` if ``self(i) > self(i+1)``.

        EXAMPLES::

            sage: S = SymmetricGroup([1,2,3])
            sage: S.one().has_descent(1)
            False
            sage: S.one().has_descent(2)
            False
            sage: s = S.simple_reflections()
            sage: x = s[1]*s[2]
            sage: x.has_descent(1, side = "right")
            False
            sage: x.has_descent(2, side = "right")
            True
            sage: x.has_descent(1, side = "left")
            True
            sage: x.has_descent(2, side = "left")
            False
            sage: S._test_has_descent()

        The symmetric group acting on a set not of the form
        `(1,\dots,n)` is also supported::

            sage: S = SymmetricGroup([2,4,1])
            sage: s = S.simple_reflections()
            sage: x = s[2]*s[4]
            sage: x.has_descent(4)
            True
            sage: S._test_has_descent()
        """
        to_gap = self._parent._domain_to_gap
        from_gap = self._parent._domain_from_gap
        if side == "right":
            self = ~self

        try:
            i1 = from_gap[to_gap[i]+1]
        except KeyError:
            return False

        return (to_gap[self(i)] > to_gap[self(i1)]) is not positive

    def matrix(self):
        """
        Returns deg x deg permutation matrix associated to the permutation
        self

        EXAMPLES::

            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.matrix()
            [0 1 0 0 0]
            [0 0 1 0 0]
            [1 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 1 0]
        """
        M = MatrixSpace(ZZ, self.n, self.n, sparse=True)
        cdef int i
        entries = {}
        for i from 0 <= i < self.n:
            entries[i, self.perm[i]] = 1
        return M(entries)

    def word_problem(self, words, display=True):
        """
        G and H are permutation groups, g in G, H is a subgroup of G
        generated by a list (words) of elements of G. If g is in H, return
        the expression for g as a word in the elements of (words).

        This function does not solve the word problem in Sage. Rather it
        pushes it over to GAP, which has optimized algorithms for the word
        problem. Essentially, this function is a wrapper for the GAP
        functions "EpimorphismFromFreeGroup" and
        "PreImagesRepresentative".

        EXAMPLE::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]], canonicalize=False)
            sage: g1, g2 = G.gens()
            sage: h = g1^2*g2*g1
            sage: h.word_problem([g1,g2], False)
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')
            sage: h.word_problem([g1,g2])
               x1^2*x2^-1*x1
               [['(1,2,3)(4,5)', 2], ['(3,4)', -1], ['(1,2,3)(4,5)', 1]]
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')
        """
        if not self._parent._has_natural_domain():
            raise NotImplementedError

        import copy
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap

        G = gap(words[0].parent())
        g = words[0].parent()(self)
        H = gap.Group(words)
        ans = G.EpimorphismFromFreeGroup().PreImagesRepresentative(g)

        l1 = str(ans)
        l2 = copy.copy(l1)
        l4 = []
        l3 = l1.split("*")
        for i in range(1,len(words)+1):
            l2 = l2.replace("x"+str(i),str(words[i-1]))

        if display:
            for i in range(len(l3)):    ## parsing the word for display
                if len(l3[i].split("^"))==2:
                    l4.append([l3[i].split("^")[0],int(l3[i].split("^")[1])])
                if len(l3[i].split("^"))==1:
                    l4.append([l3[i].split("^")[0],1])
            l5 = copy.copy(l4)
            for i in range(len(l4)):
                for j in range(1,len(words)+1):
                    l5[i][0] = l5[i][0].replace("x"+str(j),str(words[j-1]))
            print "         ",l1
            print "         ",l5
        return l1,l2

    def conjugacy_class(self):
        r"""
        Return the conjugacy class of ``self``.

        EXAMPLES::

            sage: D = DihedralGroup(5)
            sage: g = D((1,3,5,2,4))
            sage: g.conjugacy_class()
            Conjugacy class of (1,3,5,2,4) in Dihedral group of order 10 as a permutation group
        """
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        return ConjugacyClassGAP(self.parent(), self)


cdef bint is_valid_permutation(int* perm, int n):
    """
    This is used in the __init__ method.

    Returns True iff the first n elements of perm are literally a
    permutation of [0, ..., n-1].

    TESTS::

        sage: S = SymmetricGroup(10)
        sage: PermutationGroupElement([2,1],S,check=False)
        (1,2)
        sage: PermutationGroupElement([1,1],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: Invalid permutation vector: [1, 1]
        sage: PermutationGroupElement([1,-1],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: Invalid permutation vector: [1, -1]
        sage: PermutationGroupElement([1,2,3,10],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: Invalid permutation vector: [1, 2, 3, 10]
    """
    cdef int i, ix
    # make everything is in bounds
    for i from 0 <= i < n:
        if not 0 <= perm[i] < n:
            return False
    # mark hit points by sign
    for i from 0 <= i < n:
        ix = -1-perm[i] if perm[i] < 0 else perm[i]
        perm[ix] = -1-perm[ix]
    # make sure everything is hit once, and reset signs
    for i from 0 <= i < n:
        if perm[i] >= 0:
            return False
        perm[i] = -1-perm[i]

    return True

