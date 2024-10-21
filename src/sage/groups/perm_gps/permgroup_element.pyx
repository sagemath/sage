# sage_setup: distribution = sagemath-gap
r"""
Permutation group elements

AUTHORS:

- David Joyner (2006-02)

- David Joyner (2006-03): word problem method and reorganization

- Robert Bradshaw (2007-11): convert to Cython

- Sebastian Oehms (2018-11): Added :meth:`gap` as synonym to
  :meth:`_gap_` (compatibility to libgap framework, see :issue:`26750`)

- Sebastian Oehms (2019-02): Implemented :meth:`gap` properly (:issue:`27234`)

There are several ways to define a permutation group element:

-  Define a permutation group `G`, then use ``G.gens()``
   and multiplication ``*`` to construct elements.

-  Define a permutation group `G`, then use, e.g.,
   ``G([(1,2),(3,4,5)])`` to construct an element of the
   group. You could also use ``G('(1,2)(3,4,5)')``

-  Use, e.g.,
   ``PermutationGroupElement([(1,2),(3,4,5)])`` or
   ``PermutationGroupElement('(1,2)(3,4,5)')`` to make a
   permutation group element with parent `S_5`.

EXAMPLES:

We illustrate construction of permutation using several
different methods.

First we construct elements by multiplying together generators for
a group::

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
already-created group::

    sage: g = G([(1,2),(3,5,6)])
    sage: g
    (1,2)(3,5,6)
    sage: g.parent()
    Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]
    sage: g == s[0]*s[1]
    True

We can also use a string or one-line notation to specify the
permutation::

    sage: h = G('(1,2)(3,5,6)')
    sage: i = G([2,1,5,4,6,3])
    sage: g == h == i
    True

The Rubik's cube group::

    sage: f = [(17,19,24,22),(18,21,23,20),( 6,25,43,16),( 7,28,42,13),( 8,30,41,11)]
    sage: b = [(33,35,40,38),(34,37,39,36),( 3, 9,46,32),( 2,12,47,29),( 1,14,48,27)]
    sage: l = [( 9,11,16,14),(10,13,15,12),( 1,17,41,40),( 4,20,44,37),( 6,22,46,35)]
    sage: r = [(25,27,32,30),(26,29,31,28),( 3,38,43,19),( 5,36,45,21),( 8,33,48,24)]
    sage: u = [( 1, 3, 8, 6),( 2, 5, 7, 4),( 9,33,25,17),(10,34,26,18),(11,35,27,19)]
    sage: d = [(41,43,48,46),(42,45,47,44),(14,22,30,38),(15,23,31,39),(16,24,32,40)]
    sage: cube = PermutationGroup([f, b, l, r, u, d])
    sage: F, B, L, R, U, D = cube.gens()
    sage: cube.order()
    43252003274489856000
    sage: F.order()
    4

We create element of a permutation group of large degree::

    sage: G = SymmetricGroup(30)
    sage: s = G(srange(30,0,-1)); s
    (1,30)(2,29)(3,28)(4,27)(5,26)(6,25)(7,24)(8,23)(9,22)(10,21)(11,20)(12,19)(13,18)(14,17)(15,16)

TESTS:

Check that :issue:`13569` is fixed::

    sage: [g*h for g in SymmetricGroup(3) for h in AlternatingGroup(3)]
    [(), (1,2,3), (1,3,2), (1,3,2), (), (1,2,3), (1,2,3), (1,3,2), (), (2,3),
    (1,2), (1,3), (1,3), (2,3), (1,2), (1,2), (1,3), (2,3)]
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Joyner
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import copy

from libc.stdlib cimport qsort

from cysignals.memory cimport sig_malloc, sig_calloc, sig_realloc, sig_free
from cpython.list cimport *

try:
    from cypari2.gen import Gen
except ImportError:
    Gen = ()

from sage.ext.stdsage cimport HAS_DICTIONARY
from sage.interfaces.abc import GpElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.structure.coerce cimport coercion_model
from sage.structure.element import Matrix
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool

import sage.interfaces.abc
import sage.groups.old as group
import sage.structure.coerce as coerce

from sage.ext.stdsage cimport HAS_DICTIONARY
from sage.interfaces.abc import GpElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.structure.coerce cimport coercion_model
from sage.structure.element import is_Matrix
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool

from sage.libs.gap.libgap import libgap
from sage.libs.gap.gap_includes cimport (UInt, UInt2, UInt4, T_PERM2, T_PERM4,
        NEW_PERM2, TNUM_OBJ, DEG_PERM2, DEG_PERM4, CONST_ADDR_PERM2,
        CONST_ADDR_PERM4, ADDR_PERM2)
from sage.libs.gap.util cimport initialize
from sage.libs.gap.element cimport (GapElement, GapElement_List,
        GapElement_String, GapElement_Permutation, make_GapElement_Permutation)
from sage.libs.gap.gap_includes cimport Obj, GAP_ValueInt, ELM_LIST

import operator

from sage.rings.fast_arith cimport arith_llong
cdef arith_llong arith = arith_llong()
cdef extern from *:
    long long LLONG_MAX

cdef int etuple_index_cmp(const void * a, const void * b) noexcept nogil:
    return ((<int *> a)[0] > (<int *> b)[0]) - ((<int *> a)[0] < (<int *> b)[0])


def make_permgroup_element(G, x):
    r"""
    Return a :class:`PermutationGroupElement` given the permutation group
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
    r"""
    Return a :class:`PermutationGroupElement` given the permutation group
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
    G._domain_to_gap = {key: i+1 for i, key in enumerate(domain)}
    G._domain_from_gap = {i+1: key for i, key in enumerate(domain)}
    return G.element_class(x, G, check=False)


def is_PermutationGroupElement(x):
    r"""
    Return ``True`` if ``x`` is a :class:`PermutationGroupElement`.

    EXAMPLES::

        sage: p = PermutationGroupElement([(1,2),(3,4,5)])
        sage: from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
        sage: is_PermutationGroupElement(p)
        doctest:warning...
        DeprecationWarning: The function is_PermutationGroupElement is deprecated;
        use 'isinstance(..., PermutationGroupElement)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38184,
                       "The function is_PermutationGroupElement is deprecated; "
                       "use 'isinstance(..., PermutationGroupElement)' instead.")
    return isinstance(x, PermutationGroupElement)


cdef class PermutationGroupElement(MultiplicativeGroupElement):
    r"""
    An element of a permutation group.

    EXAMPLES::

        sage: G = PermutationGroup(['(1,2,3)(4,5)']); G
        Permutation Group with generators [(1,2,3)(4,5)]
        sage: g = G.random_element()
        sage: g in G
        True
        sage: g = G.gen(0); g
        (1,2,3)(4,5)
        sage: print(g)
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
    def __dealloc__(self):
        if self.perm != NULL and self.perm != self.perm_buf:
            sig_free(self.perm)

    def __init__(self, g, parent, check=True):
        r"""
        Create element of a permutation group.

        There are several ways to define a permutation group element:

        -  Define a permutation group `G`, then use
           ``G.gens()`` and multiplication ``*`` to construct
           elements.

        -  Define a permutation group `G`, then use, e.g.,
           ``G([(1,2),(3,4,5)])`` to construct an element of the
           group. You could also use ``G('(1,2)(3,4,5)')``

        -  Use, e.g., ``PermutationGroupElement([(1,2),(3,4,5)])`` or
           ``PermutationGroupElement('(1,2)(3,4,5)')`` to make a
           permutation group element with parent `S_5`.

        INPUT:

        - ``g`` -- defines element

        - ``parent`` -- defines parent group (``g`` must be in
          parent if specified, or a :exc:`TypeError` is raised)

        - ``check`` -- boolean (default: ``True``); if ``False`` assumes ``g``
          is a gap element in parent (if specified)

        EXAMPLES:

        We illustrate construction of permutation using several
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

            sage: g = G([(1,2),(3,5,6)]); g
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
        :class:`PermutationGroupElement` command. Note that the
        parent is then the full symmetric group `S_n`, where
        `n` is the largest integer that is moved by the
        permutation.

        ::

            sage: k = PermutationGroupElement('(1,2)(3,5,6)'); k
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

        Initialization from gap or pari objects::

            sage: S = SymmetricGroup(5)
            sage: S(Permutation([5,1,4,3,2]))
            (1,5,2)(3,4)
            sage: S(gp.Vecsmall([5,1,4,3,2]))                                           # needs sage.libs.pari
            (1,5,2)(3,4)
            sage: S(gap.PermList([5,1,4,3,2]))
            (1,5,2)(3,4)
            sage: S(pari.Vecsmall([5,1,4,3,2]))                                         # needs sage.libs.pari
            (1,5,2)(3,4)
            sage: S(libgap.PermList([5,1,4,3,2]))
            (1,5,2)(3,4)
            sage: S(libgap([3,1,4,2,5]))
            (1,3,4,2)
            sage: S(libgap("(1,2,4)(3,5)"))
            (1,2,4)(3,5)

        Note that the conversion from gap permutations is agnostic to domains::

            sage: S = SymmetricGroup(['a', 'b', 'c'])
            sage: S(libgap.eval("(1,3)"))
            ('a','c')
            sage: S(gap("(1,3)"))
            ('a','c')

        Though not for gap list or string::

            sage: S(libgap([1, 3, 2]))
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: S(libgap(['a', 'c', 'b']))
            ('b','c')

            sage: S(libgap("(1,2)"))
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: S(libgap("('a','c')"))
            ('a','c')

        EXAMPLES::

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k._gap_()
            (1,2)(3,5,6)
            sage: k._gap_().parent()
            C library interface to GAP

        List notation::

            sage: PermutationGroupElement([1,2,4,3,5])
            (3,4)

        TESTS::

            sage: PermutationGroupElement(())
            ()
            sage: PermutationGroupElement([()])
            ()

        We check that :issue:`16678` is fixed::

            sage: Permutations.options.display='cycle'
            sage: p = Permutation((1,2))
            sage: PermutationGroupElement(p)
            (1,2)

        Bad input::

            sage: S5 = SymmetricGroup(5)
            sage: S5((3,-1,5))
            Traceback (most recent call last):
            ...
            ValueError: invalid list of cycles to initialize a permutation
            sage: S5((1,2,6,3))
            Traceback (most recent call last):
            ...
            ValueError: invalid list of cycles to initialize a permutation

            sage: P = PermutationGroup(['(1,2,3)'])
            sage: P('(1,2)')
            Traceback (most recent call last):
            ...
            ValueError: permutation (1,2) not in Permutation Group with generators [(1,2,3)]
        """
        cdef int i, degree = parent.degree()
        self._parent = parent
        self._alloc(degree)

        if not g:
            self._set_identity()
            return

        # TODO: this is for 1 base domain (ie {1, 2, ..., n}). We should
        # also introduce shortcut for 0 base domain.
        convert = not parent._has_natural_domain()
        if isinstance(g, tuple) and not isinstance(g[0], tuple):
            self._set_list_cycles([g], convert)
        elif isinstance(g, (list, tuple)):
            if isinstance(g[0], tuple):
                self._set_list_cycles(g, convert)
            else:
                self._set_list_images(g, convert)
        elif isinstance(g, str):
            self._set_string(g)
        elif isinstance(g, PermutationGroupElement):
            self._set_permutation_group_element(g, convert or not g.parent()._has_natural_domain())
        elif isinstance(g, Gen):
            self._set_list_images(g, convert)
        elif isinstance(g, GapElement):
            if isinstance(g, GapElement_Permutation):
                self._set_libgap(g)
            elif isinstance(g, GapElement_String):
                self._set_string(g.sage())
            elif isinstance(g, GapElement_List):
                self._set_list_images(g.sage(), convert)
            else:
                raise ValueError("invalid data to initialize a permutation")
        elif isinstance(g, sage.interfaces.abc.GapElement):
            if g.IsPerm():
                self._set_list_images(g.ListPerm(), False)
            else:
                raise ValueError("invalid data to initialize a permutation")
        elif isinstance(g, GpElement):
            self._set_list_images(list(g), convert)
            return
        else:
            # one more conversion that needs an extra import
            # (with circular dependencies)
            from sage.combinat.permutation import Permutation
            if isinstance(g, Permutation):
                self._set_list_images(g._list, convert)
            else:
                raise ValueError("invalid data to initialize a permutation")

        # We do this check even if check=False because it's fast
        # (relative to other things in this function) and the
        # rest of the code assumes that self.perm specifies
        # a valid permutation (else segfaults, infinite loops may occur).
        if not is_valid_permutation(self.perm, self.n):
            raise ValueError("invalid data to initialize a permutation")

        # This is more expensive
        if check and not isinstance(self, SymmetricGroupElement):
            P = parent._libgap_()
            p = libgap.PermList([self.perm[i] + 1 for i in range(self.n)])
            if p not in P:
                raise ValueError('permutation %s not in %s' % (g, parent))

    cpdef _set_identity(self):
        r"""
        TESTS::

            sage: p = PermutationGroupElement([3,1,4,2,5])
            sage: p._set_identity()
            sage: p
            ()

            sage: S = SymmetricGroup(["a", "b", "c"])
            sage: p = S(["b", "a", "c"])
            sage: p._set_identity()
            sage: p
            ()
        """
        cdef int i
        for i in range(self.n):
            self.perm[i] = i

    cpdef _set_list_images(self, v, bint convert):
        r"""
        TESTS::

            sage: p = PermutationGroupElement([3,1,4,2,5])
            sage: p._set_list_images([1,4,2,3,5], True)
            sage: p
            (2,4,3)
            sage: p._set_list_images([1,4,2,3,5], False)
            sage: p
            (2,4,3)

            sage: S = SymmetricGroup(["a", "b", "c"])
            sage: p = S(["b", "a", "c"])
            sage: p._set_list_images(["a", "c", "b"], True)
            sage: p
            ('b','c')
            sage: p._set_list_images([1, 3, 2], False)
            sage: p
            ('b','c')
        """
        cdef int i, j, vn = len(v)
        assert(vn <= self.n)
        if convert:
            convert_dict = self._parent._domain_to_gap
            for i in range(len(v)):
                self.perm[i] = convert_dict[v[i]] - 1
        else:
            for i, j in enumerate(v):
                self.perm[i] = j - 1

        for i in range(vn, self.n):
            self.perm[i] = i

    cpdef _set_libgap(self, GapElement p):
        r"""
        TESTS::

            sage: S = SymmetricGroup(4)
            sage: S(libgap.eval('(1,2)'))
            (1,2)
            sage: S(libgap.eval('(2,3)'))
            (2,3)
            sage: S(libgap.eval('(1,4)'))
            (1,4)
            sage: S(libgap.eval('(1,2,4)(3,5)') * libgap.eval('(3,5)'))
            (1,2,4)

            sage: S(libgap.eval('(1,5)'))
            Traceback (most recent call last):
            ...
            ValueError: invalid data to initialize a permutation
        """
        cdef const UInt2* p2
        cdef const UInt4* p4
        cdef int i
        cdef UInt d

        if TNUM_OBJ(p.value) == T_PERM2:
            d = DEG_PERM2(p.value)
            if d > self.n:
                d = self.n
            else:
                for i in range(d, self.n):
                    self.perm[i] = i
            p2 = CONST_ADDR_PERM2(p.value)
            for i in range(d):
                self.perm[i] = p2[i]
        elif TNUM_OBJ(p.value) == T_PERM4:
            d = DEG_PERM4(p.value)
            if d > self.n:
                d = self.n
            else:
                for i in range(d, self.n):
                    self.perm[i] = i
            p4 = CONST_ADDR_PERM4(p.value)
            for i in range(d):
                self.perm[i] = p4[i]
        else:
            raise TypeError("not a gap permutation")

        self._libgap = p

    cpdef _set_permutation_group_element(self, PermutationGroupElement p, bint convert):
        r"""
        TESTS::

            sage: S1 = SymmetricGroup(5)
            sage: S2 = SymmetricGroup(6)
            sage: S3 = SymmetricGroup([1,3,5])
            sage: for U in [S1,S2,S3]:
            ....:     for V in [S1,S2,S3]:
            ....:         assert U(V((1,3))) == U((1,3))

            sage: S3(S2('(1,2,3)'))
            Traceback (most recent call last):
            ...
            KeyError: 2

            sage: S1(S2('(1,6)'))
            Traceback (most recent call last):
            ...
            ValueError: invalid data to initialize a permutation

            sage: S0 = SymmetricGroup(0)
            sage: S0(S1.one())
            ()
            sage: S0(S1('(1,5)'))
            Traceback (most recent call last):
            ...
            ValueError: invalid data to initialize a permutation
        """
        cdef i, j
        if convert:
            self._set_identity()
            convert_dict = self._parent._domain_to_gap
            p_domain = p._parent._domain
            for i in range(p.n):
                if p.perm[i] == i:
                    continue
                j = p.perm[i]
                i = convert_dict[p_domain[i]] - 1
                j = convert_dict[p_domain[j]] - 1
                self.perm[i] = j
        else:
            if p.n > self.n:
                for i in range(self.n, p.n):
                    if p.perm[i] != i:
                        raise ValueError("invalid data to initialize a permutation")
                for i in range(self.n):
                    self.perm[i] = p.perm[i]
            else:
                for i in range(p.n):
                    self.perm[i] = p.perm[i]
                for i in range(p.n, self.n):
                    self.perm[i] = i

    cpdef _set_list_cycles(self, c, bint convert):
        r"""
        TESTS::

            sage: p = PermutationGroupElement([3,1,4,2,5])
            sage: p._set_list_cycles([(1,4),(2,3,5)], True)
            sage: p
            (1,4)(2,3,5)
            sage: p._set_list_cycles([(1,4),(2,3,5)], False)
            sage: p
            (1,4)(2,3,5)

            sage: S = SymmetricGroup(["a", "b", "c"])
            sage: p = S(["a", "b", "c"])
            sage: p._set_list_cycles([("a", "c")], True)
            sage: p
            ('a','c')
            sage: p._set_list_cycles([(1, 3)], False)
            sage: p
            ('a','c')
        """
        cdef int i, j
        self._set_identity()
        if convert:
            convert_dict = self._parent._domain_to_gap
            for t in c:
                if len(t) <= 1:
                    continue
                for i in range(len(t) - 1):
                    j = convert_dict[t[i]] - 1
                    self.perm[j] = convert_dict[t[i+1]] - 1
                j = convert_dict[t[-1]] - 1
                self.perm[j] = convert_dict[t[0]] - 1
        else:
            for t in c:
                if len(t) <= 1:
                    continue
                for i in range(len(t) - 1):
                    j = t[i] - 1
                    if j < 0 or j >= self.n:
                        raise ValueError("invalid list of cycles to initialize a permutation")
                    self.perm[j] = t[i+1] - 1
                j = t[-1] - 1
                if j < 0 or j >= self.n:
                    raise ValueError("invalid list of cycles to initialize a permutation")
                self.perm[j] = t[0] - 1

    cpdef _set_string(self, str s):
        r"""
        TESTS::

            sage: s = SymmetricGroup(8)()
            sage: s._set_string('(4,2,1,6)(3,5)')
            sage: s
            (1,6,4,2)(3,5)
            sage: s._set_string('(4,2) (1,6)\n \n(3,5)')
            sage: s
            (1,6)(2,4)(3,5)

            sage: S = SymmetricGroup(['a', 'b', 'c', 'd', 'e', 'f'])
            sage: S("('a','c')('b','e','d')")
            ('a','c')('b','e','d')

            sage: S = SymmetricGroup(8)
            sage: S("()")
            ()
            sage: S("()(1,3)()()")
            (1,3)
            sage: S("(1,3")
            Traceback (most recent call last):
            ...
            ValueError: invalid string to initialize a permutation
            sage: S("1,3)")
            Traceback (most recent call last):
            ...
            ValueError: invalid string to initialize a permutation
            sage: S("1,3")
            Traceback (most recent call last):
            ...
            ValueError: invalid string to initialize a permutation
        """
        cdef Py_ssize_t i=0, j
        cdef int k, m
        cdef str c
        cdef list cycle
        cdef dict convert_dict = self._parent._domain_to_gap
        cdef bint convert = self._parent._has_natural_domain()

        self._set_identity()

        while i < len(s):
            while i < len(s) and s[i].isspace():
                i += 1
            if i == len(s):
                return
            j = s.find(')', i + 1)
            if s[i] != '(' or j == -1:
                raise ValueError("invalid string to initialize a permutation")
            i += 1
            if i == j:
                i = j + 1
                continue

            cycle = []
            for c in s[i:j].split(','):
                c = c.strip()
                if convert:
                    cycle.append(int(c) - 1)
                else:
                    cycle.append(convert_dict[eval(c)] - 1)
            for k in range(len(cycle) - 1):
                m = cycle[k]
                if m < 0 or m >= self.n:
                    raise ValueError("invalid list of cycles to initialize a permutation")
                self.perm[m] = cycle[k+1]
            m = cycle[-1]
            if m < 0 or m >= self.n:
                raise ValueError("invalid list of cycles to initialize a permutation")
            self.perm[m] = cycle[0]

            i = j + 1

    def __reduce__(self):
        r"""
        Return a function and its arguments needed to create this
        permutation group element.  This is used in pickling.

        EXAMPLES::

           sage: g = PermutationGroupElement([(1,2,3),(4,5)]); g
           (1,2,3)(4,5)
           sage: func, args = g.__reduce__()
           sage: func(*args)
           (1,2,3)(4,5)
        """
        return make_permgroup_element_v2, (self._parent, self.domain(), self._parent.domain())

    cdef _alloc(self, int n):
        if n < 16 and self.perm == NULL:
            self.perm = self.perm_buf
        elif self.perm == NULL:
            self.perm = <int *> sig_calloc(n, sizeof(int))
        elif n > self.n:
            self.perm = <int *> sig_realloc(self.perm, n * sizeof(int))

        self.n = n

    cdef PermutationGroupElement _new_c(self):
        cdef type t = type(self)
        cdef PermutationGroupElement other = t.__new__(t)
        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__
        other._parent = self._parent
        other._alloc(self.n)
        return other

    def _gap_(self, gap=None):
        r"""
        TESTS::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)]); g
            (1,2,3)(4,5)
            sage: g._gap_()
            (1,2,3)(4,5)
        """
        return self._libgap_()

    def _libgap_(self):
        r"""
        Return ``self`` as a libgap element.

        EXAMPLES::

            sage: S = SymmetricGroup(4)
            sage: p = S('(2,4)')
            sage: p_libgap = libgap(p)
            sage: p_libgap.Order()
            2
            sage: S(p_libgap) == p
            True

            sage: # needs sage.rings.finite_rings
            sage: P = PGU(8,2)
            sage: p, q = P.gens()
            sage: p_libgap  = p.gap()

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: P = PGU(8,2)
            sage: p, q = P.gens()
            sage: p_pexpect = gap(p)
            sage: p_libgap == p_pexpect
            True
            sage: type(p_libgap) == type(p_pexpect)
            True

        If the permutation element is built from a libgap element, it is cached
        and returned by this function::

            sage: S = SymmetricGroup(4)
            sage: p = libgap.eval('(1,3)')
            sage: libgap(S(p)) is p
            True

        Test the empty permutation::

            sage: p = SymmetricGroup(0).an_element()
            sage: p._libgap_()
            ()
        """
        if self._libgap is not None:
            return self._libgap
        initialize()

        cdef Obj res = NEW_PERM2(self.n)
        cdef UInt2* p = ADDR_PERM2(res)
        cdef UInt i
        for i in range(self.n):
            p[i] = self.perm[i]
        self._libgap = make_GapElement_Permutation(libgap, res)
        return self._libgap

    # for compatibility with sage.groups.libgap_wrapper.ElementLibGAP
    # see sage.groups.perm_gps.permgroup.PermutationGroup_generic.gap
    gap = _libgap_

    def _gap_init_(self):
        r"""
        Return a GAP string representation for this
        PermutationGroupElement.

        EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._gap_init_()
            'PermList([2, 3, 1, 5, 4])'
        """
        return 'PermList(%s)' % self._gap_list()

    def _repr_(self):
        r"""
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
        Return a latex representation of this permutation.

        EXAMPLES::

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: latex(g)
            (1,2,3)(4,5)

            sage: S = SymmetricGroup(['a', 'b'])
            sage: latex(S.gens())
            \left((\text{\texttt{a}},\text{\texttt{b}})\right)
            sage: latex(S.one())
            1
        """
        cycle_tuples = self.cycle_tuples()
        if not cycle_tuples:
            return '1'
        from sage.misc.latex import latex
        return "".join(("(" + ",".join(latex(x) for x in cycle) + ")")
                       for cycle in cycle_tuples)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th permutation cycle in the disjoint cycle
        representation of ``self``.

        INPUT:

        - ``i`` -- integer

        OUTPUT: a permutation group element

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3),(4,5)]])
            sage: g = G.gen(0)
            sage: g[0]
            (1,2,3)
            sage: g[1]
            (4,5)
        """
        return self.cycles()[i]

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare group elements ``self`` and ``other``.

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

        Verify that we fixed bug :issue:`5537`::

            sage: h = PermutationGroupElement('(1,3,2)')
            sage: k = PermutationGroupElement('(1,2,3)(4,5)')
            sage: k^2 == h, h == k^2
            (True, True)
            sage: k^6 == PermutationGroupElement('()')
            True
        """
        cdef int i
        cdef PermutationGroupElement right = <PermutationGroupElement>other
        for i in range(self.n):  # Equal parents, so self.n == other.n
            li = self.perm[i]
            ri = right.perm[i]
            if li != ri:
                return richcmp_not_equal(li, ri, op)
        return rich_to_bool(op, 0)

    def __call__(self, i):
        r"""
        Return the image of the integer `i` under this permutation.
        Alternately, if `i` is a list, tuple or string, returns the result of
        ``self`` acting on `i`.

        EXAMPLES::

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

            sage: g(x)                                                                  # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: must be in the domain or a list, tuple or string
            sage: g(3/2)
            Traceback (most recent call last):
            ...
            ValueError: must be in the domain or a list, tuple or string
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
            if isinstance(i, (int, Integer)):
                return i

            if not isinstance(i, (list, tuple, str)):
                raise ValueError("must be in the domain or a list, tuple or string")

            permuted = [i[self.perm[j]] for j from 0 <= j < self.n]
            if isinstance(i, tuple):
                permuted = tuple(permuted)
            elif isinstance(i, str):
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
        r"""
        Return the right action of ``self`` on the list ``x``. This is the
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
        r"""
        Return the right action of ``self`` on the :class:`ClonableIntArray`
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

    cpdef ETuple _act_on_etuple_on_position(self, ETuple x):
        r"""
        Return the right action of this permutation on the ETuple ``x``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: S = SymmetricGroup(6)
            sage: e = ETuple([1,2,3,4,5,6])
            sage: S("(1,4)")._act_on_etuple_on_position(e)
            (4, 2, 3, 1, 5, 6)
            sage: S("(1,2,3,4,5,6)")._act_on_etuple_on_position(e)
            (6, 1, 2, 3, 4, 5)

            sage: e = ETuple([1,2,0,0,0,6])
            sage: S("(1,4)")._act_on_etuple_on_position(e)
            (0, 2, 0, 1, 0, 6)
            sage: S("(1,2,3,4,5,6)")._act_on_etuple_on_position(e)
            (6, 1, 2, 0, 0, 0)

        It is indeed a right action::

            sage: p, q = S('(1,2,3,4,5,6)'), S('(1,2)(3,4)(5,6)')
            sage: e = ETuple([10..15])
            sage: right = lambda x, p: p._act_on_etuple_on_position(x)
            sage: right(e, p * q) == right(right(e, p), q)
            True
        """
        cdef size_t ind
        cdef ETuple result = ETuple.__new__(ETuple)

        result._length = x._length
        result._nonzero = x._nonzero
        result._data = <int*> sig_malloc(sizeof(int)*result._nonzero*2)
        for ind in range(x._nonzero):
            result._data[2*ind] = self.perm[x._data[2*ind]] # index
            result._data[2*ind + 1] = x._data[2*ind+1] # exponent
        qsort(result._data, result._nonzero, 2 * sizeof(int), etuple_index_cmp)
        return result

    cpdef _act_on_(self, x, bint self_on_left):
        r"""
        Return the result of the action of ``self`` on ``x``.

        For example, if ``x=f(z)`` is a polynomial, then this function returns
        f(sigma\*z), which is the image of f under the right action of sigma on
        the indeterminates. This is a right action since the image of
        f(sigma\*z) under tau is f(sigma\*tau\*z).

        Additionally, if ``x`` is a matrix, then sigma acts on the matrix
        by permuting the columns when acting from the right and by permuting
        the rows when acting from the left.

        INPUT:

        - ``x`` -- element of space on which permutations act

        - ``self_on_left`` -- if ``True``, this permutation acts on ``x`` from
          the left, otherwise from the right

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

            sage: M = matrix(ZZ,[[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]])
            sage: sigma * M
            [0 2 0 0 0]
            [0 0 3 0 0]
            [1 0 0 0 0]
            [0 0 0 0 5]
            [0 0 0 4 0]
            sage: (M * sigma) * tau == M * (sigma * tau)
            True
            sage: (M * sigma) * tau == (M * sigma.matrix()) * tau.matrix()
            True
            sage: (tau * sigma) * M == tau * (sigma * M)
            True
        """
        if not self_on_left:
            left = x
            if isinstance(left, Polynomial):
                if self != 1:
                    raise ValueError("%s does not act on %s" % (self,
                                                                left.parent()))
                return left
            elif isinstance(left, MPolynomial):
                R = left.parent()
                vars = R.gens()
                try:
                    sigma_x  = [vars[self(i+1)-1] for i in range(R.ngens())]
                except IndexError:
                    raise TypeError("%s does not act on %s" % (self,
                                                               left.parent()))
                return left(tuple(sigma_x))
            elif isinstance(left, Matrix):
                return left.with_permuted_columns(~self)
        else:
            if isinstance(x, Matrix):
                return x.with_permuted_rows(self)

    def __mul__(left, right):
        r"""
        TESTS::

            sage: S = SymmetricGroup(5)
            sage: P = PermutationGroup([(1,2,3),(4,5)])
            sage: Q = PermutationGroup([(1,2),(2,3),(4,5)])
            sage: prod = S('(1,2,3)(4,5)')
            sage: for P1 in [S, P,Q]:
            ....:     for P2 in [S, P,Q]:
            ....:         prod = P1('(1,2,3)') * P2('(4,5)')
            ....:         assert prod.parent() == coercion_model.common_parent(P1, P2)
            ....:         assert prod == S('(1,2,3)(4,5)')
        """
        if type(left) is type(right) and \
           (<PermutationGroupElement> left)._parent is (<PermutationGroupElement> right)._parent:
            return (<PermutationGroupElement> left)._mul_(right)

        # shortcut the case when one of them belong to a full symmetric group
        # and domains are equal
        cdef PermutationGroupElement pl, pr, prod
        cdef int i

        if isinstance(left, SymmetricGroupElement):
            if isinstance(right, PermutationGroupElement):
                pl = <PermutationGroupElement> left
                pr = <PermutationGroupElement> right
                Pleft = pl._parent
                Pright = pr._parent
                if Pleft._domain == Pright._domain:
                    prod = pl._new_c()
                    for i in range(pl.n):
                        prod.perm[i] = pr.perm[pl.perm[i]]
                    return prod

        elif isinstance(right, SymmetricGroupElement):
            pl = <PermutationGroupElement> left
            pr = <PermutationGroupElement> right
            Pleft = pl._parent
            Pright = pr._parent
            if Pleft._domain == Pright._domain:
                prod = pr._new_c()
                for i in range(pr.n):
                    prod.perm[i] = pr.perm[pl.perm[i]]
                return prod

        return coercion_model.bin_op(left, right, operator.mul)

    cpdef PermutationGroupElement _transpose_left(self, j, k):
        r"""
        Return the product of the transposition `(j, k)` with ``self``.

        EXAMPLES::

            sage: S = SymmetricGroup(5)
            sage: s = S([5, 2, 4, 3, 1])
            sage: s._transpose_left(2, 3)
            (1,5)(2,4,3)
            sage: S((2, 3)) * s
            (1,5)(2,4,3)

            sage: S = SymmetricGroup(["a", "b", "c", "d", "e"])
            sage: s = S(["e", "b", "d", "c", "a"])
            sage: s._transpose_left("b", "c")
            ('a','e')('b','d','c')
            sage: S(("b", "c")) * s
            ('a','e')('b','d','c')
        """
        if j == k:
            return self
        cdef PermutationGroupElement prod = self._new_c()
        cdef int i
        for i in range(self.n):
            prod.perm[i] = self.perm[i]
        if not self._parent._has_natural_domain():
            convert_dict = self._parent._domain_to_gap
            j = convert_dict[j]
            k = convert_dict[k]
        prod.perm[j - 1], prod.perm[k - 1] = self.perm[k - 1], self.perm[j - 1]
        return prod

    cpdef _mul_(left, _right):
        r"""
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

    cpdef PermutationGroupElement _generate_new(self, list v):
        r"""
        Generate a new permutation group element with the same parent
        as ``self`` from ``v``.

        EXAMPLES::

            sage: P = PermutationGroup([(1,2),(1,2,3,4)])
            sage: one = P.one()
            sage: one._generate_new([])
            ()
            sage: one._generate_new([4,3,2,1])
            (1,4)(2,3)
        """
        cdef PermutationGroupElement new = self._new_c()
        new._set_list_images(v, False)
        return new

    cpdef PermutationGroupElement _generate_new_GAP(self, lst_in):
        r"""
        Generate a new permutation group element with the same parent
        as ``self`` from the GAP list ``lst_in``.

        EXAMPLES::

            sage: from sage.libs.gap.libgap import libgap

            sage: P = PermutationGroup([(1,2),(1,2,3,4)])
            sage: one = P.one()
            sage: perm = libgap.eval('[]')
            sage: one._generate_new_GAP(perm)
            ()
            sage: perm = libgap.eval('[4,3,2,1]')
            sage: one._generate_new_GAP(perm)
            (1,4)(2,3)
        """
        cdef GapElement_List lst = <GapElement_List?> lst_in
        cdef Obj obj = lst.value

        cdef PermutationGroupElement new = self._new_c()
        cdef Py_ssize_t i, j, vn = len(lst)

        assert vn <= self.n

        for i in range(vn):
            j = GAP_ValueInt(ELM_LIST(obj, i+1))
            new.perm[i] = j - 1
        for i in range(vn, self.n):
            new.perm[i] = i
        return new

    def __invert__(self):
        r"""
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
        r"""
        Return this permutation in list notation compatible with the
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
        r"""
        Return a cycle string for this permutation compatible with
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

    cpdef domain(self):
        r"""
        Return the domain of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: x = G([2,1,4,3]); x
            (1,2)(3,4)
            sage: v = x.domain(); v
            [2, 1, 4, 3]
            sage: type(v[0])
            <... 'int'>
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
        if self.n == 0:
            return []

        cdef int i
        from_gap = self._parent._domain_from_gap
        return [from_gap[self.perm[i]+1] for i from 0 <= i < self.n]

    def __hash__(self):
        r"""
        Return a hash for this permutation.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: hash(G([2,1,5,3,4]))
            -1203337681           # 32-bit
            -1527414595000039889  # 64-bit

        Check that the hash looks reasonable::

            sage: s = set()
            sage: s.update(map(hash,SymmetricGroup(0)))
            sage: s.update(map(hash,SymmetricGroup(1)))
            sage: s.update(map(hash,SymmetricGroup(2)))
            sage: s.update(map(hash,SymmetricGroup(3)))
            sage: s.update(map(hash,SymmetricGroup(4)))
            sage: s.update(map(hash,SymmetricGroup(5)))
            sage: len(s) == 1 + 1 + 2 + 6 + 24 + 120
            True
        """
        cdef size_t i
        cdef long ans = self.n
        for i in range(self.n):
            ans = (ans ^ (self.perm[i])) * 1000003L
        return ans

    def tuple(self):
        r"""
        Return tuple of images of the domain under ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: s = G([2,1,5,3,4])
            sage: s.tuple()
            (2, 1, 5, 3, 4)

            sage: S = SymmetricGroup(['a', 'b'])
            sage: S.gen().tuple()
            ('b', 'a')
        """
        if self.n == 0:
            return ()

        cdef int i
        from_gap = self._parent._domain_from_gap
        return tuple([from_gap[self.perm[i]+1] for i in range(self.n)])

    def dict(self):
        r"""
        Return a dictionary associating each element of the domain with its
        image.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: g = G((1,2,3,4)); g
            (1,2,3,4)
            sage: v = g.dict(); v
            {1: 2, 2: 3, 3: 4, 4: 1}
            sage: type(v[1])
            <... 'int'>
            sage: x = G([2,1]); x
            (1,2)
            sage: x.dict()
            {1: 2, 2: 1, 3: 3, 4: 4}
        """
        from_gap = self._parent._domain_from_gap
        to_gap = self._parent._domain_to_gap
        cdef int i
        return {e:from_gap[self.perm[i-1]+1] for e,i in to_gap.iteritems()}

    def multiplicative_order(self):
        r"""
        Return the order of this group element, which is the smallest
        positive integer `n` for which `g^n = 1`.

        EXAMPLES::

            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.multiplicative_order()
            6

        :meth:`order` is just an alias for :meth:`multiplicative_order`::

            sage: s.order()
            6

        TESTS::

            sage: # needs sage.libs.pari
            sage: prod(primes(150))
            1492182350939279320058875736615841068547583863326864530410
            sage: L = [tuple(range(sum(primes(p))+1, sum(primes(p))+1+p)) for p in primes(150)]
            sage: t = PermutationGroupElement(L).multiplicative_order(); t
            1492182350939279320058875736615841068547583863326864530410
            sage: type(t)
            <class 'sage.rings.integer.Integer'>
        """
        order = None
        cdef long long order_c = 1
        cdef int cycle_len
        cdef int i, k
        cdef bint* seen = <bint *>sig_malloc(sizeof(bint) * self.n)
        for i in range(self.n):
            seen[i] = 0
        for i in range(self.n):
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
                if order_c > LLONG_MAX / (self.n - i):
                    order = Integer(order_c)
        sig_free(seen)
        return Integer(order_c) if order is None else order

    def inverse(self):
        r"""
        Return the inverse permutation.

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
        r"""
        Return the sign of ``self``, which is `(-1)^{s}`, where
        `s` is the number of swaps.

        EXAMPLES::

            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.sign()
            -1

        ALGORITHM: Only even cycles contribute to the sign, thus

        .. MATH::

            sign(sigma) = (-1)^{\sum_c len(c)-1}


        where the sum is over cycles in ``self``.
        """
        cdef int cycle_len_sum = 0
        cdef int i, k
        cdef bint* seen = <bint *>sig_malloc(sizeof(bint) * self.n)
        for i in range(self.n):
            seen[i] = 0
        for i in range(self.n):
            if seen[i] or self.perm[i] == i:
                continue
            k = self.perm[i]
            while k != i:
                seen[k] = 1
                k = self.perm[k]
                cycle_len_sum += 1
        sig_free(seen)
        return 1 - 2*(cycle_len_sum % 2) # == (-1)^cycle_len

    def orbit(self, n, bint sorted=True):
        r"""
        Return the orbit of the integer `n` under this group
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
        r"""
        Return ``self`` as a list of disjoint cycles.

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
        cdef bint* seen = <bint *>sig_malloc(sizeof(bint) * self.n)
        for i in range(self.n):
            seen[i] = 0
        for i in range(self.n):
            if seen[i] or self.perm[i] == i:
                continue
            cycle = self._new_c()
            for j in range(self.n):
                cycle.perm[j] = j
            k = cycle.perm[i] = self.perm[i]
            while k != i:
                seen[k] = 1
                next_k = cycle.perm[k] = self.perm[k]
                k = next_k
            PyList_Append(L, cycle)
        sig_free(seen)
        return L

    def cycle_tuples(self, singletons=False):
        r"""
        Return ``self`` as a list of disjoint cycles, represented as tuples
        rather than permutation group elements.

        INPUT:

        - ``singletons`` -- boolean (default: ``False``); whether or not consider the
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
        cdef bint* seen = <bint *>sig_malloc(sizeof(bint) * self.n)
        for i in range(self.n):
            seen[i] = 0
        for i in range(self.n):
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
        sig_free(seen)
        return L

    def cycle_string(self, singletons=False):
        r"""
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
        if not cycles:
            return '()'
        text = ''.join(repr(c) for c in cycles)
        return text.replace(', ', ',').replace(',)', ')')

    def cycle_type(self, singletons=True, as_list=False):
        r"""
        Return the partition that gives the cycle type of ``self`` as an element of
        its parent.

        INPUT:

        - ``g`` -- an element of the permutation group ``self.parent()``

        - ``singletons`` -- boolean depending on whether or not
          trivial cycles should be counted (default: ``True``)

        - ``as_list`` -- boolean depending on whether the cycle
          type should be returned as a :class:`list` or as a :class:`Partition`
          (default: ``False``)

        OUTPUT: a :class:`Partition`, or :class:`list` if ``is_list`` is
        ``True``, giving the cycle type of ``self``

        If speed is a concern, then ``as_list=True`` should be used.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: G = DihedralGroup(3)
            sage: [g.cycle_type() for g in G]
            [[1, 1, 1], [3], [3], [2, 1], [2, 1], [2, 1]]
            sage: PermutationGroupElement('(1,2,3)(4,5)(6,7,8)').cycle_type()
            [3, 3, 2]
            sage: G = SymmetricGroup(3); G('(1,2)').cycle_type()
            [2, 1]
            sage: G = SymmetricGroup(4); G('(1,2)').cycle_type()
            [2, 1, 1]
            sage: G = SymmetricGroup(4); G('(1,2)').cycle_type(singletons=False)
            [2]
            sage: G = SymmetricGroup(4); G('(1,2)').cycle_type(as_list=False)
            [2, 1, 1]
        """
        cycle_type = [len(c) for c in self.cycle_tuples(singletons)]
        cycle_type.sort(reverse=True)
        if as_list:
            return cycle_type
        else:
            from sage.combinat.partition import _Partitions
            return _Partitions(cycle_type)

    def has_descent(self, i, side='right', positive=False):
        r"""
        Return whether ``self`` has a left (resp. right) descent at
        position ``i``. If ``positive`` is ``True``, then test for a non
        descent instead.

        Beware that, since permutations are acting on the right, the
        meaning of descents is the reverse of the usual
        convention. Hence, ``self`` has a left descent at position
        ``i`` if ``self(i) > self(i+1)``.

        INPUT:

        - ``i`` -- an element of the index set
        - ``side`` -- ``'left'`` or ``'right'`` (default: ``'right'``)
        - ``positive`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: S = SymmetricGroup([1,2,3])
            sage: S.one().has_descent(1)
            False
            sage: S.one().has_descent(2)
            False
            sage: s = S.simple_reflections()
            sage: x = s[1]*s[2]
            sage: x.has_descent(1, side='right')
            False
            sage: x.has_descent(2, side='right')
            True
            sage: x.has_descent(1, side='left')
            True
            sage: x.has_descent(2, side='left')
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
        r"""
        Return a deg `\times` deg permutation matrix associated to the permutation
        ``self``.

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
        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(ZZ, self.n, self.n, sparse=True)
        cdef int i
        entries = {}
        for i from 0 <= i < self.n:
            entries[i, self.perm[i]] = 1
        return M(entries)

    def word_problem(self, words, display=True, as_list=False):
        r"""
        Try to solve the word problem for ``self``.

        INPUT:

        - ``words`` -- list of elements of the ambient group, generating
          a subgroup

        - ``display`` -- boolean (default: ``True``); whether to display
          additional information

        - ``as_list`` -- boolean (default: ``False``); whether to return
          the result as a list of pairs (generator, exponent)

        OUTPUT:

        - a pair of strings, both representing the same word

        or

        - a list of pairs representing the word, each pair being
          (generator as a string, exponent as an integer)

        Let `G` be the ambient permutation group, containing the given
        element `g`. Let `H` be the subgroup of `G` generated by the list
        ``words`` of elements of `G`. If `g` is in `H`, this function
        returns an expression for `g` as a word in the elements of
        ``words`` and their inverses.

        This function does not solve the word problem in Sage. Rather it
        pushes it over to GAP, which has optimized algorithms for the word
        problem. Essentially, this function is a wrapper for the GAP
        functions ``EpimorphismFromFreeGroup`` and
        ``PreImagesRepresentative``.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]], canonicalize=False)
            sage: g1, g2 = G.gens()
            sage: h = g1^2*g2*g1
            sage: h.word_problem([g1,g2], False)
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')

            sage: h.word_problem([g1,g2])
               x1^2*x2^-1*x1
               [['(1,2,3)(4,5)', 2], ['(3,4)', -1], ['(1,2,3)(4,5)', 1]]
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')

            sage: h.word_problem([g1,g2], False, as_list=True)
            [['(1,2,3)(4,5)', 2], ['(3,4)', -1], ['(1,2,3)(4,5)', 1]]

        TESTS:

        Check for :issue:`28556`::

            sage: G = SymmetricGroup(6)
            sage: g = G('(1,2,3)')
            sage: g.word_problem([g], False)
            ('x1', '(1,2,3)')
        """
        if not self._parent._has_natural_domain():
            raise NotImplementedError

        def convert_back(string):
            L = copy.copy(string)
            for i, w_i in enumerate(words):
                L = L.replace("x" + str(i + 1), str(w_i))
            return L

        g = words[0].parent()(self)
        H = libgap.Group(words)
        ans = H.EpimorphismFromFreeGroup().PreImagesRepresentative(g)

        l1 = str(ans)
        l2 = convert_back(l1)

        if display or as_list:
            l3 = l1.split("*")
            l4 = []
            for m in l3:  # parsing the word for display
                m_split = m.split("^")
                if len(m_split) == 2:
                    l4.append([m_split[0], int(m_split[1])])
                else:
                    l4.append([m_split[0], 1])
            l5 = [[convert_back(w), e] for w, e in l4]

        if display:
            print(l1)
            print(l5)

        if as_list:
            return l5
        else:
            return l1, l2


cdef class SymmetricGroupElement(PermutationGroupElement):
    r"""
    An element of the symmetric group.
    """
    def absolute_length(self):
        r"""
        Return the absolute length of ``self``.

        The absolute length is the size minus the number of its disjoint
        cycles. Alternatively, it is the length of the shortest
        expression of the element as a product of reflections.

        .. SEEALSO::

            :meth:`absolute_le`

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: [x.absolute_length() for x in S]                                      # needs sage.combinat
            [0, 2, 2, 1, 1, 1]
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self).absolute_length()

    def has_left_descent(self, i):
        r"""
        Return whether `i` is a left descent of ``self``.

        EXAMPLES::

            sage: W = SymmetricGroup(4)
            sage: w = W.from_reduced_word([1,3,2,1])
            sage: [i for i in W.index_set() if w.has_left_descent(i)]
            [1, 3]
        """
        return self.has_descent(i, side='left')


cdef bint is_valid_permutation(int* perm, int n) noexcept:
    r"""
    This is used in the ``__init__`` method.

    Return ``True`` iff the first ``n`` elements of perm are literally a
    permutation of ``[0, ..., n-1]``.

    TESTS::

        sage: S = SymmetricGroup(10)
        sage: PermutationGroupElement([2,1],S,check=False)
        (1,2)
        sage: PermutationGroupElement([1,1],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: invalid data to initialize a permutation
        sage: PermutationGroupElement([1,-1],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: invalid data to initialize a permutation
        sage: PermutationGroupElement([1,2,3,10],S,check=False)
        Traceback (most recent call last):
        ...
        ValueError: invalid data to initialize a permutation
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
