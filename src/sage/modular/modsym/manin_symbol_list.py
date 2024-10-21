# sage_setup: distribution = sagemath-schemes
"""
Manin symbol lists

There are various different classes holding lists of Manin symbols of
different types.  The hierarchy is as follows:

- :class:`ManinSymbolList`

  - :class:`ManinSymbolList_group`

    - :class:`ManinSymbolList_gamma0`
    - :class:`ManinSymbolList_gamma1`
    - :class:`ManinSymbolList_gamma_h`

  - :class:`ManinSymbolList_character`
"""
# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.modular.modsym.p1list as p1list
import sage.modular.modsym.g1list as g1list
import sage.modular.modsym.ghlist as ghlist
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.misc.persist import register_unpickle_override
from sage.structure.richcmp import richcmp_method, richcmp
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from .apply import apply_to_monomial

from sage.modular.modsym.manin_symbol import ManinSymbol


@richcmp_method
class ManinSymbolList(Parent):
    """
    Base class for lists of all Manin symbols for a given weight, group or character.
    """

    Element = ManinSymbol

    def __init__(self, weight, lst):
        """
        Constructor for a ManinSymbolList.

        INPUT:

        - ``weight`` -- the weight of the symbols

        - ``lst`` -- the list of symbols

        On construction, a ManinSymbolList constructs a dict for
        rapid determination of the index of any given symbol.

        This is a base class only; users will only directly construct
        objects in the derived classes ManinSymbolList_gamma0,
        ManinSymbolList_gamma1, ManinSymbolList_gamma_h,
        ManinSymbolList_gamma_character.  Many standard methods are
        only implemented in the derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: ManinSymbolList(6,P1List(11))
            <sage.modular.modsym.manin_symbol_list.ManinSymbolList_with_category object at ...>
        """
        self._weight = weight
        self._symbol_list = lst
        self._index = {x: i for i,x in enumerate(lst)}
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        TESTS::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6, P1List(11))
            sage: x = m((2, 3, 5)); x
            [X^2*Y^2,(3,5)]
            sage: m(x) == x
            True
        """
        if isinstance(x, ManinSymbol):
            x = x.tuple()
        return self.element_class(self, x)

    def __richcmp__(self, right, op):
        """
        Comparison function for ManinSymbolList objects.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m1 = ManinSymbolList(6,P1List(11))
            sage: m2 = ManinSymbolList(6,P1List(13))
            sage: m3 = ManinSymbolList(4,P1List(11))
            sage: m1 < m2
            True
            sage: m2 < m3
            False
            sage: m1 < m3
            False
        """
        if not isinstance(right, ManinSymbolList):
            return NotImplemented
        return richcmp((self._weight, self._symbol_list),
                       (right._weight, right._symbol_list), op)

    def symbol_list(self):
        """
        Return the list of symbols of ``self``.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6, P1List(11))
        """
        return list(self._symbol_list)  # This makes a shallow copy

    def __len__(self):
        """
        Return the length of this :class:`ManinSymbolList`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: len(m)
            12
        """
        return len(self._symbol_list)

    def apply(self, j, X):
        """
        Apply the matrix `X = [a, b; c, d]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply(10, [1,2,0,1])
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError("Only implemented in derived classes")

    def _apply_S_only_0pm1(self) -> bool:
        """
        Return ``True`` if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            True
            sage: eps = DirichletGroup(7).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            False
        """
        return False  # derived classes could overload and put True

    def apply_S(self, j):
        """
        Apply the matrix `S = [0, -1; 1, 0]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_S(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError("Only implemented in derived classes")

    def apply_I(self, j):
        """
        Apply the matrix `I = [-1, 0; 0, 1]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_I(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError("Only implemented in derived classes")

    def apply_T(self, j):
        """
        Apply the matrix `T = [0, 1; -1, -1]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_T(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError("Only implemented in derived classes")

    def apply_TT(self, j):
        """
        Apply the matrix `TT = T^2 = [-1, -1; 0, 1]` to the `j`-th
        Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_TT(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError("Only implemented in derived classes")

    def index(self, x):
        """
        Return the index of ``x`` in the list of Manin symbols.

        INPUT:

        - ``x`` -- a triple of integers `(i, u, v)` defining a valid
          Manin symbol, which need not be normalized

        OUTPUT:

        integer -- the index of the normalized Manin symbol equivalent
        to `(i, u, v)`.  If ``x`` is not in ``self``, -1 is returned.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.index(m.symbol_list()[2])
            2
            sage: S = m.symbol_list()
            sage: all(i == m.index(S[i]) for i in range(len(S)))
            True
        """
        if x in self._index:
            return self._index[x]
        x = self.normalize(x)
        try:
            return self._index[x]
        except KeyError:
            return -1

    def manin_symbol_list(self):
        """
        Return all the Manin symbols in ``self`` as a list.

        Cached for subsequent calls.

        OUTPUT:

        A list of :class:`ManinSymbol` objects, which is a copy of the
        complete list of Manin symbols.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.manin_symbol_list() # not implemented for the base class

        ::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: m.manin_symbol_list()
            [[Y^2,(0,1)],
            [Y^2,(1,0)],
            [Y^2,(1,1)],
            ...
            [X^2,(3,1)],
            [X^2,(3,2)]]
        """
        import copy
        try:
            return copy.copy(self.__manin_symbol_list)
        except AttributeError:
            self.__manin_symbol_list = [self.manin_symbol(i)
                                        for i in range(len(self))]
        return copy.copy(self.__manin_symbol_list)

    list = manin_symbol_list

    def manin_symbol(self, i):
        """
        Return the ``i``-th Manin symbol in this :class:`ManinSymbolList`.

        INPUT:

        - ``i`` -- integer; a valid index of a symbol in this list

        OUTPUT:

        :class:`ManinSymbol` -- the `i`-th Manin symbol in the list.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.manin_symbol(3) # not implemented for base class

        ::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: s = m.manin_symbol(3); s
            [Y^2,(1,2)]
            sage: type(s)
            <class 'sage.modular.modsym.manin_symbol.ManinSymbol'>
        """
        return self.element_class(self, self._symbol_list[i])

    def normalize(self, x):
        """
        Return a normalized Manin symbol from ``x``.

        To be implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.normalize((0,6,7)) # not implemented in base class
        """
        raise NotImplementedError("Only implemented in derived classes")

    def weight(self):
        """
        Return the weight of the Manin symbols in this :class:`ManinSymbolList`.

        OUTPUT: integer; the weight of the Manin symbols in the list

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: m.weight()
            4
        """
        return self._weight


class ManinSymbolList_group(ManinSymbolList):
    """
    Base class for Manin symbol lists for a given group.

    INPUT:

    - ``level`` -- integer level

    - ``weight`` -- integer weight

    - ``syms`` -- something with ``normalize`` and ``list`` methods,
      e.g. :class:`~sage.modular.modsym.p1list.P1List`

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_group
        sage: ManinSymbolList_group(11, 2, P1List(11))
        <sage.modular.modsym.manin_symbol_list.ManinSymbolList_group_with_category object at ...>
    """
    def __init__(self, level, weight, syms):
        """
        Constructor for class ManinSymbolList_group.

        INPUT:

        - ``level`` -- integer level

        - ``weight`` -- integer weight

        - ``syms`` -- something with ``normalize`` and ``list``
          methods, e.g. :class:`~sage.modular.modsym.p1list.P1List`

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_group
            sage: L = ManinSymbolList_group(11, 2, P1List(11)); L
            <sage.modular.modsym.manin_symbol_list.ManinSymbolList_group_with_category object at ...>
        """
        self.__level = level
        self.__syms = syms  # syms is anything with a normalize and list method.

        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight - 2 + 1)
             for u, v in syms.list()]
        ManinSymbolList.__init__(self, weight, L)

    def level(self):
        """
        Return the level of this :class:`ManinSymbolList`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: ManinSymbolList_gamma0(5,2).level()
            5

        ::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma1
            sage: ManinSymbolList_gamma1(51,2).level()
            51

        ::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma_h
            sage: ManinSymbolList_gamma_h(GammaH(117, [4]),2).level()
            117
        """
        return self.__level

    def apply_S(self, j):
        """
        Apply the matrix `S = [0, -1; 1, 0]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        ``(k, s)`` where k is the index of the symbol obtained by acting on the
        `j`-th symbol with `S`, and `s` is the parity of the `j`-th symbol
        (a Python ``int``, either 1 or -1).

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_S(4)
            (40, 1)
            sage: [m.apply_S(i) for i in range(len(m))]
            [(37, 1),
            (36, 1),
            (41, 1),
            (39, 1),
            (40, 1),
            (38, 1),
            (31, -1),
            (30, -1),
            (35, -1),
            (33, -1),
            (34, -1),
            (32, -1),
            ...
            (4, 1),
            (2, 1)]
        """
        i, u, v = self._symbol_list[j]
        k = self.index((self._weight - 2 - i, v, -u))
        return k, -1 if i % 2 else 1

    def _apply_S_only_0pm1(self):
        """
        Return ``True`` if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: ManinSymbolList_gamma0(5,8)._apply_S_only_0pm1()
            True
        """
        return True

    def apply_I(self, j):
        """
        Apply the matrix `I=[-1,0,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        ``(k, s)`` where k is the index of the symbol obtained by acting on the
        `j`-th symbol with `I`, and `s` is the parity of the `j`-th symbol
        (a Python ``int``, either 1 or -1)

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_I(4)
            (3, 1)
            sage: [m.apply_I(i) for i in range(10)]
            [(0, 1),
            (1, 1),
            (5, 1),
            (4, 1),
            (3, 1),
            (2, 1),
            (6, -1),
            (7, -1),
            (11, -1),
            (10, -1)]
        """
        i, u, v = self._symbol_list[j]
        k = self.index((i, -u, v))
        if i % 2 == 0:
            return k, 1
        else:
            return k, -1

    def apply_T(self, j):
        """
        Apply the matrix `T=[0,1,-1,-1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT: see documentation for apply()

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_T(4)
            [(3, 1), (9, -6), (15, 15), (21, -20), (27, 15), (33, -6), (39, 1)]
            sage: [m.apply_T(i) for i in range(10)]
            [[(5, 1), (11, -6), (17, 15), (23, -20), (29, 15), (35, -6), (41, 1)],
            [(0, 1), (6, -6), (12, 15), (18, -20), (24, 15), (30, -6), (36, 1)],
            [(4, 1), (10, -6), (16, 15), (22, -20), (28, 15), (34, -6), (40, 1)],
            [(2, 1), (8, -6), (14, 15), (20, -20), (26, 15), (32, -6), (38, 1)],
            [(3, 1), (9, -6), (15, 15), (21, -20), (27, 15), (33, -6), (39, 1)],
            [(1, 1), (7, -6), (13, 15), (19, -20), (25, 15), (31, -6), (37, 1)],
            [(5, 1), (11, -5), (17, 10), (23, -10), (29, 5), (35, -1)],
            [(0, 1), (6, -5), (12, 10), (18, -10), (24, 5), (30, -1)],
            [(4, 1), (10, -5), (16, 10), (22, -10), (28, 5), (34, -1)],
            [(2, 1), (8, -5), (14, 10), (20, -10), (26, 5), (32, -1)]]
        """
        k = self._weight
        i, u, v = self._symbol_list[j]
        u, v = self.__syms.normalize(v,-u-v)
        if (k-2) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        a = Integer(k-2-i)
        for j in range(k-2-i+1):
            m = self.index((j, u, v))
            z.append((m, s * a.binomial(j)))
            s *= -1
        return z

    def apply_TT(self, j):
        """
        Apply the matrix `TT=[-1,-1,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT: see documentation for apply()

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_TT(4)
            [(38, 1)]
            sage: [m.apply_TT(i) for i in range(10)]
            [[(37, 1)],
            [(41, 1)],
            [(39, 1)],
            [(40, 1)],
            [(38, 1)],
            [(36, 1)],
            [(31, -1), (37, 1)],
            [(35, -1), (41, 1)],
            [(33, -1), (39, 1)],
            [(34, -1), (40, 1)]]
        """
        k = self._weight
        i, u, v = self._symbol_list[j]
        u, v = self.__syms.normalize(-u-v,u)
        if (k-2-i) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        a = Integer(i)
        for j in range(i+1):
            m = self.index((k-2-i+j, u, v))
            z.append((m, s * a.binomial(j)))
            s *= -1
        return z

    def apply(self, j, m):
        r"""
        Apply the matrix `m = [a, b; c, d]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        - ``m = [a, b, c, d]`` a list of 4 integers, which defines a 2x2 matrix

        OUTPUT:

        a list of pairs `(j_i, \alpha_i)`, where each `\alpha_i` is a nonzero
        integer, `j_i` is an integer (index of the `j_i`-th Manin symbol), and
        `\sum_i \alpha_i\*x_{j_i}` is the image of the `j`-th Manin symbol under
        the right action of the matrix [a,b;c,d]. Here the right action of
        `g = [a, b; c, d]` on a Manin symbol `[P(X,Y),(u,v)]` is
        `[P(aX+bY,cX+dY),(u,v)\*g]`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply(40, [2,3,1,1])
            [(0, 729), (6, 2916), (12, 4860), (18, 4320),
             (24, 2160), (30, 576), (36, 64)]
        """
        a, b, c, d = m[0], m[1], m[2], m[3]
        i, u, v = self._symbol_list[j]
        P = apply_to_monomial(i, self._weight-2, a, b, c, d)
        m = self.index((0, u*a+v*c, u*b+v*d))
        if m == -1:
            return []
        r = len(self.__syms)
        return [(m + r*k, P[k]) for k in range(self._weight-2+1)
                if P[k] != 0]

    def normalize(self, x):
        """
        Return the normalization of the Manin symbol ``x`` with
        respect to this list.

        INPUT:

        - ``x`` -- (3-tuple of ints) a tuple defining a ManinSymbol

        OUTPUT:

        ``(i,u,v)`` -- (3-tuple of ints) another tuple defining the associated
        normalized ManinSymbol

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: [m.normalize(s.tuple()) for s in m.manin_symbol_list()][:10]
            [(0, 0, 1),
             (0, 1, 0),
             (0, 1, 1),
             (0, 1, 2),
             (0, 1, 3),
             (0, 1, 4),
             (1, 0, 1),
             (1, 1, 0),
             (1, 1, 1),
             (1, 1, 2)]
        """
        u,v = self.__syms.normalize(x[1],x[2])
        return (x[0],u,v)


class ManinSymbolList_gamma0(ManinSymbolList_group):
    r"""
    Class for Manin symbols for `\Gamma_0(N)`.

    INPUT:

    - ``level`` -- integer; the level

    - ``weight`` -- integer; the weight

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
        sage: m = ManinSymbolList_gamma0(5,2); m
        Manin Symbol List of weight 2 for Gamma0(5)
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (1,4)]
        sage: m = ManinSymbolList_gamma0(6,4); m
        Manin Symbol List of weight 4 for Gamma0(6)
        sage: len(m)
        36
    """
    def __init__(self, level, weight):
        """
        Constructor for a ModularSymbolList for Gamma_0(N).

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: M11 = ManinSymbolList_gamma0(11,2)
            sage: M11
            Manin Symbol List of weight 2 for Gamma0(11)
            sage: M11 == loads(dumps(M11))
            True
        """
        ManinSymbolList_group.__init__(self, level, weight, p1list.P1List(level))

    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
            sage: M11 = ManinSymbolList_gamma0(11,2)
            sage: str(M11)
            'Manin Symbol List of weight 2 for Gamma0(11)'
        """
        return "Manin Symbol List of weight %s for Gamma0(%s)" % (
            self.weight(), self.level())


class ManinSymbolList_gamma1(ManinSymbolList_group):
    r"""
    Class for Manin symbols for `\Gamma_1(N)`.

    INPUT:

    - ``level`` -- integer; the level

    - ``weight`` -- integer; the weight

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma1
        sage: m = ManinSymbolList_gamma1(5,2); m
        Manin Symbol List of weight 2 for Gamma1(5)
        sage: m.manin_symbol_list()
        [(0,1),
        (0,2),
        (0,3),
        ...
        (4,3),
        (4,4)]
        sage: m = ManinSymbolList_gamma1(6,4); m
        Manin Symbol List of weight 4 for Gamma1(6)
        sage: len(m)
        72
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, level, weight):
        r"""
        Constructor for a ModularSymbolList for `\Gamma_0(N)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma1
            sage: M11 = ManinSymbolList_gamma1(11,2)
            sage: M11
            Manin Symbol List of weight 2 for Gamma1(11)
        """
        ManinSymbolList_group.__init__(self, level, weight, g1list.G1list(level))

    def __repr__(self):
        """
        Return the string representation of this :class:`ManinSymbolList`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma1
            sage: M11 = ManinSymbolList_gamma1(11,4)
            sage: str(M11)
            'Manin Symbol List of weight 4 for Gamma1(11)'
        """
        return "Manin Symbol List of weight %s for Gamma1(%s)" % (
            self.weight(), self.level())


class ManinSymbolList_gamma_h(ManinSymbolList_group):
    r"""
    Class for Manin symbols for `\Gamma_H(N)`.

    INPUT:

    - ``group`` -- integer; the congruence subgroup

    - ``weight`` -- integer; the weight

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma_h
        sage: G = GammaH(117, [4])
        sage: m = ManinSymbolList_gamma_h(G,2); m
        Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(117) with H generated by [4]
        sage: m.manin_symbol_list()[100:110]
        [(1,88),
        (1,89),
        (1,90),
        (1,91),
        (1,92),
        (1,93),
        (1,94),
        (1,95),
        (1,96),
        (1,97)]
        sage: len(m.manin_symbol_list())
        2016
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, group, weight):
        r"""
        Constructor for Manin symbols for `\Gamma_H(N)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma_h
            sage: G = GammaH(117, [4])
            sage: m = ManinSymbolList_gamma_h(G,2); m
            Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(117) with H generated by [4]
        """
        self.__group = group
        ManinSymbolList_group.__init__(self, group.level(), weight, ghlist.GHlist(group))

    def group(self):
        """
        Return the group associated to ``self``.

        EXAMPLES::

            sage: ModularSymbols(GammaH(12, [5]), 2).manin_symbols().group()
            Congruence Subgroup Gamma_H(12) with H generated by [5]
        """
        return self.__group

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularSymbols(GammaH(12, [5]), 2).manin_symbols().__repr__()
            'Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(12) with H generated by [5]'
        """
        return "Manin Symbol List of weight %s for %s" % (
            self.weight(), self.group())


class ManinSymbolList_character(ManinSymbolList):
    """
    List of Manin symbols with character.

    INPUT:

    - ``character`` -- (DirichletCharacter) the Dirichlet character

    - ``weight`` -- integer; the weight

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: eps = DirichletGroup(4).gen(0)
        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
        sage: m = ManinSymbolList_character(eps,2); m
        Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (2,1)]
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, character, weight):
        """
        Constructor for :class:`ManinSymbolList_character` objects.

        INPUT:

        - ``character`` -- (DirichletCharacter) the Dirichlet character

        - ``weight`` -- integer; the weight

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.manin_symbol_list()
            [(0,1), (1,0), (1,1), (1,2), (1,3), (2,1)]
            sage: TestSuite(m).run()
        """
        self.__level = character.modulus()
        self.__P1 = p1list.P1List(self.level())

        # We make a copy of the character *only* to program around what seems
        # to be a bug in the cPickle module in some obscure case.
        # If we don't due this, then this doctest fails.
        #    sage: M = ModularSymbols(DirichletGroup(5).0)
        #    sage: loads(dumps(M)) == M

        self.__character = character.__copy__()

        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight - 2 + 1)
             for u, v in self.__P1.list()]
        self.__list = L
        ManinSymbolList.__init__(self, weight, L)

    def __repr__(self):
        """
        Standard function returning string representation.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: str(m) # indirect doctest
            'Manin Symbol List of weight 2 for Gamma1(4) with character [-1]'
        """
        return "Manin Symbol List of weight %s for Gamma1(%s) with character %s" % (
            self.weight(), self.level(), self.character()._repr_short_())

    def level(self):
        """
        Return the level of this :class:`ManinSymbolList`.

        OUTPUT: integer; the level of the symbols in this list

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,4).level()
            4
        """
        return self.__level

    def apply(self, j, m):
        """
        Apply the integer matrix `m=[a,b;c,d]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; the index of the symbol to act on

        - ``m`` -- list of integers `[a,b,c,d]` where `m = [a, b; c, d]` is the
          matrix to be applied

        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of ``self`` under the right action
        of the matrix `[a,b;c,d]`. Here the right action of
        `g = [a,b;c,d]` on a Manin symbol `[P(X,Y),(u,v)]` is by
        definition `[P(aX+bY,cX+dY),(u,v)*g]`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4)
            sage: m[6]
            [X*Y,(0,1)]
            sage: m.apply(4, [1,0,0,1])
            [(4, 1)]
            sage: m.apply(1, [-1,0,0,1])
            [(1, -1)]
        """
        a, b, c, d = m[0], m[1], m[2], m[3]
        i, u, v = self._symbol_list[j]
        P = apply_to_monomial(i, self._weight-2, a, b, c, d)
        m, s = self.index((0, u*a+v*c, u*b+v*d))
        if m == -1 or s == 0:
            return []
        r = len(self.__P1)
        return [(m + r*k, s*P[k]) for k in range(self._weight-2+1)
                if P[k] != 0]

    def apply_S(self, j):
        """
        Apply the matrix `S=[0,1;-1,0]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        ``(k, s)`` where `k` is the index of the symbol obtained by acting
        on the `j`-th symbol with `S`, and `s` is the parity of the
        `j`-th symbol.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_S(4)
            (2, -1)
            sage: [m.apply_S(i) for i in range(len(m))]
            [(1, 1), (0, -1), (4, 1), (5, -1), (2, -1), (3, 1)]
        """
        i, u, v = self._symbol_list[j]
        k, s = self.index((self._weight - 2 - i, v, -u))
        if i % 2 == 0:
            return k, s
        else:
            return k, -s

    def _apply_S_only_0pm1(self):
        """
        Return ``True`` if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            True
            sage: ManinSymbolList_character(DirichletGroup(13).0,2)._apply_S_only_0pm1()
            False
        """
        return self.__character.order() <= 2

    def apply_I(self, j):
        """
        Apply the matrix `I=[-1,0,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        ``(k, s)`` where `k` is the index of the symbol obtained by acting
        on the `j`-th symbol with `I`, and `s` is the parity of the
        `j`-th symbol.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_I(4)
            (2, -1)
            sage: [m.apply_I(i) for i in range(len(m))]
            [(0, 1), (1, -1), (4, -1), (3, -1), (2, -1), (5, 1)]
        """
        i, u, v = self._symbol_list[j]
        k, s = self.index((i, -u, v))
        if i % 2 == 0:
            return k, s
        else:
            return k, -s

    def apply_T(self, j):
        """
        Apply the matrix `T=[0,1,-1,-1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of ``self`` under the right action
        of the matrix `T`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_T(4)
            [(1, -1)]
            sage: [m.apply_T(i) for i in range(len(m))]
            [[(4, 1)], [(0, -1)], [(3, 1)], [(5, 1)], [(1, -1)], [(2, 1)]]
        """
        k = self._weight
        i, u, v = self._symbol_list[j]
        u, v, r = self.__P1.normalize_with_scalar(v,-u-v)
        r = self.__character(r)
        if (k-2) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        a = Integer(k-2-i)
        for j in range(k-2-i+1):
            m, r = self.index((j, u, v))
            z.append((m, s * r * a.binomial(j)))
            s *= -1
        return z

    def apply_TT(self, j):
        """
        Apply the matrix `TT=[-1,-1,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` -- integer; a symbol index

        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of ``self`` under the right action
        of the matrix `T^2`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_TT(4)
            [(0, 1)]
            sage: [m.apply_TT(i) for i in range(len(m))]
            [[(1, -1)], [(4, -1)], [(5, 1)], [(2, 1)], [(0, 1)], [(3, 1)]]
        """
        k = self._weight
        i, u, v = self._symbol_list[j]
        u, v, r = self.__P1.normalize_with_scalar(-u-v,u)
        r = self.__character(r)
        if (k-2-i) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        a = Integer(i)
        for j in range(i+1):
            m, r = self.index((k-2-i+j, u, v))
            z.append((m, s * r * a.binomial(j)))
            s *= -1
        return z

    def character(self):
        """
        Return the character of this :class:`ManinSymbolList_character` object.

        OUTPUT: the Dirichlet character of this Manin symbol list

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.character()
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1
        """
        return self.__character

    def index(self, x):
        """
        Return the index of a standard Manin symbol equivalent to
        ``x``, together with a scaling factor.

        INPUT:

        - ``x`` -- 3-tuple of integers defining an element of this
          list of Manin symbols, which need not be normalized

        OUTPUT:

        A pair ``(i, s)`` where ``i`` is the index of the Manin symbol
        equivalent to ``x`` and ``s`` is the scalar (an element of the
        base field).  If there is no Manin symbol equivalent to ``x``
        in the list, then ``(-1, 0)`` is returned.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4); m
            Manin Symbol List of weight 4 for Gamma1(4) with character [-1]
            sage: [m.index(s.tuple()) for s in m.manin_symbol_list()]
            [(0, 1),
            (1, 1),
            (2, 1),
            (3, 1),
            ...
            (16, 1),
            (17, 1)]
        """
        if x in self._index:
            return self._index[x], 1
        x, s = self.normalize(x)
        try:
            return self._index[x], s
        except KeyError:
            return -1, 0

    def normalize(self, x):
        """
        Return the normalization of the Manin Symbol ``x`` with
        respect to this list, together with the normalizing scalar.

        INPUT:

        - ``x`` -- 3-tuple of integers ``(i,u,v)``, defining an element of this
          list of Manin symbols, which need not be normalized

        OUTPUT:

        ``((i,u,v),s)``, where ``(i,u,v)`` is the normalized Manin symbol equivalent
        to ``x``, and ``s`` is the normalizing scalar.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4); m
            Manin Symbol List of weight 4 for Gamma1(4) with character [-1]
            sage: [m.normalize(s.tuple()) for s in m.manin_symbol_list()]
            [((0, 0, 1), 1),
             ((0, 1, 0), 1),
             ((0, 1, 1), 1),
             ...
             ((2, 1, 3), 1),
             ((2, 2, 1), 1)]
        """
        u,v,s = self.__P1.normalize_with_scalar(x[1],x[2])
        return (x[0],u,v), self.__character(s)


register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList', ManinSymbolList)
register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList_group', ManinSymbolList_group)
register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList_gamma0', ManinSymbolList_gamma0)
register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList_gamma1', ManinSymbolList_gamma1)
register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList_gamma_h', ManinSymbolList_gamma_h)
register_unpickle_override('sage.modular.modsym.manin_symbols',
                           'ManinSymbolList_character', ManinSymbolList_character)
