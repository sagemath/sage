# sage.doctest: needs sage.combinat sage.modules
r"""
Spin Crystals

These are the crystals associated with the three spin
representations: the spin representations of odd orthogonal groups
(or rather their double covers); and the `+` and `-` spin
representations of the even orthogonal groups.

We follow Kashiwara and Nakashima (Journal of Algebra 165, 1994) in
representing the elements of the spin crystal by sequences of signs
`\pm`.
"""
#TODO: Do we want the following two representations?
#
#Two other representations are available as attributes
#:meth:`Spin.internal_repn` and :meth:`Spin.signature` of the crystal element.
#
#- A numerical internal representation, an integer `n` such that if `n-1`
#  is written in binary and the `1`'s are replaced by ``-``, the `0`'s by
#  ``+``
#
#- The signature, which is a list in which ``+`` is replaced by `+1` and
#  ``-`` by `-1`.


#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                          Daniel Bump    <bump at match.stanford.edu>
#                     2019 Travis Scrimshaw <tcscrims at gmail.com>
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
#                  http://www.gnu.org/licenses/
#****************************************************************************

from cpython.object cimport Py_EQ, Py_NE, Py_LE, Py_GE, Py_LT, Py_GT
from cysignals.memory cimport sig_malloc, sig_free
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element, parent
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.tableau import Tableau
from sage.rings.integer_ring import ZZ
from sage.typeset.ascii_art import AsciiArt
from sage.typeset.unicode_art import UnicodeArt


#########################
# Type B spin
#########################

def CrystalOfSpins(ct):
    r"""
    Return the spin crystal of the given type `B`.

    This is a combinatorial model for the crystal with highest weight
    `Lambda_n` (the `n`-th fundamental weight). It has
    `2^n` elements, here called Spins. See also
    :func:`~sage.combinat.crystals.letters.CrystalOfLetters`,
    :func:`~sage.combinat.crystals.spins.CrystalOfSpinsPlus`,
    and :func:`~sage.combinat.crystals.spins.CrystalOfSpinsMinus`.

    INPUT:

    -  ``['B', n]`` -- A Cartan type `B_n`.

    EXAMPLES::

        sage: C = crystals.Spins(['B',3])
        sage: C.list()
        [+++, ++-, +-+, -++, +--, -+-, --+, ---]
        sage: C.cartan_type()
        ['B', 3]

    ::

        sage: [x.signature() for x in C]
        ['+++', '++-', '+-+', '-++', '+--', '-+-', '--+', '---']

    TESTS::

        sage: crystals.TensorProduct(C,C,generators=[[C.list()[0],C.list()[0]]]).cardinality()
        35
    """
    ct = CartanType(ct)
    if ct[0] == 'B':
        return GenericCrystalOfSpins(ct, Spin_crystal_type_B_element, "spins")
    else:
        raise NotImplementedError


#########################
# Type D spins
#########################

def CrystalOfSpinsPlus(ct):
    r"""
    Return the plus spin crystal of the given type D.

    This is the crystal with highest weight `Lambda_n` (the
    `n`-th fundamental weight).

    INPUT:

    -  ``['D', n]`` -- A Cartan type `D_n`.

    EXAMPLES::

        sage: D = crystals.SpinsPlus(['D',4])
        sage: D.list()
        [++++, ++--, +-+-, -++-, +--+, -+-+, --++, ----]

    ::

        sage: [x.signature() for x in D]
        ['++++', '++--', '+-+-', '-++-', '+--+', '-+-+', '--++', '----']

    TESTS::

        sage: TestSuite(D).run()
    """
    ct = CartanType(ct)
    if ct[0] == 'D':
        return GenericCrystalOfSpins(ct, Spin_crystal_type_D_element, "plus")
    else:
        raise NotImplementedError


def CrystalOfSpinsMinus(ct):
    r"""
    Return the minus spin crystal of the given type D.

    This is the crystal with highest weight `Lambda_{n-1}`
    (the `(n-1)`-st fundamental weight).

    INPUT:

    -  ``['D', n]`` -- A Cartan type `D_n`.

    EXAMPLES::

        sage: E = crystals.SpinsMinus(['D',4])
        sage: E.list()
        [+++-, ++-+, +-++, -+++, +---, -+--, --+-, ---+]
        sage: [x.signature() for x in E]
        ['+++-', '++-+', '+-++', '-+++', '+---', '-+--', '--+-', '---+']

    TESTS::

        sage: len(crystals.TensorProduct(E,E,generators=[[E[0],E[0]]]).list())
        35
        sage: D = crystals.SpinsPlus(['D',4])
        sage: len(crystals.TensorProduct(D,E,generators=[[D.list()[0],E.list()[0]]]).list())
        56
    """
    ct = CartanType(ct)
    if ct[0] == 'D':
        return GenericCrystalOfSpins(ct, Spin_crystal_type_D_element, "minus")
    else:
        raise NotImplementedError


class GenericCrystalOfSpins(UniqueRepresentation, Parent):
    """
    A generic crystal of spins.
    """
    def __init__(self, ct, element_class, case):
        """
        EXAMPLES::

            sage: E = crystals.SpinsMinus(['D',4])
            sage: TestSuite(E).run()
        """
        self._cartan_type = CartanType(ct)
        if case == "spins":
            self.rename("The crystal of spins for type %s" % ct)
        elif case == "plus":
            self.rename("The plus crystal of spins for type %s" % ct)
        else:
            self.rename("The minus crystal of spins for type %s" % ct)

        self.Element = element_class
        Parent.__init__(self, category=ClassicalCrystals())

        if case == "minus":
            generator = [1]*(ct[1]-1)
            generator.append(-1)
        else:
            generator = [1]*ct[1]
        self.module_generators = (self.element_class(self, tuple(generator)),)

    def _element_constructor_(self, value):
        """
        Construct an element of ``self`` from ``value``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: x = C((1,1,1)); x
            +++
            sage: y = C([1,1,1]); y
            +++
            sage: x == y
            True
        """
        return self.element_class(self, tuple(value))

    @lazy_attribute
    def _digraph_closure(self):
        """
        The transitive closure of the digraph associated to ``self``.

        EXAMPLES::

            sage: crystals.Spins(['B',4])._digraph_closure
            Transitive closure of : Digraph on 16 vertices
        """
        return self.digraph().transitive_closure()

    def lt_elements(self, x, y):
        r"""
        Return ``True`` if and only if there is a path from ``x`` to ``y``
        in the crystal graph.

        Because the crystal graph is classical, it is a directed acyclic
        graph which can be interpreted as a poset. This function implements
        the comparison function of this poset.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: x = C([1,1,1])
            sage: y = C([-1,-1,-1])
            sage: C.lt_elements(x, y)
            True
            sage: C.lt_elements(y, x)
            False
            sage: C.lt_elements(x, x)
            False
        """
        if parent(x) is not self or parent(y) is not self:
            raise ValueError("both elements must be in this crystal")
        return self._digraph_closure.has_edge(x, y)


cdef class Spin(Element):
    """
    A spin letter in the crystal of spins.

    EXAMPLES::

        sage: C = crystals.Spins(['B',3])
        sage: c = C([1,1,1])
        sage: c
        +++
        sage: c.parent()
        The crystal of spins for type ['B', 3]

        sage: D = crystals.Spins(['B',4])
        sage: a = C([1,1,1])
        sage: b = C([-1,-1,-1])
        sage: c = D([1,1,1,1])
        sage: a == a
        True
        sage: a == b
        False
        sage: b == c
        False
    """
    # cdef bint* self._value  # A + is a 0/False and a - is a 1/True

    def __init__(self, parent, tuple val):
        """
        Initialize ``self``.

        TESTS::

            sage: C = crystals.Spins(['B',3])
            sage: c = C([1,1,1])
            sage: TestSuite(c).run()
        """
        cdef int i
        self._n = parent.cartan_type().rank()
        self._value = <bint*>sig_malloc(self._n*sizeof(bint))
        for i in range(self._n):
            self._value[i] = (val[i] != 1)
        Element.__init__(self, parent)

    cdef Spin _new_c(self, bint* value):
        r"""
        Fast creation of a spin element.
        """
        cdef Spin ret = type(self).__new__(type(self))
        ret._parent = self._parent
        ret._n = self._n
        ret._value = value
        ret._hash = 0
        return ret

    def __dealloc__(self):
        """
        Deallocate ``self``.

        TESTS::

            sage: C = crystals.Spins(['B',3])
            sage: c = C([1,1,1])
            sage: del c
        """
        sig_free(self._value)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: C = crystals.Spins(['B',3])
            sage: len(set(C)) == len(set([hash(x) for x in C]))
            True
        """
        cdef int i
        if self._hash == 0:
            self._hash = hash(tuple([-1 if self._value[i] else 1 for i in range(self._n)]))
        return self._hash

    def __reduce__(self):
        r"""
        Used to pickle ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: a = C([1,-1,1])
            sage: a.__reduce__()
            (The crystal of spins for type ['B', 3], ((1, -1, 1),))
        """
        tup = tuple([-1 if self._value[i] else 1 for i in range(self._n)])
        return (self._parent, (tup,))

    cpdef _richcmp_(left, right, int op):
        """
        Return ``True`` if ``left`` compares with ``right`` based on ``op``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: x = C([1,1,1])
            sage: y = C([-1,-1,-1])
            sage: x < y
            True
            sage: x >= y
            False
            sage: x < x
            False
            sage: x <= x
            True
            sage: x != y
            True
            sage: x == y
            False
        """
        cdef Spin self, x
        cdef int i
        self = left
        x = right
        if op == Py_EQ:
            for i in range(self._n):
                if self._value[i] != x._value[i]:
                    return False
            return True
        if op == Py_NE:
            for i in range(self._n):
                if self._value[i] != x._value[i]:
                    return True
            return False
        if op == Py_LT:
            return self._parent._digraph_closure.has_edge(self, x)
        if op == Py_GT:
            return x._parent._digraph_closure.has_edge(x, self)
        if op == Py_LE:
            return self == x or self._parent._digraph_closure.has_edge(self, x)
        if op == Py_GE:
            return self == x or x._parent._digraph_closure.has_edge(x, self)
        return False

    @property
    def value(self):
        r"""
        Return ``self`` as a tuple with `+1` and `-1`.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: C([1,1,1]).value
            (1, 1, 1)
            sage: C([1,1,-1]).value
            (1, 1, -1)
        """
        cdef int i
        one = ZZ.one()
        return tuple([-one if self._value[i] else one for i in range(self._n)])

    def signature(self):
        """
        Return the signature of ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: C([1,1,1]).signature()
            '+++'
            sage: C([1,1,-1]).signature()
            '++-'
        """
        cdef int i
        cdef str sword = ""
        for i in range(self._n):
            sword += "+" if self._value[i] != 1 else "-"
        return sword

    _repr_ = signature

    def _repr_diagram(self):
        """
        Return a representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: print(b._repr_diagram())
            +
            +
            -
        """
        return '\n'.join(self.signature())

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: ascii_art(b)
            +
            +
            -
        """
        return AsciiArt(list(self.signature()))

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: unicode_art(b)
            +
            +
            -
        """
        return UnicodeArt(list(self.signature()))

    def pp(self):
        """
        Pretty print ``self`` as a column.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: b.pp()
            +
            +
            -
        """
        print(self._repr_diagram())

    def _latex_(self):
        r"""
        Gives the latex output of a spin column.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: print(b._latex_())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{1}c}\cline{1-1}
            \lr{-}\\\cline{1-1}
            \lr{+}\\\cline{1-1}
            \lr{+}\\\cline{1-1}
            \end{array}$}
            }
        """
        return Tableau([[i] for i in reversed(self.signature())])._latex_()

    def latex_standalone(self, **kwds):
        r"""
        Return Standalone picture illustrating ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: b = C([1,1,-1])
            sage: s = b.latex_standalone()
            sage: _ = s.pdf(view=False)           # long time (2s), optional - latex

        """
        from sage.misc.latex_standalone import Standalone
        if not 'standalone_config' in kwds:
            kwds['standalone_config'] = ["border=1pt"]
        return Standalone(self._latex_(), **kwds)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: [v.weight() for v in crystals.Spins(['B',3])]
            [(1/2, 1/2, 1/2), (1/2, 1/2, -1/2),
             (1/2, -1/2, 1/2), (-1/2, 1/2, 1/2),
             (1/2, -1/2, -1/2), (-1/2, 1/2, -1/2),
             (-1/2, -1/2, 1/2), (-1/2, -1/2, -1/2)]
        """
        WLR = self._parent.weight_lattice_realization()
        cdef int i
        mone = -WLR.base_ring().one()
        # The ambient space is indexed by 0,...,n-1
        return WLR._from_dict({i: mone**int(self._value[i]) / 2 for i in range(self._n)},
                              remove_zeros=False, coerce=False)

cdef class Spin_crystal_type_B_element(Spin):
    r"""
    Type B spin representation crystal element
    """
    cpdef Spin e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: [[C[m].e(i) for i in range(1,4)] for m in range(8)]
            [[None, None, None], [None, None, +++], [None, ++-, None], [+-+, None, None],
             [None, None, +-+], [+--, None, -++], [None, -+-, None], [None, None, --+]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        cdef int j
        cdef bint* ret
        if i == self._n:
            if self._value[i-1]:
                ret = <bint*>sig_malloc(self._n*sizeof(bint))
                for j in range(self._n):
                    ret[j] = self._value[j]
                ret[i-1] = False
                return self._new_c(ret)
            return None

        if self._value[i-1] and not self._value[i]:
            ret = <bint*>sig_malloc(self._n*sizeof(bint))
            for j in range(self._n):
                ret[j] = self._value[j]
            ret[i-1] = False
            ret[i] = True
            return self._new_c(ret)
        return None

    cpdef Spin f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: [[C[m].f(i) for i in range(1,4)] for m in range(8)]
            [[None, None, ++-], [None, +-+, None], [-++, None, +--], [None, None, -+-],
             [-+-, None, None], [None, --+, None], [None, None, ---], [None, None, None]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        cdef int j
        cdef bint* ret
        if i == self._n:
            if not self._value[i-1]:
                ret = <bint*>sig_malloc(self._n*sizeof(bint))
                for j in range(self._n):
                    ret[j] = self._value[j]
                ret[i-1] = True
                return self._new_c(ret)
            return None

        if self._value[i] and not self._value[i-1]:
            ret = <bint*>sig_malloc(self._n*sizeof(bint))
            for j in range(self._n):
                ret[j] = self._value[j]
            ret[i-1] = True
            ret[i] = False
            return self._new_c(ret)
        return None

    cpdef int epsilon(self, int i) noexcept:
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: [[C[m].epsilon(i) for i in range(1,4)] for m in range(8)]
            [[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0],
             [0, 0, 1], [1, 0, 1], [0, 1, 0], [0, 0, 1]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        if i == self._n:
            return self._value[i-1]
        return self._value[i-1] and not self._value[i]

    cpdef int phi(self, int i) noexcept:
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Spins(['B',3])
            sage: [[C[m].phi(i) for i in range(1,4)] for m in range(8)]
            [[0, 0, 1], [0, 1, 0], [1, 0, 1], [0, 0, 1],
             [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        if i == self._n:
            return not self._value[i-1]
        return self._value[i] and not self._value[i-1]

cdef class Spin_crystal_type_D_element(Spin):
    r"""
    Type D spin representation crystal element
    """
    cpdef Spin e(self, int i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: D = crystals.SpinsPlus(['D',4])
            sage: [[D.list()[m].e(i) for i in range(1,4)] for m in range(8)]
            [[None, None, None], [None, None, None], [None, ++--, None], [+-+-, None, None],
             [None, None, +-+-], [+--+, None, -++-], [None, -+-+, None], [None, None, None]]

        ::

            sage: E = crystals.SpinsMinus(['D',4])
            sage: [[E[m].e(i) for i in range(1,4)] for m in range(8)]
            [[None, None, None], [None, None, +++-], [None, ++-+, None], [+-++, None, None],
             [None, None, None], [+---, None, None], [None, -+--, None], [None, None, --+-]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        cdef int j
        cdef bint* ret
        if i == self._n:
            if self._value[i-1] and self._value[i-2]:
                ret = <bint*>sig_malloc(self._n*sizeof(bint))
                for j in range(self._n):
                    ret[j] = self._value[j]
                ret[i-1] = False
                ret[i-2] = False
                return self._new_c(ret)
            return None

        if self._value[i-1] and not self._value[i]:
            ret = <bint*>sig_malloc(self._n*sizeof(bint))
            for j in range(self._n):
                ret[j] = self._value[j]
            ret[i-1] = False
            ret[i] = True
            return self._new_c(ret)
        return None

    cpdef Spin f(self, int i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: D = crystals.SpinsPlus(['D',4])
            sage: [[D.list()[m].f(i) for i in range(1,4)] for m in range(8)]
            [[None, None, None], [None, +-+-, None], [-++-, None, +--+], [None, None, -+-+],
             [-+-+, None, None], [None, --++, None], [None, None, None], [None, None, None]]

        ::

            sage: E = crystals.SpinsMinus(['D',4])
            sage: [[E[m].f(i) for i in range(1,4)] for m in range(8)]
            [[None, None, ++-+], [None, +-++, None], [-+++, None, None], [None, None, None],
             [-+--, None, None], [None, --+-, None], [None, None, ---+], [None, None, None]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        cdef int j
        cdef bint* ret
        if i == self._n:
            if not self._value[i-1] and not self._value[i-2]:
                ret = <bint*>sig_malloc(self._n*sizeof(bint))
                for j in range(self._n):
                    ret[j] = self._value[j]
                ret[i-1] = True
                ret[i-2] = True
                return self._new_c(ret)
            return None

        if self._value[i] and not self._value[i-1]:
            ret = <bint*>sig_malloc(self._n*sizeof(bint))
            for j in range(self._n):
                ret[j] = self._value[j]
            ret[i-1] = True
            ret[i] = False
            return self._new_c(ret)
        return None

    cpdef int epsilon(self, int i) noexcept:
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.SpinsMinus(['D',4])
            sage: [[C[m].epsilon(i) for i in C.index_set()] for m in range(8)]
            [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0],
             [0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        if i == self._n:
            return self._value[i-1] and self._value[i-2]
        return self._value[i-1] and not self._value[i]

    cpdef int phi(self, int i) noexcept:
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.SpinsPlus(['D',4])
            sage: [[C[m].phi(i) for i in C.index_set()] for m in range(8)]
            [[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 0],
             [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]]
        """
        if i < 1 or i > self._n:
            raise ValueError("i is not in the index set")
        if i == self._n:
            return not self._value[i - 1] and not self._value[i - 2]
        return self._value[i] and not self._value[i - 1]
