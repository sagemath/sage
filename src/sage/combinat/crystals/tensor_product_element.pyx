# sage.doctest: needs sage.combinat sage.modules
"""
Tensor Products of Crystal Elements

AUTHORS:

- Anne Schilling, Nicolas Thiery (2007): Initial version
- Ben Salisbury, Travis Scrimshaw (2013): Refactored tensor products to handle
  non-regular crystals and created new subclass to take advantage of
  the regularity
- Travis Scrimshaw (2017): Cythonized element classes
- Franco Saliola (2017): Tensor products for crystal of super algebras
- Anne Schilling (2018): Tensor products for crystals of queer super algebras
"""
#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                     2017 Franco Saliola <saliola@gmail.com>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
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
#****************************************************************************

from sage.structure.parent cimport Parent

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.combinat.tableau import Tableau
from sage.rings.integer_ring import ZZ

##############################################################################
# Support classes
##############################################################################

cdef class ImmutableListWithParent(ClonableArray):
    r"""
    A class for lists having a parent.

    Specification: any subclass ``C`` should implement ``__init__`` which
    accepts the following form ``C(parent, list=list)``
    """
    def __init__(self, Parent parent, list):
        """
        Initialize ``self``.

        TESTS::

            sage: b = crystals.Tableaux(['A',2], shape=[2,1]).module_generators[0]
            sage: TestSuite(b).run()
        """
        ClonableArray.__init__(self, parent, list, check=False)

    def __setstate__(self, state):
        """
        For unpickling old pickles.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',2], shape=[2,1])
            sage: b = T.module_generators[0]
            sage: b.__setstate__([T, {'_list': list(b)}])
        """
        self._parent = state[0]
        self._list = state[1]['_list']
        self._is_immutable = True
        self._hash = 0

    cpdef _set_index(self, k, value):
        r"""
        Return a sibling of ``self`` obtained by setting the
        `k`-th entry of ``self`` to value.

        EXAMPLES::

            sage: b = crystals.Tableaux(['A',2], shape=[3]).module_generators[0]
            sage: list(b._set_index(0, 2))
            [2, 1, 1]
            sage: list(b._set_index(1, 4))
            [1, 4, 1]
        """
        cdef list l = list(self._list) # Make a (shallow) copy
        l[k] = value
        return type(self)(self._parent, list=l)

##############################################################################
# Primary classes
##############################################################################

cdef class TensorProductOfCrystalsElement(ImmutableListWithParent):
    r"""
    A class for elements of tensor products of crystals.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',3])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2))
            [1, 2]
        """
        if self._parent.options.convention == "Kashiwara":
            return repr(list(reversed(self._list)))
        return repr(self._list)

    def _latex_(self):
        r"""
        Return latex code for ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(["A",2])
            sage: D = crystals.Tableaux(["A",2], shape=[2])
            sage: E = crystals.TensorProduct(C,D)
            sage: latex(E.module_generators[0])
            1 \otimes {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{1}\\\cline{1-2}
            \end{array}$}
            }
        """
        from sage.misc.latex import latex
        if self._parent.options.convention == "Kashiwara":
            return r' \otimes '.join(latex(c) for c in reversed(self))
        return r' \otimes '.join(latex(c) for c in self)

    def latex_standalone(self, **kwds):
        r"""
        Return Standalone picture illustrating ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(["A",2])
            sage: D = crystals.Tableaux(["A",2], shape=[2])
            sage: E = crystals.TensorProduct(C,D)
            sage: x = E.module_generators[0]
            sage: s = x.latex_standalone()
            sage: _ = s.pdf(view=False)       # long time (2s), optional - latex

        """
        from sage.misc.latex_standalone import Standalone
        if not 'standalone_config' in kwds:
            kwds['standalone_config'] = ["border=1pt"]
        return Standalone(self._latex_(), **kwds)

    def _ascii_art_(self):
        """
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D',4,1],[[3,3],[2,1],[1,2]])
            sage: ascii_art(KT.module_generators[0])
              1  1  1
              2  2  2 #   1 #   1  1
              3  3  3     2
             -4 -4 -4
        """
        if self._parent.options.convention == "Kashiwara":
            lst = list(reversed(self))
        else:
            lst = self
        from sage.typeset.ascii_art import ascii_art, AsciiArt
        s = ascii_art(lst[0])
        s._baseline = s._h // 2
        ret = s
        for tableau in lst[1:]:
            s = ascii_art(tableau)
            s._baseline = s._h // 2
            ret += AsciiArt([" # "]) + s
        return ret

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D',4,1],[[3,3],[2,1],[1,2]])
            sage: unicode_art(KT.module_generators[0])
            ┌───┬───┬───┐
            │ 1 │ 1 │ 1 │
            ├───┼───┼───┤      ┌───┐
            │ 2 │ 2 │ 2 │      │ 1 │   ┌───┬───┐
            ├───┼───┼───┤    ⊗ ├───┤ ⊗ │ 1 │ 1 │
            │ 3 │ 3 │ 3 │      │ 2 │   └───┴───┘
            ├───┼───┼───┤      └───┘
            │ 4̄ │ 4̄ │ 4̄ │
            └───┴───┴───┘
        """
        if self._parent.options.convention == "Kashiwara":
            lst = list(reversed(self))
        else:
            lst = self
        from sage.typeset.unicode_art import unicode_art, UnicodeArt
        s = unicode_art(lst[0])
        s._baseline = s._h // 2
        ret = s
        for tableau in lst[1:]:
            s = unicode_art(tableau)
            s._baseline = s._h // 2
            ret += UnicodeArt([u" ⊗ "]) + s
        return ret

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[3,1])
            sage: D = crystals.Tableaux(['A',3], shape=[1])
            sage: E = crystals.Tableaux(['A',3], shape=[2,2,2])
            sage: T = crystals.TensorProduct(C,D,E)
            sage: print(T.module_generators[0]._repr_diagram())
              1  1  1 (X)   1 (X)   1  1
              2                     2  2
                                    3  3
        """
        pplist = []
        max_widths = []
        num_cols = len(self._list)
        for c in self:
            try:
                pplist.append(c._repr_diagram().split('\n'))
            except AttributeError:
                pplist.append(c._repr_().split('\n'))
            max_widths.append(max(map(len, pplist[-1])))
        num_rows = max(map(len, pplist))
        ret = ""
        for i in range(num_rows):
            if i > 0:
                ret += '\n'
            for j in range(num_cols):
                if j > 0:
                    if i == 0:
                        ret += ' (X) '
                    else:
                        ret += '     '
                if i < len(pplist[j]):
                    ret += pplist[j][i]
                    ret += ' '*(max_widths[j] - len(pplist[j][i]))
                else:
                    ret += ' '*max_widths[j]
        return ret

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',3], shape=[3,1])
            sage: D = crystals.Tableaux(['A',3], shape=[1])
            sage: E = crystals.Tableaux(['A',3], shape=[2,2,2])
            sage: T = crystals.TensorProduct(C,D,E)
            sage: T.module_generators[0].pp()
              1  1  1 (X)   1 (X)   1  1
              2                     2  2
                                    3  3
        """
        print(self._repr_diagram())

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("A3")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([2,1,3])
            sage: b2 = B.highest_weight_vector().f(1)
            sage: t = T(b2, b1)
            sage: t
            [[[1, 1, 1, 2], [2, 2], [3]], [[1, 1, 1, 1, 2], [2, 2, 4], [3]]]
            sage: t.weight()
            (-2, 1, 0, 1)

        ::

            sage: C = crystals.Letters(['A',3])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2)).weight()
            (1, 1, 0, 0)
            sage: T = crystals.Tableaux(['D',4],shape=[])
            sage: T.list()[0].weight()
            (0, 0, 0, 0)
        """
        WLR = self._parent.weight_lattice_realization()
        return WLR(sum(elt.weight() for elt in self))

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("G2")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f(2)
            sage: b2 = B.highest_weight_vector().f_string([2,2,1])
            sage: t = T(b2, b1)
            sage: [t.epsilon(i) for i in B.index_set()]
            [0, 3]
        """
        return max(self._sig(i, k) for k in range(1, len(self._list)+1))

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
            sage: B = crystals.GeneralizedYoungWalls(2,La[0]+La[1])
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([1,0])
            sage: b2 = B.highest_weight_vector().f_string([0,1])
            sage: t = T(b2, b1)
            sage: [t.phi(i) for i in B.index_set()]
            [1, 1, 4]

        TESTS:

        Check that :issue:`15462` is fixed::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: La = RootSystem(['A',2]).ambient_space().fundamental_weights()
            sage: T = crystals.TensorProduct(crystals.elementary.T(['A',2], La[1]+La[2]), B)
            sage: t = T.an_element()
            sage: t.phi(1)
            2
            sage: t.phi(2)
            2
        """
        P = self._list[-1].parent().weight_lattice_realization()
        h = P.simple_coroots()
        omega = P(self.weight()).scalar(h[i])
        return max(omega + self._sig(i, k) for k in range(1, len(self._list)+1))

    @cached_in_parent_method
    def _sig(self, i, k):
        r"""
        Return `a_i(k)` of ``self``.

        The value `a_i(k)` of a crystal `b = b_N \otimes \cdots \otimes b_1`
        is defined as:

        .. MATH::

            a_i(k) = \varepsilon_i(b_k) - \sum_{j=1}^{k-1} \langle h_i,
            \mathrm{wt}(b_j) \rangle

        where `\mathrm{wt}` is the :meth:`weight` of `b_j`.

        INPUT:

        - ``i`` -- an element of the index set

        - ``k`` -- the (1-based) index of the tensor factor of ``self``

        EXAMPLES::

            sage: B = crystals.infinity.GeneralizedYoungWalls(3)
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([0,3,1])
            sage: b2 = B.highest_weight_vector().f_string([3,2,1,0,2,3])
            sage: t = T(b1, b2)
            sage: [[t._sig(i,k) for k in range(1,len(t)+1)] for i in B.index_set()]
            [[0, -1], [0, 0], [0, 1], [1, 2]]

        TESTS:

        Check that :issue:`18469` is fixed::

            sage: E1 = crystals.elementary.B(['A',2], 1)
            sage: E2 = crystals.elementary.B(['A',2], 2)
            sage: T = crystals.TensorProduct(E1, E2)
            sage: x = T(E1.module_generators[0], E2.module_generators[0]); x
            [0, 0]
            sage: [[x._sig(i,k) for k in range(1,3)] for i in T.index_set()]
            [[-inf, 0], [0, -inf]]
            sage: x.f(1)
            [-1, 0]
            sage: x.e(1)
            [1, 0]
        """
        if k == 1:
            return self._list[-1].epsilon(i)
        ep = self._list[-k].epsilon(i)
        if ep == float("-inf"):
            return ep

        P = self._list[-1].parent().weight_lattice_realization()
        h = P.simple_coroots()
        wt = P.sum(P(self._list[-j].weight()) for j in range(1, k))
        return ep - wt.scalar(h[i])

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux("D4")
            sage: T = crystals.TensorProduct(B,B)
            sage: b1 = B.highest_weight_vector().f_string([1,4,3])
            sage: b2 = B.highest_weight_vector().f_string([2,2,3,1,4])
            sage: t = T(b2, b1)
            sage: t.e(1)
            [[[1, 1, 1, 1, 1], [2, 2, 3, -3], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
            sage: t.e(2)
            sage: t.e(3)
            [[[1, 1, 1, 1, 1, 2], [2, 2, 3, -4], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
            sage: t.e(4)
            [[[1, 1, 1, 1, 1, 2], [2, 2, 3, 4], [3]], [[1, 1, 1, 1, 2], [2, 2, 2], [3, -3]]]
        """
        N = len(self._list) + 1
        for k in range(1, N):
            if all(self._sig(i,k) > self._sig(i,j) for j in range(1, k)) and \
                   all(self._sig(i,k) >= self._sig(i,j) for j in range(k+1, N)):
                crystal = self._list[-k].e(i)
                if crystal is None:
                    return None
                return self._set_index(-k, crystal)
        return None

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: La = RootSystem(['A',3,1]).weight_lattice(extended=True).fundamental_weights()
            sage: B = crystals.GeneralizedYoungWalls(3,La[0])
            sage: T = crystals.TensorProduct(B,B,B)
            sage: b1 = B.highest_weight_vector().f_string([0,3])
            sage: b2 = B.highest_weight_vector().f_string([0])
            sage: b3 = B.highest_weight_vector()
            sage: t = T(b3, b2, b1)
            sage: t.f(0)
            [[[0]], [[0]], [[0, 3]]]
            sage: t.f(1)
            [[], [[0]], [[0, 3], [1]]]
            sage: t.f(2)
            [[], [[0]], [[0, 3, 2]]]
            sage: t.f(3)
            [[], [[0, 3]], [[0, 3]]]
        """
        N = len(self._list) + 1
        for k in range(1, N):
            if all(self._sig(i,k) >= self._sig(i,j) for j in range(1, k)) and \
                   all(self._sig(i,k) > self._sig(i,j) for j in range(k+1, N)):
                crystal = self._list[-k].f(i)
                if crystal is None:
                    return None
                return self._set_index(-k, crystal)
        return None


cdef class TensorProductOfRegularCrystalsElement(TensorProductOfCrystalsElement):
    """
    Element class for a tensor product of regular crystals.

    TESTS::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C, C)
        sage: elt = T(C(1), C(2))
        sage: from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        sage: isinstance(elt, TensorProductOfRegularCrystalsElement)
        True
    """
    def e(self, i):
        """
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(2)).e(1) == T(C(1),C(1))
            True
            sage: T(C(2),C(1)).e(1) is None
            True
            sage: T(C(2),C(2)).e(1) == T(C(1),C(2))
            True
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        k = self.position_of_first_unmatched_plus(i)
        if k is None:
            return None
        return self._set_index(k, self._list[k].e(i))

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).f(1)
            [1, 2]
            sage: T(C(1),C(2)).f(1)
            [2, 2]
            sage: T(C(2),C(1)).f(1) is None
            True
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        k = self.position_of_last_unmatched_minus(i)
        if k is None:
            return None
        return self._set_index(k, self._list[k].f(i))

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).phi(1)
            2
            sage: T(C(1),C(2)).phi(1)
            1
            sage: T(C(2),C(1)).phi(1)
            0
        """
        height = 0
        for elt in reversed(self._list):
            plus = elt.epsilon(i)
            minus = elt.phi(i)
            if height - plus < 0:
                height = minus
            else:
                height = height - plus + minus
        return height

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(1),C(1)).epsilon(1)
            0
            sage: T(C(1),C(2)).epsilon(1)
            1
            sage: T(C(2),C(1)).epsilon(1)
            0
        """
        height = 0
        for elt in self:
            minus = elt.phi(i)
            plus = elt.epsilon(i)
            if height - minus < 0:
                height = plus
            else:
                height = height - minus + plus
        return height

    cpdef position_of_last_unmatched_minus(self, i):
        """
        Return the position of the last unmatched `-` or ``None`` if
        there is no unmatched `-`.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).position_of_last_unmatched_minus(1)
            sage: T(C(1),C(2)).position_of_last_unmatched_minus(1)
            0
        """
        unmatched_minus = None
        height = 0
        cdef int j
        for j,elt in enumerate(self):
            plus = elt.epsilon(i)
            minus = elt.phi(i)
            if height - minus < 0:
                unmatched_minus = j
                height = plus
            else:
                height = height - minus + plus
        return unmatched_minus

    cpdef position_of_first_unmatched_plus(self, i):
        """
        Return the position of the first unmatched `+` or ``None`` if
        there is no unmatched `+`.

        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).position_of_first_unmatched_plus(1)
            sage: T(C(1),C(2)).position_of_first_unmatched_plus(1)
            1
        """
        unmatched_plus = None
        height = 0
        cdef int N = len(self._list) - 1
        cdef int j
        for j, elt in enumerate(reversed(self._list)):
            plus = elt.epsilon(i)
            minus = elt.phi(i)
            if height - plus < 0:
                unmatched_plus = N - j
                height = minus
            else:
                height = height - plus + minus
        return unmatched_plus

    # Legacy function
    def positions_of_unmatched_minus(self, i, dual=False, reverse=False):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_minus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_minus(1)
            [0]
        """
        cdef list unmatched_plus = []
        cdef int j
        height = 0
        if reverse:
            self = type(self)(self._parent, list(reversed(self._list)))
        if not dual:
            for j,elt in enumerate(self):
                minus = elt.phi(i)
                plus = elt.epsilon(i)
                if height-minus < 0:
                    unmatched_plus.append(j)
                    height = plus
                else:
                    height = height - minus + plus
        else:
            for j,elt in enumerate(self):
                plus = elt.epsilon(i)
                minus = elt.phi(i)
                if height-plus < 0:
                    unmatched_plus.append(j)
                    height = minus
                else:
                    height = height - plus + minus
        return unmatched_plus

    # Legacy function
    def positions_of_unmatched_plus(self, i):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',5])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_plus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_plus(1)
            [1]
        """
        cdef list L = self.positions_of_unmatched_minus(i, dual=True, reverse=True)
        L.reverse()
        cdef int N = len(self._list) - 1
        return [N - val for val in L]

cdef class CrystalOfTableauxElement(TensorProductOfRegularCrystalsElement):
    """
    Element in a crystal of tableaux.
    """
    def __init__(self, parent, *args, **options):
        """
        There are several ways to input tableaux, by rows, by columns,
        by columns, as the list of column elements, or as a sequence
        of numbers in column reading.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(columns=[[3,1],[4,2]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(list=[3,1,4,2])
            sage: t
            [[1, 2], [3, 4]]

            sage: t = T(3,1,4,2)
            sage: t
            [[1, 2], [3, 4]]

        Currently inputting the empty tableau as an empty sequence is
        broken due to a bug in the generic __call__ method (see :issue:`8648`).

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape=[])
            sage: t = T()
            sage: list(t)
            [0]

        TESTS:

        Integer types that are not a Sage ``Integer`` (such as a Python ``int``
        and typically arise from compiled code) were not converted into a
        letter. This caused certain functions to fail. This is fixed in
        :issue:`13204`::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(list=[int(3),1,4,2])
            sage: type(t[0])
            <class 'sage.combinat.crystals.letters.Crystal_of_letters_type_A_element'>
            sage: t = T(list=[3,int(1),4,2])
            sage: type(t[1])
            <class 'sage.combinat.crystals.letters.Crystal_of_letters_type_A_element'>
            sage: C = crystals.KirillovReshetikhin(['A',int(3),1], 1,1)
            sage: C[0].e(0)
            [[4]]
        """
        if len(args) == 1:
            if isinstance(args[0], Tableau):
                options['rows'] = args[0]
        if 'list' in options:
            the_list = options['list']
        elif 'rows' in options:
            rows = options['rows']
#            the_list=Tableau(rows).to_word_by_column()
            rows = Tableau(rows).conjugate()
            the_list = []
            for col in rows:
                the_list += reversed(col)
        elif 'columns' in options:
            columns = options['columns']
            the_list = []
            for col in columns:
                the_list += col
        else:
            the_list = [i for i in args]
        TensorProductOfRegularCrystalsElement.__init__(self, parent,
                                                       [parent.letters(i)
                                                        for i in the_list])

    def _repr_(self):
        """
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t._repr_()
            '[[1, 2], [3, 4]]'
        """
        return repr(self.to_tableau())

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A', 4], shape=[4,2,1])
            sage: elt = C(rows=[[1,1,1,2], [2,3], [4]])
            sage: print(elt._repr_diagram())
              1  1  1  2
              2  3
              4
        """
        return self.to_tableau()._repr_diagram()

    def pp(self):
        """
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t.pp()
            1  2
            3  4
        """
        return self.to_tableau().pp()

    def _ascii_art_(self):
        """
        Return an ascii art version of ``self``.

        EXAMPLES:

        We check that :issue:`16486` is fixed::

            sage: T = crystals.Tableaux(['B',6], shape=[1]*5)
            sage: ascii_art(T.module_generators[0])
              1
              2
              3
              4
              5
            sage: T = crystals.Tableaux(['D',4], shape=[2,1])
            sage: t = T.module_generators[0].f_string([1,2,3,4,2,2,3,4])
            sage: ascii_art(t)
              1 -2
             -3
        """
        return self.to_tableau()._ascii_art_()

    def _unicode_art_(self):
        """
        Return a unicode art version of ``self``.

        EXAMPLES::

            sage: T = crystals.Tableaux(['B',4], shape=[1]*3)
            sage: unicode_art(T.module_generators[0])
            ┌───┐
            │ 1 │
            ├───┤
            │ 2 │
            ├───┤
            │ 3 │
            └───┘
            sage: T = crystals.Tableaux(['D',4], shape=[2,1])
            sage: t = T.module_generators[0].f_string([1,2,3,4,2,2,3,4])
            sage: unicode_art(t)
            ┌───┬───┐
            │ 1 │ 2̄ │
            ├───┼───┘
            │ 3̄ │
            └───┘
        """
        if not self._list:
            return Tableau([])._unicode_art_()
        cdef list lst = self._list
        cdef list tab = [ [lst[0]] ]
        cdef int i
        for i in range(1,len(self)):
            if lst[i-1] < lst[i] or (lst[i-1].value != 0 and lst[i-1] == lst[i]):
                tab.append([lst[i]])
            else:
                tab[len(tab)-1].append(lst[i])
        for x in tab:
            x.reverse()
        return Tableau(tab).conjugate()._unicode_art_()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [4,2])
            sage: t = T(rows=[[1,1,2,3],[2,3]])
            sage: latex(t) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{1}&\lr{2}&\lr{3}\\\cline{1-4}
            \lr{2}&\lr{3}\\\cline{1-2}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        # Modified version of to_tableau() to have the entries be letters
        #   rather than their values
        if not self._list:
            return "{\\emptyset}"

        tab = [ [self[0]] ]
        for i in range(1,len(self)):
            if self[i-1] < self[i] or (self[i-1].value != 0 and self[i-1] == self[i]):
                tab.append([self[i]])
            else:
                l = len(tab)-1
                tab[l].append(self[i])
        for x in tab:
            x.reverse()
        T = Tableau(tab).conjugate()
        return tex_from_array([[letter._latex_() for letter in row] for row in T])

    def latex_standalone(self, **kwds):
        r"""
        Return Standalone picture illustrating ``self``.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [4,2])
            sage: t = T(rows=[[1,1,2,3],[2,3]])
            sage: s = t.latex_standalone()
            sage: _ = s.pdf(view=False)        # long time (2s), optional - latex

        """
        from sage.misc.latex_standalone import Standalone
        if not 'standalone_config' in kwds:
            kwds['standalone_config'] = ["border=1pt"]
        return Standalone(self._latex_(), **kwds)

    @cached_method
    def to_tableau(self):
        """
        Return the :class:`Tableau` object corresponding to ``self``.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]]).to_tableau(); t
            [[1, 2], [3, 4]]
            sage: type(t)
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
            sage: type(t[0][0])
            <class 'int'>
            sage: T = crystals.Tableaux(['D',3], shape = [1,1])
            sage: t=T(rows=[[-3],[3]]).to_tableau(); t
            [[-3], [3]]
            sage: t=T(rows=[[3],[-3]]).to_tableau(); t
            [[3], [-3]]
            sage: T = crystals.Tableaux(['B',2], shape = [1,1])
            sage: t = T(rows=[[0],[0]]).to_tableau(); t
            [[0], [0]]
        """
        if not self._list:
            return Tableau([])
        cdef list lst = self._list
        cdef list tab = [ [lst[0].value] ]
        cdef int i
        for i in range(1,len(self)):
            if lst[i-1] < lst[i] or (lst[i-1].value != 0 and lst[i-1] == lst[i]):
                tab.append([lst[i].value])
            else:
                tab[len(tab)-1].append(lst[i].value)
        for x in tab:
            x.reverse()
        return Tableau(tab).conjugate()

    def shape(self):
        r"""
        Return the shape of the tableau corresponding to ``self``.

        OUTPUT: an instance of :class:`Partition`

        .. SEEALSO::

            :meth:`to_tableau`

        EXAMPLES::

            sage: C = crystals.Tableaux(["A", 2], shape=[2,1])
            sage: x = C.an_element()
            sage: x.to_tableau().shape()
            [2, 1]
            sage: x.shape()
            [2, 1]
        """
        return self.to_tableau().shape()

    def promotion(self):
        """
        Return the result of applying promotion on ``self``.

        Promotion for type A crystals of tableaux of rectangular shape.
        This method only makes sense in type A with rectangular shapes.

        EXAMPLES::

            sage: C = crystals.Tableaux(["A",3], shape = [3,3,3])
            sage: t = C(Tableau([[1,1,1],[2,2,3],[3,4,4]]))
            sage: t
            [[1, 1, 1], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion()
            [[1, 1, 2], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion().parent()
            The crystal of tableaux of type ['A', 3] and shape(s) [[3, 3, 3]]
        """
        crystal = self._parent
        cartan_type = crystal.cartan_type()
        assert cartan_type.type() == 'A'
        return crystal(self.to_tableau().promotion(cartan_type.rank()))

    def promotion_inverse(self):
        """
        Return the result of applying inverse promotion on ``self``.

        Inverse promotion for type A crystals of tableaux of rectangular shape.
        This method only makes sense in type A with rectangular shapes.

        EXAMPLES::

            sage: C = crystals.Tableaux(["A",3], shape = [3,3,3])
            sage: t = C(Tableau([[1,1,1],[2,2,3],[3,4,4]]))
            sage: t
            [[1, 1, 1], [2, 2, 3], [3, 4, 4]]
            sage: t.promotion_inverse()
            [[1, 1, 2], [2, 3, 3], [4, 4, 4]]
            sage: t.promotion_inverse().parent()
            The crystal of tableaux of type ['A', 3] and shape(s) [[3, 3, 3]]
        """
        crystal = self._parent
        cartan_type = crystal.cartan_type()
        assert cartan_type.type() == 'A'
        return crystal(self.to_tableau().promotion_inverse(cartan_type.rank()))

cdef class InfinityCrystalOfTableauxElement(CrystalOfTableauxElement):
    def e(self,i):
        r"""
        Return the action of `\widetilde{e}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['B',3])
            sage: b = B(rows=[[1,1,1,1,1,1,1,2,0,-3,-1,-1,-1,-1],[2,2,2,2,-2,-2],[3,-3,-3]])
            sage: b.e(3).pp()
            1  1  1  1  1  1  1  2  0 -3 -1 -1 -1 -1
            2  2  2  2 -2 -2
            3  0 -3
            sage: b.e(1).pp()
            1  1  1  1  1  1  1  0 -3 -1 -1 -1 -1
            2  2  2  2 -2 -2
            3 -3 -3
        """
        if i not in self.index_set():
            raise ValueError('i is not in the index set')
        k = self.position_of_first_unmatched_plus(i)
        if k is None:
            return None
        cdef InfinityCrystalOfTableauxElement ret
        ret = <InfinityCrystalOfTableauxElement>(self._set_index(k, self._list[k].e(i)))
        if k+i > len(self._list):
            return ret
        for j in reversed(range(1, i+1)):
            if ret._list[k+i-j].value != j:
                return ret
        # We've found a column, so we need to remove it
        for j in range(i):
            ret._list.pop(k)
        return ret

    def f(self, i):
        r"""
        Return the action of `\widetilde{f}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['C',4])
            sage: b = B.highest_weight_vector()
            sage: b.f(1).pp()
            1  1  1  1  2
            2  2  2
            3  3
            4
            sage: b.f(3).pp()
            1  1  1  1  1
            2  2  2  2
            3  3  4
            4
            sage: b.f(3).f(4).pp()
            1  1  1  1  1
            2  2  2  2
            3  3 -4
            4
        """
        if i not in self.index_set():
            raise ValueError('i is not in the index set')
        k = self.position_of_last_unmatched_minus(i)
        if k is None:
            return None
        cdef InfinityCrystalOfTableauxElement ret
        ret = <InfinityCrystalOfTableauxElement>(self._set_index(k, self._list[k].f(i)))
        if k+i > len(self._list):
            return ret
        for j in reversed(range(1,i+1)):
            if self._list[k+i-j].value != j:
                return ret
        # We've found a full column, so we'll need to add a new column
        for j in range(i):
            ret._list.insert(k, self._parent.letters(j+1))
        return ret

cdef class InfinityCrystalOfTableauxElementTypeD(InfinityCrystalOfTableauxElement):
    def e(self, i):
        r"""
        Return the action of `\widetilde{e}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['D',4])
            sage: b = B.highest_weight_vector().f_string([1,4,3,1,2]); b.pp()
            1  1  1  1  2  3
            2  2  2
            3 -3
            sage: b.e(2).pp()
            1  1  1  1  2  2
            2  2  2
            3 -3
        """
        if i not in self.index_set():
            raise ValueError('i is not in the index set')
        k = self.position_of_first_unmatched_plus(i)
        if k is None:
            return None
        cdef InfinityCrystalOfTableauxElementTypeD ret
        ret = <InfinityCrystalOfTableauxElementTypeD>(self._set_index(k, self._list[k].e(i)))
        if i == self.cartan_type().rank():
            i -= 1
        if k+i > len(self._list):
            return ret
        for j in reversed(range(1, i+1)):
            if ret._list[k+i-j].value != j:
                return ret
        # We've found a column, so we need to remove it
        for j in range(i):
            ret._list.pop(k)
        return ret

    def f(self, i):
        r"""
        Return the action of `\widetilde{f}_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['D',5])
            sage: b = B.highest_weight_vector().f_string([1,4,3,1,5]); b.pp()
            1  1  1  1  1  1  2  2
            2  2  2  2  2
            3  3  3 -5
            4  5
            sage: b.f(1).pp()
            1  1  1  1  1  1  2  2  2
            2  2  2  2  2
            3  3  3 -5
            4  5
            sage: b.f(5).pp()
            1  1  1  1  1  1  2  2
            2  2  2  2  2
            3  3  3 -5
            4 -4
        """
        cdef InfinityCrystalOfTableauxElementTypeD ret
        ret = <InfinityCrystalOfTableauxElementTypeD>(InfinityCrystalOfTableauxElement.f(self, i))
        if ret._list[0].value == -self._parent.cartan_type().rank():
            # Exceptional case for f_n where we need to add a new column
            for j in range(i-1):
                ret._list.insert(0, self._parent.letters(j+1))
        return ret

#####################################################################
#  BKK crystal elements

cdef class TensorProductOfSuperCrystalsElement(TensorProductOfRegularCrystalsElement):
    r"""
    Element class for a tensor product of crystals for Lie superalgebras.

    This implements the tensor product rule for crystals of
    Lie superalgebras of [BKK2000]_.

    TESTS::

        sage: C = crystals.Letters(['A', [2, 1]])
        sage: T = tensor([C,C])
        sage: T
        Full tensor product of the crystals [The crystal of letters for type ['A', [2, 1]], The crystal of letters for type ['A', [2, 1]]]
        sage: T.cardinality()
        25
        sage: t = T.an_element(); t
        [-3, -3]
        sage: t.weight()
        (2, 0, 0, 0, 0)
    """

    def e(self, i):
        r"""
        Return `e_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A', [2, 1]])
            sage: T = tensor([C,C])
            sage: t = T(C(1),C(1))
            sage: t.e(0)
            [-1, 1]
        """
        if i > 0:
            return TensorProductOfRegularCrystalsElement.e(self, i)
        if i < 0:
            x = type(self)(self._parent, reversed(self))
            k = x.position_of_first_unmatched_plus(i)
            if k is None:
                return None
            k = len(self._list) - k - 1
            return self._set_index(k, self._list[k].e(i))
        # Otherwise i == 0
        for k,elt in enumerate(self._list):
            if elt.f(i) is not None:
                return None
            x = elt.e(i)
            if x is not None:
                return self._set_index(k, x)
        return None

    def f(self, i):
        r"""
        Return `f_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A', [2, 1]])
            sage: T = tensor([C,C])
            sage: t = T(C(1),C(1))
            sage: t.f(0)
            sage: t.f(1)
            [1, 2]
        """
        if i > 0:
            return TensorProductOfRegularCrystalsElement.f(self, i)
        if i < 0:
            x = type(self)(self._parent, reversed(self))
            k = x.position_of_last_unmatched_minus(i)
            if k is None:
                return None
            k = len(self._list) - k - 1
            return self._set_index(k, self._list[k].f(i))
        # Otherwise i == 0
        for k,elt in enumerate(self._list):
            if elt.e(i) is not None:
                return None
            x = elt.f(i)
            if x is not None:
                return self._set_index(k, x)
        return None

    # Override epsilon/phi (for now)
    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A', [2, 1]])
            sage: T = tensor([C,C])
            sage: t = T(C(1),C(1))
            sage: t.epsilon(0)
            1
        """
        string_length = 0
        x = self
        while True:
            x = x.e(i)
            if x is None:
                return string_length
            else:
                string_length += 1

    def phi(self, i):
        r"""
        Return `\varphi_i` on ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A', [2, 1]])
            sage: T = tensor([C,C])
            sage: t = T(C(1),C(1))
            sage: t.phi(0)
            0
        """
        string_length = 0
        x = self
        while True:
            x = x.f(i)
            if x is None:
                return string_length
            else:
                string_length += 1

cdef class CrystalOfBKKTableauxElement(TensorProductOfSuperCrystalsElement):
    """
    Element class for the crystal of tableaux for Lie superalgebras
    of [BKK2000]_.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: C.an_element()
            [[-2], [-1]]
        """
        return repr(self.to_tableau())

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: c._repr_diagram()
            ' -2\n -1'
        """
        return self.to_tableau()._repr_diagram()

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: c.pp()
            -2
            -1
        """
        return self.to_tableau().pp()

    def _ascii_art_(self):
        """
        Return an ascii art version of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: ascii_art(c)
            -2
            -1
        """
        return self.to_tableau()._ascii_art_()

    def _unicode_art_(self):
        """
        Return a unicode art version of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: unicode_art(c)
            ┌───┐
            │ 2̄ │
            ├───┤
            │ 1̄ │
            └───┘
        """
        return self.to_tableau()._unicode_art_()

    def _latex_(self):
        r"""
        Return the latex code of ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: latex(c)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{1}c}\cline{1-1}
            \lr{-2}\\\cline{1-1}
            \lr{-1}\\\cline{1-1}
            \end{array}$}
            }
        """
        return self.to_tableau()._latex_()

    def latex_standalone(self, **kwds):
        r"""
        Return Standalone picture illustrating ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: s = c.latex_standalone()
            sage: _ = s.pdf(view=False)        # long time (2s), optional - latex

        """
        from sage.misc.latex_standalone import Standalone
        if not 'standalone_config' in kwds:
            kwds['standalone_config'] = ["border=1pt"]
        return Standalone(self._latex_(), **kwds)

    @cached_method
    def to_tableau(self):
        """
        Return the :class:`Tableau` object corresponding to ``self``.

        EXAMPLES::

            sage: C = crystals.Tableaux(['A',[1,2]], shape=[1,1])
            sage: c = C.an_element()
            sage: c.to_tableau()
            [[-2], [-1]]
            sage: type(c.to_tableau())
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
            sage: type(c)
            <class 'sage.combinat.crystals.bkk_crystals.CrystalOfBKKTableaux_with_category.element_class'>
        """
        sh = self._parent._shape.conjugate()
        tab = [[None]*row for row in sh]
        cur = 0
        lst = list(reversed(self._list))
        for r,row_len in enumerate(sh):
            for c in reversed(range(row_len)):
                tab[r][c] = lst[cur]
                cur += 1
        return Tableau(tab).conjugate()

#####################################################################
#  Queer crystal elements

cdef class TensorProductOfQueerSuperCrystalsElement(TensorProductOfRegularCrystalsElement):
    r"""
    Element class for a tensor product of crystals for queer Lie superalgebras.

    This implements the tensor product rule for crystals of Grantcharov et al.
    [GJK+2014]_. Given crystals `B_1` and `B_2` of type `\mathfrak{q}_{n+1}`,
    we define the tensor product `b_1 \otimes b_2 \in B_1 \otimes B_2`,
    where `b_1 \in B_1` and `b_2 \in B_2`, as the following:

    In addition to the tensor product rule for type `A_n`, the tensor product
    rule for `e_{-1}` and `f_{-1}` on `b_1\otimes b_2` are given by

    .. MATH::

        \begin{aligned}
        e_{-1}(b_1\otimes b_2) &=
        \begin{cases}
        b_1 \otimes e_{-1}b_2 &
         \text{if } \operatorname{wt}(b_1)_1 = \operatorname{wt}(b_1)_2 = 0,\\
        e_{-1}b_1 \otimes b_2 & \text{otherwise},
        \end{cases}
        \\
        f_{-1}(b_1\otimes b_2) &=
        \begin{cases}
        b_1 \otimes f_{-1}b_2 &
         \text{if } \operatorname{wt}(b_1)_1 = \operatorname{wt}(b_1)_2 = 0,\\
        f_{-1}b_1 \otimes b_2 & \text{otherwise}.
        \end{cases}
        \end{aligned}

    For `1 < i \leq n`, the operators `e_{-i}` and `f_{-i}` are defined as

    .. MATH::

        e_{-i} = s_{w^{-1}_i} e_{-1} s_{w_i}, \quad
        f_{-i} = s_{w^{-1}_i} f_{-1} s_{w_i}.

    Here, `w_i = s_2 \cdots s_i s_1 \cdots s_{i-1}` and `s_i` is the reflection
    along the `i`-string in the crystal. Moreover, for `1<i\leq n`, we define
    the operators `e_{-i'}` and `f_{-i'}` as

    .. MATH::

        e_{-i'} = s_{w_0} f_{-(n+1-i)} s_{w_0}, \quad
        f_{-i'} = s_{w_0} e_{-(n+1-i)} s_{w_0},

    where `w_0` is the longest element in the symmetric group `S_{n+1}`
    generated by `s_1,\ldots,s_n`. In this implementation, we use the integers
    `-2n, \ldots, -(n+1)` to respectively denote the indices `-n', \ldots, -1'`.

    TESTS::

        sage: Q = crystals.Letters(['Q', 3])
        sage: T = tensor([Q,Q]); T
        Full tensor product of the crystals
         [The queer crystal of letters for q(3),
          The queer crystal of letters for q(3)]
        sage: T.cardinality()
        9
        sage: t = T.an_element(); t
        [1, 1]
        sage: t.weight()
        (2, 0, 0)
    """
    def e(self, i):
        r"""
        Return `e_i` on ``self``.

        EXAMPLES::

            sage: Q = crystals.Letters(['Q', 3])
            sage: T = tensor([Q,Q])
            sage: t = T(Q(1),Q(1))
            sage: t.e(-1)
            sage: t = T(Q(2),Q(1))
            sage: t.e(-1)
            [1, 1]

            sage: T = tensor([Q,Q,Q,Q])
            sage: t = T(Q(1),Q(3),Q(2),Q(1))
            sage: t.e(-2)
            [2, 2, 1, 1]
        """
        if i > 0:
            from sage.categories.regular_supercrystals import RegularSuperCrystals
            if self._parent in RegularSuperCrystals():
                return TensorProductOfRegularCrystalsElement.e(self, i)
            else:
                return TensorProductOfCrystalsElement.e(self, i)
        cdef tuple w
        cdef int k, a, l
        l = len(self._list)
        if i == -1:
            for k in range(l):
                b = self._list[k].e(i)
                if b is not None:
                    return self._set_index(k, b)
                if self._list[k].f(i) is not None:
                    # There are no (-1)-string of length > 1
                    return None
            return None
        n = self._parent.cartan_type().n
        if i < -1 and i >= -n:
            j = -i
            b = self
            for a in range(j-1, 0, -1):
                b = b.s(a)
            for a in range(j, 1, -1):
                b = b.s(a)
            b = b.e(-1)
            if b is None:
                return None
            for a in range(2, j+1):
                b = b.s(a)
            for a in range(1, j):
                b = b.s(a)
            return b
        if i < -n:
            j = -(i+n)
            w = <tuple> (self._parent._long_element())
            b = self
            for a in w:
                b = b.s(a)
            b = b.f(-(n+1-j))
            if b is None:
                return None
            for a in w:
                b = b.s(a)
            return b
        return None

    def f(self, i):
        r"""
        Return `f_i` on ``self``.

        EXAMPLES::

            sage: Q = crystals.Letters(['Q', 3])
            sage: T = tensor([Q, Q])
            sage: t = T(Q(1), Q(1))
            sage: t.f(-1)
            [2, 1]
        """
        if i > 0:
            from sage.categories.regular_supercrystals import RegularSuperCrystals
            if self._parent in RegularSuperCrystals():
                return TensorProductOfRegularCrystalsElement.f(self, i)
            else:
                return TensorProductOfCrystalsElement.f(self, i)
        cdef tuple w
        cdef int k, a, l
        l = len(self._list)
        if i == -1:
            for k in range(l):
                b = self._list[k].f(i)
                if b is not None:
                    return self._set_index(k, b)
                if self._list[k].e(i) is not None:
                    # There are no (-1)-string of length > 1
                    return None
            return None
        n = self._parent.cartan_type().n
        if i < -1 and i >= -n:
            j = -i
            b = self
            for a in range(j-1, 0, -1):
                b = b.s(a)
            for a in range(j, 1, -1):
                b = b.s(a)
            b = b.f(-1)
            if b is None:
                return None
            for a in range(2, j+1):
                b = b.s(a)
            for a in range(1, j):
                b = b.s(a)
            return b
        if i < -n:
            j = -(i+n)
            w = <tuple> (self._parent._long_element())
            b = self
            for a in w:
                b = b.s(a)
            b = b.e(-(n+1-j))
            if b is None:
                return None
            for a in w:
                b = b.s(a)
            return b
        return None

    # Override epsilon/phi (for now)
    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` on ``self``.

        EXAMPLES::

            sage: Q = crystals.Letters(['Q', 3])
            sage: T = tensor([Q, Q, Q, Q])
            sage: t = T(Q(1), Q(3), Q(2), Q(1))
            sage: t.epsilon(-2)
            1
        """
        string_length = -1
        x = self
        while x is not None:
            string_length += 1
            x = x.e(i)
        return string_length

    def phi(self, i):
        r"""
        Return `\varphi_i` on ``self``.

        EXAMPLES::

            sage: Q = crystals.Letters(['Q', 3])
            sage: T = tensor([Q, Q, Q, Q])
            sage: t = T(Q(1), Q(3), Q(2), Q(1))
            sage: t.phi(-2)
            0
            sage: t.phi(-1)
            1
        """
        string_length = -1
        x = self
        while x is not None:
            string_length += 1
            x = x.f(i)
        return string_length


cdef class InfinityQueerCrystalOfTableauxElement(TensorProductOfQueerSuperCrystalsElement):
    def __init__(self, parent, list, row_lengths=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,2,1],[3,3,3],[2,2],[1]])
            sage: t
            [[4, 4, 4, 4, 2, 1], [3, 3, 3], [2, 2], [1]]
            sage: TestSuite(t).run()
        """
        if row_lengths is None:
            row_lengths = []
        if not row_lengths and list and not isinstance(list[0], parent.letters.element_class):
            ret = []
            L = parent.letters
            row_lengths = []
            for row in list:
                ret.extend(L(val) for val in reversed(row))
                row_lengths.append(len(row))
            list = ret
        self._row_lengths = row_lengths
        super().__init__(parent, list)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B.an_element()
            sage: t
            [[4, 4, 4, 4], [3, 3, 3], [2, 2], [1]]
        """
        return repr([list(reversed(row)) for row in self.rows()])

    def _ascii_art_(self):
        """
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,2,1],[3,3,3],[2,2],[1]])
            sage: ascii_art(t)
              4  4  4  4  2  1
                 3  3  3
                    2  2
                       1
        """
        from sage.typeset.ascii_art import AsciiArt
        ret = [" "*(3*i) + "".join("%3s" % str(x) for x in reversed(row))
               for i, row in enumerate(self.rows())]
        return AsciiArt(ret)

    def _latex_(self):
        r"""
        Return latex code for ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: latex(t)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{7}c}\cline{1-7}
            \lr{4}&\lr{4}&\lr{4}&\lr{4}&\lr{4}&\lr{2}&\lr{1}\\\cline{1-7}
            &\lr{3}&\lr{3}&\lr{3}&\lr{3}\\\cline{2-5}
            &&\lr{2}&\lr{2}&\lr{1}\\\cline{3-5}
            &&&\lr{1}\\\cline{4-4}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        return tex_from_array([[None]*i + list(reversed(row))
                              for i, row in enumerate(self.rows())])

    def latex_standalone(self, **kwds):
        r"""
        Return Standalone picture illustrating ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: s = t.latex_standalone()
            sage: _ = s.pdf(view=False)        # long time (2s), optional - latex

        """
        from sage.misc.latex_standalone import Standalone
        if not 'standalone_config' in kwds:
            kwds['standalone_config'] = ["border=1pt"]
        return Standalone(self._latex_(), **kwds)

    def rows(self):
        """
        Return the list of rows of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: t.rows()
            [[1, 2, 4, 4, 4, 4, 4], [3, 3, 3, 3], [1, 2, 2], [1]]
        """
        if not self:
            return []

        cdef list ret = []
        cdef Py_ssize_t pos = 0
        for l in self._row_lengths:
            ret.append(self._list[pos:pos+l])
            pos += l
        return ret

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: t.e(1)
            [[4, 4, 4, 4, 4, 4, 2, 1], [3, 3, 3, 3, 3], [2, 2, 1, 1], [1]]
            sage: t.e(3)
            [[4, 4, 4, 4, 4, 3, 2, 1], [3, 3, 3, 3], [2, 2, 1], [1]]
            sage: t.e(-1)
        """
        ret = super().e(i)
        if ret is None:
            return None
        (<InfinityQueerCrystalOfTableauxElement> ret)._row_lengths = self._row_lengths
        if i < 0:
            i = -i
        L = self._parent.letters
        n = self._parent._cartan_type.n
        rows = ret.rows()
        row_lens = list(self._row_lengths)
        if count_leading(rows[n-i], L(i+1)) != len(rows[n-i+1]) + 1:
            for j in range(n-i+1):
                rows[j].append(L(n+1-j))
                row_lens[j] += 1
        return type(self)(self._parent, sum(rows, []), row_lens)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: t.f(1)
            [[4, 4, 4, 4, 4, 2, 2], [3, 3, 3, 3], [2, 2, 1], [1]]
            sage: t.f(3)
            sage: t.f(-1)
            [[4, 4, 4, 4, 4, 2, 2], [3, 3, 3, 3], [2, 2, 1], [1]]
        """
        ret = super(InfinityQueerCrystalOfTableauxElement, self).f(i)
        if ret is None:
            return None
        (<InfinityQueerCrystalOfTableauxElement> ret)._row_lengths = self._row_lengths
        if i < 0:
            i = -i
        L = self._parent.letters
        n = self._parent._cartan_type.n
        rows = ret.rows()
        row_lens = list(self._row_lengths)
        if count_leading(rows[n-i], L(i+1)) != len(rows[n-i+1]) + 1:
            for j in range(n-i+1):
                rows[j].pop()
                row_lens[j] -= 1
        return type(self)(self._parent, sum(rows, []), row_lens)

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: [t.epsilon(i) for i in B.index_set()]
            [-1, 1, -2, 0]
        """
        if i == -1:
            if self.e(-1) is None:
                return ZZ.zero()
            return ZZ.one()
        P = self._parent.weight_lattice_realization()
        h = P.simple_coroots()
        return self.phi(i) - P(self.weight()).scalar(h[i])

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['Q',4])
            sage: t = B([[4,4,4,4,4,2,1],[3,3,3,3],[2,2,1],[1]])
            sage: t.weight()
            (4, 2, 2, 0)
        """
        ret = super(InfinityQueerCrystalOfTableauxElement, self).weight()
        L = self._parent.letters
        n = self._parent._cartan_type.n + 1
        zero = self._parent.weight_lattice_realization().zero()
        La = self._parent.weight_lattice_realization().fundamental_weights()

        def fwt(i):
            return zero if i == n else La[i]
        ret -= sum((self._row_lengths[i] - 1 - self._row_lengths[i+1])*(fwt(n-i)-fwt(n-i-1))
                   for i in range(n-1))
        for i, l in enumerate(self._row_lengths[1:]):
            ret -= L(n-i).weight() * (l + 1)
        ret -= L(1).weight()  # From the 1 on the bottom row
        return ret

cdef Py_ssize_t count_leading(list row, letter) noexcept:
    cdef Py_ssize_t i
    for i in range(len(row)-1,-1,-1):
        if row[i] != letter:
            return len(row) - 1 - i
    return len(row)

# for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.crystals.tensor_product', 'ImmutableListWithParent', ImmutableListWithParent)
