# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Affinization Crystals
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
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

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.rings.infinity import Infinity


class AffinizationOfCrystal(UniqueRepresentation, Parent):
    r"""
    An affinization of a crystal.

    Let `\mathfrak{g}` be a Kac-Moody algebra of affine type. The
    affinization of a finite `U_q^{\prime}(\mathfrak{g})`-crystal `B`
    is the (infinite) `U_q(\mathfrak{g})`-crystal with underlying set:

    .. MATH::

        B^{\mathrm{aff}} = \{ b(m) \mid b \in B, m \in \ZZ \}

    and crystal structure determined by:

    .. MATH::

        \begin{aligned}
            e_i(b(m)) & =
            \begin{cases}
              (e_0 b)(m+1) & i = 0, \\
              (e_i b)(m)   & i \neq 0,
            \end{cases} \\
            f_i(b(m)) &=
            \begin{cases}
              (f_0 b)(m-1) & i = 0, \\
              (f_i b)(m)   & i \neq 0,
            \end{cases} \\
            \mathrm{wt}(b(m)) &= \mathrm{wt}(b) + m \delta.
        \end{aligned}

    EXAMPLES:

    We first construct a Kirillov-Reshetikhin crystal and then take it's
    corresponding affinization::

        sage: K = crystals.KirillovReshetikhin(['A',2,1], 2, 2)
        sage: A = K.affinization()

    Next we construct an affinization crystal from a tensor product of KR
    crystals::

        sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C',2,1], [[1,2],[2,1]])
        sage: A = crystals.AffinizationOf(KT)

    REFERENCES:

    - [HK2002]_ Chapter 10
    """

    def __init__(self, B):
        """
        Initialize ``self``.

        EXAMPLES:

        We skip the Stembridge axioms test since this is an abstract crystal::

            sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
            sage: TestSuite(A).run(skip='_test_stembridge_local_axioms') # long time
        """
        if not B.cartan_type().is_affine():
            raise ValueError("must be an affine crystal")
        if B.cardinality() == Infinity:
            raise ValueError("must be finite crystal")
        self._B = B
        self._cartan_type = B.cartan_type()
        Parent.__init__(self, category=(RegularCrystals(), InfiniteEnumeratedSets()))
        self.module_generators = tuple([self.element_class(self, b, 0)
                                        for b in B.module_generators])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.KirillovReshetikhin(['A',2,1], 1, 1).affinization()
            Affinization of Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1)
        """
        return "Affinization of {}".format(self._B)

    class Element(Element):
        """
        An element in an affinization crystal.
        """

        def __init__(self, parent, b, m):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: mg = A.module_generators[0]
                sage: TestSuite(mg).run()
            """
            self._b = b
            self._m = m
            Element.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: A.module_generators[0]
                [[1, 1], [2, 2]](0)
                sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C',2,1], [[1,2],[2,1]])
                sage: A = crystals.AffinizationOf(KT)
                sage: A.module_generators[0]
                [[1, 1]] (X) [[1], [2]](0)
            """
            return "{!r}({})".format(self._b, self._m)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: latex(A.module_generators[0])
                {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
                \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
                \lr{1}&\lr{1}\\\cline{1-2}
                \lr{2}&\lr{2}\\\cline{1-2}
                \end{array}$}
                } (0)
            """
            from sage.misc.latex import latex
            return latex(self._b) + "({})".format(self._m)

        def __hash__(self):
            r"""
            TESTS::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: mg = A.module_generators[0]
                sage: hash(mg) == hash(mg._b) ^^ hash(mg._m)
                True
            """
            return hash(self._b) ^ hash(self._m)

        def _richcmp_(self, other, op):
            """
            Comparison.

            TESTS::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg == mg
                True
                sage: mg == mg.f(2).e(2)
                True
                sage: KT = crystals.TensorProductOfKirillovReshetikhinTableaux(['C',2,1], [[1,2],[2,1]])
                sage: A = crystals.AffinizationOf(KT)
                sage: A(KT.module_generators[3], 1).f(0) == A.module_generators[0]
                True

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg != mg.f(2)
                True
                sage: mg != mg.f(2).e(2)
                False


                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2, 2).affinization()
                sage: S = A.subcrystal(max_depth=2)
                sage: sorted(S)
                [[[1, 1], [2, 2]](0),
                 [[1, 1], [2, 3]](0),
                 [[1, 2], [2, 3]](0),
                 [[1, 1], [3, 3]](0),
                 [[1, 1], [2, 3]](1),
                 [[1, 2], [2, 3]](1),
                 [[1, 2], [3, 3]](1),
                 [[2, 2], [3, 3]](2)]
            """
            return richcmp((self._m, self._b), (other._m, other._b), op)

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2,2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg.e(0)
                [[1, 2], [2, 3]](1)
                sage: mg.e(1)
                sage: mg.e(0).e(1)
                [[1, 1], [2, 3]](1)
            """
            bp = self._b.e(i)
            if bp is None:
                return None
            if i == 0:
                return self.__class__(self.parent(), bp, self._m+1)
            return self.__class__(self.parent(), bp, self._m)

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2,2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg.f(2)
                [[1, 1], [2, 3]](0)
                sage: mg.f(2).f(2).f(0)
                sage: mg.f_string([2,1,1])
                sage: mg.f_string([2,1])
                [[1, 2], [2, 3]](0)
                sage: mg.f_string([2,1,0])
                [[1, 1], [2, 2]](-1)
            """
            bp = self._b.f(i)
            if bp is None:
                return None
            if i == 0:
                return self.__class__(self.parent(), bp, self._m-1)
            return self.__class__(self.parent(), bp, self._m)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2,2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg.epsilon(0)
                2
                sage: mg.epsilon(1)
                0
            """
            return self._b.epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2,2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg.phi(0)
                0
                sage: mg.phi(2)
                2
            """
            return self._b.phi(i)

        def weight(self):
            r"""
            Return the weight of ``self``.

            The weight `\mathrm{wt}` of an element is:

            .. MATH::

                \mathrm{wt}\bigl( b(m) \bigr) = \mathrm{wt}(b) + m \delta,

            where `\delta` is the null root.

            EXAMPLES::

                sage: A = crystals.KirillovReshetikhin(['A',2,1], 2,2).affinization()
                sage: mg = A.module_generators[0]
                sage: mg.weight()
                -2*Lambda[0] + 2*Lambda[2]
                sage: mg.e(0).weight()
                -Lambda[1] + Lambda[2] + delta
                sage: mg.e(0).e(0).weight()
                2*Lambda[0] - 2*Lambda[1] + 2*delta
            """
            WLR = self.parent().weight_lattice_realization()
            La = WLR.fundamental_weights()
            return WLR.sum(c*La[i] for i,c in self._b.weight()) + self._m * WLR.null_root()
