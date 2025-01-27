r"""
The Pre-Lie operad

AUTHORS:

- Florent Hivert, Frédéric Chapoton (2011-2025)
"""
# ****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@lri.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from typing import NewType

from sage.categories.operads_with_basis import OperadsWithBasis
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.rooted_tree import LabelledRootedTrees
from sage.functions.other import factorial
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family

Label = NewType("Label", str)

# PreLie operad : rooted trees

# bases combinatoires : les etiquettes d'un arbre
# x.labels()

# bases combinatoires : la racine d'un arbre
# x.label()

# bases combinatoires : la largeur d'un arbre
# len(x)

# bases combinatoires : les sous-arbres t=B(t0,..,tk)
# x[0],x[1],...


class PreLieOperad(CombinatorialFreeModule):
    r"""
    An example of an operad with basis: the PreLie operad

    This class illustrates a minimal implementation of an operad with basis.
    """
    def __init__(self, R) -> None:
        """
        EXAMPLES::

            sage: A = PreLieOperad(QQ); A
            The Pre-Lie operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, LabelledRootedTrees(),
                                         latex_prefix="",
                                         sorting_key=lambda x: x.sort_key(),
                                         category=(OperadsWithBasis(R),))

    def _sort_key(self, x):
        """
        Return the key used to sort the terms.

        INPUT:

        - x -- a labelled rooted tree

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = LT([LT([],'b')], label='a')
            sage: y = LT([LT([],'d'), LT([],'e')], label='c')
            sage: A._sort_key(x)
            ((1, 'a'), (0, 'b'))
            sage: A._sort_key(y)
            ((2, 'c'), (0, 'd'), (0, 'e'))
        """
        return x.sort_key()

    def _repr_(self) -> str:
        """
        Return the string representation.

        EXAMPLES::

            sage: PreLieOperad(QQ)  # indirect doctest
            The Pre-Lie operad over Rational Field
        """
        return f"The Pre-Lie operad over {self.base_ring()}"

    def species(self):
        """
        Return the species of rooted trees.

        EXAMPLES::

            sage: f = PreLieOperad(QQ).species()
            sage: f.generating_series()[:5]
            [1, 1, 3/2, 8/3]
        """
        from sage.combinat.species.library import (
            CombinatorialSpecies,
            SetSpecies,
            SingletonSpecies,
        )
        X = SingletonSpecies()
        E = SetSpecies()
        R = CombinatorialSpecies()
        R.define(X * E(R))
        return R

    @cached_method
    def one_basis(self, letter: Label = "@"):
        """
        Return the tree with one vertex, which index the one of this operad.

        INPUT:

        - letter (default '@') -- letter used to label the vertex

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.one_basis("a")
            a[]
            sage: A.one("a")
            B[a[]]
        """
        return self.basis().keys()([], label=letter)

    def degree_on_basis(self, t):
        """
        Return the degree of a rooted tree in the Pre-Lie operad.

        This is the number of vertices.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.degree_on_basis(LT([LT([],'b')], label='a'))
            2
        """
        return t.node_number()

    def labelling_on_basis(self, t):
        """
        Put canonical labels on a tree in the Pre-Lie operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.labelling_on_basis(LT([LT([],'b')], label='a'))
            B[1[2[]]]
        """
        return self.basis()[t.canonical_labelling()]

    def unlabelling_on_basis(self, t):
        """
        Remove the labels of a tree in the Pre-Lie operad.

        The image is an element in the free pre-Lie algebra on one generator.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: A.unlabelling_on_basis(LT([LT([],'b')], label='a'))
            B[[[]]]
        """
        from sage.combinat.free_prelie_algebra import FreePreLieAlgebra
        F = FreePreLieAlgebra(self.base_ring(), '@')
        return F.basis()[t.shape()]

    def one(self, letter: Label = "@"):
        """
        I overload the one of the operad so that it can also serve as the one
        of the Hopf algebra

        MOVE TO FREE ALGEBRAS !

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.one("a")
            B[a[]]
            sage: A.one()
            B[@[]]
        """
        # super(PreLieOperad, self).one(letter) bug of cached_method
        return self.monomial(self.one_basis(letter))

    def some_elements(self):
        """
        Return some elements of the operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.some_elements()
            [B[a[]],
             B[a[b[]]],
             B[a[b[c[d[]]]]] + B[a[b[], c[d[]]]],
             B[a[b[d[]]]] + B[a[b[], d[]]]]
        """
        x = self.one("a") < self.one("b")
        y = self.one("c") < self.one("d")
        return [self.one("a"), x, x < y, x < self.one("d")]

    # procedures de composition
    # nota bene : these algorithms only work if all labels of y are distinct !

    def singleGraft(self, y, x, graftingFunction, path_prefix=()):
        """
        Return the rooted tree obtained from a rooted tree `y`, a
        rooted tree `x` and a grafting function ``graftingFunction`` from
        range(len(x)) to the set of paths in `y`.

        .. WARNING:

            These algorithms only work if all labels of `y` are distinct !

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.singleGraft(y,x,dict([[0,(0,)]]))
            a[b[d[]]]
            sage: t = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: A.singleGraft(t,s,dict([[0,(0,)],[1,(1,)]]))
            a[b[d[]], c[e[]]]
        """
        y1 = self.basis().keys()([self.singleGraft(suby, x, graftingFunction,
                                                   path_prefix + (i,))
                                  for i, suby in enumerate(y)],
                                 label=y.label())
        with y1.clone() as y2:
            for k in range(len(x)):
                if graftingFunction[k] == path_prefix:
                    y2.append(x[k])
        return y2

    def composition_on_basis_in_root(self, x, y):
        """
        Return a list of rooted trees obtained from a rooted tree `x`
        and a rooted tree `y` by the composition `x o_i y` where `i`
        is the label of the root of `x`.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.composition_on_basis_in_root(x,y)
            [a[b[], d[]], a[b[d[]]]]
        """
        return [self.singleGraft(y, x, graftingFunction)
                for graftingFunction in
                cartesian_product([list(y.paths())] * len(x))]

    def composition_on_basis_iter(self, x, y, i: Label):
        """
        Return an iterator of rooted trees obtained from a rooted tree `x`
        and a rooted tree `y` by the composition `x o_i y`.

        The composition index `i` should be a label of `x`.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: list(A.composition_on_basis_iter(x,y,'c'))
            [a[b[], d[]], a[b[d[]]]]
        """
        if i not in x.labels():
            raise ValueError("the composition index is not present")

        if x.label() == i:
            yield from self.composition_on_basis_in_root(x, y)
            return

        for j in range(len(x)):
            if i in x[j].labels():
                for sx in self.composition_on_basis_iter(x[j], y, i):
                    with x.clone() as x1:
                        x1[j] = sx
                    yield x1

    def composition_on_basis(self, x, y, i: Label):
        """
        Return the composition `x o_i y` as a sum of rooted trees,
        for rooted trees `x` and `y`.

        The composition index `i` should be a label of `x`.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: y = LT([LT([],'b')], label='a')
            sage: x = LT([LT([],'d')], label='c')
            sage: A.composition_on_basis(x,y,'a')
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
            sage: A.composition(A(x), A(y), 'd')
            B[c[a[b[]]]]
            sage: A.composition(A(x), A(y), 'c')
            B[a[b[d[]]]] + B[a[b[], d[]]]
        """
        if i not in x.labels():
            raise ValueError("the composition index is not present")
        return sum(self.basis()[t] for t in
                   self.composition_on_basis_iter(x, y, i))

    def operad_generator_basis(self, label0: Label = "0", label1: Label = "1"):
        """
        Return the generator of ``self`` as a labelled rooted tree.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.operad_generator_basis()
            0[1[]]
        """
        LT = self.basis().keys()
        return LT([LT([], label=label1)], label=label0)

    def operad_generator(self, label0: Label = "0", label1: Label = "1"):
        """
        Return the generator of ``self``.

        This is the rooted tree with two vertices.

        INPUT:

        - labels of the vertices -- optional, by default 0 for the root vertex
          and 1 for the other vertex

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: A.operad_generator()
            B[0[1[]]]
        """
        return self.monomial(self.operad_generator_basis(label0, label1))

    def operad_generators(self) -> Family:
        """
        Return the family of generators.

        EXAMPLES::

            sage: PreLieOperad(QQ).operad_generators()
            Finite family {'pre_Lie_product': B[1[2[]]]}
        """
        return Family({"pre_Lie_product": self.operad_generator(1, 2)})

    def pre_Lie_product(self, x, y):
        """
        Return the pre-Lie product of ``x`` and ``y`` inside the operad.

        EXAMPLES::

            sage: A = PreLieOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b')], label='a'))
            sage: y = A(LT([LT([],'d')], label='c'))
            sage: A.pre_Lie_product(x, y)
            B[a[b[c[d[]]]]] + B[a[b[], c[d[]]]]
        """
        return self.operad_generator().compose(x, "0").compose(y, "1")

    def operad_morphism_on_basis(self, t, cod):
        """
        Apply a morphism from the PreLie operad to the target operad.

        The target operad has to possess a method called ``pre_Lie_product``.

        The argument t (a rooted tree) should not have repeated labels.

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: RT = PL.basis().keys()
            sage: tr = RT([],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL) == PL(tr)
            True
            sage: tr = RT([RT([],label='b')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL) == PL(tr)
            True
            sage: tr = RT([RT([],label='b'),RT([],label='c')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
            sage: tr = RT([RT([RT([],label='c')],label='b'),RT([],label='e'),RT([],label='d')],label='a')
            sage: PL.operad_morphism_on_basis(tr,PL)==PL(tr)
            True
        """
        targetPreLieProduct = cod.pre_Lie_product
        width = len(t)
        if width == 0:
            return cod.one(t.label())
        elif width == 1:
            return targetPreLieProduct(cod.one(t.label()),
                                       self.operad_morphism_on_basis(t[0],
                                                                     cod))

        t_red = self.basis().keys()([t[i] for i in range(width - 1)],
                                    label=t.label())
        somme1 = targetPreLieProduct(self.operad_morphism_on_basis(t_red, cod),
                                     self.operad_morphism_on_basis(t[width - 1],
                                                                   cod))
        somme2 = cod(0)
        for j in range(width - 1):
            with t_red.clone() as tj:
                tj[j] = self.basis().keys()([], label=t_red[j].label())
                somme2 += cod.composition(self.operad_morphism_on_basis(tj, cod), targetPreLieProduct(self.operad_morphism_on_basis(t_red[j], cod), self.operad_morphism_on_basis(t[width - 1], cod)), t_red[j].label())
        return somme1 - somme2

    def corolla(self, n, x, y, N=10):
        r"""
        Evaluate the corolla with `n` leaves, with `x` in the root and `y`
        in the leaves.

        SHOULD BE A METHOD OF PRELIE-ALGEBRAS

        The result is computed up to order `N` (included).

        INPUT:

        - `n` -- an integer

        - `x`, `y` -- two tree-indexed series

        - `N` (optional) -- an integer (default: 10) truncation order

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(ZZ)
            sage: PLT = PL.basis().keys()

            sage: a = PL(PLT([PLT([],label="@")],label="@"))
            sage: PL.corolla(1,a,a)
            B[@[@[@[@[]]]]] + B[@[@[], @[@[]]]]

            sage: b = PL.one()
            sage: PL.corolla(3,b,b,4)
            B[@[@[], @[], @[]]]

            sage: PL.corolla(2,a,b,4)
            B[@[@[@[], @[]]]] + 2*B[@[@[], @[@[]]]] + B[@[@[], @[], @[]]]
        """
        PL = x.parent()
        if n + 1 > N:
            return PL.zero()

        xx = x.truncate(N + 1 - n)
        yy = y.truncate(N + 1 - n)
        PLT = PL.basis().keys()
        crln = PL(PLT([PLT([], str(i + 1)) for i in range(n)], "0"))
        resu = crln.compose(xx, "0")
        for i in range(1, n + 1):
            resu = PL.composition_truncated(resu, yy, str(i), N)
        return resu

    def sum_corolla_prelie(self, n, x, y, N=10):
        r"""
        Evaluate the sum of all corollas with up to `n` leaves, with
        `x` in the root and `y` in the leaves.

        The result is computed up to order `N` (included).

        INPUT:

        - `n` -- an integer

        - `a`, `b` -- two tree-indexed series

        - `N` (optional) -- an integer (default: 10) truncation order

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: a = PL.one("@")
            sage: b = PL.one("O")
            sage: PL.sum_corolla_prelie(3,a,b,3)
            B[@[]] + B[@[O[]]] + 1/2*B[@[O[], O[]]]
            sage: PL.sum_corolla_prelie(3,a,b,10)
            B[@[]] + B[@[O[]]] + 1/2*B[@[O[], O[]]] + 1/6*B[@[O[], O[], O[]]]
        """
        PL = x.parent()
        return x + sum(PL.corolla(i, x, y, N) / factorial(i)
                       for i in range(1, n + 1))

    def diese_product(self, x, y, N=10):
        r"""
        Evaluate the `\#` product of `x` and `y` up to order `N`

        The `\#` product is an associative product on tree-indexed
            `\#` series. This is similar to the
            Baker-Campbell-Hausdorff formula.

        The tree-indexed series `x` is in the root. The tree-indexed
        series `y` is in the leaves.

        The result is computed up to order `N` (included).

        SHOULD BE A METHOD OF PRELIE-ALGEBRAS

        INPUT:

        - `x`, `y` -- two tree-indexed series

        - `N` (optional) -- an integer (default: 10) truncation order

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: PLT = PL.basis().keys()

            sage: a = PL(PLT([PLT([],label='@')],label='@'))
            sage: PL.diese_product(a,a,2)
            2*B[@[@[]]]
            sage: PL.diese_product(a,a,4)
            2*B[@[@[]]] + B[@[@[@[@[]]]]] + B[@[@[], @[@[]]]]
            sage: b = PL.one("O")
            sage: PL.diese_product(b,a,5)
            B[O[]] + B[@[@[]]] + B[O[@[@[]]]] + 1/2*B[O[@[@[]], @[@[]]]]
        """
        PL = x.parent()
        yy = y.truncate(N + 1)
        xx = x.truncate(N + 1)
        return yy + xx + sum(PL.corolla(i, xx, yy, N) / factorial(i)
                             for i in range(1, N))

    def inverse_diese(self, x, N=10):
        r"""
        Return the inverse of `x` for the `\#` product up to order `N`

        The result is computed up to order `N` (included).

        INPUT:

        - `x` -- a tree-indexed series

        - `N` (optional) -- an integer (default: 10) truncation order

        OUTPUT:

        - a tree-indexed series

        EXAMPLES::

            sage: PL = PreLieOperad(QQ)
            sage: a = PL.one()
            sage: PL.inverse_diese(a,2)
            -B[@[]] + B[@[@[]]]
        """
        PL = x.parent()
        inverse = PL.zero()
        for i in range(1, N + 1):
            inverse -= PL.diese_product(x, inverse, i)
        return inverse

    class Element(CombinatorialFreeModule.Element):
        def __lt__(self, other):
            r"""
            Shortcut for the prelie product

            EXAMPLES::

                sage: A = PreLieOperad(QQ)
                sage: A.one("x") < A.one("x")
                B[x[x[]]]

            .. warning::

                Due to priority rules for operators, term must be put
                within parentheses inside sum, product... For example you must
                write::

                    sage: a = A.one('a'); b = A.one('b'); c = A.one('c')
                    sage: (a<b) + c
                    B[c[]] + B[a[b[]]]

                Indeed ``a<b + c`` is understood as ``a< (b + c)``::

                    sage: (a<b + c) - (a < (b + c))
                    0
            """
            parent = self.parent()
            if parent != other.parent():
                raise TypeError('not in the same parent')
            return parent.pre_Lie_product(self, other)
