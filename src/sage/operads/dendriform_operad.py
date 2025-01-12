from sage.categories.operads_with_basis import OperadsWithBasis
from sage.combinat.binary_tree import LabelledBinaryTrees
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method


class DendriformOperad(CombinatorialFreeModule):
    r"""
    The Dendriform operad

    This is an operad on the species of planar binary trees.

    EXAMPLES::

        sage: A = DendriformOperad(QQ)
        sage: LT = A.basis().keys()
        sage: x = A(LT([LT([],'b'),None], label='a'))
        sage: y = A(LT([None, LT([],'c')], label='d'))
        sage: x.compose(y, 'a')
        B[d[b[., .], c[., .]]]

        sage: z = A(LT([LT([], 'e'), LT([],'c')], label='d'))
        sage: x.compose(z,'a')
        B[d[b[., e[., .]], c[., .]]] + B[d[e[b[., .], .], c[., .]]]

    REFERENCES:

    - [Loday2001]_
    """
    def __init__(self, R) -> None:
        """
        EXAMPLES::

            sage: A = DendriformOperad(QQ); A
            The Dendriform operad over Rational Field
            sage: TestSuite(A).run()

            sage: A = DendriformOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'u'),None], label='v'))
            sage: y = A(LT([None, LT([],'j')], label='i'))
            sage: x.compose(y, 'v')
            B[i[u[., .], j[., .]]]
        """
        def key(t):
            return (t.to_dyck_word(), t.labels())
        CombinatorialFreeModule.__init__(self, R, LabelledBinaryTrees(),
                                         category=OperadsWithBasis(R),
                                         sorting_key=key)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: DendriformOperad(QQ)  # indirect doctest
            The Dendriform operad over Rational Field
        """
        return f"The Dendriform operad over {self.base_ring()}"

    def species(self):
        """
        The species of planar binary trees with labels at inner vertices

        This is the species underlying the Dendriform operad.

        EXAMPLES::

            sage: f = DendriformOperad(QQ).species()
            sage: f.generating_series()[:5]
            [1, 2, 5, 14]
        """
        from sage.combinat.species.library import (
            CombinatorialSpecies,
            EmptySetSpecies,
            SingletonSpecies,
        )
        X = SingletonSpecies()
        u = EmptySetSpecies()
        R = CombinatorialSpecies()
        R.define((u + R) * X * (u + R))
        return R

    @cached_method
    def one_basis(self, letter='@'):
        """
        Return the planar binary tree with one vertex, which indexes
        the one of this operad, as per
        :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        INPUT:

        - ``letter`` (default ``'@'``) -- letter used to label the unit

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: A.one_basis("a")
            a[., .]
        """
        return self.basis().keys()([], label=letter)

    def degree_on_basis(self, t):
        """
        Return the degree of a tree in the Dendriform operad.

        This is the number of nodes (inner vertices).

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: A.degree_on_basis(Trees([None,Trees([],label="c")],label="d"))
            2
        """
        return t.node_number()

    def labelling_on_basis(self, t):
        """
        Put canonical labels on a tree in the Dendriform operad.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: A.labelling_on_basis(Trees([None,Trees([],label="c")],label="d"))
            B[1[., 2[., .]]]
        """
        return self.basis()[t.canonical_labelling()]

    def unlabelling_on_basis(self, t):
        """
        Remove the labels of a tree in the Dendriform operad.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: A.unlabelling_on_basis(Trees([None,Trees([],label="c")],label="d"))
            B[[., [., .]]]
        """
        return self.basis()[t.shape()]

    def shuffle_on_basis_list(self, x, y):
        """
        Return the shuffle product (associative) of two planar binary
        trees as a list of planar binary trees.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: ouest = A.one_basis("o")
            sage: est = A.one_basis("e")
            sage: A.shuffle_on_basis_list(ouest,est)
            [e[o[., .], .], o[., e[., .]]]
        """
        if x.is_empty():
            return [y]
        elif y.is_empty():
            return [x]
        return self.composition_on_basis_list(self.basis().keys()([x, None], label='diese'), y, 'diese') + self.composition_on_basis_list(self.basis().keys()([None, y], label='diese'), x, 'diese')

    def composition_on_basis_list(self, x, y, i):
        """
        Return the composition of two planar binary
        trees in the dendriform operad as a list of planar binary trees.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: A.composition_on_basis_list(Trees([],label="a"), Trees([None,Trees([],label="c")],label="d"),"a")
            [d[., c[., .]]]
        """
        if i in x[0].labels():
            return [self.basis().keys()([z, x[1]], label=x.label())
                    for z in self.composition_on_basis_list(x[0], y, i)]
        elif i in x[1].labels():
            return [self.basis().keys()([x[0], z], label=x.label())
                    for z in self.composition_on_basis_list(x[1], y, i)]
        return [self.basis().keys()([z, t], label=y.label())
                for z in self.shuffle_on_basis_list(x[0], y[0])
                for t in self.shuffle_on_basis_list(y[1], x[1])]

    def composition_on_basis(self, x, y, i):
        r"""
        Return the composition `x o_i y` as a sum of planar binary trees.

        INPUT:

        - planar binary trees `x` and `y`
        - composition index `i`

        The composition index `i` should be a label of `x`.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: A.composition_on_basis(Trees([None,Trees([],label="b")],label="a"), Trees([None,Trees([],label="c")],label="d"),"a")
            B[d[., c[., b[., .]]]] + B[d[., b[c[., .], .]]]

        TESTS::

            sage: A.composition_on_basis(Trees([],label="a"), Trees([None,Trees([],label="c")],label="d"),"a")
            B[d[., c[., .]]]

            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b'),None], label='a'))
            sage: y = A(LT([None, LT([],'c')], label='d'))
            sage: x.compose(y, 'a')
            B[d[b[., .], c[., .]]]

            sage: A.composition_on_basis(Trees([],label="a"), Trees([None,Trees([],label="c")],label="d"),"e")
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
        """
        if i not in x.labels():
            raise ValueError("the composition index is not present")
        return sum(self.basis()[t] for t in
                   self.composition_on_basis_list(x, y, i))

    def pre_Lie_product(self, x, y):
        """
        Return the pre-Lie product of ``x`` and ``y``.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b'),None], label='a'))
            sage: y = A(LT([LT([],'d'),None], label='c'))
            sage: A.pre_Lie_product(x, y)
            -B[a[b[., .], c[d[., .], .]]] + B[a[c[d[., .], b[., .]], .]] + B[a[b[c[d[., .], .], .], .]]
        """
        LT = self.basis().keys()
        expr = self.monomial(LT([LT([], label="1"), None], label="0")) - self.monomial(LT([None, LT([], label="1")], label="0"))
        return expr.compose(x, "0").compose(y, "1")

    def associative_product(self, x, y):
        """
        Return the associative product of ``x`` and ``y``.

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b'),None], label='a'))
            sage: y = A(LT([LT([],'d'),None], label='c'))
            sage: A.associative_product(x, y)
            B[a[b[., .], c[d[., .], .]]] + B[c[a[b[., .], d[., .]], .]] + B[c[d[a[b[., .], .], .], .]]
        """
        LT = self.basis().keys()
        expr = self.monomial(LT([LT([], label="0"), None], label="1")) + self.monomial(LT([None, LT([], label="1")], label="0"))
        return expr.compose(x, "0").compose(y, "1")

    def chosen_product(self, x, y, name='assoc'):
        """
        Return one among different binary products.

        The possible choices are:

        "assoc", "prelie", "lie", "left_dend", "right_dend"

        EXAMPLES::

            sage: A = DendriformOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b'),None], label='a'))
            sage: y = A(LT([LT([],'d'),None], label='c'))
            sage: A.chosen_product(x, y, 'assoc')
            B[a[b[., .], c[d[., .], .]]] + B[c[a[b[., .], d[., .]], .]]
            + B[c[d[a[b[., .], .], .], .]]
            sage: A.chosen_product(x, y, 'lie')
            B[a[b[., .], c[d[., .], .]]] - B[c[d[., .], a[b[., .], .]]] - B[a[c[d[., .], b[., .]], .]] + B[c[a[b[., .], d[., .]], .]] - B[a[b[c[d[., .], .], .], .]] + B[c[d[a[b[., .], .], .], .]]

            sage: A.chosen_product(x, y, 'prelie')
            -B[a[b[., .], c[d[., .], .]]] + B[a[c[d[., .], b[., .]], .]] + B[a[b[c[d[., .], .], .], .]]
            sage: A.chosen_product(x, y, 'left_dend')
            B[c[a[b[., .], d[., .]], .]] + B[c[d[a[b[., .], .], .], .]]
            sage: A.chosen_product(x, y, 'right_dend')
            B[a[b[., .], c[d[., .], .]]]
        """
        LT = self.basis().keys()
        if name == 'assoc':
            expr = self.monomial(LT([LT([], label="0"), None], label="1"))
            expr += self.monomial(LT([None, LT([], label="1")], label="0"))
        elif name == 'prelie':
            expr = self.monomial(LT([LT([], label="1"), None], label="0"))
            expr -= self.monomial(LT([None, LT([], label="1")], label="0"))
        elif name == 'lie':
            expr = self.monomial(LT([LT([], label="0"), None], label="1"))
            expr += self.monomial(LT([None, LT([], label="1")], label="0"))
            expr -= self.monomial(LT([LT([], label="1"), None], label="0"))
            expr -= self.monomial(LT([None, LT([], label="0")], label="1"))
        elif name == 'left_dend':
            expr = self.monomial(LT([LT([], label="0"), None], label="1"))
        elif name == 'right_dend':
            expr = self.monomial(LT([None, LT([], label="1")], label="0"))
        return expr.compose(x, "0").compose(y, "1")
