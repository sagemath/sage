from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ordered_tree import LabelledOrderedTrees
from sage.misc.cachefunc import cached_method

LT = LabelledOrderedTrees()


def map_leaves(self, f):
    if self.is_empty():
        return self
    if len(self) == 0:
        return self.parent()([t.map_labels(f) for t in self],
                             label=f(self.label()))
    return self.parent()([t.map_labels(f) for t in self],
                         label=self.label())


LT.map_labels = map_leaves


class FreeOperad(CombinatorialFreeModule):
    r"""
    The free operad over any given set of generators
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = FreeOperad(QQ); A
            The Free operad over Rational Field with generators
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, LT,
                                         category=OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: FreeOperad(QQ)            # indirect doctest
            The Free operad over Rational Field with generators
        """
        msg = "The Free operad over {} with generators "
        return msg.format(self.base_ring())

    @cached_method
    def one_basis(self, letter=None):
        """
        Return the planar tree with no vertex and one leaf, which indexes
        the one of this operad.

        EXAMPLES::

            sage: A = FreeOperad(QQ)
            sage: A.one_basis("a")
            a[]
        """
        if letter is None:
            letter = '@'
        return self.basis().keys()([], label=letter)

    def composition_on_basis_as_tree(self, x, y, i):
        """
        Return the composition of two planar trees in the free operad
        as a planar tree.

        This is just grafting `y` on the leaf `i` of `x`.

        INPUT:

        - `x,y` -- two planar trees

        - `i` -- composition index

        The composition index `i` should be a leaf label of `x`.

        EXAMPLES::

            sage: A = FreeOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: one = Trees([],label='a')
            sage: t = Trees([Trees([],label='c')],label='d')
            sage: A.composition_on_basis_as_tree(one,one,'a')
            a[]
            sage: A.composition_on_basis_as_tree(one,t,'a')
            d[c[]]
            sage: A.composition_on_basis_as_tree(t,one,'c')
            d[a[]]
            sage: t1 = Trees([Trees([],label='c')],label='d')
            sage: t2 = Trees([Trees([],label='e'),Trees([],label='b')],
            ....:   label='f')
            sage: A.composition_on_basis_as_tree(t1,t2,'c')
            d[f[e[], b[]]]
            sage: A.composition_on_basis_as_tree(t2,t1,'b')
            f[e[], d[c[]]]
        """
        if x.node_number() == 1:
            return y

        with x.clone() as t:
            for j in range(len(t)):
                if i in t[j].leaf_labels():
                    t[j] = self.composition_on_basis_as_tree(t[j], y, i)
        return t

    def composition_on_basis(self, x, y, i):
        r"""
        Return the composition `x o_i y` for planar trees `x` and `y`.

        The composition index `i` should be a leaf label of `x`.

        EXAMPLES::

            sage: A = FreeOperad(QQ)
            sage: Trees = A.basis().keys()
            sage: one = Trees([],label='a')
            sage: t = Trees([Trees([],label='c')],label='d')
            sage: A.composition_on_basis(one,one,'a')
            B[a[]]
            sage: A.composition_on_basis(one,t,'a')
            B[d[c[]]]
            sage: A.composition_on_basis(t,one,'c')
            B[d[a[]]]
            sage: t1 = Trees([Trees([],label='c')],label='d')
            sage: t2 = Trees([Trees([],label='e'),Trees([],label='b')],
            ....:   label='f')
            sage: A.composition_on_basis(t1,t2,'c')
            B[d[f[e[], b[]]]]
            sage: A.composition_on_basis(t2,t1,'b')
            B[f[e[], d[c[]]]]

        TESTS::

            sage: A = FreeOperad(QQ)
            sage: Trees = A.basis().keys()
        """
        if i not in x.leaf_labels():
            raise ValueError("the composition index is not present")

        return self.basis()[self.composition_on_basis_as_tree(x, y, i)]

    def magmatic_product(self, x, y):
        """
        Return the binary magmatic product of ``x`` and ``y`` inside the operad.

        EXAMPLES::

            sage: A = FreeOperad(QQ)
            sage: LT = A.basis().keys()
            sage: x = A(LT([LT([],'b')], label='a'))
            sage: y = A(LT([LT([],'d')], label='c'))
            sage: A.magmatic_product(x, y)
            B[@[a[b[]], c[d[]]]]
        """
        LT = self.basis().keys()
        t = LT([LT([], label=0), LT([], label=1)], label='@')
        gen = self.monomial(t)
        return gen.compose(x, 0).compose(y, 1)

    def operad_morphism_on_basis(self, t, cod, fun=None):
        """
        Define a morphism from the free operad to the target operad.

        INPUT:

        - If given ``fun`` is a function associating the generator of
          ``self`` to element of the codomain ``cod``. If not given
          the generator are supposed to be indexed by string and the
          image of the generator ``gen`` is computed via ``cod.gen``.

        EXAMPLES::

            sage: A = FreeOperad(QQ)
            sage: PL = PreLieOperad(QQ)
            sage: Tr = A.basis().keys()
            sage: t = Tr([Tr([], label='a'), Tr([], label='b')],
            ....:        label='pre_Lie_product')
            sage: A.operad_morphism_on_basis(t, PL,
            ....:      PL.operad_generators().__getitem__)
            B[a[b[]]]
            sage: A.operad_morphism(PL,
            ....:      PL.operad_generators().__getitem__)(A(t))
            B[a[b[]]]
            sage: A.operad_morphism(PL,
            ....:      PL.operad_generators().__getitem__)(2*A(t))
            2*B[a[b[]]]
        """
        if fun is None:
            def fun(st):
                return getattr(cod, st)
        if len(t) == 0:
            return cod.one(t.label())
        res = fun(t.label())
        for i in range(len(t), 0, -1):
            ti = t[i - 1]  # trees starts from zero / operads from 1
            ri = self.operad_morphism_on_basis(ti, cod, fun)
            res = cod.composition(res, ri, i)
        return res
