from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.rooted_tree import LabelledRootedTrees
from sage.misc.cachefunc import cached_method


class NapOperad(CombinatorialFreeModule):
    r"""
    The Nap operad

    The word Nap stands here for 'Not Associative Permutative'

    This is an operad on the species of rooted trees.

    This class illustrates a minimal implementation of an operad with basis.

    EXAMPLES::

        sage: NAP = NapOperad(QQ)
        sage: B = NAP.basis()
        sage: LT = NAP.basis().keys()
        sage: t = LT([LT([],label='b'),LT([],label='c')], label='a')
        sage: s = LT([LT([],label='d'),LT([],label='e')], label='f')
        sage: NAP.composition(B[t],B[s],"b")
        B[a[c[], f[d[], e[]]]]
        sage: NAP.composition(B[s],B[t],"d")
        B[f[e[], a[b[], c[]]]]
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = NapOperad(QQ); A
            The Nap operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, LabelledRootedTrees(),
                                         latex_prefix="",
                                         sorting_key=lambda x: x.sort_key(),
                                         category=(OperadsWithBasis(R),))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: NapOperad(QQ)         # indirect doctest
            The Nap operad over Rational Field
        """
        return f"The Nap operad over {self.base_ring()}"

    def species(self):
        """
        Return the species of rooted trees.

        This is the species underlying the Nap operad.

        EXAMPLES::

            sage: f = NapOperad(QQ).species()
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
    def one_basis(self, letter='@'):
        """
        Return the tree with one vertex, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        INPUT:

        - ``letter`` (default ``'@'``) -- letter used to label the unit

        EXAMPLES::

            sage: A = NapOperad(QQ)
            sage: A.one_basis("a")
            a[]
        """
        return self.basis().keys()([], label=letter)

    def composition_on_basis_in_root(self, x, y):
        r"""
        Return the rooted tree
        obtained from a rooted tree `x` and a rooted tree `y` by the
        composition `x o_i y` where `i` is the root of `x`.

        EXAMPLES::

            sage: NAP = NapOperad(QQ)
            sage: LT = NAP.basis().keys()
            sage: t = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: NAP.composition_on_basis_in_root(t,s)
            f[b[], c[], d[], e[]]

        TESTS::

            sage: toto = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: titi = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: test1 = NAP.composition_on_basis_in_root(toto,titi)
        """
        return self.basis().keys()(list(x) + list(y), label=y.label())

    def composition_on_basis_as_tree(self, x, y, i):
        r"""
        Return the rooted tree obtained from a rooted tree `x`
        and a rooted tree `y` by the composition `x o_i y`.

        The composition index `i` should be a label of `x`.

        EXAMPLES::

            sage: NAP = NapOperad(QQ)
            sage: LT = NAP.basis().keys()
            sage: t = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: NAP.composition_on_basis_as_tree(t,s,"b")
            a[c[], f[d[], e[]]]

        TESTS::

            sage: NAP.composition_on_basis_as_tree(t,s,"f")
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
        """
        if i not in x.labels():
            raise ValueError("the composition index is not present")
        elif x.label() == i:
            return self.composition_on_basis_in_root(x, y)

        j = [k for k in range(len(x)) if i in x[k].labels()][0]
        with x.clone() as x1:
            x1[j] = self.composition_on_basis_as_tree(x[j], y, i)
        return x1

    def composition_on_basis(self, x, y, i):
        r"""
        Return a monomial obtained from a rooted tree `x`
        and a rooted tree `y` by the composition `x o_i y`.

        The composition index `i` should be a label of `x`.

        EXAMPLES::

            sage: NAP = NapOperad(QQ)
            sage: LT = NAP.basis().keys()
            sage: t = LT([LT([],label='b'),LT([],label='c')], label='a')
            sage: s = LT([LT([],label='d'),LT([],label='e')], label='f')
            sage: NAP.composition_on_basis(t,s,"b")
            B[a[c[], f[d[], e[]]]]

        TESTS::

            sage: NAP.composition_on_basis(t,s,"f")
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
        """
        return self.basis()[self.composition_on_basis_as_tree(x, y, i)]
