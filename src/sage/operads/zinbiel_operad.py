from sage.categories.operads_with_basis import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method


class ZinbielOperad(CombinatorialFreeModule):
    r"""
    The Zinbiel operad

    This is an operad on the species of non-empty lists.

    EXAMPLES::

        sage: Z = ZinbielOperad(QQ)
        sage: Z('abc')
        B[word: abc]
        sage: Z('abc').compose(Z('de'), 'a')
        B[word: dbce] + B[word: dbec] + B[word: debc]
        sage: Z('abc').compose(Z('de'), 'c')
        B[word: abde]

    REFERENCES:

    .. [todo]
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = ZinbielOperad(QQ); A
            The Zinbiel operad over Rational Field
            sage: TestSuite(A).run()

            sage: W = Words()
            sage: A.composition(A(W("abc")), A(W("def")), "b")
            B[word: adcef] + B[word: adecf] + B[word: adefc]
            sage: A.composition(A("abc"), A("def"), "b")
            B[word: adcef] + B[word: adecf] + B[word: adefc]
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category=OperadsWithBasis(R))

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ZinbielOperad(QQ)     # indirect doctest
            The Zinbiel operad over Rational Field
        """
        return f"The Zinbiel operad over {self.base_ring()}"

    def species(self):
        """
        Return the species of non-empty lists.

        This is the species underlying the Zinbiel operad.

        EXAMPLES::

            sage: f = ZinbielOperad(QQ).species()
            sage: f.generating_series()[:5]
            [1, 1, 1, 1]
        """
        from sage.combinat.species.library import LinearOrderSpecies
        return LinearOrderSpecies().restricted(min=1)

    def _coerce_end(self, st):
        """
        Allow for the shortcut ``A(<string>)``.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: A("abc")  # indirect doctest
            B[word: abc]
        """
        if isinstance(st, str):
            return self._from_key(st)
        raise TypeError

    def _from_key(self, k):
        """
        Return an element from a word.

        EXAMPLES::

            sage: ZinbielOperad(QQ)._from_key("abc")
            B[word: abc]
        """
        return self._element_constructor_(self.basis().keys()(k))

    @cached_method
    def one_basis(self, letter='@'):
        """
        Return the word of length one, which index the one of this operad.

        INPUT:

        - ``letter`` (default ``'@'``) -- letter used to label the unit

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def operad_generators(self):
        """
        Return the generators of the operad.

        EXAMPLES::

            sage: ZinbielOperad(QQ).operad_generators()
            Finite family {'zinbiel_product': B[word: 12]}
        """
        from sage.sets.family import Family
        return Family({"zinbiel_product": self._from_key([1, 2])})

    def composition_on_basis_list(self, x, y, i):
        r"""
        Return the composition of two words `x o_i y` as a list of
        words.

        The composition index `i` must be a label of `x`.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: W = A.basis().keys()
            sage: list(A.composition_on_basis_list(W("abc"), W("de"), "a"))
            [word: dbce, word: dbec, word: debc]

        TESTS::

            sage: A.composition_on_basis_list(W("abc"), W("de"), "f")
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
        """
        if i not in x:
            raise ValueError("the composition index is not present")
        elif x[0] == i:
            return (y[:1] + u for u in ShuffleProduct_w1w2(x[1:], y[1:]))
        return (x[:1] + u for u in self.composition_on_basis_list(x[1:], y, i))

    def composition_on_basis(self, x, y, i):
        """
        Return the composition of words as a sum of basis elements.

        INPUT:

        - x, y -- words

        - i -- the composition index

        The composition index `i` must be a label of `x`.

        EXAMPLES::

            sage: A = ZinbielOperad(QQ)
            sage: W = A.basis().keys()
            sage: A.composition_on_basis(W(["a","b","c"]), W(["d","e"]), "a")
            B[word: dbce] + B[word: dbec] + B[word: debc]

        TESTS::

            sage: A.composition_on_basis(W(["a","b","c"]), W(["d","e"]), "u")
            Traceback (most recent call last):
            ...
            ValueError: the composition index is not present
        """
        if i not in x:
            raise ValueError("the composition index is not present")
        return self.sum_of_monomials(t for t in
                                     self.composition_on_basis_list(x, y, i))
