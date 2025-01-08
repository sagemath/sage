"""
The Commutative Operad
"""
from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method


class CommutativeOperad(CombinatorialFreeModule):
    r"""
    The Commutative operad

    This is an operad on the species of non-empty sets.

    EXAMPLES::

        sage: A = CommutativeOperad(QQ)
        sage: W = A.basis().keys()
        sage: x = A(W('ab'))
        sage: y = A(W('dc'))
        sage: x.compose(y, 'a')
        B[word: bcd]

    REFERENCES:

    .. [todo]
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = CommutativeOperad(QQ); A
            The Commutative operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category=OperadsWithBasis(R))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CommutativeOperad(QQ)       # indirect doctest
            The Commutative operad over Rational Field
        """
        return f"The Commutative operad over {self.base_ring()}"

    def species(self):
        """
        Return the species of non-empty sets.

        This is the species underlying the Commutative operad.

        EXAMPLES::

            sage: f = CommutativeOperad(QQ).species()
            sage: f.generating_series()[:5]
            [1, 1/2, 1/6, 1/24]
        """
        from sage.combinat.species.library import SetSpecies
        return SetSpecies().restricted(min=1)

    def _coerce_end(self, st):
        """
        Allow for the shortcut ``A(<string>)``.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
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

            sage: CommutativeOperad(QQ)._from_key("abc")
            B[word: abc]
            sage: CommutativeOperad(QQ)._from_key("bac")
            B[word: abc]
        """
        return self._element_constructor(self.basis().keys()(sorted(k)))

    @cached_method
    def one_basis(self, letter='@'):
        """
        Return the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        INPUT:

        - ``letter`` (default ``'@'``) -- letter used to label the unit

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def degree_on_basis(self, t):
        """
        Return the degree of a word `t` in the Commutative operad.

        This is the length of the word.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.degree_on_basis(m)
            4
        """
        return t.length()

    def map_labels(self, t, f):
        """
        Map the function `f` on the word `t`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.map_labels(m,lambda u:u)
            word: 1234
        """
        return self.basis().keys()(sorted([f(u) for u in t]))

    def labelling_on_basis(self, t):
        """
        Put canonical labels on a word in the Commutative operad.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.labelling_on_basis(m)
            B[word: 1234]
        """
        B = self.basis()
        return B[B.keys()(sorted([1 + i for i in range(t.length())]))]

    def unlabelling_on_basis(self, t):
        """
        Remove the labels of a word in the Commutative operad.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.unlabelling_on_basis(m)
            B[word: 1111]
        """
        B = self.basis()
        return B[B.keys()([1 for i in range(t.length())])]

    def grafts(self, x, y, i):
        r"""
        Return the word obtained by inserting a word `y` at position `i`
        in a word `x`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.grafts(Words("acb"), Words("de"),"c")
            word: abde
        """
        if x[0] == i:
            return self.basis().keys()(sorted(y + x[1:]))
        return self.basis().keys()(sorted(x[:1] + self.grafts(x[1:], y, i)))

    def composition_on_basis(self, x, y, i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words("acb"), Words("de"),"c")
            B[word: abde]
        """
        if i not in x:
            raise ValueError("the composition index is not present")

        return self.basis()[self.grafts(x, y, i)]

    def commutative_product(self, x, y):
        """
        Compute the commutative product.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
            sage: W = A.basis().keys()
            sage: x = A(W('ab'))
            sage: y = A(W('dc'))
            sage: A.commutative_product(x, y)
            B[word: abcd]
        """
        gen = self.basis()[self.basis().keys()([0, 1])]
        return gen.compose(x, 0).compose(y, 1)

    def operad_generators(self):
        """
        Return the generators of the operad.

        EXAMPLES::

            sage: CommutativeOperad(QQ).operad_generators()
            Finite family {'commutative_product': B[word: 12]}
        """
        from sage.sets.family import Family
        return Family({"commutative_product": self.basis()[self.basis().keys()([1, 2])]})

    def operad_morphism_on_basis(self, t, codomain):
        """
        Define a morphism from the Commutative operad to the target operad

        The target operad has to possess a method called
        ``commutative_product``.

        The argument `t` should not have repeated labels.

        EXAMPLES::

            sage: A = CommutativeOperad(QQ)
        """
        targetProduct = codomain.commutative_product
        n = len(t)
        if n == 1:
            return codomain.one(t[0])
        return targetProduct(self.operad_morphism_on_basis(t[0], codomain),
                             self.operad_morphism_on_basis(t[1:], codomain))
