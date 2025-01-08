from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method


class AssociativeOperad(CombinatorialFreeModule):
    r"""
    The Associative operad
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = AssociativeOperad(QQ); A
            The Associative operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category=OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: AssociativeOperad(QQ)       # indirect doctest
            The Associative operad over Rational Field
        """
        return f"The Associative operad over {self.base_ring()}"

    def species(self):
        """
        The species of non-empty lists

        EXAMPLES::

            sage: f = AssociativeOperad(QQ).species()
            sage: f.generating_series()[:5]
            [1, 1, 1, 1]
        """
        from sage.combinat.species.linear_order_species import LinearOrderSpecies
        return LinearOrderSpecies().restricted(min=1)

    @cached_method
    def one_basis(self, letter):
        """
        Return the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def degree_on_basis(self, t):
        """
        Return the degree of a word `t` in the Associative operad.

        This is the length of the word.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.degree_on_basis(m)
            4
        """
        return t.length()

    def map_labels(self, t, f):
        """
        Map the function `f` on the word `t`.

        INPUT:

        - t -- the index of a basis element

        - f -- a map that can be applied to the labels of t

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.map_labels(m,lambda u:u)
            word: 4321
        """
        return self.basis().keys()([f(u) for u in t])

    def labelling_on_basis(self, t):
        """
        Put canonical labels on a word in the Associative operad.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.labelling_on_basis(m)
            B[word: 1234]
        """
        B = self.basis()
        return B[B.keys()([1 + i for i in range(t.length())])]

    def unlabelling_on_basis(self, t):
        """
        Removes the labels of a tree in the Associative operad.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.unlabelling_on_basis(m)
            B[word: 1111]
        """
        B = self.basis()
        return B[B.keys()([1 for i in range(t.length())])]

    def grafts(self, x, y, i):
        """
        Insert a word y at position i in a word x and return a word

        This is the composition of the set-theoretic Associative operad.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.grafts(Words("acb"), Words("de"),"c")
            word: adeb
        """
        if x[0] == i:
            return y + x[1:]
        return x[:1] + self.grafts(x[1:], y, i)

    def composition_on_basis(self, x, y, i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words("acb"), Words("de"),"c")
            B[word: adeb]
        """
        if i not in x:
            raise ValueError("the composition index is not present")
        return self.basis()[self.grafts(x, y, i)]

    def associative_product(self, x, y):
        """
        Return the associative product of ``x`` and ``y``.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: W = A.basis().keys()
            sage: x = A(W('ab'))
            sage: y = A(W('dc'))
            sage: A.associative_product(x, y)
            B[word: abdc]
        """
        gen = self.basis()[self.basis().keys()([0, 1])]
        return gen.compose(x, 0).compose(y, 1)

    def operad_generators(self):
        """
        EXAMPLES::

            sage: AssociativeOperad(QQ).operad_generators()
            Finite family {'associative_product': B[word: 12]}
        """
        from sage.sets.family import Family
        return Family({"associative_product": self.basis()[self.basis().keys()([1, 2])]})

    def operad_morphism_on_basis(self, t, codomain):
        """
        Defines a morphism from the Associative operad to the target operad

        The target operad has to possess a method called ``associative_product``.

        The argument should not have repeated labels.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: D = DendriformOperad(QQ)
            sage: A.operad_morphism_on_basis(A.one_basis('a'), D)
            B[a[., .]]
        """
        targetProduct = codomain.associative_product
        n = len(t)
        if n == 1:
            return codomain.one(t[0])
        return targetProduct(self.operad_morphism_on_basis(t[0], codomain),
                             self.operad_morphism_on_basis(t[1:], codomain))
