# sage.doctest: needs sage.combinat
r"""
Monoids
"""

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method


class Monoid_class(Parent):
    def __init__(self, names, category=None) -> None:
        r"""
        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid_class
            sage: Monoid_class(('a','b'))
            <sage.monoids.monoid.Monoid_class_with_category object at ...>

        TESTS::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: TestSuite(F).run()
        """
        from sage.categories.monoids import Monoids
        if category is None:
            cat = Monoids().FinitelyGeneratedAsMagma()
        else:
            cat = category & Monoids().FinitelyGeneratedAsMagma()
        Parent.__init__(self, base=self, names=names, category=cat)

    @cached_method
    def gens(self) -> tuple:
        r"""
        Return the generators for ``self``.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.gens()
            (a, b, c, d, e)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))

    def monoid_generators(self):
        r"""
        Return the generators for ``self``.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.monoid_generators()
            Family (a, b, c, d, e)
        """
        from sage.sets.family import Family
        return Family(self.gens())
