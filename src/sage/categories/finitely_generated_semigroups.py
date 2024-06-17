# sage_setup: distribution = sagemath-categories
r"""
Finitely generated semigroups
"""
# ****************************************************************************
#  Copyright (C) 2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import itertools
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.semigroups import Semigroups
from sage.categories.enumerated_sets import EnumeratedSets


class FinitelyGeneratedSemigroups(CategoryWithAxiom):
    r"""
    The category of finitely generated (multiplicative) semigroups.

    A :class:`finitely generated semigroup <Semigroups>` is a
    :class:`semigroup <Semigroups>` endowed with a distinguished
    finite set of generators (see
    :meth:`FinitelyGeneratedSemigroups.ParentMethods.semigroup_generators`).
    This makes it into an :class:`enumerated set <EnumeratedSets>`.

    EXAMPLES::

        sage: C = Semigroups().FinitelyGenerated(); C
        Category of finitely generated semigroups
        sage: C.super_categories()
        [Category of semigroups,
         Category of finitely generated magmas,
         Category of enumerated sets]
        sage: sorted(C.axioms())
        ['Associative', 'Enumerated', 'FinitelyGeneratedAsMagma']
        sage: C.example()
        An example of a semigroup: the free semigroup generated
        by ('a', 'b', 'c', 'd')

    TESTS::

        sage: TestSuite(C).run()
    """

    _base_category_class_and_axiom = (Semigroups, "FinitelyGeneratedAsMagma")

    @cached_method
    def extra_super_categories(self):
        r"""
        State that a finitely generated semigroup is endowed with a
        default enumeration.

        EXAMPLES::

            sage: Semigroups().FinitelyGenerated().extra_super_categories()
            [Category of enumerated sets]
        """
        return [EnumeratedSets()]

    def example(self):
        r"""
        EXAMPLES::

            sage: Semigroups().FinitelyGenerated().example()
            An example of a semigroup: the free semigroup generated
            by ('a', 'b', 'c', 'd')
        """
        return Semigroups().example("free")

    class ParentMethods:

        @abstract_method
        def semigroup_generators(self):
            r"""
            Return distinguished semigroup generators for ``self``.

            OUTPUT: a finite family

            This method should be implemented by all semigroups in
            :class:`FinitelyGeneratedSemigroups`.

            EXAMPLES::

                sage: S = FiniteSemigroups().example()
                sage: S.semigroup_generators()
                Family ('a', 'b', 'c', 'd')
            """

        # TODO: update transitive ideal

        def succ_generators(self, side='twosided'):
            r"""
            Return the successor function of the ``side``-sided Cayley
            graph of ``self``.

            This is a function that maps an element of ``self`` to all
            the products of ``x`` by a generator of this semigroup,
            where the product is taken on the left, right, or both
            sides.

            INPUT:

            - ``side`` -- ``'left'``, ``'right'``, or ``'twosided'`` (default)

            .. TODO:: Design choice:

               - find a better name for this method
               - should we return a set? a family?

            EXAMPLES::

                sage: S = FiniteSemigroups().example()
                sage: S.succ_generators("left" )(S('ca'))
                ('ac', 'bca', 'ca', 'dca')
                sage: S.succ_generators("right")(S('ca'))
                ('ca', 'cab', 'ca', 'cad')
                sage: S.succ_generators("twosided" )(S('ca'))
                ('ac', 'bca', 'ca', 'dca', 'ca', 'cab', 'ca', 'cad')
            """
            left = (side == "left" or side == "twosided")
            right = (side == "right" or side == "twosided")
            generators = self.semigroup_generators()
            return lambda x: (tuple(g * x for g in generators) if left else ()) + (tuple(x * g for g in generators) if right else ())

        def __iter__(self):
            r"""
            Return an iterator over the elements of ``self``.

            This brute force implementation recursively multiplies
            together the distinguished semigroup generators.

            .. SEEALSO:: :meth:`semigroup_generators`

            EXAMPLES::

                sage: S = FiniteSemigroups().example(alphabet=('x','y'))
                sage: it = S.__iter__()
                sage: sorted(it)
                ['x', 'xy', 'y', 'yx']
            """
            from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
            return iter(RecursivelyEnumeratedSet(self.semigroup_generators(),
                                                 self.succ_generators(side='right'),
                                                 enumeration='breadth'))

        def ideal(self, gens, side='twosided'):
            r"""
            Return the ``side``-sided ideal generated by ``gens``.

            This brute force implementation recursively multiplies the
            elements of ``gens`` by the distinguished generators of
            this semigroup.

            .. SEEALSO:: :meth:`semigroup_generators`

            INPUT:

            - ``gens`` -- list (or iterable) of elements of ``self``
            - ``side`` -- [default: "twosided"] "left", "right" or "twosided"

            EXAMPLES::

                sage: S = FiniteSemigroups().example()
                sage: sorted(S.ideal([S('cab')], side='left'))
                ['abc', 'abcd', 'abdc', 'acb', 'acbd', 'acdb', 'adbc',
                 'adcb', 'bac', 'bacd', 'badc', 'bca', 'bcad', 'bcda',
                 'bdac', 'bdca', 'cab', 'cabd', 'cadb', 'cba', 'cbad',
                 'cbda', 'cdab', 'cdba', 'dabc', 'dacb', 'dbac', 'dbca',
                 'dcab', 'dcba']
                sage: list(S.ideal([S('cab')], side='right'))
                ['cab', 'cabd']
                sage: sorted(S.ideal([S('cab')], side='twosided'))
                ['abc', 'abcd', 'abdc', 'acb', 'acbd', 'acdb', 'adbc',
                 'adcb', 'bac', 'bacd', 'badc', 'bca', 'bcad', 'bcda',
                 'bdac', 'bdca', 'cab', 'cabd', 'cadb', 'cba', 'cbad',
                 'cbda', 'cdab', 'cdba', 'dabc', 'dacb', 'dbac', 'dbca',
                 'dcab', 'dcba']
                sage: sorted(S.ideal([S('cab')]))
                ['abc', 'abcd', 'abdc', 'acb', 'acbd', 'acdb', 'adbc',
                 'adcb', 'bac', 'bacd', 'badc', 'bca', 'bcad', 'bcda',
                 'bdac', 'bdca', 'cab', 'cabd', 'cadb', 'cba', 'cbad',
                 'cbda', 'cdab', 'cdba', 'dabc', 'dacb', 'dbac', 'dbca',
                 'dcab', 'dcba']
            """
            from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
            return RecursivelyEnumeratedSet(gens,
                                            self.succ_generators(side=side))

    class Finite(CategoryWithAxiom):

        class ParentMethods:
            def some_elements(self):
                r"""
                Return an iterable containing some elements of the semigroup.

                OUTPUT: the ten first elements of the semigroup, if they exist

                EXAMPLES::

                    sage: S = FiniteSemigroups().example(alphabet=('x','y'))
                    sage: sorted(S.some_elements())
                    ['x', 'xy', 'y', 'yx']
                    sage: S = FiniteSemigroups().example(alphabet=('x','y','z'))
                    sage: X = S.some_elements()
                    sage: len(X)
                    10
                    sage: all(x in S for x in X)
                    True
                """
                return list(itertools.islice(self, 10))
