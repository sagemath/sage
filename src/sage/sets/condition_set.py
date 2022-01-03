r"""
Subsets of a Universe Defined by Predicates
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent, Set_generic
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.misc import _stable_uniq
from sage.structure.element import Expression
from .set import Set, Set_base, Set_boolean_operators, Set_add_sub_operators


class ConditionSet(Set_generic, Set_base, Set_boolean_operators, Set_add_sub_operators,
                   UniqueRepresentation):
    r"""
    Set of elements of a universe that satisfy given predicates

    INPUT:

    - ``universe`` -- a set

    - ``*predicates`` -- callables

    - ``vars`` or ``names`` -- (default: inferred from ``predicates`` if any predicate is
      an element of a :class:`~sage.symbolic.callable.CallableSymbolicExpressionRing_class`)
      variables or names of variables

    - ``category`` -- (default: inferred from ``universe``) a category

    EXAMPLES::

        sage: Evens = ConditionSet(ZZ, is_even); Evens
        { x ∈ Integer Ring : <function is_even at 0x...>(x) }
        sage: 2 in Evens
        True
        sage: 3 in Evens
        False
        sage: 2.0 in Evens
        True

        sage: Odds = ConditionSet(ZZ, is_odd); Odds
        { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
        sage: EvensAndOdds = Evens | Odds; EvensAndOdds
        Set-theoretic union of
         { x ∈ Integer Ring : <function is_even at 0x...>(x) } and
         { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
        sage: 5 in EvensAndOdds
        True
        sage: 7/2 in EvensAndOdds
        False

        sage: var('y')
        y
        sage: SmallOdds = ConditionSet(ZZ, is_odd, abs(y) <= 11, vars=[y]); SmallOdds
        { y ∈ Integer Ring : abs(y) <= 11, <function is_odd at 0x...>(y) }

        sage: P = polytopes.cube(); P
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: P.rename("P")
        sage: P_inter_B = ConditionSet(P, lambda x: x.norm() < 1.2); P_inter_B
        { x ∈ P : <function <lambda> at 0x...>(x) }
        sage: vector([1, 0, 0]) in P_inter_B
        True
        sage: vector([1, 1, 1]) in P_inter_B
        False

        sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 1.2; predicate
        (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 1.20000000000000
        sage: P_inter_B_again = ConditionSet(P, predicate); P_inter_B_again
        { (x, y, z) ∈ P : sqrt(x^2 + y^2 + z^2) < 1.20000000000000 }
        sage: vector([1, 0, 0]) in P_inter_B_again
        True
        sage: vector([1, 1, 1]) in P_inter_B_again
        False

    Iterating over subsets determined by predicates::

        sage: Odds = ConditionSet(ZZ, is_odd); Odds
        { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
        sage: list(Odds.iterator_range(stop=6))
        [1, -1, 3, -3, 5, -5]

        sage: R = IntegerModRing(8)
        sage: R_primes = ConditionSet(R, is_prime); R_primes
        { x ∈ Ring of integers modulo 8 : <function is_prime at 0x...>(x) }
        sage: R_primes.is_finite()
        True
        sage: list(R_primes)
        [2, 6]

    Using ``ConditionSet`` without predicates provides a way of attaching variable names
    to a set::

        sage: Z3 = ConditionSet(ZZ^3, vars=['x', 'y', 'z']); Z3
        { (x, y, z) ∈ Ambient free module of rank 3 over the principal ideal domain Integer Ring }
        sage: Z3.variable_names()
        ('x', 'y', 'z')
        sage: Z3.arguments()
        (x, y, z)

        sage: Q4.<a, b, c, d> = ConditionSet(QQ^4); Q4
        { (a, b, c, d) ∈ Vector space of dimension 4 over Rational Field }
        sage: Q4.variable_names()
        ('a', 'b', 'c', 'd')
        sage: Q4.arguments()
        (a, b, c, d)

    TESTS::

        sage: TestSuite(P_inter_B).run(skip='_test_pickling')  # cannot pickle lambdas
        sage: TestSuite(P_inter_B_again).run()
    """
    @staticmethod
    def __classcall_private__(cls, universe, *predicates, vars=None, names=None, category=None):
        r"""
        Normalize init arguments.

        TESTS::

            sage: ConditionSet(ZZ, names=["x"]) is ConditionSet(ZZ, names=x)
            True
            sage: ConditionSet(RR, x > 0, names=x) is ConditionSet(RR, (x > 0).function(x))
            True
        """
        if category is None:
            category = Sets()
        if isinstance(universe, Parent):
            if universe in Sets().Finite():
                category &= Sets().Finite()
            if universe in EnumeratedSets():
                category &= EnumeratedSets()

        if vars is not None:
            if names is not None:
                raise ValueError('cannot use names and vars at the same time; they are aliases')
            names, vars = vars, None

        if names is not None:
            names = normalize_names(-1, names)

        callable_symbolic_predicates = []
        other_predicates = []

        for predicate in predicates:
            if isinstance(predicate, Expression) and predicate.is_callable():
                if names is None:
                    names = tuple(str(var) for var in predicate.args())
                elif len(names) != len(predicate.args()):
                    raise TypeError('mismatch in number of arguments')
                if vars is None:
                    vars = predicate.args()
                callable_symbolic_predicates.append(predicate)
            elif isinstance(predicate, Expression):
                if names is None:
                    raise TypeError('use callable symbolic expressions or provide variable names')
                if vars is None:
                    from sage.symbolic.ring import SR
                    vars = tuple(SR.var(name) for name in names)
                callable_symbolic_predicates.append(predicate.function(*vars))
            else:
                other_predicates.append(predicate)

        predicates = list(_stable_uniq(callable_symbolic_predicates + other_predicates))

        if not other_predicates and not callable_symbolic_predicates:
            if names is None and category is None:
                # No conditions, no variable names, no category, just use Set.
                return Set(universe)

        if any(predicate.args() != vars
               for predicate in callable_symbolic_predicates):
            # TODO: Implement safe renaming of the arguments of a callable symbolic expressions
            raise NotImplementedError('all callable symbolic expressions must use the same arguments')

        if names is None:
            names = ("x",)
        return super().__classcall__(cls, universe, *predicates,
                                     names=names, category=category)

    def __init__(self, universe, *predicates, names=None, category=None):
        r"""
        TESTS::

            sage: Evens = ConditionSet(ZZ, is_even); Evens
            { x ∈ Integer Ring : <function is_even at 0x...>(x) }
            sage: TestSuite(Evens).run()
        """
        self._universe = universe
        self._predicates = predicates
        facade = None
        if isinstance(universe, Parent):
            facade = universe
        super().__init__(facade=facade, category=category,
                         names=names, normalize=False) # names already normalized by classcall

    def _first_ngens(self, n):
        r"""
        Return the list of variables.

        This is useful only for the use of Sage preparser::

            sage: preparse("Q3.<x,y,z> = ConditionSet(QQ^3)")
            "Q3 = ConditionSet(QQ**Integer(3), names=('x', 'y', 'z',)); (x, y, z,) = Q3._first_ngens(3)"

        """
        return self.arguments()

    def _repr_(self):
        """
        Print representation of this set.

        EXAMPLES::

            sage: var('t') # parameter
            t
            sage: ZeroDimButNotNullary = ConditionSet(ZZ^0, t > 0, vars=("q"))
            sage: ZeroDimButNotNullary._repr_()
            '{ q ∈ Ambient free module of rank 0
                   over the principal ideal domain Integer Ring : t > 0 }'
        """
        s = "{ "
        names = self.variable_names()
        comma_sep_names = ", ".join(str(name) for name in names)
        if len(names) == 1:
            s += f"{comma_sep_names}"
        else:
            s += f"({comma_sep_names})"
        universe = self._universe
        s += f" ∈ {universe}"
        sep = " : "
        for predicate in self._predicates:
            s += sep + self._repr_condition(predicate)
            sep = ", "
        s += " }"
        return s

    @cached_method
    def _repr_condition(self, predicate):
        """
        Format the predicate, applied to the arguments.

        EXAMPLES::

            sage: Evens = ConditionSet(ZZ, is_even)
            sage: Evens._repr_condition(is_even)
            '<function is_even at 0x...>(x)'
            sage: BigSin = ConditionSet(RR, sin(x) > 0.9, vars=[x])
            sage: BigSin._repr_condition(BigSin._predicates[0])
            'sin(x) > 0.900000000000000'
            sage: var('t') # parameter
            t
            sage: ZeroDimButNotNullary = ConditionSet(ZZ^0, t > 0, vars=("q"))
            sage: ZeroDimButNotNullary._repr_condition(ZeroDimButNotNullary._predicates[0])
            't > 0'
        """
        if isinstance(predicate, Expression) and predicate.is_callable():
            args = self.arguments()
            if len(args) == 1:
                args = args[0]
            condition = self._call_predicate(predicate, args)
            return str(condition)
        comma_sep_names = ", ".join(str(name)
                                    for name in self.variable_names())
        return f"{predicate}({comma_sep_names})"

    @cached_method
    def arguments(self):
        """
        Return the variables of ``self`` as elements of the symbolic ring.

        EXAMPLES::

            sage: Odds = ConditionSet(ZZ, is_odd); Odds
            { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
            sage: args = Odds.arguments(); args
            (x,)
            sage: args[0].parent()
            Symbolic Ring
        """
        from sage.symbolic.ring import SR
        return SR.var(self.variable_names())

    def _element_constructor_(self, *args, **kwds):
        """
        Construct an element of the set.

        This element constructor raises an error if the element does not
        satisfy the predicates.

        EXAMPLES::

            sage: Evens = ConditionSet(ZZ, is_even); Evens
            { x ∈ Integer Ring : <function is_even at 0x...>(x) }
            sage: element_two = Evens(2r); element_two
            2
            sage: element_two.parent()
            Integer Ring
            sage: element_too = Evens(2.0); element_too
            2
            sage: element_too.parent()
            Integer Ring
            sage: Evens(3)
            Traceback (most recent call last):
            ...
            ValueError: 3 does not satisfy the condition

        """
        try:
            universe_element_constructor = self._universe._element_constructor_
        except AttributeError:
            if len(args) != 1 or kwds:
                raise ValueError('element constructor only takes 1 argument')
            element = args[0]
            if element not in self._universe:
                raise ValueError(f'{element} is not an element of the universe')
        else:
            element = universe_element_constructor(*args, **kwds)
        if not all(self._call_predicate(predicate, element)
                   for predicate in self._predicates):
            raise ValueError(f'{element} does not satisfy the condition')
        return element

    def _call_predicate(self, predicate, element):
        r"""
        Call ``predicate`` on an ``element`` of the universe of ``self``.

        TESTS::

            sage: TripleDigits = ZZ^3
            sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 12; predicate
            (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 12
            sage: SmallTriples = ConditionSet(ZZ^3, predicate); SmallTriples
            { (x, y, z) ∈ Ambient free module of rank 3 over the principal
                          ideal domain Integer Ring : sqrt(x^2 + y^2 + z^2) < 12 }
            sage: predicate = SmallTriples._predicates[0]
            sage: element = TripleDigits((1, 2, 3))
            sage: SmallTriples._call_predicate(predicate, element)
            sqrt(14) < 12

            sage: var('t')
            t
            sage: TinyUniverse = ZZ^0
            sage: Nullary = ConditionSet(TinyUniverse, t > 0, vars=())
            sage: predicate = Nullary._predicates[0]
            sage: element = TinyUniverse(0)
            sage: Nullary._call_predicate(predicate, element)
            t > 0
        """
        if isinstance(predicate, Expression) and predicate.is_callable():
            if len(predicate.arguments()) != 1:
                return predicate(*element)
        return predicate(element)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        This may raise ``NotImplementedError``.

        TESTS::

            sage: TripleDigits = ZZ^3
            sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 12; predicate
            (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 12
            sage: SmallTriples = ConditionSet(ZZ^3, predicate); SmallTriples
            { (x, y, z) ∈ Ambient free module of rank 3 over the principal
                          ideal domain Integer Ring : sqrt(x^2 + y^2 + z^2) < 12 }
            sage: SmallTriples.an_element()  # indirect doctest
            (1, 0, 0)
        """
        for element in self._universe.some_elements():
            if element in self:
                return element
        raise NotImplementedError

    def ambient(self):
        r"""
        Return the universe of ``self``.

        EXAMPLES::

            sage: Evens = ConditionSet(ZZ, is_even); Evens
            { x ∈ Integer Ring : <function is_even at 0x...>(x) }
            sage: Evens.ambient()
            Integer Ring
        """
        return self._universe

    @cached_method
    def _sympy_(self):
        r"""
        Return an instance of a subclass of SymPy ``Set`` corresponding to ``self``.

        EXAMPLES::

            sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 12; predicate
            (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 12
            sage: SmallTriples = ConditionSet(ZZ^3, predicate); SmallTriples
            { (x, y, z) ∈ Ambient free module of rank 3 over the principal
                          ideal domain Integer Ring : sqrt(x^2 + y^2 + z^2) < 12 }
            sage: ST = SmallTriples._sympy_(); ST
            ConditionSet((x, y, z), sqrt(x**2 + y**2 + z**2) < 12,
                         ProductSet(Integers, Integers, Integers))
            sage: (1, 3, 5) in ST
            True
            sage: (5, 7, 9) in ST
            False

            sage: Interval = ConditionSet(RR, x >= -7, x <= 4, vars=[x]); Interval
            { x ∈ Real Field with 53 bits of precision : x >= -7, x <= 4 }
            sage: Interval._sympy_()
            ConditionSet(x, (x >= -7) & (x <= 4), SageSet(Real Field with 53 bits of precision))

        If a predicate is not symbolic, we fall back to creating a wrapper::

            sage: Evens = ConditionSet(ZZ, is_even); Evens
            { x ∈ Integer Ring : <function is_even at 0x...>(x) }
            sage: Evens._sympy_()
            SageSet({ x ∈ Integer Ring : <function is_even at 0x...>(x) })
        """
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        import sympy

        args = self.arguments()
        single_arg = len(args) == 1
        if single_arg:
            args = args[0]

        try:
            conditions = [self._call_predicate(predicate, args)
                          for predicate in self._predicates]

            sym = tuple(x._sympy_() for x in self.arguments())
            if single_arg:
                sym = sym[0]
            result = sympy.ConditionSet(sym,
                                        sympy.And(*[condition._sympy_()
                                                    for condition in conditions]),
                                        base_set=self._universe._sympy_())
            result._sage_object = self
            return result
        except TypeError:
            # Fall back to creating a wrapper
            return super()._sympy_()

    def intersection(self, X):
        r"""
        Return the intersection of ``self`` and ``X``.

        EXAMPLES::

            sage: in_small_oblong(x, y) = x^2 + 3 * y^2 <= 42
            sage: SmallOblongUniverse = ConditionSet(QQ^2, in_small_oblong)
            sage: SmallOblongUniverse
            { (x, y) ∈ Vector space of dimension 2 over Rational Field : x^2 + 3*y^2 <= 42 }
            sage: parity_check(x, y) = abs(sin(pi/2*(x + y))) < 1/1000
            sage: EvenUniverse = ConditionSet(ZZ^2, parity_check); EvenUniverse
            { (x, y) ∈ Ambient free module of rank 2 over the principal ideal
                       domain Integer Ring : abs(sin(1/2*pi*x + 1/2*pi*y)) < (1/1000) }
            sage: SmallOblongUniverse & EvenUniverse
            { (x, y) ∈ Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 1] : x^2 + 3*y^2 <= 42, abs(sin(1/2*pi*x + 1/2*pi*y)) < (1/1000) }

        Combining two ``ConditionSet``s with different formal variables works correctly.
        The formal variables of the intersection are taken from ``self``::

            sage: SmallMirrorUniverse = ConditionSet(QQ^2, in_small_oblong, vars=(y, x))
            sage: SmallMirrorUniverse
            { (y, x) ∈ Vector space of dimension 2 over Rational Field : 3*x^2 + y^2 <= 42 }
            sage: SmallOblongUniverse & SmallMirrorUniverse
            { (x, y) ∈ Vector space of dimension 2 over Rational Field : x^2 + 3*y^2 <= 42 }
            sage: SmallMirrorUniverse & SmallOblongUniverse
            { (y, x) ∈ Vector space of dimension 2 over Rational Field : 3*x^2 + y^2 <= 42 }
        """
        if isinstance(X, ConditionSet):
            return ConditionSet(self.ambient().intersection(X.ambient()),
                                *(self._predicates + X._predicates),
                                vars=self.arguments())
        return super().intersection(X)

    def __iter__(self):
        r"""
        Iterate over ``self``.

        TESTS::

            sage: Odds = ConditionSet(ZZ, is_odd); Odds
            { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
            sage: list(Odds.iterator_range(stop=6))
            [1, -1, 3, -3, 5, -5]

        """
        for x in self._universe:
            if x in self:
                yield x
