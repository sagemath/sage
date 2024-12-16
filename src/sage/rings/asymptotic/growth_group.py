r"""
(Asymptotic) Growth Groups

This module provides support for (asymptotic) growth groups.

Such groups are equipped with a partial order: the elements can be
seen as functions, and the behavior as their argument (or arguments)
gets large (tend to `\infty`) is compared.

Growth groups are used for the calculations done in the
:doc:`asymptotic ring <asymptotic_ring>`. There, take a look at the
:ref:`informal definition <asymptotic_ring_definition>`, where
examples of growth groups and elements are given as well.


.. _growth_group_description:

Description of Growth Groups
============================

Many growth groups can be described by a string, which can also be used to
create them. For example, the string ``'x^QQ * log(x)^ZZ * QQ^y * y^QQ'``
represents a growth group with the following properties:

- It is a growth group in the two variables `x` and `y`.

- Its elements are of the form

  .. MATH::

      x^r \cdot \log(x)^s \cdot a^y \cdot y^q

  for `r\in\QQ`, `s\in\ZZ`, `a\in\QQ` and `q\in\QQ`.

- The order is with respect to `x\to\infty` and `y\to\infty` independently
  of each other.

- To compare such elements, they are split into parts belonging to
  only one variable. In the example above,

  .. MATH::

      x^{r_1} \cdot \log(x)^{s_1} \leq x^{r_2} \cdot \log(x)^{s_2}

  if `(r_1, s_1) \leq (r_2, s_2)` lexicographically. This reflects the fact
  that elements `x^r` are larger than elements `\log(x)^s` as `x\to\infty`.
  The factors belonging to the variable `y` are compared analogously.

  The results of these comparisons are then put together using the
  :wikipedia:`product order <Product_order>`, i.e., `\leq` if each component
  satisfies `\leq`.


Each description string consists of ordered factors---yes, this means
``*`` is noncommutative---of strings describing "elementary" growth
groups (see the examples below). As stated in the example above, these
factors are split by their variable; factors with the same variable are
grouped. Reading such factors from left to right determines the order:
Comparing elements of two factors (growth groups) `L` and `R`, then all
elements of `L` are considered to be larger than each element of `R`.


.. _growth_group_creating:

Creating a Growth Group
=======================

For many purposes the factory ``GrowthGroup`` (see
:class:`GrowthGroupFactory`) is the most convenient way to generate a
growth group.
::

    sage: from sage.rings.asymptotic.growth_group import GrowthGroup

Here are some examples::

    sage: GrowthGroup('z^ZZ')
    Growth Group z^ZZ
    sage: M = GrowthGroup('z^QQ'); M
    Growth Group z^QQ

Each of these two generated groups is a :class:`MonomialGrowthGroup`,
whose elements are powers of a fixed symbol (above ``'z'``).
For the order of the elements it is assumed that `z\to\infty`.

.. NOTE::

    Growth groups where the variable tend to some value distinct from
    `\infty` are not yet implemented.

To create elements of `M`, a generator can be used::

    sage: z = M.gen()
    sage: z^(3/5)
    z^(3/5)

Strings can also be parsed::

    sage: M('z^7')
    z^7

Similarly, we can construct logarithmic factors by::

    sage: GrowthGroup('log(z)^QQ')
    Growth Group log(z)^QQ

which again creates a
:class:`MonomialGrowthGroup`. An :class:`ExponentialGrowthGroup` is generated in the same way. Our factory gives
::

    sage: E = GrowthGroup('(QQ_+)^z'); E
    Growth Group QQ^z

and a typical element looks like this::

    sage: E.an_element()
    (1/2)^z

More complex groups are created in a similar fashion. For example
::

    sage: C = GrowthGroup('(QQ_+)^z * z^QQ * log(z)^QQ'); C
    Growth Group QQ^z * z^QQ * log(z)^QQ

This contains elements of the form
::

    sage: C.an_element()
    (1/2)^z*z^(1/2)*log(z)^(1/2)

The group `C` itself is a Cartesian product; to be precise a
:class:`~sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct`. We
can see its factors::

    sage: C.cartesian_factors()
    (Growth Group QQ^z, Growth Group z^QQ, Growth Group log(z)^QQ)

Multivariate constructions are also possible::

    sage: GrowthGroup('x^QQ * y^QQ')
    Growth Group x^QQ * y^QQ

This gives a
:class:`~sage.rings.asymptotic.growth_group_cartesian.MultivariateProduct`.

Both these Cartesian products are derived from the class
:class:`~sage.rings.asymptotic.growth_group_cartesian.GenericProduct`. Moreover
all growth groups have the abstract base class
:class:`GenericGrowthGroup` in common.

Some Examples
^^^^^^^^^^^^^

::

    sage: from sage.rings.asymptotic.growth_group import GrowthGroup
    sage: G_x = GrowthGroup('x^ZZ'); G_x
    Growth Group x^ZZ
    sage: G_xy = GrowthGroup('x^ZZ * y^ZZ'); G_xy
    Growth Group x^ZZ * y^ZZ
    sage: G_xy.an_element()
    x*y
    sage: x = G_xy('x'); y = G_xy('y')
    sage: x^2
    x^2
    sage: elem = x^21*y^21; elem^2
    x^42*y^42

A monomial growth group itself is totally ordered, all elements
are comparable. However, this does **not** hold for Cartesian
products::

    sage: e1 = x^2*y; e2 = x*y^2
    sage: e1 <= e2 or e2 <= e1
    False

In terms of uniqueness, we have the following behaviour::

    sage: GrowthGroup('x^ZZ * y^ZZ') is GrowthGroup('y^ZZ * x^ZZ')
    True

The above is ``True`` since the order of the factors does not play a role here; they use different variables. But when using the same variable, it plays a role::

    sage: GrowthGroup('x^ZZ * log(x)^ZZ') is GrowthGroup('log(x)^ZZ * x^ZZ')
    False

In this case the components are ordered lexicographically, which
means that in the second growth group, ``log(x)`` is assumed to
grow faster than ``x`` (which is nonsense, mathematically). See
:class:`CartesianProduct <sage.rings.asymptotic.growth_group_cartesian.CartesianProductFactory>`
for more details or see :ref:`above <growth_group_description>`
for a more extensive description.

Short notation also allows the construction of more complicated
growth groups::

    sage: G = GrowthGroup('(QQ_+)^x * x^ZZ * log(x)^QQ * y^QQ')
    sage: G.an_element()
    (1/2)^x*x*log(x)^(1/2)*y^(1/2)
    sage: x, y = var('x y')
    sage: G(2^x * log(x) * y^(1/2)) * G(x^(-5) * 5^x * y^(1/3))
    10^x*x^(-5)*log(x)*y^(5/6)

AUTHORS:

- Benjamin Hackl (2015)
- Daniel Krenn (2015)
- Clemens Heuberger (2016)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from collections import namedtuple

from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.asymptotic.growth_group_cartesian', 'CartesianProductGrowthGroups')

from sage.categories.pushout import ConstructionFunctor
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import (CachedRepresentation,
                                                  UniqueRepresentation)
from sage.structure.richcmp import richcmp_by_eq_and_lt
import sage.rings.abc
from .misc import WithLocals


class Variable(CachedRepresentation, SageObject):
    r"""
    A class managing the variable of a growth group.

    INPUT:

    - ``var`` -- an object whose representation string is used as the
      variable. It has to be a valid Python identifier. ``var`` can
      also be a tuple (or other iterable) of such objects.

    - ``repr`` -- (default: ``None``) if specified, then this string
      will be displayed instead of ``var``. Use this to get
      e.g. ``log(x)^ZZ``: ``var`` is then used to specify the variable `x`.

    - ``latex_name`` -- (default: ``None``) if specified, then this string
      will be used as LaTeX-representation of ``var``

    - ``ignore`` -- (default: ``None``) a tuple (or other iterable)
      of strings which are not variables

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import Variable
        sage: v = Variable('x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x1'); repr(v), v.variable_names()
        ('x1', ('x1',))
        sage: v = Variable('x_42'); repr(v), v.variable_names()
        ('x_42', ('x_42',))
        sage: v = Variable(' x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x '); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable(''); repr(v), v.variable_names()
        ('', ())

    ::

        sage: v = Variable(('x', 'y')); repr(v), v.variable_names()
        ('x, y', ('x', 'y'))
        sage: v = Variable(('x', 'log(y)')); repr(v), v.variable_names()
        ('x, log(y)', ('x', 'y'))
        sage: v = Variable(('x', 'log(x)')); repr(v), v.variable_names()
        Traceback (most recent call last):
        ...
        ValueError: Variable names ('x', 'x') are not pairwise distinct.

    ::

        sage: v = Variable('log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))
        sage: v = Variable('log(log(x))'); repr(v), v.variable_names()
        ('log(log(x))', ('x',))

    ::

        sage: v = Variable('x', repr='log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))

    ::

        sage: v = Variable('e^x', ignore=('e',)); repr(v), v.variable_names()
        ('e^x', ('x',))

    ::

        sage: v = Variable('(e^n)', ignore=('e',)); repr(v), v.variable_names()
        ('e^n', ('n',))
        sage: v = Variable('(e^(n*log(n)))', ignore=('e',)); repr(v), v.variable_names()
        ('e^(n*log(n))', ('n',))
    """
    def __init__(self, var, repr=None, latex_name=None, ignore=None):
        r"""
        See :class:`Variable` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')
            blub
            sage: Variable('blub') is Variable('blub')
            True

        ::

            sage: Variable('(:-)')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: : !!! -
            sage: Variable('(:-)', repr='icecream')
            Traceback (most recent call last):
            ...
            ValueError: ':-' is not a valid name for a variable.

        Check :issue:`26452`::

            sage: Variable(('w',),
            ....:          repr='w^(Number Field in i with defining polynomial x^2 + 1) * log(w)^ZZ')
            w^(Number Field in i with defining polynomial x^2 + 1) * log(w)^ZZ
        """
        from sage.symbolic.ring import isidentifier
        from .misc import split_str_by_op

        if not isinstance(var, (list, tuple)):
            var = (var,)
        var = tuple(''.join(split_str_by_op(str(v), None)) for v in var)  # we strip off parentheses

        if ignore is None:
            ignore = tuple()

        from sage.misc.latex import latex
        from sage.symbolic.ring import SR

        if repr is None:
            var_bases = tuple(i for i in sum(iter(
                self.extract_variable_names(v)
                if not isidentifier(v) else (v,)
                for v in var), tuple()) if i not in ignore)
            var_repr = ', '.join(var)
            if latex_name is None:
                latex_name = ', '.join(latex(SR(v)) for v in var if v)
        else:
            for v in var:
                if not isidentifier(v):
                    raise ValueError("'%s' is not a valid name for a variable." % (v,))
            var_bases = var
            var_repr = str(repr).strip()
            if latex_name is None:
                try:
                    latex_name = latex(SR(var_repr))
                except (TypeError, ValueError):
                    latex_name = latex(var_repr)

        if len(var_bases) != len(set(var_bases)):
            raise ValueError('Variable names %s are not pairwise distinct.' %
                             (var_bases,))

        self.var_bases = var_bases
        self.var_repr = var_repr

        self.latex_name = latex_name

    def __hash__(self):
        r"""
        Return the hash of this variable.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: hash(Variable('blub'))  # random
            -123456789
        """
        return hash((self.var_repr,) + self.var_bases)

    def __eq__(self, other):
        r"""
        Compare whether this variable equals ``other``.

        INPUT:

        - ``other`` -- another variable

        OUTPUT: boolean

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') == Variable('x')
            True
            sage: Variable('x') == Variable('y')
            False
        """
        return self.var_repr == other.var_repr and self.var_bases == other.var_bases

    def __ne__(self, other):
        r"""
        Return whether this variable does not equal ``other``.

        INPUT:

        - ``other`` -- another variable

        OUTPUT: boolean

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') != Variable('x')
            False
            sage: Variable('x') != Variable('y')
            True
        """
        return not self == other

    def _repr_(self):
        r"""
        Return a representation string of this variable.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')  # indirect doctest
            blub
        """
        return self.var_repr

    def _latex_(self):
        r"""
        Return a LaTeX-representation string of this variable.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: latex(Variable('x'))  # indirect doctest
            x
            sage: latex(Variable('x1'))  # indirect doctest
            x_{1}
            sage: latex(Variable('x_1'))  # indirect doctest
            x_{1}
        """
        return self.latex_name

    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').variable_names()
            ('x',)
            sage: Variable('log(x)').variable_names()
            ('x',)
        """
        return self.var_bases

    def is_monomial(self):
        r"""
        Return whether this is a monomial variable.

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').is_monomial()
            True
            sage: Variable('log(x)').is_monomial()
            False
        """
        return len(self.var_bases) == 1 and self.var_bases[0] == self.var_repr

    @staticmethod
    def extract_variable_names(s):
        r"""
        Determine the name of the variable for the given string.

        INPUT:

        - ``s`` -- string

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable.extract_variable_names('')
            ()
            sage: Variable.extract_variable_names('x')
            ('x',)
            sage: Variable.extract_variable_names('exp(x)')
            ('x',)
            sage: Variable.extract_variable_names('sin(cos(ln(x)))')
            ('x',)

        ::

            sage: Variable.extract_variable_names('log(77w)')
            ('w',)
            sage: Variable.extract_variable_names('log(x')
            Traceback (most recent call last):
            ...
            TypeError: Bad function call: log(x !!!
            sage: Variable.extract_variable_names('x)')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: x) !!!
            sage: Variable.extract_variable_names('log)x(')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: log) !!! x(
            sage: Variable.extract_variable_names('log(x)+y')
            ('x', 'y')
            sage: Variable.extract_variable_names('icecream(summer)')
            ('summer',)

        ::

            sage: Variable.extract_variable_names('a + b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a+b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a +b')
            ('a', 'b')
            sage: Variable.extract_variable_names('+a')
            ('a',)
            sage: Variable.extract_variable_names('a+')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: a+ !!!
            sage: Variable.extract_variable_names('b!')
            ('b',)
            sage: Variable.extract_variable_names('-a')
            ('a',)
            sage: Variable.extract_variable_names('a*b')
            ('a', 'b')
            sage: Variable.extract_variable_names('2^q')
            ('q',)
            sage: Variable.extract_variable_names('77')
            ()

        ::

            sage: Variable.extract_variable_names('a + (b + c) + d')
            ('a', 'b', 'c', 'd')
        """
        from sage.symbolic.ring import SR
        if s == '':
            return ()
        return tuple(str(s) for s in SR(s).variables())

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this variable.

        INPUT:

        - ``rules`` -- dictionary

        OUTPUT: an object

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x^2')._substitute_({'x': SR.var('z')})
            z^2
            sage: _.parent()
            Symbolic Ring

        ::

            sage: Variable('1/x')._substitute_({'x': 'z'})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute in 1/x in
            <class 'sage.rings.asymptotic.growth_group.Variable'>.
            > *previous* TypeError: unsupported operand parent(s) for /:
            'Integer Ring' and '<class 'str'>'
            sage: Variable('1/x')._substitute_({'x': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in 1/x in
            <class 'sage.rings.asymptotic.growth_group.Variable'>.
            > *previous* ZeroDivisionError: rational division by zero
        """
        from sage.misc.sage_eval import sage_eval
        try:
            return sage_eval(self.var_repr, locals=rules)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)


class PartialConversionValueError(ValueError):
    r"""
    A special :python:`ValueError<library/exceptions.html#exceptions.ValueError>`
    which is raised when (partial) conversion fails.

    INPUT:

    - ``element`` -- a :class:`PartialConversionElement`

    The remaining argument passed on to
    :python:`ValueError<library/exceptions.html#exceptions.ValueError>`.
    """
    def __init__(self, element, *args, **kwds):
        r"""
        See :exc:`PartialConversionValueError` for more information.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import PartialConversionValueError, PartialConversionElement, GrowthGroup
            sage: raise PartialConversionValueError(
            ....:     PartialConversionElement(GrowthGroup('QQ^n'), -2), 'wrong value')
            Traceback (most recent call last):
            ...
            PartialConversionValueError: wrong value
        """
        super().__init__(*args, **kwds)
        self.element = element


class PartialConversionElement(SageObject):
    r"""
    A not converted element of a growth group.

    INPUT:

    - ``growth_group`` -- a group group

    - ``raw_element`` -- an object

    A :class:`PartialConversionElement` is an element ``growth_group(raw_element)``
    which usually appears in conjunction with :exc:`PartialConversionValueError`.
    In this case, it was to possible to create that element, although
    the conversion went partially well in the sense that a ``raw_element``
    (e.g. an exponent for :class:`MonomialGrowthElement` or a base for
    :class:`ExponentialGrowthElement`) could be extracted.

    Its main purpose is to carry data used during the creation of
    elements of
    :mod:`cartesian products of growth groups <sage.rings.asymptotic.growth_group_cartesian>`.
    """
    def __init__(self, growth_group, raw_element):
        r"""
        See :class:`PartialConversionElement` for more information.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import PartialConversionElement, GrowthGroup
            sage: PartialConversionElement(GrowthGroup('(QQ_+)^n'), -2)
            element with parameter -2 (Integer Ring) in Growth Group QQ^n
        """
        self.growth_group = growth_group
        self.raw_element = raw_element

    def _repr_(self):
        r"""
        Return a representation string of this partial conversion element.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import PartialConversionElement, GrowthGroup
            sage: PartialConversionElement(GrowthGroup('(QQ_+)^n'), -42)  # indirect doctest
            element with parameter -42 (Integer Ring) in Growth Group QQ^n
        """
        from sage.structure.element import parent
        return 'element with parameter {} ({}) in {}'.format(self.raw_element,
                                                             parent(self.raw_element),
                                                             self.growth_group)

    def split(self):
        r"""
        Split the contained ``raw_element`` according to the growth group's
        :meth:`GrowthGroup._split_raw_element_`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup, PartialConversionValueError
            sage: E = ExponentialGrowthGroup(ZZ, 'x')
            sage: try:
            ....:     E((-2)^x)
            ....: except PartialConversionValueError as e:
            ....:     e.element.split()
            (2^x, element with parameter -1 (<class 'int'>) in Growth Group ZZ^x)

        TESTS::

            sage: try:
            ....:     E((2/3)^x)
            ....: except PartialConversionValueError as e:
            ....:     e.element.split()
            Traceback (most recent call last):
            ...
            ValueError: cannot split element with parameter 2/3 (Symbolic Ring) in Growth Group ZZ^x
            > *previous* PartialConversionValueError: 2/3 (Rational Field) is not in Integer Ring
            >> *previous* TypeError: no conversion of this rational to integer
        """
        raw_here, raw_other = self.growth_group._split_raw_element_(self.raw_element)
        try:
            here = self.growth_group.element_class(self.growth_group, raw_here)
        except PartialConversionValueError as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ValueError('cannot split {}'.format(self)), e)

        other = PartialConversionElement(self.growth_group, raw_other)
        return here, other

    def is_compatible(self, other):
        r"""
        Wrapper to :meth:`GenericGrowthGroup.is_compatible`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup, ExponentialNonGrowthGroup, PartialConversionElement
            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: Q = ExponentialGrowthGroup(QQ, 'n')
            sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'n')
            sage: PartialConversionElement(Q, -42/5).is_compatible(UU)
            True
        """
        return self.growth_group.is_compatible(other)


# The following function is used in the classes GenericGrowthElement and
# GenericProduct.Element as a method.
def _is_lt_one_(self):
    r"""
    Return whether this element is less than `1`.

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('x^ZZ'); x = G(x)
        sage: (x^42).is_lt_one()  # indirect doctest
        False
        sage: (x^(-42)).is_lt_one()  # indirect doctest
        True
    """
    one = self.parent().one()
    return self <= one and self != one


# The following function is used in the classes GenericGrowthElement and
# GenericProduct.Element as a method.
def _log_(self, base=None):
    r"""
    Return the logarithm of this element.

    INPUT:

    - ``base`` -- the base of the logarithm. If ``None``
      (default value) is used, the natural logarithm is taken.

    OUTPUT: a growth element

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('x^ZZ * log(x)^ZZ')
        sage: x, = G.gens_monomial()
        sage: log(x)  # indirect doctest
        log(x)
        sage: log(x^5)  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: When calculating log(x^5) a factor 5 != 1 appeared,
        which is not contained in Growth Group x^ZZ * log(x)^ZZ.

    ::

        sage: G = GrowthGroup('(QQ_+)^x * x^ZZ')
        sage: x, = G.gens_monomial()
        sage: el = x.rpow(2); el
        2^x
        sage: log(el)  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: When calculating log(2^x) a factor log(2) != 1
        appeared, which is not contained in Growth Group QQ^x * x^ZZ.
        sage: log(el, base=2)  # indirect doctest
        x

    ::

        sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        sage: x = GenericGrowthGroup(ZZ).an_element()
        sage: log(x)  # indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: Cannot determine logarithmized factorization of
        GenericGrowthElement(1) in abstract base class.

    ::

        sage: x = GrowthGroup('x^ZZ').an_element()
        sage: log(x)  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot build log(x) since log(x) is not in
        Growth Group x^ZZ.

    TESTS::

        sage: G = GrowthGroup("(e^x)^QQ * x^ZZ")
        sage: x, = G.gens_monomial()
        sage: log(exp(x))  # indirect doctest
        x
        sage: G.one().log()  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: log(1) is zero, which is not contained in
        Growth Group (e^x)^QQ * x^ZZ.

    ::

        sage: G = GrowthGroup("(e^x)^ZZ * x^ZZ")
        sage: x, = G.gens_monomial()
        sage: log(exp(x))  # indirect doctest
        x
        sage: G.one().log()  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: log(1) is zero, which is not contained in
        Growth Group (e^x)^ZZ * x^ZZ.

    ::

        sage: G = GrowthGroup('(QQ_+)^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ')
        sage: x, y = G.gens_monomial()
        sage: (x * y).log()  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: Calculating log(x*y) results in a sum,
        which is not contained in
        Growth Group QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ.
    """
    from .misc import log_string

    log_factor = self.log_factor(base=base)
    if not log_factor:
        raise ArithmeticError('%s is zero, '
                              'which is not contained in %s.' %
                              (log_string(self, base), self.parent()))

    if len(log_factor) != 1:
        raise ArithmeticError('Calculating %s results in a sum, '
                              'which is not contained in %s.' %
                              (log_string(self, base), self.parent()))
    g, c = log_factor[0]
    if c != 1:
        raise ArithmeticError('When calculating %s a factor %s != 1 '
                              'appeared, which is not contained in %s.' %
                              (log_string(self, base), c, self.parent()))
    return g


# The following function is used in the classes GenericGrowthElement and
# GenericProduct.Element as a method.
def _log_factor_(self, base=None, locals=None):
    r"""
    Return the logarithm of the factorization of this
    element.

    INPUT:

    - ``base`` -- the base of the logarithm. If ``None``
      (default value) is used, the natural logarithm is taken.

    - ``locals`` -- dictionary which may contain the following keys and values:

      - ``'log'`` -- value: a function. If not used, then the usual
        :class:`log <sage.functions.log.Function_log>` is taken.

    OUTPUT:

    A tuple of pairs, where the first entry is a growth
    element and the second a multiplicative coefficient.

    ALGORITHM:

        This function factors the given element and calculates
        the logarithm of each of these factors.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('(QQ_+)^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ')
        sage: x, y = G.gens_monomial()
        sage: (x * y).log_factor()  # indirect doctest
        ((log(x), 1), (log(y), 1))
        sage: (x^123).log_factor()  # indirect doctest
        ((log(x), 123),)
        sage: (G('2^x') * x^2).log_factor(base=2)  # indirect doctest
        ((x, 1), (log(x), 2/log(2)))

    ::

        sage: G(1).log_factor()
        ()

    ::

        sage: log(x).log_factor()  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot build log(log(x)) since log(log(x)) is
        not in Growth Group QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ.

    .. SEEALSO::

        :meth:`~sage.rings.asymptotic.growth_group.GenericGrowthElement.factors`,
        :meth:`~sage.rings.asymptotic.growth_group.GenericGrowthElement.log`.

    TESTS::

        sage: G = GrowthGroup("(e^x)^ZZ * x^ZZ * log(x)^ZZ")
        sage: x, = G.gens_monomial()
        sage: (exp(x) * x).log_factor()  # indirect doctest
        ((x, 1), (log(x), 1))
    """
    log_factor = self._log_factor_(base=base, locals=locals)

    for g, c in log_factor:
        if hasattr(g, 'parent') and \
           isinstance(g.parent(), GenericGrowthGroup):
            continue
        from .misc import log_string
        raise ArithmeticError('Cannot build %s since %s '
                              'is not in %s.' % (log_string(self, base),
                                                 g, self.parent()))

    return log_factor


# The following function is used in the classes GenericGrowthElement and
# GenericProduct.Element as a method.
def _rpow_(self, base):
    r"""
    Calculate the power of ``base`` to this element.

    INPUT:

    - ``base`` -- an element

    OUTPUT: a growth element

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('(QQ_+)^x * x^ZZ')
        sage: x = G('x')
        sage: x.rpow(2)  # indirect doctest
        2^x
        sage: x.rpow(1/2)  # indirect doctest
        (1/2)^x

    ::

        sage: x.rpow(0)  # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: 0 is not an allowed base for calculating the power to x.
        sage: (x^2).rpow(2)  # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot construct 2^(x^2) in Growth Group QQ^x * x^ZZ
        > *previous* TypeError: unsupported operand parent(s) for *:
        'Growth Group QQ^x * x^ZZ' and 'Growth Group ZZ^(x^2)'

    ::

        sage: G = GrowthGroup('QQ^(x*log(x)) * x^ZZ * log(x)^ZZ')
        sage: x = G('x')
        sage: (x * log(x)).rpow(2)  # indirect doctest
        2^(x*log(x))

    ::

        sage: n = GrowthGroup('(QQ_+)^n * n^QQ')('n')
        sage: n.rpow(2)
        2^n
        sage: _.parent()
        Growth Group QQ^n * n^QQ

    ::

        sage: n = GrowthGroup('QQ^n * n^QQ')('n')
        sage: n.rpow(-2)
        2^n*(-1)^n

    TESTS::

        sage: SCR = SR.subring(no_variables=True)
        sage: G = GrowthGroup('QQ^x * x^ZZ')
        sage: x = G('x')
        sage: x.rpow(SCR(5))
        5^x
        sage: _.parent()
        Growth Group (Symbolic Constants Subring)^x * x^ZZ * Signs^x

    ::

        sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
        sage: UU = RootsOfUnityGroup()
        sage: asymptotic_expansions.SingularityAnalysis(
        ....:     'n', UU(-1), alpha=2, beta=1, precision=5,
        ....:     normalized=False)
        n*log(n)*(-1)^n + (euler_gamma - 1)*n*(-1)^n + log(n)*(-1)^n
        + (euler_gamma + 1/2)*(-1)^n + O(n^(-1)*(-1)^n)
        sage: _.parent()
        Asymptotic Ring <n^ZZ * log(n)^ZZ * UU^n> over Symbolic Constants Subring
    """
    if base == 0:
        raise ValueError('%s is not an allowed base for calculating the '
                         'power to %s.' % (base, self))

    var = str(self)

    try:
        element = self._rpow_element_(base)
    except ValueError:
        if base == 'e':
            from sage.rings.integer_ring import ZZ
            from .misc import repr_op
            MM = MonomialGrowthGroup(ZZ, repr_op('e', '^', var),
                                     ignore_variables=('e',))
            element = MM(raw_element=ZZ(1))
        else:
            EEUU = ExponentialGrowthGroup.factory(base.parent(), var)
            try:
                factors = EEUU.cartesian_factors()
            except AttributeError:
                factors = (EEUU,)
            if len(factors) == 1:
                EE, = factors
                element = EE(raw_element=base)
            else:
                EE, UU = factors
                try:
                    element = EE(raw_element=base)
                except PartialConversionValueError as e:
                    element = EEUU._convert_factors_([e.element])

    try:
        return self.parent().one() * element
    except (TypeError, ValueError) as e:
        from .misc import combine_exceptions, repr_op
        raise combine_exceptions(
            ArithmeticError('Cannot construct %s in %s' %
                            (repr_op(base, '^', var), self.parent())), e)


class GenericGrowthElement(MultiplicativeGroupElement):
    r"""
    A basic implementation of a generic growth element.

    Growth elements form a group by multiplication, and (some of) the
    elements can be compared to each other, i.e., all elements form a
    poset.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`

    - ``raw_element`` -- an element from the base of the parent

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import (GenericGrowthGroup,
        ....:                                                 GenericGrowthElement)
        sage: G = GenericGrowthGroup(ZZ)
        sage: g = GenericGrowthElement(G, 42); g
        GenericGrowthElement(42)
        sage: g.parent()
        Growth Group Generic(ZZ)
        sage: G(raw_element=42) == g
        True
    """

    def __init__(self, parent, raw_element):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)
            GenericGrowthElement(42)

        TESTS::

            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthElement
            sage: GenericGrowthElement(None, 0)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        super().__init__(parent=parent)

        try:
            self._raw_element_ = parent.base()(raw_element)
        except (TypeError, ValueError) as e:
            from .misc import combine_exceptions
            from sage.structure.element import parent as parent_function
            raise combine_exceptions(
                PartialConversionValueError(
                    PartialConversionElement(parent, raw_element),
                    '{} ({}) is not in {}'.format(raw_element,
                                                  parent_function(raw_element),
                                                  parent.base())),
                e)

        self._check_()

    def _check_(self):
        r"""
        Perform an additional check at the end of :meth:`__init__`.

        No check is performed for this class.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)  # indirect doctest
            GenericGrowthElement(42)
        """
        pass

    def _repr_(self):
        r"""
        A representation string for this generic element.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)  # indirect doctest
            GenericGrowthElement(42)
            sage: H = GenericGrowthGroup(ZZ, 'h')
            sage: H(raw_element=42)  # indirect doctest
            GenericGrowthElement(42, h)
        """
        vars = ', '.join(self.parent()._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'GenericGrowthElement(%s%s)' % (self._raw_element_, vars)

    def __hash__(self):
        r"""
        Return the hash of this element.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: hash(G(raw_element=42))  # random
            5656565656565656
        """
        return hash((self.parent(), self._raw_element_))

    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic elements.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`

        OUTPUT:

        A :class:`GenericGrowthElement` representing the product with
        ``other``.

        .. NOTE::

            Inherited classes must override this.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: g = G.an_element()
            sage: g*g
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')

    def __invert__(self):
        r"""
        Return the inverse of this growth element.

        OUTPUT: an instance of :class:`GenericGrowthElement`

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: ~G.an_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of GenericGrowthElement(1) not implemented
            (in this abstract method).
            sage: G.an_element()^7
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
            sage: P = GrowthGroup('x^ZZ')
            sage: ~P.an_element()
            x^(-1)
        """
        raise NotImplementedError('Inversion of %s not implemented '
                                  '(in this abstract method).' % (self,))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this :class:`GenericGrowthElement` is equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`

        OUTPUT: boolean

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G.an_element() == G.an_element()
            True
            sage: G(raw_element=42) == G(raw_element=7)
            False

        ::

            sage: G_ZZ = GenericGrowthGroup(ZZ)
            sage: G_QQ = GenericGrowthGroup(QQ)
            sage: G_ZZ(raw_element=1) == G_QQ(raw_element=1)
            True

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P_ZZ = GrowthGroup('x^ZZ')
            sage: P_QQ = GrowthGroup('x^QQ')
            sage: P_ZZ.gen() == P_QQ.gen()
            True
            sage: ~P_ZZ.gen() == P_ZZ.gen()
            False
            sage: ~P_ZZ(1) == P_ZZ(1)
            True

        TESTS::

            sage: P = GrowthGroup('x^ZZ')
            sage: e1 = P(raw_element=1)
            sage: e1 == P.gen()
            True
            sage: e2 = e1^4
            sage: e2 == e1^2*e1*e1
            True
            sage: e2 == e1
            False

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: G.one() != G(1)
            False
            sage: G.one() != G.one()
            False
            sage: G(1) != G(1)
            False
        """
        return self._raw_element_ == other._raw_element_

    def _lt_(self, other):
        r"""
        Return whether this :class:`GenericGrowthElement` is at most (less
        than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`

        OUTPUT: boolean

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P_ZZ = GrowthGroup('x^ZZ')
            sage: P_QQ = GrowthGroup('x^QQ')
            sage: P_ZZ.gen() <= P_QQ.gen()^2
            True
            sage: ~P_ZZ.gen() <= P_ZZ.gen()
            True

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: e1 = G(raw_element=1); e2 = G(raw_element=2)
            sage: e1 <= e2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')

    log = _log_
    log_factor = _log_factor_

    def _log_factor_(self, base=None, locals=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(QQ)
            sage: G.an_element().log_factor()  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot determine logarithmized factorization of
            GenericGrowthElement(1/2) in abstract base class.
        """
        raise NotImplementedError('Cannot determine logarithmized factorization '
                                  'of %s in abstract base class.' % (self,))

    rpow = _rpow_

    def _rpow_element_(self, base):
        r"""
        Return an element which is the power of ``base`` to this
        element.

        INPUT:

        - ``base`` -- an element

        OUTPUT: nothing since a :exc:`ValueError` is raised in this generic method

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('(QQ_+)^x')
            sage: x = G(raw_element=3)
            sage: x._rpow_element_(2) is None
            Traceback (most recent call last):
            ...
            ValueError: Cannot compute 2 to the generic element 3^x.
        """
        raise ValueError('Cannot compute %s to the generic element %s.' %
                         (base, self))

    def factors(self):
        r"""
        Return the atomic factors of this growth element. An atomic factor
        cannot be split further.

        OUTPUT: a tuple of growth elements

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: G.an_element().factors()
            (x,)
        """
        return (self,)

    is_lt_one = _is_lt_one_

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this generic growth element.

        INPUT:

        - ``rules`` -- dictionary

        OUTPUT: an object

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)._substitute_({})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute in GenericGrowthElement(42) in
            Growth Group Generic(ZZ).
            > *previous* TypeError: Cannot substitute in the abstract base class
            Growth Group Generic(ZZ).
        """
        from .misc import substitute_raise_exception
        substitute_raise_exception(self, TypeError(
            'Cannot substitute in the abstract '
            'base class %s.' % (self.parent(),)))

    def variable_names(self):
        r"""
        Return the names of the variables of this growth element.

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('m^QQ')
            sage: G('m^2').variable_names()
            ('m',)
            sage: G('m^0').variable_names()
            ()

        ::

            sage: G = GrowthGroup('QQ^m')
            sage: G('2^m').variable_names()
            ('m',)
            sage: G('1^m').variable_names()
            ()

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(QQ)
            sage: G(raw_element=2).variable_names()
            Traceback (most recent call last):
            ...
            AttributeError: 'GenericGrowthGroup_with_category.element_class' object
            has no attribute 'is_one'...
        """
        if self.is_one():
            return tuple()
        else:
            return self.parent().variable_names()

    def _singularity_analysis_(self, var, zeta, precision):
        r"""
        Perform singularity analysis on this growth element.

        INPUT:

        - ``var`` -- string denoting the variable

        - ``zeta`` -- a number

        - ``precision`` -- integer

        OUTPUT:

        An asymptotic expansion for `[z^n] f` where `n` is ``var``
        and `f` has this growth element as a singular expansion
        in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this element
        is a growth element in `T`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=2)._singularity_analysis_('n', 2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedError: singularity analysis of GenericGrowthElement(2)
            not implemented
        """
        raise NotImplementedError('singularity analysis of {} '
                                  'not implemented '.format(self))

    def _find_minimum_(self, valid_from):
        r"""
        Find the minimum of this growth element over the range implied by ``valid_from``.

        INPUT:

        - ``valid_from`` -- dictionary describing the range of the minimization:
          the keys are names of variables and the range is the intersection over
          the ranges where the absolute value of the variable designated by the
          key is at least the corresponding value

        OUTPUT: a number

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup

            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)._find_minimum_(valid_from={'m': 10})
            Traceback (most recent call last):
            ...
            NotImplementedError: find minimum for GenericGrowthElement(42) not implemented
        """
        raise NotImplementedError(f'find minimum for {self} not implemented')


class GenericGrowthGroup(UniqueRepresentation, Parent, WithLocals):
    r"""
    A basic implementation for growth groups.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``)

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    - ``ignore_variables`` -- (default: ``None``) a tuple (or other
      iterable) of strings. The specified names are not considered as
      variables.

    .. NOTE::

        This class should be derived for concrete implementations.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        sage: G = GenericGrowthGroup(ZZ); G
        Growth Group Generic(ZZ)

    .. SEEALSO::

        :class:`MonomialGrowthGroup`,
        :class:`ExponentialGrowthGroup`
    """
    # TODO: implement some sort of 'assume', where basic assumptions
    # for the variables can be stored. --> within the Cartesian product

    # enable the category framework for elements
    Element = GenericGrowthElement

    # set everything up to determine category
    from sage.categories.sets_cat import Sets
    from sage.categories.posets import Posets
    from sage.categories.magmas import Magmas
    from sage.categories.additive_magmas import AdditiveMagmas

    _determine_category_subcategory_mapping_ = [
        (Sets(), Sets(), True),
        (Posets(), Posets(), False)]

    _determine_category_axiom_mapping_ = []

    @staticmethod
    def __classcall__(cls, base, var=None, category=None, ignore_variables=None):
        r"""
        Normalize the input in order to ensure a unique
        representation.

        For more information see :class:`GenericGrowthGroup`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: P1 = MonomialGrowthGroup(ZZ, 'x')
            sage: P2 = MonomialGrowthGroup(ZZ, ZZ['x'].gen())
            sage: P3 = MonomialGrowthGroup(ZZ, SR.var('x'))
            sage: P1 is P2 and P2 is P3
            True
            sage: P5 = MonomialGrowthGroup(ZZ, 'x ')
            sage: P1 is P5
            True

        ::

            sage: L1 = MonomialGrowthGroup(QQ, log(x))
            sage: L2 = MonomialGrowthGroup(QQ, 'log(x)')
            sage: L1 is L2
            True

        Test determining of the category (:class:`GenericGrowthGroup`)::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(ZZ, 'x').category()  # indirect doctest
            Category of posets
            sage: GenericGrowthGroup(ZZ, 'x', category=Groups()).category()  # indirect doctest
            Category of groups

        Test determining of the category (:class:`MonomialGrowthGroup`)::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: MonomialGrowthGroup(ZZ, 'x').category()  # indirect doctest
            Join of Category of commutative groups and Category of posets
            sage: MonomialGrowthGroup(ZZ, 'x', category=Monoids()).category()  # indirect doctest
            Category of monoids
            sage: W = Words([0, 1])
            sage: W.category()
            Category of sets
            sage: MonomialGrowthGroup(W, 'x').category()  # indirect doctest
            Category of sets

        Test determining of the category (:class:`ExponentialGrowthGroup`)::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(ZZ, 'x').category()  # indirect doctest
            Join of Category of commutative monoids and Category of posets
            sage: ExponentialGrowthGroup(QQ, 'x').category()  # indirect doctest
            Join of Category of commutative groups and Category of posets
            sage: ExponentialGrowthGroup(ZZ, 'x', category=Groups()).category()  # indirect doctest
            Category of groups
            sage: ExponentialGrowthGroup(QQ, 'x', category=Monoids()).category()  # indirect doctest
            Category of monoids

        ::

            sage: MonomialGrowthGroup(AsymptoticRing('z^ZZ', QQ), 'x')
            Traceback (most recent call last):
            ...
            TypeError: Asymptotic Ring <z^ZZ> over Rational Field is not a valid base.
        """
        from .asymptotic_ring import AsymptoticRing
        if not isinstance(base, Parent) or \
           isinstance(base, AsymptoticRing):
            raise TypeError('%s is not a valid base.' % (base,))

        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var, ignore=ignore_variables)

        if category is None:
            from .misc import transform_category
            category = transform_category(
                base.category(),
                cls._determine_category_subcategory_mapping_,
                cls._determine_category_axiom_mapping_,
                initial_category=cls._initial_category_(base))

        return super().__classcall__(cls, base, var, category)

    @staticmethod
    def _initial_category_(base):
        r"""
        Return a category with which creating the actual category
        of this growth group starts.

        INPUT:

        - ``base`` -- a SageMath parent

        OUTPUT: a category or ``None``

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup._initial_category_(ZZ)
            Category of posets
            sage: GenericGrowthGroup._initial_category_(QQ)
            Category of posets
            sage: GenericGrowthGroup._initial_category_(SR) is None
            True
        """
        from sage.categories.posets import Posets
        # The following block can be removed once #19269 is fixed.
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        if base is ZZ or base is QQ or \
                isinstance(base, PolynomialRing_general) and \
                (base.base_ring() is ZZ or base.base_ring() is QQ):
            return Posets()
        else:
            return None

    def __init__(self, base, var, category):
        r"""
        See :class:`GenericGrowthGroup` for more information.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(ZZ).category()
            Category of posets

        ::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: MonomialGrowthGroup(ZZ, 'x')
            Growth Group x^ZZ
            sage: MonomialGrowthGroup(QQ, SR.var('n'))
            Growth Group n^QQ
            sage: MonomialGrowthGroup(ZZ, ZZ['y'].gen())
            Growth Group y^ZZ
            sage: MonomialGrowthGroup(QQ, 'log(x)')
            Growth Group log(x)^QQ

        ::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(QQ, 'x')
            Growth Group QQ^x
            sage: assume(SR.an_element() > 0)
            sage: ExponentialGrowthGroup(SR, ZZ['y'].gen())
            Growth Group SR^y
            sage: forget()

        TESTS::

            sage: G = GenericGrowthGroup(ZZ)
            sage: G.is_parent_of(G(raw_element=42))
            True
            sage: G2 = GenericGrowthGroup(ZZ, category=FiniteGroups() & Posets())
            sage: G2.category()
            Join of Category of finite groups and Category of finite posets

        ::

            sage: G = GenericGrowthGroup('42')
            Traceback (most recent call last):
            ...
            TypeError: 42 is not a valid base.

        ::

            sage: MonomialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.
            sage: MonomialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.

        ::

            sage: ExponentialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.
            sage: ExponentialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.
        """
        self._var_ = var
        super().__init__(category=category, base=base)

    def _repr_short_(self):
        r"""
        A short representation string of this abstract growth group.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(QQ)._repr_short_()
            'Generic(QQ)'
            sage: GenericGrowthGroup(QQ)
            Growth Group Generic(QQ)
            sage: GenericGrowthGroup(QQ, ('a', 'b'))
            Growth Group Generic(QQ, a, b)
        """
        from .misc import parent_to_repr_short
        vars = ', '.join(self._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'Generic(%s%s)' % (parent_to_repr_short(self.base()), vars)

    def _repr_(self, condense=False):
        r"""
        A representations string of this growth group.

        INPUT:

        - ``condense`` -- boolean (default: ``False``); if set, then a shorter
          output is returned, e.g. the prefix-string ``Growth Group`` is not
          shown in this case

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ')  # indirect doctest
            Growth Group x^ZZ
            sage: GrowthGroup('log(x)^QQ')  # indirect doctest
            Growth Group log(x)^QQ

        TESTS::

            sage: GrowthGroup('log(x)^QQ')._repr_(condense=True)
            'log(x)^QQ'
        """
        pre = 'Growth Group ' if not condense else ''
        return '%s%s' % (pre, self._repr_short_())

    def __hash__(self):
        r"""
        Return the hash of this group.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import (GenericGrowthGroup,
            ....:                                                 GrowthGroup)
            sage: hash(GenericGrowthGroup(ZZ))  # random
            4242424242424242

        ::

            sage: P = GrowthGroup('x^ZZ')
            sage: hash(P)  # random
            -1234567890123456789

        ::

            sage: P = GrowthGroup('QQ^x')
            sage: hash(P)  # random
            -1234567890123456789
        """
        return hash((self.__class__, self.base(), self._var_))

    def _an_element_(self):
        r"""
        Return an element of this growth group.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import (GenericGrowthGroup,
            ....:                                                 GrowthGroup)
            sage: GenericGrowthGroup(ZZ).an_element()  # indirect doctest
            GenericGrowthElement(1)
            sage: GrowthGroup('z^ZZ').an_element()  # indirect doctest
            z
            sage: GrowthGroup('log(z)^QQ').an_element()  # indirect doctest
            log(z)^(1/2)
            sage: GrowthGroup('(QQ_+)^(x*log(x))').an_element()  # indirect doctest
            (1/2)^(x*log(x))
            sage: GrowthGroup('QQ^(x*log(x))').an_element()  # indirect doctest
            (1/2)^(x*log(x))*(-1)^(x*log(x))
        """
        return self.element_class(self, self.base().an_element())

    def some_elements(self):
        r"""
        Return some elements of this growth group.

        See :class:`TestSuite` for a typical use case.

        OUTPUT: an iterator

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: tuple(GrowthGroup('z^ZZ').some_elements())
            (1, z, z^(-1), z^2, z^(-2), z^3, z^(-3),
             z^4, z^(-4), z^5, z^(-5), ...)
            sage: tuple(GrowthGroup('z^QQ').some_elements())
            (z^(1/2), z^(-1/2), z^2, z^(-2),
             1, z, z^(-1), z^42,
             z^(2/3), z^(-2/3), z^(3/2), z^(-3/2),
             z^(4/5), z^(-4/5), z^(5/4), z^(-5/4), ...)
        """
        return iter(self.element_class(self, e)
                    for e in self.base().some_elements())

    def _create_element_in_extension_(self, raw_element):
        r"""
        Create an element in an extension of this growth group which
        is chosen according to the input ``raw_element``.

        INPUT:

        - ``raw_element`` -- the element data

        OUTPUT: an element

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ')
            sage: G._create_element_in_extension_(3).parent()
            Growth Group z^ZZ
            sage: G._create_element_in_extension_(1/2).parent()
            Growth Group z^QQ
        """
        if raw_element.parent() is self.base():
            parent = self
        else:
            parent = self._underlying_class()(raw_element.parent(), self._var_,
                                              category=self.category())
        return parent(raw_element=raw_element)

    def le(self, left, right):
        r"""
        Return whether the growth of ``left`` is at most (less than or
        equal to) the growth of ``right``.

        INPUT:

        - ``left`` -- an element

        - ``right`` -- an element

        OUTPUT: boolean

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: x = G.gen()
            sage: G.le(x, x^2)
            True
            sage: G.le(x^2, x)
            False
            sage: G.le(x^0, 1)
            True
        """
        return self(left) <= self(right)

    def _element_constructor_(self, data, raw_element=None):
        r"""
        Convert a given object to this growth group.

        INPUT:

        - ``data`` -- an object representing the element to be initialized

        - ``raw_element`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no conversion
          is performed).

        OUTPUT: an element of this growth group

        .. NOTE::

            Either ``data`` or ``raw_element`` has to be given. If
            ``raw_element`` is specified, then no positional argument
            may be passed.

            This method calls :meth:`_convert_`, which does the actual
            conversion from ``data``.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G_ZZ = GenericGrowthGroup(ZZ)
            sage: z = G_ZZ(raw_element=42); z  # indirect doctest
            GenericGrowthElement(42)
            sage: z is G_ZZ(z)  # indirect doctest
            True

        ::

            sage: G_QQ = GenericGrowthGroup(QQ)
            sage: q = G_QQ(raw_element=42)  # indirect doctest
            sage: q is z
            False
            sage: G_ZZ(q)  # indirect doctest
            GenericGrowthElement(42)
            sage: G_QQ(z)  # indirect doctest
            GenericGrowthElement(42)
            sage: q is G_ZZ(q)  # indirect doctest
            False

        ::

            sage: G_ZZ()
            Traceback (most recent call last):
            ...
            ValueError: No input specified. Cannot continue.
            sage: G_ZZ('blub')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: blub is not in Growth Group Generic(ZZ).
            sage: G_ZZ('x', raw_element=42)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: input is ambiguous: x as well as raw_element=42 are specified

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: x = GrowthGroup('x^ZZ')(raw_element=1)  # indirect doctest
            sage: G_y = GrowthGroup('y^ZZ')
            sage: G_y(x)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: x is not in Growth Group y^ZZ.

        ::

            sage: GrowthGroup('QQ^x')(GrowthGroup('ZZ^x')('2^x'))
            2^x

        ::

            sage: from sage.rings.asymptotic.growth_group import PartialConversionValueError
            sage: G = GrowthGroup('(QQ_+)^n')
            sage: n = SR.var('n')
            sage: try:
            ....:     G((-1/42)^n)
            ....: except PartialConversionValueError as e:
            ....:     a = e.element
            sage: a
            element with parameter -1/42 (Rational Field) in Growth Group QQ^n
            sage: G(a)
            Traceback (most recent call last):
            ...
            PartialConversionValueError: no conversion of
            element with parameter -1/42 (Rational Field)
            in Growth Group QQ^n:
            this was already unsuccessful earlier
            sage: UU = GrowthGroup('UU^n')
            sage: b, c = a.split()
            sage: G(b), UU(c)
            ((1/42)^n, (-1)^n)
        """
        from .misc import combine_exceptions

        if raw_element is None:
            if isinstance(data, int) and data == 0:
                raise ValueError('No input specified. Cannot continue.')

            elif isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                if self._var_ != data.parent()._var_:
                    raise ValueError('%s is not in %s.' % (data, self))
                raw_element = data._raw_element_

            elif isinstance(data, self.Element):
                if self._var_ == data.parent()._var_:
                    try:
                        raw_element = self.base()(data._raw_element_)
                    except (TypeError, ValueError) as e:
                        raise combine_exceptions(
                            ValueError('%s is not in %s.' % (data, self)), e)

            elif isinstance(data, GenericGrowthElement):
                if data.is_one():
                    return self.one()

            elif isinstance(data, PartialConversionElement):
                if data.growth_group is self:
                    raise PartialConversionValueError(
                        data,
                        'no conversion of {}: this was already unsuccessful '
                        'earlier'.format(data))
                if not data.is_compatible(self):
                    raise TypeError(
                        'cannot (partially) convert {} because its '
                        'growth group {} is not compatible to this '
                        'growth group {}'.format(data.raw_element, data.growth_group, self))
                raw_element = data.raw_element

            else:
                raw_element = self._convert_(data)

            if raw_element is None:
                raise ValueError('%s is not in %s.' % (data, self))
        elif not isinstance(data, int) or data != 0:
            raise ValueError('input is ambiguous: '
                             '%s as well as raw_element=%s '
                             'are specified' % (data, raw_element))

        return self.element_class(self, raw_element)

    def _convert_(self, data):
        r"""
        Convert ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            This method always returns ``None`` in this abstract base
            class, and should be overridden in inherited class.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G._convert_('icecream') is None
            True
        """
        pass

    def is_compatible(self, other):
        r"""
        Return whether this growth group is compatible with ``other`` meaning
        that both are of the same type and have the same variables, but
        maybe a different base.

        INPUT:

        - ``other`` -- a growth group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup, ExponentialNonGrowthGroup
            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: EQ = ExponentialGrowthGroup(QQ, 'n')
            sage: EZ = ExponentialGrowthGroup(ZZ, 'n')
            sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'n')
            sage: for a in (EQ, EZ, UU):
            ....:     for b in (EQ, EZ, UU):
            ....:         print('{} is {}compatible with {}'.format(
            ....:             a, '' if a.is_compatible(b) else 'not ', b))
            Growth Group QQ^n is compatible with Growth Group QQ^n
            Growth Group QQ^n is compatible with Growth Group ZZ^n
            Growth Group QQ^n is compatible with Growth Group UU^n
            Growth Group ZZ^n is compatible with Growth Group QQ^n
            Growth Group ZZ^n is compatible with Growth Group ZZ^n
            Growth Group ZZ^n is compatible with Growth Group UU^n
            Growth Group UU^n is not compatible with Growth Group QQ^n
            Growth Group UU^n is not compatible with Growth Group ZZ^n
            Growth Group UU^n is compatible with Growth Group UU^n
        """
        return isinstance(other, self._underlying_class()) and self._var_ == other._var_

    @staticmethod
    def _split_raw_element_(raw_element):
        r"""
        Split ``raw_element`` in a part convertible to this growth group
        and a part which needs to be converted by some compatible growth group.

        INPUT:

        - ``raw_element`` -- an object

        OUTPUT: a pair of objects

        .. NOTE::

            This method is called by
            :meth:`~sage.rings.asymptotic.growth_group.PartialConversionElement.split`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G._split_raw_element_(0)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented in concrete realizations
        """
        raise NotImplementedError('only implemented in concrete realizations')

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent

        OUTPUT: boolean

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: G_QQ = GrowthGroup('x^QQ')
            sage: G_ZZ.has_coerce_map_from(G_QQ)  # indirect doctest
            False
            sage: G_QQ.has_coerce_map_from(G_ZZ)  # indirect doctest
            True

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P_x_ZZ = GrowthGroup('x^ZZ')
            sage: P_x_QQ = GrowthGroup('x^QQ')
            sage: P_x_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            True
            sage: P_y_ZZ = GrowthGroup('y^ZZ')
            sage: P_y_ZZ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            False
            sage: P_x_ZZ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False
            sage: P_y_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P_x_ZZ = GrowthGroup('ZZ^x')
            sage: P_x_QQ = GrowthGroup('QQ^x')
            sage: P_x_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            True
            sage: P_y_ZZ = GrowthGroup('ZZ^y')
            sage: P_y_ZZ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            False
            sage: P_x_ZZ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False
            sage: P_y_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False

        ::

            sage: GrowthGroup('x^QQ').has_coerce_map_from(GrowthGroup('QQ^x'))  # indirect doctest
            False
        """
        if isinstance(S, self._underlying_class()) and self._var_ == S._var_:
            if self.base().has_coerce_map_from(S.base()):
                return True

    def _pushout_(self, other):
        r"""
        Construct the pushout of this and the other growth group. This is called by
        :func:`sage.categories.pushout.pushout`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.categories.pushout import pushout
            sage: cm = sage.structure.element.get_coercion_model()
            sage: A = GrowthGroup('(QQ_+)^x')
            sage: B = GrowthGroup('y^ZZ')

        When using growth groups with disjoint variable lists, then a
        pushout can be constructed::

            sage: A._pushout_(B)
            Growth Group QQ^x * y^ZZ
            sage: cm.common_parent(A, B)
            Growth Group QQ^x * y^ZZ

        In general, growth groups of the same variable cannot be
        combined automatically, since there is no order relation between the two factors::

            sage: C = GrowthGroup('x^QQ')
            sage: cm.common_parent(A, C)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Growth Group QQ^x' and 'Growth Group x^QQ'

        However, combining is possible if the factors with the same variable
        overlap::

            sage: cm.common_parent(GrowthGroup('x^ZZ * log(x)^ZZ'), GrowthGroup('exp(x)^ZZ * x^ZZ'))
            Growth Group exp(x)^ZZ * x^ZZ * log(x)^ZZ
            sage: cm.common_parent(GrowthGroup('x^ZZ * log(x)^ZZ'), GrowthGroup('y^ZZ * x^ZZ'))
            Growth Group x^ZZ * log(x)^ZZ * y^ZZ

        ::

            sage: cm.common_parent(GrowthGroup('x^ZZ'), GrowthGroup('y^ZZ'))
            Growth Group x^ZZ * y^ZZ

        ::

            sage: cm.record_exceptions()
            sage: cm.common_parent(GrowthGroup('x^ZZ'), GrowthGroup('y^ZZ'))
            Growth Group x^ZZ * y^ZZ
            sage: sage.structure.element.coercion_traceback()  # not tested
        """
        if not isinstance(other, GenericGrowthGroup) and \
           not (other.construction() is not None and
                isinstance(other.construction()[0], AbstractGrowthGroupFunctor)):
            return

        if set(self.variable_names()).isdisjoint(set(other.variable_names())):
            from sage.categories.cartesian_product import cartesian_product
            return cartesian_product([self, other])

    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth
        group.

        OUTPUT: an empty tuple

        .. NOTE::

            A generator is called monomial generator if the variable
            of the underlying growth group is a valid identifier. For
            example, ``x^ZZ`` has ``x`` as a monomial generator,
            while ``log(x)^ZZ`` or ``icecream(x)^ZZ`` do not have
            monomial generators.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(ZZ).gens_monomial()
            ()

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^QQ').gens_monomial()
            (x,)
            sage: GrowthGroup('QQ^x').gens_monomial()
            ()
        """
        return tuple()

    def gens(self):
        r"""
        Return a tuple of all generators of this growth group.

        OUTPUT: a tuple whose entries are growth elements

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: P.gens()
            (x,)
            sage: GrowthGroup('log(x)^ZZ').gens()
            (log(x),)
        """
        return (self(raw_element=self.base().one()),)

    def gen(self, n=0):
        r"""
        Return the `n`-th generator (as a group) of this growth group.

        INPUT:

        - ``n`` -- (default: `0`)

        OUTPUT: a :class:`MonomialGrowthElement`

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: P.gen()
            x

        ::

            sage: P = GrowthGroup('(QQ_+)^x')
            sage: P.gen()
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.gens()[n]

    def ngens(self):
        r"""
        Return the number of generators (as a group) of this growth group.

        OUTPUT: a Python integer

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: P.ngens()
            1
            sage: GrowthGroup('log(x)^ZZ').ngens()
            1

        ::

            sage: P = GrowthGroup('(QQ_+)^x')
            sage: P.ngens()
            0
        """
        return len(self.gens())

    def variable_names(self):
        r"""
        Return the names of the variables of this growth group.

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(ZZ).variable_names()
            ()

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').variable_names()
            ('x',)
            sage: GrowthGroup('log(x)^ZZ').variable_names()
            ('x',)

        ::

            sage: GrowthGroup('(QQ_+)^x').variable_names()
            ('x',)
            sage: GrowthGroup('(QQ_+)^(x*log(x))').variable_names()
            ('x',)
        """
        return self._var_.variable_names()

    CartesianProduct = CartesianProductGrowthGroups

    def extended_by_non_growth_group(self):
        r"""
        Extend to a cartesian product of this growth group
        and a suitable non growth group.

        OUTPUT: a group group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('(QQ_+)^x').extended_by_non_growth_group()
            Growth Group QQ^x * Signs^x
            sage: GrowthGroup('(RR_+)^x').extended_by_non_growth_group()
            Growth Group RR^x * Signs^x
            sage: GrowthGroup('(RIF_+)^x').extended_by_non_growth_group()
            Growth Group RIF^x * Signs^x
            sage: GrowthGroup('(RBF_+)^x').extended_by_non_growth_group()
            Growth Group RBF^x * Signs^x
            sage: GrowthGroup('(CC_+)^x').extended_by_non_growth_group()
            Growth Group CC^x * UU_RR^x
            sage: GrowthGroup('(CIF_+)^x').extended_by_non_growth_group()
            Growth Group CIF^x * UU_RIF^x
            sage: GrowthGroup('(CBF_+)^x').extended_by_non_growth_group()
            Growth Group CBF^x * UU_RBF^x
        """
        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product((self, self.non_growth_group()))

    def non_growth_group(self):
        r"""
        Return a non-growth group compatible with this growth group.

        OUTPUT: a group group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(ZZ, 'n').non_growth_group()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented in concrete realizations
        """
        raise NotImplementedError('only implemented in concrete realizations')


class AbstractGrowthGroupFunctor(ConstructionFunctor):
    r"""
    A base class for the functors constructing growth groups.

    INPUT:

    - ``var`` -- string or list of strings (or anything else :class:`Variable`
      accepts)

    - ``domain`` -- a category

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: GrowthGroup('z^QQ').construction()[0]  # indirect doctest
        MonomialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`ExponentialGrowthGroupFunctor`,
        :class:`MonomialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.
    """

    _functor_name = 'AbstractGrowthGroup'

    rank = 13

    def __init__(self, var, domain):
        r"""
        See :class:`AbstractGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import AbstractGrowthGroupFunctor
            sage: AbstractGrowthGroupFunctor('x', Groups())
            AbstractGrowthGroup[x]
        """
        from sage.categories.monoids import Monoids
        from sage.categories.posets import Posets

        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var)
        self.var = var
        super().__init__(domain, Monoids() & Posets())

    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('(QQ_+)^t').construction()[0]  # indirect doctest
            ExponentialGrowthGroup[t]
        """
        return '%s[%s]' % (self._functor_name, self.var)

    def merge(self, other):
        r"""
        Merge this functor with ``other`` of possible.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: a functor or ``None``

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('(QQ_+)^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F.merge(F)
            ExponentialGrowthGroup[t]
            sage: F.merge(G) is None
            True
        """
        if self == other:
            return self

    def __eq__(self, other):
        r"""
        Return whether this functor is equal to ``other``.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('QQ^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F == F
            True
            sage: F == G
            False
        """
        return type(self) is type(other) and self.var == other.var

    def __ne__(self, other):
        r"""
        Return whether this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('QQ^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F != F
            False
            sage: F != G
            True
        """
        return not self == other


class DecreasingGrowthElementError(ValueError):
    r"""
    A special :python:`ValueError<library/exceptions.html#exceptions.ValueError>`
    which is raised when a growth element is less than one.

    INPUT:

    - ``element`` -- a :class:`GenericGrowthElement`

    The remaining arguments are passed on to
    :python:`ValueError<library/exceptions.html#exceptions.ValueError>`.
    """
    def __init__(self, element, *args, **kwds):
        r"""
        See :exc:`DecreasingGrowthElementError` for more information.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import DecreasingGrowthElementError, GenericGrowthElement, MonomialGrowthGroup
            sage: raise DecreasingGrowthElementError(
            ....:     GenericGrowthElement(MonomialGrowthGroup(QQ, 'x'), 1/2), 'wrong value')
            Traceback (most recent call last):
            ...
            DecreasingGrowthElementError: wrong value
        """
        super().__init__(*args, **kwds)
        self.element = element


class MonomialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of monomial growth elements.

    INPUT:

    - ``parent`` -- a :class:`MonomialGrowthGroup`

    - ``raw_element`` -- an element from the base ring of the parent

      This ``raw_element`` is the exponent of the created monomial
      growth element.

    A monomial growth element represents a term of the type
    `\operatorname{variable}^{\operatorname{exponent}}`. The multiplication
    corresponds to the addition of the exponents.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
        sage: P = MonomialGrowthGroup(ZZ, 'x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """

    @property
    def exponent(self):
        r"""
        The exponent of this growth element.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: P(x^42).exponent
            42
        """
        return self._raw_element_

    def _repr_(self, latex=False):
        r"""
        A representation string for this monomial growth element.

        INPUT:

        - ``latex`` -- boolean (default: ``False``); if set, then LaTeX-output
          is returned

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^QQ')
            sage: P(1)._repr_()
            '1'
            sage: P(x^5)  # indirect doctest
            x^5
            sage: P(x^(1/2))  # indirect doctest
            x^(1/2)

        TESTS::

            sage: P(x^-1)  # indirect doctest
            x^(-1)
            sage: P(x^-42)  # indirect doctest
            x^(-42)

        ::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: from sage.rings.asymptotic.growth_group import MonomialNonGrowthGroup
            sage: J = MonomialNonGrowthGroup(ImaginaryGroup(ZZ), 'n')
            sage: J.an_element()  # indirect doctest
            n^I
        """
        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
        else:
            f = repr

        from sage.symbolic.ring import isidentifier
        from sage.rings.integer_ring import ZZ
        from .misc import repr_op

        var = f(self.parent()._var_)
        if self.exponent.is_zero():
            return '1'
        elif self.exponent == 1:
            return var
        elif latex:
            return repr_op(var, '^', latex=True) + \
                '{' + latex_repr(self.exponent)._latex_() + '}'
        elif self.exponent in ZZ and self.exponent > 0 \
                or isidentifier(str(self.exponent)):
            return repr_op(var, '^') + str(self.exponent)
        else:
            return repr_op(var, '^') + '(' + str(self.exponent) + ')'

    def _latex_(self):
        r"""
        A LaTeX-representation string for this monomial growth element.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^QQ')
            sage: latex(P(1))  # indirect doctest
            1
            sage: latex(P(x^5))  # indirect doctest
            x^{5}
            sage: latex(P(x^(1/2)))  # indirect doctest
            x^{\frac{1}{2}}

        ::

            sage: latex(P(x^-1))  # indirect doctest
            x^{-1}
            sage: latex(P(x^-42))  # indirect doctest
            x^{-42}
        """
        return self._repr_(latex=True)

    def _mul_(self, other):
        r"""
        Multiply this monomial growth element with another.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`

        OUTPUT: the product as a :class:`MonomialGrowthElement`

        .. NOTE::

            Two monomial growth elements are multiplied by adding
            their exponents.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: a = P(x^2)
            sage: b = P(x^3)
            sage: c = a._mul_(b); c
            x^5
            sage: c == a*b
            True
            sage: a*b*a  # indirect doctest
            x^7
        """
        return self.parent()(raw_element=self.exponent + other.exponent)

    def __invert__(self):
        r"""
        Return the multiplicative inverse of this monomial growth element.

        OUTPUT: the multiplicative inverse as a :class:`MonomialGrowthElement`

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            x^(-2)
            sage: e2 == ~e1
            True
            sage: Q = GrowthGroup('x^NN'); Q
            Growth Group x^(Non negative integer semiring)
            sage: e3 = ~Q('x'); e3
            x^(-1)
            sage: e3.parent()
            Growth Group x^ZZ
        """
        return self.parent()._create_element_in_extension_(-self.exponent)

    def __pow__(self, exponent):
        r"""
        Calculate the power of this growth element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- a number. This can be anything that is a
          valid right hand side of ``*`` with elements of the
          parent's base.

        OUTPUT: the result of this exponentiation, a :class:`MonomialGrowthElement`

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: x = P.gen()
            sage: a = x^7; a
            x^7
            sage: a^(1/2)
            x^(7/2)
            sage: (a^(1/2)).parent()
            Growth Group x^QQ
            sage: a^(1/7)
            x
            sage: (a^(1/7)).parent()
            Growth Group x^QQ
            sage: P = GrowthGroup('x^QQ')
            sage: b = P.gen()^(7/2); b
            x^(7/2)
            sage: b^12
            x^42
        """
        from .misc import strip_symbolic
        return self.parent()._create_element_in_extension_(
            self.exponent * strip_symbolic(exponent))

    def _log_factor_(self, base=None, locals=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None`` (default value) is
          used, the natural logarithm is taken.

        - ``locals`` -- dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^QQ')
            sage: G('x').log_factor()  # indirect doctest
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot build log(x) since log(x) is not in
            Growth Group x^QQ.

        ::

            sage: G = GrowthGroup('exp(x)^ZZ * x^ZZ')
            sage: log(G('exp(x)'), base=2)
            Traceback (most recent call last):
            ...
            ArithmeticError: When calculating log(exp(x), base=2) a factor
            1/log(2) != 1 appeared, which is not contained in
            Growth Group exp(x)^ZZ * x^ZZ.

        ::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: L.<log2> = ZZ[]
            sage: G = MonomialGrowthGroup(L, 'x')
            sage: G(raw_element=log2)._log_factor_(base=2)
            (('log(x)', log2/log(2)),)
            sage: G(raw_element=log2)._log_factor_(base=2,
            ....:       locals={'log': lambda z: log2 if z == 2 else log(z)})
            (('log(x)', 1),)
        """
        if self.is_one():
            return tuple()
        coefficient = self.exponent

        var = str(self.parent()._var_)

        from .misc import split_str_by_op
        split = split_str_by_op(var, '^')
        if len(split) == 2:
            b, e = split
            if base is None and b == 'e' or \
               base is not None and b == str(base):
                return ((e, coefficient),)

        if var.startswith('exp('):
            assert (var[-1] == ')')
            v = var[4:-1]
        else:
            v = 'log(%s)' % (var,)

        if base is not None:
            log = self.parent().locals(locals)['log']
            coefficient = coefficient / log(base)
        return ((v, coefficient),)

    def _rpow_element_(self, base, locals=None):
        r"""
        Return an element which is the power of ``base`` to this
        element.

        INPUT:

        - ``base`` -- an element

        - ``locals`` -- dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT: a growth element

        .. NOTE::

            The parent of the result can be different from the parent
            of this element.

        A :exc:`ValueError` is raised if the calculation is not possible
        within this method. (Then the calling method should take care
        of the calculation.)

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: x = G('x')
            sage: x._rpow_element_(2)
            Traceback (most recent call last):
            ...
            ValueError: Variable x is not a log of something.

        The previous example does not work since the result would not
        live in a monomial growth group. When using
        :meth:`~GenericGrowthElement.rpow`, this
        case is handled by the calling method :meth:`_rpow_`.

        ::

            sage: G = GrowthGroup('log(x)^ZZ')
            sage: lx = G(raw_element=1); lx
            log(x)
            sage: rp = lx._rpow_element_('e'); rp
            x
            sage: rp.parent()
            Growth Group x^ZZ

        ::

            sage: G = GrowthGroup('log(x)^SR')
            sage: lx = G('log(x)')
            sage: lx._rpow_element_(2)
            x^(log(2))
        """
        var = str(self.parent()._var_)
        if not (var.startswith('log(') and self.exponent.is_one()):
            raise ValueError('Variable %s is not a log of something.' % (var,))
        new_var = var[4:-1]
        if base == 'e':
            from sage.rings.integer_ring import ZZ
            M = MonomialGrowthGroup(ZZ, new_var)
            return M(raw_element=ZZ(1))
        else:
            log = self.parent().locals(locals)['log']
            new_exponent = log(base)
            M = MonomialGrowthGroup(new_exponent.parent(), new_var)
            return M(raw_element=new_exponent)

    def _lt_(self, other):
        r"""
        Return whether this :class:`MonomialGrowthElement` is
        less than ``other``.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`

        OUTPUT: boolean

        .. NOTE::

            This function compares two instances of
            :class:`MonomialGrowthElement`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P_ZZ = GrowthGroup('x^ZZ')
            sage: P_QQ = GrowthGroup('x^QQ')
            sage: P_ZZ.gen() <= P_QQ.gen()^2  # indirect doctest
            True
        """
        return self.exponent < other.exponent

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this monomial growth element.

        INPUT:

        - ``rules`` -- dictionary; the neutral element of the group is replaced
          by the value to key ``'_one_'``

        OUTPUT: an object

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: G(x^42)._substitute_({'x': SR.var('z')})
            z^42
            sage: _.parent()
            Symbolic Ring
            sage: G(x^3)._substitute_({'x': 2})
            8
            sage: _.parent()
            Integer Ring
            sage: G(1 / x)._substitute_({'x': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in x^(-1) in Growth Group x^ZZ.
            > *previous* ZeroDivisionError: rational division by zero
            sage: G(1)._substitute_({'_one_': 'one'})
            'one'
        """
        if self.is_one():
            return rules['_one_']
        try:
            return self.parent()._var_._substitute_(rules) ** self.exponent
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)

    def _singularity_analysis_(self, var, zeta, precision):
        r"""
        Perform singularity analysis on this monomial growth element.

        INPUT:

        - ``var`` -- string denoting the variable

        - ``zeta`` -- a number

        - ``precision`` -- integer

        OUTPUT:

        An asymptotic expansion for `[z^n] f` where `n` is ``var``
        and `f` has this growth element as a singular expansion
        in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this element
        is a growth element in `T`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^QQ')
            sage: G(x^(1/2))._singularity_analysis_('n', 2, precision=2)
            1/sqrt(pi)*(1/2)^n*n^(-1/2) - 1/8/sqrt(pi)*(1/2)^n*n^(-3/2)
            + O((1/2)^n*n^(-5/2))
            sage: G = GrowthGroup('log(x)^QQ')
            sage: G(log(x))._singularity_analysis_('n', 1, precision=5)
            n^(-1) + O(n^(-3))
            sage: G(log(x)^2)._singularity_analysis_('n', 2, precision=3)
            2*(1/2)^n*n^(-1)*log(n) + 2*euler_gamma*(1/2)^n*n^(-1)
            + O((1/2)^n*n^(-2)*log(n)^2)

        TESTS::

            sage: G(log(x)^(1/2))._singularity_analysis_('n', 2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedError: singularity analysis of log(x)^(1/2)
            not implemented since exponent 1/2 is not an integer
            sage: G = GrowthGroup('log(log(x))^QQ')
            sage: G(log(log(x))^(1/2))._singularity_analysis_('n', 2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedError: singularity analysis of log(log(x))^(1/2)
            not implemented
        """
        from sage.rings.integer_ring import ZZ

        if self.parent()._var_.is_monomial():
            from sage.rings.asymptotic.asymptotic_expansion_generators import \
                asymptotic_expansions
            return asymptotic_expansions.SingularityAnalysis(
                var=var, zeta=zeta, alpha=self.exponent, beta=0, delta=0,
                precision=precision)
        elif self.parent().gens_logarithmic():
            if self.exponent not in ZZ:
                raise NotImplementedError(
                    'singularity analysis of {} not implemented '
                    'since exponent {} is not an integer'.format(
                        self, self.exponent))
            from sage.rings.asymptotic.asymptotic_expansion_generators import \
                asymptotic_expansions
            return asymptotic_expansions.SingularityAnalysis(
                var=var, zeta=zeta, alpha=0, beta=ZZ(self.exponent), delta=0,
                precision=precision, normalized=False)
        else:
            raise NotImplementedError(
                'singularity analysis of {} not implemented'.format(self))

    def _find_minimum_(self, valid_from):
        r"""
        Find the minimum of this growth element over the range implied by ``valid_from``.

        INPUT:

        - ``valid_from`` -- dictionary describing the range of the minimization:
          the keys are names of variables and the range is the intersection over
          the ranges where the absolute value of the variable designated by the
          key is at least the corresponding value

        OUTPUT: a number

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup, GrowthGroup

            sage: G = MonomialGrowthGroup(QQ, 'x')
            sage: G('x^3')._find_minimum_(valid_from={'x': 10})
            1000
            sage: e1 = G(raw_element=3); e2 = G(raw_element=2)
            sage: e3 = e2 / e1
            sage: e3
            x^(-1)
            sage: e3._find_minimum_(valid_from={'x': 5})
            Traceback (most recent call last):
            ...
            DecreasingGrowthElementError: the growth of x^(-1) is less than one
            sage: H = GrowthGroup('log(x)^ZZ')
            sage: l1 = H(raw_element=2)
            sage: l1._find_minimum_(valid_from={'x': 5})
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimum of log(x)^2 is not implemented
            sage: I = GrowthGroup('log(log(x))^ZZ')
            sage: l2 = I(raw_element=5)
            sage: l2._find_minimum_(valid_from={'x': 5})
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimum of log(log(x))^5 is not implemented
        """
        if not self.parent().gens_monomial():
            raise NotImplementedError(f'Minimum of {self} is not implemented')
        if self.is_lt_one():
            raise DecreasingGrowthElementError(self, f'the growth of {self} is less than one')
        elif self.is_one():
            return 1
        assert self.variable_names(), f'{self.variable_names()} is empty'
        return valid_from[self.variable_names()[0]] ** self.exponent


class MonomialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with powers of a fixed object/symbol.

    The elements :class:`MonomialGrowthElement` of this group represent powers
    of a fixed base; the group law is the multiplication, which corresponds
    to the addition of the exponents of the monomials.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``)

      As monomials are represented by this group, the elements in
      ``base`` are the exponents of these monomials.

    - ``var`` -- an object

      The string representation of ``var`` acts as a base of the
      monomials represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
        sage: P = MonomialGrowthGroup(ZZ, 'x'); P
        Growth Group x^ZZ
        sage: MonomialGrowthGroup(ZZ, log(SR.var('y')))
        Growth Group log(y)^ZZ

    .. SEEALSO::

        :class:`GenericGrowthGroup`

    TESTS::

        sage: L1 = MonomialGrowthGroup(QQ, log(x))
        sage: L2 = MonomialGrowthGroup(QQ, 'log(x)')
        sage: L1 is L2
        True
    """

    # enable the category framework for elements
    Element = MonomialGrowthElement

    # set everything up to determine category
    from sage.categories.sets_cat import Sets
    from sage.categories.posets import Posets
    from sage.categories.magmas import Magmas
    from sage.categories.additive_magmas import AdditiveMagmas

    _determine_category_subcategory_mapping_ = [
        (Sets(), Sets(), True),
        (Posets(), Posets(), False),
        (AdditiveMagmas(), Magmas(), False)]

    _determine_category_axiom_mapping_ = [
        ('AdditiveAssociative', 'Associative', False),
        ('AdditiveUnital', 'Unital', False),
        ('AdditiveInverse', 'Inverse', False),
        ('AdditiveCommutative', 'Commutative', False)]

    def _repr_short_(self):
        r"""
        A short representation string of this monomial growth group.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: MonomialGrowthGroup(ZZ, 'a')  # indirect doctest
            Growth Group a^ZZ


        TESTS::

            sage: MonomialGrowthGroup(ZZ, 'a')._repr_short_()
            'a^ZZ'
            sage: MonomialGrowthGroup(QQ, 'a')._repr_short_()
            'a^QQ'
            sage: MonomialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'a^QQ[x]'
        """
        from .misc import parent_to_repr_short, repr_op
        return repr_op(self._var_, '^', parent_to_repr_short(self.base()))

    def _convert_(self, data):
        r"""
        Convert ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('x')  # indirect doctest
            x

        ::

            sage: P(x)  # indirect doctest
            x
            sage: P(x^-333)  # indirect doctest
            x^(-333)
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: log(x)^2 is not in Growth Group x^ZZ.

        ::

            sage: PR.<x> = ZZ[]; x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: P(x^2)  # indirect doctest
            x^2

        ::

            sage: PSR.<x> = ZZ[[]]
            sage: P(x^42)  # indirect doctest
            x^42
            sage: P(x^12 + O(x^17))
            Traceback (most recent call last):
            ...
            ValueError: x^12 + O(x^17) is not in Growth Group x^ZZ.

        ::

            sage: R.<w,x> = ZZ[]
            sage: P(x^4242)  # indirect doctest
            x^4242
            sage: P(w^4242)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: w^4242 is not in Growth Group x^ZZ.

        ::

            sage: PSR.<w,x> = ZZ[[]]
            sage: P(x^7)  # indirect doctest
            x^7
            sage: P(w^7)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: w^7 is not in Growth Group x^ZZ.

        ::

            sage: P('x^7')
            x^7
            sage: P('1/x')
            x^(-1)
            sage: P('x^(-2)')
            x^(-2)
            sage: P('x^-2')
            x^(-2)

        ::

            sage: P('1')
            1

        ::

            sage: GrowthGroup('x^QQ')(GrowthGroup('x^ZZ')(1))
            1
        """
        if data == 1 or data == '1':
            return self.base().zero()
        var = repr(self._var_)
        if str(data) == var:
            return self.base().one()

        try:
            P = data.parent()
        except AttributeError:
            if var not in str(data):
                return  # this has to end here
            from sage.symbolic.ring import SR
            return self._convert_(SR(data))

        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_base import \
            MPolynomialRing_base
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        import operator
        if isinstance(P, sage.rings.abc.SymbolicRing):
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == var:
                    return exponent
        elif isinstance(P, (PolynomialRing_general, MPolynomialRing_base)):
            if data.is_monomial() and len(data.variables()) == 1:
                if var == str(data.variables()[0]):
                    return data.degree()
        elif isinstance(P, PowerSeriesRing_generic):
            if hasattr(data, 'variables') and len(data.variables()) == 1:
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    if var == str(data.variables()[0]):
                        return data.degree()
            elif len(P.variable_names()) == 1 and \
                            var == str(data.variable()[0]):
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    return data.degree()

    @staticmethod
    def _split_raw_element_(raw_element):
        r"""
        Split ``raw_element`` in a part convertible to this growth group
        and a part which needs to be converted by some compatible growth group.

        For this monomial growth group the two parts are
        real and imaginary part of ``raw_element``.

        INPUT:

        - ``raw_element`` -- an object

        OUTPUT: a pair of objects

        .. NOTE::

            This method is called by
            :meth:`~sage.rings.asymptotic.growth_group.PartialConversionElement.split`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('n^ZZ')
            sage: G._split_raw_element_(3 + 4*I)
            (3, 4)
        """
        from sage.functions.other import real, imag
        return real(raw_element), imag(raw_element)

    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth
        group.

        OUTPUT: a tuple containing elements of this growth group

        .. NOTE::

            A generator is called monomial generator if the variable
            of the underlying growth group is a valid identifier. For
            example, ``x^ZZ`` has ``x`` as a monomial generator,
            while ``log(x)^ZZ`` or ``icecream(x)^ZZ`` do not have
            monomial generators.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').gens_monomial()
            (x,)
            sage: GrowthGroup('log(x)^QQ').gens_monomial()
            ()
        """
        if not self._var_.is_monomial():
            return tuple()
        return (self(raw_element=self.base().one()),)

    def gens_logarithmic(self):
        r"""
        Return a tuple containing logarithmic generators of this growth
        group.

        OUTPUT: a tuple containing elements of this growth group

        .. NOTE::

            A generator is called logarithmic generator if the variable
            of the underlying growth group is the logarithm of a valid
            identifier. For
            example, ``x^ZZ`` has no logarithmic generator,
            while ``log(x)^ZZ`` has ``log(x)`` as
            logarithmic generator.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').gens_logarithmic()
            ()
            sage: GrowthGroup('log(x)^QQ').gens_logarithmic()
            (log(x),)
        """
        if str(self.gen()) == "log({})".format(self.variable_name()):
            return (self(raw_element=self.base().one()),)
        else:
            return tuple()

    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is a
        :class:`monomial construction functor <MonomialGrowthGroupFunctor>`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').construction()
            (MonomialGrowthGroup[x], Integer Ring)
        """
        return MonomialGrowthGroupFunctor(self._var_), self.base()

    @classmethod
    def factory(cls,
                base, var,
                extend_by_non_growth_group=False,
                return_factors=False,
                **kwds):
        r"""
        Create a monomial growth group.

        INPUT:

        - ``base``, ``var``, keywords -- use in the initialization of the
          exponential growth group; see :class:`MonomialGrowthGroup`
          for details.

        - ``extend_by_non_growth_group`` -- boolean (default: ``False``); if
          set, then the growth group consists of two parts, one part dealing
          with the absolute values of the bases and one for their arguments.

        - ``return_factors`` -- boolean (default: ``False``); if set,
          then a tuple of the (cartesian) factors of this growth group
          is returned

        OUTPUT: a growth group or tuple of growth groups

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: MonomialGrowthGroup.factory(ZZ, 'n')
            Growth Group n^ZZ
            sage: MonomialGrowthGroup.factory(ImaginaryGroup(ZZ), 'n')
            Growth Group n^(ZZ*I)

        TESTS::

            sage: MonomialGrowthGroup.factory(ZZ, 'n', return_factors=True)
            (Growth Group n^ZZ,)
            sage: MonomialGrowthGroup.factory(ZZ, 'n', extend_by_non_growth_group=True)
            Growth Group n^ZZ * n^(ZZ*I)
            sage: MonomialGrowthGroup.factory(ZZ, 'n', return_factors=True,
            ....:                             extend_by_non_growth_group=True)
            (Growth Group n^ZZ, Growth Group n^(ZZ*I))
        """
        from sage.categories.cartesian_product import cartesian_product
        from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup

        if isinstance(base, ImaginaryGroup):
            groups = (cls._non_growth_group_class_(base, var, **kwds),)
        elif extend_by_non_growth_group:
            M = cls(base, var, **kwds)
            groups = (M, M.non_growth_group())
        else:
            groups = (cls(base, var, **kwds),)

        if return_factors:
            return tuple(groups)
        else:
            return cartesian_product(groups)

    def non_growth_group(self):
        r"""
        Return a non-growth group
        (with an imaginary group as base)
        compatible with this monomial growth group.

        OUTPUT: a group group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('n^ZZ').non_growth_group()
            Growth Group n^(ZZ*I)
        """
        from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
        J = ImaginaryGroup(self.base())
        return self._non_growth_group_class_(J, self._var_)


class MonomialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`monomial growth groups <MonomialGrowthGroup>`.

    INPUT:

    - ``var`` -- string or list of strings (or anything else
      :class:`Variable` accepts)

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, MonomialGrowthGroupFunctor
        sage: GrowthGroup('z^QQ').construction()[0]
        MonomialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AbstractGrowthGroupFunctor`,
        :class:`ExponentialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, MonomialGrowthGroupFunctor
        sage: cm = sage.structure.element.get_coercion_model()
        sage: A = GrowthGroup('x^QQ')
        sage: B = MonomialGrowthGroupFunctor('x')(ZZ['t'])
        sage: cm.common_parent(A, B)
        Growth Group x^QQ[t]
    """

    _functor_name = 'MonomialGrowthGroup'

    def __init__(self, var):
        r"""
        See :class:`MonomialGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroupFunctor
            sage: MonomialGrowthGroupFunctor('x')
            MonomialGrowthGroup[x]
        """
        from sage.categories.commutative_additive_monoids import CommutativeAdditiveMonoids
        super().__init__(var, CommutativeAdditiveMonoids())

    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` -- anything :class:`MonomialGrowthGroup` accepts

        OUTPUT: a monomial growth group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('z^QQ').construction()
            sage: F(R)  # indirect doctest
            Growth Group z^QQ
        """
        return MonomialGrowthGroup(base, self.var)


class ExponentialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of exponential growth elements.

    INPUT:

    - ``parent`` -- an :class:`ExponentialGrowthGroup`

    - ``raw_element`` -- an element from the base ring of the parent

      This ``raw_element`` is the base of the created exponential
      growth element.

    An exponential growth element represents a term of the type
    `\operatorname{base}^{\operatorname{variable}}`. The multiplication
    corresponds to the multiplication of the bases.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: P = GrowthGroup('(ZZ_+)^x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        2^x
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P(1)) and P.le(P(1), e2)
        True
    """

    def _check_(self):
        r"""
        Perform an additional check at the end of :meth:`__init__`.

        This check is whether the base of this
        exponential growth group is positive.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('(QQ_+)^x')
            sage: P(raw_element=2/3)  # indirect doctest
            (2/3)^x
            sage: P(raw_element=-2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            PartialConversionValueError: base -2/3 (Rational Field) must be positive
        """
        if not self.base > 0:
            from sage.structure.element import parent
            raise PartialConversionValueError(
                PartialConversionElement(self.parent(), self.base),
                'base {} ({}) must be positive'.format(self.base, parent(self.base)))

    @property
    def base(self):
        r"""
        The base of this exponential growth element.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('(ZZ_+)^x')
            sage: P(42^x).base
            42
        """
        return self._raw_element_

    def _repr_(self, latex=False):
        r"""
        A representation string for this exponential growth element.

        INPUT:

        - ``latex`` -- boolean (default: ``False``); if set, then
          LaTeX-output is returned

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('QQ^x')
            sage: P(1)._repr_()
            '1'
            sage: P(5^x)  # indirect doctest
            5^x
            sage: P((1/2)^x)  # indirect doctest
            (1/2)^x

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'x')
            sage: UU((-1)^x)  # indirect doctest
            (-1)^x

        ::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: G = ExponentialGrowthGroup(ZZ['x'], 'y'); G
            Growth Group ZZ[x]^y
            sage: G('(1-x)^y')
            (-x + 1)^y
            sage: G('(1+x)^y')
            (x + 1)^y
        """
        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
        else:
            f = repr

        from .misc import repr_op

        var = f(self.parent()._var_)
        if self.base.is_one():
            return '1'
        if latex:
            return repr_op(latex_repr(self.base)._latex_(), '^', latex=True) + \
                '{' + latex_repr(var)._latex_() + '}'
        else:
            return repr_op(str(self.base), '^', var)

    def _latex_(self):
        r"""
        A LaTeX-representation string for this exponential growth element.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('QQ^x')
            sage: latex(P(1))
            1
            sage: latex(P(5^x))  # indirect doctest
            5^{x}
            sage: latex(P((1/2)^x))  # indirect doctest
            \left(\frac{1}{2}\right)^{x}

        ::

            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'x')
            sage: latex(UU((-1)^x))  # indirect doctest
            \left(-1\right)^{x}

        ::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: G = ExponentialGrowthGroup(ZZ['x'], 'y'); G
            Growth Group ZZ[x]^y
            sage: latex(G('(1-x)^y'))
            \left(-x + 1\right)^{y}
            sage: latex(G('(1+x)^y'))
            \left(x + 1\right)^{y}
        """
        return self._repr_(latex=True)

    def _mul_(self, other):
        r"""
        Multiply this exponential growth element with another.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`

        OUTPUT: the product as a :class:`ExponentialGrowthElement`

        .. NOTE::

            Two exponential growth elements are multiplied by
            multiplying their bases.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('ZZ^x')
            sage: a = P(2^x)
            sage: b = P(3^x)
            sage: c = a._mul_(b); c
            6^x
            sage: c == a*b
            True
            sage: a*b*a  # indirect doctest
            12^x
        """
        return self.parent()(raw_element=self.base * other.base)

    def __invert__(self):
        r"""
        Return the multiplicative inverse of this exponential growth element.

        OUTPUT: the multiplicative inverse as a :class:`ExponentialGrowthElement`

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('(ZZ_+)^x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            (1/2)^x
            sage: e2 == ~e1
            True
            sage: e2.parent()
            Growth Group QQ^x

        ::

            sage: (~P(raw_element=1)).parent()
            Growth Group QQ^x

        ::

            sage: UU = GrowthGroup('UU^n')
            sage: zeta = UU.an_element(); zeta
            (-1)^n
            sage: ~zeta
            (-1)^n
        """
        return self.parent()._create_element_in_extension_(~self.base)

    def __pow__(self, exponent):
        r"""
        Calculate the power of this growth element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- a number. This can be anything that is valid to be
          on the right hand side of ``*`` with an elements of the
          parent's base.

        OUTPUT: the result of this exponentiation as an :class:`ExponentialGrowthElement`

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('(ZZ_+)^x')
            sage: assume(SR.an_element() > 0)
            sage: a = P(7^x); a
            7^x
            sage: b = a^(1/2); b
            sqrt(7)^x
            sage: b.parent()
            Growth Group SR^x
            sage: b^12
            117649^x
            sage: forget()

        TESTS::

             sage: SCR = SR.subring(no_variables=True)
             sage: G = GrowthGroup('QQ^x * x^ZZ'); G
             Growth Group QQ^x * x^ZZ * Signs^x
             sage: x = G('x')
             sage: x^SCR(1)
             x
             sage: _.parent()
             Growth Group QQ^x * x^ZZ * Signs^x
        """
        from .misc import strip_symbolic
        return self.parent()._create_element_in_extension_(
            self.base ** strip_symbolic(exponent))

    def _log_factor_(self, base=None, locals=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second is a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('(QQ_+)^x')
            sage: G('4^x').log_factor(base=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot build log(4^x, base=2) since x is not in
            Growth Group QQ^x.

        ::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: L.<log2> = ZZ[]
            sage: G = ExponentialGrowthGroup(L, 'x')
            sage: G(raw_element=2)._log_factor_()
            (('x', log(2)),)
            sage: G(raw_element=2)._log_factor_(
            ....:       locals={'log': lambda z, base: log2 if z == 2 else log(z)})
            (('x', log2),)
        """
        if self.is_one():
            return tuple()
        b = self.base
        if base is None and hasattr(b, 'is_monomial') and b.is_monomial() and \
                        b.variable_name() == 'e':
            coefficient = b.valuation()
        elif base is None and str(b) == 'e':
            coefficient = self.parent().base().one()
        else:
            log = self.parent().locals(locals)['log']
            coefficient = log(b, base=base)

        return ((str(self.parent()._var_), coefficient),)

    def _lt_(self, other):
        r"""
        Return whether this :class:`ExponentialGrowthElement` is
        less than ``other``.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`

        OUTPUT: boolean

        .. NOTE::

            This function compares two instances of
            :class:`ExponentialGrowthElement`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: assume(SR.an_element() > 0)
            sage: P_ZZ = GrowthGroup('(ZZ_+)^x')
            sage: P_SR = GrowthGroup('(SR_+)^x')
            sage: P_ZZ(2^x) <= P_SR(sqrt(3)^x)^2  # indirect doctest
            True
            sage: forget()

        Check that :issue:`19999` is fixed::

            sage: P_ZZ_UU = GrowthGroup('ZZ^x * UU^x')
            sage: P_ZZ_UU((-2)^x) <= P_ZZ_UU(2^x) or P_ZZ_UU(2^x) <= P_ZZ_UU((-2)^x)
            False
        """
        return bool(self.base < other.base)

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this exponential growth element.

        INPUT:

        - ``rules`` -- dictionary; the neutral element of the group is replaced
          by the value to key ``'_one_'``

        OUTPUT: an object

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('(QQ_+)^x')
            sage: G((1/2)^x)._substitute_({'x': SR.var('z')})
            (1/2)^z
            sage: _.parent()
            Symbolic Ring
            sage: G((1/2)^x)._substitute_({'x': 2})
            1/4
            sage: _.parent()
            Rational Field
            sage: G(1)._substitute_({'_one_': 'one'})
            'one'
        """
        if self.is_one():
            return rules['_one_']
        try:
            return self.base ** self.parent()._var_._substitute_(rules)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)


class ExponentialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with expressions involving a fixed
    variable/symbol as the exponent.

    The elements :class:`ExponentialGrowthElement` of this group
    represent exponential functions with bases from a fixed base
    ring; the group law is the multiplication.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``)

      As exponential expressions are represented by this group,
      the elements in ``base`` are the bases of these exponentials.

    - ``var`` -- an object

      The string representation of ``var`` acts as an exponent of the
      elements represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
        sage: P = ExponentialGrowthGroup(QQ, 'x'); P
        Growth Group QQ^x

    .. SEEALSO::

        :class:`GenericGrowthGroup`
    """

    # enable the category framework for elements
    Element = ExponentialGrowthElement

    # set everything up to determine category
    from sage.categories.sets_cat import Sets
    from sage.categories.posets import Posets
    from sage.categories.magmas import Magmas
    from sage.categories.groups import Groups
    from sage.categories.division_rings import DivisionRings

    _determine_category_subcategory_mapping_ = [
        (Sets(), Sets(), True),
        (Posets(), Posets(), False),
        (Magmas(), Magmas(), False),
        (DivisionRings(), Groups(), False)]

    _determine_category_axiom_mapping_ = [
        ('Associative', 'Associative', False),
        ('Unital', 'Unital', False),
        ('Inverse', 'Inverse', False),
        ('Commutative', 'Commutative', False)]

    def __init__(self, base, *args, **kwds):
        r"""
        See :class:`ExponentialGrowthGroup` for more information.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(SR.subring(no_variables=True), 't')  # indirect doctest
            Growth Group (Symbolic Constants Subring)^t
            sage: ExponentialGrowthGroup(SR, 't')  # indirect doctest
            doctest:warning
            ...
            RuntimeWarning: When using the Exponential Growth Group SR^t,
            make assumptions on the used symbolic elements.
            In particular, use something like 'assume(SR.an_element() > 0)'
            to make coercions work properly.
            Growth Group SR^t
            sage: assume(SR.an_element() > 0)
            sage: ExponentialGrowthGroup(SR, 't')  # indirect doctest
            Growth Group SR^t
            sage: forget()
        """
        from warnings import warn

        super().__init__(base, *args, **kwds)
        if isinstance(base, sage.rings.abc.SymbolicRing) and not self._an_element_base_() > 0:
            warn("When using the Exponential {}, make "
                 "assumptions on the used symbolic elements.\n"
                 "In particular, use something like "
                 "'assume(SR.an_element() > 0)' to make "
                 "coercions work properly.".format(self),
                 RuntimeWarning, 2)

    def _repr_short_(self):
        r"""
        A short representation string of this exponential growth group.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(QQ, 'a')  # indirect doctest
            Growth Group QQ^a


        TESTS::

            sage: ExponentialGrowthGroup(QQ, 'a')._repr_short_()
            'QQ^a'
            sage: ExponentialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'QQ[x]^a'
        """
        from .misc import parent_to_repr_short, repr_op
        return repr_op(parent_to_repr_short(self.base()), '^', self._var_)

    def _convert_(self, data):
        r"""
        Convert given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('(QQ_+)^x')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('2^x')  # indirect doctest
            2^x

        ::

            sage: P(2^x)  # indirect doctest
            2^x
            sage: P((-333)^x)  # indirect doctest
            Traceback (most recent call last):
            ...
            PartialConversionValueError: base -333 (Rational Field) must be positive
            sage: P(0)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 0 is not in Growth Group QQ^x.

        ::

            sage: P('7^x')
            7^x
            sage: P('(-2)^x')
            Traceback (most recent call last):
            ...
            PartialConversionValueError: base -2 (Rational Field) must be positive

        ::

            sage: P = GrowthGroup('SR^x')
            sage: P(sqrt(3)^x)
            sqrt(3)^x
            sage: P((3^(1/3))^x)
            (3^(1/3))^x
            sage: P(e^x)
            e^x
            sage: P(exp(2*x))
            (e^2)^x

        ::

            sage: GrowthGroup('QQ^x')(GrowthGroup('ZZ^x')(1))
            1

        ::

            sage: E = GrowthGroup('(QQ_+)^x * UU^x')
            sage: E((-333)^x)  # indirect doctest
            333^x*(-1)^x
            sage: E('(-2)^x')
            2^x*(-1)^x
        """
        if data == '1' or isinstance(data, int) and data == 1:
            return self.base().one()
        var = repr(self._var_)
        try:
            P = data.parent()
        except AttributeError:
            if data == 1:
                return self.base().one()
            s = str(data)
            if var not in s:
                return  # this has to end here

            elif s.endswith('^' + var):
                return self.base()(s.replace('^' + var, '')
                                   .replace('(', '').replace(')', ''))
            else:
                return  # end of parsing

        import operator
        from sage.functions.log import Function_exp
        from sage.symbolic.operators import mul_vararg

        if isinstance(P, sage.rings.abc.SymbolicRing):
            op = data.operator()
            if op == operator.pow:
                base, exponent = data.operands()
                if str(exponent) == var:
                    return base
                elif exponent.operator() == mul_vararg:
                    return base ** (exponent / P(var))
            elif isinstance(op, Function_exp):
                from sage.functions.log import exp
                base = exp(1)
                exponent = data.operands()[0]
                if str(exponent) == var:
                    return base
                elif exponent.operator() == mul_vararg:
                    return base ** (exponent / P(var))

        elif data == 1:  # can be expensive, so let's put it at the end
            return self.base().one()

    @staticmethod
    def _split_raw_element_(base):
        r"""
        Split ``raw_element`` in a part convertible to this growth group
        and a part which needs to be converted by some compatible growth group.

        For this exponential growth group the two parts are
        absolute value and argument of ``raw_element``.

        INPUT:

        - ``raw_element`` -- an object

        OUTPUT: a pair of objects

        .. NOTE::

            This method is called by
            :meth:`~sage.rings.asymptotic.growth_group.PartialConversionElement.split`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('(QQ_+)^n')
            sage: G._split_raw_element_(ZZ(-2))
            (2, -1)
            sage: G._split_raw_element_(ZZ(2))
            (2, None)
            sage: G._split_raw_element_(QQ(-2/3))
            (2/3, -1)
            sage: G._split_raw_element_(QQ(2/3))
            (2/3, None)
            sage: G._split_raw_element_(AA(-3/4))
            (3/4, -1)
            sage: G._split_raw_element_(AA(3/4))
            (3/4, None)
            sage: G._split_raw_element_(RR(-3.14))
            (3.14000000000000, -1)
            sage: G._split_raw_element_(RR(3.14))
            (3.14000000000000, None)
           sage: G._split_raw_element_(RIF(-3.14))
            (3.1400000000000002?, -1)
            sage: G._split_raw_element_(RIF(3.14))
            (3.1400000000000002?, None)
            sage: G._split_raw_element_(RBF(-3.14))
            ([3.140000000000000 +/- 1.25e-16], -1)
            sage: G._split_raw_element_(RBF(3.14))
            ([3.140000000000000 +/- 1.25e-16], None)
            sage: G._split_raw_element_(CC(-3.14))
            (3.14000000000000, -1.00000000000000)
            sage: G._split_raw_element_(CC(3.14))
            (3.14000000000000, 1.00000000000000)
            sage: G._split_raw_element_(CC(1+I))
            (1.41421356237310, 0.707106781186547 + 0.707106781186547*I)
            sage: G._split_raw_element_(CC(I))
            (1.00000000000000, 1.00000000000000*I)
            sage: G._split_raw_element_(CIF(1+I))
            (1.414213562373095?, 0.707106781186548? + 0.707106781186548?*I)
            sage: G._split_raw_element_(CBF(1+I))
            ([1.414213562373095 +/- 2.99e-16],
             [0.707106781186548 +/- 6.50e-16] + [0.707106781186548 +/- 6.50e-16]*I)

            sage: G._split_raw_element_(SR(-2/3))
            (2/3, -1)

            sage: G._split_raw_element_(x)
            Traceback (most recent call last):
            ...
            ValueError: cannot split x (Symbolic Ring) into abs and arg
            sage: assume(x > 0)
            sage: G._split_raw_element_(x)
            (x, None)
            sage: forget()
            sage: assume(x < 0)
            sage: G._split_raw_element_(x)
            (-x, -1)
            sage: forget()
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.qqbar import AA
        from sage.structure.element import parent

        P = base.parent()
        if isinstance(P, sage.rings.abc.SymbolicRing):
            try:
                base = base.pyobject()
            except TypeError:
                pass
            else:
                P = base.parent()

        if P in (ZZ, QQ, AA) or isinstance(P, (sage.rings.abc.SymbolicRing,
                                               sage.rings.abc.RealField,
                                               sage.rings.abc.RealIntervalField,
                                               sage.rings.abc.RealBallField)):
            if base > 0:
                return base, None
            if base < 0:
                return -base, -1
        elif isinstance(P, (sage.rings.abc.ComplexField,
                            sage.rings.abc.ComplexIntervalField,
                            sage.rings.abc.ComplexBallField)):
            size = abs(base)
            direction = base / size
            return size, direction

        raise ValueError('cannot split {} ({}) into '
                         'abs and arg'.format(base, parent(base)))

    def _an_element_(self):
        r"""
        Return an element of this exponential growth group.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(SR, 'n').an_element()  # indirect doctest
            Traceback (most recent call last):
            ...
            PartialConversionValueError: base abs(some_variable) (Symbolic Ring) must be positive

            sage: assume(SR.an_element() > 0)
            sage: ExponentialGrowthGroup(SR, 'n').an_element()  # indirect doctest
            some_variable^n
            sage: forget()
            sage: ExponentialGrowthGroup(SR.subring(no_variables=True), 'n').an_element()  # indirect doctest
            (pi*e)^n
        """
        return self.element_class(self, self._an_element_base_())

    def _an_element_base_(self):
        r"""
        Return a base for :meth:`_an_element_` of this exponential growth group.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup(SR, 'n')._an_element_base_()
            abs(some_variable)
            sage: ExponentialGrowthGroup(SR.subring(no_variables=True), 'n')._an_element_base_()
            pi*e
        """
        e = self.base().an_element()
        return e if e > 0 else abs(e)

    def some_elements(self):
        r"""
        Return some elements of this exponential growth group.

        See :class:`TestSuite` for a typical use case.

        OUTPUT: an iterator

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: tuple(GrowthGroup('(QQ_+)^z').some_elements())
            ((1/2)^z, 2^z, 1, 42^z, (2/3)^z, (3/2)^z, ...)
        """
        return iter(self.element_class(self, e)
                    for e in self.base().some_elements() if e > 0)

    def gens(self):
        r"""
        Return a tuple of all generators of this exponential growth
        group.

        OUTPUT: an empty tuple

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: E = GrowthGroup('(ZZ_+)^x')
            sage: E.gens()
            ()
        """
        return tuple()

    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is an
        :class:`exponential construction functor <ExponentialGrowthGroupFunctor>`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('(QQ_+)^x').construction()
            (ExponentialGrowthGroup[x], Rational Field)
        """
        return ExponentialGrowthGroupFunctor(self._var_), self.base()

    @classmethod
    def factory(cls,
                base, var,
                extend_by_non_growth_group=True,
                return_factors=False,
                **kwds):
        r"""
        Create an exponential growth group.

        This factory takes care of the splitting of the bases into their
        absolute values and arguments.

        INPUT:

        - ``base``, ``var``, keywords -- use in the initialization of the
          exponential growth group; see :class:`ExponentialGrowthGroup`
          for details.

        - ``extend_by_non_growth_group`` -- boolean (default: ``True``); if
          set, then the growth group consists of two parts, one part dealing
          with the absolute values of the bases and one for their arguments.

        - ``return_factors`` -- boolean (default: ``False``); if set,
          then a tuple of the (cartesian) factors of this growth group
          is returned

        OUTPUT: a growth group or tuple of growth groups

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            sage: ExponentialGrowthGroup.factory(QQ, 'n')
            Growth Group QQ^n * Signs^n

        TESTS::

            sage: ExponentialGrowthGroup.factory(QQ, 'n', return_factors=True)
            (Growth Group QQ^n, Growth Group Signs^n)
            sage: ExponentialGrowthGroup.factory(QQ, 'n', extend_by_non_growth_group=False)
            Growth Group QQ^n
            sage: from sage.groups.misc_gps.argument_groups import ArgumentGroup
            sage: UU = ArgumentGroup(exponents=QQ)
            sage: ExponentialGrowthGroup.factory(UU, 'n')
            Growth Group UU^n

            sage: ExponentialGrowthGroup.factory(CC, 'n')
            Growth Group RR^n * UU_RR^n
            sage: ExponentialGrowthGroup.factory(CyclotomicField(3), 'n')
            Growth Group (Algebraic Real Field)^n * (Arg_(Cyclotomic Field of order 3 and degree 2))^n
        """
        from sage.categories.cartesian_product import cartesian_product
        from sage.groups.misc_gps.argument_groups import AbstractArgumentGroup
        from sage.groups.misc_gps.argument_groups import ArgumentGroup
        from sage.rings.number_field.number_field import NumberField_cyclotomic
        from sage.rings.qqbar import QQbar, AA

        if isinstance(base, AbstractArgumentGroup):
            groups = (cls._non_growth_group_class_(base, var, **kwds),)
        elif extend_by_non_growth_group:
            if base == QQbar or isinstance(base, NumberField_cyclotomic):
                EE = cls(AA, var, **kwds)
                UU = cls._non_growth_group_class_(
                    ArgumentGroup(domain=base), var)
                groups = (EE, UU)
            elif isinstance(base, (sage.rings.abc.ComplexField,
                                   sage.rings.abc.ComplexIntervalField,
                                   sage.rings.abc.ComplexBallField)):
                EE = cls(base._real_field(), var, **kwds)
                UU = cls._non_growth_group_class_(
                    ArgumentGroup(exponents=base._real_field()), var)
                groups = (EE, UU)
            else:
                EE = cls(base, var, **kwds)
                groups = (EE, EE.non_growth_group())
        else:
            groups = (cls(base, var, **kwds),)

        if return_factors:
            return tuple(groups)
        else:
            return cartesian_product(groups)

    def non_growth_group(self):
        r"""
        Return a non-growth group
        (with an argument group, e.g. roots of unity, as base)
        compatible with this exponential growth group.

        OUTPUT: a group group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('(QQ_+)^x').non_growth_group()
            Growth Group Signs^x
            sage: GrowthGroup('(RR_+)^x').non_growth_group()
            Growth Group Signs^x
            sage: GrowthGroup('(RIF_+)^x').non_growth_group()
            Growth Group Signs^x
            sage: GrowthGroup('(RBF_+)^x').non_growth_group()
            Growth Group Signs^x
            sage: GrowthGroup('(CC_+)^x').non_growth_group()
            Growth Group UU_RR^x
            sage: GrowthGroup('(CIF_+)^x').non_growth_group()
            Growth Group UU_RIF^x
            sage: GrowthGroup('(CBF_+)^x').non_growth_group()
            Growth Group UU_RBF^x
        """
        from sage.groups.misc_gps.argument_groups import ArgumentGroup
        UU = ArgumentGroup(domain=self.base())
        return self._non_growth_group_class_(UU, self._var_)


class ExponentialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`exponential growth groups <ExponentialGrowthGroup>`.

    INPUT:

    - ``var`` -- string or list of strings (or anything else
      :class:`Variable` accepts)

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, ExponentialGrowthGroupFunctor
        sage: GrowthGroup('(QQ_+)^z').construction()[0]
        ExponentialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AbstractGrowthGroupFunctor`,
        :class:`MonomialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, ExponentialGrowthGroupFunctor
        sage: cm = sage.structure.element.get_coercion_model()
        sage: A = GrowthGroup('(QQ_+)^x')
        sage: B = ExponentialGrowthGroupFunctor('x')(ZZ['t'])
        sage: cm.common_parent(A, B)
        Growth Group QQ[t]^x
    """

    _functor_name = 'ExponentialGrowthGroup'

    def __init__(self, var):
        r"""
        See :class:`ExponentialGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroupFunctor
            sage: ExponentialGrowthGroupFunctor('x')
            ExponentialGrowthGroup[x]
        """
        from sage.categories.monoids import Monoids
        super().__init__(var, Monoids())

    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` -- anything :class:`ExponentialGrowthGroup` accepts

        OUTPUT: an exponential growth group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('(QQ_+)^z').construction()
            sage: F(R)  # indirect doctest
            Growth Group QQ^z
        """
        return ExponentialGrowthGroup(base, self.var)


class GenericNonGrowthElement(GenericGrowthElement):
    r"""
    An element of :class:`GenericNonGrowthGroup`.
    """

    def _lt_(self, other):
        r"""
        Return ``False`` as elements are not comparable.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'n')
            sage: UU(raw_element=-1) < UU(raw_element=1)
            False
            sage: UU(raw_element=-1) > UU(raw_element=1)
            False

            sage: from sage.rings.asymptotic.growth_group import MonomialNonGrowthGroup
            sage: MM = MonomialNonGrowthGroup(RootsOfUnityGroup(), 'n')
            sage: MM(raw_element=-1) < MM(raw_element=1)
            False
            sage: MM(raw_element=-1) > MM(raw_element=1)
            False
        """
        return False


class GenericNonGrowthGroup(GenericGrowthGroup):
    r"""
    A (abstract) growth group whose elements are all of the same growth `1`.

    See :class:`ExponentialNonGrowthGroup` for a concrete
    realization.
    """

    @staticmethod
    def _initial_category_(base):
        r"""
        Return a category with which creating the actual category
        of this growth group starts.

        INPUT:

        - ``base`` -- a SageMath parent

        OUTPUT: always the category of posets

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: ExponentialNonGrowthGroup._initial_category_(ZZ)
            Category of posets
            sage: ExponentialNonGrowthGroup._initial_category_(QQ)
            Category of posets
            sage: ExponentialNonGrowthGroup._initial_category_(SR)
            Category of posets
        """
        from sage.categories.posets import Posets
        return Posets()


class ExponentialNonGrowthElement(GenericNonGrowthElement,
                                  ExponentialGrowthElement):
    r"""
    An element of :class:`ExponentialNonGrowthGroup`.
    """

    def _check_(self):
        r"""
        Perform an additional check at the end of :meth:`__init__`.

        No check is performed for this class.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: UU = RootsOfUnityGroup()
            sage: EE = ExponentialNonGrowthGroup(UU, 'n')
            sage: EE(raw_element=UU(-1))  # indirect doctest
            (-1)^n
        """
        pass


class ExponentialNonGrowthGroup(GenericNonGrowthGroup,
                                ExponentialGrowthGroup):
    r"""
    A growth group whose base is an
    :mod:`argument group <sage.groups.misc_gps.argument_groups>`.

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
        sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
        sage: UU = ExponentialNonGrowthGroup(RootsOfUnityGroup(), 'n')
        sage: UU(raw_element=-1)
        (-1)^n

    TESTS::

        sage: UU(raw_element=int(-1))
        (-1)^n

    ::

        sage: UU.category()
        Join of Category of commutative groups and Category of posets
    """

    Element = ExponentialNonGrowthElement

    def _an_element_base_(self):
        r"""
        Return a base for :meth:`_an_element_` of this exponential non growth group.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: from sage.rings.asymptotic.growth_group import ExponentialNonGrowthGroup
            sage: ExponentialNonGrowthGroup(SignGroup(), 'n').an_element()  # indirect doctest
            (-1)^n
        """
        return self.base().an_element()

    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is an
        :class:`ExponentialNonGrowthGroupFunctor`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('UU^x').construction()
            (ExponentialNonGrowthGroup[x], Group of Roots of Unity)
        """
        return ExponentialNonGrowthGroupFunctor(self._var_), self.base()


ExponentialGrowthGroup._non_growth_group_class_ = ExponentialNonGrowthGroup


class ExponentialNonGrowthGroupFunctor(ExponentialGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`ExponentialNonGrowthGroup`.
    """

    _functor_name = 'ExponentialNonGrowthGroup'

    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` -- anything :class:`ExponentialNonGrowthGroup` accepts

        OUTPUT: an exponential argument growth group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('UU^z').construction()
            sage: F(R)  # indirect doctest
            Growth Group UU^z
        """
        return ExponentialNonGrowthGroup(base, self.var)


class MonomialNonGrowthElement(GenericNonGrowthElement,
                               MonomialGrowthElement):
    r"""
    An element of :class:`MonomialNonGrowthGroup`.
    """
    pass


class MonomialNonGrowthGroup(GenericNonGrowthGroup,
                             MonomialGrowthGroup):
    r"""
    A growth group whose base is an
    :mod:`imaginary group <sage.groups.misc_gps.imaginary_groups>`.

    EXAMPLES::

        sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
        sage: from sage.rings.asymptotic.growth_group import MonomialNonGrowthGroup
        sage: J = MonomialNonGrowthGroup(ImaginaryGroup(ZZ), 'n')
        sage: J.an_element()
        n^I

    TESTS::

        sage: J.category()
        Join of Category of commutative groups and Category of posets
    """

    Element = MonomialNonGrowthElement

    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is an
        :class:`MonomialNonGrowthGroupFunctor`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^(QQ*I)').construction()
            (MonomialNonGrowthGroup[x], Imaginary Group over Rational Field)
        """
        return MonomialNonGrowthGroupFunctor(self._var_), self.base()


MonomialGrowthGroup._non_growth_group_class_ = MonomialNonGrowthGroup


class MonomialNonGrowthGroupFunctor(MonomialGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`MonomialNonGrowthGroup`.
    """

    _functor_name = 'MonomialNonGrowthGroup'

    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` -- anything :class:`MonomialNonGrowthGroup` accepts

        OUTPUT: an exponential argument growth group

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('n^(ZZ*I)').construction()
            sage: F(R)  # indirect doctest
            Growth Group n^(ZZ*I)
        """
        return MonomialNonGrowthGroup(base, self.var)


GrowthGroupFactor = namedtuple('GrowthGroupFactor',
                               ['cls', 'base', 'var',
                                'extend_by_non_growth_group'])


class GrowthGroupFactory(UniqueFactory):
    r"""
    A factory creating asymptotic growth groups.

    INPUT:

    - ``specification`` -- string

    - keyword arguments are passed on to the growth group
      constructor.
      If the keyword ``ignore_variables`` is not specified, then
      ``ignore_variables=('e',)`` (to ignore ``e`` as a variable name)
      is used.

    OUTPUT: an asymptotic growth group

    .. NOTE::

        An instance of this factory is available as ``GrowthGroup``.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: GrowthGroup('x^ZZ')
        Growth Group x^ZZ
        sage: GrowthGroup('log(x)^QQ')
        Growth Group log(x)^QQ

    This factory can also be used to construct Cartesian products
    of growth groups::

        sage: GrowthGroup('x^ZZ * y^ZZ')
        Growth Group x^ZZ * y^ZZ
        sage: GrowthGroup('x^ZZ * log(x)^ZZ')
        Growth Group x^ZZ * log(x)^ZZ
        sage: GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ')
        Growth Group x^ZZ * log(x)^ZZ * y^QQ
        sage: GrowthGroup('(QQ_+)^x * x^ZZ * y^QQ * (QQ_+)^z')
        Growth Group QQ^x * x^ZZ * y^QQ * QQ^z
        sage: GrowthGroup('QQ^x * x^ZZ * y^QQ * QQ^z')
        Growth Group QQ^x * x^ZZ * Signs^x * y^QQ * QQ^z * Signs^z
        sage: GrowthGroup('exp(x)^ZZ * x^ZZ')
        Growth Group exp(x)^ZZ * x^ZZ
        sage: GrowthGroup('(e^x)^ZZ * x^ZZ')
        Growth Group (e^x)^ZZ * x^ZZ

    ::

        sage: GrowthGroup('QQ^n * n^ZZ')
        Growth Group QQ^n * n^ZZ * Signs^n
        sage: GrowthGroup('(QQ_+)^n * n^ZZ * UU^n')
        Growth Group QQ^n * n^ZZ * UU^n
        sage: GrowthGroup('(QQ_+)^n * n^ZZ')
        Growth Group QQ^n * n^ZZ

    ::

        sage: GrowthGroup('n^(ZZ)')
        Growth Group n^ZZ
        sage: GrowthGroup('n^(ZZ[I])')
        Growth Group n^ZZ * n^(ZZ*I)
        sage: GrowthGroup('n^(I*ZZ)')
        Growth Group n^(ZZ*I)
        sage: GrowthGroup('n^(ZZ*I)')
        Growth Group n^(ZZ*I)

    TESTS::

        sage: G = GrowthGroup('(e^(n*log(n)))^ZZ')
        sage: G, G._var_
        (Growth Group (e^(n*log(n)))^ZZ, e^(n*log(n)))
        sage: G = GrowthGroup('(e^n)^ZZ')
        sage: G, G._var_
        (Growth Group (e^n)^ZZ, e^n)
        sage: G = GrowthGroup('(e^(n*log(n)))^ZZ * (e^n)^ZZ * n^ZZ * log(n)^ZZ')
        sage: G, tuple(F._var_ for F in G.cartesian_factors())
        (Growth Group (e^(n*log(n)))^ZZ * (e^n)^ZZ * n^ZZ * log(n)^ZZ,
         (e^(n*log(n)), e^n, n, log(n)))

    ::

        sage: GrowthGroup('m^(ZZ[I]) * log(m)^(ZZ[I]) * n^(ZZ[I])')
        Growth Group m^ZZ * log(m)^ZZ * m^(ZZ*I) * log(m)^(ZZ*I) * n^ZZ * n^(ZZ*I)

    ::

        sage: TestSuite(GrowthGroup('x^ZZ')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_inverse() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(GrowthGroup('(QQ_+)^y')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_inverse() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(GrowthGroup('x^QQ * log(x)^ZZ')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_inverse() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """

    def create_key_and_extra_args(self, specification, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup.create_key_and_extra_args('asdf')
            Traceback (most recent call last):
            ...
            ValueError: 'asdf' is not a valid substring of 'asdf' describing a growth group.
            sage: GrowthGroup('as^df')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 'as^df' is not a valid substring of as^df
            describing a growth group.
            > *previous* ValueError: Cannot create a parent out of 'as'.
            >> *previous* ValueError: unknown specification as
            >> *and* SyntaxError: ... (<string>, line 1)
            > *and* ValueError: Cannot create a parent out of 'df'.
            >> *previous* ValueError: unknown specification df
            >> *and* NameError: name 'df' is not defined
            sage: GrowthGroup('x^y^z')
            Traceback (most recent call last):
            ...
            ValueError: 'x^y^z' is an ambiguous substring of
            a growth group description of 'x^y^z'.
            Use parentheses to make it unique.
            sage: GrowthGroup('(x^y)^z')
            Traceback (most recent call last):
            ...
            ValueError: '(x^y)^z' is not a valid substring of (x^y)^z
            describing a growth group.
            > *previous* ValueError: Cannot create a parent out of 'x^y'.
            >> *previous* ValueError: unknown specification x^y
            >> *and* NameError: name 'x' is not defined
            > *and* ValueError: Cannot create a parent out of 'z'.
            >> *previous* ValueError: unknown specification z
            >> *and* NameError: name 'z' is not defined
            sage: GrowthGroup('x^(y^z)')
            Traceback (most recent call last):
            ...
            ValueError: 'x^(y^z)' is not a valid substring of x^(y^z)
            describing a growth group.
            > *previous* ValueError: Cannot create a parent out of 'x'.
            >> *previous* ValueError: unknown specification x
            >> *and* NameError: name 'x' is not defined
            > *and* ValueError: Cannot create a parent out of 'y^z'.
            >> *previous* ValueError: unknown specification y^z
            >> *and* NameError: name 'y' is not defined

        ::

            sage: GrowthGroup('n^(I*ZZ)')
            Growth Group n^(ZZ*I)
            sage: GrowthGroup('n^(I  *   ZZ)')
            Growth Group n^(ZZ*I)
        """
        from .misc import repr_short_to_parent, split_str_by_op
        from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup

        kwds.setdefault('ignore_variables', ('e',))

        sfactors = split_str_by_op(
            ' '.join(specification.split()).replace('**', '^'), '*')

        def remove_parentheses(s):
            while s.startswith('(') and s.endswith(')'):
                s = s[1:-1].strip()
            return s

        def has_l_property(s, properties, invert=False):
            for p in properties:
                if s.startswith(p):
                    return s[len(p):].strip(), not invert
            return s, invert

        def has_r_property(s, properties, invert=False):
            for p in properties:
                if s.endswith(p):
                    return s[:-len(p)].strip(), not invert
            return s, invert

        factors = []

        for factor in sfactors:
            if '^' not in factor:
                raise ValueError("'{}' is not a valid substring of '{}' describing "
                                 "a growth group.".format(factor, specification))

            split = split_str_by_op(factor, '^')
            if len(split) != 2:
                raise ValueError("'{}' is an ambiguous substring of a growth group "
                                 "description of '{}'. Use parentheses to make it "
                                 "unique.".format(factor, ' * '.join(sfactors)))

            b, e = split
            b = remove_parentheses(b)
            e = remove_parentheses(e)

            b, extend_B_by_non_growth_group = has_r_property(
                b, ['_+'], invert=True)
            e, extend_E_by_non_growth_group = has_r_property(
                e, ['[I]', '[i]'], invert=False)
            e, l_E_only_imaginary_group = has_l_property(e, ['I*', 'I *'])
            e, r_E_only_imaginary_group = has_r_property(e, ['*I', '* I'])
            E_only_imaginary_group = l_E_only_imaginary_group or r_E_only_imaginary_group
            if E_only_imaginary_group and extend_E_by_non_growth_group:
                raise ValueError("'{}' is not a valid substring of '{}' describing "
                                 "a growth group.".format(factor, specification))

            try:
                B = repr_short_to_parent(b)
            except ValueError as exc:
                exc_b = exc
                exc_b.__traceback__ = None
                B = None
            try:
                E = repr_short_to_parent(e)
            except ValueError as exc:
                exc_e = exc
                exc_e.__traceback__ = None
                E = None

            if B is None and E is None:
                from .misc import combine_exceptions
                raise combine_exceptions(
                    ValueError("'{}' is not a valid substring of {} describing "
                               "a growth group.".format(factor, ' * '.join(sfactors))),
                    exc_b, exc_e)
            elif B is None and E is not None:
                if E_only_imaginary_group:
                    E = ImaginaryGroup(E)
                factors.append(GrowthGroupFactor(
                    cls=MonomialGrowthGroup,
                    base=E,
                    var=b,
                    extend_by_non_growth_group=extend_E_by_non_growth_group))
            elif B is not None and E is None:
                factors.append(GrowthGroupFactor(
                    cls=ExponentialGrowthGroup,
                    base=B,
                    var=e,
                    extend_by_non_growth_group=extend_B_by_non_growth_group))
            else:
                raise ValueError("'{}' is an ambiguous substring of a growth group "
                                 "description of '{}'.".format(factor, ' * '.join(factors)))

        return tuple(factors), kwds

    def create_object(self, version, factors, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^n')  # indirect doctest
            Growth Group QQ^n * Signs^n
        """
        groups = []
        non_growth_groups = []
        for factor in factors:
            grps = factor.cls.factory(
                factor.base,
                factor.var,
                extend_by_non_growth_group=factor.extend_by_non_growth_group,
                return_factors=True,
                **kwds)
            for grp in grps:
                if isinstance(grp, GenericNonGrowthGroup):
                    non_growth_groups.append(grp)
                else:
                    groups.append(grp)
        groups.extend(non_growth_groups)

        if len(groups) == 1:
            return groups[0]

        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product(groups)


GrowthGroup = GrowthGroupFactory("sage.rings.asymptotic.growth_group.GrowthGroup")
r"""
A factory for growth groups.
This is an instance of :class:`GrowthGroupFactory` whose documentation
provides more details.
"""
