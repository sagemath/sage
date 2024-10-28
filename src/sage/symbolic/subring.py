r"""
Subrings of the Symbolic Ring

Subrings of the symbolic ring can be created via the
:meth:`~sage.symbolic.ring.SymbolicRing.subring` method of
``SR``. This will call :class:`SymbolicSubring <SymbolicSubringFactory>`
of this module.

The following kinds of subrings are supported:

- A symbolic subring of expressions, whose variables are contained in
  a given set of symbolic variables (see
  :class:`SymbolicSubringAcceptingVars`). E.g.
  ::

      sage: SR.subring(accepting_variables=('a', 'b'))
      Symbolic Subring accepting the variables a, b

- A symbolic subring of expressions, whose variables are disjoint to a
  given set of symbolic variables (see
  :class:`SymbolicSubringRejectingVars`). E.g.
  ::

      sage: SR.subring(rejecting_variables=('r', 's'))
      Symbolic Subring rejecting the variables r, s

- The subring of symbolic constants (see
  :class:`SymbolicConstantsSubring`). E.g.
  ::

      sage: SR.subring(no_variables=True)
      Symbolic Constants Subring


TESTS:

In the following we have a couple of tests to see whether the coercion
framework works properly::

    sage: from sage.symbolic.subring import SymbolicSubring
    sage: V = var('a, r, x')
    sage: A = SymbolicSubring(accepting_variables=(a,)); A
    Symbolic Subring accepting the variable a
    sage: R = SymbolicSubring(rejecting_variables=(r,)); R
    Symbolic Subring rejecting the variable r
    sage: C = SymbolicSubring(no_variables=True); C
    Symbolic Constants Subring

::

    sage: sage.categories.pushout.pushout(A, R)
    Symbolic Subring rejecting the variable r
    sage: sage.categories.pushout.pushout(R, C)
    Symbolic Subring rejecting the variable r
    sage: sage.categories.pushout.pushout(C, A)
    Symbolic Subring accepting the variable a
    sage: sage.categories.pushout.pushout(A, SR)
    Symbolic Ring
    sage: sage.categories.pushout.pushout(R, SR)
    Symbolic Ring
    sage: sage.categories.pushout.pushout(C, SR)
    Symbolic Ring

::

    sage: cm = sage.structure.element.get_coercion_model()
    sage: cm.common_parent(A, R)
    Symbolic Subring rejecting the variable r
    sage: cm.common_parent(R, C)
    Symbolic Subring rejecting the variable r
    sage: cm.common_parent(C, A)
    Symbolic Subring accepting the variable a
    sage: cm.common_parent(A, SR)
    Symbolic Ring
    sage: cm.common_parent(R, SR)
    Symbolic Ring
    sage: cm.common_parent(C, SR)
    Symbolic Ring


AUTHORS:

- Daniel Krenn (2015)


Classes and Methods
===================
"""

# ****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc
from .ring import SymbolicRing, SR
from sage.categories.pushout import ConstructionFunctor
from sage.structure.factory import UniqueFactory


class SymbolicSubringFactory(UniqueFactory):
    r"""
    A factory creating a symbolic subring.

    INPUT:

    Specify one of the following keywords to create a subring.

    - ``accepting_variables`` -- (default: ``None``) a tuple or other
      iterable of variables. If specified, then a symbolic subring of
      expressions in only these variables is created.

    - ``rejecting_variables`` -- (default: ``None``) a tuple or other
      iterable of variables. If specified, then a symbolic subring of
      expressions in variables distinct to these variables is
      created.

    - ``no_variables`` -- boolean (default: ``False``); if set,
      then a symbolic subring of constant expressions (i.e.,
      expressions without a variable) is created.

    EXAMPLES::

        sage: from sage.symbolic.subring import SymbolicSubring
        sage: V = var('a, b, c, r, s, t, x, y, z')

    ::

        sage: A = SymbolicSubring(accepting_variables=(a, b, c)); A
        Symbolic Subring accepting the variables a, b, c
        sage: tuple((v, v in A) for v in V)
        ((a, True), (b, True), (c, True),
         (r, False), (s, False), (t, False),
         (x, False), (y, False), (z, False))

    ::

        sage: R = SymbolicSubring(rejecting_variables=(r, s, t)); R
        Symbolic Subring rejecting the variables r, s, t
        sage: tuple((v, v in R) for v in V)
        ((a, True), (b, True), (c, True),
         (r, False), (s, False), (t, False),
         (x, True), (y, True), (z, True))

    ::

        sage: C = SymbolicSubring(no_variables=True); C
        Symbolic Constants Subring
        sage: tuple((v, v in C) for v in V)
        ((a, False), (b, False), (c, False),
         (r, False), (s, False), (t, False),
         (x, False), (y, False), (z, False))

    TESTS::

        sage: SymbolicSubring(accepting_variables=tuple()) is C
        True

    ::

        sage: SymbolicSubring(rejecting_variables=tuple()) is SR
        True
    """
    def create_key_and_extra_args(
            self, accepting_variables=None, rejecting_variables=None,
            no_variables=False, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        See :class:`SymbolicSubringFactory` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring.create_key_and_extra_args()
            Traceback (most recent call last):
            ...
            ValueError: cannot create a symbolic subring since nothing is specified
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     accepting_variables=('a',), rejecting_variables=('r',))
            Traceback (most recent call last):
            ...
            ValueError: cannot create a symbolic subring since input is ambiguous
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     accepting_variables=('a',), no_variables=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot create a symbolic subring since input is ambiguous
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     rejecting_variables=('r',), no_variables=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot create a symbolic subring since input is ambiguous
        """
        if accepting_variables is None and \
           rejecting_variables is None and \
           not no_variables:
            raise ValueError('cannot create a symbolic subring '
                             'since nothing is specified')
        if accepting_variables is not None and rejecting_variables is not None or \
           rejecting_variables is not None and no_variables or \
           no_variables and accepting_variables is not None:
            raise ValueError('cannot create a symbolic subring '
                             'since input is ambiguous')

        if accepting_variables is not None:
            vars = tuple(accepting_variables)
            if vars:
                cls = SymbolicSubringAcceptingVars
            else:
                cls = SymbolicConstantsSubring
        elif rejecting_variables is not None:
            vars = tuple(rejecting_variables)
            cls = SymbolicSubringRejectingVars
        elif no_variables:
            vars = tuple()
            cls = SymbolicConstantsSubring

        vars = tuple(sorted(iter(SR(v) for v in vars), key=str))
        return (cls, vars), kwds

    def create_object(self, version, key, **kwds):
        r"""
        Create an object from the given arguments.

        See :class:`SymbolicSubringFactory` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=tuple()) is SR  # indirect doctest
            True
        """
        cls, vars = key
        if cls is SymbolicSubringRejectingVars and not vars:
            return SR
        return cls(vars, **kwds)


SymbolicSubring = SymbolicSubringFactory("SymbolicSubring")


class GenericSymbolicSubring(SymbolicRing):

    def __init__(self, vars):
        r"""
        An abstract base class for a symbolic subring.

        INPUT:

        - ``vars`` -- tuple of symbolic variables

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
            sage: SymbolicSubring(no_variables=True)  # indirect doctest
            Symbolic Constants Subring
            sage: SymbolicSubring(rejecting_variables=tuple())  # indirect doctest
            Symbolic Ring

        ::

            sage: SR.subring(accepting_variables=(0, pi, sqrt(2), 'zzz', I))
            Traceback (most recent call last):
            ...
            ValueError: Invalid variables: 0, I, pi, sqrt(2)
        """
        super().__init__()
        self._vars_ = set(vars)
        if not all(v.is_symbol() for v in self._vars_):
            raise ValueError('Invalid variables: {}'.format(
                ', '.join(str(v) for v in sorted(self._vars_, key=str)
                          if not v.is_symbol())))

    def _repr_variables_(self):
        r"""
        Return a representation string of the variables.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=tuple())._repr_variables_()
            'no variable'
            sage: SymbolicSubring(accepting_variables=('a',))._repr_variables_()
            'the variable a'
            sage: SymbolicSubring(accepting_variables=('a', 'b'))._repr_variables_()
            'the variables a, b'
        """
        if not self._vars_:
            s = 'no variable'
        elif len(self._vars_) == 1:
            s = 'the variable '
        else:
            s = 'the variables '
        return s + ', '.join(str(v) for v in sorted(self._vars_, key=str))

    def has_valid_variable(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.subring import GenericSymbolicSubring
            sage: GenericSymbolicSubring(vars=tuple()).has_valid_variable(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not implemented in this abstract base class
        """
        raise NotImplementedError('Not implemented in this abstract base class')

    def _element_constructor_(self, x):
        r"""
        Create the element of this subring specified by the input ``x``.

        INPUT:

        - ``x`` -- an object

        OUTPUT: an element of this symbolic subring

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(accepting_variables=('a',))
            sage: S('a')  # indirect doctest
            a
            sage: _.parent()
            Symbolic Subring accepting the variable a
            sage: S('x')  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: x is not contained in Symbolic Subring accepting the variable a
        """
        expression = super()._element_constructor_(x)
        assert expression.parent() is self
        if not all(self.has_valid_variable(var)
                   for var in expression.variables()):
            raise TypeError('%s is not contained in %s' % (x, self))
        return expression

    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent

        OUTPUT: boolean or ``None``

        TESTS::

            sage: from sage.symbolic.subring import GenericSymbolicSubring
            sage: GenericSymbolicSubring(vars=tuple()).has_coerce_map_from(SR)  # indirect doctest  # not tested see #19231
            False

        ::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: C = SymbolicSubring(no_variables=True)
            sage: C.has_coerce_map_from(ZZ)  # indirect doctest
            True
            sage: C.has_coerce_map_from(QQ)  # indirect doctest
            True
            sage: C.has_coerce_map_from(RR)  # indirect doctest
            True
            sage: C.has_coerce_map_from(RIF)  # indirect doctest
            True
            sage: C.has_coerce_map_from(CC)  # indirect doctest
            True
            sage: C.has_coerce_map_from(CIF)  # indirect doctest
            True
            sage: C.has_coerce_map_from(AA)  # indirect doctest
            True
            sage: C.has_coerce_map_from(QQbar)  # indirect doctest
            True
            sage: C.has_coerce_map_from(SR)  # indirect doctest
            False
        """
        from sage.rings.infinity import InfinityRing
        from sage.rings.qqbar import AA, QQbar
        from sage.rings.real_lazy import RLF, CLF

        if isinstance(P, type):
            return SR._coerce_map_from_(P)

        if RLF.has_coerce_map_from(P) or \
           CLF.has_coerce_map_from(P) or \
           AA.has_coerce_map_from(P) or \
           QQbar.has_coerce_map_from(P):
            return True

        if (P is InfinityRing or
            isinstance(P, (sage.rings.abc.RealIntervalField,
                           sage.rings.abc.ComplexIntervalField))):
            return True

        if P._is_numerical():
            return P not in (RLF, CLF, AA, QQbar)

    def __eq__(self, other):
        """
        Compare two symbolic subrings.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: B = SymbolicSubring(accepting_variables=('b',))
            sage: AB = SymbolicSubring(accepting_variables=('a', 'b'))
            sage: A == A
            True
            sage: A == B
            False
            sage: A == AB
            False
        """
        if not isinstance(other, GenericSymbolicSubring):
            return False
        return self._vars_ == other._vars_

    def __ne__(self, other):
        """
        Check whether ``self`` and ``other`` are not equal.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: B = SymbolicSubring(accepting_variables=('b',))
            sage: AB = SymbolicSubring(accepting_variables=('a', 'b'))
            sage: A != A
            False
            sage: A != B
            True
            sage: A != AB
            True
        """
        return not self == other

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: B = SymbolicSubring(accepting_variables=('b',))
            sage: hash(A) == hash(A)
            True
            sage: hash(A) == hash(B)
            False
        """
        return hash(tuple(sorted(self._vars_)))


class GenericSymbolicSubringFunctor(ConstructionFunctor):
    r"""
    A base class for the functors constructing symbolic subrings.

    INPUT:

    - ``vars`` -- tuple, set, or other iterable of symbolic variables

    EXAMPLES::

        sage: from sage.symbolic.subring import SymbolicSubring
        sage: SymbolicSubring(no_variables=True).construction()[0]  # indirect doctest
        Subring<accepting no variable>

    .. SEEALSO::

        :class:`sage.categories.pushout.ConstructionFunctor`.
    """

    _functor_name = 'GenericSymbolicSubringFunctor'

    rank = 11

    # The symbolic subring construction returns an object admitting a
    # coercion map into the original, not vice versa.
    coercion_reversed = True

    _repr_type_ = 'generic'

    def __init__(self, vars):
        r"""
        See :class:`GenericSymbolicSubringFunctor` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).construction()[0]  # indirect doctest
            Subring<accepting a>
        """
        self.vars = set(vars)
        from sage.categories.rings import Rings
        super().__init__(Rings(), Rings())

    def _repr_variables_(self):
        r"""
        Return a representation string of the variables.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F._repr_variables_()
            'a'
        """
        return ', '.join(str(v) for v in sorted(self.vars, key=str))

    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
            sage: SymbolicSubring(no_variables=True)  # indirect doctest
            Symbolic Constants Subring
        """
        return 'Subring<%s%s%s>' % (
            self._repr_type_, ' ' if self._repr_type_ else '',
            self._repr_variables_() if self.vars else 'no variable')

    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: a functor or ``None``

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F.merge(F) is F
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

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F == F
            True
        """
        return type(self) is type(other) and self.vars == other.vars

    def __ne__(self, other):
        r"""
        Return whether this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F != F
            False
        """
        return not self == other


class SymbolicSubringAcceptingVars(GenericSymbolicSubring):
    r"""
    The symbolic subring consisting of symbolic expressions in the given variables.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
        """
        return 'Symbolic Subring accepting %s' % \
            (self._repr_variables_())

    def has_valid_variable(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(accepting_variables=('a',))
            sage: S.has_valid_variable('a')
            True
            sage: S.has_valid_variable('r')
            False
            sage: S.has_valid_variable('x')
            False
        """
        return SR(variable) in self._vars_

    def construction(self):
        r"""
        Return the functorial construction of this symbolic subring.

        OUTPUT:

        A tuple whose first entry is a construction functor and its second
        is the symbolic ring.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).construction()
            (Subring<accepting a>, Symbolic Ring)
        """
        return (SymbolicSubringAcceptingVarsFunctor(self._vars_), SR)

    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent

        OUTPUT: boolean or ``None``

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: AB = SymbolicSubring(accepting_variables=('a', 'b'))
            sage: A.has_coerce_map_from(AB)  # indirect doctest
            False
            sage: AB.has_coerce_map_from(A)  # indirect doctest
            True
        """
        if isinstance(P, SymbolicSubringAcceptingVars):
            return self._vars_ >= P._vars_
        return super()._coerce_map_from_(P)

    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT: a symbolic expression

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).an_element()
            a
            sage: _.parent()
            Symbolic Subring accepting the variable a
        """
        return self(sorted(self._vars_, key=str)[0])


class SymbolicSubringAcceptingVarsFunctor(GenericSymbolicSubringFunctor):

    _functor_name = 'SymbolicSubringAcceptingVarsFunctor'

    _repr_type_ = 'accepting'

    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: a functor or ``None``

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: G = SymbolicSubring(rejecting_variables=('r',)).construction()[0]
            sage: F.merge(F) is F
            True
            sage: F.merge(G) is G
            True
        """
        if self == other:
            return self
        elif type(self) is type(other):
            return type(self)(self.vars | other.vars)
        elif isinstance(other, SymbolicSubringRejectingVarsFunctor):
            if not (self.vars & other.vars):
                return other

    def _apply_functor(self, R):
        """
        Apply this functor to the given symbolic ring `R`.

        INPUT:

        - ``R`` -- a symbolic ring

        OUTPUT: a subring of ``R``

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F, R = SymbolicSubring(accepting_variables=('a',)).construction()
            sage: F(R)  # indirect doctest
            Symbolic Subring accepting the variable a

        TESTS::

            sage: F(F(R))
            Traceback (most recent call last):
            ...
            NotImplementedError: This functor can only be applied on the
            symbolic ring but Symbolic Subring accepting the variable a given.
        """
        if R is not SR:
            raise NotImplementedError('This functor can only be applied on '
                                      'the symbolic ring but %s given.' % (R,))
        return SymbolicSubring(accepting_variables=self.vars)


class SymbolicSubringRejectingVars(GenericSymbolicSubring):
    r"""
    The symbolic subring consisting of symbolic expressions whose variables
    are none of the given variables.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
        """
        return 'Symbolic Subring rejecting %s' % \
            (self._repr_variables_())

    def has_valid_variable(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(rejecting_variables=('r',))
            sage: S.has_valid_variable('a')
            True
            sage: S.has_valid_variable('r')
            False
            sage: S.has_valid_variable('x')
            True
        """
        return SR(variable) not in self._vars_

    def construction(self):
        r"""
        Return the functorial construction of this symbolic subring.

        OUTPUT:

        A tuple whose first entry is a construction functor and its second
        is the symbolic ring.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',)).construction()
            (Subring<rejecting r>, Symbolic Ring)
        """
        return (SymbolicSubringRejectingVarsFunctor(self._vars_), SR)

    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent

        OUTPUT: boolean or ``None``

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: R = SymbolicSubring(rejecting_variables=('r',))
            sage: RS = SymbolicSubring(rejecting_variables=('r', 's'))
            sage: RS.has_coerce_map_from(R)  # indirect doctest
            False
            sage: R.has_coerce_map_from(RS)  # indirect doctest
            True
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: R.has_coerce_map_from(A)
            True
            sage: A.has_coerce_map_from(R)
            False
        """
        if isinstance(P, SymbolicSubringRejectingVars):
            return self._vars_ <= P._vars_
        elif isinstance(P, SymbolicSubringAcceptingVars):
            return not (self._vars_ & P._vars_)
        return super()._coerce_map_from_(P)

    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT: a symbolic expression

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',)).an_element()
            some_variable
            sage: _.parent()
            Symbolic Subring rejecting the variable r
            sage: SymbolicSubring(rejecting_variables=('some_variable',)).an_element()
            some_some_variable
            sage: _.parent()
            Symbolic Subring rejecting the variable some_variable
            sage: SymbolicSubring(rejecting_variables=('some_some_variable',)).an_element()
            some_variable
            sage: _.parent()
            Symbolic Subring rejecting the variable some_some_variable
            sage: SymbolicSubring(rejecting_variables=('some_variable','some_some_variable')).an_element()
            some_some_some_variable
            sage: _.parent()
            Symbolic Subring rejecting the variables some_some_variable, some_variable
        """
        v = SR.an_element()
        while not self.has_valid_variable(v):
            v = SR('some_' + str(v))
        return self(v)


class SymbolicSubringRejectingVarsFunctor(GenericSymbolicSubringFunctor):

    _functor_name = 'SymbolicSubringRejectingVarsFunctor'

    _repr_type_ = 'rejecting'

    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor

        OUTPUT: a functor or ``None``

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: G = SymbolicSubring(rejecting_variables=('r',)).construction()[0]
            sage: G.merge(G) is G
            True
            sage: G.merge(F) is G
            True
        """
        if self == other:
            return self
        elif type(self) is type(other):
            return type(self)(self.vars & other.vars)
        elif isinstance(other, SymbolicSubringAcceptingVarsFunctor):
            if not (self.vars & other.vars):
                return self

    def _apply_functor(self, R):
        """
        Apply this functor to the given symbolic ring `R`.

        INPUT:

        - ``R`` -- a symbolic ring

        OUTPUT: a subring of ``R``

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F, R = SymbolicSubring(rejecting_variables=('r',)).construction()
            sage: F(R)  # indirect doctest
            Symbolic Subring rejecting the variable r

        TESTS::

            sage: F(F(R))
            Traceback (most recent call last):
            ...
            NotImplementedError: This functor can only be applied on the
            symbolic ring but Symbolic Subring rejecting the variable r given.
        """
        if R is not SR:
            raise NotImplementedError('This functor can only be applied on '
                                      'the symbolic ring but %s given.' % (R,))
        return SymbolicSubring(rejecting_variables=self.vars)


class SymbolicConstantsSubring(SymbolicSubringAcceptingVars):
    r"""
    The symbolic subring consisting of symbolic constants.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT: string

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(no_variables=True)  # indirect doctest
            Symbolic Constants Subring
        """
        return 'Symbolic Constants Subring'

    def has_valid_variable(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(no_variables=True)
            sage: S.has_valid_variable('a')
            False
            sage: S.has_valid_variable('r')
            False
            sage: S.has_valid_variable('x')
            False
        """
        return False

    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT: a symbolic expression

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(no_variables=True).an_element()
            I*pi*e
            sage: _.parent()
            Symbolic Constants Subring
        """
        return self(SR('I') * SR('pi') * SR('e'))
