# sage_setup: distribution = sagemath-objects
"""
Coerce actions
"""

#*****************************************************************************
#     Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

from cpython.long cimport *
from cpython.number cimport *
from cysignals.signals cimport sig_check

from sage.structure.coerce cimport coercion_model
from sage.structure.element cimport parent, Element, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.coerce_exceptions import CoercionException
from sage.categories.action cimport InverseAction, PrecomposedAction
from sage.arith.long cimport integer_check_long


cdef _record_exception():
    coercion_model._record_exception()

cdef inline an_element(R):
    if isinstance(R, Parent):
        return R.an_element()
    else:
        for x in ([(1, 2)], "abc", 10.5, 10):
            try:
                return R(x)
            except Exception:
                pass


# In the code below, I take the convention that g is acting on a.

cdef class GenericAction(Action):

    cdef _codomain

    def __init__(self, Parent G, S, is_left, bint check=True):
        """
        TESTS:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import ActedUponAction, GenericAction
            sage: M = MatrixSpace(ZZ, 2)                                                # needs sage.modules
            sage: ActedUponAction(M, Cusps, True)                                       # needs sage.modular sage.modules
            Left action
             by Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
             on Set P^1(QQ) of all cusps

            sage: Z6 = Zmod(6)
            sage: GenericAction(QQ, Z6, True)
            Traceback (most recent call last):
            ...
            NotImplementedError: action for <class 'sage.structure.coerce_actions.GenericAction'> not implemented

        This will break if we tried to use it::

            sage: GenericAction(QQ, Z6, True, check=False)
            Left action by Rational Field on Ring of integers modulo 6
        """
        Action.__init__(self, G, S, is_left, operator.mul)
        if check:
            res = self.act(G.an_element(), S.an_element())
            if res is None:
                raise CoercionException
            _codomain = parent(res)

    def codomain(self):
        """
        Return the "codomain" of this action, i.e. the Parent in which the
        result elements live. Typically, this should be the same as the
        acted upon set.

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains, for
        otherwise they could be garbage collected, giving rise to random
        errors (see :issue:`18157`). ::


            sage: M = MatrixSpace(ZZ, 2)                                                # needs sage.modules
            sage: A = sage.structure.coerce_actions.ActedUponAction(M, Cusps, True)     # needs sage.modular sage.modules
            sage: A.codomain()                                                          # needs sage.modular sage.modules
            Set P^1(QQ) of all cusps

            sage: # needs sage.groups
            sage: S3 = SymmetricGroup(3)
            sage: QQxyz = QQ['x,y,z']
            sage: A = sage.structure.coerce_actions.ActOnAction(S3, QQxyz, False)
            sage: A.codomain()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if self._codomain is None:
            self._codomain = parent(self.act(an_element(self.G),
                                               an_element(self.underlying_set())))
        return self._codomain


cdef class ActOnAction(GenericAction):
    """
    Class for actions defined via the _act_on_ method.
    """
    cpdef _act_(self, g, x):
        """
        TESTS::

            sage: # needs sage.groups
            sage: G = SymmetricGroup(3)
            sage: R.<x,y,z> = QQ[]
            sage: A = sage.structure.coerce_actions.ActOnAction(G, R, False)
            sage: A(x^2 + y - z, G((1,2)))
            y^2 + x - z
            sage: A(x + 2*y + 3*z, G((1,3,2)))
            2*x + 3*y + z
            sage: type(A)
            <... 'sage.structure.coerce_actions.ActOnAction'>
        """
        return (<Element>g)._act_on_(x, self._is_left)


cdef class ActedUponAction(GenericAction):
    """
    Class for actions defined via the _acted_upon_ method.
    """
    cpdef _act_(self, g, x):
        """
        TESTS::

            sage: # needs sage.modular sage.modules
            sage: M = MatrixSpace(ZZ, 2)
            sage: A = sage.structure.coerce_actions.ActedUponAction(M, Cusps, True)
            sage: A.act(matrix(ZZ, 2, [1,0,2,-1]), Cusp(1,2))
            Infinity
            sage: A(matrix(ZZ, 2, [1,0,2,-1]), Cusp(1,2))
            Infinity
            sage: type(A)
            <... 'sage.structure.coerce_actions.ActedUponAction'>
        """
        return (<Element>x)._acted_upon_(g, not self._is_left)


def detect_element_action(Parent X, Y, bint X_on_left, X_el=None, Y_el=None):
    r"""
    Return an action of X on Y as defined by elements of X, if any.

    EXAMPLES:

    Note that coerce actions should only be used inside of the coercion
    model. For this test, we need to strongly reference the domains,
    for otherwise they could be garbage collected, giving rise to
    random errors (see :issue:`18157`). ::

        sage: from sage.structure.coerce_actions import detect_element_action
        sage: ZZx = ZZ['x']
        sage: M = MatrixSpace(ZZ, 2)                                                    # needs sage.modules
        sage: detect_element_action(ZZx, ZZ, False)
        Left scalar multiplication by Integer Ring
         on Univariate Polynomial Ring in x over Integer Ring
        sage: detect_element_action(ZZx, QQ, True)
        Right scalar multiplication by Rational Field
         on Univariate Polynomial Ring in x over Integer Ring
        sage: detect_element_action(Cusps, M, False)                                    # needs sage.modular sage.modules
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
         on Set P^1(QQ) of all cusps
        sage: detect_element_action(Cusps, M, True),                                    # needs sage.modular sage.modules
        (None,)
        sage: detect_element_action(ZZ, QQ, True),
        (None,)

    TESTS:

    This test checks that the issue in :issue:`7718` has been fixed::

        sage: class MyParent(Parent):
        ....:     def an_element(self):
        ....:         pass
        ....:
        sage: A = MyParent()
        sage: detect_element_action(A, ZZ, True)
        Traceback (most recent call last):
        ...
        RuntimeError: an_element() for <__main__.MyParent object at ...> returned None
    """
    cdef Element x

    # sample elements x and y
    if X_el is None or (parent(X_el) is not X):
        x = an_element(X)
    else:
        x = X_el
    if x is None:
        raise RuntimeError("an_element() for %s returned None" % X)
    if Y_el is None or (parent(Y_el) is not Y):
        y = an_element(Y)
    else:
        y = Y_el
    if y is None:
        if isinstance(Y, Parent):
            raise RuntimeError("an_element() for %s returned None" % Y)
        else:
            return # don't know how to make elements of this type...

    # element x defining _lmul_ or _rmul_
    if isinstance(x, ModuleElement) and isinstance(y, Element):
        try:
            return (RightModuleAction if X_on_left else LeftModuleAction)(Y, X, y, x)
        except CoercionException as msg:
            _record_exception()

    # element x defining _act_on_
    try:
        if x._act_on_(y, X_on_left) is not None:
            return ActOnAction(X, Y, X_on_left, False)
    except CoercionException:
        _record_exception()

    # element x defining _acted_upon_
    if isinstance(Y, Parent):
        try:
            if x._acted_upon_(y, X_on_left) is not None:
                return ActedUponAction(Y, X, not X_on_left, False)
        except CoercionException:
            _record_exception()


cdef class ModuleAction(Action):
    """
    Module action.

    .. SEEALSO::

        This is an abstract class, one must actually instantiate a
        :class:`LeftModuleAction` or a :class:`RightModuleAction`.

    INPUT:

    - ``G`` -- the actor, an instance of :class:`~sage.structure.parent.Parent`
    - ``S`` -- the object that is acted upon
    - ``g`` -- (optional) an element of ``G``
    - ``a`` -- (optional) an element of ``S``
    - ``check`` -- if ``True`` (default), then there will be no consistency tests
      performed on sample elements

    NOTE:

    By default, the sample elements of ``S`` and ``G`` are obtained from
    :meth:`~sage.structure.parent.Parent.an_element`, which relies on the
    implementation of an ``_an_element_()`` method. This is not always
    available. But usually, the action is only needed when one already
    *has* two elements. Hence, by :issue:`14249`, the coercion model will
    pass these two elements to the :class:`ModuleAction` constructor.

    The actual action is implemented by the ``_rmul_`` or ``_lmul_``
    function on its elements. We must, however, be very particular about
    what we feed into these functions, because they operate under the
    assumption that the inputs lie exactly in the base ring and may
    segfault otherwise. Thus we handle all possible base extensions
    manually here.
    """
    def __init__(self, G, S, g=None, a=None, check=True):
        """
        This creates an action of an element of a module by an element of its
        base ring. The simplest example to keep in mind is R acting on the
        polynomial ring R[x].

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::


            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: ZZx = ZZ['x']
            sage: QQx = QQ['x']
            sage: QQy = QQ['y']
            sage: ZZxy = ZZ['x']['y']
            sage: LeftModuleAction(ZZ, ZZx)
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: LeftModuleAction(ZZ, QQx)
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Rational Field
            sage: LeftModuleAction(QQ, ZZx)
            Left scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: LeftModuleAction(QQ, ZZxy)
            Left scalar multiplication by Rational Field
             on Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring

        The following tests against a problem that was relevant during work on
        :issue:`9944`::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<x> = PolynomialRing(ZZ, sparse=True)
            sage: 1/R.0
            1/x
            sage: 1/S.0
            1/x

        If there is a coercion from ``G`` to ``S``, we do not create
        the module action of ``G`` on the pushout of ``G`` and ``S``::

            sage: # needs sage.modules
            sage: G = PolynomialRing(QQ, "x")
            sage: S = PolynomialRing(MatrixSpace(QQ, 2), "x")
            sage: G.gen() * S.gen()
            [1 0]
            [0 1]*x^2

        Contrast the previous example with the following, where we
        have no coercion from ``G`` to ``S``::

            sage: # needs sage.modules
            sage: S = PolynomialRing(MatrixSpace(QQ, 2), "y")
            sage: G.gen() * S.gen()
            [x 0]
            [0 x]*y
        """
        Action.__init__(self, G, S, not isinstance(self, RightModuleAction), operator.mul)
        if not isinstance(G, Parent):
            # only let Parents act
            raise TypeError("Actor must be a parent.")
        base = S.base()
        if base is S or base is None:
            # The right thing to do is a normal multiplication
            raise CoercionException("Best viewed as standard multiplication")
        # Objects are implemented with the assumption that
        # _lmul_/_rmul_ are given an element of the base ring
        if G is not base:
            # first we try the easy case of coercing G to the base ring of S
            self.connecting = base._internal_coerce_map_from(G)
            if self.connecting is None:
                if S._internal_coerce_map_from(G) is not None:
                    raise CoercionException("Best viewed as standard coercion multiplication.")
                # otherwise, we try and find a base extension
                from sage.categories.pushout import pushout
                # this may raise a type error, which we propagate
                self.extended_base = pushout(G, S)
                # make sure the pushout actually gave a correct base extension of S
                if self.extended_base.base() != pushout(G, base):
                    raise CoercionException("Actor must be coercible into base.")
                else:
                    self.connecting = self.extended_base.base()._internal_coerce_map_from(G)
                    if self.connecting is None:
                        # this may happen if G is, say, int rather than a parent
                        # TODO: let python types be valid actions
                        raise CoercionException("Missing connecting morphism")

        # Don't waste time if our connecting morphisms happen to be the identity.
        if self.connecting is not None and self.connecting.codomain() is G:
            self.connecting = None

        if self.extended_base is not None and self.extended_base is S:
            self.extended_base = None

        # At this point, we can assert it is safe to call _Xmul_
        the_ring = G if self.connecting is None else self.connecting.codomain()
        the_set = S if self.extended_base is None else self.extended_base
        assert the_ring is the_set.base(), "BUG in coercion model\n    Apparently there are two versions of\n        %s\n    in the cache." % the_ring

        if not check:
            return
        if g is None:
            g = G.an_element()
        if parent(g) is not G:
            raise CoercionException("The parent of %s is not %s but %s" % (g, G, parent(g)))
        if a is None:
            a = S.an_element()
        if parent(a) is not S:
            raise CoercionException("The parent of %s is not %s but %s" % (a, S, parent(a)))
        if not isinstance(g, Element) or not isinstance(a, ModuleElement):
            raise CoercionException("not an Element acting on a ModuleElement")
        res = self.act(g, a)
        if parent(res) is not the_set:
            # In particular we will raise an error if res is None
            raise CoercionException("Result is None or has wrong parent.")

    def __reduce__(self):
        """
        Used in pickling.

        TESTS:

        Check that this action can be pickled (:issue:`29031`)::

            sage: A = ZZ['x'].get_action(QQ, self_on_left=False, op=operator.mul)
            sage: loads(dumps(A)) is not None
            True
        """
        return (type(self), (self.G, self.underlying_set()))

    def _repr_name_(self):
        """
        The default name of this action type, which is has a sane default.

        EXAMPLES:

        Note that coerce actions should only be used inside of the
        coercion model. For this test, we need to strongly reference the
        domains, for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
            sage: ZZx = ZZ['x']
            sage: A = LeftModuleAction(ZZ, ZZx); A
            Left scalar multiplication by Integer Ring
             on Univariate Polynomial Ring in x over Integer Ring
            sage: A._repr_name_()
            'scalar multiplication'

            sage: GF5 = GF(5)
            sage: GF5t = GF5[['t']]
            sage: RightModuleAction(GF5, GF5t)
            Right scalar multiplication by Finite Field of size 5
             on Power Series Ring in t over Finite Field of size 5
        """
        return "scalar multiplication"

    def codomain(self):
        """
        The codomain of self, which may or may not be equal to the domain.

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: ZZxyz = ZZ['x,y,z']
            sage: A = LeftModuleAction(QQ, ZZxyz)
            sage: A.codomain()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if self.extended_base is not None:
            return self.extended_base
        return self.underlying_set()

    def domain(self):
        """
        The domain of self, which is the module that is being acted on.

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: ZZxyz = ZZ['x,y,z']
            sage: A = LeftModuleAction(QQ, ZZxyz)
            sage: A.domain()
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        return self.underlying_set()

    def __invert__(self):
        """
        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains, for
        otherwise they could be garbage collected, giving rise to random
        errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import RightModuleAction

            sage: ZZx = ZZ['x']
            sage: x = ZZx.gen()
            sage: QQx = QQ['x']
            sage: A = ~RightModuleAction(QQ, QQx); A
            Right inverse action by Rational Field on Univariate Polynomial Ring in x over Rational Field
            sage: A(x, 2)
            1/2*x

            sage: A = ~RightModuleAction(QQ, ZZx); A
            Right inverse action by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: A.codomain()
            Univariate Polynomial Ring in x over Rational Field
            sage: A(x, 2)
            1/2*x

            sage: A = ~RightModuleAction(ZZ, ZZx); A
            Right inverse action by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: A.codomain()
            Univariate Polynomial Ring in x over Rational Field
            sage: A(x, 2)
            1/2*x

            sage: GF5x = GF(5)['x']
            sage: A = ~RightModuleAction(ZZ, GF5x); A
            Right inverse action by Finite Field of size 5
            on Univariate Polynomial Ring in x over Finite Field of size 5
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 5
            sage: A(x, 2)
            3*x

            sage: GF5xy = GF5x['y']
            sage: A = ~RightModuleAction(ZZ, GF5xy); A
            Right inverse action by Finite Field of size 5
            on Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Finite Field of size 5
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 5

            sage: ZZy = ZZ['y']
            sage: ZZxyzw = ZZx['y']['z']['w']
            sage: A = ~RightModuleAction(ZZy, ZZxyzw); A
            Right inverse action by Fraction Field of Univariate Polynomial Ring in y
                over Univariate Polynomial Ring in x
                over Integer Ring
            on Univariate Polynomial Ring in w
                over Univariate Polynomial Ring in z
                over Univariate Polynomial Ring in y
                over Univariate Polynomial Ring in x
                over Integer Ring
            with precomposition on right by Conversion via FractionFieldElement map:
              From: Univariate Polynomial Ring in y over Integer Ring
              To:   Fraction Field of Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring

        TESTS:

        See :issue:`19521`::

            sage: # needs sage.symbolic
            sage: Q.<y> = SR.subring(no_variables=True)[[]]
            sage: (y / 1).parent()
            Power Series Ring in y over Symbolic Constants Subring
            sage: R.<x> = SR.subring(no_variables=True)[]
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(x, 1, operator.truediv)
            Action discovered.
                Right inverse action by Symbolic Constants Subring
                 on Univariate Polynomial Ring in x over Symbolic Constants Subring
                with precomposition on right by Conversion via _symbolic_ method map:
                  From: Integer Ring
                  To:   Symbolic Constants Subring
            Result lives in Univariate Polynomial Ring in x over Symbolic Constants Subring
            Univariate Polynomial Ring in x over Symbolic Constants Subring
        """
        K = self.G._pseudo_fraction_field()
        if K is self.G:
            return InverseAction(self)
        else:
            # Try to find a suitable ring between G and R in which to compute
            # the inverse.
            from sage.categories.pushout import construction_tower, pushout
            R = self.G if self.connecting is None else self.connecting.codomain()
            K = pushout(self.G._pseudo_fraction_field(), R)
            if K is None:
                K = R._pseudo_fraction_field()

            # Suppose we have parents of a construction tower
            #   A -> B -> C <- D <- E -> F <- G -> H,
            # where the arrows specify the direction of coercion (i.e.
            # A -> B means A coerces to B). Note that B = FA(A), C = FB(B), ...
            # As we want to find a "smallest" parent with some properties,
            # we need the order
            #   A, B, E, D, C, G, F, H
            # for our search. Thus the elements connected with a <- have to
            # be reversed. See code below.
            tower = []
            reversed_part = []
            for Fi, Ri in reversed(construction_tower(K)):
                if Fi is not None and Fi.coercion_reversed:
                    reversed_part.append((Fi, Ri))
                else:
                    tower.append((Fi, Ri))
                    if reversed_part:
                        tower += reversed(reversed_part)
                        reversed_part = []
            assert(not reversed_part)
            for _, Ri in tower:
                if not Ri.has_coerce_map_from(self.G):
                    continue
                Ki = Ri._pseudo_fraction_field()
                if Ki is Ri:
                    K = Ki
                    break
            module_action = type(self)(K, self.codomain())
            connecting = K.coerce_map_from(self.G)
            left, right = (connecting, None) if self._is_left else (None, connecting)
            return PrecomposedAction(InverseAction(module_action), left, right)


cdef class LeftModuleAction(ModuleAction):
    cpdef _act_(self, g, a):
        """
        A left module action is an action that takes the ring element as the
        first argument (the left side) and the module element as the second
        argument (the right side).

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: R.<x> = QQ['x']
            sage: A = LeftModuleAction(ZZ, R)
            sage: A(5, x+1)
            5*x + 5
            sage: R.<x> = ZZ['x']
            sage: A = LeftModuleAction(QQ, R)
            sage: A(1/2, x+1)
            1/2*x + 1/2
            sage: A(1/2, x+1)
            1/2*x + 1/2
        """
        # The way we are called, we know for sure that a is a
        # ModuleElement and that g is an Element.
        # If we change a or g, we need to explicitly check this,
        # which explains the <?> casts below.
        if self.connecting is not None:
            g = <Element?>self.connecting._call_(g)
        if self.extended_base is not None:
            a = <ModuleElement?>self.extended_base(a)
        return (<ModuleElement>a)._rmul_(<Element>g)  # g * a


cdef class RightModuleAction(ModuleAction):
    cpdef _act_(self, g, a):
        """
        A right module action is an action that takes the module element as the
        first argument (the left side) and the ring element as the second
        argument (the right side).

        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the domains,
        for otherwise they could be garbage collected, giving rise to
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import RightModuleAction
            sage: R.<x> = QQ['x']
            sage: A = RightModuleAction(ZZ, R)
            sage: A(x+5, 2)
            2*x + 10
            sage: A(x+5, 2)
            2*x + 10
        """
        # The way we are called, we know for sure that a is a
        # ModuleElement and that g is an Element.
        # If we change a or g, we need to explicitly check this,
        # which explains the <?> casts below.
        if self.connecting is not None:
            g = <Element?>self.connecting._call_(g)
        if self.extended_base is not None:
            a = <ModuleElement?>self.extended_base(a)
        return (<ModuleElement>a)._lmul_(<Element>g)  # a * g


cdef class IntegerAction(Action):
    """
    Abstract base class representing some action by integers on
    something. Here, "integer" is defined loosely in the "duck typing"
    sense.

    INPUT:

    - ``Z`` -- a type or parent representing integers

    For the other arguments, see :class:`Action`.

    .. NOTE::

        This class is used internally in Sage's coercion model. Outside
        of the coercion model, special precautions are needed to prevent
        domains of the action from being garbage collected.
    """
    def __init__(self, Z, S, is_left, op):
        if isinstance(Z, type):
            from sage.sets.pythonclass import Set_PythonType
            Z = Set_PythonType(Z)
        super().__init__(Z, S, is_left, op)

    def __reduce__(self):
        """
        Used in pickling.

        TESTS:

        Check that this action can be pickled (:issue:`29031`)::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: act = IntegerMulAction(ZZ, CDF)                                       # needs sage.rings.complex_double
            sage: loads(dumps(act)) is not None                                         # needs sage.rings.complex_double
            True
        """
        # All base classes must take the signature
        #   (Z, S, is_left)
        return (type(self), (self.G, self.underlying_set(), self._is_left))

    def __invert__(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: act = IntegerMulAction(ZZ, CDF)                                       # needs sage.rings.complex_double
            sage: ~act                                                                  # needs sage.rings.complex_double
            Traceback (most recent call last):
            ...
            TypeError: actions by ZZ cannot be inverted
        """
        raise TypeError("actions by ZZ cannot be inverted")


cdef class IntegerMulAction(IntegerAction):
    r"""
    Implement the action `n \cdot a = a + a + ... + a` via repeated
    doubling.

    Both addition and negation must be defined on the set `M`.

    INPUT:

    - ``Z`` -- a type or parent representing integers

    - ``M`` -- a ``ZZ``-module

    - ``m`` -- (optional) an element of ``M``

    EXAMPLES::

        sage: from sage.structure.coerce_actions import IntegerMulAction
        sage: R.<x> = QQ['x']
        sage: act = IntegerMulAction(ZZ, R)
        sage: act(5, x)
        5*x
        sage: act(0, x)
        0
        sage: act(-3, x-1)
        -3*x + 3
    """
    def __init__(self, Z, M, is_left=True, m=None):
        if m is None:
            m = M.an_element()
        test = m + (-m)  # make sure addition and negation is allowed
        super().__init__(Z, M, is_left, operator.mul)

    cpdef _act_(self, nn, a):
        """
        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the field
        ``GF(101)``, for otherwise it could be garbage collected, giving rise
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: GF101 = GF(101)
            sage: act = IntegerMulAction(ZZ, GF101)
            sage: act(3, 9)
            27
            sage: act(3^689, 9)
            42
            sage: 3^689 * mod(9, 101)
            42

        TESTS:

        Use round off error to verify this is doing actual repeated addition
        instead of just multiplying::

            sage: act = IntegerMulAction(ZZ, RR)
            sage: act(49, 1/49) == 49*RR(1/49)
            False

        This used to hang before :issue:`17844`::

            sage: E = EllipticCurve(GF(5), [4,0])                                       # needs sage.schemes
            sage: P = E.random_element()                                                # needs sage.schemes
            sage: (-2^63)*P                                                             # needs sage.schemes
            (0 : 1 : 0)

        Check that large multiplications can be interrupted::

            sage: # needs sage.schemes
            sage: P = E([2,1,1])
            sage: alarm(0.001); 2^(10^8) * P
            Traceback (most recent call last):
            ...
            AlarmInterrupt

        Verify that cysignals correctly detects that the above
        exception has been handled::

            sage: # needs sage.schemes
            sage: from cysignals.tests import print_sig_occurred
            sage: print_sig_occurred()
            No current exception
        """
        cdef int err = 0
        cdef long n_long

        if integer_check_long(nn, &n_long, &err) and not err:
            return fast_mul_long(a, n_long)

        return fast_mul(a, nn)

    def _repr_name_(self):
        """
        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the field
        ``GF(5)``, for otherwise it could be garbage collected, giving rise
        random errors (see :issue:`18157`). ::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: GF5 = GF(5)
            sage: IntegerMulAction(ZZ, GF5)
            Left Integer Multiplication by Integer Ring
            on Finite Field of size 5
        """
        return "Integer Multiplication"


cdef class IntegerPowAction(IntegerAction):
    r"""
    The right action ``a ^ n = a * a * ... * a`` where `n` is an
    integer.

    The action is implemented using the ``_pow_int`` method on elements.

    INPUT:

    - ``Z`` -- a type or parent representing integers

    - ``M`` -- a parent whose elements implement ``_pow_int``

    - ``m`` -- (optional) an element of ``M``

    EXAMPLES::

        sage: from sage.structure.coerce_actions import IntegerPowAction
        sage: R.<x> = LaurentSeriesRing(QQ)
        sage: act = IntegerPowAction(ZZ, R)
        sage: act(x, 5)
        x^5
        sage: act(x, -2)
        x^-2
        sage: act(x, int(5))
        x^5

    TESTS::

        sage: IntegerPowAction(ZZ, R, True)
        Traceback (most recent call last):
        ...
        ValueError: powering must be a right action
        sage: IntegerPowAction(ZZ, QQ^3)                                                # needs sage.modules
        Traceback (most recent call last):
        ...
        TypeError: no integer powering action defined on Vector space of dimension 3 over Rational Field

    ::

        sage: var('x,y')                                                                # needs sage.symbolic
        (x, y)
        sage: RDF('-2.3')^(x+y^3+sin(x))                                                # needs sage.symbolic
        (-2.3)^(y^3 + x + sin(x))
        sage: RDF('-2.3')^x                                                             # needs sage.symbolic
        (-2.3)^x
    """
    def __init__(self, Z, M, is_left=False, m=None):
        if is_left:
            raise ValueError("powering must be a right action")
        if m is None:
            m = M.an_element()
        try:
            # Check that there is a _pow_int() method
            m._pow_int
        except AttributeError:
            raise TypeError(f"no integer powering action defined on {M}")
        super().__init__(Z, M, False, operator.pow)

    cpdef _act_(self, n, a):
        """
        EXAMPLES:

        Note that coerce actions should only be used inside of the coercion
        model. For this test, we need to strongly reference the field
        ``GF(101)``::

            sage: from sage.structure.coerce_actions import IntegerPowAction
            sage: GF101 = GF(101)
            sage: act = IntegerPowAction(ZZ, GF101)
            sage: act(3, 100)
            1
            sage: act(3, -1)
            34
            sage: act(3, 1000000000000000000000000000000000000000000001)
            3
        """
        cdef Element e = <Element>a
        cdef long value = 0
        cdef int err
        integer_check_long(n, &value, &err)
        if not err:
            return e._pow_long(value)
        return e._pow_int(n)

    def _repr_name_(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerPowAction
            sage: IntegerPowAction(ZZ, QQ)
            Right Integer Powering by Integer Ring on Rational Field
        """
        return "Integer Powering"


cdef inline fast_mul(a, n):
    if n < 0:
        n = -n
        a = -a
    pow2a = a
    while n & 1 == 0:
        sig_check()
        pow2a += pow2a
        n = n >> 1
    sum = pow2a
    n = n >> 1
    while n != 0:
        sig_check()
        pow2a += pow2a
        if n & 1:
            sum += pow2a
        n = n >> 1
    return sum

cdef inline fast_mul_long(a, long s):
    # It's important to change the signed s to an unsigned n,
    # since -LONG_MIN = LONG_MIN.  See Issue #17844.
    cdef unsigned long n
    if s < 0:
        n = -s
        a = -a
    else:
        n = s
    if n < 4:
        if n == 0:
            p = parent(a)
            try:
                return p.zero()
            except AttributeError:
                return p(0)
        elif n == 1: return a
        elif n == 2: return a+a
        elif n == 3: return a+a+a
    pow2a = a
    while n & 1 == 0:
        pow2a += pow2a
        n = n >> 1
    sum = pow2a
    n = n >> 1
    while n != 0:
        pow2a += pow2a
        if n & 1:
            sum += pow2a
        n = n >> 1
    return sum
