r"""
Group, ring, etc. actions on objects

The terminology and notation used is suggestive of groups acting on sets,
but this framework can be used for modules, algebras, etc.

A group action `G \times S \rightarrow S` is a functor from `G` to Sets.

.. WARNING::

    An :class:`Action` object only keeps a weak reference to the underlying
    set which is acted upon. This decision was made in :issue:`715` in order
    to allow garbage collection within the coercion framework (this is where
    actions are mainly used) and avoid memory leaks.

    ::

        sage: from sage.categories.action import Action
        sage: class P: pass
        sage: A = Action(P(),P())
        sage: import gc
        sage: _ = gc.collect()
        sage: A
        <repr(<sage.categories.action.Action at 0x...>) failed:
         RuntimeError: This action acted on a set that became garbage collected>

    To avoid garbage collection of the underlying set, it is sufficient to
    create a strong reference to it before the action is created.

    ::

        sage: _ = gc.collect()
        sage: from sage.categories.action import Action
        sage: class P: pass
        sage: q = P()
        sage: A = Action(P(),q)
        sage: gc.collect()
        0
        sage: A
        Left action by <__main__.P ... at ...> on <__main__.P ... at ...>

    Once an :class:`Action` is stored in a coercion cache, the cache layer
    also weakens its reference to the actor ``G`` (:issue:`27358`); actions
    explicitly registered via :meth:`~sage.structure.parent.Parent.register_action`
    opt out via :meth:`Action._pin_actor` so the registered parent's
    ``_action_list`` keeps the actor alive.

AUTHOR:

- Robert Bradshaw: initial version
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.tuple cimport PyTuple_GET_ITEM

from sage.categories.functor cimport Functor
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map
from sage.structure.element cimport parent
from sage.structure.parent cimport Parent

from sage.categories import homset
from weakref import ref


cdef inline category(x):
    try:
        return x.category()
    except AttributeError:
        from sage.categories.objects import Objects
        return Objects()


cdef class Action(Functor):
    """
    The action of ``G`` on ``S``.

    INPUT:

    - ``G`` -- a parent or Python type

    - ``S`` -- a parent or Python type

    - ``is_left`` -- boolean (default: ``True``); whether elements of
      ``G`` are on the left

    - ``op`` -- (default: ``None``) operation. This is not used by
      :class:`Action` itself, but other classes may use it

    TESTS:

    Check that repeated binary operations involving large numbers of
    parents do not leak those parents into memory (:issue:`27358`)::

        sage: import gc, operator
        sage: for p in primes(2000):
        ....:     F = GF(p)
        ....:     M = MatrixSpace(F, 2)
        ....:     V = F^2
        ....:     _ = M.one() * V.zero()
        sage: _ = gc.collect(); _ = gc.collect()
        sage: n = len(list(primes(2000)))
        sage: len([v for v in gc.get_objects() if type(v) is type(V)]) < n
        True
    """
    def __init__(self, G, S, is_left=True, op=None):
        from sage.categories.objects import Objects
        # Pass the singleton :class:`Objects` category as the Functor
        # domain rather than ``Groupoid(G)``: the latter would pin
        # ``G`` in memory forever through the
        # :class:`UniqueRepresentation` cache key of ``Groupoid(G)``
        # (see :issue:`27358`).
        Functor.__init__(self, Objects(), category(S))
        # Hold ``G`` strongly by default; the coercion cache layer
        # (see :meth:`_make_actor_weak`) explicitly weakens this
        # reference once it has stored the action, so cached actions
        # do not pin their acting parents (:issue:`27358`).  Actions
        # registered via :meth:`Parent.register_action` opt out via
        # :meth:`_pin_actor`.  ``S`` is always held weakly
        # (:issue:`715`).
        self._G = G
        self._G_weakref = None
        self._pinned = False
        self.US = ref(S)
        self._is_left = is_left
        self.op = op

    @property
    def G(self):
        """
        The actor of this action.

        EXAMPLES::

            sage: from sage.categories.action import Action
            sage: class P: pass
            sage: g = P()
            sage: s = P()
            sage: A = Action(g, s)
            sage: A.G is g
            True

        For cached actions the actor is held only weakly (see
        :issue:`27358`); accessing :attr:`G` after the actor has been
        garbage-collected raises :exc:`RuntimeError`.
        """
        return self._actor()

    cdef _actor(self):
        """
        Return the actor parent ``G``.

        Prefers the strong slot ``_G``; falls back to the weak slot
        ``_G_weakref`` set by :meth:`_make_actor_weak`.  Raises
        :exc:`RuntimeError` if the actor has been garbage-collected.
        """
        cdef G
        if self._G is not None:
            return self._G
        if self._G_weakref is not None:
            G = self._G_weakref()
            if G is not None:
                return G
        raise RuntimeError("This action acted with an actor that became garbage collected")

    cpdef _pin_actor(self):
        """
        Mark this action so that :meth:`_make_actor_weak` becomes a
        no-op.

        Called by :meth:`Parent.register_action` because the parent's
        ``_action_list`` is the only thing keeping the actor alive for
        a registered action; weakening would let the actor be
        garbage-collected while the parent is still using the action
        (see :issue:`27358`).
        """
        self._pinned = True

    cpdef _make_actor_weak(self):
        """
        Convert the actor reference from strong to weak.

        No-op if the action has been pinned via :meth:`_pin_actor`, if
        the actor is already weakly held, or if the actor is not
        weak-referenceable (some Python built-in types).  Used by the
        coercion cache to avoid pinning actor parents through cached
        actions (:issue:`27358`).
        """
        if self._pinned or self._G is None:
            return
        try:
            self._G_weakref = ref(self._G)
        except TypeError:
            # Not weak-referenceable; keep the strong reference.
            return
        self._G = None

    def _apply_functor(self, x):
        return self(x)

    def __reduce__(self):
        """
        Used in pickling.

        .. WARNING::

            If you change the signature of the ``__init__`` for a subclass,
            you must override this method as well.

        TESTS:

        Check that this action can be pickled (:issue:`29031`)::

            sage: P = QQ['x']
            sage: R = (ZZ['x'])['y']
            sage: A = R.get_action(P, operator.mul, True)
            sage: loads(dumps(A)) is not None
            True
        """
        return (type(self), (self.G, self.underlying_set(), self._is_left, self.op))

    def __call__(self, *args):
        """
        Let this action act.

        For a left action, ``action(a, b)`` lets ``a`` act on ``b``.
        For a right action, ``action(a, b)`` lets ``b`` act on ``a``.

        If needed, ``a`` and ``b`` are converted to the correct parent.

        .. SEEALSO::

            :meth:`act` which lets you pass the acting and acted-on
            elements directly.

        EXAMPLES::

            sage: R.<x> = ZZ []
            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: A = IntegerMulAction(ZZ, R, True)   # Left action
            sage: A(5, x)
            5*x
            sage: A(int(5), x)
            5*x
            sage: A(x, 5)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial
            sage: A = IntegerMulAction(ZZ, R, False)  # Right action
            sage: A(x, 5)
            5*x
            sage: A(x, int(5))
            5*x
            sage: A(5, x)
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial
        """
        if len(args) == 2:
            # Normal case, called with (g, x) or (x, g) as arguments
            #
            # For a left action:  g = args[0] and x = args[1]
            # For a right action: g = args[1] and x = args[0]
            g = <object>PyTuple_GET_ITEM(args, 1 - self._is_left)
            x = <object>PyTuple_GET_ITEM(args, self._is_left)
            return self._act_convert(g, x)
        elif len(args) == 1:
            g = <object>PyTuple_GET_ITEM(args, 0)
            G = self._actor()
            if g in G:
                return ActionEndomorphism(self, G(g))
            elif g == G:
                return self.underlying_set()
            else:
                raise TypeError("%s not an element of %s" % (g, G))
        else:
            raise TypeError("actions should be called with 1 or 2 arguments")

    cdef _act_convert(self, g, x):
        """
        Let ``g`` act on ``x`` under this action, converting ``g``
        and ``x`` to the correct parents first.
        """
        cdef G = self._actor()
        U = self.underlying_set()
        if parent(g) is not G:
            g = G(g)
        if parent(x) is not U:
            x = U(x)
        return self._act_(g, x)

    cpdef _act_(self, g, x):
        """
        Let ``g`` act on ``x`` under this action.

        Regardless of whether this is a left or right action, the acting
        element comes first.

        INPUT:

        - ``g`` -- an object with parent ``self.G``

        - ``x`` -- an object with parent ``self.US()``

        .. WARNING::

            This is meant to be a fast internal function, so the
            conditions on the input are not checked!
        """
        raise NotImplementedError(f"action for {type(self)} not implemented")

    def act(self, g, x):
        """
        This is a consistent interface for acting on ``x`` by ``g``,
        regardless of whether it's a left or right action.

        If needed, ``g`` and ``x`` are converted to the correct parent.

        EXAMPLES::

            sage: R.<x> = ZZ []
            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: A = IntegerMulAction(ZZ, R, True)   # Left action
            sage: A.act(5, x)
            5*x
            sage: A.act(int(5), x)
            5*x
            sage: A = IntegerMulAction(ZZ, R, False)  # Right action
            sage: A.act(5, x)
            5*x
            sage: A.act(int(5), x)
            5*x
        """
        return self._act_convert(g, x)

    def __invert__(self):
        return InverseAction(self)

    def is_left(self):
        return self._is_left

    def _repr_(self):
        side = "Left" if self._is_left else "Right"
        return "%s %s by %r on %r" % (side, self._repr_name_(), self.G,
                                      self.underlying_set())

    def _repr_name_(self):
        return "action"

    def actor(self):
        return self.G

    cdef underlying_set(self):
        """
        The set on which the actor acts (it is not necessarily the codomain of
        the action).

        .. NOTE::

            Since this is a cdef'ed method, we can only provide an indirect doctest.

        EXAMPLES::

            sage: P = QQ['x']
            sage: R = (ZZ['x'])['y']
            sage: A = R.get_action(P,operator.mul,True)
            sage: A                 # indirect doctest
            Right scalar multiplication
             by Univariate Polynomial Ring in x over Rational Field
             on Univariate Polynomial Ring in y over
                 Univariate Polynomial Ring in x over Integer Ring

        In this example, the underlying set is the ring ``R``. This is the same
        as the left domain, which is different from the codomain of the action::

            sage: A.codomain()
            Univariate Polynomial Ring in y
             over Univariate Polynomial Ring in x over Rational Field
            sage: A.codomain() == R
            False
            sage: A.left_domain() is R
            True

        By :issue:`715`, there is only a weak reference to the underlying set.
        Hence, the underlying set may be garbage collected, even when the
        action is still alive. This may result in a runtime error, as follows::

            sage: from sage.categories.action import Action
            sage: class P: pass
            sage: p = P()
            sage: q = P()
            sage: A = Action(p,q)
            sage: A
            Left action by <__main__.P ... at ...> on <__main__.P ... at ...>
            sage: del q
            sage: import gc
            sage: _ = gc.collect()
            sage: A
            <repr(<sage.categories.action.Action at 0x...>) failed:
             RuntimeError: This action acted on a set that became garbage collected>
        """
        S = self.US()
        if S is None:
            raise RuntimeError("This action acted on a set that became garbage collected")
        return S

    def codomain(self):
        return self.underlying_set()

    def domain(self):
        return self.underlying_set()

    def left_domain(self):
        if self._is_left:
            return self.G
        return self.domain()

    def right_domain(self):
        if self._is_left:
            return self.domain()
        return self.G

    def operation(self):
        return self.op


cdef class InverseAction(Action):
    """
    An action that acts as the inverse of the given action.

    EXAMPLES::

        sage: V = QQ^3                                                                  # needs sage.modules
        sage: v = V((1, 2, 3))                                                          # needs sage.modules
        sage: cm = get_coercion_model()

        sage: # needs sage.modules
        sage: a = cm.get_action(V, QQ, operator.mul)
        sage: a
        Right scalar multiplication by Rational Field
         on Vector space of dimension 3 over Rational Field
        sage: ~a
        Right inverse action by Rational Field
         on Vector space of dimension 3 over Rational Field
        sage: (~a)(v, 1/3)
        (3, 6, 9)

        sage: # needs sage.modules
        sage: b = cm.get_action(QQ, V, operator.mul)
        sage: b
        Left scalar multiplication by Rational Field
         on Vector space of dimension 3 over Rational Field
        sage: ~b
        Left inverse action by Rational Field
         on Vector space of dimension 3 over Rational Field
        sage: (~b)(1/3, v)
        (3, 6, 9)

        sage: c = cm.get_action(ZZ, list, operator.mul)
        sage: c
        Left action by Integer Ring on <... 'list'>
        sage: ~c
        Traceback (most recent call last):
        ...
        TypeError: no inverse defined for Left action by Integer Ring on <... 'list'>

    TESTS:

        sage: x = polygen(QQ,'x')
        sage: a = 2*x^2+2; a
        2*x^2 + 2
        sage: a / 2
        x^2 + 1
        sage: a /= 2
        sage: a
        x^2 + 1
    """
    def __init__(self, Action action):
        G = action.G
        try:
            from sage.groups.group import Group
            # We must be in the case that parent(~a) == parent(a)
            # so we can invert in _call_ code below.
            if (isinstance(G, Group) and G.is_multiplicative()) or G.is_field():
                Action.__init__(self, G, action.underlying_set(), action._is_left)
                self._action = action
                return
        except (AttributeError, NotImplementedError):
            pass
        raise TypeError(f"no inverse defined for {action!r}")

    def __reduce__(self):
        """
        Used in pickling.

        TESTS:

        Check that this action can be pickled (:issue:`29031`)::

            sage: # needs sage.modules
            sage: V = QQ^3
            sage: v = V((1, 2, 3))
            sage: cm = get_coercion_model()
            sage: a = cm.get_action(V, QQ, operator.mul)
            sage: loads(dumps(~a)) is not None
            True
        """
        return (type(self), (self._action,))

    cpdef _act_(self, g, x):
        if self.S_precomposition is not None:
            x = self.S_precomposition(x)
        return self._action._act_(~g, x)

    def codomain(self):
        return self._action.codomain()

    def __invert__(self):
        return self._action

    def _repr_name_(self):
        return "inverse action"


cdef class PrecomposedAction(Action):
    """
    A precomposed action first applies given maps, and then applying an action
    to the return values of the maps.

    EXAMPLES:

    We demonstrate that an example discussed on :issue:`14711` did not become a
    problem::

        sage: # needs sage.libs.flint sage.modular
        sage: E = ModularSymbols(11).2
        sage: s = E.modular_symbol_rep()
        sage: del E,s
        sage: import gc
        sage: _ = gc.collect()
        sage: E = ModularSymbols(11).2
        sage: v = E.manin_symbol_rep()
        sage: c,x = v[0]
        sage: y = x.modular_symbol_rep()
        sage: coercion_model.get_action(QQ, parent(y), op=operator.mul)
        Left scalar multiplication by Rational Field
         on Abelian Group of all Formal Finite Sums over Rational Field
         with precomposition on right by Coercion map:
          From: Abelian Group of all Formal Finite Sums over Integer Ring
          To:   Abelian Group of all Formal Finite Sums over Rational Field
    """
    def __init__(self, Action action, Map left_precomposition, Map right_precomposition):
        left = action.left_domain()
        right = action.right_domain()
        US = action.underlying_set()
        cdef Parent lco, rco
        if left_precomposition is not None:
            lco = left_precomposition._codomain
            if lco is not left:
                left_precomposition = homset.Hom(lco, left).natural_map() * left_precomposition
            left = left_precomposition.domain()
        if right_precomposition is not None:
            rco = right_precomposition._codomain
            if rco is not right:
                right_precomposition = homset.Hom(rco, right).natural_map() * right_precomposition
            right = right_precomposition.domain()
        if action._is_left:
            Action.__init__(self, left, US, 1)
        else:
            Action.__init__(self, right, US, 0)
        self._action = action
        if self._is_left:
            self.G_precomposition = left_precomposition
            self.S_precomposition = right_precomposition
        else:
            self.G_precomposition = right_precomposition
            self.S_precomposition = left_precomposition

    def __reduce__(self):
        """
        Used in pickling.

        TESTS:

        Check that this action can be pickled (:issue:`29031`)::

            sage: # needs sage.libs.flint sage.modular
            sage: E = ModularSymbols(11).2
            sage: v = E.manin_symbol_rep()
            sage: c,x = v[0]
            sage: y = x.modular_symbol_rep()
            sage: act = coercion_model.get_action(QQ, parent(y), op=operator.mul)
            sage: loads(dumps(act)) is not None
            True
        """
        return (type(self), (self._action, self.G_precomposition, self.S_precomposition))

    cpdef _act_(self, g, x):
        if self.G_precomposition is not None:
            g = self.G_precomposition._call_(g)
        if self.S_precomposition is not None:
            x = self.S_precomposition._call_(x)
        return self._action._act_(g, x)

    def domain(self):
        if self.S_precomposition is not None:
            return self.S_precomposition.domain()
        return self._action.domain()

    def codomain(self):
        return self._action.codomain()

    @property
    def left_precomposition(self):
        """
        The left map to precompose with, or ``None`` if there is no
        left precomposition map.
        """
        if self._is_left:
            return self.G_precomposition
        return self.S_precomposition

    @property
    def right_precomposition(self):
        """
        The right map to precompose with, or ``None`` if there is no
        right precomposition map.
        """
        if self._is_left:
            return self.S_precomposition
        return self.G_precomposition

    def __invert__(self):
        return PrecomposedAction(~self._action, self.left_precomposition, self.right_precomposition)

    def _repr_(self):
        s = repr(self._action)
        if self.left_precomposition is not None:
            s += "\nwith precomposition on left by %s" % self.left_precomposition._default_repr_()
        if self.right_precomposition is not None:
            s += "\nwith precomposition on right by %s" % self.right_precomposition._default_repr_()
        return s


cdef class ActionEndomorphism(Morphism):
    """
    The endomorphism defined by the action of one element.

    EXAMPLES::

        sage: A = ZZ['x'].get_action(QQ, self_on_left=False, op=operator.mul)
        sage: A
        Left scalar multiplication by Rational Field
         on Univariate Polynomial Ring in x over Integer Ring
        sage: A(1/2)
        Action of 1/2 on Univariate Polynomial Ring in x over Integer Ring
        under Left scalar multiplication by Rational Field on Univariate
        Polynomial Ring in x over Integer Ring.
    """
    def __init__(self, Action action, g):
        Morphism.__init__(self, homset.Hom(action.underlying_set(),
                                           action.underlying_set()))
        self._action = action
        self._g = g

    cdef dict _extra_slots(self):
        """
        Helper for pickling and copying.

        TESTS::

            sage: P.<x> = ZZ[]
            sage: A = P.get_action(QQ, self_on_left=False, op=operator.mul)
            sage: phi = A(1/2)
            sage: psi = copy(phi)  # indirect doctest
            sage: psi
            Action of 1/2 on Univariate Polynomial Ring in x over
            Integer Ring under Left scalar multiplication by Rational
            Field on Univariate Polynomial Ring in x over Integer Ring.
            sage: psi(x) == phi(x)
            True
        """
        slots = Morphism._extra_slots(self)
        slots['_action'] = self._action
        slots['_g'] = self._g
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for pickling and copying.

        TESTS::

            sage: P.<x> = ZZ[]
            sage: A = P.get_action(QQ, self_on_left=False, op=operator.mul)
            sage: phi = A(1/2)
            sage: psi = copy(phi)  # indirect doctest
            sage: psi
            Action of 1/2 on Univariate Polynomial Ring in x over
            Integer Ring under Left scalar multiplication by Rational
            Field on Univariate Polynomial Ring in x over Integer Ring.
            sage: psi(x) == phi(x)
            True
        """
        self._action = _slots['_action']
        self._g = _slots['_g']
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        return self._action._act_(self._g, x)

    def _repr_(self):
        return "Action of %s on %s under %s." % (self._g,
                                                 self._action.underlying_set(),
                                                 self._action)

    def __mul__(left, right):
        cdef ActionEndomorphism left_c, right_c
        if isinstance(left, ActionEndomorphism) and isinstance(right, ActionEndomorphism):
            left_c = left
            right_c = right
            if left_c._action is right_c._action:
                if left_c._action._is_left:
                    return ActionEndomorphism(left_c._action, left_c._g * right_c._g)
                return ActionEndomorphism(left_c._action, right_c._g * left_c._g)
        return Morphism.__mul__(left, right)

    def __invert__(self):
        inv_g = ~self._g
        if parent(inv_g) is parent(self._g):
            return ActionEndomorphism(self._action, inv_g)
        return (~self._action)(self._g)
