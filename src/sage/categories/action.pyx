r"""
Group, ring, etc. actions on objects.

The terminology and notation used is suggestive of groups acting on sets,
but this framework can be used for modules, algebras, etc.

A group action $G \times S \rightarrow S$ is a functor from $G$ to Sets.

.. WARNING::

    An :class:`Action` object only keeps a weak reference to the underlying set
    which is acted upon. This decision was made in :trac:`715` in order to
    allow garbage collection within the coercion framework (this is where
    actions are mainly used) and avoid memory leaks.

    ::

        sage: from sage.categories.action import Action
        sage: class P: pass
        sage: A = Action(P(),P())
        sage: import gc
        sage: _ = gc.collect()
        sage: A
        <repr(<sage.categories.action.Action at 0x...>) failed: RuntimeError: This action acted on a set that became garbage collected>

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
        Left action by <__main__.P instance at ...> on <__main__.P instance at ...>

AUTHOR:

- Robert Bradshaw: initial version
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functor cimport Functor
from morphism cimport Morphism
from map cimport Map
from sage.structure.parent cimport Parent

import homset
import sage.structure.element
from weakref import ref
from sage.misc.constant_function import ConstantFunction


cdef inline category(x):
    try:
        return x.category()
    except AttributeError:
        import sage.categories.all
        return sage.categories.all.Objects()

cdef class Action(Functor):

    def __init__(self, G, S, bint is_left = 1, op=None):
        from groupoid import Groupoid
        Functor.__init__(self, Groupoid(G), category(S))
        self.G = G
        self.US = ref(S)
        self._is_left = is_left
        self.op = op

    def _apply_functor(self, x):
        return self(x)

    def __call__(self, *args):
        if len(args) == 1:
            g = args[0]
            if g in self.G:
                return ActionEndomorphism(self, self.G(g))
            elif g == self.G:
                return self.underlying_set()
            else:
                raise TypeError("%s not an element of %s" % (g, self.G))
        elif len(args) == 2:
            if self._is_left:
                return self._call_(self.G(args[0]), self.underlying_set()(args[1]))
            else:
                return self._call_(self.underlying_set()(args[0]), self.G(args[1]))

    cpdef _call_(self, a, b):
        raise NotImplementedError("Action not implemented.")

    def act(self, g, a):
        """
        This is a consistent interface for acting on a by g,
        regardless of whether it's a left or right action.
        """
        if self._is_left:
            return self._call_(g, a)
        else:
            return self._call_(a, g)

    def __invert__(self):
        return InverseAction(self)

    def is_left(self):
        return self._is_left

    def _repr_(self):
        side = "Left" if self._is_left else "Right"
        return "%s %s by %r on %r"%(side, self._repr_name_(), self.G,
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
            Right scalar multiplication by Univariate Polynomial Ring in x over
            Rational Field on Univariate Polynomial Ring in y over Univariate
            Polynomial Ring in x over Integer Ring

        In this example, the underlying set is the ring ``R``. This is the same
        as the left domain, which is different from the codomain of the action::

            sage: A.codomain()
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: A.codomain() == R
            False
            sage: A.left_domain() is R
            True

        By :trac:`715`, there is only a weak reference to the underlying set.
        Hence, the underlying set may be garbage collected, even when the
        action is still alive. This may result in a runtime error, as follows::

            sage: from sage.categories.action import Action
            sage: class P: pass
            sage: p = P()
            sage: q = P()
            sage: A = Action(p,q)
            sage: A
            Left action by <__main__.P instance at ...> on <__main__.P instance at ...>
            sage: del q
            sage: import gc
            sage: _ = gc.collect()
            sage: A
            <repr(<sage.categories.action.Action at 0x...>) failed: RuntimeError: This action acted on a set that became garbage collected>
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
        else:
            return self.domain()

    def right_domain(self):
        if self._is_left:
            return self.domain()
        else:
            return self.G

    def operation(self):
        return self.op


cdef class InverseAction(Action):
    """
    An action that acts as the inverse of the given action.

    TESTS:

    This illustrates a shortcoming in the current coercion model.
    See the comments in _call_ below::

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
            from sage.groups.group import is_Group
            # We must be in the case that parent(~a) == parent(a)
            # so we can invert in call_c code below.
            if (is_Group(G) and G.is_multiplicative()) or G.is_field():
                Action.__init__(self, G, action.underlying_set(), action._is_left)
                self._action = action
                return
            else:
                K = G._pseudo_fraction_field()
                Action.__init__(self, K, action.underlying_set(), action._is_left)
                self._action = action
                return
        except (AttributeError, NotImplementedError):
            pass
        raise TypeError("No inverse defined for %r." % action)

    cpdef _call_(self, a, b):
        if self._action._is_left:
            if self.S_precomposition is not None:
                b = self.S_precomposition(b)
            return self._action._call_(~a, b)
        else:
            if self.S_precomposition is not None:
                a = self.S_precomposition(a)
            return self._action._call_(a, ~b)

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

    We demonstrate that an example discussed on :trac:`14711` did not become a
    problem::

        sage: E = ModularSymbols(11).2
        sage: s = E.modular_symbol_rep()
        sage: del E,s
        sage: import gc
        sage: _ = gc.collect()
        sage: E = ModularSymbols(11).2
        sage: v = E.manin_symbol_rep()
        sage: c,x = v[0]
        sage: y = x.modular_symbol_rep()
        sage: A = y.parent().get_action(QQ, self_on_left=False, op=operator.mul)
        sage: A
        Left scalar multiplication by Rational Field on Abelian Group of all
        Formal Finite Sums over Rational Field
        with precomposition on right by Conversion map:
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
        self.left_precomposition = left_precomposition
        self.right_precomposition = right_precomposition

    cpdef _call_(self, a, b):
        if self.left_precomposition is not None:
            a = self.left_precomposition._call_(a)
        if self.right_precomposition is not None:
            b = self.right_precomposition._call_(b)
        return self._action._call_(a, b)

    def domain(self):
        if self._is_left and self.right_precomposition is not None:
            return self.right_precomposition.domain()
        elif not self._is_left and self.left_precomposition is not None:
            return self.left_precomposition.domain()
        else:
            return self._action.domain()

    def codomain(self):
        return self._action.codomain()

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
        Left scalar multiplication by Rational Field on Univariate Polynomial
        Ring in x over Integer Ring
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

    cdef dict _extra_slots(self, dict _slots):
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
        _slots['_action'] = self._action
        _slots['_g'] = self._g
        return Morphism._extra_slots(self, _slots)

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
        if self._action._is_left:
            return self._action._call_(self._g, x)
        else:
            return self._action._call_(x, self._g)

    def _repr_(self):
        return "Action of %s on %s under %s."%(self._g,
                                               self._action.underlying_set(), self._action)

    def __mul__(left, right):
        cdef ActionEndomorphism left_c, right_c
        if isinstance(left, ActionEndomorphism) and isinstance(right, ActionEndomorphism):
            left_c = left
            right_c = right
            if left_c._action is right_c._action:
                if left_c._action._is_left:
                    return ActionEndomorphism(left_c._action, left_c._g * right_c._g)
                else:
                    return ActionEndomorphism(left_c._action, right_c._g * left_c._g)
        return Morphism.__mul__(left, right)

    def __invert__(self):
            inv_g = ~self._g
            if sage.structure.element.parent(inv_g) is sage.structure.element.parent(self._g):
                return ActionEndomorphism(self._action, inv_g)
            else:
                return (~self._action)(self._g)


