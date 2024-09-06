r"""
Discrete dynamical systems
==========================

A *discrete dynamical system* (henceforth *DDS*) is a pair
`(X, \phi)` of a set `X` and a map `\phi : X \to X`.
(This is one of several things known as a "discrete dynamical
system" in mathematics.)

This file implements the following classes for discrete
dynamical systems:

- :class:`DiscreteDynamicalSystem`: general discrete dynamical
  system, as above.
  Inherit from this class if the ground set of your DDS is
  infinite or large enough that you want to avoid it getting
  stored as a list.
  See the doc of this class for further details.

- :class:`FiniteDynamicalSystem`: finite discrete dynamical
  system.
  This can be instantiated by calling
  :class:`DiscreteDynamicalSystem` with the parameter
  ``is_finite`` set to ``True``.

- :class:`InvertibleDiscreteDynamicalSystem`: invertible
  discrete dynamical system.
  This implements an ``inverse_evolution`` method for `\phi^{-1}`
  (the default implementation simply applies `\phi` over and
  over until the original value is revisited; the last value
  before that is then taken to be the result).
  This can be instantiated by calling
  :class:`DiscreteDynamicalSystem` with the parameter
  ``inverse`` provided.

- :class:`InvertibleFiniteDynamicalSystem`: invertible
  finite discrete dynamical system.
  This can be instantiated by calling
  :class:`DiscreteDynamicalSystem` with the parameter
  ``is_finite`` set to ``True`` and the parameter
  ``inverse`` provided.

.. TODO::

    - Implement some more functionality for homomesy and
      invariance testing:
      Checking invariance on a sublist;
      computing the first `k` entries of an orbit (useful
      when orbits can be too large);
      orbits_iterator (for when there are too many orbits to
      list);
      etc.

    - Further examples for non-auto functionality: e.g.,
      infection on a chessboard; Conway's game of life.

    - Subclasses for DDSes whose ground set is an enumerated set.
      Should we have those?

    - Implement caching for orbits (can be useful: some DDSes
      have a complicated evolution that shouldn't be recomputed
      every time).
      Does this require a whole new class?

    - Further functionality for non-invertible DDSes:
      is_recurrent, recurrent_entries, idempotent_power, etc.

    - Wrap (some of) the cyclic_sieving_phenomenon.py methods
      (:mod:`sage.combinat.cyclic_sieving_phenomenon`).

    - Interact with sage.dynamics. This requires someone who
      knows the latter part of the Sage library well.
"""
# ****************************************************************************
#       Copyright (C) 2018 Darij Grinberg <darijgrinberg@gmail.com>,
#                     2018 Tom Roby <tomrobyuconn@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.sets_cat import Sets
from sage.structure.sage_object import SageObject
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall


class DiscreteDynamicalSystem(SageObject, metaclass=ClasscallMetaclass):
    r"""
    A discrete dynamical system.

    A *discrete dynamical system* (henceforth *DDS*) is a
    pair `(X, \phi)` of a set `X` and a map `\phi : X \to X`.
    The set `X` is called the *ground set* of the DDS;
    the map `\phi` is called the *evolution* of the DDS;
    the inverse map `\phi^{-1}` (when it exists) is called the
    *inverse evolution* of the DDS.

    The DDS is called *finite* if `X` is finite.
    The DDS is called *invertible* if the map `\phi` is
    invertible.

    Given a DDS `(X, \phi)`, we can study

    * its orbits (i.e., the lists
      `(s, \phi(s), \phi^2(s), \phi^3(s), \ldots)` for `s \in X`),

    * its invariants (i.e., maps `f : X \to Y` satisfying
      `f \circ \phi = f`),

    * its cycles (i.e., lists `(u_1, u_2, \ldots, u_k)` of elements
      of `X` such that `\phi(u_i) = u_{i+1}` for each `i \leq k`,
      where we set `u_{k+1} = u_1`),

    * its homomesies (i.e., maps `h : X \to A` to a
      `\QQ`-vector space `A` such that the average of the values
      of `h` on each cycle is the same),

    and various other features.
    (Some of these require `X` to be finite or at least to have
    finite orbits.)

    INPUT:

    - ``X`` -- set, list, tuple, or another iterable, or
      ``None`` (default: ``None``); the ground set for the DDS.
      This can be ``None`` (in which case Sage will not know
      the ground set, but can still apply evolution to any
      elements that are provided to it).
      Make sure to set the ``create_tuple`` argument to
      ``True`` if the ``X`` you provide is an iterator or
      a list, as otherwise your ``X`` would be exposed
      (and thus subject to mutation or exhaustion).

    - ``phi`` -- function, or callable that acts like a
      function; the evolution of the DDS

    - ``cache_orbits`` -- boolean (default: ``False``);
      whether or not the orbits should be cached once they
      are computed.
      This currently does nothing, as we are not caching
      orbits yet.

    - ``create_tuple`` -- boolean (default: ``False``);
      whether or not the input ``X`` should be translated
      into a tuple. Set this to ``True`` to prevent
      mutation if ``X`` is a list, and to prevent
      exhaustion if ``X`` is an iterator.

    - ``inverse`` -- function, or callable that acts like a
      function, or boolean or ``None`` (default: ``None``);
      the inverse evolution of the DDS, if the DDS is
      invertible. Set this to ``None`` or ``False``
      if the DDS is not invertible (or you don't want Sage
      to treat it as such).
      Alternatively, by setting this argument to ``True``,
      you can signal that the DDS is invertible
      without providing the inverse evolution. (In this case,
      Sage will compute the inverse, assuming the orbits
      to be finite.)

    - ``is_finite`` -- boolean or ``None`` (default: ``None``);
      whether the DDS is finite. The default option ``None``
      leaves this to Sage to decide.
      Only set this to ``True`` if you provide the ground
      set ``X``.

    EXAMPLES:

    The following discrete dynamical system is neither
    finite nor invertible::

        sage: D = DiscreteDynamicalSystem(NN, lambda x: x + 2)
        sage: D.ground_set()
        Non negative integer semiring
        sage: D.evolution()(5)
        7
        sage: D.evolution_power(7)(5)
        19
        sage: D.evolution_power(0)(5)
        5

    The necessity of ``create_tuple=True``::

        sage: X = [0, 1, 2, 3, 4]
        sage: D_wrong = DiscreteDynamicalSystem(X, lambda x: (x**3) % 5)
        sage: D_right = DiscreteDynamicalSystem(X, lambda x: (x**3) % 5, create_tuple=True)
        sage: X[4] = 666 # evil
        sage: D_wrong.ground_set()
        [0, 1, 2, 3, 666]
        sage: D_right.ground_set()
        (0, 1, 2, 3, 4)

    Here is an invertible (but infinite) discrete dynamical
    system whose orbits are finite::

        sage: D = DiscreteDynamicalSystem(NN, lambda x: (x + 2 if x % 6 < 4 else x - 4), inverse=True)
        sage: D.ground_set()
        Non negative integer semiring
        sage: D.evolution()(5)
        1
        sage: D.evolution()(1)
        3
        sage: D.evolution()(3)
        5
        sage: D.evolution_power(2)(5)
        3
        sage: D.evolution_power(3)(5)
        5
        sage: D.evolution_power(-2)(5)
        1
        sage: D.inverse_evolution()(4)
        2
        sage: D.orbit(3)
        [3, 5, 1]

    Setting the ``inverse`` parameter to ``None`` or ``False``
    would give the same system without the functionality that
    relies on invertibility::

        sage: D = DiscreteDynamicalSystem(NN, lambda x: (x + 2 if x % 6 < 4 else x - 4), inverse=False)
        sage: D.ground_set()
        Non negative integer semiring
        sage: D.evolution()(5)
        1
        sage: D.inverse_evolution()(4)
        Traceback (most recent call last):
        ...
        AttributeError: 'DiscreteDynamicalSystem' object has no attribute 'inverse_evolution'...
        sage: D.orbit(3)
        [3, 5, 1]

        sage: D = DiscreteDynamicalSystem(NN, lambda x: (x + 2 if x % 6 < 4 else x - 4), inverse=None)
        sage: D.ground_set()
        Non negative integer semiring
        sage: D.evolution()(5)
        1
        sage: D.inverse_evolution()(4)
        Traceback (most recent call last):
        ...
        AttributeError: 'DiscreteDynamicalSystem' object has no attribute 'inverse_evolution'...
        sage: D.orbit(3)
        [3, 5, 1]

    Next, let us try out a finite non-invertible DDS::

        sage: D = DiscreteDynamicalSystem(tuple(range(13)), lambda x: (x**2) % 13)
        sage: D.ground_set()
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
        sage: D.evolution()(4)
        3
        sage: D.orbit(4)
        [4, 3, 9]
        sage: D.orbit(1)
        [1]
        sage: D.orbit(3)
        [3, 9]

    Note that the finiteness is automatically being inferred here,
    since the (finite) tuple ``tuple(range(13))`` has been
    provided as the ground set.

    Finally, here is a finite invertible DDS::

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (s[0],)
        sage: D = DiscreteDynamicalSystem(Y, rot, inverse=True)
        sage: D.evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
        sage: D.inverse_evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 0, 1, 1, 0, 1, 0, 0)
        sage: sorted(D.orbit_lengths())
        [2, 4, 8, 8, 8, 8, 8, 8, 8, 8]

    We could have just as well provided its inverse explicitly::

        sage: rot7 = lambda s: (s[-1],) + s[:-1]
        sage: D = DiscreteDynamicalSystem(Y, rot, inverse=rot7)
        sage: D.evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
        sage: D.inverse_evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 0, 1, 1, 0, 1, 0, 0)
    """
    @staticmethod
    def __classcall_private__(cls, X, phi, cache_orbits=False, create_tuple=False, inverse=None, is_finite=None):
        """
        Return the correct object based on input.

        The main purpose of this method is to decide which
        subclass the object will belong to based on the
        ``inverse`` and ``is_finite`` arguments.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(NN, lambda x: x + 1)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.DiscreteDynamicalSystem'>

            sage: f1 = lambda x: (x + 2 if x % 6 < 4 else x - 4)

            sage: D = DiscreteDynamicalSystem(NN, f1, inverse=False)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.DiscreteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem(NN, f1, inverse=None)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.DiscreteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem(NN, f1, inverse=True)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.InvertibleDiscreteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem(NN, lambda x: x + 1, is_finite=False)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.DiscreteDynamicalSystem'>

            sage: f2 = lambda x: (x + 1 if x < 3 else 1)
            sage: f3 = lambda x: (x - 1 if x > 3 else 3)

            sage: D = DiscreteDynamicalSystem([1, 2, 3], f2)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.FiniteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem([1, 2, 3], f2, inverse=True)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.InvertibleFiniteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem([1, 2, 3], f2, inverse=f3)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.InvertibleFiniteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem([1, 2, 3], f2, inverse=False)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.FiniteDynamicalSystem'>

            sage: D = DiscreteDynamicalSystem([1, 2, 3], f2, inverse=None)
            sage: parent(D)
            <class 'sage.dynamics.finite_dynamical_system.FiniteDynamicalSystem'>
        """
        if is_finite is None:
            is_finite = (X in Sets().Finite() or isinstance(X, (list,tuple,set,frozenset)))
        if inverse:
            if inverse is True: # invertibility claimed, but inverse not provided
                # This is how the input for these subclasses work
                inverse = None
            ret_cls = (InvertibleFiniteDynamicalSystem if is_finite
                   else InvertibleDiscreteDynamicalSystem)
            return ret_cls(X, phi, cache_orbits=cache_orbits,
                           create_tuple=create_tuple, inverse=inverse)
        if is_finite:
            return FiniteDynamicalSystem(X, phi, cache_orbits=cache_orbits,
                                         create_tuple=create_tuple)
        return typecall(cls, X, phi, cache_orbits=cache_orbits, create_tuple=create_tuple)

    def __init__(self, X, phi, cache_orbits=False, create_tuple=False):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: (3 if x == 4 else 1), create_tuple=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: (3 if x == 4 else 1), create_tuple=True, is_finite=False)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem(NN, lambda x: (3 if x == 4 else 1))
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem(None, lambda x: (3 if x == 4 else 1))
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: x, create_tuple=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
        """
        if create_tuple:
            X = tuple(X)
        self._X = X
        self._phi = phi
        self._cache_orbits = cache_orbits

    def ground_set(self):
        r"""
        Return the ground set of ``self``.

        This will return ``None`` if no ground set was
        provided in the construction of ``self``.

        .. WARNING::

            Unless ``self`` has been constructed with the
            ``create_tuple`` parameter set to ``True``,
            this method will return whatever ground set was
            provided to the constructor.
            In particular, if a list was provided, then this
            precise list will be returned; mutating this list
            will then corrupt ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: (3 if x == 4 else 1), create_tuple=True)
            sage: D.ground_set()
            (1, 3, 4)
        """
        return self._X

    def evolution(self):
        r"""
        Return the evolution of ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: (3 if x == 4 else 1), create_tuple=True)
            sage: ev = D.evolution()
            sage: ev(1)
            1
            sage: ev(4)
            3
        """
        return self._phi

    def evolution_power(self, n):
        r"""
        Return the `n`-th power (with respect to composition)
        of the evolution of ``self``.

        This requires `n` to be a nonnegative integer.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(range(10), lambda x: (x + 3) % 10, create_tuple=True)
            sage: ev3 = D.evolution_power(3)
            sage: ev3(1)
            0
            sage: ev3(2)
            1
            sage: ev0 = D.evolution_power(0)
            sage: ev0(1)
            1
            sage: ev0(2)
            2
            sage: D.evolution_power(-1)
            Traceback (most recent call last):
            ...
            ValueError: the n-th power of evolution is only defined for nonnegative integers n
        """
        from sage.rings.semirings.non_negative_integer_semiring import NN
        if n not in NN:
            raise ValueError("the n-th power of evolution is only defined for nonnegative integers n")
        ev = self.evolution()

        def evn(x):
            y = x
            for _ in range(n):
                y = ev(y)
            return y
        return evn

    def __iter__(self):
        r"""
        Iterate over the ground set of ``self``.

        This assumes that an iterable ground set of ``self``
        has been provided.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(NN, lambda x: x + 2)
            sage: D[:3]
            [0, 1, 2]
        """
        return iter(self._X)

    def __getitem__(self, i):
        r"""
        Accessor for the ``i``-th entry of the ground set of
        ``self``, assuming that the latter ground set was
        provided as a list or tuple when ``self`` was constructed.

        EXAMPLES::

            sage: X = (5, 6, 7, 8)
            sage: D = DiscreteDynamicalSystem(X, lambda x: (x**3) % 4 + 5)
            sage: X[3]
            8
            sage: X[2]
            7
        """
        return self._X[i]

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(NN, lambda x: x + 2)
            sage: D # indirect doctest
            A discrete dynamical system with ground set
             Non negative integer semiring
            sage: D = DiscreteDynamicalSystem(None, lambda x: x + 2)
            sage: D # indirect doctest
            A discrete dynamical system with unspecified ground set
        """
        if self._X is None:
            return "A discrete dynamical system with unspecified ground set"
        return "A discrete dynamical system with ground set " \
               + repr(self._X)

    def orbit(self, x, preperiod=False):
        r"""
        Return the orbit of the element ``x`` of the ground set
        of ``self`` under the evolution `\phi` of ``self``.

        This orbit is a list beginning with ``x`` and ending
        with the last element that is not a repetition of a
        previous element.
        If the orbit is infinite, then this method does not
        terminate!

        If the optional argument ``preperiod`` is set to
        ``True``, then this method returns a pair ``(o, k)``,
        where ``o`` is the orbit of ``self``, while `k` is the
        smallest nonnegative integer such that
        `\phi^k(x) \in \left\{ \phi^i(x) \mid i > k \right\}`.

        The orbit of the element ``x`` is also called the
        "rho" of ``x``, due to its shape when it is depicted
        as a directed graph.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(11)), lambda x: (x ** 2) % 11)
            sage: D.orbit(6)
            [6, 3, 9, 4, 5]
            sage: D.orbit(6, preperiod=True)
            ([6, 3, 9, 4, 5], 1)
            sage: D.orbit(3)
            [3, 9, 4, 5]
            sage: D.orbit(3, preperiod=True)
            ([3, 9, 4, 5], 0)
            sage: D.orbit(9)
            [9, 4, 5, 3]
            sage: D.orbit(0)
            [0]
        """
        orb = [x]
        phi = self._phi
        curr = phi(x)
        while curr not in orb:
            orb.append(curr)
            curr = phi(curr)
        if not preperiod:
            return orb
        return (orb, orb.index(curr))

    def is_homomesic(self, h, average=None, find_average=False, elements=None):
        r"""
        Check if ``h`` (a map from the ground set of ``self`` to
        a `\QQ`-vector space) is homomesic with respect to ``self``.

        If the optional argument ``average`` is provided, then
        this also checks that the averages are equal to ``average``.

        If the optional argument ``find_average`` is set to
        ``True``, then this method returns the average of ``h``
        in case ``h`` is homomesic (instead of returning ``True``).

        If the optional argument ``elements`` (an iterable of
        elements of the ground set of ``self``) is provided, then
        this method only checks homomesy for the cycles in the
        orbits of the elements given in the list ``elements``.
        Note that ``elements`` must be provided if the ground set of
        ``self`` is infinite (or cannot be iterated through for any
        other reason), since there is no way to check all the cycles
        in this case.

        This method will fail to terminate if any element of
        ``elements`` has an infinite orbit.

        Let us recall the definition of homomesy:
        Let `(X, \phi)` be a DDS.
        A *cycle* of `(X, \phi)` is a finite list
        `u = (u_1, u_2, \ldots, u_k)` of elements of `X` such
        that `\phi(u_i) = u_{i+1}` for each `i \leq k`, where we
        set `u_{k+1} = u_1`.
        Note that any element of `X` whose orbit is finite has a
        cycle in its orbit.
        Now, let `h` be a map from `X` to a `\QQ`-vector space `A`.
        If `u = (u_1, u_2, \ldots, u_k)` is any cycle of
        `(X, \phi)`, then the *average* of `h` on this cycle is
        defined to be the element
        `(h(u_1) + h(u_2) + \cdots + h(u_k)) / k` of `A`.
        We say that `h` is *homomesic* (with respect to the DDS
        `(X, \phi)`) if and only if the averages of `h` on all
        cycles of `(X, \phi)` are equal.

        EXAMPLES::

            sage: W = Words(2, 5)
            sage: F = DiscreteDynamicalSystem(W, lambda x: x[1:] + Word([x[0]]), is_finite=True, inverse=True)
            sage: F.is_homomesic(lambda w: sum(w))
            False
            sage: F.is_homomesic(lambda w: 1, average=1)
            True
            sage: F.is_homomesic(lambda w: 1, average=0)
            False
            sage: F.is_homomesic(lambda w: 1)
            True
            sage: F.is_homomesic(lambda w: 1, find_average=True)
            1
            sage: F.is_homomesic(lambda w: w[0] - w[1], average=0)
            True
            sage: F.is_homomesic(lambda w: w[0] - w[1], find_average=True)
            0

        Now, let us check homomesy restricted to specific cycles::

            sage: F = finite_dynamical_systems.bitstring_rotation(7)
            sage: descents = lambda x: sum(1 for i in range(6) if x[i] > x[i+1])
            sage: F.is_homomesic(descents)
            False
            sage: F.is_homomesic(descents, elements=[(1, 0, 1, 0, 0, 0, 0), (1, 0, 0, 1, 0, 0, 0)])
            True
            sage: F.is_homomesic(descents, elements=[(1, 0, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0, 0)])
            False
            sage: F.is_homomesic(descents, elements=[(1, 0, 1, 0, 0, 0, 0)])
            True
            sage: F.is_homomesic(descents, elements=[])
            True

        And here is a non-invertible finite dynamical system::

            sage: F = finite_dynamical_systems.one_line([9, 1, 1, 6, 5, 4, 5, 5, 1])
            sage: F.is_homomesic(lambda i: i)
            True
            sage: F.is_homomesic(lambda i: i % 2)
            False
            sage: F.is_homomesic(lambda i: i % 2, elements=[2, 9, 7])
            True
            sage: F.is_homomesic(lambda i: i % 2, elements=[2, 9, 4])
            False
            sage: F.is_homomesic(lambda i: i % 2, elements=[2, 9, 5, 7, 8, 2])
            True
        """
        from sage.rings.rational_field import QQ
        orbavgs = [] # This will be the list of all averages on cycles.
        if elements is None:
            # The user has not provided elements, so we need to
            # check all cycles of the DDS.
            for cyc in self.cycles():
                l = len(cyc)
                avg = ~(QQ(l)) * sum(h(i) for i in cyc)
                if avg not in orbavgs:
                    if orbavgs:
                        return False
                    orbavgs.append(avg)
        else:
            # Checking only the cycles of the elements provided
            # by the user.
            for element in elements:
                (orb, ix) = self.orbit(element, preperiod=True)
                cyc = orb[ix:] # the cycle in the orbit of element
                l = len(cyc)
                avg = ~(QQ(l)) * sum(h(i) for i in cyc)
                if avg not in orbavgs:
                    if orbavgs:
                        return False
                    orbavgs.append(avg)
        if not orbavgs:
            return True
        if average is None:
            if find_average:
                return orbavgs[0]
            return True
        return orbavgs[0] == average

class InvertibleDiscreteDynamicalSystem(DiscreteDynamicalSystem):
    r"""
    An invertible discrete dynamical system.

    A *discrete dynamical system* (henceforth *DDS*) is a
    pair `(X, \phi)` of a set `X` and a map `\phi : X \to X`.
    This set `X` is called the *ground set* of the DDS, while
    the map `\phi` is called the *evolution* of the DDS.
    An *invertible DDS* is a DDS `(X, \phi)` whose evolution
    `\phi` is invertible.
    In that case, `\phi^{-1}` is called the *inverse evolution*
    of the DDS.

    See :class:`DiscreteDynamicalSystem` for details.

    INPUT:

    - ``X`` -- set, list, tuple, or another iterable, or
      ``None``; the ground set for the DDS. This can be
      ``None`` in case of a
      :class:`DiscreteDynamicalSystem` or a
      :class:`InvertibleDiscreteDynamicalSystem`.
      Make sure to set the ``create_tuple`` argument to
      ``True`` if you provide an iterator or a list for
      ``X``, as otherwise the input would be exposed.

    - ``phi`` -- function, or callable that acts like a
      function; the evolution of the DDS

    - ``inverse`` -- function, or callable that acts like a
      function; the inverse evolution of the DDS. (A
      default implementation is implemented when this
      argument is not provided; but it assumes the orbits
      to be finite.)

    - ``cache_orbits`` -- boolean (default: ``False``);
      whether or not the orbits should be cached once they
      are computed.

    - ``create_tuple`` -- boolean (default: ``False``);
      whether or not the input ``X`` should be translated
      into a tuple (set this to ``True`` to prevent
      mutation if ``X`` is a list, and to prevent
      exhaustion if ``X`` is an iterator).

    EXAMPLES::

        sage: from sage.dynamics.finite_dynamical_system import InvertibleDiscreteDynamicalSystem
        sage: D = InvertibleDiscreteDynamicalSystem(NN, lambda x: (x + 2 if x % 4 < 2 else x - 2))
        sage: D.ground_set()
        Non negative integer semiring
        sage: D.evolution()(5)
        7
        sage: D.evolution()(6)
        4
        sage: D.evolution()(4)
        6
        sage: D.inverse_evolution()(4)
        6

    The necessity of ``create_tuple=True``::

        sage: X = [0, 1, 2, 3, 4]
        sage: D_wrong = InvertibleDiscreteDynamicalSystem(X, lambda x: (x**3) % 5)
        sage: D_right = InvertibleDiscreteDynamicalSystem(X, lambda x: (x**3) % 5, create_tuple=True)
        sage: X[4] = 666 # evil
        sage: D_wrong.ground_set()
        [0, 1, 2, 3, 666]
        sage: D_right.ground_set()
        (0, 1, 2, 3, 4)
    """
    def __init__(self, X, phi, inverse=None, cache_orbits=False, create_tuple=False):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: x, create_tuple=True, inverse=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x: x, create_tuple=True, is_finite=False, inverse=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem(NN, lambda x: x, inverse=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
            sage: D = DiscreteDynamicalSystem(None, lambda x: x, inverse=True)
            sage: TestSuite(D).run(skip ='_test_pickling') # indirect doctest
        """
        if create_tuple:
            X = tuple(X)
        self._X = X
        self._phi = phi
        if inverse is None:
            self._inverse = self.inverse_evolution_default
        else:
            self._inverse = inverse
        self._cache_orbits = cache_orbits

    def evolution_power(self, n):
        r"""
        Return the `n`-th power (with respect to composition)
        of the evolution of ``self``.

        This requires `n` to be an integer.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(range(10), lambda x: (x + 3) % 10, create_tuple=True, inverse=True)
            sage: ev3 = D.evolution_power(3)
            sage: ev3(1)
            0
            sage: ev3(2)
            1
            sage: ev0 = D.evolution_power(0)
            sage: ev0(1)
            1
            sage: ev0(2)
            2
            sage: evm1 = D.evolution_power(-1)
            sage: evm1(1)
            8
            sage: evm1(2)
            9
            sage: evm2 = D.evolution_power(-2)
            sage: evm2(1)
            5
            sage: evm2(2)
            6
        """
        from sage.rings.integer_ring import ZZ
        if n not in ZZ:
            raise ValueError("the n-th power of evolution is only defined for integers n")
        if n >= 0:
            ev = self.evolution()
        else:
            ev = self.inverse_evolution()
            n = -n

        def evn(x):
            y = x
            for _ in range(n):
                y = ev(y)
            return y
        return evn

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(NN, lambda x: (x + 2 if x % 4 < 2 else x - 2), inverse=True)
            sage: D # indirect doctest
            An invertible discrete dynamical system with ground set
             Non negative integer semiring
            sage: D = DiscreteDynamicalSystem(None, lambda x: x + 2, inverse=True)
            sage: D # indirect doctest
            An invertible discrete dynamical system with unspecified ground set
        """
        if self._X is None:
            return "An invertible discrete dynamical system with unspecified ground set"
        return "An invertible discrete dynamical system with ground set " \
               + repr(self._X)

    def inverse_evolution(self):
        r"""
        Return the inverse evolution of ``self`` (as a map
        from the ground set of ``self`` to itself).

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(8)), lambda x: (x + 2) % 8, inverse=True)
            sage: D.inverse_evolution()(1)
            7
            sage: D.inverse_evolution()(3)
            1

            sage: D = DiscreteDynamicalSystem(ZZ, lambda x: (x + 2) % 8, inverse=True)
            sage: D.inverse_evolution()(1)
            7
            sage: D.inverse_evolution()(3)
            1
        """
        return self._inverse

    def verify_inverse_evolution(self, x=None):
        r"""
        Verify that the composition of evolution and
        inverse evolution on ``self`` is the identity
        (both ways).

        The optional argument ``x``, if provided, restricts
        the testing to the element ``x`` only.
        Otherwise, all elements of the ground set are
        tested (if they can be enumerated).

        This is mostly used to check the correctness of
        self-implemented inverse evolution methods.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(8)), lambda x: (x + 2) % 8, inverse=True)
            sage: D.verify_inverse_evolution()
            True
            sage: D.verify_inverse_evolution(3)
            True
            sage: fake_inverse = lambda x: x
            sage: D = DiscreteDynamicalSystem(tuple(range(8)), lambda x: (x + 2) % 8, inverse=fake_inverse)
            sage: D.verify_inverse_evolution()
            False
            sage: D.verify_inverse_evolution(3)
            False
        """
        ev = self.evolution()
        iev = self.inverse_evolution()
        if x is None:
            els = self.ground_set()
        else:
            els = [x]
        for y in els:
            if ev(iev(y)) != y:
                return False
            if iev(ev(y)) != y:
                return False
        return True

    def orbit(self, x, preperiod=False):
        r"""
        Return the orbit of the element ``x`` of the ground set
        of ``self``.

        This orbit is a list beginning with ``x`` and ending
        with the last element until ``x`` reappears.
        If ``x`` never reappears, then this will not terminate!

        If the optional argument ``preperiod`` is set to
        ``True``, then this method returns a pair ``(o, k)``,
        where ``o`` is the orbit of ``self``, while `k` is the
        smallest nonnegative integer such that
        `\phi^k(x) \in \left\{ \phi^i(x) \mid i > k \right\}`.
        Note that `k` is necessarily `0`, since the DDS
        ``self`` is invertible!

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(8)), lambda x: (x + 2) % 8, inverse=True)
            sage: D.ground_set()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: D.orbit(2)
            [2, 4, 6, 0]
            sage: D.orbit(5)
            [5, 7, 1, 3]
            sage: D.orbit(5, preperiod=True)
            ([5, 7, 1, 3], 0)

            sage: D = DiscreteDynamicalSystem(ZZ, lambda x: (x + 2) % 8, inverse=True)
            sage: D.ground_set()
            Integer Ring
            sage: D.orbit(2)
            [2, 4, 6, 0]
            sage: D.orbit(5)
            [5, 7, 1, 3]
            sage: D.orbit(5, preperiod=True)
            ([5, 7, 1, 3], 0)
        """
        orb = [x]
        phi = self._phi
        curr = phi(x)
        while curr != x:
            orb.append(curr)
            curr = phi(curr)
        if not preperiod:
            return orb
        return (orb, 0)

    def inverse_evolution_default(self, x):
        r"""
        Return the inverse evolution of ``self``, applied
        to the element ``x`` of the ground set of ``self``.

        This is the default implementation, assuming that
        the orbit of ``x`` is finite.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(8)), lambda x: (x + 2) % 8, inverse=True)
            sage: D.inverse_evolution_default(1)
            7
            sage: D.inverse_evolution_default(3)
            1
        """
        return self.orbit(x)[-1]

class FiniteDynamicalSystem(DiscreteDynamicalSystem):
    r"""
    A finite discrete dynamical system.

    A *finite discrete dynamical system* (henceforth *FDDS*) is a
    pair `(X, \phi)` of a finite set `X` and a map `\phi : X \to X`.
    This set `X` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.

    The ground set `X` should always be provided as an
    iterable when defining a :class:`FiniteDynamicalSystem`.

    EXAMPLES::

        sage: from sage.dynamics.finite_dynamical_system import FiniteDynamicalSystem
        sage: D = FiniteDynamicalSystem(tuple(range(11)), lambda x: (x**2) % 11)
        sage: D.ground_set()
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        sage: D.evolution()(4)
        5
        sage: D.orbit(4)
        [4, 5, 3, 9]
        sage: D.orbit(1)
        [1]
        sage: D.orbit(2)
        [2, 4, 5, 3, 9]

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (0,)
        sage: D = FiniteDynamicalSystem(Y, rot)
        sage: D.evolution()((1, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
    """
    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(11)), lambda x: (x**2) % 11)
            sage: D
            A finite discrete dynamical system with ground set
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        """
        return "A finite discrete dynamical system with ground set " \
               + repr(self._X)

    def is_invariant(self, f):
        r"""
        Check if ``f`` is an invariant of ``self``.

        Let `(X, \phi)` be a discrete dynamical system.
        Let `Y` be any set. Let `f : X \to Y` be any map.
        Then, we say that `f` is an *invariant* of `(X, \phi)`
        if and only if `f \circ \phi = f`.

        EXAMPLES::

            sage: W = Words(2, 5)
            sage: F = DiscreteDynamicalSystem(W, lambda x: x[1:] + Word([x[0]]), is_finite=True)
            sage: F.is_invariant(lambda w: sum(w))
            True
            sage: F.is_invariant(lambda w: 1)
            True
            sage: F.is_invariant(lambda w: w[0] - w[1])
            False
            sage: F.is_invariant(lambda w: sum(i**2 for i in w))
            True

        Invariants and non-invariants of a permutation::

            sage: F = finite_dynamical_systems.permutation([3, 4, 5, 6, 1, 2])
            sage: F.is_invariant(lambda i: i % 2)
            True
            sage: F.is_invariant(lambda i: i % 3)
            False
            sage: F.is_invariant(lambda i: i > 1)
            False
            sage: F.is_invariant(lambda i: i % 2 == 0)
            True
        """
        phi = self._phi
        return all(f(phi(i)) == f(i) for i in self._X)

    def cycles(self):
        r"""
        Return a list of all cycles of ``self``, up to
        cyclic rotation.

        We recall the definition of cycles:
        Let `(X, \phi)` be a DDS.
        A *cycle* of `(X, \phi)` is a finite list
        `u = (u_1, u_2, \ldots, u_k)` of elements of `X` such
        that `\phi(u_i) = u_{i+1}` for each `i \leq k`, where we
        set `u_{k+1} = u_1`.
        Note that any element of `X` whose orbit is finite has a
        cycle in its orbit.

        EXAMPLES::

            sage: BS = finite_dynamical_systems.bulgarian_solitaire
            sage: BS(8).cycles()
            [[[4, 3, 1], [3, 3, 2], [3, 2, 2, 1], [4, 2, 1, 1]],
             [[4, 2, 2], [3, 3, 1, 1]]]
            sage: BS(6).cycles()
            [[[3, 2, 1]]]

            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 2) % 6)
            sage: D.cycles()
            [[5, 1, 3], [4, 0, 2]]
            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x ** 2) % 6)
            sage: D.cycles()
            [[1], [4], [3], [0]]
            sage: D = DiscreteDynamicalSystem(tuple(range(11)), lambda x: (x ** 2 - 1) % 11)
            sage: D.cycles()
            [[10, 0], [8], [4]]

            sage: F = finite_dynamical_systems.one_line([4, 7, 2, 6, 2, 10, 9, 11, 5, 6, 12, 12, 12, 6])
            sage: F.cycles()
            [[6, 10], [12], [9, 5, 2, 7]]
        """
        l = list(self)
        cycs = []
        while l:
            start = l[-1]
            (orb, ix) = self.orbit(start, preperiod=True)
            if orb[ix] in l:
                # This means we've actually found a new cycle,
                # not just a new path to an old cycle.
                cycs.append(orb[ix:])
            for j in orb:
                try:
                    l.remove(j)
                except ValueError:
                    # Here we break out of the for-loop, because
                    # if ``j`` has already been removed from
                    # ``l``, then all later elements of the orbit
                    # must have been removed from ``l`` as well
                    # (indeed, the set of elements that have been
                    # removed from ``l`` is closed under ``phi``).
                    break
        return cycs

class InvertibleFiniteDynamicalSystem(InvertibleDiscreteDynamicalSystem, FiniteDynamicalSystem):
    r"""
    An invertible finite discrete dynamical system.

    A *finite discrete dynamical system* (henceforth *FDDS*) is a
    pair `(X, \phi)` of a finite set `X` and a map `\phi : X \to X`.
    This set `X` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.
    An FDDS `(X, \phi)` is called *invertible* if the map `\phi`
    is invertible; in this case, `\phi^{-1}` is called the
    *inverse evolution* of the FDDS.

    The ground set `X` should always be provided as an
    iterable when defining a :class:`FiniteDynamicalSystem`.

    EXAMPLES::

        sage: from sage.dynamics.finite_dynamical_system import InvertibleFiniteDynamicalSystem
        sage: D = InvertibleFiniteDynamicalSystem(tuple(range(5)), lambda x: (x + 2) % 5)
        sage: D.ground_set()
        (0, 1, 2, 3, 4)
        sage: D.evolution()(4)
        1
        sage: D.orbits()
        [[4, 1, 3, 0, 2]]
        sage: D.inverse_evolution()(2)
        0
        sage: D.inverse_evolution()(1)
        4
        sage: D.evolution_power(-1)(1)
        4
        sage: D.evolution_power(-2)(1)
        2

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (s[0],)
        sage: D = InvertibleFiniteDynamicalSystem(Y, rot)
        sage: D.evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
        sage: D.inverse_evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 0, 1, 1, 0, 1, 0, 0)
        sage: sorted(D.orbit_lengths())
        [2, 4, 8, 8, 8, 8, 8, 8, 8, 8]
    """
    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(5)), lambda x: (x + 2) % 5, inverse=True)
            sage: D
            An invertible finite discrete dynamical system with ground set
            (0, 1, 2, 3, 4)
        """
        return "An invertible finite discrete dynamical system with ground set " \
               + repr(self._X)

    def orbits(self):
        r"""
        Return a list of all orbits of ``self``, up to
        cyclic rotation.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 2) % 6, inverse=True)
            sage: D.orbits()
            [[5, 1, 3], [4, 0, 2]]
            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 3) % 6, inverse=True)
            sage: D.orbits()
            [[5, 2], [4, 1], [3, 0]]
        """
        phi = self._phi
        l = list(self)
        orbs = []
        while l:
            start = l.pop()
            orb = [start]
            curr = phi(start)
            while curr != start:
                l.remove(curr)
                orb.append(curr)
                curr = phi(curr)
            orbs.append(orb)
        return orbs

    def cycles(self):
        r"""
        Return a list of all cycles of ``self``, up to
        cyclic rotation.

        We recall the definition of cycles:
        Let `(X, \phi)` be a DDS.
        A *cycle* of `(X, \phi)` is a finite list
        `u = (u_1, u_2, \ldots, u_k)` of elements of `X` such
        that `\phi(u_i) = u_{i+1}` for each `i \leq k`, where we
        set `u_{k+1} = u_1`.
        Note that any element of `X` whose orbit is finite has a
        cycle in its orbit.

        Since ``self`` is invertible, the cycles of ``self``
        are the same as its orbits.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 2) % 6, inverse=True)
            sage: D.cycles()
            [[5, 1, 3], [4, 0, 2]]
            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 3) % 6, inverse=True)
            sage: D.cycles()
            [[5, 2], [4, 1], [3, 0]]
        """
        return self.orbits()

    def orbit_lengths(self):
        r"""
        Return a list of the lengths of all orbits of
        ``self``.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 2) % 6, inverse=True)
            sage: D.orbit_lengths()
            [3, 3]
            sage: D = DiscreteDynamicalSystem(tuple(range(6)), lambda x: (x + 3) % 6, inverse=True)
            sage: D.orbit_lengths()
            [2, 2, 2]
        """
        return [len(orb) for orb in self.orbits()]
