"""
Soliton Cellular Automata

AUTHORS:

- Travis Scrimshaw (2017-06-30): Initial version
- Travis Scrimshaw (2018-02-03): Periodic version
"""

# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
from sage.combinat.root_system.cartan_type import CartanType
from sage.typeset.ascii_art import ascii_art
from sage.rings.integer_ring import ZZ


class SolitonCellularAutomata(SageObject):
    r"""
    Soliton cellular automata.

    Fix an affine Lie algebra `\mathfrak{g}` with index `I` and
    classical index set `I_0`. Fix some `r \in I_0`. A *soliton
    cellular automaton* (SCA) is a discrete (non-linear) dynamical
    system given as follows. The *states* are given by elements of
    a semi-infinite tensor product of Kirillov-Reshetihkin crystals
    `B^{r,1}`, where only a finite number of factors are not the
    maximal element `u`, which we will call the *vacuum*. The *time
    evolution* `T_s` is defined by

    .. MATH::

        R(p \otimes u_s) = u_s \otimes T_s(p),

    where `p = \cdots \otimes p_3 \otimes p_2 \otimes p_1 \otimes p_0`
    is a state and `u_s` is the maximal element of `B^{r,s}`.
    In more detail, we have `R(p_i \otimes u^{(i)}) =
    u^{(i+1)} \otimes \widetilde{p}_i` with `u^{(0)} = u_s` and
    `T_s(p) = \cdots \otimes \widetilde{p}_1 \otimes \widetilde{p}_0`.
    This is well-defined since `R(u \otimes u_s) = u_s \otimes u`
    and `u^{(k)} = u_s` for all `k \gg 1`.

    INPUT:

    - ``initial_state`` -- the list of elements, can also be a string
      when ``vacuum`` is 1 and ``n`` is `\mathfrak{sl}_n`
    - ``cartan_type`` -- (default: 2) the value ``n``, for `\mathfrak{sl}_n`,
      or a Cartan type
    - ``r`` -- (default: 1) the node index `r`; typically this
      corresponds to the height of the vacuum element

    EXAMPLES:

    We first create an example in `\mathfrak{sl}_4` (type `A_3`)::

        sage: B = SolitonCellularAutomata('3411111122411112223', 4)
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          34......224....2223
          evoltuions: []
          current state:
          34......224....2223

    We then apply a standard evolution::

        sage: B.evolve()
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          34......224....2223
          evoltuions: [(1, 19)]
          current state:
          .................34.....224...2223....

    Next, we apply a smaller carrier evolution. Note that the soliton
    of size 4 moves only 3 steps::

        sage: B.evolve(3)
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          34......224....2223
          evoltuions: [(1, 19), (1, 3)]
          current state:
          ...............34....224...2223.......

    We can also use carriers corresponding to non-vacuum indices.
    In these cases, the carrier might not return to its initial
    state, which results in a message being displayed about
    the resulting state of the carrier::

        sage: B.evolve(carrier_capacity=7, carrier_index=3)
        Last carrier:
          1  1  1  1  1  1  1
          2  2  2  2  2  3  3
          3  3  3  3  3  4  4
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          34......224....2223
          evoltuions: [(1, 19), (1, 3), (3, 7)]
          current state:
          .....................23....222....2223.......

        sage: B.evolve(carrier_capacity=3, carrier_index=2)
        Last carrier:
          1  1  1
          2  2  3
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          34......224....2223
          evoltuions: [(1, 19), (1, 3), (3, 7), (2, 3)]
          current state:
          .......................22.....223...2222........

    To summarize our current evolutions, we can use :meth:`print_states`::

        sage: B.print_states(5)
        t: 0
             .............................34......224....2223
        t: 1
             ...........................34.....224...2223....
        t: 2
             .........................34....224...2223.......
        t: 3
             ........................23....222....2223.......
        t: 4
             .......................22.....223...2222........

    To run the SCA further under the standard evolutions, one can use
    :meth:`print_states` or :meth:`latex_states`::

        sage: B.print_states(15)
        t: 0
             ................................................34......224....2223
        t: 1
             ..............................................34.....224...2223....
        t: 2
             ............................................34....224...2223.......
        t: 3
             ...........................................23....222....2223.......
        t: 4
             ..........................................22.....223...2222........
        t: 5
             ........................................22....223..2222............
        t: 6
             ......................................22..2223..222................
        t: 7
             ..................................2222..23...222...................
        t: 8
             ..............................2222....23..222......................
        t: 9
             ..........................2222......23.222.........................
        t: 10
              ......................2222.......223.22............................
        t: 11
              ..................2222........223..22..............................
        t: 12
              ..............2222.........223...22................................
        t: 13
              ..........2222..........223....22..................................
        t: 14
              ......2222...........223.....22....................................

    Next, we use `r = 2` in type `A_3`. Here, we give the data as lists of
    values corresponding to the entries of the column of height 2 from
    the largest entry to smallest. Our columns are drawn in French
    convention::

        sage: B = SolitonCellularAutomata([[4,1],[4,1],[2,1],[2,1],[2,1],[2,1],[3,1],[3,1],[3,2]], 4, 2)

    We perform 3 evolutions and obtain the following::

        sage: B.evolve(number=3)
        sage: B
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 2
          initial state:
          44    333
          11....112
          evoltuions: [(2, 9), (2, 9), (2, 9)]
          current state:
             44 333
          ...11.112.........

    We construct Example 2.9 from [LS2017]_::

        sage: B = SolitonCellularAutomata([[2],[-3],[1],[1],[1],[4],[0],[-2],
        ....:   [1],[1],[1],[1],[3],[-4],[-3],[-3],[1]], ['D',5,2])
        sage: B.print_states(10)
        t: 0                                    _     _     ___
             ..................................23...402....3433.
        t: 1                                  _    _    ___
             ................................23..402...3433.....
        t: 2                                _   _   ___
             ..............................23.402..3433.........
        t: 3                              _  _  ___
             ...........................243.02.3433.............
        t: 4                           _  __ __
             .......................2403..42333.................
        t: 5                       _   ___ _
             ...................2403...44243....................
        t: 6                   _    ___  _
             ...............2403....442.43......................
        t: 7               _     ___   _
             ...........2403.....442..43........................
        t: 8           _      ___    _
             .......2403......442...43..........................
        t: 9       _       ___     _
             ...2403.......442....43............................

    Example 3.4 from [LS2017]_::

        sage: B = SolitonCellularAutomata([['E'],[1],[1],[1],[3],[0],
        ....: [1],[1],[1],[1],[2],[-3],[-1],[1]], ['D',4,2])
        sage: B.print_states(10)
        t: 0                                                      __
             ..........................................E...30....231.
        t: 1                                                  __
             .........................................E..30..231.....
        t: 2                                            _  _
             ........................................E303.21.........
        t: 3                                       _ _
             ....................................303E2.22............
        t: 4                                   _    _
             ................................303E...222..............
        t: 5                               _       _
             ............................303E......12................
        t: 6                           _         _
             ........................303E........1.2.................
        t: 7                       _           _
             ....................303E..........1..2..................
        t: 8                   _             _
             ................303E............1...2...................
        t: 9               _               _
             ............303E..............1....2....................

    Example 3.12 from [LS2017]_::

        sage: B = SolitonCellularAutomata([[-1,3,2],[3,2,1],[3,2,1],[-3,2,1],
        ....:   [-2,-3,1]], ['B',3,1], 3)
        sage: B.print_states(6)
                                           -1    -3-2
        t: 0                                3     2-3
              . . . . . . . . . . . . . . . 2 . . 1 1
                                         -1-3-2
        t: 1                              3 2-3
              . . . . . . . . . . . . . . 2 1 1 . . .
                                     -3-1
        t: 2                          2-2
              . . . . . . . . . . . . 1-3 . . . . . .
                               -3-1  -3
        t: 3                    2-2   2
              . . . . . . . . . 1 3 . 1 . . . . . . .
                         -3-1      -3
        t: 4              2-2       2
              . . . . . . 1 3 . . . 1 . . . . . . . .
                   -3-1          -3
        t: 5        2-2           2
              . . . 1 3 . . . . . 1 . . . . . . . . .

    Example 4.12 from [LS2017]_::

        sage: K = crystals.KirillovReshetikhin(['E',6,1], 1,1, 'KR')
        sage: u = K.module_generators[0]
        sage: x = u.f_string([1,3,4,5])
        sage: y = u.f_string([1,3,4,2,5,6])
        sage: a = u.f_string([1,3,4,2])
        sage: B = SolitonCellularAutomata([a, u,u,u, x,y], ['E',6,1], 1)
        sage: B
        Soliton cellular automata of type ['E', 6, 1] and vacuum = 1
          initial state:
              (-2, 5)          .          .          . (-5, 2, 6)(-2, -6, 4)
          evoltuions: []
          current state:
              (-2, 5)          .          .          . (-5, 2, 6)(-2, -6, 4)
        sage: B.print_states(8)
        t: 0 ...
        t: 7
                          .       (-2, 5)(-2, -5, 4, 6) ... (-6, 2) ...
    """
    def __init__(self, initial_state, cartan_type=2, vacuum=1):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = SolitonCellularAutomata('3411111122411112223', 4)
            sage: TestSuite(B).run()
        """
        if cartan_type in ZZ:
            cartan_type = CartanType(['A',cartan_type-1,1])
        else:
            cartan_type = CartanType(cartan_type)
        self._cartan_type = cartan_type
        self._vacuum = vacuum
        K = KirillovReshetikhinTableaux(self._cartan_type, self._vacuum, 1)
        try:
            # FIXME: the maximal_vector() does not work in type E and F
            self._vacuum_elt = K.maximal_vector()
        except (ValueError, TypeError, AttributeError):
            self._vacuum_elt = K.module_generators[0]

        if isinstance(initial_state, str):
            # We consider things 1-9
            initial_state = [[ZZ(x) if x != '.' else ZZ.one()] for x in initial_state]
        try:
            KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type,
                                                             [[vacuum, len(st)//vacuum]
                                                              for st in initial_state])
            self._states = [KRT(pathlist=initial_state)]
        except TypeError:
            KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type,
                                                             [[vacuum, 1]
                                                              for st in initial_state])
            self._states = [KRT(*initial_state)]

        self._evolutions = []
        self._initial_carrier = []
        self._nballs = len(self._states[0])

    def __eq__(self, other):
        """
        Check equality.

        Two SCAs are equal when they have the same initial state
        and evolutions.

        TESTS::

            sage: B1 = SolitonCellularAutomata('34112223', 4)
            sage: B2 = SolitonCellularAutomata('34112223', 4)
            sage: B1 == B2
            True
            sage: B1.evolve()
            sage: B1 == B2
            False
            sage: B2.evolve()
            sage: B1 == B2
            True
            sage: B1.evolve(5)
            sage: B2.evolve(6)
            sage: B1 == B2
            False
        """
        return (isinstance(other, SolitonCellularAutomata)
                and self._states[0] == other._states[0]
                and self._evolutions == other._evolutions)

    def __ne__(self, other):
        """
        Check non equality.

        TESTS::

            sage: B1 = SolitonCellularAutomata('34112223', 4)
            sage: B2 = SolitonCellularAutomata('34112223', 4)
            sage: B1 != B2
            False
            sage: B1.evolve()
            sage: B1 != B2
            True
            sage: B2.evolve()
            sage: B1 != B2
            False
            sage: B1.evolve(5)
            sage: B2.evolve(6)
            sage: B1 != B2
            True
        """
        return not (self == other)

    # Evolution functions
    # -------------------

    def evolve(self, carrier_capacity=None, carrier_index=None, number=None):
        r"""
        Evolve ``self``.

        Time evolution `T_s` of a SCA state `p` is determined by

        .. MATH::

            u_{r,s} \otimes T_s(p) = R(p \otimes u_{r,s}),

        where `u_{r,s}` is the maximal element of `B^{r,s}`.

        INPUT:

        - ``carrier_capacity`` -- (default: the number of balls in
          the system) the size `s` of carrier

        - ``carrier_index`` -- (default: the vacuum index) the index `r`
          of the carrier

        - ``number`` -- (optional) the number of times to perform
          the evolutions

        To perform multiple evolutions of the SCA, ``carrier_capacity``
        and ``carrier_index`` may be lists of the same length.

        EXAMPLES::

            sage: B = SolitonCellularAutomata('3411111122411112223', 4)
            sage: for k in range(10):
            ....:     B.evolve()
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(1, 19), (1, 19), (1, 19), (1, 19), (1, 19),
                           (1, 19), (1, 19), (1, 19), (1, 19), (1, 19)]
              current state:
              ......2344.......222....23...............................

            sage: B.reset()
            sage: B.evolve(number=10); B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(1, 19), (1, 19), (1, 19), (1, 19), (1, 19),
                           (1, 19), (1, 19), (1, 19), (1, 19), (1, 19)]
              current state:
              ......2344.......222....23...............................

            sage: B.reset()
            sage: B.evolve(carrier_capacity=[1,2,3,4,5,6,7,8,9,10]); B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
                           (1, 6), (1, 7), (1, 8), (1, 9), (1, 10)]
              current state:
              ........2344....222..23..............................

            sage: B.reset()
            sage: B.evolve(carrier_index=[1,2,3])
            Last carrier:
              1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
              2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3  4  4
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(1, 19), (2, 19), (3, 19)]
              current state:
              ..................................22......223...2222.....

            sage: B.reset()
            sage: B.evolve(carrier_capacity=[1,2,3], carrier_index=[1,2,3])
            Last carrier:
              1  1
              3  4
            Last carrier:
              1  1  1
              2  2  3
              3  3  4
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(1, 1), (2, 2), (3, 3)]
              current state:
              .....22.......223....2222..

            sage: B.reset()
            sage: B.evolve(1, 2, number=3)
            Last carrier:
              1
              3
            Last carrier:
              1
              4
            Last carrier:
              1
              3
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: [(2, 1), (2, 1), (2, 1)]
              current state:
              .24......222.....2222.
        """
        if isinstance(carrier_capacity, (list, tuple)):
            if not isinstance(carrier_index, (list, tuple)):
                carrier_index = [carrier_index] * len(carrier_capacity)
            if len(carrier_index) != len(carrier_capacity):
                raise ValueError("carrier_index and carrier_capacity"
                                 " must have the same length")
            for i, r in zip(carrier_capacity, carrier_index):
                self.evolve(i, r)
            return
        if isinstance(carrier_index, (list, tuple)):
            # carrier_capacity must be not be a list/tuple if given
            for r in carrier_index:
                self.evolve(carrier_capacity, r)
            return

        if carrier_capacity is None:
            carrier_capacity = self._nballs
        if carrier_index is None:
            carrier_index = self._vacuum

        if number is not None:
            for _ in range(number):
                self.evolve(carrier_capacity, carrier_index)
            return

        passed = False
        K = KirillovReshetikhinTableaux(self._cartan_type, carrier_index, carrier_capacity)
        try:
            # FIXME: the maximal_vector() does not work in type E and F
            empty_carrier = K.maximal_vector()
        except (ValueError, TypeError, AttributeError):
            empty_carrier = K.module_generators[0]
        self._initial_carrier.append(empty_carrier)
        carrier_factor = (carrier_index, carrier_capacity)
        last_final_carrier = empty_carrier
        state = self._states[-1]
        dims = state.parent().dims
        while not passed:
            KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type,
                                                             dims + (carrier_factor,))
            elt = KRT(*(list(state) + [empty_carrier]))
            RC = RiggedConfigurations(self._cartan_type, (carrier_factor,) + dims)
            elt2 = RC(*elt.to_rigged_configuration()).to_tensor_product_of_kirillov_reshetikhin_tableaux()
            # Back to an empty carrier or we are not getting any better
            if elt2[0] == empty_carrier or elt2[0] == last_final_carrier:
                passed = True
                KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type, dims)
                self._states.append(KRT(*elt2[1:]))
                self._evolutions.append(carrier_factor)
                if elt2[0] != empty_carrier:
                    print("Last carrier:")
                    print(ascii_art(last_final_carrier))
            else:
                # We need to add more vacuum states
                last_final_carrier = elt2[0]
                dims = tuple([(self._vacuum, 1)]*carrier_capacity) + dims

    def state_evolution(self, num):
        """
        Return a list of the carrier values at state ``num`` evolving to
        the next state.

        If ``num`` is greater than the number of states, this performs
        the standard evolution `T_k`, where `k` is the number of balls
        in the system.

        .. SEEALSO::

            :meth:`print_state_evolution`, :meth:`latex_state_evolution`

        EXAMPLES::

            sage: B = SolitonCellularAutomata('1113123', 3)
            sage: B.evolve(3)
            sage: B.state_evolution(0)
            [[[1, 1, 1]],
             [[1, 1, 1]],
             [[1, 1, 1]],
             [[1, 1, 3]],
             [[1, 1, 2]],
             [[1, 2, 3]],
             [[1, 1, 3]],
             [[1, 1, 1]]]
            sage: B.state_evolution(2)
            [[[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 3]],
             [[1, 1, 1, 1, 1, 3, 3]],
             [[1, 1, 1, 1, 1, 1, 3]],
             [[1, 1, 1, 1, 1, 1, 2]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]],
             [[1, 1, 1, 1, 1, 1, 1]]]
        """
        if num + 2 > len(self._states):
            for _ in range(num + 2 - len(self._states)):
                self.evolve()

        carrier = KirillovReshetikhinTableaux(self._cartan_type, *self._evolutions[num])
        num_factors = len(self._states[num+1])
        vacuum = self._vacuum_elt
        state = [vacuum]*(num_factors - len(self._states[num])) + list(self._states[num])
        final = []
        u = [self._initial_carrier[num]]
        # Assume every element has the same parent
        R = state[0].parent().R_matrix(carrier)
        for elt in reversed(state):
            up, eltp = R(R.domain()(elt, u[0]))
            u.insert(0, up)
            final.insert(0, eltp)
        return u

    def reset(self):
        r"""
        Reset ``self`` back to the initial state.

        EXAMPLES::

            sage: B = SolitonCellularAutomata('34111111224', 4)
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224
              evoltuions: []
              current state:
              34......224
            sage: B.evolve()
            sage: B.evolve()
            sage: B.evolve()
            sage: B.evolve()
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224
              evoltuions: [(1, 11), (1, 11), (1, 11), (1, 11)]
              current state:
              ...34..224............
            sage: B.reset()
            sage: B
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224
              evoltuions: []
              current state:
              34......224
        """
        self._states = [self._states[0]]
        self._evolutions = []
        self._initial_carrier = []

    # Output functions
    # ----------------

    def _column_repr(self, b, vacuum_letter=None):
        """
        Return a string representation of the column ``b``.

        EXAMPLES::

            sage: B = SolitonCellularAutomata([[-2,1],[2,1],[-3,1],[-3,2]], ['D',4,2], 2)
            sage: K = crystals.KirillovReshetikhin(['D',4,2], 2,1, 'KR')
            sage: B._column_repr(K(-2,1))
            -2
             1
            sage: B._column_repr(K.module_generator())
            2
            1
            sage: B._column_repr(K.module_generator(), 'x')
            x
        """
        if vacuum_letter is not None and b == self._vacuum_elt:
            return ascii_art(vacuum_letter)
        if self._vacuum_elt.parent()._tableau_height == 1:
            s = str(b[0])
            return ascii_art(s if s[0] != '-' else '_\n' + s[1:])
        letter_str = [str(letter) for letter in b]
        max_width = max(len(s) for s in letter_str)
        return ascii_art('\n'.join(' '*(max_width-len(s)) + s for s in letter_str))

    def _repr_state(self, state, vacuum_letter='.'):
        """
        Return a string representation of ``state``.

        EXAMPLES::

            sage: B = SolitonCellularAutomata('3411111122411112223', 4)
            sage: B.evolve(number=10)
            sage: print(B._repr_state(B._states[0]))
            34......224....2223
            sage: print(B._repr_state(B._states[-1], '_'))
            ______2344_______222____23_______________________________
        """
        output = [self._column_repr(b, vacuum_letter) for b in state]
        max_width = max(cell.width() for cell in output)
        return sum((ascii_art(' '*(max_width-b.width())) + b for b in output),
                    ascii_art(''))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SolitonCellularAutomata('3411111122411112223', 4)
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
              initial state:
              34......224....2223
              evoltuions: []
              current state:
              34......224....2223
            sage: SolitonCellularAutomata([[4,1],[2,1],[2,1],[3,1],[3,2]], 4, 2)
            Soliton cellular automata of type ['A', 3, 1] and vacuum = 2
              initial state:
              4  33
              1..12
              evoltuions: []
              current state:
              4  33
              1..12
            sage: SolitonCellularAutomata([[4,1],[2,1],[2,1],[3,1],[3,2]], ['C',4,1], 2)
            Soliton cellular automata of type ['C', 4, 1] and vacuum = 2
              initial state:
              4  33
              1..12
              evoltuions: []
              current state:
              4  33
              1..12
            sage: SolitonCellularAutomata([[4,3],[2,1],[-3,1],[-3,2]], ['B',4,1], 2)
            Soliton cellular automata of type ['B', 4, 1] and vacuum = 2
              initial state:
               4  -3-3
               3 . 1 2
              evoltuions: []
              current state:
               4  -3-3
               3 . 1 2
        """
        ret = "Soliton cellular automata of type {} and vacuum = {}\n".format(self._cartan_type, self._vacuum)
        ret += "  initial state:\n{}\n  evoltuions: {}\n  current state:\n{}".format(
            ascii_art('  ') + self._repr_state(self._states[0]),
            self._evolutions,
            ascii_art('  ') + self._repr_state(self._states[-1])
            )
        return ret

    def print_state(self, num=None, vacuum_letter='.', remove_trailing_vacuums=False):
        """
        Print the state ``num``.

        INPUT:

        - ``num`` -- (default: the current state) the state to print
        - ``vacuum_letter`` -- (default: ``'.'``) the letter to print
          for the vacuum
        - ``remove_trailing_vacuums`` -- boolean (default: ``False``); if ``True``
          then this does not print the vacuum letters at the right end
          of the state

        EXAMPLES::

            sage: B = SolitonCellularAutomata('3411111122411112223', 4)
            sage: B.print_state()
            34......224....2223
            sage: B.evolve(number=2)
            sage: B.print_state(vacuum_letter=',')
            ,,,,,,,,,,,,,,,34,,,,224,,2223,,,,,,,,
            sage: B.print_state(10, '_')
            ______2344_______222____23_______________________________
            sage: B.print_state(10, '_', True)
            ______2344_______222____23
        """
        if num is None:
            num = len(self._states) - 1
        if num + 1 > len(self._states):
            for _ in range(num + 1 - len(self._states)):
                self.evolve()
        state = self._states[num]
        if remove_trailing_vacuums:
            pos = len(state) - 1
            # The pos goes negative if and only if the state consists
            #   entirely of vacuum elements.
            while pos >= 0 and state[pos] == self._vacuum_elt:
                pos -= 1
            state = state[:pos+1]
        print(self._repr_state(state, vacuum_letter))

    def print_states(self, num=None, vacuum_letter='.'):
        r"""
        Print the first ``num`` states of ``self``.

        .. NOTE::

            If the number of states computed for ``self`` is less than
            ``num``, then this evolves the system using the default
            time evolution.

        INPUT:

        - ``num`` -- the number of states to print

        EXAMPLES::

            sage: B = SolitonCellularAutomata([[2],[-1],[1],[1],[1],[1],[2],[2],[3],
            ....:   [-2],[1],[1],[2],[-1],[1],[1],[1],[1],[1],[1],[2],[3],[3],[-3],[-2]],
            ....:   ['C',3,1])
            sage: B.print_states(7)
            t: 0                           _       _   _         __
                 .........................21....2232..21......23332
            t: 1                        _      _    _       __
                 ......................21...2232...21....23332.....
            t: 2                     _     _     _     __
                 ...................21..2232....21..23332..........
            t: 3                  _    _      _   __
                 ...............221..232...2231..332...............
            t: 4              _     _    _     __
                 ...........221...232.2231....332..................
            t: 5          _     __          __
                 .......221...2321223......332.....................
            t: 6      _    __            __
                 ..2221...321..223......332........................

            sage: B = SolitonCellularAutomata([[2],[1],[1],[1],[3],[-2],[1],[1],
            ....:   [1],[2],[2],[-3],[1],[1],[1],[1],[1],[1],[2],[3],[3],[-3]],
            ....:   ['B',3,1])
            sage: B.print_states(9, ' ')
            t: 0                            _     _         _
                                       2   32   223      2333
            t: 1                          _    _        _
                                      2  32  223     2333
            t: 2                        _   _       _
                                     2 32 223    2333
            t: 3                      _  _      _
                                   23 2223   2333
            t: 4                     __     _
                                 23 213  2333
            t: 5                _  _    _
                             2233 222 333
            t: 6            _    _   _
                         2233  23223 3
            t: 7        _     _     _
                     2233   232 23  3
            t: 8    _      _       _
                 2233    232  23   3

            sage: B = SolitonCellularAutomata([[2],[-2],[1],[1],[1],[1],[2],[0],[-3],
            ....:   [1],[1],[1],[1],[1],[2],[2],[3],[-3],], ['D',4,2])
            sage: B.print_states(10)
            t: 0                                      _      _        _
                 ....................................22....203.....2233
            t: 1                                    _     _       _
                 ..................................22...203....2233....
            t: 2                                  _    _      _
                 ................................22..203...2233........
            t: 3                                _   _     _
                 ..............................22.203..2233............
            t: 4                              _  _    _
                 ............................22203.2233................
            t: 5                            _ _   _
                 ........................220223.233....................
            t: 6                        _   _  _
                 ....................2202.223.33.......................
            t: 7                    _    _   _
                 ................2202..223..33.........................
            t: 8                _     _    _
                 ............2202...223...33...........................
            t: 9            _      _     _
                 ........2202....223....33.............................

        Example 4.13 from [Yamada2007]_::

            sage: B = SolitonCellularAutomata([[3],[3],[1],[1],[1],[1],[2],[2],[2]], ['D',4,3])
            sage: B.print_states(15)
            t: 0
                 ....................................33....222
            t: 1
                 ..................................33...222...
            t: 2
                 ................................33..222......
            t: 3
                 ..............................33.222.........
            t: 4
                 ............................33222............
            t: 5
                 ..........................3022...............
            t: 6                          _
                 ........................332..................
            t: 7                        _
                 ......................03.....................
            t: 8                     _
                 ....................3E.......................
            t: 9                   _
                 .................21..........................
            t: 10
                  ..............20E............................
            t: 11              _
                  ...........233...............................
            t: 12
                  ........2302.................................
            t: 13
                  .....23322...................................
            t: 14
                  ..233.22.....................................

        Example 4.14 from [Yamada2007]_::

            sage: B = SolitonCellularAutomata([[3],[1],[1],[1],[2],[3],[1],[1],[1],[2],[3],[3]], ['D',4,3])
            sage: B.print_states(15)
            t: 0
                 ....................................3...23...233
            t: 1
                 ...................................3..23..233...
            t: 2
                 ..................................3.23.233......
            t: 3
                 .................................323233.........
            t: 4
                 ................................0033............
            t: 5                                _
                 ..............................313...............
            t: 6
                 ...........................30E.3................
            t: 7                           _
                 ........................333...3.................
            t: 8
                 .....................3302....3..................
            t: 9
                 ..................33322.....3...................
            t: 10
                  ...............333.22......3....................
            t: 11
                  ............333..22.......3.....................
            t: 12
                  .........333...22........3......................
            t: 13
                  ......333....22.........3.......................
            t: 14
                  ...333.....22..........3........................
        """
        if num is None:
            num = len(self._states)
        if num > len(self._states):
            for _ in range(num - len(self._states)):
                self.evolve()

        vacuum = self._vacuum_elt
        num_factors = len(self._states[num-1])
        for i,state in enumerate(self._states[:num]):
            state = [vacuum]*(num_factors - len(state)) + list(state)
            output = [self._column_repr(b, vacuum_letter) for b in state]
            max_width = max(b.width() for b in output)
            start = ascii_art("t: %s \n" % i)
            start._baseline = -1
            print(start
                  + sum((ascii_art(' '*(max_width-b.width())) + b for b in output),
                        ascii_art('')))

    def latex_states(self, num=None, as_array=True, box_width='5pt'):
        r"""
        Return a latex version of the states.

        INPUT:

        - ``num`` -- the number of states
        - ``as_array`` -- boolean (default: ``True``); if ``True``, then the
          states are placed inside of an array; if ``False``, then the states
          are given as a word
        - ``box_width`` -- (default: ``'5pt'``) the width of the ``.`` used
          to represent the vacuum state when ``as_array`` is ``True``

        If ``as_array`` is ``False``, then the vacuum element is printed
        in a gray color. If ``as_array`` is ``True``, then the vacuum
        is given as ``.``

        Use the ``box_width`` to help create more even spacing when
        a column in the output contains only vacuum elements.

        EXAMPLES::

            sage: B = SolitonCellularAutomata('411122', 4)
            sage: B.latex_states(8)
            {\arraycolsep=0.5pt \begin{array}{c|ccccccccccccccccccc}
            t = 0 & \cdots & ... & \makebox[5pt]{.} & 4 & \makebox[5pt]{.}
             & \makebox[5pt]{.} & \makebox[5pt]{.} & 2 & 2 \\
            t = 1 & \cdots & ... & 4 & \makebox[5pt]{.} & \makebox[5pt]{.} & 2 & 2 & ... \\
            t = 2 & \cdots & ... & 4 & \makebox[5pt]{.} & 2 & 2 & ... \\
            t = 3 & \cdots & ... & 4 & 2 & 2 & ... \\
            t = 4 & \cdots & ... & 2 & 4 & 2 & ... \\
            t = 5 & \cdots & ... & 2 & 4 & \makebox[5pt]{.} & 2 & ... \\
            t = 6 & \cdots & ... & 2 & 4 & \makebox[5pt]{.} & \makebox[5pt]{.}
             & 2 & ... \\
            t = 7 & \cdots & \makebox[5pt]{.} & 2 & 4 & \makebox[5pt]{.}
             & \makebox[5pt]{.} & \makebox[5pt]{.} & 2 & ... \\
            \end{array}}

            sage: B = SolitonCellularAutomata('511122', 5)
            sage: B.latex_states(8, as_array=False)
            {\begin{array}{c|c}
            t = 0 & \cdots ... {\color{gray} 1} 5 {\color{gray} 1}
             {\color{gray} 1} {\color{gray} 1} 2 2 \\
            t = 1 & \cdots ... 5 {\color{gray} 1} {\color{gray} 1} 2 2 ... \\
            t = 2 & \cdots ... 5 {\color{gray} 1} 2 2 ... \\
            t = 3 & \cdots ... 5 2 2 ... \\
            t = 4 & \cdots ... 2 5 2 ... \\
            t = 5 & \cdots ... 2 5 {\color{gray} 1} 2 ... \\
            t = 6 & \cdots ... 2 5 {\color{gray} 1} {\color{gray} 1} 2 ... \\
            t = 7 & \cdots {\color{gray} 1} 2 5 {\color{gray} 1}
             {\color{gray} 1} {\color{gray} 1} 2 ... \\
            \end{array}}
        """
        from sage.misc.latex import latex, LatexExpr
        if not as_array:
            latex.add_package_to_preamble_if_available('xcolor')

        if num is None:
            num = len(self._states)
        if num > len(self._states):
            for _ in range(num - len(self._states)):
                self.evolve()

        vacuum = self._vacuum_elt

        def compact_repr(b):
            if as_array and b == vacuum:
                return "\\makebox[%s]{.}" % box_width

            if b.parent()._tableau_height == 1:
                temp = latex(b[0])
            else:
                temp = "\\begin{array}{@{}c@{}}"  # No padding around columns
                temp += r"\\".join(latex(letter) for letter in reversed(b))
                temp += "\\end{array}"

            if b == vacuum:
                return "{\\color{gray} %s}" % temp
            return temp # "\\makebox[%s]{$%s$}"%(box_width, temp)

        num_factors = len(self._states[num-1])
        if as_array:
            ret = "{\\arraycolsep=0.5pt \\begin{array}"
            ret += "{c|c%s}\n" % ('c'*num_factors)
        else:
            ret = "{\\begin{array}"
            ret += "{c|c}\n"
        for i,state in enumerate(self._states[:num]):
            state = [vacuum]*(num_factors-len(state)) + list(state)
            if as_array:
                ret += "t = %s & \\cdots & %s \\\\\n" % (i, r" & ".join(compact_repr(b) for b in state))
            else:
                ret += "t = %s & \\cdots %s \\\\\n" % (i, r" ".join(compact_repr(b) for b in state))
        ret += "\\end{array}}\n"
        return LatexExpr(ret)

    def print_state_evolution(self, num):
        r"""
        Print the evolution process of the state ``num``.

        .. SEEALSO::

            :meth:`state_evolution`, :meth:`latex_state_evolution`

        EXAMPLES::

            sage: B = SolitonCellularAutomata('1113123', 3)
            sage: B.evolve(3)
            sage: B.evolve(3)
            sage: B.print_state_evolution(0)
                  1         1         1         3         1         2         3
                  |         |         |         |         |         |         |
            111 --+-- 111 --+-- 111 --+-- 113 --+-- 112 --+-- 123 --+-- 113 --+-- 111
                  |         |         |         |         |         |         |
                  1         1         3         2         3         1         1
            sage: B.print_state_evolution(1)
                  1         1         3         2         3         1         1
                  |         |         |         |         |         |         |
            111 --+-- 113 --+-- 133 --+-- 123 --+-- 113 --+-- 111 --+-- 111 --+-- 111
                  |         |         |         |         |         |         |
                  3         3         2         1         1         1         1
        """
        u = self.state_evolution(num)  # Also evolves as necessary
        final = self._states[num+1]
        vacuum = self._vacuum_elt
        state = [vacuum]*(len(final) - len(self._states[num])) + list(self._states[num])
        carrier = KirillovReshetikhinTableaux(self._cartan_type, *self._evolutions[num])

        def simple_repr(x):
            return ''.join(repr(x).strip('[]').split(', '))

        def carrier_repr(x):
            if carrier._tableau_height == 1:
                return sum((ascii_art(repr(b)) if repr(b)[0] != '-'
                            else ascii_art("_" + '\n' + repr(b)[1:])
                            for b in x),
                           ascii_art(''))
            return ascii_art(''.join(repr(x).strip('[]').split(', ')))

        def cross_repr(i):
            ret = ascii_art(
"""
{!s:^7}
   |
 --+--
   |
{!s:^7}
""".format(simple_repr(state[i]), simple_repr(final[i])))
            ret._baseline = 2
            return ret
        art = sum((cross_repr(i)
                 + carrier_repr(u[i+1])
                 for i in range(len(state))), ascii_art(''))
        print(ascii_art(carrier_repr(u[0])) + art)

    def latex_state_evolution(self, num, scale=1):
        r"""
        Return a latex version of the evolution process of
        the state ``num``.

        .. SEEALSO::

            :meth:`state_evolution`, :meth:`print_state_evolution`

        EXAMPLES::

            sage: B = SolitonCellularAutomata('113123', 3)
            sage: B.evolve(3)
            sage: B.latex_state_evolution(0)
            \begin{tikzpicture}[scale=1]
            \node (i0) at (0.0,0.9) {$1$};
            \node (i1) at (2.48,0.9) {$1$};
            \node (i2) at (4.96,0.9) {$3$};
            ...
            \draw[->] (i5) -- (t5);
            \draw[->] (u6) -- (u5);
            \end{tikzpicture}
            sage: B.latex_state_evolution(1)
            \begin{tikzpicture}[scale=1]
            ...
            \end{tikzpicture}
        """
        from sage.graphs.graph_latex import setup_latex_preamble
        from sage.misc.latex import LatexExpr
        setup_latex_preamble()
        u = self.state_evolution(num) # Also evolves as necessary
        final = self._states[num+1]
        vacuum = self._vacuum_elt
        initial = [vacuum]*(len(final) - len(self._states[num])) + list(self._states[num])
        cs = len(u[0]) * 0.08 + 1 # carrier scaling

        def simple_repr(x):
            return ''.join(repr(x).strip('[]').split(', '))
        ret = '\\begin{{tikzpicture}}[scale={}]\n'.format(scale)
        for i,val in enumerate(initial):
            ret += '\\node (i{}) at ({},0.9) {{${}$}};\n'.format(i, 2*i*cs, simple_repr(val))
        for i,val in enumerate(final):
            ret += '\\node (t{}) at ({},-1) {{${}$}};\n'.format(i, 2*i*cs, simple_repr(val))
        for i,val in enumerate(u):
            ret += '\\node (u{}) at ({},0) {{${}$}};\n'.format(i, (2*i-1)*cs, simple_repr(val))
        for i in range(len(initial)):
            ret += '\\draw[->] (i{}) -- (t{});\n'.format(i, i)
            ret += '\\draw[->] (u{}) -- (u{});\n'.format(i+1, i)
        ret += '\\end{tikzpicture}'
        return LatexExpr(ret)


class PeriodicSolitonCellularAutomata(SolitonCellularAutomata):
    r"""
    A periodic soliton cellular automata.

    Fix some `r \in I_0`. A *periodic soliton cellular automata* is a
    :class:`SolitonCellularAutomata` with a state being a fixed number of
    tensor factors `p = p_{\ell} \otimes \cdots \otimes p_1 \otimes p_0`
    and the *time evolution* `T_s` is defined by

    .. MATH::

        R(p \otimes u) = u \otimes T_s(p),

    for some element `u \in B^{r,s}`.

    INPUT:

    - ``initial_state`` -- the list of elements, can also be a string
      when ``vacuum`` is 1 and ``n`` is `\mathfrak{sl}_n`
    - ``cartan_type`` -- (default: 2) the value ``n``, for `\mathfrak{sl}_n`,
      or a Cartan type
    - ``r`` -- (default: 1) the node index `r`; typically this
      corresponds to the height of the vacuum element

    EXAMPLES:

    The construction and usage is the same as for
    :class:`SolitonCellularAutomata`::

        sage: P = PeriodicSolitonCellularAutomata('1123334111241111423111411123112', 4)
        sage: P.evolve()
        sage: P
        Soliton cellular automata of type ['A', 3, 1] and vacuum = 1
          initial state:
          ..23334...24....423...4...23..2
          evoltuions: [(1, 31)]
          current state:
          34......24....243....4.223.233.
        sage: P.evolve(carrier_capacity=2)
        sage: P.evolve(carrier_index=2)
        sage: P.evolve(carrier_index=2, carrier_capacity=3)
        sage: P.print_states(10)
        t: 0
             ..23334...24....423...4...23..2
        t: 1
             34......24....243....4.223.233.
        t: 2
             ......24....24.3....4223.2333.4
        t: 3
             .....34....34.2..223234.24...3.
        t: 4
             ....34...23..242223.4..33....4.
        t: 5
             ..34.2223.224.3....4.33.....4..
        t: 6
             34223...24...3....433......4.22
        t: 7
             23....24....3...343....222434..
        t: 8
             ....24.....3..34.322244...3..23
        t: 9
             ..24.....332442342.......3.23..

    Using `r = 2` in type `A_3^{(1)}`::

        sage: initial = [[2,1],[2,1],[4,1],[2,1],[2,1],[2,1],[3,1],[3,1],[3,2]]
        sage: P = PeriodicSolitonCellularAutomata(initial, 4, 2)
        sage: P.print_states(10)
        t: 0   4   333
             ..1...112
        t: 1  4 333
             .1.112...
        t: 2 433     3
             112.....1
        t: 3 3    334
             2....111.
        t: 4   334   3
             ..111...2
        t: 5 34     33
             11.....21
        t: 6     3334
             ....1121.
        t: 7  333  4
             .112..1..
        t: 8 3    4 33
             2....1.11
        t: 9    3433
             ...1112..

    We do some examples in other types::

        sage: initial = [[1],[2],[2],[1],[1],[1],[3],[1],['E'],[1],[1]]
        sage: P = PeriodicSolitonCellularAutomata(initial, ['D',4,3])
        sage: P.print_states(10)
        t: 0
             .22...3.E..
        t: 1
             2....3.E..2
        t: 2
             ....3.E.22.
        t: 3
             ...3.E22...
        t: 4
             ..32E2.....
        t: 5
             .00.2......
        t: 6 _
             22.2.......
        t: 7
             2.2......3E
        t: 8
             .2.....30.2
        t: 9
             2....332.2.

        sage: P = PeriodicSolitonCellularAutomata([[3],[2],[1],[1],[-2]], ['C',2,1])
        sage: P.print_state_evolution(0)
                3           2           1           1          -2
            _   |           |           |       _   |      __   |       _
        11112 --+-- 11112 --+-- 11111 --+-- 11112 --+-- 11122 --+-- 11112
                |           |           |           |           |
                2           1          -2          -2           1

    REFERENCES:

    - [KTT2006]_
    - [KS2006]_
    - [YT2002]_
    - [YYT2003]_
    """
    def evolve(self, carrier_capacity=None, carrier_index=None, number=None):
        r"""
        Evolve ``self``.

        Time evolution `T_s` of a SCA state `p` is determined by

        .. MATH::

            u \otimes T_s(p) = R(p \otimes u),

        where `u` is some element in `B^{r,s}`.

        INPUT:

        - ``carrier_capacity`` -- (default: the number of balls in
          the system) the size `s` of carrier

        - ``carrier_index`` -- (default: the vacuum index) the index `r`
          of the carrier

        - ``number`` -- (optional) the number of times to perform
          the evolutions

        To perform multiple evolutions of the SCA, ``carrier_capacity``
        and ``carrier_index`` may be lists of the same length.

        .. WARNING::

            Time evolution is only guaranteed to result in a solution
            when the ``carrier_index`` is the defining `r` of the SCA.
            If no solution is found, then this will raise an error.

        EXAMPLES::

            sage: P = PeriodicSolitonCellularAutomata('12411133214131221122', 4)
            sage: P.evolve()
            sage: P.print_state(0)
            .24...332.4.3.22..22
            sage: P.print_state(1)
            4...33.2.42322..22..
            sage: P.evolve(carrier_capacity=2)
            sage: P.print_state(2)
            ..33.22.4232..22...4
            sage: P.evolve(carrier_capacity=[1,3,1,2])
            sage: P.evolve(1, number=3)
            sage: P.print_states(10)
            t: 0
                 .24...332.4.3.22..22
            t: 1
                 4...33.2.42322..22..
            t: 2
                 ..33.22.4232..22...4
            t: 3
                 .33.22.4232..22...4.
            t: 4
                 3222..43.2.22....4.3
            t: 5
                 222..43.2.22....4.33
            t: 6
                 2...4322.2.....43322
            t: 7
                 ...4322.2.....433222
            t: 8
                 ..4322.2.....433222.
            t: 9
                 .4322.2.....433222..

            sage: P = PeriodicSolitonCellularAutomata('12411132121', 4)
            sage: P.evolve(carrier_index=2, carrier_capacity=3)
            sage: P.state_evolution(0)
            [[[1, 1, 1], [2, 2, 4]],
             [[1, 1, 2], [2, 2, 4]],
             [[1, 1, 3], [2, 2, 4]],
             [[1, 1, 1], [2, 2, 3]],
             [[1, 1, 1], [2, 2, 3]],
             [[1, 1, 1], [2, 2, 3]],
             [[1, 1, 2], [2, 2, 3]],
             [[1, 1, 1], [2, 2, 2]],
             [[1, 1, 1], [2, 2, 2]],
             [[1, 1, 1], [2, 2, 2]],
             [[1, 1, 1], [2, 2, 4]],
             [[1, 1, 1], [2, 2, 4]]]
        """
        if isinstance(carrier_capacity, (list, tuple)):
            if not isinstance(carrier_index, (list, tuple)):
                carrier_index = [carrier_index] * len(carrier_capacity)
            if len(carrier_index) != len(carrier_capacity):
                raise ValueError("carrier_index and carrier_capacity"
                                 " must have the same length")
            for i, r in zip(carrier_capacity, carrier_index):
                self.evolve(i, r)
            return
        if isinstance(carrier_index, (list, tuple)):
            # carrier_capacity must be not be a list/tuple if given
            for r in carrier_index:
                self.evolve(carrier_capacity, r)
            return

        if carrier_capacity is None:
            carrier_capacity = self._nballs
        if carrier_index is None:
            carrier_index = self._vacuum

        if number is not None:
            for _ in range(number):
                self.evolve(carrier_capacity, carrier_index)
            return
        if carrier_capacity is None:
            carrier_capacity = self._nballs
        if carrier_index is None:
            carrier_index = self._vacuum

        K = KirillovReshetikhinTableaux(self._cartan_type, carrier_index, carrier_capacity)
        carrier_factor = (carrier_index, carrier_capacity)
        state = self._states[-1]
        dims = state.parent().dims
        for carrier in K:
            KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type,
                                                             dims + (carrier_factor,))
            elt = KRT(*(list(state) + [carrier]))
            RC = RiggedConfigurations(self._cartan_type, (carrier_factor,) + dims)
            elt2 = RC(*elt.to_rigged_configuration()).to_tensor_product_of_kirillov_reshetikhin_tableaux()
            # Back to an empty carrier or we are not getting any better
            if elt2[0] == carrier:
                KRT = TensorProductOfKirillovReshetikhinTableaux(self._cartan_type, dims)
                self._states.append(KRT(*elt2[1:]))
                self._evolutions.append(carrier_factor)
                self._initial_carrier.append(carrier)
                break
        else:
            raise ValueError("cannot find solution to time evolution")

    def __eq__(self, other):
        """
        Check equality.

        Two periodic SCAs are equal when they have the same initial
        state and evolutions.

        TESTS::

            sage: P1 = PeriodicSolitonCellularAutomata('34112223', 4)
            sage: P2 = PeriodicSolitonCellularAutomata('34112223', 4)
            sage: P1 == P2
            True
            sage: P1.evolve()
            sage: P1 == P2
            False
            sage: P2.evolve()
            sage: P1 == P2
            True
            sage: P1.evolve(5)
            sage: P2.evolve(6)
            sage: P1 == P2
            False

            sage: P = PeriodicSolitonCellularAutomata('34112223', 4)
            sage: B = SolitonCellularAutomata('34112223', 4)
            sage: P == B
            False
            sage: B == P
            False
        """
        return (isinstance(other, PeriodicSolitonCellularAutomata)
                and SolitonCellularAutomata.__eq__(self, other))
