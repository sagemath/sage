r"""
Cartan types

.. TODO::

   Why does sphinx complain if I use sections here?

Introduction

Loosely speaking, Dynkin diagrams (or equivalently Cartan matrices)
are graphs which are used to classify root systems, Coxeter and Weyl
groups, Lie algebras, Lie groups, crystals, etc. up to an
isomorphism. *Cartan types* are a standard set of names for those
Dynkin diagrams (see :wikipedia:`Dynkin_diagram`).

Let us consider, for example, the Cartan type `A_4`::

    sage: T = CartanType(['A', 4]); T
    ['A', 4]

It is the name of the following Dynkin diagram::

    sage: DynkinDiagram(T)                                                              # needs sage.graphs
    O---O---O---O
    1   2   3   4
    A4

.. NOTE::

    For convenience, the following shortcuts are available::

        sage: DynkinDiagram(['A',4])                                                    # needs sage.graphs
        O---O---O---O
        1   2   3   4
        A4
        sage: DynkinDiagram('A4')                                                       # needs sage.graphs
        O---O---O---O
        1   2   3   4
        A4
        sage: T.dynkin_diagram()                                                        # needs sage.graphs
        O---O---O---O
        1   2   3   4
        A4

See :class:`~sage.combinat.root_system.dynkin_diagram.DynkinDiagram`
for how to further manipulate Dynkin diagrams.

From this data (the *Cartan datum*), one can construct the associated
root system::

    sage: RootSystem(T)
    Root system of type ['A', 4]

The associated Weyl group of `A_n` is the symmetric group `S_{n+1}`::

    sage: W = WeylGroup(T); W                                                           # needs sage.libs.gap
    Weyl Group of type ['A', 4] (as a matrix group acting on the ambient space)
    sage: W.cardinality()                                                               # needs sage.libs.gap
    120

while the Lie algebra is `sl_{n+1}`, and the Lie group `SL_{n+1}`
(TODO: illustrate this once this is implemented).

One may also construct crystals associated to various Dynkin diagrams.
For example::

    sage: C = crystals.Letters(T); C                                                    # needs sage.combinat
    The crystal of letters for type ['A', 4]
    sage: C.list()                                                                      # needs sage.combinat
    [1, 2, 3, 4, 5]

    sage: C = crystals.Tableaux(T, shape=[2]); C                                        # needs sage.combinat
    The crystal of tableaux of type ['A', 4] and shape(s) [[2]]
    sage: C.cardinality()                                                               # needs sage.combinat
    15

Here is a sample of all the finite irreducible crystallographic Cartan
types::

    sage: CartanType.samples(finite=True, crystallographic=True)
    [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
     ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2]]

One can also get latex representations of the crystallographic Cartan
types and their corresponding Dynkin diagrams::

    sage: [latex(ct) for ct in CartanType.samples(crystallographic=True)]
    [A_{1}, A_{5}, B_{1}, B_{5}, C_{1}, C_{5}, D_{2}, D_{3}, D_{5},
     E_6, E_7, E_8, F_4, G_2,
     A_{1}^{(1)}, A_{5}^{(1)}, B_{1}^{(1)}, B_{5}^{(1)},
     C_{1}^{(1)}, C_{5}^{(1)}, D_{3}^{(1)}, D_{5}^{(1)},
     E_6^{(1)}, E_7^{(1)}, E_8^{(1)}, F_4^{(1)}, G_2^{(1)},
     BC_{1}^{(2)}, BC_{5}^{(2)},
     B_{5}^{(1)\vee}, C_{4}^{(1)\vee}, F_4^{(1)\vee},
     G_2^{(1)\vee}, BC_{1}^{(2)\vee}, BC_{5}^{(2)\vee}]
    sage: view([DynkinDiagram(ct) for ct in CartanType.samples(crystallographic=True)]) # not tested

Non-crystallographic Cartan types are also partially supported::

    sage: CartanType.samples(finite=True, crystallographic=False)
    [['I', 5], ['H', 3], ['H', 4]]

In Sage, a Cartan type is used as a database of type-specific
information and algorithms (see e.g. :mod:`sage.combinat.root_system.type_A`).
This database includes how to construct the Dynkin diagram, the ambient space
for the root system (see :wikipedia:`Root_system`), and further
mathematical properties::

    sage: T.is_finite(), T.is_simply_laced(), T.is_affine(), T.is_crystallographic()
    (True, True, False, True)

In particular, a Sage Cartan type is endowed with a fixed choice of
labels for the nodes of the Dynkin diagram. This choice follows the
conventions of Nicolas Bourbaki, Lie Groups and Lie Algebras: Chapter 4-6,
Elements of Mathematics, Springer (2002). ISBN 978-3540426509. For example::

    sage: T = CartanType(['D', 4])
    sage: DynkinDiagram(T)                                                              # needs sage.graphs
        O 4
        |
        |
    O---O---O
    1   2   3
    D4

    sage: E6 = CartanType(['E',6])
    sage: DynkinDiagram(E6)                                                             # needs sage.graphs
            O 2
            |
            |
    O---O---O---O---O
    1   3   4   5   6
    E6

.. NOTE::

    The direction of the arrows is the **opposite** (i.e. the transpose)
    of Bourbaki's convention, but agrees with Kac's.

    For example, in type `C_2`, we have::

        sage: C2 = DynkinDiagram(['C',2]); C2                                           # needs sage.graphs
        O=<=O
        1   2
        C2
        sage: C2.cartan_matrix()                                                        # needs sage.graphs
        [ 2 -2]
        [-1  2]

    However Bourbaki would have the Cartan matrix as:

    .. MATH::

        \begin{bmatrix}
        2 & -1 \\
        -2 & 2
        \end{bmatrix}.

If desired, other node labelling conventions can be achieved. For
example the Kac labelling for type `E_6` can be obtained via::

    sage: E6.relabel({1:1,2:6,3:2,4:3,5:4,6:5}).dynkin_diagram()                        # needs sage.graphs
            O 6
            |
            |
    O---O---O---O---O
    1   2   3   4   5
    E6 relabelled by {1: 1, 2: 6, 3: 2, 4: 3, 5: 4, 6: 5}

Contributions implementing other conventions are very welcome.

Another option is to build from scratch a new Dynkin diagram.  The
architecture has been designed to make it fairly easy to add other
labelling conventions. In particular, we strived at choosing type free
algorithms whenever possible, so in principle most features should
remain available even with custom Cartan types. This has not been used
much yet, so some rough corners certainly remain.

Here, we construct the hyperbolic example of Exercise 4.9 p. 57 of
Kac, Infinite Dimensional Lie Algebras. We start with an empty Dynkin
diagram, and add a couple nodes::

    sage: g = DynkinDiagram()                                                           # needs sage.graphs
    sage: g.add_vertices([1,2,3])                                                       # needs sage.graphs

Note that the diagonal of the Cartan matrix is already initialized::

    sage: g.cartan_matrix()                                                             # needs sage.graphs
    [2 0 0]
    [0 2 0]
    [0 0 2]

Then we add a couple edges::

    sage: g.add_edge(1,2,2)                                                             # needs sage.graphs
    sage: g.add_edge(1,3)                                                               # needs sage.graphs
    sage: g.add_edge(2,3)                                                               # needs sage.graphs

and we get the desired Cartan matrix::

    sage: g.cartan_matrix()                                                             # needs sage.graphs
    [2 0 0]
    [0 2 0]
    [0 0 2]

Oops, the Cartan matrix did not change! This is because it is cached
for efficiency (see :class:`cached_method`). In general, a Dynkin
diagram should not be modified after having been used.

.. WARNING:: this is not checked currently

.. TODO:: add a method :meth:`set_mutable` as, say, for matrices

Here, we can work around this by clearing the cache::

    sage: delattr(g, 'cartan_matrix')                                                   # needs sage.graphs

Now we get the desired Cartan matrix::

    sage: g.cartan_matrix()                                                             # needs sage.graphs
    [ 2 -1 -1]
    [-2  2 -1]
    [-1 -1  2]

Note that backward edges have been automatically added::

    sage: g.edges(sort=True)                                                            # needs sage.graphs
    [(1, 2, 2), (1, 3, 1), (2, 1, 1), (2, 3, 1), (3, 1, 1), (3, 2, 1)]

.. rubric:: Reducible Cartan types

Reducible Cartan types can be specified by passing a sequence
or list of irreducible Cartan types::

    sage: CartanType(['A',2],['B',2])
    A2xB2
    sage: CartanType([['A',2],['B',2]])
    A2xB2
    sage: CartanType(['A',2],['B',2]).is_reducible()
    True

or using the following short hand notation::

    sage: CartanType("A2xB2")
    A2xB2
    sage: CartanType("A2","B2") == CartanType("A2xB2")
    True

.. rubric:: Degenerate cases

When possible, type `I_n` is automatically converted to the isomorphic
crystallographic Cartan types (any reason not to do so?)::

    sage: CartanType(["I",1])
    A1xA1
    sage: CartanType(["I",3])
    ['A', 2]
    sage: CartanType(["I",4])
    ['C', 2]
    sage: CartanType(["I",6])
    ['G', 2]

The Dynkin diagrams for types `B_1`, `C_1`, `D_2`, and `D_3` are
isomorphic to that for `A_1`, `A_1`, `A_1 \times A_1`, and `A_3`,
respectively. However their natural ambient space realizations (stemming
from the corresponding infinite families of Lie groups) are different.
Therefore, the Cartan types are considered as distinct::

    sage: CartanType(['B',1])
    ['B', 1]
    sage: CartanType(['C',1])
    ['C', 1]
    sage: CartanType(['D',2])
    ['D', 2]
    sage: CartanType(['D',3])
    ['D', 3]

.. rubric:: Affine Cartan types

For affine types, we use the usual conventions for affine Coxeter
groups: each affine type is either untwisted (that is arise from the
natural affinisation of a finite Cartan type)::

    sage: CartanType(["A", 4, 1]).dynkin_diagram()                                      # needs sage.graphs
    0
    O-----------+
    |           |
    |           |
    O---O---O---O
    1   2   3   4
    A4~
    sage: CartanType(["B", 4, 1]).dynkin_diagram()                                      # needs sage.graphs
        O 0
        |
        |
    O---O---O=>=O
    1   2   3   4
    B4~

or dual thereof::

    sage: CartanType(["B", 4, 1]).dual().dynkin_diagram()                               # needs sage.graphs
        O 0
        |
        |
    O---O---O=<=O
    1   2   3   4
    B4~*

or is of type `\widetilde{BC}_n` (which yields an irreducible, but
nonreduced root system)::

    sage: CartanType(["BC", 4, 2]).dynkin_diagram()                                     # needs sage.graphs
    O=<=O---O---O=<=O
    0   1   2   3   4
    BC4~

This includes the two degenerate cases::

    sage: CartanType(["A", 1, 1]).dynkin_diagram()                                      # needs sage.graphs
    O<=>O
    0   1
    A1~
    sage: CartanType(["BC", 1, 2]).dynkin_diagram()                                     # needs sage.graphs
      4
    O=<=O
    0   1
    BC1~

For the user convenience, Kac's notations for twisted affine types are
automatically translated into the previous ones::

    sage: # needs sage.graphs
    sage: CartanType(["A", 9, 2])
    ['B', 5, 1]^*
    sage: CartanType(["A", 9, 2]).dynkin_diagram()
        O 0
        |
        |
    O---O---O---O=<=O
    1   2   3   4   5
    B5~*
    sage: CartanType(["A", 10, 2]).dynkin_diagram()
    O=<=O---O---O---O=<=O
    0   1   2   3   4   5
    BC5~
    sage: CartanType(["D", 5, 2]).dynkin_diagram()
    O=<=O---O---O=>=O
    0   1   2   3   4
    C4~*
    sage: CartanType(["D", 4, 3]).dynkin_diagram()
      3
    O=>=O---O
    2   1   0
    G2~* relabelled by {0: 0, 1: 2, 2: 1}
    sage: CartanType(["E", 6, 2]).dynkin_diagram()
    O---O---O=<=O---O
    0   1   2   3   4
    F4~*

Additionally one can set the notation option to use Kac's notation::

    sage: # needs sage.graphs
    sage: CartanType.options['notation'] = 'Kac'
    sage: CartanType(["A", 9, 2])
    ['A', 9, 2]
    sage: CartanType(["A", 9, 2]).dynkin_diagram()
        O 0
        |
        |
    O---O---O---O=<=O
    1   2   3   4   5
    A9^2
    sage: CartanType(["A", 10, 2]).dynkin_diagram()
    O=<=O---O---O---O=<=O
    0   1   2   3   4   5
    A10^2
    sage: CartanType(["D", 5, 2]).dynkin_diagram()
    O=<=O---O---O=>=O
    0   1   2   3   4
    D5^2
    sage: CartanType(["D", 4, 3]).dynkin_diagram()
      3
    O=>=O---O
    2   1   0
    D4^3
    sage: CartanType(["E", 6, 2]).dynkin_diagram()
    O---O---O=<=O---O
    0   1   2   3   4
    E6^2
    sage: CartanType.options['notation'] = 'BC'

.. rubric:: Infinite Cartan types

There are minimal implementations of the Cartan types `A_{\infty}` and
`A_{+\infty}`.  In sage `oo` is the same as `+Infinity`, so `NN` and `ZZ` are
used to differentiate between the `A_{+\infty}` and `A_{\infty}` root systems::

    sage: CartanType(['A', NN])
    ['A', NN]
    sage: print(CartanType(['A', NN]).ascii_art())
    O---O---O---O---O---O---O---..
    0   1   2   3   4   5   6
    sage: CartanType(['A', ZZ])
    ['A', ZZ]
    sage: print(CartanType(['A', ZZ]).ascii_art())
    ..---O---O---O---O---O---O---O---..
        -3  -2  -1   0   1   2   3

There are also the following shorthands::

    sage: CartanType("Aoo")
    ['A', ZZ]
    sage: CartanType("A+oo")
    ['A', NN]

.. rubric:: Abstract classes for Cartan types

- :class:`CartanType_abstract`
- :class:`CartanType_crystallographic`
- :class:`CartanType_simply_laced`
- :class:`CartanType_simple`
- :class:`CartanType_finite`
- :class:`CartanType_affine` (see also :ref:`sage.combinat.root_system.type_affine`)
- :obj:`sage.combinat.root_system.cartan_type.CartanType`
- :ref:`sage.combinat.root_system.type_dual`
- :ref:`sage.combinat.root_system.type_reducible`
- :ref:`sage.combinat.root_system.type_relabel`

Concrete classes for Cartan types

- :class:`CartanType_standard`
- :class:`CartanType_standard_finite`
- :class:`CartanType_standard_affine`
- :class:`CartanType_standard_untwisted_affine`

Type specific data

The data essentially consists of a description of the Dynkin/Coxeter
diagram and, when relevant, of the natural embedding of the root
system in a Euclidean space. Everything else is reconstructed from
this data.

- :ref:`sage.combinat.root_system.type_A`
- :ref:`sage.combinat.root_system.type_B`
- :ref:`sage.combinat.root_system.type_C`
- :ref:`sage.combinat.root_system.type_D`
- :ref:`sage.combinat.root_system.type_E`
- :ref:`sage.combinat.root_system.type_F`
- :ref:`sage.combinat.root_system.type_G`
- :ref:`sage.combinat.root_system.type_H`
- :ref:`sage.combinat.root_system.type_I`
- :ref:`sage.combinat.root_system.type_super_A`
- :ref:`sage.combinat.root_system.type_Q`
- :ref:`sage.combinat.root_system.type_A_affine`
- :ref:`sage.combinat.root_system.type_B_affine`
- :ref:`sage.combinat.root_system.type_C_affine`
- :ref:`sage.combinat.root_system.type_D_affine`
- :ref:`sage.combinat.root_system.type_E_affine`
- :ref:`sage.combinat.root_system.type_F_affine`
- :ref:`sage.combinat.root_system.type_G_affine`
- :ref:`sage.combinat.root_system.type_BC_affine`
- :ref:`sage.combinat.root_system.type_A_infinity`

.. TODO:: Should those indexes come before the introduction?
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_import import LazyImport
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.global_options import GlobalOptions
from sage.sets.family import Family

# TODO:
# Implement the Kac conventions by relabeling/dual/... of the above
# Implement Coxeter diagrams for non crystallographic


# Intention: we want simultaneously CartanType to be a factory for
# the various subtypes of CartanType_abstract, as in:
#     CartanType(["A",4,1])
# and to behaves as a "module" for some extra utilities:
#     CartanType.samples()
#
# Implementation: CartanType is the unique instance of this class
# CartanTypeFactory. Is there a better/more standard way to do it?

class CartanTypeFactory(SageObject):

    def __call__(self, *args):
        """
        Construct a Cartan type object.

        INPUT:

        - ``[letter, rank]`` -- letter is one of 'A', 'B', 'C', 'D', 'E', 'F', 'G'
          and rank is an integer or a pair of integers

        - ``[letter, rank, twist]`` -- letter is one of 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'BC'
          and rank and twist are integers

        - ``str`` -- string

        - ``object`` -- a Cartan type, or an object with a Cartan type method

        EXAMPLES:

        We construct the Cartan type `D_4`::

            sage: d4 = CartanType(['D',4])
            sage: d4
            ['D', 4]

        or, for short::

            sage: CartanType("D4")
            ['D', 4]

        .. SEEALSO:: :func:`~sage.combinat.root_system.cartan_type.CartanType`

        TESTS:

        Check that this is compatible with :class:`CartanTypeFolded`::

            sage: fct = CartanType(['C', 4, 1]).as_folding()
            sage: CartanType(fct)
            ['C', 4, 1]

        Check that :issue:`13774` is fixed::

            sage: CT = CartanType([['A',2]])
            sage: CT.is_irreducible()
            True
            sage: CT.cartan_matrix()                                                    # needs sage.graphs
            [ 2 -1]
            [-1  2]
            sage: CT = CartanType(['A2'])
            sage: CT.is_irreducible()
            True
            sage: CartanType('A2')
            ['A', 2]

        Check that we can pass any Cartan type as a single element list::

            sage: CT = CartanType(['A2', 'A2', 'A2'])
            sage: CartanType([CT])
            A2xA2xA2

            sage: CT = CartanType('A2').relabel({1:-1, 2:-2})
            sage: CartanType([CT])
            ['A', 2] relabelled by {1: -1, 2: -2}

        Check the errors from :issue:`20973`::

            sage: CartanType(['A',-1])
            Traceback (most recent call last):
            ...
            ValueError: ['A', -1] is not a valid Cartan type

        Check that unicode is handled properly (:issue:`23323`)::

            sage: CartanType(u"A3")
            ['A', 3]
        """
        if len(args) == 1:
            t = args[0]
        else:
            t = args

        if isinstance(t, (CartanType_abstract, SuperCartanType_standard)):
            return t
        if hasattr(t, "cartan_type"):
            return t.cartan_type()

        if len(t) == 1: # Fix for trac #13774
            t = t[0]

        # We need to make another check
        if isinstance(t, (CartanType_abstract, SuperCartanType_standard)):
            return t

        from sage.rings.semirings.non_negative_integer_semiring import NN
        if isinstance(t, str):
            if "x" in t:
                from . import type_reducible
                return type_reducible.CartanType([CartanType(u) for u in t.split("x")])
            elif t[-1] == "*":
                return CartanType(t[:-1]).dual()
            elif t[-1] == "~":
                return CartanType(t[:-1]).affine()
            elif t in ["Aoo", "A∞"]:
                return CartanType(['A', Infinity])
            elif t == "A+oo":
                from . import type_A_infinity
                return type_A_infinity.CartanType(NN)
            else:
                return CartanType([t[0], eval(t[1:])])

        t = list(t)
        if isinstance(t[0], str) and t[1] in [Infinity, ZZ, NN]:
            letter, n = t[0], t[1]
            if letter == 'A':
                from . import type_A_infinity
                if t[1] == NN:
                    return type_A_infinity.CartanType(NN)
                else:
                    return type_A_infinity.CartanType(ZZ)

        if isinstance(t[0], str) and t[1] in ZZ and t[1] >= 0:
            letter, n = t[0], t[1]
            if len(t) == 2:
                if letter == "A":
                    if n >= 0:
                        from . import type_A
                        return type_A.CartanType(n)
                if letter == "B":
                    if n >= 1:
                        from . import type_B
                        return type_B.CartanType(n)
                if letter == "C":
                    if n >= 1:
                        from . import type_C
                        return type_C.CartanType(n)
                if letter == "D":
                    from . import type_D
                    if n >= 2:
                        return type_D.CartanType(n)
                if letter == "E":
                    if n >= 6 and n <= 8:
                        from . import type_E
                        return type_E.CartanType(n)
                if letter == "F":
                    if n == 4:
                        from . import type_F
                        return type_F.CartanType()
                if letter == "G":
                    if n == 2:
                        from . import type_G
                        return type_G.CartanType()
                if letter == "H":
                    if n in [3, 4]:
                        from . import type_H
                        return type_H.CartanType(n)
                if letter == "I":
                    if n == 1:
                        return CartanType([["A", 1], ["A", 1]])
                    if n == 3:
                        return CartanType(["A", 2])
                    if n == 4:
                        return CartanType(["C", 2])
                    if n == 6:
                        return CartanType(["G", 2])
                    if n >= 1:
                        from . import type_I
                        return type_I.CartanType(n)
                if letter == "Q":
                    if n >= 1:
                        from . import type_Q
                        return type_Q.CartanType(n)

            if len(t) == 3:
                if t[2] == 1: # Untwisted affine
                    if letter == "A":
                        if n >= 1:
                            from . import type_A_affine
                            return type_A_affine.CartanType(n)
                    if letter == "B":
                        if n >= 1:
                            from . import type_B_affine
                            return type_B_affine.CartanType(n)
                    if letter == "C":
                        if n >= 1:
                            from . import type_C_affine
                            return type_C_affine.CartanType(n)
                    if letter == "D":
                        from . import type_D_affine
                        if n >= 3:
                            return type_D_affine.CartanType(n)
                    if letter == "E":
                        if n >= 6 and n <= 8:
                            from . import type_E_affine
                            return type_E_affine.CartanType(n)
                    if letter == "F":
                        if n == 4:
                            from . import type_F_affine
                            return type_F_affine.CartanType()
                    if letter == "G":
                        if n == 2:
                            from . import type_G_affine
                            return type_G_affine.CartanType()
                if t[2] in [2,3]:
                    if letter == "BC" and t[2] == 2:
                        if n >= 1:
                            from . import type_BC_affine
                            return type_BC_affine.CartanType(n)
                    if letter == "A" and t[2] == 2:
                        if n % 2 == 0: # Kac' A_2n^(2)
                            return CartanType(["BC", ZZ(n//2), 2])
                        else:        # Kac' A_2n-1^(2)
                            return CartanType(["B", ZZ((n+1)//2), 1]).dual()
                    if letter == "D" and t[2] == 2:
                        return CartanType(["C", n-1, 1]).dual()
                    if letter == "D" and t[2] == 3 and n == 4:
                        return CartanType(["G", 2, 1]).dual().relabel([0,2,1])
                    if letter == "E" and t[2] == 2 and n == 6:
                        return CartanType(["F", 4, 1]).dual()
            raise ValueError("%s is not a valid Cartan type" % t)

        if isinstance(t[0], str) and isinstance(t[1], (list, tuple)):
            letter, n = t[0], t[1]
            if len(t) == 2 and len(n) == 2:
                from . import type_super_A
                return type_super_A.CartanType(n[0], n[1])
            raise ValueError("%s is not a valid super Cartan type" % t)

        # As the Cartan type has not been recognised try subtypes - but check
        # for the error noted in trac:???
        from . import type_reducible
        try:
            return type_reducible.CartanType([ CartanType(subtype) for subtype in t ])
        except (SyntaxError, ValueError):
            raise ValueError("%s is not a valid Cartan type" % t)

    def _repr_(self):
        """
        EXAMPLES::

            sage: CartanType    # indirect doctest
            CartanType
        """
        return "CartanType"

    def samples(self, finite=None, affine=None, crystallographic=None):
        """
        Return a sample of the available Cartan types.

        INPUT:

        - ``finite`` -- boolean or ``None`` (default: ``None``)

        - ``affine`` -- boolean or ``None`` (default: ``None``)

        - ``crystallographic`` -- boolean or ``None`` (default: ``None``)

        The sample contains all the exceptional finite and affine
        Cartan types, as well as typical representatives of the
        infinite families.

        EXAMPLES::

            sage: CartanType.samples()
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4],
             ['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['BC', 1, 2], ['BC', 5, 2],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*, ['BC', 1, 2]^*, ['BC', 5, 2]^*]

        The finite, affine and crystallographic options allow
        respectively for restricting to (non) finite, (non) affine,
        and (non) crystallographic Cartan types::

            sage: CartanType.samples(finite=True)
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4]]

            sage: CartanType.samples(affine=True)
            [['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['BC', 1, 2], ['BC', 5, 2],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*, ['BC', 1, 2]^*, ['BC', 5, 2]^*]

            sage: CartanType.samples(crystallographic=True)
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2],
             ['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['BC', 1, 2], ['BC', 5, 2],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*, ['BC', 1, 2]^*, ['BC', 5, 2]^*]

            sage: CartanType.samples(crystallographic=False)
            [['I', 5], ['H', 3], ['H', 4]]

        .. TODO:: add some reducible Cartan types (suggestions?)

        TESTS::

            sage: for ct in CartanType.samples(): TestSuite(ct).run()
        """
        result = self._samples()
        if crystallographic is not None:
            result = [t for t in result if t.is_crystallographic() == crystallographic ]
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        return result

    @cached_method
    def _samples(self):
        """
        Return a sample of all implemented Cartan types.

        .. NOTE:: This is intended to be used through :meth:`samples`.

        EXAMPLES::

            sage: CartanType._samples()
            [['A', 1], ['A', 5], ['B', 1], ['B', 5], ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
             ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4],
             ['A', 1, 1], ['A', 5, 1], ['B', 1, 1], ['B', 5, 1],
             ['C', 1, 1], ['C', 5, 1], ['D', 3, 1], ['D', 5, 1],
             ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['BC', 1, 2], ['BC', 5, 2],
             ['B', 5, 1]^*, ['C', 4, 1]^*, ['F', 4, 1]^*, ['G', 2, 1]^*, ['BC', 1, 2]^*, ['BC', 5, 2]^*]
        """
        finite_crystallographic = [CartanType(t)
            for t in [['A', 1], ['A', 5], ['B', 1], ['B', 5],
                      ['C', 1], ['C', 5], ['D', 2], ['D', 3], ['D', 5],
                      ["E", 6], ["E", 7], ["E", 8],
                      ["F", 4],
                      ["G", 2]]]

        # Support for hand constructed Dynkin diagrams as Cartan types is not yet ready enough for including an example here.
        # from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
        # g = DynkinDiagram_class.an_instance()
        return finite_crystallographic + \
            [CartanType(t) for t in [["I", 5], ["H", 3], ["H", 4]]] + \
            [t.affine() for t in finite_crystallographic if t.is_irreducible()] + \
            [CartanType(t) for t in [["BC", 1, 2], ["BC", 5, 2]]] + \
            [CartanType(t).dual() for t in [["B", 5, 1], ["C", 4, 1],
                                            ["F", 4, 1], ["G", 2, 1],
                                            ["BC", 1, 2], ["BC", 5, 2]]]  # + \
            # [ g ]

    _colors = {1: 'blue', -1: 'blue',
               2: 'red', -2: 'red',
               3: 'green', -3: 'green',
               4: 'cyan', -4: 'cyan',
               5: 'magenta', -5: 'magenta',
               6: 'yellow', -6: 'yellow'}

    @classmethod
    def color(cls, i):
        """
        Default color scheme for the vertices of a Dynkin diagram (and associated objects).

        EXAMPLES::

            sage: CartanType.color(1)
            'blue'
            sage: CartanType.color(2)
            'red'
            sage: CartanType.color(3)
            'green'

        The default color is black::

            sage: CartanType.color(0)
            'black'

        Negative indices get the same color as their positive counterparts::

            sage: CartanType.color(-1)
            'blue'
            sage: CartanType.color(-2)
            'red'
            sage: CartanType.color(-3)
            'green'
        """
        return cls._colors.get(i, 'black')

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for Cartan types. If no parameters
        are set, then the function returns a copy of the options dictionary.

        The ``options`` to partitions can be accessed as the method
        :obj:`CartanType.options` of
        :class:`CartanType <CartanTypeFactory>`.

        @OPTIONS@

        EXAMPLES::

            sage: ct = CartanType(['D',5,2]); ct
            ['C', 4, 1]^*
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            O=<=O---O---O=>=O
            0   1   2   3   4
            C4~*
            sage: latex(ct)
            C_{4}^{(1)\vee}
            sage: CartanType.options(dual_str='#', dual_latex='\\ast',)
            sage: ct
            ['C', 4, 1]^#
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            O=<=O---O---O=>=O
            0   1   2   3   4
            C4~#
            sage: latex(ct)
            C_{4}^{(1)\ast}
            sage: CartanType.options(notation='kac', mark_special_node='both')
            sage: ct
            ['D', 5, 2]
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            @=<=O---O---O=>=O
            0   1   2   3   4
            D5^2
            sage: latex(ct)
            D_{5}^{(2)}

        For type `A_{2n}^{(2)\dagger}`, the dual string/latex options are
        automatically overridden::

            sage: dct = CartanType(['A',8,2]).dual(); dct
            ['A', 8, 2]^+
            sage: latex(dct)
            A_{8}^{(2)\dagger}
            sage: dct.dynkin_diagram()                                                  # needs sage.graphs
            @=>=O---O---O=>=O
            0   1   2   3   4
            A8^2+
            sage: CartanType.options._reset()
        """
        NAME = 'CartanType'
        module = 'sage.combinat.root_system.cartan_type'
        option_class = 'CartanTypeFactory'
        notation = dict(default='Stembridge',
                      description='Specifies which notation Cartan types should use when printed',
                      values=dict(Stembridge="use Stembridge's notation",
                                  Kac="use Kac's notation"),
                      case_sensitive=False,
                      alias=dict(BC='Stembridge', tilde='Stembridge', twisted='Kac'))
        dual_str = dict(default='*',
                      description='The string used for dual Cartan types when printing',
                      checker=lambda char: isinstance(char, str))
        dual_latex = dict(default='\\vee',
                        description='The latex used for dual CartanTypes when latexing',
                        checker=lambda char: isinstance(char, str))
        mark_special_node = dict(default='none',
                               description="Make the special nodes",
                               values=dict(none="no markup", latex="only in latex",
                                           printing="only in printing", both="both in latex and printing"),
                               case_sensitive=False)
        special_node_str = dict(default='@',
                              description="The string used to indicate which node is special when printing",
                              checker=lambda char: isinstance(char, str))
        marked_node_str = dict(default='X',
                             description="The string used to indicate a marked node when printing",
                             checker=lambda char: isinstance(char, str))
        latex_relabel = dict(default=True,
                           description="Indicate in the latex output if a Cartan type has been relabelled",
                           checker=lambda x: isinstance(x, bool))
        latex_marked = dict(default=True,
                          description="Indicate in the latex output if a Cartan type has been marked",
                          checker=lambda x: isinstance(x, bool))


CartanType = CartanTypeFactory()
CartanType.__doc__ = __doc__


class CartanType_abstract:
    r"""
    Abstract class for Cartan types

    Subclasses should implement:

    - :meth:`dynkin_diagram()`

    - :meth:`cartan_matrix()`

    - :meth:`is_finite()`

    - :meth:`is_affine()`

    - :meth:`is_irreducible()`
    """

    def type(self):
        r"""
        Return the type of ``self``, or ``None`` if unknown.

        This method should be overridden in any subclass.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.type() is None
            True
        """
        return None

    def _add_abstract_superclass(self, classes):
        """
        Add abstract super-classes to the class of ``self``.

        INPUT:

        - ``classes`` -- an abstract class or tuple thereof

        EXAMPLES::

            sage: C = CartanType(["A",3,1])
            sage: class MyCartanType:
            ....:     def my_method(self):
            ....:         return 'I am here!'
            sage: C._add_abstract_superclass(MyCartanType)
            sage: C.__class__
            <class 'sage.combinat.root_system.type_A_affine.CartanType_with_superclass_with_superclass'>
            sage: C.__class__.__bases__
            (<class 'sage.combinat.root_system.type_A_affine.CartanType_with_superclass'>,
             <class ...__main__.MyCartanType...>)
            sage: C.my_method()
            'I am here!'

        .. TODO:: Generalize to :class:`SageObject`?
        """
        from sage.structure.dynamic_class import dynamic_class
        assert isinstance(classes, (tuple, type))
        if not isinstance(classes, tuple):
            classes = (classes,)
        bases = (self.__class__,) + classes
        self.__class__ = dynamic_class(self.__class__.__name__+"_with_superclass", bases)

    def _ascii_art_node(self, label):
        """
        Return the ascii art for the node labelled by ``label``.

        EXAMPLES::

            sage: CartanType(['A',3])._ascii_art_node(2)
            'O'
        """
        return "O"

    def _latex_draw_node(self, x, y, label, position='below=4pt', fill='white'):
        r"""
        Draw (possibly marked [crossed out]) circular node ``i`` at the
        position ``(x,y)`` with node label ``label`` .

        - ``position`` -- position of the label relative to the node
        - ``anchor`` -- (optional) the anchor point for the label

        EXAMPLES::

            sage: CartanType(['A',3])._latex_draw_node(0, 0, 1)
            '\\draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};\n'
        """
        return "\\draw[fill={}] ({} cm, {} cm) circle (.25cm) node[{}]{{${}$}};\n".format(
                fill, x, y, position, label)

    def _latex_draw_arrow_tip(self, x, y, rot=0):
        r"""
        Draw an arrow tip at the point ``(x, y)`` rotated by ``rot``.

        INPUT:

        - ``(x, y)`` -- the coordinates of a point, in cm

        - ``rot`` -- an angle, in degrees

        This is an internal function used to assist drawing the Dynkin
        diagrams. See e.g. :meth:`~sage.combinat.root_system.type_B.CartanType._latex_dynkin_diagram`.

        EXAMPLES::

            sage: CartanType(['B',2])._latex_draw_arrow_tip(1, 0, 180)
            '\\draw[shift={(1, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);\n'
        """
        return "\\draw[shift={(%s, %s)}, rotate=%s] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);\n" % (x, y, rot)

    @abstract_method
    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter or
        Dynkin diagram.

        EXAMPLES::

            sage: CartanType(['A', 4]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            5
            sage: CartanType(['I', 8]).rank()
            2
        """
        #return len(self.index_set())

    @abstract_method
    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter or
        Dynkin diagram.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).index_set()
            (0, 1, 2, 3)
            sage: CartanType(['D', 4]).index_set()
            (1, 2, 3, 4)
            sage: CartanType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 6, 2]).index_set()
            (0, 1, 2, 3)
            sage: CartanType(['D', 6, 2]).index_set()
            (0, 1, 2, 3, 4, 5)
            sage: CartanType(['E', 6, 1]).index_set()
            (0, 1, 2, 3, 4, 5, 6)
            sage: CartanType(['E', 6, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 2, 2]).index_set()
            (0, 1)
            sage: CartanType(['G', 2, 1]).index_set()
            (0, 1, 2)
            sage: CartanType(['F', 4, 1]).index_set()
            (0, 1, 2, 3, 4)
        """

    # This coloring scheme is used for crystal graphs and will eventually
    # be used for Coxeter groups etc. (experimental feature)
    _index_set_coloring = {1:"blue", 2:"red", 3:"green"}

    @abstract_method(optional=True)
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram for ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: CartanType(['B',3]).coxeter_diagram()
            Graph on 3 vertices
            sage: CartanType(['A',3]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 3)]
            sage: CartanType(['B',3]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 4)]
            sage: CartanType(['G',2]).coxeter_diagram().edges(sort=True)
            [(1, 2, 6)]
            sage: CartanType(['F',4]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
        """

    @cached_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix for ``self``.

        EXAMPLES::

            sage: CartanType(['A', 4]).coxeter_matrix()                                 # needs sage.graphs
            [1 3 2 2]
            [3 1 3 2]
            [2 3 1 3]
            [2 2 3 1]
        """
        from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
        return CoxeterMatrix(self)

    def coxeter_type(self):
        """
        Return the Coxeter type for ``self``.

        EXAMPLES::

            sage: CartanType(['A', 4]).coxeter_type()
            Coxeter type of ['A', 4]
        """
        from sage.combinat.root_system.coxeter_type import CoxeterType
        return CoxeterType(self)

    def dual(self):
        """
        Return the dual Cartan type, possibly just as a formal dual.

        EXAMPLES::

            sage: CartanType(['A',3]).dual()
            ['A', 3]
            sage: CartanType(["B", 3]).dual()
            ['C', 3]
            sage: CartanType(['C',2]).dual()
            ['B', 2]
            sage: CartanType(['D',4]).dual()
            ['D', 4]
            sage: CartanType(['E',8]).dual()
            ['E', 8]
            sage: CartanType(['F',4]).dual()
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
        """
        from . import type_dual
        return type_dual.CartanType(self)

    def relabel(self, relabelling):
        """
        Return a relabelled copy of this Cartan type.

        INPUT:

        - ``relabelling`` -- a function (or a list or dictionary)

        OUTPUT:

        an isomorphic Cartan type obtained by relabelling the nodes of
        the Dynkin diagram. Namely, the node with label ``i`` is
        relabelled ``f(i)`` (or, by ``f[i]`` if ``f`` is a list or
        dictionary).

        EXAMPLES::

           sage: CartanType(['F',4]).relabel({ 1:4, 2:3, 3:2, 4:1 }).dynkin_diagram()   # needs sage.graphs
           O---O=>=O---O
           4   3   2   1
           F4 relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
        """
        from . import type_relabel
        return type_relabel.CartanType(self, relabelling)

    def subtype(self, index_set):
        """
        Return a subtype of ``self`` given by ``index_set``.

        A subtype can be considered the Dynkin diagram induced from
        the Dynkin diagram of ``self`` by ``index_set``.

        EXAMPLES::

            sage: ct = CartanType(['A',6,2])
            sage: ct.dynkin_diagram()                                                   # needs sage.graphs
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: ct.subtype([1,2,3])                                                   # needs sage.graphs
            ['C', 3]
        """
        return self.cartan_matrix().subtype(index_set).cartan_type()

    def marked_nodes(self, marked_nodes):
        """
        Return a Cartan type with the nodes ``marked_nodes`` marked.

        INPUT:

        - ``marked_nodes`` -- list of nodes to mark

        EXAMPLES::

            sage: CartanType(['F',4]).marked_nodes([1, 3]).dynkin_diagram()             # needs sage.graphs
            X---O=>=X---O
            1   2   3   4
            F4 with nodes (1, 3) marked
        """
        if not marked_nodes:
            return self
        from . import type_marked
        return type_marked.CartanType(self, marked_nodes)

    def is_reducible(self):
        """
        Report whether the root system is reducible (i.e. not simple), that
        is whether it can be factored as a product of root systems.

        EXAMPLES::

            sage: CartanType("A2xB3").is_reducible()
            True
            sage: CartanType(['A',2]).is_reducible()
            False
        """
        return not self.is_irreducible()

    def is_irreducible(self):
        """
        Report whether this Cartan type is irreducible (i.e. simple). This
        should be overridden in any subclass.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_irreducible()
            False
        """
        return False

    def is_atomic(self):
        r"""
        This method is usually equivalent to :meth:`is_reducible`,
        except for the Cartan type `D_2`.

        `D_2` is not a standard Cartan type. It is equivalent to type
        `A_1 \times A_1` which is reducible; however the isomorphism
        from its ambient space (for the orthogonal group of degree 4)
        to that of `A_1 \times A_1` is non trivial, and it is useful to
        have it.

        From a programming point of view its implementation is more
        similar to the irreducible types, and so the method
        :meth:`is_atomic()` is supplied.

        EXAMPLES::

            sage: CartanType("D2").is_atomic()
            True
            sage: CartanType("D2").is_irreducible()
            False

        TESTS::

            sage: all( T.is_irreducible() == T.is_atomic() for T in CartanType.samples() if T != CartanType("D2"))
            True
        """
        return self.is_irreducible()

    def is_compound(self):
        """
        A short hand for not :meth:`is_atomic`.

        TESTS::

            sage: all( T.is_compound() == (not T.is_atomic()) for T in CartanType.samples())
            True
        """
        return not self.is_atomic()

    @abstract_method
    def is_finite(self):
        """
        Return whether this Cartan type is finite.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method is_finite at ...>

        ::

            sage: CartanType(['A',4]).is_finite()
            True
            sage: CartanType(['A',4, 1]).is_finite()
            False
        """

    @abstract_method
    def is_affine(self):
        """
        Return whether ``self`` is affine.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """

    def is_crystallographic(self):
        """
        Return whether this Cartan type is crystallographic.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_crystallographic() ] for t in CartanType.samples(finite=True) ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], True],
             [['C', 1], True], [['C', 5], True],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], True], [['G', 2], True],
             [['I', 5], False], [['H', 3], False], [['H', 4], False]]
        """
        return False

    def is_simply_laced(self):
        """
        Return whether this Cartan type is simply laced.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_simply_laced() ] for t in CartanType.samples() ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], False],
             [['C', 1], True], [['C', 5], False],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], False], [['G', 2], False], [['I', 5], False],
             [['H', 3], False], [['H', 4], False],
             [['A', 1, 1], False], [['A', 5, 1], True],
             [['B', 1, 1], False], [['B', 5, 1], False],
             [['C', 1, 1], False], [['C', 5, 1], False],
             [['D', 3, 1], True], [['D', 5, 1], True],
             [['E', 6, 1], True], [['E', 7, 1], True], [['E', 8, 1], True],
             [['F', 4, 1], False], [['G', 2, 1], False],
             [['BC', 1, 2], False], [['BC', 5, 2], False],
             [['B', 5, 1]^*, False], [['C', 4, 1]^*, False],
             [['F', 4, 1]^*, False], [['G', 2, 1]^*, False],
             [['BC', 1, 2]^*, False], [['BC', 5, 2]^*, False]]
        """
        return False

    def is_implemented(self):
        """
        Check whether the Cartan datum for ``self`` is actually implemented.

        EXAMPLES::

            sage: CartanType(["A",4,1]).is_implemented()                                # needs sage.graphs
            True
            sage: CartanType(['H',3]).is_implemented()
            True
        """
        try:
            self.coxeter_diagram()
            return True
        except Exception:
            return False

    def root_system(self):
        """
        Return the root system associated to ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).root_system()
            Root system of type ['A', 4]
        """
        from sage.combinat.root_system.root_system import RootSystem
        return RootSystem(self)

    def as_folding(self, folding_of=None, sigma=None):
        r"""
        Return ``self`` realized as a folded Cartan type.

        For finite and affine types, this is realized by the Dynkin
        diagram foldings:

        .. MATH::

            \begin{array}{ccl}
            C_n^{(1)}, A_{2n}^{(2)}, A_{2n}^{(2)\dagger}, D_{n+1}^{(2)}
            & \hookrightarrow & A_{2n-1}^{(1)}, \\
            A_{2n-1}^{(2)}, B_n^{(1)} & \hookrightarrow & D_{n+1}^{(1)}, \\
            E_6^{(2)}, F_4^{(1)} & \hookrightarrow & E_6^{(1)}, \\
            D_4^{(3)}, G_2^{(1)} & \hookrightarrow & D_4^{(1)}, \\
            C_n & \hookrightarrow & A_{2n-1}, \\
            B_n & \hookrightarrow & D_{n+1}, \\
            F_4 & \hookrightarrow & E_6, \\
            G_2 & \hookrightarrow & D_4.
            \end{array}

        For general types, this returns ``self`` as a folded type of ``self``
        with `\sigma` as the identity map.

        For more information on these foldings and folded Cartan types, see
        :class:`sage.combinat.root_system.type_folded.CartanTypeFolded`.

        If the optional inputs ``folding_of`` and ``sigma`` are specified, then
        this returns the folded Cartan type of ``self`` in ``folding_of`` given
        by the automorphism ``sigma``.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).as_folding()
            ['B', 3, 1] as a folding of  ['D', 4, 1]
            sage: CartanType(['F', 4]).as_folding()
            ['F', 4] as a folding of  ['E', 6]
            sage: CartanType(['BC', 3, 2]).as_folding()
            ['BC', 3, 2] as a folding of  ['A', 5, 1]
            sage: CartanType(['D', 4, 3]).as_folding()
            ['G', 2, 1]^* relabelled by {0: 0, 1: 2, 2: 1} as a folding of ['D', 4, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        if folding_of is None and sigma is None:
            return self._default_folded_cartan_type()
        if folding_of is None or sigma is None:
            raise ValueError("Both folding_of and sigma must be given")
        return CartanTypeFolded(self, folding_of, sigma)

    def _default_folded_cartan_type(self):
        r"""
        Return the default folded Cartan type.

        In general, this just returns ``self`` in ``self`` with `\sigma` as
        the identity map.

        EXAMPLES::

            sage: D = CartanMatrix([[2, -3], [-2, 2]]).dynkin_diagram()                 # needs sage.graphs
            sage: D._default_folded_cartan_type()                                       # needs sage.graphs
            Dynkin diagram of rank 2 as a folding of  Dynkin diagram of rank 2
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        return CartanTypeFolded(self, self, [[i] for i in self.index_set()])

    options = CartanType.options


class CartanType_crystallographic(CartanType_abstract):
    """
    An abstract class for crystallographic Cartan types.
    """
    # The default value should really be lambda x:x, but sphinx does
    # not like it currently (see #14553); since this is an abstract method
    # the value won't actually be used, so we put a fake instead.
    @abstract_method(optional=True)
    def ascii_art(self, label='lambda x: x', node=None):
        r"""
        Return an ascii art representation of the Dynkin diagram.

        INPUT:

        - ``label`` -- (default: the identity) a relabeling function
          for the nodes
        - ``node`` -- (optional) a function which returns
          the character for a node

        EXAMPLES::

            sage: cartan_type = CartanType(['B',5,1])
            sage: print(cartan_type.ascii_art())
                O 0
                |
                |
            O---O---O---O=>=O
            1   2   3   4   5

        The label option is useful to visualize various statistics on
        the nodes of the Dynkin diagram::

            sage: a = cartan_type.col_annihilator(); a                                  # needs sage.graphs
            Finite family {0: 1, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2}
            sage: print(CartanType(['B',5,1]).ascii_art(label=a.__getitem__))           # needs sage.graphs
                O 1
                |
                |
            O---O---O---O=>=O
            1   2   2   2   2
        """

    # The default value of label should really be lambda i:i, but sphinx does
    # not like it currently (see #14553); since this is an abstract method
    # the value won't actually be used, so we put a fake instead.
    @abstract_method(optional=True)
    def _latex_dynkin_diagram(self, label='lambda i: i',
                              node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        INPUT:

        - ``label`` -- (default: the identity) a relabeling function
          for the nodes

        - ``node`` -- (optional) a function which returns the latex for a node

        - ``node_dist`` -- (default: 2) the distance between nodes in cm

        EXAMPLES::

            sage: latex(CartanType(['A',4]).dynkin_diagram())  # indirect doctest       # needs sage.graphs
            \begin{tikzpicture}[scale=0.5]
            \draw (-1,0) node[anchor=east] {$A_{4}$};
            \draw (0 cm,0) -- (6 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            \end{tikzpicture}
        """

    @abstract_method
    def dynkin_diagram(self):
        """
        Return the Dynkin diagram associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).dynkin_diagram()                                  # needs sage.graphs
            O---O---O---O
            1   2   3   4
            A4

        .. NOTE::

            Derived subclasses should typically implement this as a cached
            method.
        """

    @cached_method
    def cartan_matrix(self):
        """
        Return the Cartan matrix associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).cartan_matrix()                                   # needs sage.graphs
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        from sage.combinat.root_system.cartan_matrix import CartanMatrix
        return CartanMatrix(self.dynkin_diagram())

    def coxeter_diagram(self):
        """
        Return the Coxeter diagram for ``self``.

        This implementation constructs it from the Dynkin diagram.

        .. SEEALSO:: :meth:`CartanType_abstract.coxeter_diagram`

        EXAMPLES::

            sage: # needs sage.graphs
            sage: CartanType(['A',3]).coxeter_diagram()
            Graph on 3 vertices
            sage: CartanType(['A',3]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 3)]
            sage: CartanType(['B',3]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 4)]
            sage: CartanType(['G',2]).coxeter_diagram().edges(sort=True)
            [(1, 2, 6)]
            sage: CartanType(['F',4]).coxeter_diagram().edges(sort=True)
            [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
            sage: CartanType(['A',2,2]).coxeter_diagram().edges(sort=True)
            [(0, 1, +Infinity)]
        """
        return self.dynkin_diagram().coxeter_diagram()

    def is_crystallographic(self):
        """
        Implement :meth:`CartanType_abstract.is_crystallographic`
        by returning ``True``.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_crystallographic()
            True
        """
        return True

    @cached_method
    def symmetrizer(self):
        """
        Return the symmetrizer of the Cartan matrix of ``self``.

        A Cartan matrix `M` is symmetrizable if there exists a non
        trivial diagonal matrix `D` such that `DM` is a symmetric
        matrix, that is `DM = M^tD`. In that case, `D` is unique, up
        to a scalar factor for each connected component of the Dynkin
        diagram.

        This method computes the unique minimal such `D` with positive
        integral coefficients. If `D` exists, it is returned as a
        family. Otherwise ``None`` is returned.

        The coefficients are coerced to ``base_ring``.

        EXAMPLES::

            sage: CartanType(["B",5]).symmetrizer()                                     # needs sage.graphs
            Finite family {1: 2, 2: 2, 3: 2, 4: 2, 5: 1}

        Here is a neat trick to visualize it better::

            sage: T = CartanType(["B",5])
            sage: print(T.ascii_art(T.symmetrizer().__getitem__))                       # needs sage.graphs
            O---O---O---O=>=O
            2   2   2   2   1

            sage: T = CartanType(["BC",5, 2])
            sage: print(T.ascii_art(T.symmetrizer().__getitem__))                       # needs sage.graphs
            O=<=O---O---O---O=<=O
            1   2   2   2   2   4

        Here is the symmetrizer of some reducible Cartan types::

            sage: T = CartanType(["D", 2])
            sage: print(T.ascii_art(T.symmetrizer().__getitem__))                       # needs sage.graphs
            O   O
            1   1

            sage: T = CartanType(["B",5],["BC",5, 2])
            sage: print(T.ascii_art(T.symmetrizer().__getitem__))                       # needs sage.graphs
            O---O---O---O=>=O
            2   2   2   2   1
            O=<=O---O---O---O=<=O
            1   2   2   2   2   4

        Property: up to an overall scalar factor, this gives the norm
        of the simple roots in the ambient space::

            sage: T = CartanType(["C",5])
            sage: print(T.ascii_art(T.symmetrizer().__getitem__))                       # needs sage.graphs
            O---O---O---O=<=O
            1   1   1   1   2

            sage: alpha = RootSystem(T).ambient_space().simple_roots()
            sage: print(T.ascii_art(lambda i: alpha[i].scalar(alpha[i])))
            O---O---O---O=<=O
            2   2   2   2   4
        """
        from sage.matrix.constructor import matrix, diagonal_matrix
        m = self.cartan_matrix()
        n = m.nrows()
        M = matrix(ZZ, n, n*n, sparse=True)
        for (i,j) in m.nonzero_positions():
            M[i, n * i + j] = m[i,j]
            M[j, n * i + j] -= m[j,i]
        kern = M.integer_kernel()
        c = len(self.dynkin_diagram().connected_components(sort=False))
        if kern.dimension() < c:
            # the Cartan matrix is not symmetrizable
            return None
        assert kern.dimension() == c
        # Now the basis contains one vector v per connected component
        # C of the Dynkin diagram, or equivalently diagonal block of
        # the Cartan matrix. The support of v is exactly that
        # connected component, and it symmetrizes the corresponding
        # diagonal block of the Cartan matrix. We sum all those vectors.
        D = sum(kern.basis())
        assert diagonal_matrix(D) * m == m.transpose() * diagonal_matrix(D)
        I = self.index_set()
        return Family( dict( (I[i], D[i]) for i in range(n) ) )

    def index_set_bipartition(self):
        r"""
        Return a bipartition `\{L,R\}` of the vertices of the Dynkin diagram.

        For `i` and `j` both in `L` (or both in `R`), the simple
        reflections `s_i` and `s_j` commute.

        Of course, the Dynkin diagram should be bipartite. This is
        always the case for all finite types.

        EXAMPLES::

            sage: CartanType(['A',5]).index_set_bipartition()                           # needs sage.graphs
            ({1, 3, 5}, {2, 4})

            sage: CartanType(['A',2,1]).index_set_bipartition()                         # needs sage.graphs
            Traceback (most recent call last):
            ...
            ValueError: the Dynkin diagram must be bipartite
        """
        from sage.graphs.graph import Graph
        G = Graph(self.dynkin_diagram())
        if not G.is_bipartite():
            raise ValueError("the Dynkin diagram must be bipartite")
        return G.bipartite_sets()


class CartanType_simply_laced(CartanType_crystallographic):
    """
    An abstract class for simply laced Cartan types.
    """

    def is_simply_laced(self):
        """
        Return whether ``self`` is simply laced, which is ``True``.

        EXAMPLES::

            sage: CartanType(['A',3,1]).is_simply_laced()
            True
            sage: CartanType(['A',2]).is_simply_laced()
            True
        """
        return True

    def dual(self):
        """
        Simply laced Cartan types are self-dual, so return ``self``.

        EXAMPLES::

            sage: CartanType(["A", 3]).dual()
            ['A', 3]
            sage: CartanType(["A", 3, 1]).dual()
            ['A', 3, 1]
            sage: CartanType(["D", 3]).dual()
            ['D', 3]
            sage: CartanType(["D", 4, 1]).dual()
            ['D', 4, 1]
            sage: CartanType(["E", 6]).dual()
            ['E', 6]
            sage: CartanType(["E", 6, 1]).dual()
            ['E', 6, 1]
        """
        return self


class CartanType_simple(CartanType_abstract):
    """
    An abstract class for simple Cartan types.
    """

    def is_irreducible(self):
        """
        Return whether ``self`` is irreducible, which is ``True``.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_irreducible()
            True
        """
        return True


class CartanType_finite(CartanType_abstract):
    """
    An abstract class for simple affine Cartan types.
    """

    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_finite()
            True
        """
        return True

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_affine()
            False
        """
        return False


class CartanType_affine(CartanType_simple, CartanType_crystallographic):
    """
    An abstract class for simple affine Cartan types
    """

    AmbientSpace = LazyImport('sage.combinat.root_system.type_affine', 'AmbientSpace')

    def _ascii_art_node(self, label):
        """
        Return the ascii art for the node labeled by ``label``.

        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: CartanType.options(mark_special_node='both')
            sage: ct._ascii_art_node(0)
            '@'
            sage: CartanType.options._reset()
        """
        if (label == self.special_node()
                and self.options('mark_special_node') in ['printing', 'both']):
            return self.options('special_node_str')
        return super()._ascii_art_node(label)

    def _latex_draw_node(self, x, y, label, position='below=4pt'):
        r"""
        Draw (possibly marked [crossed out]) circular node ``i`` at the
        position ``(x,y)`` with node label ``label`` .

        - ``position`` -- position of the label relative to the node
        - ``anchor`` -- (optional) the anchor point for the label

        EXAMPLES::

            sage: CartanType.options(mark_special_node='both')
            sage: CartanType(['A',3,1])._latex_draw_node(0, 0, 0)
            '\\draw[fill=black] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$0$};\n'
            sage: CartanType.options._reset()
        """
        if label == self.special_node() and self.options('mark_special_node') in ['latex', 'both']:
            fill = 'black'
        else:
            fill = 'white'
        return super()._latex_draw_node(x, y, label, position, fill)

    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_finite()
            False
        """
        return False

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """
        return True

    def is_untwisted_affine(self):
        """
        Return whether ``self`` is untwisted affine.

        A Cartan type is untwisted affine if it is the canonical
        affine extension of some finite type. Every affine type is
        either untwisted affine, dual thereof, or of type ``BC``.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_untwisted_affine()
            True
            sage: CartanType(['A', 3, 1]).dual().is_untwisted_affine()  # this one is self dual!
            True
            sage: CartanType(['B', 3, 1]).dual().is_untwisted_affine()
            False
            sage: CartanType(['BC', 3, 2]).is_untwisted_affine()
            False
        """
        return False

    @abstract_method
    def special_node(self):
        r"""
        Return a special node of the Dynkin diagram.

        A *special* node is a node of the Dynkin diagram such that
        pruning it yields a Dynkin diagram for the associated
        classical type (see :meth:`classical`).

        This method returns the label of some special node. This is
        usually `0` in the standard conventions.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).special_node()
            0

        The choice is guaranteed to be consistent with the indexing of
        the nodes of the classical Dynkin diagram::

            sage: CartanType(['A', 3, 1]).index_set()
            (0, 1, 2, 3)
            sage: CartanType(['A', 3, 1]).classical().index_set()
            (1, 2, 3)
        """

    @cached_method
    def special_nodes(self):
        r"""
        Return the set of special nodes of the affine Dynkin diagram.

        EXAMPLES::

            sage: # needs sage.graphs sage.groups
            sage: CartanType(['A',3,1]).special_nodes()
            (0, 1, 2, 3)
            sage: CartanType(['C',2,1]).special_nodes()
            (0, 2)
            sage: CartanType(['D',4,1]).special_nodes()
            (0, 1, 3, 4)
            sage: CartanType(['E',6,1]).special_nodes()
            (0, 1, 6)
            sage: CartanType(['D',3,2]).special_nodes()
            (0, 2)
            sage: CartanType(['A',4,2]).special_nodes()
            (0,)
        """
        return tuple(sorted(self.dynkin_diagram().automorphism_group(edge_labels=True).orbit(self.special_node())))

    @abstract_method
    def classical(self):
        r"""
        Return the classical Cartan type associated with this affine Cartan type.

        EXAMPLES::

            sage: CartanType(['A', 1, 1]).classical()
            ['A', 1]
            sage: CartanType(['A', 3, 1]).classical()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).classical()
            ['B', 3]

            sage: CartanType(['A', 2, 2]).classical()
            ['C', 1]
            sage: CartanType(['BC', 1, 2]).classical()
            ['C', 1]
            sage: CartanType(['A', 4, 2]).classical()
            ['C', 2]
            sage: CartanType(['BC', 2, 2]).classical()
            ['C', 2]
            sage: CartanType(['A', 10, 2]).classical()
            ['C', 5]
            sage: CartanType(['BC', 5, 2]).classical()
            ['C', 5]

            sage: CartanType(['D', 5, 2]).classical()
            ['B', 4]
            sage: CartanType(['E', 6, 1]).classical()
            ['E', 6]
            sage: CartanType(['G', 2, 1]).classical()
            ['G', 2]
            sage: CartanType(['E', 6, 2]).classical()
            ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}
            sage: CartanType(['D', 4, 3]).classical()
            ['G', 2]

        We check that :meth:`classical`,
        :meth:`sage.combinat.root_system.cartan_type.CartanType_crystallographic.dynkin_diagram`,
        and :meth:`special_node` are consistent::

            sage: for ct in CartanType.samples(affine=True):                            # needs sage.graphs
            ....:     g1 = ct.classical().dynkin_diagram()
            ....:     g2 = ct.dynkin_diagram()
            ....:     g2.delete_vertex(ct.special_node())
            ....:     assert g1.vertices(sort=True) == g2.vertices(sort=True)
            ....:     assert g1.edges(sort=True) == g2.edges(sort=True)
        """

    @abstract_method
    def basic_untwisted(self):
        r"""
        Return the basic untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic untwisted type is `X_n`.
        In other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: CartanType(['A', 1, 1]).basic_untwisted()
            ['A', 1]
            sage: CartanType(['A', 3, 1]).basic_untwisted()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).basic_untwisted()
            ['B', 3]
            sage: CartanType(['E', 6, 1]).basic_untwisted()
            ['E', 6]
            sage: CartanType(['G', 2, 1]).basic_untwisted()
            ['G', 2]

            sage: CartanType(['A', 2, 2]).basic_untwisted()
            ['A', 2]
            sage: CartanType(['A', 4, 2]).basic_untwisted()
            ['A', 4]
            sage: CartanType(['A', 11, 2]).basic_untwisted()
            ['A', 11]
            sage: CartanType(['D', 5, 2]).basic_untwisted()
            ['D', 5]
            sage: CartanType(['E', 6, 2]).basic_untwisted()
            ['E', 6]
            sage: CartanType(['D', 4, 3]).basic_untwisted()
            ['D', 4]
        """

    def row_annihilator(self, m=None):
        r"""
        Return the unique minimal non trivial annihilating linear
        combination of `\alpha_0, \alpha_1, \ldots, \alpha_n` with
        nonnegative coefficients (or alternatively, the unique minimal
        non trivial annihilating linear combination of the rows of the
        Cartan matrix with nonnegative coefficients).

        Throw an error if the existence of uniqueness does not hold

        The optional argument ``m`` is for internal use only.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: RootSystem(['C',2,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 1, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 1, 2: 2, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().acheck()
            Finite family {0: 1, 1: 2, 2: 3, 3: 2, 4: 1}
            sage: RootSystem(['BC',4,2]).cartan_type().acheck()
            Finite family {0: 1, 1: 2, 2: 2, 3: 2, 4: 2}

        ``acheck`` is a shortcut for row_annihilator::

            sage: RootSystem(['BC',4,2]).cartan_type().row_annihilator()                # needs sage.graphs
            Finite family {0: 1, 1: 2, 2: 2, 3: 2, 4: 2}

        FIXME:

        - The current implementation assumes that the Cartan matrix
          is indexed by `[0,1,...]`, in the same order as the index set.
        - This really should be a method of :class:`CartanMatrix`.
        """
        if m is None:
            m = self.cartan_matrix()
        if self.index_set() != tuple(range(m.ncols())):
            raise NotImplementedError("the Cartan matrix currently must be indexed by [0,1,...,{}]".format(m.ncols()))
        annihilator_basis = m.integer_kernel().gens()
        if len(annihilator_basis) != 1:
            raise ValueError("the kernel is not 1 dimensional")
        assert (all(coef > 0 for coef in annihilator_basis[0]))

        return Family(dict((i,annihilator_basis[0][i])for i in self.index_set()))

    acheck = row_annihilator

    def col_annihilator(self):
        r"""
        Return the unique minimal non trivial annihilating linear
        combination of `\alpha^\vee_0, \alpha^\vee, \ldots, \alpha^\vee` with
        nonnegative coefficients (or alternatively, the unique minimal
        non trivial annihilating linear combination of the columns of the
        Cartan matrix with nonnegative coefficients).

        Throw an error if the existence or uniqueness does not hold

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: RootSystem(['C',2,1]).cartan_type().a()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().a()
            Finite family {0: 1, 1: 1, 2: 2, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().a()
            Finite family {0: 1, 1: 2, 2: 3, 3: 4, 4: 2}
            sage: RootSystem(['BC',4,2]).cartan_type().a()
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}

        ``a`` is a shortcut for col_annihilator::

            sage: RootSystem(['BC',4,2]).cartan_type().col_annihilator()                # needs sage.graphs
            Finite family {0: 2, 1: 2, 2: 2, 3: 2, 4: 1}
        """
        return self.row_annihilator(self.cartan_matrix().transpose())

    a = col_annihilator

    def c(self):
        r"""
        Return the family (c_i)_i of integer coefficients defined by
        `c_i=max(1, a_i/a^vee_i)` (see e.g. [FSS07]_ p. 3)

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: RootSystem(['C',2,1]).cartan_type().c()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: RootSystem(['D',4,1]).cartan_type().c()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
            sage: RootSystem(['F',4,1]).cartan_type().c()
            Finite family {0: 1, 1: 1, 2: 1, 3: 2, 4: 2}
            sage: RootSystem(['BC',4,2]).cartan_type().c()
            Finite family {0: 2, 1: 1, 2: 1, 3: 1, 4: 1}

        TESTS::

            sage: CartanType(["B", 3, 1]).c().map(parent)                               # needs sage.graphs
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Integer Ring}

        REFERENCES:

        .. [FSS07] \G. Fourier, A. Schilling, and M. Shimozono,
           *Demazure structure inside Kirillov-Reshetikhin crystals*,
           J. Algebra, Vol. 309, (2007), p. 386-404
           :arxiv:`math/0605451`
        """
        a = self.a()
        acheck = self.acheck()
        return Family(dict((i, max(ZZ(1), a[i] // acheck[i]))
                           for i in self.index_set()))

    def translation_factors(self):
        r"""
        Return the translation factors for ``self``.

        Those are the smallest factors `t_i` such that the translation
        by `t_i \alpha_i` maps the fundamental polygon to another
        polygon in the alcove picture.

        OUTPUT:

        a dictionary from ``self.index_set()`` to `\ZZ`
        (or `\QQ` for affine type `BC`)

        Those coefficients are all `1` for dual untwisted, and in
        particular for simply laced. They coincide with the usual
        `c_i` coefficients (see :meth:`c`) for untwisted and dual
        thereof. See the discussion below for affine type `BC`.

        .. NOTE::

            One usually realizes the alcove picture in the coweight
            lattice, with translations by coroots; in that case, one will
            use the translation factors for the dual Cartan type.

        FIXME: the current implementation assumes that the Cartan
        matrix is indexed by `[0,1,...]`, in the same order as the
        index set.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: CartanType(['C',2,1]).translation_factors()
            Finite family {0: 1, 1: 2, 2: 1}
            sage: CartanType(['C',2,1]).dual().translation_factors()
            Finite family {0: 1, 1: 1, 2: 1}
            sage: CartanType(['D',4,1]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
            sage: CartanType(['F',4,1]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 2, 4: 2}
            sage: CartanType(['BC',4,2]).translation_factors()
            Finite family {0: 1, 1: 1, 2: 1, 3: 1, 4: 1/2}

        We proceed with systematic tests taken from MuPAD-Combinat's
        testsuite::

            sage: # needs sage.graphs
            sage: list(CartanType(["A", 1, 1]).translation_factors())
            [1, 1]
            sage: list(CartanType(["A", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["B", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 2]
            sage: list(CartanType(["C", 5, 1]).translation_factors())
            [1, 2, 2, 2, 2, 1]
            sage: list(CartanType(["D", 5, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 6, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 7, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["E", 8, 1]).translation_factors())
            [1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["F", 4, 1]).translation_factors())
            [1, 1, 1, 2, 2]
            sage: list(CartanType(["G", 2, 1]).translation_factors())
            [1, 3, 1]
            sage: list(CartanType(["A", 2, 2]).translation_factors())
            [1, 1/2]
            sage: list(CartanType(["A", 2, 2]).dual().translation_factors())
            [1/2, 1]
            sage: list(CartanType(["A", 10, 2]).translation_factors())
            [1, 1, 1, 1, 1, 1/2]
            sage: list(CartanType(["A", 10, 2]).dual().translation_factors())
            [1/2, 1, 1, 1, 1, 1]
            sage: list(CartanType(["A", 9, 2]).translation_factors())
            [1, 1, 1, 1, 1, 1]
            sage: list(CartanType(["D", 5, 2]).translation_factors())
            [1, 1, 1, 1, 1]
            sage: list(CartanType(["D", 4, 3]).translation_factors())
            [1, 1, 1]
            sage: list(CartanType(["E", 6, 2]).translation_factors())
            [1, 1, 1, 1, 1]

        We conclude with a discussion of the appropriate value for
        affine type `BC`. Let us consider the alcove picture realized
        in the weight lattice. It is obtained by taking the level-`1`
        affine hyperplane in the weight lattice, and projecting it
        along `\Lambda_0`::

            sage: R = RootSystem(["BC",2,2])
            sage: alpha = R.weight_space().simple_roots()                               # needs sage.graphs
            sage: alphacheck = R.coroot_space().simple_roots()
            sage: Lambda = R.weight_space().fundamental_weights()

        Here are the levels of the fundamental weights::

            sage: Lambda[0].level(), Lambda[1].level(), Lambda[2].level()               # needs sage.graphs
            (1, 2, 2)

        So the "center" of the fundamental polygon at level `1` is::

            sage: O = Lambda[0]
            sage: O.level()                                                             # needs sage.graphs
            1

        We take the projection `\omega_1` at level `0` of `\Lambda_1`
        as unit vector on the `x`-axis, and the projection `\omega_2`
        at level 0 of `\Lambda_2` as unit vector of the `y`-axis::

            sage: omega1 = Lambda[1] - 2*Lambda[0]
            sage: omega2 = Lambda[2] - 2*Lambda[0]
            sage: omega1.level(), omega2.level()                                        # needs sage.graphs
            (0, 0)

        The projections of the simple roots can be read off::

            sage: alpha[0]                                                              # needs sage.graphs
            2*Lambda[0] - Lambda[1]
            sage: alpha[1]                                                              # needs sage.graphs
            -2*Lambda[0] + 2*Lambda[1] - Lambda[2]
            sage: alpha[2]                                                              # needs sage.graphs
            -2*Lambda[1] + 2*Lambda[2]

        Namely `\alpha_0 = -\omega_1`, `\alpha_1 = 2\omega_1 -
        \omega_2` and `\alpha_2 = -2 \omega_1 + 2 \omega_2`.

        The reflection hyperplane defined by `\alpha_0^\vee` goes
        through the points `O+1/2 \omega_1` and `O+1/2 \omega_2`::

            sage: (O+(1/2)*omega1).scalar(alphacheck[0])
            0
            sage: (O+(1/2)*omega2).scalar(alphacheck[0])
            0

        Hence, the fundamental alcove is the triangle `(O, O+1/2
        \omega_1, O+1/2 \omega_2)`. By successive reflections, one can
        tile the full plane. This induces a tiling of the full plane
        by translates of the fundamental polygon.

        .. TODO::

            Add the picture here, once root system plots in the
            weight lattice will be implemented. In the mean time, the
            reader may look up the dual picture on Figure 2 of [HST09]_
            which was produced by MuPAD-Combinat.

        From this picture, one can read that translations by
        `\alpha_0`, `\alpha_1`, and `1/2\alpha_2` map the fundamental
        polygon to translates of it in the alcove picture, and are
        smallest with this property. Hence, the translation factors
        for affine type `BC` are `t_0=1, t_1=1, t_2=1/2`::

            sage: CartanType(['BC',2,2]).translation_factors()                          # needs sage.graphs
            Finite family {0: 1, 1: 1, 2: 1/2}

        TESTS::

            sage: CartanType(["B", 3, 1]).translation_factors().map(parent)             # needs sage.graphs
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Integer Ring}
            sage: CartanType(["BC", 3, 2]).translation_factors().map(parent)            # needs sage.graphs
            Finite family {0: Integer Ring, 1: Integer Ring, 2: Integer Ring, 3: Rational Field}

        REFERENCES:

        .. [HST09] \F. Hivert, A. Schilling, and N. M. Thiery,
           *Hecke group algebras as quotients of affine Hecke
           algebras at level 0*, JCT A, Vol. 116, (2009) p. 844-863
           :arxiv:`0804.3781`
        """
        a = self.a()
        acheck = self.acheck()
        if set([1/ZZ(2), 2]).issubset( set(a[i]/acheck[i] for i in self.index_set()) ):
            # The test above and the formula below are rather meaningless
            # But they detect properly type BC or dual and return the correct value
            return Family(dict((i, min(ZZ(1), a[i] / acheck[i]))
                               for i in self.index_set()))

        else:
            return self.c()

    def _test_dual_classical(self, **options):
        r"""
        Test whether the special node of the dual is still the same and whether
        the methods dual and classical commute.

        TESTS::

            sage: C = CartanType(['A',2,2])
            sage: C._test_dual_classical()
        """
        tester = self._tester(**options)
        tester.assertEqual(self.classical().dual(), self.dual().classical())
        tester.assertEqual(self.special_node(), self.dual().special_node())

    def other_affinization(self):
        """
        Return the other affinization of the same classical type.

        EXAMPLES::

            sage: CartanType(["A", 3, 1]).other_affinization()
            ['A', 3, 1]
            sage: CartanType(["B", 3, 1]).other_affinization()
            ['C', 3, 1]^*
            sage: CartanType(["C", 3, 1]).dual().other_affinization()
            ['B', 3, 1]

        Is this what we want?::

            sage: CartanType(["BC", 3, 2]).dual().other_affinization()
            ['B', 3, 1]
        """
        if self.is_untwisted_affine():
            result = self.classical().dual().affine().dual()
        else:
            result = self.dual().classical().dual().affine()
        assert result.classical() is self.classical()
        return result

##############################################################################
# Concrete base classes


class CartanType_standard(UniqueRepresentation, SageObject):
    # Technical methods
    def _repr_(self, compact=False):
        """
        TESTS::

            sage: ct = CartanType(['A',3])
            sage: repr(ct)
            "['A', 3]"
            sage: ct._repr_(compact=True)
            'A3'
        """
        format = '%s%s' if compact else "['%s', %s]"
        return format % (self.letter, self.n)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',4]))
            2
            sage: len(CartanType(['A',4,1]))
            3
        """
        return 3 if self.is_affine() else 2

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['B', 3])
            sage: t[0]
            'B'
            sage: t[1]
            3
            sage: t[2]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i == 1:
            return self.n
        else:
            raise IndexError("index out of range")


class CartanType_standard_finite(CartanType_standard, CartanType_finite):
    """
    A concrete base class for the finite standard Cartan types.

    This includes for example `A_3`, `D_4`, or `E_8`.

     TESTS::

         sage: ct1 = CartanType(['A',4])
         sage: ct2 = CartanType(['A',4])
         sage: ct3 = CartanType(['A',5])
         sage: ct1 == ct2
         True
         sage: ct1 != ct3
         True
    """

    def __init__(self, letter, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])

        TESTS::

            sage: TestSuite(ct).run(verbose = True)
            running ._test_category() . . . pass
            running ._test_new() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
        """
#         assert(t[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])
#         assert(t[1] in ZZ and t[1] >= 0)
#         if t[0] in ['B', 'C']:
#             assert(t[1] >= 2)
#         if t[0] == 'D':
#             assert(t[1] >= 3)
#         if t[0] == 'E':
#             assert(t[1] <= 8)
#         if t[0] == 'F':
#             assert(t[1] <= 4)
#         if t[0] == 'G':
#             assert(t[1] <= 2)
#         if t[0] == 'H':
#             assert(t[1] <= 4)
        self.letter = letter
        self.n = n

    def __reduce__(self):
        """
        TESTS::

            sage: T = CartanType(['D', 4])
            sage: T.__reduce__()
            (CartanType, ('D', 4))
            sage: T == loads(dumps(T))
            True
        """
        return (CartanType, (self.letter, self.n))

    def __hash__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',2])
            sage: hash(ct) #random
            -5684143898951441983
        """
        return hash((self.n, self.letter))

    # mathematical methods

    def index_set(self):
        r"""
        Implement :meth:`CartanType_abstract.index_set`.

        The index set for all standard finite Cartan types is of the form
        `\{1, \ldots, n\}`. (See :mod:`~sage.combinat.root_system.type_I`
        for a slight abuse of this).

        EXAMPLES::

            sage: CartanType(['A', 5]).index_set()
            (1, 2, 3, 4, 5)
        """
        return tuple(range(1,self.n+1))

    def rank(self):
        """
        Return the rank of ``self`` which for type `X_n` is `n`.

        EXAMPLES::

            sage: CartanType(['A', 3]).rank()
            3
            sage: CartanType(['B', 3]).rank()
            3
            sage: CartanType(['C', 3]).rank()
            3
            sage: CartanType(['D', 4]).rank()
            4
            sage: CartanType(['E', 6]).rank()
            6
        """
        return self.n

    def affine(self):
        """
        Return the corresponding untwisted affine Cartan type.

        EXAMPLES::

            sage: CartanType(['A',3]).affine()
            ['A', 3, 1]
        """
        return CartanType([self.letter, self.n, 1])

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        The Coxeter number is the order of a Coxeter element of the
        corresponding Weyl group.

        See Bourbaki, Lie Groups and Lie Algebras V.6.1 or
        :wikipedia:`Coxeter_element` for more information.

        EXAMPLES::

            sage: CartanType(['A',4]).coxeter_number()
            5
            sage: CartanType(['B',4]).coxeter_number()
            8
            sage: CartanType(['C',4]).coxeter_number()
            8
        """
        return sum(self.affine().a())

    def dual_coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A',4]).dual_coxeter_number()
            5
            sage: CartanType(['B',4]).dual_coxeter_number()
            7
            sage: CartanType(['C',4]).dual_coxeter_number()
            5
        """
        return sum(self.affine().acheck())

    def type(self):
        """
        Return the type of ``self``.

        EXAMPLES::

            sage: CartanType(['A', 4]).type()
            'A'
            sage: CartanType(['A', 4, 1]).type()
            'A'
        """
        return self.letter

    @cached_method
    def opposition_automorphism(self):
        r"""
        Return the opposition automorphism.

        The *opposition automorphism* is the automorphism
        `i \mapsto i^*` of the vertices Dynkin diagram such that,
        for `w_0` the longest element of the Weyl group, and any
        simple root `\alpha_i`, one has `\alpha_{i^*} = -w_0(\alpha_i)`.

        The automorphism is returned as a :class:`Family`.

        EXAMPLES::

            sage: ct = CartanType(['A', 5])
            sage: ct.opposition_automorphism()                                          # needs sage.libs.gap
            Finite family {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}

            sage: ct = CartanType(['D', 4])
            sage: ct.opposition_automorphism()                                          # needs sage.libs.gap
            Finite family {1: 1, 2: 2, 3: 3, 4: 4}

            sage: ct = CartanType(['D', 5])
            sage: ct.opposition_automorphism()                                          # needs sage.libs.gap
            Finite family {1: 1, 2: 2, 3: 3, 4: 5, 5: 4}

            sage: ct = CartanType(['C', 4])
            sage: ct.opposition_automorphism()                                          # needs sage.libs.gap
            Finite family {1: 1, 2: 2, 3: 3, 4: 4}
        """
        Q = self.root_system().root_lattice()
        W = Q.weyl_group()
        w0 = W.long_element()
        alpha = Q.simple_roots()
        d = {i: (w0.action(alpha[i])).leading_support() for i in self.index_set()}
        return Family(d)

##########################################################################


class CartanType_standard_affine(CartanType_standard, CartanType_affine):
    r"""
    A concrete class for affine simple Cartan types.
    """

    def __init__(self, letter, n, affine=1):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: TestSuite(ct).run()

        TESTS::

            sage: ct1 = CartanType(['A',3, 1])
            sage: ct2 = CartanType(['B',3, 1])
            sage: ct3 = CartanType(['A',3])
            sage: ct1 == ct1
            True
            sage: ct1 == ct2
            False
            sage: ct1 == ct3
            False
        """
        assert (letter in ['A', 'B', 'C', 'BC', 'D', 'E', 'F', 'G'])
        self.letter = letter
        self.n = n
        self.affine = affine

    def _repr_(self, compact=False):
        """
        TESTS::

            sage: ct = CartanType(['A',3, 1])
            sage: repr(ct)
            "['A', 3, 1]"
            sage: ct._repr_(compact=True)
            'A3~'
        """
        letter = self.letter
        n = self.n
        aff = self.affine
        if self.options('notation') == "Kac":
            if letter == 'BC':
                letter = 'A'
                n *= 2
            if compact:
                return '%s%s^%s' % (letter, n, aff)
        if compact:
            return '%s%s~' % (letter, n)
        else:
            return "['%s', %s, %s]" % (letter, n, aff)

    def __reduce__(self):
        """
        TESTS::

            sage: T = CartanType(['D', 4, 1])
            sage: T.__reduce__()
            (CartanType, ('D', 4, 1))
            sage: T == loads(dumps(T))
            True
        """
        return (CartanType, (self.letter, self.n, self.affine))

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['A', 3, 1])
            sage: t[0]
            'A'
            sage: t[1]
            3
            sage: t[2]
            1
            sage: t[3]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i == 1:
            return self.n
        elif i == 2:
            return self.affine
        else:
            raise IndexError("index out of range")

    def rank(self):
        """
        Return the rank of ``self`` which for type `X_n^{(1)}` is `n + 1`.

        EXAMPLES::

            sage: CartanType(['A', 4, 1]).rank()
            5
            sage: CartanType(['B', 4, 1]).rank()
            5
            sage: CartanType(['C', 3, 1]).rank()
            4
            sage: CartanType(['D', 4, 1]).rank()
            5
            sage: CartanType(['E', 6, 1]).rank()
            7
            sage: CartanType(['E', 7, 1]).rank()
            8
            sage: CartanType(['F', 4, 1]).rank()
            5
            sage: CartanType(['G', 2, 1]).rank()
            3
            sage: CartanType(['A', 2, 2]).rank()
            2
            sage: CartanType(['A', 6, 2]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            5
            sage: CartanType(['D', 5, 2]).rank()
            5
            sage: CartanType(['E', 6, 2]).rank()
            5
            sage: CartanType(['D', 4, 3]).rank()
            3
        """
        return self.n+1

    def index_set(self):
        r"""
        Implement :meth:`CartanType_abstract.index_set`.

        The index set for all standard affine Cartan types is of the form
        `\{0, \ldots, n\}`.

        EXAMPLES::

            sage: CartanType(['A', 5, 1]).index_set()
            (0, 1, 2, 3, 4, 5)
        """
        return tuple(range(self.n+1))

    def special_node(self):
        r"""
        Implement :meth:`CartanType_abstract.special_node`.

        With the standard labelling conventions, `0` is always a
        special node.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).special_node()
            0
        """
        return 0

    def type(self):
        """
        Return the type of ``self``.

        EXAMPLES::

            sage: CartanType(['A', 4, 1]).type()
            'A'
        """
        return self.letter

##########################################################################


class CartanType_standard_untwisted_affine(CartanType_standard_affine):
    r"""
    A concrete class for the standard untwisted affine Cartan types.
    """

    def classical(self):
        r"""
        Return the classical Cartan type associated with ``self``.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).classical()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).classical()
            ['B', 3]
            sage: CartanType(['C', 3, 1]).classical()
            ['C', 3]
            sage: CartanType(['D', 4, 1]).classical()
            ['D', 4]
            sage: CartanType(['E', 6, 1]).classical()
            ['E', 6]
            sage: CartanType(['F', 4, 1]).classical()
            ['F', 4]
            sage: CartanType(['G', 2, 1]).classical()
            ['G', 2]
        """
        return CartanType([self.letter,self.n])

    def basic_untwisted(self):
        r"""
        Return the basic_untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic_untwisted type is `X_n`. In
        other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: CartanType(['A', 1, 1]).basic_untwisted()
            ['A', 1]
            sage: CartanType(['A', 3, 1]).basic_untwisted()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).basic_untwisted()
            ['B', 3]
            sage: CartanType(['E', 6, 1]).basic_untwisted()
            ['E', 6]
            sage: CartanType(['G', 2, 1]).basic_untwisted()
            ['G', 2]
        """
        return self.classical()

    def is_untwisted_affine(self):
        """
        Implement :meth:`CartanType_affine.is_untwisted_affine` by
        returning ``True``.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).is_untwisted_affine()
            True
        """
        return True

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['B',4,1]))
            B_{4}^{(1)}
            sage: latex(CartanType(['C',4,1]))
            C_{4}^{(1)}
            sage: latex(CartanType(['D',4,1]))
            D_{4}^{(1)}
            sage: latex(CartanType(['F',4,1]))
            F_4^{(1)}
            sage: latex(CartanType(['G',2,1]))
            G_2^{(1)}
        """
        return self.classical()._latex_()+"^{(1)}"

##########################################################################


class CartanType_decorator(UniqueRepresentation, SageObject, CartanType_abstract):
    """
    Concrete base class for Cartan types that decorate another Cartan type.
    """

    def __init__(self, ct):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: TestSuite(ct).run()
        """
        self._type = ct

    def is_irreducible(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.is_irreducible()
            True
        """
        return self._type.is_irreducible()

    def is_finite(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.is_finite()
            True
        """
        return self._type.is_finite()

    def is_crystallographic(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.is_crystallographic()
            True
        """
        return self._type.is_crystallographic()

    def is_affine(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.is_affine()
            False
        """
        return self._type.is_affine()

    def rank(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G', 2]).relabel({1:2,2:1})
            sage: ct.rank()
            2
        """
        return self._type.rank()

    def index_set(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4, 1]).dual()
           sage: ct.index_set()
           (0, 1, 2, 3, 4)
        """
        return self._type.index_set()

##############################################################################
# Base concrete class for superalgebras


class SuperCartanType_standard(UniqueRepresentation, SageObject):
    # Technical methods
    def _repr_(self, compact=False):
        """
        TESTS::

            sage: ct = CartanType(['A', [3,2]])
            sage: repr(ct)
            "['A', [3, 2]]"
            sage: ct._repr_(compact=True)
            'A3|2'
        """
        formatstr = '%s%s|%s' if compact else "['%s', [%s, %s]]"
        return formatstr % (self.letter, self.m, self.n)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',[4,3]]))
            2
        """
        return 2

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['A', [3,6]])
            sage: t[0]
            'A'
            sage: t[1]
            [3, 6]
            sage: t[2]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i == 1:
            return [self.m, self.n]
        else:
            raise IndexError("index out of range")

    options = CartanType.options
