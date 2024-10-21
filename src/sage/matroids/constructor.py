# sage_setup: distribution = sagemath-modules
r"""
Matroid construction

Theory
======

Matroids are combinatorial structures that capture the abstract properties
of (linear/algebraic/...) dependence. Formally, a matroid is a pair
`M = (E, I)` of a finite set `E`, the *groundset*, and a collection of
subsets `I`, the independent sets, subject to the following axioms:

* `I` contains the empty set
* If `X` is a set in `I`, then each subset of `X` is in `I`
* If two subsets `X`, `Y` are in `I`, and `|X| > |Y|`, then there exists
  `x \in X - Y` such that `Y + \{x\}` is in `I`.

See the :wikipedia:`Wikipedia article on matroids <Matroid>` for more theory
and examples. Matroids can be obtained from many types of mathematical
structures, and Sage supports a number of them.

There are two main entry points to Sage's matroid functionality. The object
:class:`matroids. <sage.matroids.matroids_catalog>` contains a number of
constructors for well-known matroids. The function
:func:`Matroid() <sage.matroids.constructor.Matroid>` allows you to define
your own matroids from a variety of sources. We briefly introduce both below;
follow the links for more comprehensive documentation.

Each matroid object in Sage comes with a number of built-in operations. An
overview can be found in the documentation of
:mod:`the abstract matroid class <sage.matroids.matroid>`.

Built-in matroids
=================

For built-in matroids, do the following:

* Within a Sage session, type ``matroids.`` (Do not press :kbd:`Enter`,
  and do not forget the final period ".")
* Hit :kbd:`Tab`.

You will see a list of methods which will construct matroids. For example::

   sage: M = matroids.Wheel(4)
   sage: M.is_connected()
   True

or::

   sage: U36 = matroids.Uniform(3, 6)
   sage: U36.equals(U36.dual())
   True

A number of special matroids are collected under a ``catalog`` submenu.
To see which, type ``matroids.catalog.<tab>`` as above::

    sage: F7 = matroids.catalog.Fano()
    sage: len(F7.nonspanning_circuits())
    7

Constructing matroids
=====================

To define your own matroid, use the function
:func:`Matroid() <sage.matroids.constructor.Matroid>`. This function attempts
to interpret its arguments to create an appropriate matroid. The input
arguments are documented in detail
:func:`below <sage.matroids.constructor.Matroid>`.

EXAMPLES::

   sage: A = Matrix(GF(2), [[1, 0, 0, 0, 1, 1, 1],
   ....:                    [0, 1, 0, 1, 0, 1, 1],
   ....:                    [0, 0, 1, 1, 1, 0, 1]])
   sage: M = Matroid(A)
   sage: M.is_isomorphic(matroids.catalog.Fano())
   True

   sage: M = Matroid(graphs.PetersenGraph())                                            # needs sage.graphs
   sage: M.rank()                                                                       # needs sage.graphs
   9

AUTHORS:

- Rudi Pendavingh, Michael Welsh, Stefan van Zwam (2013-04-01): initial
  version

Functions
=========
"""

# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from itertools import combinations
from sage.combinat.posets.lattices import FiniteLatticePoset
from sage.matrix.constructor import matrix
from sage.structure.element import Matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.rings.finite_rings.finite_field_base import FiniteField
import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from sage.matroids.rank_matroid import RankMatroid
from sage.matroids.circuits_matroid import CircuitsMatroid
from sage.matroids.flats_matroid import FlatsMatroid
from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
from sage.matroids.basis_matroid import BasisMatroid
from sage.matroids.linear_matroid import LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from sage.matroids.graphic_matroid import GraphicMatroid
import sage.matroids.utilities


def Matroid(groundset=None, data=None, **kwds):
    r"""
    Construct a matroid.

    Matroids are combinatorial structures that capture the abstract properties
    of (linear/algebraic/...) dependence. Formally, a matroid is a pair
    `M = (E, I)` of a finite set `E`, the *groundset*, and a collection of
    subsets `I`, the independent sets, subject to the following axioms:

    * `I` contains the empty set
    * If `X` is a set in `I`, then each subset of `X` is in `I`
    * If two subsets `X`, `Y` are in `I`, and `|X| > |Y|`, then there exists
      `x \in X - Y` such that `Y + \{x\}` is in `I`.

    See the :wikipedia:`Wikipedia article on matroids <Matroid>` for more
    theory and examples. Matroids can be obtained from many types of
    mathematical structures, and Sage supports a number of them.

    There are two main entry points to Sage's matroid functionality. For
    built-in matroids, do the following:

    * Within a Sage session, type "matroids." (Do not press :kbd:`Enter`, and do
      not forget the final period ".")
    * Hit :kbd:`Tab`.

    You will see a list of methods which will construct matroids. For
    example::

        sage: F7 = matroids.catalog.Fano()
        sage: len(F7.nonspanning_circuits())
        7

    or::

        sage: U36 = matroids.Uniform(3, 6)
        sage: U36.equals(U36.dual())
        True

    To define your own matroid, use the function ``Matroid()``.
    This function attempts to interpret its arguments to create an appropriate
    matroid. The following named arguments are supported:

    INPUT:

    - ``groundset`` -- (optional) the groundset of the matroid; if not
      provided, the function attempts to determine a groundset from the data

    Exactly one of the following inputs must be given (where ``data``
    must be a positional argument and anything else must be a keyword
    argument):

    - ``data`` -- a graph or a matrix or a RevLex-Index string or a list
      of independent sets containing all bases or a matroid
    - ``bases`` -- the list of bases (maximal independent sets) of the
      matroid
    - ``independent_sets`` -- the list of independent sets of the matroid
    - ``circuits`` -- the list of circuits of the matroid
    - ``nonspanning_circuits`` -- the list of nonspanning circuits of the
      matroid
    - ``flats`` -- the dictionary, list, or lattice of flats of the matroid
    - ``graph`` -- a graph, whose edges form the elements of the matroid
    - ``matrix`` -- a matrix representation of the matroid
    - ``reduced_matrix`` -- a reduced representation of the matroid: if
      ``reduced_matrix = A``
      then the matroid is represented by `[I\ \ A]` where `I` is an
      appropriately sized identity matrix
    - ``morphism`` -- a morphism representation of the matroid
    - ``reduced_morphism`` -- a reduced morphism representation of the matroid
    - ``rank_function`` -- a function that computes the rank of each subset;
      can only be provided together with a groundset
    - ``circuit_closures`` -- either a list of tuples ``(k, C)`` with ``C``
      the closure of a circuit, and ``k`` the rank of ``C``, or a dictionary
      ``D`` with ``D[k]`` the set of closures of rank-``k`` circuits
    - ``revlex`` -- the encoding as a string of ``0`` and ``*`` symbols;
      used by [Mat2012]_ and explained in [MMIB2012]_
    - ``matroid`` -- an object that is already a matroid; useful only with the
      ``regular`` option

    Further options:

    - ``regular`` -- boolean (default: ``False``); if ``True``,
      output a
      :class:`RegularMatroid <sage.matroids.linear_matroid.RegularMatroid>`
      instance such that, *if* the input defines a valid regular matroid, then
      the output represents this matroid. Note that this option can be
      combined with any type of input.
    - ``ring`` -- any ring. If provided, and the input is a ``matrix`` or
      ``reduced_matrix``, output will be a linear matroid over the ring or
      field ``ring``.
    - ``field`` -- any field. Same as ``ring``, but only fields are allowed
    - ``check`` -- boolean (default: ``True``); if ``True`` and
      ``regular`` is ``True``, the output is checked to make sure it is a valid
      regular matroid

    .. WARNING::

        Except for regular matroids, the input is not checked for validity. If
        your data does not correspond to an actual matroid, the behavior of
        the methods is undefined and may cause strange errors. To ensure you
        have a matroid, run
        :meth:`M.is_valid() <sage.matroids.matroid.Matroid.is_valid>`.

    .. NOTE::

        The ``Matroid()`` method will return instances of type
        :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`,
        :class:`CircuitsMatroid <sage.matroids.circuits_matroid.CircuitsMatroid>`,
        :class:`FlatsMatroid <sage.matroids.flats_matroid.FlatsMatroid>`,
        :class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`,
        :class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`,
        :class:`BinaryMatroid <sage.matroids.linear_matroid.LinearMatroid>`,
        :class:`TernaryMatroid <sage.matroids.linear_matroid.LinearMatroid>`,
        :class:`QuaternaryMatroid <sage.matroids.linear_matroid.LinearMatroid>`,
        :class:`RegularMatroid <sage.matroids.linear_matroid.LinearMatroid>`, or
        :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>`. To
        import these classes (and other useful functions) directly into Sage's
        main namespace, type::

            sage: from sage.matroids.advanced import *

        See :mod:`sage.matroids.advanced <sage.matroids.advanced>`.

    EXAMPLES:

    Note that in these examples we will often use the fact that strings are
    iterable in these examples. So we type ``'abcd'`` to denote the list
    ``['a', 'b', 'c', 'd']``.

    #.  List of bases:

        All of the following inputs are allowed, and equivalent::

            sage: M1 = Matroid(groundset='abcd', bases=['ab', 'ac', 'ad',
            ....:                                       'bc', 'bd', 'cd'])
            sage: M2 = Matroid(bases=['ab', 'ac', 'ad', 'bc', 'bd', 'cd'])
            sage: M3 = Matroid(['ab', 'ac', 'ad', 'bc', 'bd', 'cd'])
            sage: M4 = Matroid('abcd', ['ab', 'ac', 'ad', 'bc', 'bd', 'cd'])
            sage: M5 = Matroid('abcd', bases=[['a', 'b'], ['a', 'c'],
            ....:                             ['a', 'd'], ['b', 'c'],
            ....:                             ['b', 'd'], ['c', 'd']])
            sage: M1 == M2
            True
            sage: M1 == M3
            True
            sage: M1 == M4
            True
            sage: M1 == M5
            True

        We do not check if the provided input forms an actual matroid::

            sage: M1 = Matroid(groundset='abcd', bases=['ab', 'cd'])
            sage: M1.full_rank()
            2
            sage: M1.is_valid()
            False

        Bases may be repeated::

            sage: M1 = Matroid(['ab', 'ac'])
            sage: M2 = Matroid(['ab', 'ac', 'ab'])
            sage: M1 == M2
            True

    #.  List of independent sets:

        ::

            sage: M1 = Matroid(groundset='abcd',
            ....:              independent_sets=['', 'a', 'b', 'c', 'd', 'ab',
            ....:                               'ac', 'ad', 'bc', 'bd', 'cd'])

        We only require that the list of independent sets contains each basis
        of the matroid; omissions of smaller independent sets and
        repetitions are allowed::

            sage: M1 = Matroid(bases=['ab', 'ac'])
            sage: M2 = Matroid(independent_sets=['a', 'ab', 'b', 'ab', 'a',
            ....:                                'b', 'ac'])
            sage: M1 == M2
            True

    #.  List of circuits:

        ::

            sage: M1 = Matroid(groundset='abc', circuits=['bc'])

        A matroid specified by a list of circuits gets converted to a
        :class:`CircuitsMatroid <sage.matroids.circuits_matroid.CircuitsMatroid>`
        internally::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M2 = CircuitsMatroid(Matroid(bases=['ab', 'ac']))
            sage: M1 == M2
            True

            sage: M = Matroid(groundset='abcd', circuits=['abc', 'abd', 'acd',
            ....:                                         'bcd'])
            sage: type(M)
            <class 'sage.matroids.circuits_matroid.CircuitsMatroid'>

        Strange things can happen if the input does not satisfy the circuit
        axioms, and these can be caught by the
        :meth:`is_valid() <sage.matroids.circuits_matroid.CircuitsMatroid.is_valid>`
        method. So please check that your input makes sense!

        ::

            sage: M = Matroid('abcd', circuits=['ab', 'acd'])
            sage: M.is_valid()
            False

    #.  Flats:

        Given a dictionary of flats indexed by their rank, we get a
        :class:`FlatsMatroid <sage.matroids.flats_matroid.FlatsMatroid>`::

            sage: M = Matroid(flats={0: [''], 1: ['a', 'b'], 2: ['ab']})
            sage: M.is_isomorphic(matroids.Uniform(2, 2)) and M.is_valid()
            True
            sage: type(M)
            <class 'sage.matroids.flats_matroid.FlatsMatroid'>

        If instead we simply provide a list of flats, then the class computes
        and stores the lattice of flats upon definition. This can be
        time-consuming, but after it's done we benefit from some faster methods
        (e.g., :meth:`is_valid() <sage.matroids.flats_matroid.FlatsMatroid.is_valid>`)::

            sage: M = Matroid(flats=['', 'a', 'b', 'ab'])
            sage: for i in range(M.rank() + 1):  # print flats by rank
            ....:     print(f'{i}: {sorted([sorted(F) for F in M.flats(i)], key=str)}')
            0: [[]]
            1: [['a'], ['b']]
            2: [['a', 'b']]
            sage: M.is_valid()
            True
            sage: type(M)
            <class 'sage.matroids.flats_matroid.FlatsMatroid'>

        Finally, we can also directly provide a lattice of flats::

            sage: from sage.combinat.posets.lattices import LatticePoset
            sage: flats = [frozenset(F) for F in powerset('ab')]
            sage: L_M = LatticePoset((flats, lambda x, y: x < y))
            sage: M = Matroid(L_M)
            sage: M.is_isomorphic(matroids.Uniform(2, 2)) and M.is_valid()
            True
            sage: type(M)
            <class 'sage.matroids.flats_matroid.FlatsMatroid'>

    #.  Graph:

        Sage has great support for graphs, see :mod:`sage.graphs.graph`.

        ::

            sage: G = graphs.PetersenGraph()                                            # needs sage.graphs
            sage: Matroid(G)                                                            # needs sage.graphs
            Graphic matroid of rank 9 on 15 elements

        If each edge has a unique label, then those are used as the ground set
        labels::

            sage: G = Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')])                    # needs sage.graphs
            sage: M = Matroid(G)                                                        # needs sage.graphs
            sage: sorted(M.groundset())                                                 # needs sage.graphs
            ['a', 'b', 'c']

        If there are parallel edges, then integers are used for the ground set.
        If there are no edges in parallel, and is not a complete list of labels,
        or the labels are not unique, then vertex tuples are used::

            sage: # needs sage.graphs
            sage: G = Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'b')])
            sage: M = Matroid(G)
            sage: sorted(M.groundset())
            [(0, 1), (0, 2), (1, 2)]
            sage: H = Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'b'), (1, 2, 'c')],
            ....:           multiedges=True)
            sage: N = Matroid(H)
            sage: sorted(N.groundset())
            [0, 1, 2, 3]

        The GraphicMatroid object forces its graph to be connected. If a
        disconnected graph is used as input, it will connect the components::

            sage: # needs sage.graphs
            sage: G1 = graphs.CycleGraph(3); G2 = graphs.DiamondGraph()
            sage: G = G1.disjoint_union(G2)
            sage: M = Matroid(G); M
            Graphic matroid of rank 5 on 8 elements
            sage: M.graph()
            Looped multi-graph on 6 vertices
            sage: M.graph().is_connected()
            True
            sage: M.is_connected()
            False


        If the keyword ``regular`` is set to ``True``, the output will instead
        be an instance of :class:`RegularMatroid`.

        ::

            sage: G = Graph([(0, 1), (0, 2), (1, 2)])                                   # needs sage.graphs
            sage: M = Matroid(G, regular=True); M                                       # needs sage.graphs
            Regular matroid of rank 2 on 3 elements with 3 bases

        Note: if a groundset is specified, we assume it is in the same order
        as
        :meth:`G.edge_iterator() <sage.graphs.generic_graph.GenericGraph.edge_iterator>`
        provides::

            sage: G = Graph([(0, 1), (0, 2), (0, 2), (1, 2)], multiedges=True)          # needs sage.graphs
            sage: M = Matroid('abcd', G)                                                # needs sage.graphs
            sage: M.rank(['b', 'c'])                                                    # needs sage.graphs
            1

        As before,
        if no edge labels are present and the graph is simple, we use the
        tuples ``(i, j)`` of endpoints. If that fails, we simply use a list
        ``[0..m-1]`` ::

            sage: G = Graph([(0, 1), (0, 2), (1, 2)])                                   # needs sage.graphs
            sage: M = Matroid(G, regular=True)                                          # needs sage.graphs
            sage: sorted(M.groundset())                                                 # needs sage.graphs
            [(0, 1), (0, 2), (1, 2)]

            sage: G = Graph([(0, 1), (0, 2), (0, 2), (1, 2)], multiedges=True)          # needs sage.graphs
            sage: M = Matroid(G, regular=True)                                          # needs sage.graphs
            sage: sorted(M.groundset())                                                 # needs sage.graphs
            [0, 1, 2, 3]

        When the ``graph`` keyword is used, a variety of inputs can be
        converted to a graph automatically. The following uses a graph6 string
        (see the :class:`Graph <sage.graphs.graph.Graph>` method's
        documentation)::

            sage: Matroid(graph=':I`AKGsaOs`cI]Gb~')                                    # needs sage.graphs
            Graphic matroid of rank 9 on 17 elements

        However, this method is no more clever than ``Graph()``::

            sage: Matroid(graph=41/2)                                                   # needs sage.graphs
            Traceback (most recent call last):
            ...
            ValueError: This input cannot be turned into a graph

    #.  Matrix:

        The basic input is a
        :mod:`Sage matrix <sage.matrix.constructor>`::

            sage: A = Matrix(GF(2), [[1, 0, 0, 1, 1, 0],
            ....:                    [0, 1, 0, 1, 0, 1],
            ....:                    [0, 0, 1, 0, 1, 1]])
            sage: M = Matroid(matrix=A)
            sage: M.is_isomorphic(matroids.CompleteGraphic(4))                          # needs sage.graphs
            True

        Various shortcuts are possible::

            sage: M1 = Matroid(matrix=[[1, 0, 0, 1, 1, 0],
            ....:                      [0, 1, 0, 1, 0, 1],
            ....:                      [0, 0, 1, 0, 1, 1]], ring=GF(2))
            sage: M2 = Matroid(reduced_matrix=[[1, 1, 0],
            ....:                              [1, 0, 1],
            ....:                              [0, 1, 1]], ring=GF(2))
            sage: M3 = Matroid(groundset=[0, 1, 2, 3, 4, 5],
            ....:              matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:              ring=GF(2))
            sage: A = Matrix(GF(2), [[1, 1, 0], [1, 0, 1], [0, 1, 1]])
            sage: M4 = Matroid([0, 1, 2, 3, 4, 5], A)
            sage: M1 == M2
            True
            sage: M1 == M3
            True
            sage: M1 == M4
            True

        However, with unnamed arguments the input has to be a ``Matrix``
        instance, or the function will try to interpret it as a set of bases::

            sage: Matroid([0, 1, 2], [[1, 0, 1], [0, 1, 1]])
            Traceback (most recent call last):
            ...
            ValueError: basis has wrong cardinality

        If the groundset size equals number of rows plus number of columns, an
        identity matrix is prepended. Otherwise the groundset size must equal
        the number of columns::

            sage: A = Matrix(GF(2), [[1, 1, 0], [1, 0, 1], [0, 1, 1]])
            sage: M = Matroid([0, 1, 2], A)
            sage: N = Matroid([0, 1, 2, 3, 4, 5], A)
            sage: M.rank()
            2
            sage: N.rank()
            3

        We automatically create an optimized subclass, if available::

            sage: Matroid([0, 1, 2, 3, 4, 5],
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(2))
            Binary matroid of rank 3 on 6 elements, type (2, 7)
            sage: Matroid([0, 1, 2, 3, 4, 5],
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(3))
            Ternary matroid of rank 3 on 6 elements, type 0-
            sage: Matroid([0, 1, 2, 3, 4, 5],                                           # needs sage.rings.finite_rings
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(4, 'x'))
            Quaternary matroid of rank 3 on 6 elements
            sage: Matroid([0, 1, 2, 3, 4, 5],                                           # needs sage.graphs
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(2), regular=True)
            Regular matroid of rank 3 on 6 elements with 16 bases

        Otherwise the generic LinearMatroid class is used::

            sage: Matroid([0, 1, 2, 3, 4, 5],
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(83))
            Linear matroid of rank 3 on 6 elements represented over the Finite
            Field of size 83

        An integer matrix is automatically converted to a matrix over `\QQ`.
        If you really want integers, you can specify the ring explicitly::

            sage: A = Matrix([[1, 1, 0], [1, 0, 1], [0, 1, -1]])
            sage: A.base_ring()
            Integer Ring
            sage: M = Matroid([0, 1, 2, 3, 4, 5], A)
            sage: M.base_ring()
            Rational Field
            sage: M = Matroid([0, 1, 2, 3, 4, 5], A, ring=ZZ)
            sage: M.base_ring()
            Integer Ring

        A morphism representation of a :class:`LinearMatroid` can also be used as
        input::

            sage: M = matroids.catalog.Fano()
            sage: A = M.representation(order=True); A
            Generic morphism:
              From: Free module generated by {'a', 'b', 'c', 'd', 'e', 'f', 'g'} over
                    Finite Field of size 2
              To:   Free module generated by {0, 1, 2} over Finite Field of size 2
            sage: A._unicode_art_matrix()
              a b c d e f g
            0⎛1 0 0 0 1 1 1⎞
            1⎜0 1 0 1 0 1 1⎟
            2⎝0 0 1 1 1 0 1⎠
            sage: N = Matroid(A); N
            Binary matroid of rank 3 on 7 elements, type (3, 0)
            sage: N.groundset()
            frozenset({'a', 'b', 'c', 'd', 'e', 'f', 'g'})
            sage: M == N
            True

        The keywords ``morphism`` and ``reduced_morphism`` are also available::

            sage: M = matroids.catalog.RelaxedNonFano("abcdefg")
            sage: A = M.representation(order=True, reduced=True); A
            Generic morphism:
              From: Free module generated by {'d', 'e', 'f', 'g'} over
                    Finite Field in w of size 2^2
              To:   Free module generated by {'a', 'b', 'c'} over
                    Finite Field in w of size 2^2
            sage: A._unicode_art_matrix()
              d e f g
            a⎛1 1 0 1⎞
            b⎜1 0 1 1⎟
            c⎝0 1 w 1⎠
            sage: N = Matroid(reduced_morphism=A); N
            Quaternary matroid of rank 3 on 7 elements
            sage: N.groundset()
            frozenset({'a', 'b', 'c', 'd', 'e', 'f', 'g'})
            sage: M == N
            True

    #.  Rank function:

        Any function mapping subsets to integers can be used as input::

            sage: def f(X):
            ....:     return min(len(X), 2)
            sage: M = Matroid('abcd', rank_function=f)
            sage: M
            Matroid of rank 2 on 4 elements
            sage: M.is_isomorphic(matroids.Uniform(2, 4))
            True

    #.  Circuit closures:

        This is often a really concise way to specify a matroid. The usual way
        is a dictionary of lists::

            sage: M = Matroid(circuit_closures={3: ['edfg', 'acdg', 'bcfg',
            ....:     'cefh', 'afgh', 'abce', 'abdf', 'begh', 'bcdh', 'adeh'],
            ....:     4: ['abcdefgh']})
            sage: M.equals(matroids.catalog.P8())
            True

        You can also input tuples `(k, X)` where `X` is the closure of a
        circuit, and `k` the rank of `X`::

            sage: M = Matroid(circuit_closures=[(2, 'abd'), (3, 'abcdef'),
            ....:                               (2, 'bce')])
            sage: M.equals(matroids.catalog.Q6())                                       # needs sage.rings.finite_rings
            True

    #.  RevLex-Index:

        This requires the ``groundset`` to be given and also needs a
        additional keyword argument ``rank`` to specify the rank of the
        matroid::

            sage: M = Matroid("abcdef", "000000******0**", rank=4); M
            Matroid of rank 4 on 6 elements with 8 bases
            sage: list(M.bases())
            [frozenset({'a', 'b', 'd', 'f'}),
             frozenset({'a', 'c', 'd', 'f'}),
             frozenset({'b', 'c', 'd', 'f'}),
             frozenset({'a', 'b', 'e', 'f'}),
             frozenset({'a', 'c', 'e', 'f'}),
             frozenset({'b', 'c', 'e', 'f'}),
             frozenset({'b', 'd', 'e', 'f'}),
             frozenset({'c', 'd', 'e', 'f'})]

        Only the ``0`` symbols really matter, any symbol can be used
        instead of ``*``:

            sage: Matroid("abcdefg", revlex='0++++++++0++++0+++++0+--++----+--++', rank=4)
            Matroid of rank 4 on 7 elements with 31 bases

        It is checked that the input makes sense (but not that it
        defines a matroid)::

            sage: Matroid("abcdef", "000000******0**")
            Traceback (most recent call last):
            ...
            TypeError: for RevLex-Index, the rank needs to be specified
            sage: Matroid("abcdef", "000000******0**", rank=3)
            Traceback (most recent call last):
            ...
            ValueError: expected string of length 20 (6 choose 3), got 15
            sage: M = Matroid("abcdef", "*0000000000000*", rank=4); M
            Matroid of rank 4 on 6 elements with 2 bases
            sage: M.is_valid()
            False

    #.  Matroid:

        Most of the time, the matroid itself is returned::

            sage: M = matroids.catalog.Fano()
            sage: N = Matroid(M)
            sage: N is M
            True

        But it can be useful with the ``regular`` option::

            sage: M = Matroid(circuit_closures={2:['adb', 'bec', 'cfa',
            ....:                                  'def'], 3:['abcdef']})
            sage: N = Matroid(M, regular=True); N                                       # needs sage.graphs
            Regular matroid of rank 3 on 6 elements with 16 bases
            sage: M == N                                                                # needs sage.graphs
            False
            sage: M.is_isomorphic(N)                                                    # needs sage.graphs
            True
            sage: Matrix(N)  # random                                                   # needs sage.graphs
            [1 0 0 1 1 0]
            [0 1 0 1 1 1]
            [0 0 1 0 1 1]

    The ``regular`` option::

        sage: M = Matroid(reduced_matrix=[[1, 1, 0],                                    # needs sage.graphs
        ....:                             [1, 0, 1],
        ....:                             [0, 1, 1]], regular=True); M
        Regular matroid of rank 3 on 6 elements with 16 bases

        sage: M.is_isomorphic(matroids.CompleteGraphic(4))                              # needs sage.graphs
        True

    By default we check if the resulting matroid is actually regular. To
    increase speed, this check can be skipped::

        sage: M = matroids.catalog.Fano()
        sage: N = Matroid(M, regular=True)                                              # needs sage.graphs
        Traceback (most recent call last):
        ...
        ValueError: input is not a valid regular matroid
        sage: N = Matroid(M, regular=True, check=False); N                              # needs sage.graphs
        Regular matroid of rank 3 on 7 elements with 32 bases

        sage: N.is_valid()                                                              # needs sage.graphs
        False

    Sometimes the output is regular, but represents a different matroid
    from the one you intended::

        sage: M = Matroid(Matrix(GF(3), [[1, 0, 1, 1], [0, 1, 1, 2]]))
        sage: N = Matroid(Matrix(GF(3), [[1, 0, 1, 1], [0, 1, 1, 2]]),                  # needs sage.graphs
        ....:             regular=True)
        sage: N.is_valid()                                                              # needs sage.graphs
        True
        sage: N.is_isomorphic(M)                                                        # needs sage.graphs
        False

    TESTS::

        sage: Matroid()
        Traceback (most recent call last):
        ...
        TypeError: no input data given for Matroid()
        sage: Matroid("abc", bases=["abc"], foo='bar')
        Traceback (most recent call last):
        ...
        TypeError: ...Matroid() got an unexpected keyword argument 'foo'
        sage: Matroid(data=["x"], matrix=Matrix(1,1))
        Traceback (most recent call last):
        ...
        TypeError: ...Matroid() got an unexpected keyword argument 'matrix'
        sage: Matroid(bases=["x"], matrix=Matrix(1,1))
        Traceback (most recent call last):
        ...
        TypeError: ...Matroid() got an unexpected keyword argument 'matrix'
        sage: Matroid(Matrix(1,1), ring=ZZ, field=QQ)
        Traceback (most recent call last):
        ...
        TypeError: ...Matroid() got an unexpected keyword argument 'ring'
        sage: Matroid(rank_function=lambda X: len(X))
        Traceback (most recent call last):
        ...
        TypeError: for rank functions, the groundset needs to be specified
        sage: Matroid(matroid='rubbish')
        Traceback (most recent call last):
        ...
        TypeError: input 'rubbish' is not a matroid
    """
    # process options
    want_regular = kwds.pop('regular', False)
    check = kwds.pop('check', True)

    base_ring = None
    if 'field' in kwds:
        base_ring = kwds.pop('field')
        if check and base_ring not in Fields():
            raise TypeError("{} is not a field".format(base_ring))
    elif 'ring' in kwds:
        base_ring = kwds.pop('ring')
        if check and base_ring not in Rings():
            raise TypeError("{} is not a ring".format(base_ring))

    # "key" is the kind of data we got
    key = None
    if data is None:
        for k in ['bases', 'independent_sets', 'circuits',
                  'nonspanning_circuits', 'flats', 'graph', 'matrix',
                  'reduced_matrix', 'morphism', 'reduced_morphism',
                  'rank_function', 'revlex', 'circuit_closures', 'matroid']:
            if k in kwds:
                data = kwds.pop(k)
                key = k
                break
        else:
            # Assume that the single positional argument was actually
            # the data (instead of the groundset)
            data = groundset
            groundset = None

    if key is None:
        try:
            from sage.graphs.graph import Graph
        except ImportError:
            Graph = ()
        if isinstance(data, Graph):
            key = 'graph'
        elif isinstance(data, Matrix) or (
             isinstance(data, tuple) and isinstance(data[0], Matrix)):
            key = 'matrix'
        elif isinstance(data, sage.modules.with_basis.morphism.ModuleMorphism) or (
             isinstance(data, tuple) and
             isinstance(data[0], sage.modules.with_basis.morphism.ModuleMorphism)):
            key = 'morphism'
        elif isinstance(data, sage.matroids.matroid.Matroid):
            key = 'matroid'
        elif isinstance(data, str):
            key = 'revlex'
        elif isinstance(data, dict) or isinstance(data, FiniteLatticePoset):
            key = 'flats'
        elif data is None:
            raise TypeError("no input data given for Matroid()")
        else:
            key = 'independent_sets'

    # Bases:
    if key == 'bases':
        if groundset is None:
            groundset = set()
            for B in data:
                groundset.update(B)
        M = BasisMatroid(groundset=groundset, bases=data)

    # Independent sets:
    elif key == 'independent_sets':
        # Convert to list of bases first
        rk = -1
        bases = []
        for I in data:
            if len(I) == rk:
                bases.append(I)
            elif len(I) > rk:
                bases = [I]
                rk = len(I)
        if groundset is None:
            groundset = set()
            for B in bases:
                groundset.update(B)
        M = BasisMatroid(groundset=groundset, bases=bases)

    # Circuits:
    elif key == 'circuits':
        # Convert to list of bases first
        # Determine groundset (note that this cannot detect coloops)
        if groundset is None:
            groundset = set()
            for C in data:
                groundset.update(C)
        M = CircuitsMatroid(groundset=groundset, circuits=data)

    # Nonspanning circuits:
    elif key == 'nonspanning_circuits':
        try:
            rk = kwds.pop("rank")
        except TypeError:
            raise TypeError("the rank needs to be specified alongside the " +
                            "nonspanning circuits")
        # Determine groundset (note that this cannot detect coloops)
        if groundset is None:
            groundset = set()
            for C in data:
                groundset.update(C)
        # Construct the basis matroid of appropriate rank. Note: slow!
        B = []  # bases
        for b in combinations(groundset, rk):
            flag = True
            for C in data:
                if set(b) >= set(C):
                    flag = False
                    break
            if flag:
                B += [list(b)]
        # convert to circuits matroid defined by non-spanning circuits
        M = CircuitsMatroid(
            BasisMatroid(groundset=groundset, bases=B),
            nsc_defined=True
        )

    # Flats
    elif key == 'flats':
        # Determine groundset
        if groundset is None:
            groundset = set()
            if isinstance(data, dict):
                for i in data:
                    for F in data[i]:
                        groundset.update(F)
            else:  # iterable of flats (including lattice)
                for F in data:
                    groundset.update(F)
        M = FlatsMatroid(groundset=groundset, flats=data)

    # Graphs:
    elif key == 'graph':
        from sage.graphs.graph import Graph

        if isinstance(data, sage.graphs.generic_graph.GenericGraph):
            G = data
        else:
            G = Graph(data)
        # Decide on the groundset
        m = G.num_edges()
        if groundset is None:
            # 1. Attempt to use edge labels.
            sl = G.edge_labels()
            if len(sl) == len(set(sl)):
                groundset = sl
                # 2. If simple, use vertex tuples
            elif not G.has_multiple_edges():
                groundset = [(i, j) for i, j, k in G.edge_iterator()]
            else:
                # 3. Use numbers
                groundset = list(range(m))
        if want_regular:
            # Construct the incidence matrix
            # NOTE: we are not using Sage's built-in method because
            # 1) we would need to fix the loops anyway
            # 2) Sage will sort the columns, making it impossible to keep labels!
            V = G.vertices(sort=True)
            n = G.num_verts()
            A = matrix(ZZ, n, m, 0)
            mm = 0
            for i, j, k in G.edge_iterator():
                A[V.index(i), mm] = -1
                A[V.index(j), mm] += 1  # So loops get 0
                mm += 1
            M = RegularMatroid(matrix=A, groundset=groundset)
            want_regular = False  # Save some time, since result is already regular
        else:
            M = GraphicMatroid(G, groundset=groundset)

    # Matrices:
    elif key in ['matrix', 'reduced_matrix', 'morphism', 'reduced_morphism']:
        A = data
        is_reduced = (key == 'reduced_matrix' or key == 'reduced_morphism')
        if isinstance(data, tuple):
            A = data[0]
            if key == 'matrix' or key == 'reduced_matrix':
                if groundset is None:
                    groundset = data[1]
                    if is_reduced:
                        groundset += data[2]
        if key == 'morphism' or key == 'reduced_morphism':
            if groundset is None:
                groundset = list(A.domain().basis().keys())
                if is_reduced:
                    groundset = list(A.codomain().basis().keys()) + groundset
            A = A.matrix()

        # Fix the representation
        if not isinstance(A, Matrix):
            if base_ring is not None:
                A = matrix(base_ring, A)
            else:
                A = matrix(A)

        # Fix the ring
        if base_ring is not None:
            if A.base_ring() is not base_ring:
                A = A.change_ring(base_ring)
        elif A.base_ring() is ZZ and not want_regular:  # Usually a rational matrix is intended, we presume.
            A = A.change_ring(QQ)
            base_ring = QQ
        else:
            base_ring = A.base_ring()

        # Check groundset
        if groundset is not None:
            if not is_reduced:
                if len(groundset) == A.ncols():
                    pass
                elif len(groundset) == A.nrows() + A.ncols():
                    is_reduced = True
                else:
                    raise ValueError("groundset size does not correspond to matrix size")
            elif is_reduced:
                if len(groundset) == A.nrows() + A.ncols():
                    pass
                else:
                    raise ValueError("groundset size does not correspond to matrix size")

        if is_reduced:
            kw = dict(groundset=groundset, reduced_matrix=A)
        else:
            kw = dict(groundset=groundset, matrix=A)

        if isinstance(base_ring, FiniteField):
            q = base_ring.order()
        else:
            q = 0

        if q == 2:
            M = BinaryMatroid(**kw)
        elif q == 3:
            M = TernaryMatroid(**kw)
        elif q == 4:
            M = QuaternaryMatroid(**kw)
        else:
            M = LinearMatroid(ring=base_ring, **kw)

    # Rank functions:
    elif key == 'rank_function':
        if groundset is None:
            raise TypeError('for rank functions, the groundset needs to be specified')
        M = RankMatroid(groundset=groundset, rank_function=data)

    # RevLex-Index:
    elif key == "revlex":
        if groundset is None:
            raise TypeError('for RevLex-Index, the groundset needs to be specified')
        try:
            rk = kwds.pop("rank")
        except KeyError:
            raise TypeError('for RevLex-Index, the rank needs to be specified')

        groundset = tuple(groundset)
        data = tuple(data)
        rk = int(rk)
        N = len(groundset)

        def revlex_sort_key(s):
            return tuple(reversed(s))
        subsets = sorted(combinations(range(N), rk), key=revlex_sort_key)
        if len(data) != len(subsets):
            raise ValueError("expected string of length %s (%s choose %s), got %s" %
                (len(subsets), N, rk, len(data)))
        bases = []
        for i, x in enumerate(data):
            if x != '0':
                bases.append([groundset[c] for c in subsets[i]])
        M = BasisMatroid(groundset=groundset, bases=bases)

    # Circuit closures:
    elif key == 'circuit_closures':
        if isinstance(data, dict):
            CC = data
        else:
            # Convert to dictionary
            CC = {}
            for X in data:
                if X[0] not in CC:
                    CC[X[0]] = []
                CC[X[0]].append(X[1])

        if groundset is None:
            groundset = set()
            for X in CC.values():
                for Y in X:
                    groundset.update(Y)

        M = CircuitClosuresMatroid(groundset=groundset, circuit_closures=CC)

    # Matroids:
    elif key == 'matroid':
        if not isinstance(data, sage.matroids.matroid.Matroid):
            raise TypeError("input {!r} is not a matroid".format(data))
        M = data

    else:
        raise AssertionError("unknown key %r" % key)

    # All keywords should be used
    for k in kwds:
        raise TypeError("Matroid() got an unexpected keyword argument '{}'".format(k))

    if want_regular:
        M = sage.matroids.utilities.make_regular_matroid_from_matroid(M)
        if check and not M.is_valid():
            raise ValueError('input is not a valid regular matroid')

    return M
