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

* Within a Sage session, type ``matroids.`` (Do not press "Enter", and do not
  forget the final period ".")
* Hit "tab".

You will see a list of methods which will construct matroids. For example::

   sage: M = matroids.Wheel(4)
   sage: M.is_connected()
   True

or::

   sage: U36 = matroids.Uniform(3, 6)
   sage: U36.equals(U36.dual())
   True

A number of special matroids are collected under a ``named_matroids`` submenu.
To see which, type ``matroids.named_matroids.<tab>`` as above::

    sage: F7 = matroids.named_matroids.Fano()
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
   sage: M.is_isomorphic(matroids.named_matroids.Fano())
   True

   sage: M = Matroid(graphs.PetersenGraph())
   sage: M.rank()
   9

AUTHORS:

- Rudi Pendavingh, Michael Welsh, Stefan van Zwam (2013-04-01): initial
  version

Methods
=======
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from itertools import combinations
import sage.matrix.matrix_space as matrix_space
from sage.matrix.constructor import Matrix
from sage.graphs.all import Graph, graphs
import sage.matrix.matrix
from sage.rings.all import ZZ, QQ, FiniteField, GF
import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from .minor_matroid import MinorMatroid
from .dual_matroid import DualMatroid
from .rank_matroid import RankMatroid
from .circuit_closures_matroid import CircuitClosuresMatroid
from .basis_matroid import BasisMatroid
from .linear_matroid import LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
import sage.matroids.utilities
from networkx import NetworkXError


def Matroid(*args, **kwds):
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

    * Within a Sage session, type "matroids." (Do not press "Enter", and do
      not forget the final period ".")
    * Hit "tab".

    You will see a list of methods which will construct matroids. For
    example::

        sage: F7 = matroids.named_matroids.Fano()
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

    - ``groundset`` -- If provided, the groundset of the matroid.
      If not provided, the function attempts to determine a groundset from the
      data.
    - ``bases`` -- The list of bases (maximal independent sets) of the
      matroid.
    - ``independent_sets`` -- The list of independent sets of the matroid.
    - ``circuits`` -- The list of circuits of the matroid.
    - ``graph`` -- A graph, whose edges form the elements of the matroid.
    - ``matrix`` -- A matrix representation of the matroid.
    - ``reduced_matrix`` -- A reduced representation of the matroid: if
      ``reduced_matrix = A``
      then the matroid is represented by `[I\ \ A]` where `I` is an
      appropriately sized identity matrix.
    - ``rank_function`` -- A function that computes the rank of each subset.
      Can only be provided together with a groundset.
    - ``circuit_closures`` -- Either a list of tuples ``(k, C)`` with ``C``
      the closure of a circuit, and ``k`` the rank of ``C``, or a dictionary
      ``D`` with ``D[k]`` the set of closures of rank-``k`` circuits.
    - ``matroid`` -- An object that is already a matroid. Useful only with the
      ``regular`` option.

    Up to two unnamed arguments are allowed.

    - One unnamed argument, no named arguments other than ``regular`` -- the
      input should be either a graph, or a matrix, or a list of independent
      sets containing all bases, or a matroid.
    - Two unnamed arguments: the first is the groundset, the second a graph,
      or a matrix, or a list of independent sets containing all bases, or a
      matroid.
    - One unnamed argument, at least one named argument: the unnamed argument
      is the groundset, the named argument is as above (but must be different
      from ``groundset``).

    The examples section details how each of the input types deals with
    explicit or implicit groundset arguments.

    OPTIONS:

    - ``regular`` -- (default: ``False``) boolean. If ``True``,
      output a
      :class:`RegularMatroid <sage.matroids.linear_matroid.RegularMatroid>`
      instance such that, *if* the input defines a valid regular matroid, then
      the output represents this matroid. Note that this option can be
      combined with any type of input.
    - ``ring`` -- any ring. If provided, and the input is a ``matrix`` or
      ``reduced_matrix``, output will be a linear matroid over the ring or
      field ``ring``.
    - ``field`` -- any field. Same as ``ring``, but only fields are allowed.
    - ``check`` -- (default: ``True``) boolean. If ``True`` and
      ``regular`` is true, the output is checked to make sure it is a valid
      regular matroid.

    .. WARNING::

        Except for regular matroids, the input is not checked for validity. If
        your data does not correspond to an actual matroid, the behavior of
        the methods is undefined and may cause strange errors. To ensure you
        have a matroid, run
        :meth:`M.is_valid() <sage.matroids.matroid.Matroid.is_valid>`.

    .. NOTE::

        The ``Matroid()`` method will return instances of type
        :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`,
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

    #. List of bases:

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

    #. List of independent sets:

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

    #. List of circuits:

        ::

            sage: M1 = Matroid(groundset='abc', circuits=['bc'])
            sage: M2 = Matroid(bases=['ab', 'ac'])
            sage: M1 == M2
            True

        A matroid specified by a list of circuits gets converted to a
        :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`
        internally::

            sage: M = Matroid(groundset='abcd', circuits=['abc', 'abd', 'acd',
            ....:                                         'bcd'])
            sage: type(M)
            <type 'sage.matroids.basis_matroid.BasisMatroid'>

        Strange things can happen if the input does not satisfy the circuit
        axioms, and these are not always caught by the
        :meth:`is_valid() <sage.matroids.matroid.Matroid.is_valid>` method. So
        always check whether your input makes sense!

        ::

            sage: M = Matroid('abcd', circuits=['ab', 'acd'])
            sage: M.is_valid()
            True
            sage: [sorted(C) for C in M.circuits()]
            [['a']]



    #. Graph:

        Sage has great support for graphs, see :mod:`sage.graphs.graph`.

        ::

            sage: G = graphs.PetersenGraph()
            sage: Matroid(G)
            Regular matroid of rank 9 on 15 elements with 2000 bases


        Note: if a groundset is specified, we assume it is in the same order
        as
        :meth:`G.edge_iterator() <sage.graphs.generic_graph.GenericGraph.edge_iterator>`
        provides::

            sage: G = Graph([(0, 1), (0, 2), (0, 2), (1, 2)],multiedges=True)
            sage: M = Matroid('abcd', G)
            sage: M.rank(['b', 'c'])
            1

        If no groundset is provided, we attempt to use the edge labels::

            sage: G = Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')])
            sage: M = Matroid(G)
            sage: sorted(M.groundset())
            ['a', 'b', 'c']

        If no edge labels are present and the graph is simple, we use the
        tuples ``(i, j)`` of endpoints. If that fails, we simply use a list
        ``[0..m-1]`` ::

            sage: G = Graph([(0, 1), (0, 2), (1, 2)])
            sage: M = Matroid(G)
            sage: sorted(M.groundset())
            [(0, 1), (0, 2), (1, 2)]

            sage: G = Graph([(0, 1), (0, 2), (0, 2), (1, 2)],multiedges=True)
            sage: M = Matroid(G)
            sage: sorted(M.groundset())
            [0, 1, 2, 3]

        When the ``graph`` keyword is used, a variety of inputs can be
        converted to a graph automatically. The following uses a graph6 string
        (see the :class:`Graph <sage.graphs.graph.Graph>` method's
        documentation)::

            sage: Matroid(graph=':I`AKGsaOs`cI]Gb~')
            Regular matroid of rank 9 on 17 elements with 4004 bases

        However, this method is no more clever than ``Graph()``::

            sage: Matroid(graph=41/2)
            Traceback (most recent call last):
            ...
            ValueError: input does not seem to represent a graph.


    #. Matrix:

        The basic input is a
        :mod:`Sage matrix <sage.matrix.constructor>`::

            sage: A = Matrix(GF(2), [[1, 0, 0, 1, 1, 0],
            ....:                    [0, 1, 0, 1, 0, 1],
            ....:                    [0, 0, 1, 0, 1, 1]])
            sage: M = Matroid(matrix=A)
            sage: M.is_isomorphic(matroids.CompleteGraphic(4))
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
            ValueError: basis has wrong cardinality.


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
            sage: Matroid([0, 1, 2, 3, 4, 5],
            ....:         matrix=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
            ....:         field=GF(4, 'x'))
            Quaternary matroid of rank 3 on 6 elements
            sage: Matroid([0, 1, 2, 3, 4, 5],
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

    #. Rank function:

        Any function mapping subsets to integers can be used as input::

            sage: def f(X):
            ....:     return min(len(X), 2)
            ....:
            sage: M = Matroid('abcd', rank_function=f)
            sage: M
            Matroid of rank 2 on 4 elements
            sage: M.is_isomorphic(matroids.Uniform(2, 4))
            True

    #. Circuit closures:

        This is often a really concise way to specify a matroid. The usual way
        is a dictionary of lists::

            sage: M = Matroid(circuit_closures={3: ['edfg', 'acdg', 'bcfg',
            ....:     'cefh', 'afgh', 'abce', 'abdf', 'begh', 'bcdh', 'adeh'],
            ....:     4: ['abcdefgh']})
            sage: M.equals(matroids.named_matroids.P8())
            True

        You can also input tuples `(k, X)` where `X` is the closure of a
        circuit, and `k` the rank of `X`::

            sage: M = Matroid(circuit_closures=[(2, 'abd'), (3, 'abcdef'),
            ....:                               (2, 'bce')])
            sage: M.equals(matroids.named_matroids.Q6())
            True

    #. Matroid:

        Most of the time, the matroid itself is returned::

            sage: M = matroids.named_matroids.Fano()
            sage: N = Matroid(M)
            sage: N is M
            True

        But it can be useful with the ``regular`` option::

            sage: M = Matroid(circuit_closures={2:['adb', 'bec', 'cfa',
            ....:                                  'def'], 3:['abcdef']})
            sage: N = Matroid(M, regular=True)
            sage: N
            Regular matroid of rank 3 on 6 elements with 16 bases
            sage: Matrix(N)
            [1 0 0 1 1 0]
            [0 1 0 1 1 1]
            [0 0 1 0 1 1]

    The ``regular`` option:

        ::

            sage: M = Matroid(reduced_matrix=[[1, 1, 0],
            ....:                             [1, 0, 1],
            ....:                             [0, 1, 1]], regular=True)
            sage: M
            Regular matroid of rank 3 on 6 elements with 16 bases

            sage: M.is_isomorphic(matroids.CompleteGraphic(4))
            True

        By default we check if the resulting matroid is actually regular. To
        increase speed, this check can be skipped::

            sage: M = matroids.named_matroids.Fano()
            sage: N = Matroid(M, regular=True)
            Traceback (most recent call last):
            ...
            ValueError: input does not correspond to a valid regular matroid.
            sage: N = Matroid(M, regular=True, check=False)
            sage: N
            Regular matroid of rank 3 on 7 elements with 32 bases

            sage: N.is_valid()
            False

        Sometimes the output is regular, but represents a different matroid
        from the one you intended::

            sage: M = Matroid(Matrix(GF(3), [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: N = Matroid(Matrix(GF(3), [[1, 0, 1, 1], [0, 1, 1, 2]]),
            ....:             regular=True)
            sage: N.is_valid()
            True
            sage: N.is_isomorphic(M)
            False

    """
    # These are the valid arguments:
    inputS = set(['groundset', 'bases', 'independent_sets', 'circuits', 'graph', 'matrix', 'reduced_matrix', 'rank_function', 'circuit_closures', 'matroid'])

    # process options
    if 'regular' in kwds:
        want_regular = kwds['regular']
        kwds.pop('regular')
    else:
        want_regular = False
    if 'check' in kwds:
        want_check = kwds['check']
        kwds.pop('check')
    else:
        want_check = True

    base_ring = None
    have_field = False
    if 'field' in kwds:
        base_ring = kwds['field']
        kwds.pop('field')
        have_field = True
        if 'ring' in kwds:
            raise ValueError("only one of ring and field can be specified.")
        try:
            if not base_ring.is_field():
                raise TypeError("specified ``field`` is not a field.")
        except AttributeError:
            raise TypeError("specified ``field`` is not a field.")
    if 'ring' in kwds:
        base_ring = kwds['ring']
        kwds.pop('ring')
        try:
            if not base_ring.is_ring():
                raise TypeError("specified ``ring`` is not a ring.")
        except AttributeError:
            raise TypeError("specified ``ring`` is not a ring.")

    # Process unnamed arguments
    if len(args) > 0:
        if 'groundset' in kwds:
            raise ValueError('when using unnamed arguments, groundset must be the first unnamed argument or be implicit.')
        if len(args) > 2:
            raise ValueError('at most two unnamed arguments are allowed.')
        if len(args) == 2 or len(kwds) > 0:
            # First argument should be the groundset
            kwds['groundset'] = args[0]
        # Check for unnamed data
        dataindex = -1
        if len(args) == 2:
            dataindex = 1
        if len(args) == 1 and len(kwds) == 0:
            dataindex = 0
        if dataindex > -1:
            # One unnamed argument, no named arguments
            if isinstance(args[dataindex], sage.graphs.graph.Graph):
                kwds['graph'] = args[dataindex]
            elif isinstance(args[dataindex], sage.matrix.matrix.Matrix):
                kwds['matrix'] = args[dataindex]
            elif isinstance(args[dataindex], sage.matroids.matroid.Matroid):
                kwds['matroid'] = args[dataindex]
            else:
                kwds['independent_sets'] = args[dataindex]

    # Check for multiple types of input
    if len(set(kwds).difference(inputS)) > 0:
        raise ValueError("unknown input argument")
    if ('grondset' in kwds and len(kwds) != 2) or ('groundset' not in kwds and len(kwds) > 1):
        raise ValueError("only one type of input may be specified.")

    # Bases:
    if 'bases' in kwds:
        if 'groundset' not in kwds:
            gs = set()
            for B in kwds['bases']:
                gs.update(B)
            kwds['groundset'] = gs
        M = BasisMatroid(groundset=kwds['groundset'], bases=kwds['bases'])

    # Independent sets:
    if 'independent_sets' in kwds:
        # Convert to list of bases first
        rk = -1
        bases = []
        for I in kwds['independent_sets']:
            if len(I) == rk:
                bases.append(I)
            elif len(I) > rk:
                bases = [I]
                rk = len(I)
        if 'groundset' not in kwds:
            gs = set()
            for B in bases:
                gs.update(B)
            kwds['groundset'] = gs
        M = BasisMatroid(groundset=kwds['groundset'], bases=bases)

    # Circuits:
    if 'circuits' in kwds:
        # Convert to list of bases first
        # Determine groundset (note that this cannot detect coloops)
        if 'groundset' not in kwds:
            gs = set()
            for C in kwds['circuits']:
                gs.update(C)
            kwds['groundset'] = gs
        # determine the rank by computing a basis
        B = set(kwds['groundset'])
        for C in kwds['circuits']:
            I = B.intersection(C)
            if len(I) >= len(C):
                B.discard(I.pop())
        rk = len(B)
        # Construct the basis matroid of appropriate rank. Note: slow!
        BB = [frozenset(B) for B in combinations(kwds['groundset'], rk) if not any([frozenset(C).issubset(B) for C in kwds['circuits']])]
        M = BasisMatroid(groundset=kwds['groundset'], bases=BB)

    # Graphs:
    if 'graph' in kwds:
        # Construct the incidence matrix
        # NOTE: we are not using Sage's built-in method because
        # 1) we would need to fix the loops anyway
        # 2) Sage will sort the columns, making it impossible to keep labels!
        G = kwds['graph']
        if not isinstance(G, sage.graphs.generic_graph.GenericGraph):
            try:
                G = Graph(G)
            except (ValueError, TypeError, NetworkXError):
                raise ValueError("input does not seem to represent a graph.")
        V = G.vertices()
        E = G.edges()
        n = G.num_verts()
        m = G.num_edges()
        A = Matrix(ZZ, n, m, 0)
        mm = 0
        for i, j, k in G.edge_iterator():
            A[V.index(i), mm] = -1
            A[V.index(j), mm] += 1  # So loops get 0
            mm += 1
        # Decide on the groundset
        if 'groundset' not in kwds:
            # 1. Attempt to use edge labels.
            sl = G.edge_labels()
            if len(sl) == len(set(sl)):
                kwds['groundset'] = sl
                # 2. If simple, use vertex tuples
            elif not G.has_multiple_edges():
                kwds['groundset'] = [(i, j) for i, j, k in G.edge_iterator()]
            else:
                # 3. Use numbers
                kwds['groundset'] = range(m)
        M = RegularMatroid(matrix=A, groundset=kwds['groundset'])
        want_regular = False  # Save some time, since result is already regular

    # Matrices:
    if 'matrix' in kwds or 'reduced_matrix' in kwds:
        if 'matrix' in kwds:
            A = kwds['matrix']
        if 'reduced_matrix' in kwds:
            A = kwds['reduced_matrix']

        # Fix the representation
        if not isinstance(A, sage.matrix.matrix.Matrix):
            try:
                if base_ring is not None:
                    A = Matrix(base_ring, A)
                else:
                    A = Matrix(A)
            except ValueError:
                raise ValueError("input does not seem to contain a matrix.")

        # Fix the ring
        if base_ring is not None:
            if A.base_ring() != base_ring:
                A = A.change_ring(base_ring)
        elif A.base_ring() == ZZ and not want_regular:  # Usually a rational matrix is intended, we presume.
            A = A.change_ring(QQ)
            base_ring = QQ
        else:
            base_ring = A.base_ring()

        # Determine groundset:
        if 'matrix' in kwds:
            if 'groundset' in kwds:
                if len(kwds['groundset']) == A.nrows() + A.ncols():
                    kwds['reduced_matrix'] = A
                    kwds.pop('matrix')
                else:
                    if len(kwds['groundset']) != A.ncols():
                        raise ValueError("groundset size does not correspond to matrix size.")
            else:
                kwds['groundset'] = range(A.ncols())
        if 'reduced_matrix' in kwds:
            if 'groundset' in kwds:
                if len(kwds['groundset']) != A.nrows() + A.ncols():
                    raise ValueError("groundset size does not correspond to matrix size.")
            else:
                kwds['groundset'] = range(A.nrows() + A.ncols())
        if 'matrix' in kwds:
            if base_ring == GF(2):
                M = BinaryMatroid(groundset=kwds['groundset'], matrix=A)
            elif base_ring == GF(3):
                M = TernaryMatroid(groundset=kwds['groundset'], matrix=A)
            elif base_ring.is_field() and base_ring.order() == 4:  # GF(4) can have different generators.
                M = QuaternaryMatroid(groundset=kwds['groundset'], matrix=A)
            else:
                M = LinearMatroid(groundset=kwds['groundset'], matrix=A, ring=base_ring)

        if 'reduced_matrix' in kwds:
            if A.base_ring() == GF(2):
                M = BinaryMatroid(groundset=kwds['groundset'], reduced_matrix=A)
            elif A.base_ring() == GF(3):
                M = TernaryMatroid(groundset=kwds['groundset'], reduced_matrix=A)
            elif A.base_ring().is_field() and A.base_ring().order() == 4:  # GF(4) can have different generators.
                M = QuaternaryMatroid(groundset=kwds['groundset'], reduced_matrix=A)
            else:
                M = LinearMatroid(groundset=kwds['groundset'], reduced_matrix=A, ring=base_ring)

    # Rank functions:
    if 'rank_function' in kwds:
        if 'groundset' not in kwds:
            raise ValueError('for rank functions, groundset needs to be specified.')
        M = RankMatroid(groundset=kwds['groundset'], rank_function=kwds['rank_function'])

    # Circuit closures:
    if 'circuit_closures' in kwds:
        if 'groundset' not in kwds:
            E = set()
            if isinstance(kwds['circuit_closures'], dict):
                for X in kwds['circuit_closures'].itervalues():
                    for Y in X:
                        E.update(Y)
            else:
                for X in kwds['circuit_closures']:
                    E.update(X[1])
        else:
            E = kwds['groundset']
        if not isinstance(kwds['circuit_closures'], dict):
            # Convert to dictionary
            CC = {}
            for X in kwds['circuit_closures']:
                if X[0] not in CC:
                    CC[X[0]] = []
                CC[X[0]].append(X[1])
        else:
            CC = kwds['circuit_closures']
        M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)

    # Matroids:
    if 'matroid' in kwds:
        M = kwds['matroid']
        if not isinstance(M, sage.matroids.matroid.Matroid):
            raise ValueError("input does not appear to be of Matroid type.")
    # Regular option:
    if want_regular:
        M = sage.matroids.utilities.make_regular_matroid_from_matroid(M)
        if want_check:
            if not M.is_valid():
                    raise ValueError('input does not correspond to a valid regular matroid.')
    return M
