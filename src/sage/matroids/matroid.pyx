# -*- coding: utf-8 -*-
r"""
The abstract Matroid class

Matroids are combinatorial structures that capture the abstract properties
of (linear/algebraic/...) dependence. See the :wikipedia:`Wikipedia article on
matroids <Matroid>` for theory and examples. In Sage, various types of
matroids are supported:
:class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`,
:class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`,
:class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`
(and some specialized subclasses),
:class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>`.
To construct them, use the function
:func:`Matroid() <sage.matroids.constructor.Matroid>`.

All these classes share a common interface, which includes the following
methods (organized by category). Note that most subclasses (notably
:mod:`LinearMatroids <sage.matroids.linear_matroid>`) will implement
additional functionality (e.g. linear extensions).

- Ground set:
    - :meth:`groundset() <sage.matroids.matroid.Matroid.groundset>`
    - :meth:`size() <sage.matroids.matroid.Matroid.size>`

- Rank, bases, circuits, closure
    - :meth:`rank() <sage.matroids.matroid.Matroid.rank>`
    - :meth:`full_rank() <sage.matroids.matroid.Matroid.full_rank>`
    - :meth:`basis() <sage.matroids.matroid.Matroid.basis>`
    - :meth:`max_independent() <sage.matroids.matroid.Matroid.max_independent>`
    - :meth:`circuit() <sage.matroids.matroid.Matroid.circuit>`
    - :meth:`fundamental_circuit() <sage.matroids.matroid.Matroid.fundamental_circuit>`
    - :meth:`closure() <sage.matroids.matroid.Matroid.closure>`
    - :meth:`augment() <sage.matroids.matroid.Matroid.augment>`

    - :meth:`corank() <sage.matroids.matroid.Matroid.corank>`
    - :meth:`full_corank() <sage.matroids.matroid.Matroid.full_corank>`
    - :meth:`cobasis() <sage.matroids.matroid.Matroid.cobasis>`
    - :meth:`max_coindependent() <sage.matroids.matroid.Matroid.max_coindependent>`
    - :meth:`cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`
    - :meth:`fundamental_cocircuit() <sage.matroids.matroid.Matroid.fundamental_cocircuit>`
    - :meth:`coclosure() <sage.matroids.matroid.Matroid.coclosure>`

    - :meth:`is_independent() <sage.matroids.matroid.Matroid.is_independent>`
    - :meth:`is_dependent() <sage.matroids.matroid.Matroid.is_dependent>`
    - :meth:`is_basis() <sage.matroids.matroid.Matroid.is_basis>`
    - :meth:`is_circuit() <sage.matroids.matroid.Matroid.is_circuit>`
    - :meth:`is_closed() <sage.matroids.matroid.Matroid.is_closed>`

    - :meth:`is_coindependent() <sage.matroids.matroid.Matroid.is_coindependent>`
    - :meth:`is_codependent() <sage.matroids.matroid.Matroid.is_codependent>`
    - :meth:`is_cobasis() <sage.matroids.matroid.Matroid.is_cobasis>`
    - :meth:`is_cocircuit() <sage.matroids.matroid.Matroid.is_cocircuit>`
    - :meth:`is_coclosed() <sage.matroids.matroid.Matroid.is_coclosed>`

- Verification
    - :meth:`is_valid() <sage.matroids.matroid.Matroid.is_valid>`

- Enumeration
    - :meth:`circuits() <sage.matroids.matroid.Matroid.circuits>`
    - :meth:`nonspanning_circuits() <sage.matroids.matroid.Matroid.nonspanning_circuits>`
    - :meth:`cocircuits() <sage.matroids.matroid.Matroid.cocircuits>`
    - :meth:`noncospanning_cocircuits() <sage.matroids.matroid.Matroid.noncospanning_cocircuits>`
    - :meth:`circuit_closures() <sage.matroids.matroid.Matroid.circuit_closures>`
    - :meth:`nonspanning_circuit_closures() <sage.matroids.matroid.Matroid.nonspanning_circuit_closures>`
    - :meth:`bases() <sage.matroids.matroid.Matroid.bases>`
    - :meth:`independent_r_sets() <sage.matroids.matroid.Matroid.independent_r_sets>`
    - :meth:`nonbases() <sage.matroids.matroid.Matroid.nonbases>`
    - :meth:`dependent_r_sets() <sage.matroids.matroid.Matroid.dependent_r_sets>`
    - :meth:`flats() <sage.matroids.matroid.Matroid.flats>`
    - :meth:`coflats() <sage.matroids.matroid.Matroid.coflats>`
    - :meth:`hyperplanes() <sage.matroids.matroid.Matroid.hyperplanes>`
    - :meth:`f_vector() <sage.matroids.matroid.Matroid.f_vector>`
    - :meth:`broken_circuits() <sage.matroids.matroid.Matroid.broken_circuits>`
    - :meth:`no_broken_circuits_sets() <sage.matroids.matroid.Matroid.no_broken_circuits_sets>`

- Comparison
    - :meth:`is_isomorphic() <sage.matroids.matroid.Matroid.is_isomorphic>`
    - :meth:`equals() <sage.matroids.matroid.Matroid.equals>`
    - :meth:`is_isomorphism() <sage.matroids.matroid.Matroid.is_isomorphism>`

- Minors, duality, truncation
    - :meth:`minor() <sage.matroids.matroid.Matroid.minor>`
    - :meth:`contract() <sage.matroids.matroid.Matroid.contract>`
    - :meth:`delete() <sage.matroids.matroid.Matroid.delete>`
    - :meth:`dual() <sage.matroids.matroid.Matroid.dual>`
    - :meth:`truncation() <sage.matroids.matroid.Matroid.truncation>`
    - :meth:`has_minor() <sage.matroids.matroid.Matroid.has_minor>`
    - :meth:`has_line_minor() <sage.matroids.matroid.Matroid.has_line_minor>`


- Extension
    - :meth:`extension() <sage.matroids.matroid.Matroid.extension>`
    - :meth:`coextension() <sage.matroids.matroid.Matroid.coextension>`
    - :meth:`modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`
    - :meth:`linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`
    - :meth:`extensions() <sage.matroids.matroid.Matroid.extensions>`
    - :meth:`coextensions() <sage.matroids.matroid.Matroid.coextensions>`

- Connectivity, simplicity
    - :meth:`loops() <sage.matroids.matroid.Matroid.loops>`
    - :meth:`coloops() <sage.matroids.matroid.Matroid.coloops>`
    - :meth:`simplify() <sage.matroids.matroid.Matroid.simplify>`
    - :meth:`cosimplify() <sage.matroids.matroid.Matroid.cosimplify>`
    - :meth:`is_simple() <sage.matroids.matroid.Matroid.is_simple>`
    - :meth:`is_cosimple() <sage.matroids.matroid.Matroid.is_cosimple>`
    - :meth:`components() <sage.matroids.matroid.Matroid.components>`
    - :meth:`is_connected() <sage.matroids.matroid.Matroid.is_connected>`
    - :meth:`is_3connected() <sage.matroids.matroid.Matroid.is_3connected>`

- Optimization
    - :meth:`max_weight_independent() <sage.matroids.matroid.Matroid.max_weight_independent>`
    - :meth:`max_weight_coindependent() <sage.matroids.matroid.Matroid.max_weight_coindependent>`
    - :meth:`intersection() <sage.matroids.matroid.Matroid.intersection>`

- Invariants
    - :meth:`tutte_polynomial() <sage.matroids.matroid.Matroid.tutte_polynomial>`
    - :meth:`flat_cover() <sage.matroids.matroid.Matroid.flat_cover>`
    
- Visualization
    - :meth:`show() <sage.matroids.matroid.Matroid.show>`
    - :meth:`plot() <sage.matroids.matroid.Matroid.plot>`

- Misc
    - :meth:`broken_circuit_complex() <sage.matroids.matroid.Matroid.broken_circuit_complex>`
    - :meth:`matroid_polytope() <sage.matroids.matroid.Matroid.matroid_polytope>`
    - :meth:`independence_matroid_polytope() <sage.matroids.matroid.Matroid.independence_matroid_polytope>`

In addition to these, all methods provided by
:class:`SageObject <sage.structure.sage_object.SageObject>` are available,
notably :meth:`save() <sage.structure.sage_object.SageObject.save>` and
:meth:`rename() <sage.structure.sage_object.SageObject.rename>`.

Advanced usage
==============
Many methods (such as ``M.rank()``) have a companion method whose name starts
with an underscore (such as ``M._rank()``). The method with the underscore
does not do any checks on its input. For instance, it may assume of its input
that

- It is a subset of the groundset. The interface is compatible with Python's
  ``frozenset`` type.
- It is a list of things, supports iteration, and recursively these rules
  apply to its members.

Using the underscored version could improve the speed of code a little, but
will generate more cryptic error messages when presented with wrong input.
In some instances, no error might occur and a nonsensical answer returned.

A subclass should always override the underscored method, if available, and as
a rule leave the regular method alone.

These underscored methods are not documented in the reference manual. To see
them, within Sage you can create a matroid ``M`` and type ``M._<tab>``. Then
``M._rank?`` followed by ``<tab>`` will bring up the documentation string of
the ``_rank()`` method.

Creating new Matroid subclasses
===============================
Many mathematical objects give rise to matroids, and not all are available
through the provided code. For incidental use, the
:mod:`RankMatroid <sage.matroids.rank_matroid>` subclass may suffice. If you
regularly use matroids based on a new data type, you can write a subclass of
``Matroid``. You only need to override the ``__init__``, ``_rank()`` and
``groundset()`` methods to get a fully working class.

EXAMPLE:

In a partition matroid, a subset is independent if it has at most one
element from each partition. The following is a very basic implementation,
in which the partition is specified as a list of lists::

    sage: import sage.matroids.matroid
    sage: class PartitionMatroid(sage.matroids.matroid.Matroid):
    ....:     def __init__(self, partition):
    ....:         self.partition = partition
    ....:         E = set()
    ....:         for P in partition:
    ....:             E.update(P)
    ....:         self.E = frozenset(E)
    ....:     def groundset(self):
    ....:         return self.E
    ....:     def _rank(self, X):
    ....:         X2 = set(X)
    ....:         used_indices = set()
    ....:         rk = 0
    ....:         while len(X2) > 0:
    ....:             e = X2.pop()
    ....:             for i in range(len(self.partition)):
    ....:                 if e in self.partition[i]:
    ....:                     if i not in used_indices:
    ....:                         used_indices.add(i)
    ....:                         rk = rk + 1
    ....:                     break
    ....:         return rk
    ....:
    sage: M = PartitionMatroid([[1, 2], [3, 4, 5], [6, 7]])
    sage: M.full_rank()
    3
    sage: M.tutte_polynomial(var('x'), var('y'))
    x^2*y^2 + 2*x*y^3 + y^4 + x^3 + 3*x^2*y + 3*x*y^2 + y^3

.. NOTE::

    The abstract base class has no idea about the data used to represent the
    matroid. Hence some methods need to be customized to function properly.

    Necessary:

    - ``def __init__(self, ...)``
    - ``def groundset(self)``
    - ``def _rank(self, X)``

    Representation:

    - ``def _repr_(self)``

    Comparison:

    - ``def __hash__(self)``
    - ``def __eq__(self, other)``
    - ``def __ne__(self, other)``

    In Cythonized classes, use ``__richcmp__()`` instead of ``__eq__()``,
    ``__ne__()``.

    Copying, loading, saving:

    - ``def __copy__(self)``
    - ``def __deepcopy__(self, memo={})``
    - ``def __reduce__(self)``

    See, for instance, :class:`rank_matroid <sage.matroids.rank_matroid>` or
    :class:`circuit_closures_matroid <sage.matroids.circuit_closures_matroid>`
    for sample implementations of these.

.. NOTE::

    The example provided does not check its input at all. You may want to make
    sure the input data are not corrupt.

Some examples
=============

EXAMPLES:

Construction::

    sage: M = Matroid(Matrix(QQ, [[1, 0, 0, 0, 1, 1, 1],
    ....:                         [0, 1, 0, 1, 0, 1, 1],
    ....:                         [0, 0, 1, 1, 1, 0, 1]]))
    sage: sorted(M.groundset())
    [0, 1, 2, 3, 4, 5, 6]
    sage: M.rank([0, 1, 2])
    3
    sage: M.rank([0, 1, 5])
    2

Minors::

    sage: M = Matroid(Matrix(QQ, [[1, 0, 0, 0, 1, 1, 1],
    ....:                         [0, 1, 0, 1, 0, 1, 1],
    ....:                         [0, 0, 1, 1, 1, 0, 1]]))
    sage: N = M / [2] \ [3, 4]
    sage: sorted(N.groundset())
    [0, 1, 5, 6]
    sage: N.full_rank()
    2

Testing. Note that the abstract base class does not support pickling::

    sage: M = sage.matroids.matroid.Matroid()
    sage: TestSuite(M).run(skip="_test_pickling")

REFERENCES
==========

..  [BC79] R. E. Bixby, W. H. Cunningham, Matroids, Graphs, and 3-Connectivity. In Graph theory and related topics (Proc. Conf., Univ. Waterloo, Waterloo, ON, 1977), 91-103
..  [CMO11] C. Chun, D. Mayhew, J. Oxley, A chain theorem for internally 4-connected binary matroids. J. Combin. Theory Ser. B 101 (2011), 141-189.
..  [CMO12] C. Chun, D. Mayhew, J. Oxley,  Towards a splitter theorem for internally 4-connected binary matroids. J. Combin. Theory Ser. B 102 (2012), 688-700.
..  [GG12] Jim Geelen and Bert Gerards, Characterizing graphic matroids by a system of linear equations, submitted, 2012. Preprint: http://www.gerardsbase.nl/papers/geelen_gerards=testing-graphicness%5B2013%5D.pdf
..  [GR01] C.Godsil and G.Royle, Algebraic Graph Theory. Graduate Texts in Mathematics, Springer, 2001.
..  [Hlineny] Petr Hlineny, "Equivalence-free exhaustive generation of matroid representations", Discrete Applied Mathematics 154 (2006), pp. 1210-1222.
..  [Hochstaettler] Winfried Hochstaettler, "About the Tic-Tac-Toe Matroid", preprint.
..  [Lyons] R. Lyons, Determinantal probability measures. Publications Mathematiques de l'Institut des Hautes Etudes Scientifiques 98(1)  (2003), pp. 167-212.
..  [Oxley1] James Oxley, "Matroid theory", Oxford University Press, 1992.
..  [Oxley] James Oxley, "Matroid Theory, Second Edition". Oxford University Press, 2011.
..  [Pen12] R. Pendavingh, On the evaluation at (-i, i) of the Tutte polynomial of a binary matroid. Preprint: http://arxiv.org/abs/1203.0910

AUTHORS:

- Michael Welsh (2013-04-03): Changed flats() to use SetSystem
- Michael Welsh (2013-04-01): Added is_3connected(), using naive algorithm
- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version

Methods
=======
"""
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com >
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz >
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject
from itertools import combinations, permutations
from set_system cimport SetSystem
from sage.combinat.subset import Subsets
from sage.misc.misc import subsets

from utilities import newlabel, sanitize_contractions_deletions
from sage.rings.all import ZZ
from sage.numerical.mip import MixedIntegerLinearProgram


# On some systems, macros "minor()" and "major()" are defined in system header
# files. This will undefine those:
cdef extern from "minorfix.h":
    pass

cdef class Matroid(SageObject):
    r"""
    The abstract matroid class, from which all matroids are derived. Do not
    use this class directly!

    To implement a subclass, the least you should do is implement the
    ``__init__()``, ``_rank()`` and ``groundset()`` methods. See the source of
    :mod:`rank_matroid.py <sage.matroids.rank_matroid>` for a bare-bones
    example of this.

    EXAMPLES:

    In a partition matroid, a subset is independent if it has at most one
    element from each partition. The following is a very basic implementation,
    in which the partition is specified as a list of lists::

        sage: class PartitionMatroid(sage.matroids.matroid.Matroid):
        ....:     def __init__(self, partition):
        ....:         self.partition = partition
        ....:         E = set()
        ....:         for P in partition:
        ....:             E.update(P)
        ....:         self.E = frozenset(E)
        ....:     def groundset(self):
        ....:         return self.E
        ....:     def _rank(self, X):
        ....:         X2 = set(X)
        ....:         used_indices = set()
        ....:         rk = 0
        ....:         while len(X2) > 0:
        ....:             e = X2.pop()
        ....:             for i in range(len(self.partition)):
        ....:                 if e in self.partition[i]:
        ....:                     if i not in used_indices:
        ....:                         used_indices.add(i)
        ....:                         rk = rk + 1
        ....:                     break
        ....:         return rk
        ....:
        sage: M = PartitionMatroid([[1, 2], [3, 4, 5], [6, 7]])
        sage: M.full_rank()
        3
        sage: M.tutte_polynomial(var('x'), var('y'))
        x^2*y^2 + 2*x*y^3 + y^4 + x^3 + 3*x^2*y + 3*x*y^2 + y^3

    .. NOTE::

        The abstract base class has no idea about the data used to represent
        the matroid. Hence some methods need to be customized to function
        properly.

        Necessary:

        - ``def __init__(self, ...)``
        - ``def groundset(self)``
        - ``def _rank(self, X)``

        Representation:

        - ``def _repr_(self)``

        Comparison:

        - ``def __hash__(self)``
        - ``def __eq__(self, other)``
        - ``def __ne__(self, other)``

        In Cythonized classes, use ``__richcmp__()`` instead of ``__eq__()``,
        ``__ne__()``.

        Copying, loading, saving:

        - ``def __copy__(self)``
        - ``def __deepcopy__(self, memo={})``
        - ``def __reduce__(self)``

        See, for instance,
        :mod:`rank_matroid.py <sage.matroids.rank_matroid>` or
        :mod:`circuit_closures_matroid.pyx <sage.matroids.circuit_closures_matroid>`
        for sample implementations of these.

    .. NOTE::

        Many methods (such as ``M.rank()``) have a companion method whose name
        starts with an underscore (such as ``M._rank()``). The method with the
        underscore does not do any checks on its input. For instance, it may
        assume of its input that

        - Any input that should be a subset of the groundset, is one. The
          interface is compatible with Python's ``frozenset`` type.
        - Any input that should be a list of things, supports iteration, and
          recursively these rules apply to its members.

        Using the underscored version could improve the speed of code a
        little, but will generate more cryptic error messages when presented
        with wrong input. In some instances, no error might occur and a
        nonsensical answer returned.

        A subclass should always override the underscored method, if
        available, and as a rule leave the regular method alone.

    """

    # virtual methods

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        .. NOTE::

            Subclasses should implement this method. The return type should be
            frozenset or any type with compatible interface.

        EXAMPLES::

            sage: M = sage.matroids.matroid.Matroid()
            sage: M.groundset()
            Traceback (most recent call last):
            ...
            NotImplementedError: subclasses need to implement this.

        """
        raise NotImplementedError("subclasses need to implement this.")

    cpdef _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to
        have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        Integer.

        .. NOTE::

            Subclasses should implement this method.

        EXAMPLES::

            sage: M = sage.matroids.matroid.Matroid()
            sage: M._rank([0, 1, 2])
            Traceback (most recent call last):
            ...
            NotImplementedError: subclasses need to implement this.

        """
        raise NotImplementedError("subclasses need to implement this.")

    # internal methods, assuming verified input

    # for better efficiency, its best to override the following methods in
    # each derived class

    cpdef _max_independent(self, X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._max_independent(set(['a', 'c', 'd', 'e', 'f'])))
            ['a', 'd', 'e', 'f']

        """
        res = set([])
        r = 0
        for e in X:
            res.add(e)
            if self._rank(res) > r:
                r += 1
            else:
                res.discard(e)
        return frozenset(res)

    cpdef _circuit(self, X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.
        A ``ValueError`` is raised if the set contains no circuit.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                             set(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                                       set(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set.

        """
        Z = set(X)
        if self._is_independent(X):
            raise ValueError("no circuit in independent set.")
        l = len(X) - 1
        for x in X:
            Z.discard(x)
            if self._rank(Z) == l:
                Z.add(x)
            else:
                l -= 1
        return frozenset(Z)

    cpdef _fundamental_circuit(self, B, e):
        r"""
        Return the `B`-fundamental circuit using `e`.

        Internal version that does no input checking.

        INPUT:

        - ``B`` -- a basis of the matroid.
        - ``e`` -- an element not in ``B``.

        OUTPUT:

        A set of elements.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._fundamental_circuit(frozenset('defg'), 'c'))
            ['c', 'd', 'e', 'f']
        """
        return self._circuit(B.union([e]))

    cpdef _closure(self, X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._closure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        """
        X = set(X)
        Y = self.groundset().difference(X)
        r = self._rank(X)
        for y in Y:
            X.add(y)
            if self._rank(X) > r:
                X.discard(y)
        return frozenset(X)

    cpdef _corank(self, X):
        """
        Return the corank of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._corank(set(['a', 'e', 'g', 'd', 'h']))
            4

        """
        return len(X) + self._rank(self.groundset().difference(X)) - self.full_rank()

    cpdef _max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._max_coindependent(set(['a', 'c', 'd', 'e', 'f'])))
            ['a', 'c', 'd', 'e']

        """
        res = set([])
        r = 0
        for e in X:
            res.add(e)
            if self._corank(res) > r:
                r += 1
            else:
                res.discard(e)
        return frozenset(res)

    cpdef _cocircuit(self, X):
        """
        Return a minimal codependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.
        A ``ValueError`` is raised if the set contains no cocircuit.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(sage.matroids.matroid.Matroid._cocircuit(M,
            ....:                             set(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._cocircuit(M,
            ....:                                       set(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no cocircuit in coindependent set.

        """
        Z = set(X)
        if self._is_coindependent(X):
            raise ValueError("no cocircuit in coindependent set.")
        l = len(X) - 1
        for x in X:
            Z.discard(x)
            if self._corank(Z) == l:
                Z.add(x)
            else:
                l -= 1
        return frozenset(Z)

    cpdef _fundamental_cocircuit(self, B, e):
        r"""
        Return the `B`-fundamental circuit using `e`.

        Internal version that does no input checking.

        INPUT:

        - ``B`` -- a basis of the matroid.
        - ``e`` -- an element of ``B``.

        OUTPUT:

        A set of elements.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._fundamental_cocircuit('abch', 'c'))
            ['c', 'd', 'e', 'f']
        """
        return self._cocircuit(self.groundset().difference(B).union([e]))

    cpdef _coclosure(self, X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._coclosure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        """
        X = set(X)
        Y = self.groundset().difference(X)
        r = self._corank(X)
        for y in Y:
            X.add(y)
            if self._corank(X) > r:
                X.discard(y)
        return frozenset(X)

    cpdef _augment(self, X, Y):
        r"""
        Return a maximal subset `I` of `Y` such that `r(X + I)=r(X) + r(I)`.

        This version of ``augment`` does no type checking. In particular,
        ``Y`` is assumed to be disjoint from ``X``.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.
        - ``Y`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``, and disjoint from ``X``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._augment(set(['a']), set(['e', 'f', 'g', 'h'])))
            ['e', 'g', 'h']

        """
        X = set(X)
        res = set([])
        r = self._rank(X)
        for e in Y:
            X.add(e)
            if self.rank(X) > r:
                r += 1
                res.add(e)
        return frozenset(res)

    # override the following methods for even better efficiency

    cpdef _is_independent(self, X):
        """
        Test if input is independent.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_independent(set(['a', 'b', 'c']))
            True
            sage: M._is_independent(set(['a', 'b', 'c', 'd']))
            False

        """
        return len(X) == self._rank(X)

    cpdef _is_basis(self, X):
        """
        Test if input is a basis.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        .. WARNING::

            This method assumes that ``X`` has the right size to be a basis,
            i.e. ``len(X) == self.full_rank()``. Otherwise its behavior is
            undefined.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_basis(set(['a', 'b', 'c', 'e']))
            True
            sage: M._is_basis(set(['a', 'b', 'c', 'd']))
            False

        If ``X`` does not have the right size, behavior can become
        unpredictable::

            sage: M._is_basis(set(['a', 'b', 'c']))
            True

        """
        return self._is_independent(X)

    cpdef _is_circuit(self, X):
        """
        Test if input is a circuit.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_circuit(set(['a', 'b', 'c', 'd']))
            True
            sage: M._is_circuit(set(['a', 'b', 'c', 'e']))
            False
            sage: M._is_circuit(set(['a', 'b', 'c', 'd', 'e']))
            False

        """
        l = len(X) - 1
        if self._rank(X) != l:
            return False
        Z = set(X)
        for x in X:
            Z.discard(x)
            if self._rank(frozenset(Z)) < l:
                return False
            Z.add(x)
        return True

    cpdef _is_closed(self, X):
        """
        Test if input is a closed set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_closed(set(['a', 'b', 'c', 'd']))
            True
            sage: M._is_closed(set(['a', 'b', 'c', 'e']))
            False

        """
        X = set(X)
        Y = self.groundset().difference(X)
        r = self._rank(X)
        for y in Y:
            X.add(y)
            if self._rank(frozenset(X)) == r:
                return False
            X.discard(y)
        return True

    cpdef _is_coindependent(self, X):
        """
        Test if input is coindependent.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_coindependent(set(['a', 'b', 'c']))
            True
            sage: M._is_coindependent(set(['a', 'b', 'c', 'd']))
            False

        """
        return self._corank(X) == len(X)

    cpdef _is_cobasis(self, X):
        """
        Test if input is a cobasis.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. WARNING::

            This method assumes ``X`` has the size of a cobasis, i.e.
            ``len(X) == self.full_corank()``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_cobasis(set(['a', 'b', 'c', 'e']))
            True
            sage: M._is_cobasis(set(['a', 'b', 'c', 'd']))
            False
            sage: M._is_cobasis(set(['a', 'b', 'c']))
            False

        """
        return self._is_basis(self.groundset().difference(X))

    cpdef _is_cocircuit(self, X):
        """
        Test if input is a cocircuit.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_cocircuit(set(['a', 'b', 'c', 'd']))
            True
            sage: M._is_cocircuit(set(['a', 'b', 'c', 'e']))
            False
            sage: M._is_cocircuit(set(['a', 'b', 'c', 'd', 'e']))
            False

        """
        l = len(X) - 1
        if self._corank(X) != l:
            return False
        Z = set(X)
        for x in X:
            Z.discard(x)
            if self._corank(frozenset(Z)) < l:
                return False
            Z.add(x)
        return True

    cpdef _is_coclosed(self, X):
        """
        Test if input is a coclosed set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_coclosed(set(['a', 'b', 'c', 'd']))
            True
            sage: M._is_coclosed(set(['a', 'b', 'c', 'e']))
            False

        """
        X = set(X)
        Y = self.groundset().difference(X)
        r = self._corank(X)
        for y in Y:
            X.add(y)
            if self._corank(frozenset(X)) == r:
                return False
            X.discard(y)
        return True

    cpdef _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.
        - ``deletions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        OUTPUT:

        A matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))
            M / {'a'} \ {'b', 'c'}, where M is Vamos: Matroid of rank 4 on 8
            elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
            {'e', 'f', 'g', 'h'}, {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'}},
            4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        import minor_matroid
        return minor_matroid.MinorMatroid(self, contractions, deletions)

    cpdef _has_minor(self, N):
        """
        Test if matroid has the specified minor.

        INPUT:

        - ``N`` -- An instance of a ``Matroid`` object.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._has_minor(matroids.Whirl(3))
            False
            sage: M._has_minor(matroids.Uniform(2, 4))
            True

        .. TODO::

            This important method can (and should) be optimized considerably.
            See [Hlineny]_ p.1219 for hints to that end.
        """
        if self is N:
            return True
        rd = self.full_rank() - N.full_rank()
        cd = self.full_corank() - N.full_corank()
        if rd < 0 or cd < 0:
            return False
        YY = self.dual().independent_r_sets(cd)
        for X in self.independent_r_sets(rd):
            for Y in YY:
                if X.isdisjoint(Y):
                    if N._is_isomorphic(self._minor(contractions=X, deletions=Y)):
                        return True
        return False

    cpdef _line_length(self, F):
        """
        Compute the length of the line specified through flat ``F``.

        This is the number of elements in `si(M / F)`.

        INPUT:

        - ``F`` -- a subset of the groundset, assumed to be a closed set of
          rank `r(M) - 2`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: M._line_length(['d'])
            5
            sage: M = matroids.named_matroids.NonPappus()
            sage: M._line_length(['d'])
            6
        """
        return len(self.minor(contractions=F).simplify())

    cpdef _extension(self, element, hyperplanes):
        """
        Extend the matroid by a new element.

        The result is a matroid on ``self.groundset() + {element}``, where
        ``element`` is contained in exactly the hyperplanes of ``self``
        specified by ``hyperplanes``.

        INPUT:

        - ``element`` -- a hashable object not in ``self.groundset()``.
        - ``hyperplanes`` -- the set of hyperplanes of a linear subclass of
          ``self``.

        OUTPUT:

        A matroid.

        EXAMPLES::

            sage: M = matroids.Uniform(3, 6)
            sage: H = [frozenset([0, 1])]
            sage: N = M._extension(6, H)
            sage: N
            Matroid of rank 3 on 7 elements with 34 bases
            sage: [sorted(C) for C in N.circuits() if len(C) == 3]
            [[0, 1, 6]]
        """
        import basis_matroid
        return basis_matroid.BasisMatroid(self)._extension(element, hyperplanes)

    # ** user-facing methods **

    # Display:
    # cpdef _latex_(self):
    #     return "\\texttt{Matroid on groundset }" + latex(self.groundset())

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sage.matroids.matroid.Matroid._repr_(M)
            'Matroid of rank 4 on 8 elements'
        """
        S = "Matroid of rank "
        S = S + str(self.rank()) + " on " + str(self.size()) + " elements"
        return S

    # cpdef show(self):
    # Show either the graph, or the matrix with labels, or the lattice,
    # or (in rank 3) the geometric diagram.
    #     raise NotImplementedError

    def __len__(self):
        """
        Return the size of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: len(M)
            8
        """
        return self.size()

    cpdef size(self):
        """
        Return the size of the groundset.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.size()
            8
        """
        if self._stored_size == 0:
            self._stored_size = len(self.groundset())
        return self._stored_size

    # User-visible methods

    cpdef rank(self, X=None):
        r"""
        Return the rank of ``X``.

        The *rank* of a subset `X` is the size of the largest independent set
        contained in `X`.

        If ``X`` is ``None``, the rank of the groundset is returned.

        INPUT:

        - ``X`` -- (default: ``None``) a subset of the groundset

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.rank()
            3
            sage: M.rank(['a', 'b', 'f'])
            2
            sage: M.rank(['a', 'b', 'x'])
            Traceback (most recent call last):
            ...
            ValueError: input is not a subset of the groundset.

        """
        if X is None:
            return self.full_rank()
        else:
            if not self.groundset().issuperset(X):
                raise ValueError("input is not a subset of the groundset.")
        return self._rank(frozenset(X))

    cpdef full_rank(self):
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.full_rank()
            4
            sage: M.dual().full_rank()
            4
        """
        if self._stored_full_rank == 0:
            self._stored_full_rank = self._rank(self.groundset())
        return self._stored_full_rank

    cpdef basis(self):
        r"""
        Return an arbitrary basis of the matroid.

        A *basis* is an inclusionwise maximal independent set.

        .. NOTE::

            The output of this method can change in between calls.

        OUTPUT:

        Set of elements.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: B = M.basis()
            sage: M.is_basis(B)
            True
            sage: len(B)
            3
            sage: M.rank(B)
            3
            sage: M.full_rank()
            3
        """
        return self._max_independent(self.groundset())

    cpdef max_independent(self, X):
        """
        Compute a maximal independent subset of ``X``.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Subset of ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.max_independent(['a', 'c', 'd', 'e', 'f']))
            ['a', 'd', 'e', 'f']
            sage: M.max_independent(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input is not a subset of the groundset.")
        return self._max_independent(frozenset(X))

    cpdef circuit(self, X=None):
        """
        Return a circuit.

        A *circuit* of a matroid is an inclusionwise minimal dependent subset.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``, or ``None``. In the
          latter case, a circuit contained in the full groundset is returned.

        OUTPUT:

        Set of elements.

        - If ``X`` is not ``None``, the output is a circuit contained in ``X``
          if such a circuit exists. Otherwise an error is raised.
        - If ``X`` is ``None``, the output is a circuit contained in
          ``self.groundset()`` if such a circuit exists. Otherwise an error is
          raised.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.circuit(['a', 'c', 'd', 'e', 'f']))
            ['c', 'd', 'e', 'f']
            sage: sorted(M.circuit(['a', 'c', 'd']))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set.
            sage: M.circuit(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
            sage: sorted(M.circuit())
            ['a', 'b', 'c', 'd']
        """
        if X is None:
            X = self.groundset()
        elif not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._circuit(frozenset(X))

    cpdef fundamental_circuit(self, B, e):
        r"""
        Return the `B`-fundamental circuit using `e`.

        If `B` is a basis, and `e` an element not in `B`, then the
        `B`-*fundamental circuit* using `e` is the unique matroid circuit
        contained in `B\cup e`.

        INPUT:

        - ``B`` -- a basis of the matroid.
        - ``e`` -- an element not in ``B``.

        OUTPUT:

        A set of elements.

        .. SEEALSO::

            :meth:`M.circuit() <Matroid.circuit>`,
            :meth:`M.basis() <Matroid.basis>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.fundamental_circuit('defg', 'c'))
            ['c', 'd', 'e', 'f']
        """
        B = frozenset(B)
        if not self.is_basis(B):
            raise ValueError("input B is not a basis of the matroid.")
        if not e in self.groundset():
            raise ValueError("input e is not an element of the groundset.")
        return self._fundamental_circuit(B, e)

    cpdef closure(self, X):
        """
        Return the closure of a set.

        A set is *closed* if adding any extra element to it will increase the
        rank of the set. The *closure* of a set is the smallest closed set
        containing it.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Set of elements containing ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.closure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
            sage: M.closure(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.

        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._closure(frozenset(X))

    cpdef k_closure(self, X, k):
        r"""
        Return the ``k``-closure of ``X``.

        A set `S` is `k`-*closed* if the closure of any `k` element subsets
        is contained in `S`. The `k`-*closure* of a set `X` is the smallest
        `k`-closed set containing `X`.

        INPUT:

        - ``X`` -- a subset of ``self.groundset()``
        - ``k`` -- a positive integer

        EXAMPLES::

            sage: m = matrix([[1,2,5,2], [0,2,1,0]])
            sage: M = Matroid(m)
            sage: sorted(M.k_closure({1,3}, 2))
            [0, 1, 2, 3]
            sage: sorted(M.k_closure({0,1}, 1))
            [0, 1, 3]
            sage: sorted(M.k_closure({1,2}, 1))
            [1, 2]

            sage: Q = RootSystem(['D',4]).root_lattice()
            sage: m = matrix(map(lambda x: x.to_vector(), Q.positive_roots()))
            sage: m = m.transpose(); m
            [1 0 0 0 1 0 0 0 1 1 1 1]
            [0 1 0 0 1 1 1 1 1 1 1 2]
            [0 0 1 0 0 0 1 1 1 0 1 1]
            [0 0 0 1 0 1 0 1 0 1 1 1]
            sage: M = Matroid(m)
            sage: sorted(M.k_closure({0,2,3,11}, 3))
            [0, 2, 3, 11]
            sage: sorted(M.k_closure({0,2,3,11}, 4))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            sage: sorted(M.k_closure({0,1}, 4))
            [0, 1, 4]
        """
        cdef int cur
        cdef frozenset S, cl
        cur = 0
        S = frozenset(X)
        while cur != len(S):
            cur = len(S)
            cl = frozenset([])
            for T in Subsets(S, min(k,cur)):
                cl = cl.union(self.closure(T))
            S = cl
        return S

    cpdef augment(self, X, Y=None):
        r"""
        Return a maximal subset `I` of `Y - X` such that
        `r(X + I) = r(X) + r(I)`.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.
        - ``Y`` -- a subset of ``self.groundset()``. If ``Y`` is ``None``,
          the full groundset is used.

        OUTPUT:

        A subset of `Y - X`.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.augment(['a'], ['e', 'f', 'g', 'h']))
            ['e', 'g', 'h']
            sage: sorted(M.augment(['a']))
            ['b', 'c', 'e']
            sage: sorted(M.augment(['x']))
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
            sage: sorted(M.augment(['a'], ['x']))
            Traceback (most recent call last):
            ...
            ValueError: input Y is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        if Y is None:
            Y = self.groundset()
        else:
            if not self.groundset().issuperset(Y):
                raise ValueError("input Y is not a subset of the groundset.")

        return self._augment(frozenset(X), frozenset(Y).difference(X))

    cpdef corank(self, X=None):
        r"""
        Return the corank of ``X``, or the corank of the groundset if ``X`` is
        ``None``.

        The *corank* of a set `X` is the rank of `X` in the dual matroid.

        If ``X`` is ``None``, the corank of the groundset is returned.

        INPUT:

        - ``X`` -- (default: ``None``) a subset of the groundset.

        OUTPUT:

        Integer.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.rank() <sage.matroids.matroid.Matroid.rank>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.corank()
            4
            sage: M.corank('cdeg')
            3
            sage: M.rank(['a', 'b', 'x'])
            Traceback (most recent call last):
            ...
            ValueError: input is not a subset of the groundset.
        """
        if X is None:
            X = self.groundset()
        else:
            if not self.groundset().issuperset(X):
                raise ValueError("input is not a subset of the groundset.")
        return self._corank(frozenset(X))

    cpdef full_corank(self):
        """
        Return the corank of the matroid.

        The *corank* of the matroid equals the rank of the dual matroid. It is
        given by ``M.size() - M.full_rank()``.

        OUTPUT:

        Integer.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.full_rank() <sage.matroids.matroid.Matroid.full_rank>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.full_corank()
            4
        """
        return self.size() - self.full_rank()

    cpdef cobasis(self):
        """
        Return an arbitrary cobasis of the matroid.

        A *cobasis* is the complement of a basis. A cobasis is
        a basis of the dual matroid.

        .. NOTE::

            Output can change between calls.

        OUTPUT:

        A set of elements.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.full_rank() <sage.matroids.matroid.Matroid.full_rank>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: B = M.cobasis()
            sage: M.is_cobasis(B)
            True
            sage: len(B)
            6
            sage: M.corank(B)
            6
            sage: M.full_corank()
            6

        """
        return self.max_coindependent(self.groundset())

    cpdef max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        A set is *coindependent* if it is independent in the dual matroid.
        A set is coindependent if and only if the complement is *spanning*
        (i.e. contains a basis of the matroid).

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        A subset of ``X``.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.max_independent() <sage.matroids.matroid.Matroid.max_independent>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.max_coindependent(['a', 'c', 'd', 'e', 'f']))
            ['a', 'c', 'd', 'e']
            sage: M.max_coindependent(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input is not a subset of the groundset.")
        return self._max_coindependent(frozenset(X))

    cpdef coclosure(self, X):
        """
        Return the coclosure of a set.

        A set is *coclosed* if it is closed in the dual matroid. The
        *coclosure* of `X` is the smallest coclosed set containing `X`.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        A set of elements containing ``X``.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.coclosure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
            sage: M.coclosure(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._coclosure(frozenset(X))

    cpdef cocircuit(self, X=None):
        """
        Return a cocircuit.

        A *cocircuit* is an inclusionwise minimal subset that is dependent in
        the dual matroid.

        INPUT:

        - ``X`` -- (default: ``None``) a subset of ``self.groundset()``.

        OUTPUT:

        A set of elements.

        - If ``X`` is not ``None``, the output is a cocircuit contained in
          ``X`` if such a cocircuit exists. Otherwise an error is raised.
        - If ``X`` is ``None``, the output is a cocircuit contained in
          ``self.groundset()`` if such a cocircuit exists. Otherwise an error
          is raised.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.cocircuit(['a', 'c', 'd', 'e', 'f']))
            ['c', 'd', 'e', 'f']
            sage: sorted(M.cocircuit(['a', 'c', 'd']))
            Traceback (most recent call last):
            ...
            ValueError: no cocircuit in coindependent set.
            sage: M.cocircuit(['x'])
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
            sage: sorted(M.cocircuit())
            ['e', 'f', 'g', 'h']
        """
        if X is None:
            X = self.groundset()
        elif not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._cocircuit(frozenset(X))

    cpdef fundamental_cocircuit(self, B, e):
        r"""
        Return the `B`-fundamental cocircuit using `e`.

        If `B` is a basis, and `e` an element of `B`, then the
        `B`-*fundamental cocircuit* using `e` is the unique matroid cocircuit
        that intersects `B` only in `e`.

        This is equal to
        ``M.dual().fundamental_circuit(M.groundset().difference(B), e)``.

        INPUT:

        - ``B`` -- a basis of the matroid.
        - ``e`` -- an element of ``B``.

        OUTPUT:

        A set of elements.

        .. SEEALSO::

            :meth:`M.cocircuit() <Matroid.cocircuit>`,
            :meth:`M.basis() <Matroid.basis>`,
            :meth:`M.fundamental_circuit() <Matroid.fundamental_circuit>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M.fundamental_cocircuit('abch', 'c'))
            ['c', 'd', 'e', 'f']
        """
        B = frozenset(B)
        if not self.is_basis(B):
            raise ValueError("input B is not a basis of the matroid.")
        if not e in B:
            raise ValueError("input e is not an element of B.")
        return self._fundamental_cocircuit(B, e)

    cpdef loops(self):
        r"""
        Return the set of loops of the matroid.

        A *loop* is a one-element dependent set.

        OUTPUT:

        A set of elements.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.loops()
            frozenset()
            sage: (M / ['a', 'b']).loops()
            frozenset({'f'})
        """
        return self._closure(set())

    cpdef is_independent(self, X):
        r"""
        Check if a subset is independent in the matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_independent('abc')
            True
            sage: M.is_independent('abcd')
            False
            sage: M.is_independent('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_independent(frozenset(X))

    cpdef is_dependent(self, X):
        r"""
        Check if a subset is dependent in the matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_dependent('abc')
            False
            sage: M.is_dependent('abcd')
            True
            sage: M.is_dependent('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return not self._is_independent(frozenset(X))

    cpdef is_basis(self, X):
        r"""
        Check if a subset is a basis of the matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_basis('abc')
            False
            sage: M.is_basis('abce')
            True
            sage: M.is_basis('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        if len(X) != self.full_rank():
            return False
        return self._is_basis(frozenset(X))

    cpdef is_closed(self, X):
        r"""
        Test if a subset is a closed set of the matroid.

        A set is *closed* if adding any element to it will increase the rank
        of the set.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_closed('abc')
            False
            sage: M.is_closed('abcd')
            True
            sage: M.is_closed('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_closed(frozenset(X))

    cpdef is_subset_k_closed(self, X, k):
        r"""
        Test if ``X`` is a ``k``-closed set of the matroid.

        A set `S` is `k`-*closed* if the closure of any `k` element subsets
        is contained in `S`.

        INPUT:

        - ``X`` -- a subset of ``self.groundset()``
        - ``k`` -- a positive integer

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.k_closure() <sage.matroids.matroid.Matroid.k_closure>`

        EXAMPLES::

            sage: m = matrix([[1,2,5,2], [0,2,1,0]])
            sage: M = Matroid(m)
            sage: M.is_subset_k_closed({1,3}, 2)
            False
            sage: M.is_subset_k_closed({0,1}, 1)
            False
            sage: M.is_subset_k_closed({1,2}, 1)
            True

            sage: Q = RootSystem(['D',4]).root_lattice()
            sage: m = matrix(map(lambda x: x.to_vector(), Q.positive_roots()))
            sage: m = m.transpose(); m
            [1 0 0 0 1 0 0 0 1 1 1 1]
            [0 1 0 0 1 1 1 1 1 1 1 2]
            [0 0 1 0 0 0 1 1 1 0 1 1]
            [0 0 0 1 0 1 0 1 0 1 1 1]
            sage: M = Matroid(m)
            sage: M.is_subset_k_closed({0,2,3,11}, 3)
            True
            sage: M.is_subset_k_closed({0,2,3,11}, 4)
            False
            sage: M.is_subset_k_closed({0,1}, 4)
            False
            sage: M.is_subset_k_closed({0,1,4}, 4)
            True
        """
        if len(X) < k:
            return self.is_closed(X)

        cdef frozenset cl
        for T in Subsets(X, k):
            cl = self.closure(T)
            if not cl.issubset(T):
                return False
        return True

    cpdef is_circuit(self, X):
        r"""
        Test if a subset is a circuit of the matroid.

        A *circuit* is an inclusionwise minimal dependent subset.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_circuit('abc')
            False
            sage: M.is_circuit('abcd')
            True
            sage: M.is_circuit('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_circuit(frozenset(X))

    cpdef coloops(self):
        r"""
        Return the set of coloops of the matroid.

        A *coloop* is a one-element cocircuit. It is a loop of the dual of the
        matroid.

        OUTPUT:

        A set of elements.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.loops() <sage.matroids.matroid.Matroid.loops>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano().dual()
            sage: M.coloops()
            frozenset()
            sage: (M \ ['a', 'b']).coloops()
            frozenset({'f'})
        """
        return self._coclosure(set())

    cpdef is_coindependent(self, X):
        r"""
        Check if a subset is coindependent in the matroid.

        A set is *coindependent* if it is independent in the dual matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.is_independent() <sage.matroids.matroid.Matroid.is_independent>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_coindependent('abc')
            True
            sage: M.is_coindependent('abcd')
            False
            sage: M.is_coindependent('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_coindependent(frozenset(X))

    cpdef is_codependent(self, X):
        r"""
        Check if a subset is codependent in the matroid.

        A set is *codependent* if it is dependent in the dual of the matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.is_dependent() <sage.matroids.matroid.Matroid.is_dependent>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_codependent('abc')
            False
            sage: M.is_codependent('abcd')
            True
            sage: M.is_codependent('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return not self._is_coindependent(X)

    cpdef is_cobasis(self, X):
        r"""
        Check if a subset is a cobasis of the matroid.

        A *cobasis* is the complement of a basis. It is a basis of the dual
        matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.is_basis() <sage.matroids.matroid.Matroid.is_basis>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_cobasis('abc')
            False
            sage: M.is_cobasis('abce')
            True
            sage: M.is_cobasis('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        if len(X) != self.full_corank():
            return False
        return self._is_cobasis(frozenset(X))

    cpdef is_cocircuit(self, X):
        r"""
        Test if a subset is a cocircuit of the matroid.

        A *cocircuit* is an inclusionwise minimal subset that is dependent in
        the dual matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.is_circuit() <sage.matroids.matroid.Matroid.is_circuit>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_cocircuit('abc')
            False
            sage: M.is_cocircuit('abcd')
            True
            sage: M.is_cocircuit('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_cocircuit(frozenset(X))

    cpdef is_coclosed(self, X):
        r"""
        Test if a subset is a coclosed set of the matroid.

        A set is *coclosed* if it is a closed set of the dual matroid.

        INPUT:

        - ``X`` -- A subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.is_closed() <sage.matroids.matroid.Matroid.is_closed>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_coclosed('abc')
            False
            sage: M.is_coclosed('abcd')
            True
            sage: M.is_coclosed('abcx')
            Traceback (most recent call last):
            ...
            ValueError: input X is not a subset of the groundset.
        """
        if not self.groundset().issuperset(X):
            raise ValueError("input X is not a subset of the groundset.")
        return self._is_coclosed(frozenset(X))

    # verification

    cpdef is_valid(self):
        r"""
        Test if the data obey the matroid axioms.

        The default implementation checks the (disproportionately slow) rank
        axioms. If `r` is the rank function of a matroid, we check, for all
        pairs `X, Y` of subsets,

        * `0 \leq r(X) \leq |X|`
        * If `X \subseteq Y` then `r(X) \leq r(Y)`
        * `r(X\cup Y) + r(X\cap Y) \leq r(X) + r(Y)`

        Certain subclasses may check other axioms instead.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.is_valid()
            True

        The following is the 'Escher matroid' by Brylawski and Kelly. See
        Example 1.5.5 in [Oxley]_ ::

            sage: M = Matroid(circuit_closures={2: [[1, 2, 3], [1, 4, 5]],
            ....: 3: [[1, 2, 3, 4, 5], [1, 2, 3, 6, 7], [1, 4, 5, 6, 7]]})
            sage: M.is_valid()
            False
        """
        E = list(self.groundset())
        for i in xrange(0, len(E) + 1):
            for X in combinations(E, i):
                XX = frozenset(X)
                rX = self._rank(XX)
                if rX > i:
                    return False
                for j in xrange(i, len(E) + 1):
                    for Y in combinations(E, j):
                        YY = frozenset(Y)
                        rY = self._rank(YY)
                        if XX.issubset(YY) and rX > rY:
                            return False
                        if (self._rank(XX.union(YY)) +
                           self._rank(XX.intersection(YY)) > rX + rY):
                            return False
        return True

    # enumeration

    cpdef circuits(self):
        """
        Return the list of circuits of the matroid.

        OUTPUT:

        An iterable containing all circuits.

        .. SEEALSO::

            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted([sorted(C) for C in M.circuits()])
            [['a', 'b', 'c', 'g'], ['a', 'b', 'd', 'e'], ['a', 'b', 'f'],
            ['a', 'c', 'd', 'f'], ['a', 'c', 'e'], ['a', 'd', 'g'],
            ['a', 'e', 'f', 'g'], ['b', 'c', 'd'], ['b', 'c', 'e', 'f'],
            ['b', 'd', 'f', 'g'], ['b', 'e', 'g'], ['c', 'd', 'e', 'g'],
            ['c', 'f', 'g'], ['d', 'e', 'f']]
        """
        C = set()
        for B in self.bases():
            C.update([self._circuit(B.union(set([e])))
                      for e in self.groundset().difference(B)])
        return list(C)

    cpdef nonspanning_circuits(self):
        """
        Return the list of nonspanning circuits of the matroid.

        A *nonspanning circuit* is a circuit whose rank is strictly smaller
        than the rank of the matroid.

        OUTPUT:

        An iterable containing all nonspanning circuits.

        .. SEEALSO::

            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`M.rank() <sage.matroids.matroid.Matroid.rank>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted([sorted(C) for C in M.nonspanning_circuits()])
            [['a', 'b', 'f'], ['a', 'c', 'e'], ['a', 'd', 'g'],
            ['b', 'c', 'd'], ['b', 'e', 'g'], ['c', 'f', 'g'],
            ['d', 'e', 'f']]
        """
        C = set()
        for N in self.nonbases():
            if self._rank(N) == self.full_rank() - 1:
                C.add(self._circuit(N))
        return list(C)

    cpdef cocircuits(self):
        """
        Return the list of cocircuits of the matroid.

        OUTPUT:

        An iterable containing all cocircuits.

        .. SEEALSO::

            :meth:`M.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted([sorted(C) for C in M.cocircuits()])
            [['a', 'b', 'c', 'g'], ['a', 'b', 'd', 'e'], ['a', 'c', 'd', 'f'],
            ['a', 'e', 'f', 'g'], ['b', 'c', 'e', 'f'], ['b', 'd', 'f', 'g'],
            ['c', 'd', 'e', 'g']]
        """
        C = set()
        for B in self.bases():
            C.update([self._cocircuit(self.groundset().difference(B).union(set([e]))) for e in B])
        return list(C)

    cpdef noncospanning_cocircuits(self):
        """
        Return the list of noncospanning cocircuits of the matroid.

        A *noncospanning cocircuit* is a cocircuit whose corank is strictly
        smaller than the corank of the matroid.

        OUTPUT:

        An iterable containing all nonspanning circuits.

        .. SEEALSO::

            :meth:`M.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`,
            :meth:`M.corank() <sage.matroids.matroid.Matroid.corank>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano().dual()
            sage: sorted([sorted(C) for C in M.noncospanning_cocircuits()])
            [['a', 'b', 'f'], ['a', 'c', 'e'], ['a', 'd', 'g'],
            ['b', 'c', 'd'], ['b', 'e', 'g'], ['c', 'f', 'g'],
            ['d', 'e', 'f']]
        """
        return self.dual().nonspanning_circuits()

    cpdef circuit_closures(self):
        """
        Return the list of closures of circuits of the matroid.

        A *circuit closure* is a closed set containing a circuit.

        OUTPUT:

        A dictionary containing the circuit closures of the matroid, indexed
        by their ranks.

        .. SEEALSO::

            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`M.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: CC = M.circuit_closures()
            sage: len(CC[2])
            7
            sage: len(CC[3])
            1
            sage: len(CC[1])
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: [sorted(X) for X in CC[3]]
            [['a', 'b', 'c', 'd', 'e', 'f', 'g']]
        """
        CC = [set([]) for r in xrange(self.rank() + 1)]
        for C in self.circuits():
            CC[len(C) - 1].add(self.closure(C))
        CC = dict([(r, CC[r]) for r in xrange(self.rank() + 1)
                  if len(CC[r]) > 0])
        return CC

    cpdef nonspanning_circuit_closures(self):
        """
        Return the list of closures of nonspanning circuits of the matroid.

        A *nonspanning circuit closure* is a closed set containing a
        nonspanning circuit.

        OUTPUT:

        A dictionary containing the nonspanning circuit closures of the
        matroid, indexed by their ranks.

        .. SEEALSO::

            :meth:`M.nonspanning_circuits() <sage.matroids.matroid.Matroid.nonspanning_circuits>`,
            :meth:`M.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: CC = M.nonspanning_circuit_closures()
            sage: len(CC[2])
            7
            sage: len(CC[3])
            Traceback (most recent call last):
            ...
            KeyError: 3
        """
        CC = [set([]) for r in xrange(self.rank() + 1)]
        for C in self.nonspanning_circuits():
            CC[len(C) - 1].add(self.closure(C))
        CC = dict([(r, CC[r]) for r in xrange(self.rank() + 1)
                   if len(CC[r]) > 0])
        return CC

    cpdef nonbases(self):
        r"""
        Return the list of nonbases of the matroid.

        A *nonbasis* is a set with cardinality ``self.full_rank()`` that is
        not a basis.

        OUTPUT:

        An iterable containing the nonbases of the matroid.

        .. SEEALSO::

            :meth:`M.basis() <sage.matroids.matroid.Matroid.basis>`

        EXAMPLES::

            sage: M = matroids.Uniform(2, 4)
            sage: list(M.nonbases())
            []
            sage: [sorted(X) for X in matroids.named_matroids.P6().nonbases()]
            [['a', 'b', 'c']]

        ALGORITHM:

        Test all subsets of the groundset of cardinality ``self.full_rank()``
        """
        cdef SetSystem res
        res = SetSystem(list(self.groundset()))
        for X in combinations(self.groundset(), self.full_rank()):
            if self._rank(X) < len(X):
                res.append(X)
        return res

    cpdef dependent_r_sets(self, long r):
        r"""
        Return the list of dependent subsets of fixed size.

        INPUT:

        - ``r`` -- a nonnegative integer.

        OUTPUT:

        An iterable containing all dependent subsets of size ``r``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.dependent_r_sets(3)
            []
            sage: sorted([sorted(X) for X in
            ....: matroids.named_matroids.Vamos().dependent_r_sets(4)])
            [['a', 'b', 'c', 'd'], ['a', 'b', 'e', 'f'], ['a', 'b', 'g', 'h'],
            ['c', 'd', 'e', 'f'], ['e', 'f', 'g', 'h']]

        ALGORITHM:

        Test all subsets of the groundset of cardinality ``r``
        """
        res = []
        for X in combinations(self.groundset(), r):
            X = frozenset(X)
            if self._rank(X) < len(X):
                res.append(X)
        return res

    cpdef bases(self):
        r"""
        Return the list of bases of the matroid.

        A *basis* is a maximal independent set.

        OUTPUT:

        An iterable containing all bases of the matroid.

        EXAMPLES::

            sage: M = matroids.Uniform(2, 4)
            sage: sorted([sorted(X) for X in M.bases()])
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]

        ALGORITHM:

        Test all subsets of the groundset of cardinality ``self.full_rank()``
        """
        cdef SetSystem res
        res = SetSystem(list(self.groundset()))
        for X in combinations(self.groundset(), self.full_rank()):
            if self._rank(X) == len(X):
                res.append(X)
        return res

    cpdef independent_r_sets(self, long r):
        r"""
        Return the list of size-``r`` independent subsets of the matroid.

        INPUT:

        - ``r`` -- a nonnegative integer.

        OUTPUT:

        An iterable containing all independent subsets of the matroid of
        cardinality ``r``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: M.independent_r_sets(4)
            []
            sage: sorted(sorted(M.independent_r_sets(3))[0])
            ['a', 'c', 'e']

        ALGORITHM:

        Test all subsets of the groundset of cardinality ``r``
        """
        res = []
        for X in combinations(self.groundset(), r):
            X = frozenset(X)
            if self._rank(X) == len(X):
                res.append(X)
        return res

    cpdef _extend_flags(self, flags):
        r"""
        Recursion for the ``self._flags(r)`` method.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: F = M._flags(1)
            sage: sorted(M._extend_flags(F)) == sorted(M._flags(2))
            True
        """
        newflags = []
        E = set(self.groundset())
        for [F, B, X] in flags:
            while len(X) > 0:
                x = min(X)  # TODO: Alert: this sort of thing will break in
                            # Python 3, I think. --SvZ
                newbase = B | set([x])
                newflat = self._closure(newbase)
                newX = X - newflat
                if max(newbase) == x:
                    if min(newflat - F) == x:
                        newflags.append([newflat, newbase, newX])
                X = newX
        return newflags

    cpdef _flags(self, r):
        r"""
        Compute rank-``r`` flats, with extra information for more speed.

        Used in the :meth:`flats() <sage.matroids.matroid.Matroid.flats>`
        method.

        INPUT:

        - ``r`` -- A natural number.

        OUTPUT:

        A list of triples `(F, B, X)` with `F` a flat of rank ``r``,
        `B` a lexicographically least basis of `F`, and `X` the remaining
        elements of the groundset that can potentially be lex-least extensions
        of the flat.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted([M._flags(1)])
            [[[frozenset({'a'}), {'a'}, frozenset({'b', 'c', 'd', 'e', 'f', 'g'})],
              [frozenset({'b'}), {'b'}, frozenset({'c', 'd', 'e', 'f', 'g'})],
              [frozenset({'c'}), {'c'}, frozenset({'d', 'e', 'f', 'g'})],
              [frozenset({'d'}), {'d'}, frozenset({'e', 'f', 'g'})],
              [frozenset({'e'}), {'e'}, frozenset({'f', 'g'})],
              [frozenset({'f'}), {'f'}, frozenset({'g'})],
              [frozenset({'g'}), {'g'}, frozenset()]]]
        """
        if r < 0:
            return []
        loops = self._closure(set())
        flags = [[loops, set(), self.groundset() - loops]]
        for r in xrange(r):
            flags = self._extend_flags(flags)
        return flags

    cpdef flats(self, r):
        r"""
        Return the collection of flats of the matroid of specified rank.

        A *flat* is a closed set.

        INPUT:

        - ``r`` -- A natural number.

        OUTPUT:

        An iterable containing all flats of rank ``r``.

        .. SEEALSO::

            :meth:`M.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted([sorted(F) for F in M.flats(2)])
            [['a', 'b', 'f'], ['a', 'c', 'e'], ['a', 'd', 'g'],
            ['b', 'c', 'd'], ['b', 'e', 'g'], ['c', 'f', 'g'],
            ['d', 'e', 'f']]
        """
        return SetSystem(list(self.groundset()),
                         subsets=[f[0] for f in self._flags(r)])

    cpdef coflats(self, r):
        r"""
        Return the collection of coflats of the matroid of specified corank.

        A *coflat* is a coclosed set.

        INPUT:

        - ``r`` -- a nonnegative integer.

        OUTPUT:

        An iterable containing all coflats of corank ``r``.

        .. SEEALSO::

            :meth:`M.coclosure() <sage.matroids.matroid.Matroid.coclosure>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Q6()
            sage: sorted([sorted(F) for F in M.coflats(2)])
            [['a', 'b'], ['a', 'c'], ['a', 'd', 'f'], ['a', 'e'], ['b', 'c'],
            ['b', 'd'], ['b', 'e'], ['b', 'f'], ['c', 'd'], ['c', 'e', 'f'],
            ['d', 'e']]
        """
        return self.dual().flats(r)

    def lattice_of_flats(self):
        """
        Return the lattice of flats of the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.lattice_of_flats()
            Finite lattice containing 16 elements
        """
        from sage.combinat.posets.lattices import LatticePoset
        F = [X for i in range(self.rank() + 1)
             for X in self.flats(i)]
        return LatticePoset((F, lambda x, y: x < y))

    cpdef hyperplanes(self):
        """
        Return the set of hyperplanes of the matroid.

        A *hyperplane* is a flat of rank ``self.full_rank() - 1``. A *flat* is
        a closed set.

        OUTPUT:

        An iterable containing all hyperplanes of the matroid.

        .. SEEALSO::

            :meth:`M.flats() <sage.matroids.matroid.Matroid.flats>`

        EXAMPLES::

            sage: M = matroids.Uniform(2, 3)
            sage: sorted([sorted(F) for F in M.hyperplanes()])
            [[0], [1], [2]]

        """
        return self.flats(self.full_rank() - 1)

    cpdef f_vector(self):
        r"""
        Return the `f`-vector of the matroid.

        The `f`-*vector* is a vector `(f_0, ..., f_r)`, where `f_i` is the
        number of flats of rank `i`, and `r` is the rank of the matroid.

        OUTPUT:

        List of integers.

        EXAMPLES::

            sage: M = matroids.named_matroids.BetsyRoss()
            sage: M.f_vector()
            [1, 11, 20, 1]

        """
        loops = self._closure(set())
        flags = [[loops, set(), self.groundset() - loops]]
        f_vec = [1]
        for r in xrange(self.full_rank()):
            flags = self._extend_flags(flags)
            f_vec.append(len(flags))
        return f_vec

    cpdef broken_circuits(self, ordering=None):
        r"""
        Return the list of broken circuits of ``self``.

        A *broken circuit* `B` for is a subset of the ground set under
        some total ordering `<` such that `B \cup \{ u \}` is a circuit
        and `u < b` for all `b \in B`.

        INPUT:

        - ``ordering`` -- a total ordering of the groundset given as a list

        EXAMPLES::

            sage: M = Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]])
            sage: M.broken_circuits()
            frozenset({frozenset({2, 3}), frozenset({4, 5}), frozenset({2, 4, 5})})
            sage: M.broken_circuits([5,4,3,2,1])
            frozenset({frozenset({1, 2}), frozenset({3, 4}), frozenset({1, 2, 4})})

        ::

            sage: M = Matroid(circuits=[[1,2,3], [1,4,5], [2,3,4,5]])
            sage: M.broken_circuits([5,4,3,2,1])
            frozenset({frozenset({1, 2}), frozenset({1, 4}), frozenset({2, 3, 4})})
        """
        if ordering is None:
            ordering = sorted(self.groundset())
        else:
            orderset = frozenset(ordering)
            if len(orderset) != len(self.groundset()) or orderset != self.groundset():
                raise ValueError("not an ordering of the ground set")
            ordering = list(ordering)
        ret = set()
        for C in self.circuits():
            for k in ordering:
                if k in C:
                    ret.add(C.difference([k]))
                    break
        return frozenset(ret)

    cpdef no_broken_circuits_sets(self, ordering=None):
        r"""
        Return the no broken circuits (NBC) sets of ``self``.

        An NBC set is a subset `A` of the ground set under some total
        ordering `<` such that `A` contains no broken circuit.

        INPUT:

        - ``ordering`` -- a total ordering of the groundset given as a list

        EXAMPLES::

            sage: M = Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets())
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 4), (1, 3, 5), (1, 2, 5), (1, 2, 4)}
            sage: SimplicialComplex(M.no_broken_circuits_sets([5,4,3,2,1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 4, 5), (2, 3, 5), (1, 3, 5), (2, 4, 5)}

        ::

            sage: M = Matroid(circuits=[[1,2,3], [1,4,5], [2,3,4,5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets([5,4,3,2,1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(2, 3, 5), (1, 3, 5), (2, 4, 5), (3, 4, 5)}
        """
        ret = []
        BC = self.broken_circuits(ordering)
        for r in range(self.rank() + 1):
            for I in self.independent_r_sets(r):
                add = True
                for b in BC:
                    if b.issubset(I):
                        add = False
                        break
                if add:
                    ret.append(I)
        return ret

    # polytopes

    def matroid_polytope(self):
        r"""
        Return the matroid polytope of ``self``.

        This is defined as the convex hull of the vertices

        .. MATH::

            e_B = \sum_{i \in B} e_i

        over all bases `B` of the matroid. Here `e_i` are the standard
        basis vectors of `\RR^n`. An arbitrary labelling of the
        groundset by `{0,\dots,n-1}` is chosen.

        .. SEEALSO::

            :meth:`independence_matroid_polytope`

        EXAMPLES::

            sage: M = matroids.Whirl(4)
            sage: P = M.matroid_polytope(); P
            A 7-dimensional polyhedron in ZZ^8 defined as the convex hull
            of 46 vertices

            sage: M = matroids.named_matroids.NonFano()
            sage: M.matroid_polytope()
            A 6-dimensional polyhedron in ZZ^7 defined as the convex hull
            of 29 vertices

        REFERENCE:

        .. [DLHK2007] J. A. De Loera, D. C. Haws, M. Köppe, Ehrhart polynomials
          of matroid polytopes and polymatroids. Discrete & Computational
          Geometry, Volume 42, Issue 4. :arxiv:`0710.4346`,
          :doi:`10.1007/s00454-008-9120-8`
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.modules.free_module import FreeModule
        n = self.size()
        vector_e = FreeModule(ZZ, n).basis()
        convert = {ind: i for i, ind in enumerate(self.groundset())}
        vertices = [sum(vector_e[convert[i]] for i in B)
                    for B in self.bases()]
        return Polyhedron(vertices)

    def independence_matroid_polytope(self):
        """
        Return the independence matroid polytope of ``self``.

        This is defined as the convex hull of the vertices

        .. MATH::

            \sum_{i \in I} e_i

        over all independent sets `I` of the matroid. Here `e_i` are
        the standard basis vectors of `\RR^n`. An arbitrary labelling
        of the groundset by `{0,\dots,n-1}` is chosen.

        .. SEEALSO::

            :meth:`matroid_polytope`

        EXAMPLES::

            sage: M = matroids.Whirl(4)
            sage: M.independence_matroid_polytope()
            A 8-dimensional polyhedron in ZZ^8 defined as the convex hull
            of 135 vertices

            sage: M = matroids.named_matroids.NonFano()
            sage: M.independence_matroid_polytope()
            A 7-dimensional polyhedron in ZZ^7 defined as the convex hull
            of 58 vertices

        REFERENCE:

        [DLHK2007]_
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.modules.free_module import FreeModule
        n = self.size()
        ambient = FreeModule(ZZ, n)
        vector_e = ambient.basis()
        convert = {ind: i for i, ind in enumerate(self.groundset())}
        vertices = [ambient.sum(vector_e[convert[i]] for i in IS)
                    for r in range(self.full_rank() + 1)
                    for IS in self.independent_r_sets(r)]
        return Polyhedron(vertices)

    # isomorphism and equality

    cpdef is_isomorphic(self, other):
        r"""
        Test matroid isomorphism.

        Two matroids `M` and `N` are *isomorphic* if there is a bijection `f`
        from the groundset of `M` to the groundset of `N` such that a subset
        `X` is independent in `M` if and only if `f(X)` is independent in `N`.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.CompleteGraphic(4)
            sage: M1.is_isomorphic(M2)
            True
            sage: G3 = graphs.CompleteGraph(4)
            sage: M1.is_isomorphic(G3)
            Traceback (most recent call last):
            ...
            TypeError: can only test for isomorphism between matroids.


            sage: M1 = matroids.named_matroids.Fano()
            sage: M2 = matroids.named_matroids.NonFano()
            sage: M1.is_isomorphic(M2)
            False
        """
        if not isinstance(other, Matroid):
            raise TypeError("can only test for isomorphism between matroids.")
        return self._is_isomorphic(other)

    cpdef _is_isomorphic(self, other):
        """
        Test if ``self`` is isomorphic to ``other``.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.CompleteGraphic(4)
            sage: M1._is_isomorphic(M2)
            True

            sage: M1 = matroids.named_matroids.Fano()
            sage: M2 = matroids.named_matroids.NonFano()
            sage: M1._is_isomorphic(M2)
            False

        """
        if self is other:
            return True
        return (self.full_rank() == other.full_rank() and self.nonbases()._isomorphism(other.nonbases()) is not None)

    cpdef isomorphism(self, other):
        r"""
        Return a matroid isomorphism.

        Two matroids `M` and `N` are *isomorphic* if there is a bijection `f`
        from the groundset of `M` to the groundset of `N` such that a subset
        `X` is independent in `M` if and only if `f(X)` is independent in `N`.
        This method returns one isomorphism `f` from self to other, if such an isomorphism exists.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        A dictionary, or ``None``.

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.CompleteGraphic(4)
            sage: morphism=M1.isomorphism(M2)
            sage: M1.is_isomorphism(M2, morphism)
            True
            sage: G3 = graphs.CompleteGraph(4)
            sage: M1.isomorphism(G3)
            Traceback (most recent call last):
            ...
            TypeError: can only give isomorphism between matroids.

            sage: M1 = matroids.named_matroids.Fano()
            sage: M2 = matroids.named_matroids.NonFano()
            sage: M1.isomorphism(M2) is not None
            False
        """
        if not isinstance(other, Matroid):
            raise TypeError("can only give isomorphism between matroids.")
        return self._isomorphism(other)

    cpdef _isomorphism(self, other):
        """
        Return isomorphism from ``self`` to ``other``, if such an isomorphism exists.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        A dictionary, or ``None``

        EXAMPLES::

            sage: M1 = matroids.Wheel(3)
            sage: M2 = matroids.CompleteGraphic(4)
            sage: morphism=M1.isomorphism(M2)
            sage: M1.is_isomorphism(M2, morphism)
            True
            sage: M1 = matroids.named_matroids.Fano()
            sage: M2 = matroids.named_matroids.NonFano()
            sage: M1.isomorphism(M2) is not None
            False
        """
        if self is other:
            return {e:e for e in self.groundset()}
        if self.full_rank() == other.full_rank():
            return self.nonbases()._isomorphism(other.nonbases())
        else:
            return None
    
    cpdef equals(self, other):
        """
        Test for matroid equality.

        Two matroids `M` and `N` are *equal* if they have the same groundset
        and a subset `X` is independent in `M` if and only if it is
        independent in `N`.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        Boolean.

        .. NOTE::

            This method tests abstract matroid equality. The ``==`` operator
            takes a more restricted view: ``M == N`` returns ``True`` only if

            #. the internal representations are of the same type,
            #. those representations are equivalent (for an appropriate
               meaning of "equivalent" in that class), and
            #. ``M.equals(N)``.

        EXAMPLES:

        A :class:`BinaryMatroid <sage.matroids.linear_matroid.BinaryMatroid>`
        and :class:`BasisMatroid <sage.matroids.basis_matroid.BasisMatroid>`
        use different representations of the matroid internally, so ``==``
        yields ``False``, even if the matroids are equal::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.named_matroids.Fano()
            sage: M
            Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
            sage: M1 = BasisMatroid(M)
            sage: M2 = Matroid(groundset='abcdefg', reduced_matrix=[
            ....:      [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1]], field=GF(2))
            sage: M.equals(M1)
            True
            sage: M.equals(M2)
            True
            sage: M == M1
            False
            sage: M == M2
            True

        :class:`LinearMatroid <sage.matroids.linear_matroid.LinearMatroid>`
        instances ``M`` and ``N`` satisfy ``M == N`` if the representations
        are equivalent up to row operations and column scaling::

            sage: M1 = LinearMatroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = LinearMatroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 3]]))
            sage: M3 = LinearMatroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[2, 6, 1, 0], [6, 1, 0, 1]]))
            sage: M1.equals(M2)
            True
            sage: M1.equals(M3)
            True
            sage: M1 == M2
            False
            sage: M1 == M3
            True
        """
        if self is other:
            return True
        if not isinstance(other, Matroid):
            raise TypeError("can only test for isomorphism between matroids.")
        if (not self.size() == other.size() or len(other.groundset().difference(self.groundset())) > 0):
            return False
        morphism = {e: e for e in self.groundset()}
        return self._is_isomorphism(other, morphism)

    cpdef is_isomorphism(self, other, morphism):
        """
        Test if a provided morphism induces a matroid isomorphism.

        A *morphism* is a map from the groundset of ``self`` to the groundset
        of ``other``.

        INPUT:

        - ``other`` -- A matroid.
        - ``morphism`` -- A map. Can be, for instance,
          a dictionary, function, or permutation.

        OUTPUT:

        Boolean.

        ..SEEALSO::

            :meth:`M.is_isomorphism() <sage.matroids.matroid.Matroid.is_isomorphism>`

        .. NOTE::

            If you know the input is valid, consider using the faster method
            ``self._is_isomorphism``.

        EXAMPLES:

        ::

            sage: M = matroids.named_matroids.Pappus()
            sage: N = matroids.named_matroids.NonPappus()
            sage: N.is_isomorphism(M, {e:e for e in M.groundset()})
            False

            sage: M = matroids.named_matroids.Fano() \ ['g']
            sage: N = matroids.Wheel(3)
            sage: morphism = {'a':0, 'b':1, 'c': 2, 'd':4, 'e':5, 'f':3}
            sage: M.is_isomorphism(N, morphism)
            True

        A morphism can be specified as a dictionary (above), a permutation,
        a function, and many other types of maps::

            sage: M = matroids.named_matroids.Fano()
            sage: P = PermutationGroup([[('a', 'b', 'c'),
            ....:                        ('d', 'e', 'f'), ('g')]]).gen()
            sage: M.is_isomorphism(M, P)
            True

            sage: M = matroids.named_matroids.Pappus()
            sage: N = matroids.named_matroids.NonPappus()
            sage: def f(x):
            ....:     return x
            ....:
            sage: N.is_isomorphism(M, f)
            False
            sage: N.is_isomorphism(N, f)
            True

        There is extensive checking for inappropriate input::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.is_isomorphism(graphs.CompleteGraph(4), lambda x:x)
            Traceback (most recent call last):
            ...
            TypeError: can only test for isomorphism between matroids.

            sage: M = matroids.CompleteGraphic(4)
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5]
            sage: M.is_isomorphism(M, {0: 1, 1: 2, 2: 3})
            Traceback (most recent call last):
            ...
            ValueError: domain of morphism does not contain groundset of this
            matroid.

            sage: M = matroids.CompleteGraphic(4)
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5]
            sage: M.is_isomorphism(M, {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1})
            Traceback (most recent call last):
            ...
            ValueError: range of morphism does not contain groundset of other
            matroid.

            sage: M = matroids.CompleteGraphic(3)
            sage: N = Matroid(bases=['ab', 'ac', 'bc'])
            sage: f = [0, 1, 2]
            sage: g = {'a': 0, 'b': 1, 'c': 2}
            sage: N.is_isomorphism(M, f)
            Traceback (most recent call last):
            ...
            ValueError: the morphism argument does not seem to be an
            isomorphism.

            sage: N.is_isomorphism(M, g)
            True
        """
        from copy import copy
        if not isinstance(other, Matroid):
            raise TypeError("can only test for isomorphism between matroids.")
        if not isinstance(morphism, dict):
            mf = {}
            try:
                for e in self.groundset():
                    mf[e] = morphism[e]
            except (IndexError, TypeError, ValueError):
                try:
                    for e in self.groundset():
                        mf[e] = morphism(e)
                except (TypeError, ValueError):
                    raise ValueError("the morphism argument does not seem to be an isomorphism.")
        else:
            mf = morphism
            if len(self.groundset().difference(mf.keys())) > 0:
                raise ValueError("domain of morphism does not contain groundset of this matroid.")
        if len(other.groundset().difference([mf[e] for e in self.groundset()])) > 0:
            raise ValueError("range of morphism does not contain groundset of other matroid.")
        if self is other:
            return self._is_isomorphism(copy(other), mf)
        return self._is_isomorphism(other, mf)

    cpdef _is_isomorphism(self, other, morphism):
        """
        Version of is_isomorphism() that does no type checking.

        INPUT:

        - ``other`` -- A matroid instance.
        - ``morphism`` -- a dictionary mapping the groundset of ``self`` to
          the groundset of ``other``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: N = matroids.named_matroids.NonPappus()
            sage: N._is_isomorphism(M, {e:e for e in M.groundset()})
            False

            sage: M = matroids.named_matroids.Fano() \ ['g']
            sage: N = matroids.Wheel(3)
            sage: morphism = {'a':0, 'b':1, 'c': 2, 'd':4, 'e':5, 'f':3}
            sage: M._is_isomorphism(N, morphism)
            True
        """
        import basis_exchange_matroid
        import basis_matroid
        sf = basis_matroid.BasisMatroid(self)
        if not isinstance(other, basis_exchange_matroid.BasisExchangeMatroid):
            ot = basis_matroid.BasisMatroid(other)
        else:
            ot = other
        return sf._is_isomorphism(ot, morphism)

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to ``__richcmp__`` (in Cython) and
            ``__cmp__`` or ``__eq__``/``__ne__`` (in Python). If you override
            one, you should (and in Cython: MUST) override the other!

        EXAMPLES::

            sage: M = matroids.named_matroids.BetsyRoss()
            sage: N = matroids.named_matroids.BetsyRoss()
            sage: hash(M) == hash(N)
            True
            sage: O = matroids.named_matroids.TicTacToe()
            sage: hash(M) == hash(O)
            False
        """
        return hash((self.groundset(), self.full_rank()))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. The default implementation below, testing matroid
        equality, should be overridden by subclasses.

        .. TODO::

            In a user guide, write about "pitfalls": testing something like
            ``M in S`` could yield ``False``, even if ``N.equals(M)`` is
            ``True`` for some ``N`` in ``S``.

        .. WARNING::

            This method is linked to ``__hash__``. If you override one, you
            MUST override the other!

        .. SEEALSO::

            :meth:`M.equals() <sage.matroids.matroid.Matroid.equals>`

        EXAMPLES::

            sage: M1 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 2]]))
            sage: M2 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[1, 0, 1, 1], [0, 1, 1, 3]]))
            sage: M3 = Matroid(groundset='abcd', matrix=Matrix(GF(7),
            ....:                               [[2, 6, 1, 0], [6, 1, 0, 1]]))
            sage: M1.equals(M2)
            True
            sage: M1.equals(M3)
            True
            sage: M1 != M2  # indirect doctest
            True
            sage: M1 == M3  # indirect doctest
            True
        """
        import basis_matroid
        if op in [0, 1, 4, 5]:  # <, <=, >, >=
            return NotImplemented
        if left.__class__ != right.__class__:
            return NotImplemented
        # The above test can be tricked. One could:
        # * spoof the __class__ attribute
        # * call the method with two things that are not Matroid instances, as in
        #   sage.matroids.matroid.Matroid.__richcmp__(p, q, 2)
        # Non-abstract subclasses should just call isinstance on both left and right.
        if hash(left) != hash(right):
            if op == 2:  # ==
                return False
            if op == 3:  # !=
                return True

        res = (basis_matroid.BasisMatroid(left) == basis_matroid.BasisMatroid(right))   # Default implementation
        if op == 2:  # ==
            return res
        if op == 3:  # !=
            return not res

    # Minors and duality

    cpdef minor(self, contractions=None, deletions=None):
        r"""
        Return the minor of ``self`` obtained by contracting, respectively
        deleting, the element(s) of ``contractions`` and ``deletions``.

        A *minor* of a matroid is a matroid obtained by repeatedly removing
        elements in one of two ways: either
        :meth:`contract <sage.matroids.matroid.Matroid.contract>` or
        :meth:`delete <sage.matroids.matroid.Matroid.delete>` them. It can be
        shown that the final matroid does not depend on the order in which
        elements are removed.

        INPUT:

        - ``contractions`` -- (default: ``None``) an element or set of
          elements to be contracted.
        - ``deletions`` -- (default: ``None``) an element or set of elements
          to be deleted.

        OUTPUT:

        A matroid.

        .. NOTE::

            The output is either of the same type as ``self``, or an instance
            of
            :class:`MinorMatroid <sage.matroids.minor_matroid.MinorMatroid>`.

        .. SEEALSO::

            :meth:`M.contract() <sage.matroids.matroid.Matroid.contract>`,
            :meth:`M.delete() <sage.matroids.matroid.Matroid.delete>`

        EXAMPLES:

        ::

            sage: M = matroids.Wheel(4)
            sage: N = M.minor(contractions=[7], deletions=[0])
            sage: N.is_isomorphic(matroids.Wheel(3))
            True

        The sets of contractions and deletions need not be independent,
        respectively coindependent::

            sage: M = matroids.named_matroids.Fano()
            sage: M.rank('abf')
            2
            sage: M.minor(contractions='abf')
            Binary matroid of rank 1 on 4 elements, type (1, 0)

        However, they need to be subsets of the groundset, and disjoint::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.minor('abc', 'defg')
            M / {'a', 'b', 'c'} \ {'d', 'e', 'f', 'g'}, where M is Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
            {'e', 'f', 'g', 'h'}, {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'}},
            4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}

            sage: M.minor('defgh', 'abc')
            M / {'d', 'e', 'f', 'g'} \ {'a', 'b', 'c', 'h'}, where M is Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
            {'e', 'f', 'g', 'h'}, {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'}},
            4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}

            sage: M.minor([1, 2, 3], 'efg')
            Traceback (most recent call last):
            ...
            ValueError: input contractions is not a subset of the groundset.
            sage: M.minor('efg', [1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: input deletions is not a subset of the groundset.
            sage: M.minor('ade', 'efg')
            Traceback (most recent call last):
            ...
            ValueError: contraction and deletion sets are not disjoint.


        .. WARNING::

            There can be ambiguity if elements of the groundset are themselves
            iterable, and their elements are in the groundset. The main
            example of this is when an element is a string. See the
            documentation of the methods
            :meth:`contract() <sage.matroids.matroid.Matroid.contract>` and
            :meth:`delete() <sage.matroids.matroid.Matroid.delete>` for an
            example of this.
        """
        try:
            if contractions in self.groundset():
                contractions = [contractions]
            # Else we expect to have to enumerate over the characters in the string.
        except TypeError:
            pass
        if (not isinstance(contractions, str) and not hasattr(contractions, '__iter__') and not contractions is None):
            contractions = [contractions]
        try:
            if deletions in self.groundset():
                deletions = [deletions]
            # Else we expect to have to enumerate over the characters in the string.
        except TypeError:
            pass
        if (not isinstance(deletions, str) and not hasattr(deletions, '__iter__') and not deletions is None):
            deletions = [deletions]
        conset, delset = sanitize_contractions_deletions(self, contractions, deletions)
        return self._minor(conset, delset)

    cpdef contract(self, X):
        r"""
        Contract elements.

        If `e` is a non-loop element, then the matroid `M / e` is a matroid
        on groundset `E(M) - e`. A set `X` is independent in `M / e` if and
        only if `X \cup e` is independent in `M`. If `e` is a loop then
        contracting `e` is the same as deleting `e`. We say that `M / e` is
        the matroid obtained from `M` by *contracting* `e`. Contracting an
        element in `M` is the same as deleting an element in the dual of `M`.

        When contracting a set, the elements of that set are contracted one by
        one. It can be shown that the resulting matroid does not depend on the
        order of the contractions.

        Sage supports the shortcut notation ``M / X`` for ``M.contract(X)``.

        INPUT:

        - ``X`` -- Either a single element of the groundset, or a collection
          of elements.

        OUTPUT:

        The matroid obtained by contracting the element(s) in ``X``.

        .. SEEALSO::

            :meth:`M.delete() <sage.matroids.matroid.Matroid.delete>`
            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`
            :meth:`M.minor() <sage.matroids.matroid.Matroid.minor>`

        EXAMPLES:

        ::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: M.contract(['a', 'c'])
            Binary matroid of rank 1 on 5 elements, type (1, 0)
            sage: M.contract(['a']) == M / ['a']
            True

        One can use a single element, rather than a set::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.contract(1) == M.contract([1])
            True
            sage: M / 1
            Regular matroid of rank 2 on 5 elements with 8 bases

        Note that one can iterate over strings::

            sage: M = matroids.named_matroids.Fano()
            sage: M / 'abc'
            Binary matroid of rank 0 on 4 elements, type (0, 0)

        The following is therefore ambiguous. Sage will contract the single
        element::

            sage: M = Matroid(groundset=['a', 'b', 'c', 'abc'],
            ....:             bases=[['a', 'b', 'c'], ['a', 'b', 'abc']])
            sage: sorted((M / 'abc').groundset())
            ['a', 'b', 'c']
        """
        return self.minor(contractions=X)

    def __div__(self, X):
        r"""
        Shorthand for ``self.contract(X)``.

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.contract(1) == M / 1  # indirect doctest
            True
        """
        # Shorthand: M / X
        return self.contract(X)

    cpdef delete(self, X):
        r"""
        Delete elements.

        If `e` is an element, then the matroid `M \setminus e` is a matroid
        on groundset `E(M) - e`. A set `X` is independent in `M \setminus e`
        if and only if `X` is independent in `M`. We say that `M \setminus e`
        is the matroid obtained from `M` by *deleting* `e`.

        When deleting a set, the elements of that set are deleted one by
        one. It can be shown that the resulting matroid does not depend on the
        order of the deletions.

        Sage supports the shortcut notation ``M \ X`` for ``M.delete(X)``.

        INPUT:

        - ``X`` -- Either a single element of the groundset, or a collection
          of elements.

        OUTPUT:

        The matroid obtained by deleting the element(s) in ``X``.

        .. SEEALSO::

            :meth:`M.contract() <sage.matroids.matroid.Matroid.contract>`
            :meth:`M.minor() <sage.matroids.matroid.Matroid.minor>`

        EXAMPLES:

        ::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: M.delete(['a', 'c'])
            Binary matroid of rank 3 on 5 elements, type (1, 6)
            sage: M.delete(['a']) == M \ ['a']
            True

        One can use a single element, rather than a set::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.delete(1) == M.delete([1])
            True
            sage: M \ 1
            Regular matroid of rank 3 on 5 elements with 8 bases

        Note that one can iterate over strings::

            sage: M = matroids.named_matroids.Fano()
            sage: M \ 'abc'
            Binary matroid of rank 3 on 4 elements, type (0, 5)

        The following is therefore ambiguous. Sage will delete the single
        element::

            sage: M = Matroid(groundset=['a', 'b', 'c', 'abc'],
            ....:             bases=[['a', 'b', 'c'], ['a', 'b', 'abc']])
            sage: sorted((M \ 'abc').groundset())
            ['a', 'b', 'c']
        """
        return self.minor(deletions=X)

    cpdef _backslash_(self, X):
        r"""
        Shorthand for ``self.delete(X)``.

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.delete(1) == M \ 1  # indirect doctest
            True
        """
        return self.delete(X)

    cpdef dual(self):
        """
        Return the dual of the matroid.

        Let `M` be a matroid with ground set `E`. If `B` is the set of bases
        of `M`, then the set `\{E - b : b \in B\}` is the set of bases of
        another matroid, the *dual* of `M`.

        .. NOTE::

            This function wraps ``self`` in a ``DualMatroid`` object. For more
            efficiency, subclasses that can, should override this method.

        OUTPUT:

        The dual matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: N = M.dual()
            sage: N.rank()
            6
            sage: N
            Dual of 'Pappus: Matroid of rank 3 on 9 elements with
            circuit-closures
            {2: {{'a', 'b', 'c'}, {'a', 'f', 'h'}, {'c', 'e', 'g'},
            {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'d', 'e', 'f'},
            {'a', 'e', 'i'}, {'b', 'd', 'i'}, {'g', 'h', 'i'}},
            3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}'
        """
        import dual_matroid
        return dual_matroid.DualMatroid(self)

    cpdef truncation(self):
        """
        Return a rank-1 truncation of the matroid.

        Let `M` be a matroid of rank `r`. The *truncation* of `M` is the
        matroid obtained by declaring all subsets of size `r` dependent. It
        can be obtained by adding an element freely to the span of the matroid
        and then contracting that element.

        OUTPUT:

        A matroid.

        .. SEEALSO::

            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`,
            :meth:`M.contract() <sage.matroids.matroid.Matroid.contract>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: N = M.truncation()
            sage: N.is_isomorphic(matroids.Uniform(2, 7))
            True
        """
        if self.full_rank() == 0:
            return None
        l = newlabel(self.groundset())
        return self._extension(l, [])._minor(contractions=frozenset([l]),
                                             deletions=frozenset([]))

    cpdef has_minor(self, N):
        """
        Check if ``self`` has a minor isomorphic to ``N``.

        INPUT:

        - ``N`` -- A matroid.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.minor() <sage.matroids.matroid.Matroid.minor>`,
            :meth:`M.is_isomorphic() <sage.matroids.matroid.Matroid.is_isomorphic>`

        .. TODO::

            This important method can (and should) be optimized considerably.
            See [Hlineny]_ p.1219 for hints to that end.

        EXAMPLES::

            sage: M = matroids.Whirl(3)
            sage: matroids.named_matroids.Fano().has_minor(M)
            False
            sage: matroids.named_matroids.NonFano().has_minor(M)
            True
        """
        if not isinstance(N, Matroid):
            raise ValueError("N must be a matroid.")
        return self._has_minor(N)

    cpdef has_line_minor(self, k, hyperlines=None):
        """
        Test if the matroid has a `U_{2, k}`-minor.

        The matroid `U_{2, k}` is a matroid on `k` elements in which every
        subset of at most 2 elements is independent, and every subset of more
        than two elements is dependent.

        The optional argument ``hyperlines`` restricts the search space: this
        method returns ``False`` if `si(M/F)` is isomorphic to `U_{2, l}` with
        `l \geq k` for some `F` in ``hyperlines``, and ``True`` otherwise.

        INPUT:

        - ``k`` -- the length of the line minor
        - ``hyperlines`` -- (default: ``None``) a set of flats of codimension
          2. Defaults to the set of all flats of codimension 2.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`Matroid.has_minor`

        EXAMPLES::

            sage: M = matroids.named_matroids.N1()
            sage: M.has_line_minor(4)
            True
            sage: M.has_line_minor(5)
            False
            sage: M.has_line_minor(k=4, hyperlines=[['a', 'b', 'c']])
            False
            sage: M.has_line_minor(k=4, hyperlines=[['a', 'b', 'c'],
            ....:                                   ['a', 'b', 'd' ]])
            True

        """
        if self.full_rank() < 2:
            return False
        if self.full_corank() < k - 2:
            return False
        if hyperlines is None:
            hyperlines = self.flats(self.full_rank() - 2)
        else:
            hyperlines = [frozenset(X) for X in hyperlines]
            allright = True
            for X in hyperlines:
                if not X.issubset(self.groundset()):
                    raise ValueError("input sets need to be subset of the groundset.")
                if not self._rank(X) == self.full_rank() - 2:
                    raise ValueError("input sets need to have rank 2 less than the rank of the matroid.")
            # Note that we don't check if the sets are flats, because loops
            # get simplified away anyway.
        return self._has_line_minor(k, hyperlines)

    cpdef _has_line_minor(self, k, hyperlines):
        """
        Test if the matroid has a `U_{2, k}`-minor.

        Internal version that does no input checking.

        INPUT:

        - ``k`` -- the length of the line minor
        - ``hyperlines`` -- (default: None) a set of flats of codimension 2.
          The flats are assumed to be ``frozenset`` compatible.

        OUTPUT:

        Boolean. ``False`` if `si(M/F)` is isomorphic to `U_{2, l}` with
        `l \geq k` for some `F` in ``hyperlines``. ``True``, otherwise.

        EXAMPLES::

            sage: M = matroids.named_matroids.NonPappus()
            sage: M._has_line_minor(5, M.flats(1))
            True
        """
        if self.full_rank() < 2:
            return False
        if self.full_corank() < k - 2:
            return False
        for F in hyperlines:
            if self._line_length(F) >= k:
                return True
        return False

    # extensions

    cpdef extension(self, element=None, subsets=None):
        r"""
        Return an extension of the matroid.

        An *extension* of `M` by an element `e` is a matroid `M'` such that
        `M' \setminus e = M`. The element ``element`` is placed such that it
        lies in the :meth:`closure <sage.matroids.matroid.Matroid.closure>` of
        each set in ``subsets``, and otherwise as freely as possible. More
        precisely, the extension is defined by the
        :meth:`modular cut <sage.matroids.matroid.Matroid.modular_cut>`
        generated by the sets in ``subsets``.

        INPUT:

        - ``element`` -- (default: ``None``) the label of the new element. If
          not specified, a new label will be generated automatically.
        - ``subsets`` -- (default: ``None``) a set of subsets of the matroid.
          The extension should be such that the new element is in the span of
          each of these. If not specified, the element is assumed to be in the
          span of the full groundset.

        OUTPUT:

        A matroid.

        .. NOTE::

            Internally, sage uses the notion of a *linear subclass* for
            matroid extension. If ``subsets`` already consists of a linear
            subclass (i.e. the set of hyperplanes of a modular cut) then the
            faster method ``M._extension()`` can be used.

        .. SEEALSO::

            :meth:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`,
            :meth:`M.modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`,
            :meth:`M.coextension() <sage.matroids.matroid.Matroid.coextension>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`

        EXAMPLES:

        First we add an element in general position::

            sage: M = matroids.Uniform(3, 6)
            sage: N = M.extension(6)
            sage: N.is_isomorphic(matroids.Uniform(3, 7))
            True

        Next we add one inside the span of a specified hyperplane::

            sage: M = matroids.Uniform(3, 6)
            sage: H = [frozenset([0, 1])]
            sage: N = M.extension(6, H)
            sage: N
            Matroid of rank 3 on 7 elements with 34 bases
            sage: [sorted(C) for C in N.circuits() if len(C) == 3]
            [[0, 1, 6]]

        Putting an element in parallel with another::

            sage: M = matroids.named_matroids.Fano()
            sage: N = M.extension('z', ['c'])
            sage: N.rank('cz')
            1
        """
        r = self.full_rank() - 1
        if element is None:
            element = newlabel(self.groundset())
        elif element in self.groundset():
            raise ValueError("cannot extend by element already in groundset")
        if subsets is None:
            hyperplanes = []
        else:
            hyperplanes = [H for H in self.modular_cut(subsets) if self._rank(H) == r]
        return self._extension(element, hyperplanes)

    cpdef coextension(self, element=None, subsets=None):
        r"""
        Return a coextension of the matroid.

        A *coextension* of `M` by an element `e` is a matroid `M'` such that
        `M' / e = M`. The element ``element`` is placed such that it
        lies in the
        :meth:`coclosure <sage.matroids.matroid.Matroid.coclosure>` of
        each set in ``subsets``, and otherwise as freely as possible.

        This is the dual method of
        :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`. See
        the documentation there for more details.

        INPUT:

        - ``element`` -- (default: ``None``) the label of the new element. If
          not specified, a new label will be generated automatically.
        - ``subsets`` -- (default: ``None``) a set of subsets of the matroid.
          The coextension should be such that the new element is in the cospan
          of each of these. If not specified, the element is assumed to be in
          the cospan of the full groundset.

        OUTPUT:

        A matroid.

        .. SEEALSO::

            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`,
            :meth:`M.coextensions() <sage.matroids.matroid.Matroid.coextensions>`,
            :meth:`M.modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`,
            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`

        EXAMPLES:

        Add an element in general position::

            sage: M = matroids.Uniform(3, 6)
            sage: N = M.coextension(6)
            sage: N.is_isomorphic(matroids.Uniform(4, 7))
            True

        Add one inside the span of a specified hyperplane::

            sage: M = matroids.Uniform(3, 6)
            sage: H = [frozenset([0, 1])]
            sage: N = M.coextension(6, H)
            sage: N
            Matroid of rank 4 on 7 elements with 34 bases
            sage: [sorted(C) for C in N.cocircuits() if len(C) == 3]
            [[0, 1, 6]]

        Put an element in series with another::

            sage: M = matroids.named_matroids.Fano()
            sage: N = M.coextension('z', ['c'])
            sage: N.corank('cz')
            1
        """
        return self.dual().extension(element, subsets).dual()

    cpdef modular_cut(self, subsets):
        """
        Compute the modular cut generated by ``subsets``.

        A *modular cut* is a collection `C` of flats such that

        - If `F \in C` and `F'` is a flat containing `F`, then `F' \in C`
        - If `F_1, F_2 \in C` form a modular pair of flats, then
          `F_1\cap F_2 \in C`.

        A *flat* is a closed set, a *modular pair* is a pair `F_1, F_2` of
        flats with `r(F_1) + r(F_2) = r(F_1\cup F_2) + r(F_1\cap F_2)`,
        where `r` is the rank function of the matroid.

        The modular cut *generated by* ``subsets`` is the smallest modular cut
        `C` for which closure`(S) \in C` for all `S` in ``subsets``.

        There is a one-to-one correspondence between the modular cuts of a
        matroid and the single-element extensions of the matroid. See [Oxley]_
        Section 7.2 for more information.

        .. NOTE::

            Sage uses linear subclasses, rather than modular cuts, internally
            for matroid extension. A linear subclass is the set of hyperplanes
            (flats of rank `r(M) - 1`) of a modular cut. It determines the
            modular cut uniquely (see [Oxley]_ Section 7.2).

        INPUT:

        - ``subsets`` -- A collection of subsets of the groundset.

        OUTPUT:

        A collection of subsets.

        .. SEEALSO::

            :meth:`M.flats() <sage.matroids.matroid.Matroid.flats>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`

        EXAMPLES:

        Any extension of the Vamos matroid where the new point is placed on
        the lines through elements `\{a, b\}` and through `\{c, d\}` is an
        extension by a loop::

            sage: M = matroids.named_matroids.Vamos()
            sage: frozenset([]) in M.modular_cut(['ab', 'cd'])
            True

        In any extension of the matroid `S_8 \setminus h`, a point on the
        lines through `\{c, g\}` and `\{a, e\}` also is on the line through
        `\{b, f\}`::

            sage: M = matroids.named_matroids.S8()
            sage: N = M \ 'h'
            sage: frozenset('bf') in N.modular_cut(['cg', 'ae'])
            True

        The modular cut of the full groundset is equal to just the groundset::

            sage: M = matroids.named_matroids.Fano()
            sage: M.modular_cut([M.groundset()]).difference(
            ....:                               [frozenset(M.groundset())])
            set()
        """
        final_list = set()
        temp_list = set([self.closure(X) for X in subsets])  # Checks validity
        while len(temp_list) > 0:
            F = temp_list.pop()
            r = self._rank(F)
            # Check modular pairs
            for FF in final_list:
                H = FF.intersection(F)
                rH = self._rank(H)
                if rH < r:
                    if rH + self._rank(FF.union(F)) == self._rank(FF) + r:
                        if not H in final_list:
                            temp_list.add(H)
            # Check upper closure (going just one level up)
            if r < self.full_rank() - 1:
                for e in self.groundset().difference(F):
                    FF = self.closure(F.union([e]))
                    if self._rank(FF) > r and not FF in final_list:
                        temp_list.add(FF)
            final_list.add(F)
        return final_list

    cpdef linear_subclasses(self, line_length=None, subsets=None):
        """
        Return an iterable set of linear subclasses of the matroid.

        A *linear subclass* is a set of hyperplanes (i.e. closed sets of rank
        `r(M) - 1`) with the following property:

        - If `H_1` and `H_2` are members, and `r(H_1 \cap H_2) = r(M) - 2`,
          then any hyperplane `H_3` containing `H_1 \cap H_2` is a member too.

        A linear subclass is the set of hyperplanes of a
        :meth:`modular cut <sage.matroids.matroid.Matroid.modular_cut>` and
        uniquely determines the modular cut. Hence the collection of linear
        subclasses is in 1-to-1 correspondence with the collection of
        single-element extensions of a matroid. See [Oxley]_, section 7.2.

        INPUT:

        - ``line_length`` -- (default: ``None``) a natural number. If given,
          restricts the output to modular cuts that generate an extension by
          `e` that does not contain a minor `N` isomorphic to `U_{2, k}`,
          where ``k > line_length``, and such that `e \in E(N)`.
        - ``subsets`` -- (default: ``None``) a collection of subsets of the
          ground set. If given, restricts the output to linear subclasses such
          that each hyperplane contains an element of ``subsets``.

        OUTPUT:

        An iterable collection of linear subclasses.

        .. NOTE::

            The ``line_length`` argument only checks for lines using the new
            element of the corresponding extension. It is still possible that
            a long line exists by contracting the new element!

        .. SEEALSO::

            :meth:`M.flats() <sage.matroids.matroid.Matroid.flats>`,
            :meth:`M.modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`,
            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: len(list(M.linear_subclasses()))
            16
            sage: len(list(M.linear_subclasses(line_length=3)))
            8
            sage: len(list(M.linear_subclasses(subsets=[{'a', 'b'}])))
            5

        The following matroid has an extension by element `e` such that
        contracting `e` creates a 6-point line, but no 6-point line minor uses
        `e`. Consequently, this method returns the modular cut, but the
        :meth:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`
        method doesn't return the corresponding extension::

            sage: M = Matroid(circuit_closures={2: ['abc', 'def'],
            ....:                               3: ['abcdef']})
            sage: len(list(M.extensions('g', line_length=5)))
            43
            sage: len(list(M.linear_subclasses(line_length=5)))
            44
        """
        import extension
        return extension.LinearSubclasses(self, line_length=line_length, subsets=subsets)

    cpdef extensions(self, element=None, line_length=None, subsets=None):
        r"""
        Return an iterable set of single-element extensions of the matroid.

        An *extension* of a matroid `M` by element `e` is a matroid `M'` such
        that `M' \setminus e = M`. By default, this method returns an iterable
        containing all extensions, but it can be restricted in two ways. If
        ``line_length`` is specified, the output is restricted to those
        matroids not containing a line minor of length `k` greater than
        ``line_length``. If ``subsets`` is specified, then the output is
        restricted to those matroids for which the new element lies in the
        :meth:`closure <sage.matroids.matroid.Matroid.closure>` of each
        member of ``subsets``.

        INPUT:

        - ``element`` -- (optional) the name of the newly added element in
          each extension.
        - ``line_length`` -- (optional) a natural number. If given, restricts
          the output to extensions that do not contain a `U_{2, k}` minor
          where ``k > line_length``.
        - ``subsets`` -- (optional) a collection of subsets of the ground set.
          If given, restricts the output to extensions where the new element
          is contained in all hyperplanes that contain an element of
          ``subsets``.

        OUTPUT:

        An iterable containing matroids.

        .. NOTE::

            The extension by a loop will always occur.
            The extension by a coloop will never occur.

        .. SEEALSO::

            :meth:`M.extension() <sage.matroids.matroid.Matroid.extension>`,
            :meth:`M.modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`,
            :meth:`M.coextensions() <sage.matroids.matroid.Matroid.coextensions>`

        EXAMPLES::

            sage: M = matroids.named_matroids.P8()
            sage: len(list(M.extensions()))
            1705
            sage: len(list(M.extensions(line_length=4)))
            41
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            sage: len(list(M.extensions(subsets=[{'a', 'b'}], line_length=4)))
            5

        """
        import extension
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")
        return extension.MatroidExtensions(self, element, line_length=line_length, subsets=subsets)  # return enumerator

    def coextensions(self, element=None, coline_length=None, subsets=None):
        r"""
        Return an iterable set of single-element coextensions of the matroid.

        A *coextension* of a matroid `M` by element `e` is a matroid `M'` such
        that `M' / e = M`. By default, this method returns an iterable
        containing all coextensions, but it can be restricted in two ways. If
        ``coline_length`` is specified, the output is restricted to those
        matroids not containing a coline minor of length `k` greater than
        ``coline_length``. If ``subsets`` is specified, then the output is
        restricted to those matroids for which the new element lies in the
        :meth:`coclosure <sage.matroids.matroid.Matroid.coclosure>` of each
        member of ``subsets``.

        This method is dual to
        :meth:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`.

        INPUT:

        - ``element`` -- (optional) the name of the newly added element in
          each coextension.
        - ``coline_length`` -- (optional) a natural number. If given,
          restricts the output to coextensions that do not contain a
          `U_{k - 2, k}` minor where ``k > coline_length``.
        - ``subsets`` -- (optional) a collection of subsets of the ground set.
          If given, restricts the output to extensions where the new element
          is contained in all cohyperplanes that contain an element of
          ``subsets``.

        OUTPUT:

        An iterable containing matroids.

        .. NOTE::

            The coextension by a coloop will always occur.
            The extension by a loop will never occur.

        .. SEEALSO::

            :meth:`M.coextension() <sage.matroids.matroid.Matroid.coextension>`,
            :meth:`M.modular_cut() <sage.matroids.matroid.Matroid.modular_cut>`,
            :meth:`M.linear_subclasses() <sage.matroids.matroid.Matroid.linear_subclasses>`,
            :mod:`sage.matroids.extension <sage.matroids.extension>`,
            :meth:`M.extensions() <sage.matroids.matroid.Matroid.extensions>`,
            :meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`

        EXAMPLES::

            sage: M = matroids.named_matroids.P8()
            sage: len(list(M.coextensions()))
            1705
            sage: len(list(M.coextensions(coline_length=4)))
            41
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            sage: len(list(M.coextensions(subsets=[{'a', 'b'}], coline_length=4)))
            5

        """
        it = self.dual().extensions(element, coline_length, subsets)
        return [N.dual() for N in it]
        # for M in it:
        #     yield M.dual() # gives segfault in 5.1

    # connectivity

    cpdef simplify(self):
        r"""
        Return the simplification of the matroid.

        A matroid is *simple* if it contains no circuits of length 1 or 2.
        The *simplification* of a matroid is obtained by deleting all loops
        (circuits of length 1) and deleting all but one element from each
        parallel class (a closed set of rank 1, that is, each pair in it forms
        a circuit of length 2).

        OUTPUT:

        A matroid.

        .. SEEALSO::

            :meth:`M.is_simple() <sage.matroids.matroid.Matroid.is_simple>`,
            :meth:`M.loops() <sage.matroids.matroid.Matroid.loops>`,
            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`M.cosimplify() <sage.matroids.matroid.Matroid.cosimplify>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano().contract('a')
            sage: M.size() - M.simplify().size()
            3

        """
        E = set(self.groundset())
        E.difference_update(self._closure(frozenset([])))  # groundset minus loops
        res = set([])

        while len(E) > 0:
            e = E.pop()
            res.add(e)
            E.difference_update(self._closure(frozenset([e])))
        return self._minor(contractions=frozenset([]),
                           deletions=self.groundset().difference(res))

    cpdef cosimplify(self):
        r"""
        Return the cosimplification of the matroid.

        A matroid is *cosimple* if it contains no cocircuits of length 1 or 2.
        The *cosimplification* of a matroid is obtained by contracting
        all coloops (cocircuits of length 1) and contracting all but one
        element from each series class (a coclosed set of rank 1, that is,
        each pair in it forms a cocircuit of length 2).

        OUTPUT:

        A matroid.

        .. SEEALSO::

            :meth:`M.is_cosimple() <sage.matroids.matroid.Matroid.is_cosimple>`,
            :meth:`M.coloops() <sage.matroids.matroid.Matroid.coloops>`,
            :meth:`M.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`,
            :meth:`M.simplify() <sage.matroids.matroid.Matroid.simplify>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano().dual().delete('a')
            sage: M.cosimplify().size()
            3

        """
        E = set(self.groundset())
        E.difference_update(self._coclosure(frozenset([])))  # groundset minus coloops
        res = set([])

        while len(E) > 0:
            e = E.pop()
            res.add(e)
            E.difference_update(self._coclosure(frozenset([e])))
        return self._minor(contractions=self.groundset().difference(res),
                           deletions=frozenset([]))

    cpdef is_simple(self):
        """
        Test if the matroid is simple.

        A matroid is *simple* if it contains no circuits of length 1 or 2.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.is_cosimple() <sage.matroids.matroid.Matroid.is_cosimple>`,
            :meth:`M.loops() <sage.matroids.matroid.Matroid.loops>`,
            :meth:`M.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`M.simplify() <sage.matroids.matroid.Matroid.simplify>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.is_simple()
            True
            sage: N = M / 'a'
            sage: N.is_simple()
            False
        """
        if len(self._closure(frozenset())) > 0:
            return False
        for e in self.groundset():
            if len(self._closure(frozenset([e]))) > 1:
                return False
        return True

    cpdef is_cosimple(self):
        """
        Test if the matroid is cosimple.

        A matroid is *cosimple* if it contains no cocircuits of length 1 or 2.

        Dual method of
        :meth:`M.is_simple() <sage.matroids.matroid.Matroid.is_simple>`.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.is_simple() <sage.matroids.matroid.Matroid.is_simple>`,
            :meth:`M.coloops() <sage.matroids.matroid.Matroid.coloops>`,
            :meth:`M.cocircuit() <sage.matroids.matroid.Matroid.cocircuit>`,
            :meth:`M.cosimplify() <sage.matroids.matroid.Matroid.cosimplify>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano().dual()
            sage: M.is_cosimple()
            True
            sage: N = M \ 'a'
            sage: N.is_cosimple()
            False
        """
        if len(self._coclosure(frozenset())) > 0:
            return False
        for e in self.groundset():
            if len(self._coclosure(frozenset([e]))) > 1:
                return False
        return True

    cpdef components(self):
        """
        Return a list of the components of the matroid.

        A *component* is an inclusionwise maximal connected subset of the
        matroid. A subset is *connected* if the matroid resulting from
        deleting the complement of that subset is
        :meth:`connected <sage.matroids.matroid.Matroid.is_connected>`.

        OUTPUT:

        A list of subsets.

        .. SEEALSO::

            :meth:`M.is_connected() <sage.matroids.matroid.Matroid.is_connected>`,
            :meth:`M.delete() <sage.matroids.matroid.Matroid.delete>`

        EXAMPLES::

            sage: from sage.matroids.advanced import setprint
            sage: M = Matroid(ring=QQ, matrix=[[1, 0, 0, 1, 1, 0],
            ....:                              [0, 1, 0, 1, 2, 0],
            ....:                              [0, 0, 1, 0, 0, 1]])
            sage: setprint(M.components())
            [{0, 1, 3, 4}, {2, 5}]
        """
        B = self.basis()
        components = [frozenset([e]) for e in self.groundset()]
        for e in self.groundset() - B:
            C = self.circuit(B | set([e]))
            components2 = []
            for comp in components:
                if len(C & comp) != 0:
                    C = C | comp
                else:
                    components2.append(comp)
            components2.append(C)
            components = components2
        return components

    cpdef is_connected(self):
        """
        Test if the matroid is connected.

        A *separation* in a matroid is a partition `(X, Y)` of the
        groundset with `X, Y` nonempty and `r(X) + r(Y) = r(X\cup Y)`.
        A matroid is *connected* if it has no separations.

        OUTPUT:

        Boolean.

        .. SEEALSO::

            :meth:`M.components() <sage.matroids.matroid.Matroid.components>`,
            :meth:`M.is_3connected() <sage.matroids.matroid.Matroid.is_3connected>`

        EXAMPLES::

            sage: M = Matroid(ring=QQ, matrix=[[1, 0, 0, 1, 1, 0],
            ....:                              [0, 1, 0, 1, 2, 0],
            ....:                              [0, 0, 1, 0, 0, 1]])
            sage: M.is_connected()
            False
            sage: matroids.named_matroids.Pappus().is_connected()
            True

        """
        return len(self.components()) == 1

    cpdef is_3connected(self):
        """
        Test if the matroid is 3-connected.

        A `k`-*separation* in a matroid is a partition `(X, Y)` of the
        groundset with `|X| \geq k, |Y| \geq k` and `r(X) + r(Y) - r(M) < k`.
        A matroid is `k`-*connected* if it has no `l`-separations for `l < k`.

        OUTPUT:

        Boolean.

        .. TODO::

            Implement this using the efficient algorithm from [BC79]_.

        .. SEEALSO::

            :meth:`M.is_connected() <sage.matroids.matroid.Matroid.is_connected>`

        EXAMPLES::

            sage: matroids.Uniform(2, 3).is_3connected()
            True
            sage: M = Matroid(ring=QQ, matrix=[[1, 0, 0, 1, 1, 0],
            ....:                              [0, 1, 0, 1, 2, 0],
            ....:                              [0, 0, 1, 0, 0, 1]])
            sage: M.is_3connected()
            False
            sage: N = Matroid(circuit_closures={2: ['abc', 'cdef'],
            ....:                               3: ['abcdef']},
            ....:             groundset='abcdef')
            sage: N.is_3connected()
            False
            sage: matroids.named_matroids.BetsyRoss().is_3connected()
            True

        ALGORITHM:

        Test all subsets `X` to see if `r(X) + r(E - X) - r(E)` does not equal
        `1`.

        """
        groundset = self.groundset()
        size = len(groundset)
        if not self.is_connected():
            return False
        if (size < 4):
            return True  # vacuously true
        r = self.full_rank()
        part_sizes = xrange(2, (size / 2) + 1)  # all possible partition sizes
        for part in part_sizes:
            subs = Subsets(groundset, part)
            for X in subs:
                Y = groundset.difference(X)
                if (self.rank(X) + self.rank(Y) - r == 1):
                    return False
        return True

    # matroid k-closed

    cpdef is_k_closed(self, k):
        r"""
        Return if ``self`` is a ``k``-closed matroid.

        We say a matroid is `k`-closed if all `k`-closed subsets
        are closed in ``M``.

        EXAMPLES::

            sage: PR = RootSystem(['A',4]).root_lattice().positive_roots()
            sage: m = matrix(map(lambda x: x.to_vector(), PR)).transpose()
            sage: M = Matroid(m)
            sage: M.is_k_closed(3)
            True
            sage: M.is_k_closed(4)
            True

            sage: PR = RootSystem(['D',4]).root_lattice().positive_roots()
            sage: m = matrix(map(lambda x: x.to_vector(), PR)).transpose()
            sage: M = Matroid(m)
            sage: M.is_k_closed(3)
            False
            sage: M.is_k_closed(4)
            True
        """
        for S in subsets(self.groundset()):
            if self.is_subset_k_closed(S, k) and not self.is_closed(S):
                return False
        return True

    # optimization

    cpdef max_weight_independent(self, X=None, weights=None):
        r"""
        Return a maximum-weight independent set contained in a subset.

        The *weight* of a subset ``S`` is ``sum(weights(e) for e in S)``.

        INPUT:

        - ``X`` -- (default: ``None``) an iterable with a subset of
          ``self.groundset()``.  If ``X`` is ``None``, the method returns a
          maximum-weight independent subset of the groundset.
        - ``weights`` -- a dictionary or function mapping the elements of
          ``X`` to nonnegative weights.

        OUTPUT:

        A subset of ``X``.

        ALGORITHM:

        The greedy algorithm. Sort the elements of ``X`` by decreasing
        weight, then greedily select elements if they are independent of
        all that was selected before.

        EXAMPLES::

            sage: from sage.matroids.advanced import setprint
            sage: M = matroids.named_matroids.Fano()
            sage: X = M.max_weight_independent()
            sage: M.is_basis(X)
            True

            sage: wt = {'a': 1, 'b': 2, 'c': 2, 'd': 1/2, 'e': 1,
            ....:       'f': 2, 'g': 2}
            sage: setprint(M.max_weight_independent(weights=wt))
            {'b', 'f', 'g'}
            sage: def wt(x):
            ....:   return x
            ....:
            sage: M = matroids.Uniform(2, 8)
            sage: setprint(M.max_weight_independent(weights=wt))
            {6, 7}
            sage: setprint(M.max_weight_independent())
            {0, 1}
            sage: M.max_weight_coindependent(X=[], weights={})
            frozenset()
        """
        if X is None:
            X = self.groundset()
        else:
            if not self.groundset().issuperset(X):
                raise ValueError("input is not a subset of the groundset.")
        if len(X) == 0:
            return frozenset([])
        if weights is None:
            Y = list(X)
        else:
            wt = []
            try:
                for e in X:
                    wt.append((weights[e], e))
            except (IndexError, TypeError, ValueError):
                try:
                    wt = []
                    for e in X:
                        wt.append((weights(e), e))
                except (TypeError, ValueError):
                    raise TypeError("the weights argument does not seem to be a collection of weights for the set X.")
            wt = sorted(wt, reverse=True)
            if wt[-1][1] < 0:
                raise ValueError("nonnegative weights were expected.")
            Y = [e for (w, e) in wt]
        res = set([])
        r = 0
        for e in Y:
            res.add(e)
            if self._rank(res) > r:
                r += 1
            else:
                res.discard(e)
        return frozenset(res)

    cpdef max_weight_coindependent(self, X=None, weights=None):
        r"""
        Return a maximum-weight coindependent set contained in ``X``.

        The *weight* of a subset ``S`` is ``sum(weights(e) for e in S)``.

        INPUT:

        - ``X`` -- (default: ``None``) an iterable with a subset of
          ``self.groundset()``.  If ``X`` is ``None``, the method returns a
          maximum-weight coindependent subset of the groundset.
        - ``weights`` -- a dictionary or function mapping the elements of
          ``X`` to nonnegative weights.

        OUTPUT:

        A subset of ``X``.

        ALGORITHM:

            The greedy algorithm. Sort the elements of ``X`` by decreasing
            weight, then greedily select elements if they are coindependent of
            all that was selected before.

        EXAMPLES::

            sage: from sage.matroids.advanced import setprint
            sage: M = matroids.named_matroids.Fano()
            sage: X = M.max_weight_coindependent()
            sage: M.is_cobasis(X)
            True

            sage: wt = {'a': 1, 'b': 2, 'c': 2, 'd': 1/2, 'e': 1, 'f': 2,
            ....:       'g': 2}
            sage: setprint(M.max_weight_coindependent(weights=wt))
            {'b', 'c', 'f', 'g'}
            sage: def wt(x):
            ....:   return x
            ....:
            sage: M = matroids.Uniform(2, 8)
            sage: setprint(M.max_weight_coindependent(weights=wt))
            {2, 3, 4, 5, 6, 7}
            sage: setprint(M.max_weight_coindependent())
            {0, 1, 2, 3, 4, 5}
            sage: M.max_weight_coindependent(X=[], weights={})
            frozenset()
        """
        if X is None:
            X = self.groundset()
        else:
            if not self.groundset().issuperset(X):
                raise ValueError("input is not a subset of the groundset.")
        if len(X) == 0:
            return frozenset([])
        if weights is None:
            Y = list(X)
        else:
            wt = []
            try:
                for e in X:
                    wt.append((weights[e], e))
            except (IndexError, TypeError, ValueError):
                try:
                    wt = []
                    for e in X:
                        wt.append((weights(e), e))
                except (TypeError, ValueError):
                    raise TypeError("the weights argument does not seem to be a collection of weights for the set X.")
            wt = sorted(wt, reverse=True)
            if wt[-1][1] < 0:
                raise ValueError("nonnegative weights were expected.")
            Y = [e for (w, e) in wt]
        res = set([])
        r = 0
        for e in Y:
            res.add(e)
            if self._corank(res) > r:
                r += 1
            else:
                res.discard(e)
        return frozenset(res)

    cpdef intersection(self, other, weights=None):
        r"""
        Return a maximum-weight common independent set.

        A *common independent set* of matroids `M` and `N` with the same
        groundset `E` is a subset of `E` that is independent both in `M` and
        `N`. The *weight* of a subset ``S`` is ``sum(weights(e) for e in S)``.

        INPUT:

        - ``other`` -- a second matroid with the same groundset as this
          matroid.
        - ``weights`` -- (default: ``None``) a dictionary which specifies a
          weight for each element of the common groundset. Defaults to the
          all-1 weight function.

        OUTPUT:

        A subset of the groundset.

        EXAMPLES::

            sage: M = matroids.named_matroids.T12()
            sage: N = matroids.named_matroids.ExtendedTernaryGolayCode()
            sage: w = {'a':30, 'b':10, 'c':11, 'd':20, 'e':70, 'f':21, 'g':90,
            ....:      'h':12, 'i':80, 'j':13, 'k':40, 'l':21}
            sage: Y = M.intersection(N, w)
            sage: sorted(Y)
            ['a', 'd', 'e', 'g', 'i', 'k']
            sage: sum([w[y] for y in Y])
            330
            sage: M = matroids.named_matroids.Fano()
            sage: N = matroids.Uniform(4, 7)
            sage: M.intersection(N)
            Traceback (most recent call last):
            ...
            ValueError: matroid intersection requires equal groundsets.
        """
        if not isinstance(other, Matroid):
            raise TypeError("can only intersect two matroids.")
        if not self.groundset() == other.groundset():
            raise ValueError("matroid intersection requires equal groundsets.")
        if weights is None:
            wt = {e: 1 for e in self.groundset()}
        else:
            wt = {}
            try:
                for e in self.groundset():
                    wt[e] = weights[e]
            except (IndexError, TypeError, ValueError):
                try:
                    wt = {}
                    for e in self.groundset():
                        wt[e] = weights(e)
                except (TypeError, ValueError):
                    raise TypeError("the weights argument does not seem to be a collection of weights for the groundset.")
        return self._intersection(other, wt)

    cpdef _intersection(self, other, weights):
        r"""
        Return a maximum-weight common independent.

        INPUT:

        - ``other`` -- a second matroid with the same groundset as this
          matroid.
        - ``weights`` -- a dictionary which must specify a weight for each
          element of the common groundset.

        OUTPUT:

        A subset of the groundset.

        .. NOTE::

            This is the unguarded version of the method
            :meth:`<sage.matroids.matroid.Matroid.intersection>`, which does
            not test if the input is well-formed.

        EXAMPLES::

            sage: M = matroids.named_matroids.T12()
            sage: N = matroids.named_matroids.ExtendedTernaryGolayCode()
            sage: w = {'a':30, 'b':10, 'c':11, 'd':20, 'e':70, 'f':21, 'g':90,
            ....:      'h':12, 'i':80, 'j':13, 'k':40, 'l':21}
            sage: Y = M._intersection(N, w)
            sage: sorted(Y)
            ['a', 'd', 'e', 'g', 'i', 'k']
            sage: sum([w[y] for y in Y])
            330
        """
        Y = set()
        U = self._intersection_augmentation(other, weights, Y)
        while U is not None and sum([weights[x] for x in U - Y]) > sum([weights[y] for y in U.intersection(Y)]):
            Y = Y.symmetric_difference(U)
            U = self._intersection_augmentation(other, weights, Y)
        return Y

    cpdef _intersection_augmentation(self, other, weights, Y):
        r"""
        Return a an augmenting set for the matroid intersection problem.

        INPUT:

        - ``other`` -- a matroid with the same ground set as ``self``.
        - ``weights`` -- a dictionary specifying a weight for each element of
          the ground set
        - ``Y`` -- an extremal common independent set of ``self`` and
          ``other`` of size `k`. That is, a common independent set of maximum
          weight among common independent sets of size `k`.

        OUTPUT:

        A set ``U`` such that the symmetric difference of ``Y`` and ``U``
        is extremal and has `k + 1` elements; or ``None``, if there is no
        common independent ste of size `k + 1`.

        .. NOTE::

            This is an unchecked method. In particular, if the given ``Y`` is
            not extremal, the algorithm will not terminate. The usage is to
            run the algorithm on its own output.

        EXAMPLES::

            sage: M = matroids.named_matroids.T12()
            sage: N = matroids.named_matroids.ExtendedTernaryGolayCode()
            sage: w = {'a':30, 'b':10, 'c':11, 'd':20, 'e':70, 'f':21, 'g':90,
            ....:      'h':12, 'i':80, 'j':13, 'k':40, 'l':21}
            sage: Y = M.intersection(N, w)
            sage: sorted(Y)
            ['a', 'd', 'e', 'g', 'i', 'k']
            sage: M._intersection_augmentation(N, w, Y) is None
            True
        """
        X = self.groundset() - Y

        X1 = frozenset([x for x in X if self._is_independent(Y.union([x]))])
        X2 = frozenset([x for x in X if other._is_independent(Y.union([x]))])

        w = {x: -weights[x] for x in X1}
        predecessor = {x: None for x in X1}
        out_neighbors = {x: set() for x in X2}

        todo = set(X1)
        next_layer = set()
        while todo:
            while todo:
                u = todo.pop()
                m = w[u]
                if u in Y:
                    if u not in out_neighbors:
                        out_neighbors[u] = X - self._closure(Y - set([u]))
                    for x in out_neighbors[u]:
                        m2 = m - weights[x]
                        if not x in w or w[x] > m2:
                            predecessor[x] = u
                            w[x] = m2
                            next_layer.add(x)
                else:
                    if u not in out_neighbors:
                        out_neighbors[u] = other._circuit(Y.union([u])) - set([u])  # have to make sure that u is not in X2 for this
                    for y in out_neighbors[u]:
                        m2 = m + weights[y]
                        if not y in w or w[y] > m2:
                            predecessor[y] = u
                            w[y] = m2
                            next_layer.add(y)
            todo = next_layer
            next_layer = set()
        X3 = X2.intersection(w)
        if not X3:
            return None
        s = min([w[x] for x in X3])
        for u in X3:
            if w[u] == s:
                path = set([u])
                while predecessor[u] is not None:
                    u = predecessor[u]
                    path.add(u)
                return path

    # invariants

    cpdef _internal(self, B):
        """
        Return the set of internally active elements of a basis `B`.

        An element `e` is *internally active* if it is the smallest element in
        the `B`-fundamental cocircuit using `e`. Smallest is interpreted as
        the output of the built-in ``min`` function on the subset.

        The `B`-fundamental cocircuit using `e` is the unique cocircuit
        intersecting basis `B` in exactly element `e`.

        INPUT:

        - ``B`` -- a basis of the matroid, assumed to have Python's
          ``frozenset`` interface.

        OUTPUT:

        A subset of ``B``.

        .. SEEALSO::

            :meth:`M.tutte_polynomial() <sage.matroids.matroid.Matroid.tutte_polynomial>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted(M._internal({'a', 'b', 'c'}))
            ['a', 'b', 'c']
            sage: sorted(M._internal({'e', 'f', 'g'}))
            []
        """
        N = self.groundset() - B
        A = set()
        for e in B:
            if min(self._cocircuit(N | set([e]))) == e:
                A.add(e)
        return A

    cpdef _external(self, B):
        """
        Return the set of externally active elements of a basis `B`.

        An element `e` is *externally active* if it is the smallest element in
        the `B`-fundamental circuit using `e`. Smallest is interpreted as
        the output of the built-in ``min`` function on the subset.

        The `B`-fundamental circuit using `e` is the unique circuit contained
        in `B + e`.

        INPUT:

        - ``B`` -- a basis of the matroid, assumed to have Python's
          ``frozenset`` interface.

        OUTPUT:

        A subset of ``self.groundset() - B``.

        .. SEEALSO::

            :meth:`M.tutte_polynomial() <sage.matroids.matroid.Matroid.tutte_polynomial>`

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: sorted(M._external({'a', 'b', 'c'}))
            []
            sage: sorted(M._external({'e', 'f', 'g'}))
            ['a', 'b', 'c', 'd']

        """
        N = self.groundset() - B
        A = set()
        for e in N:
            if min(self._circuit(B | set([e]))) == e:
                A.add(e)
        return A

    cpdef tutte_polynomial(self, x=None, y=None):
        """
        Return the Tutte polynomial of the matroid.

        The *Tutte polynomial* of a matroid is the polynomial

        .. MATH::

            T(x, y) = \sum_{A \subseteq E} (x - 1)^{r(E) - r(A)} (y - 1)^{r^*(E) - r^*(E\setminus A)},

        where `E` is the groundset of the matroid, `r` is the rank function,
        and `r^*` is the corank function. Tutte defined his polynomial
        differently:

        .. MATH::

            T(x, y)=\sum_{B} x^i(B) y^e(B),

        where the sum ranges over all bases of the matroid, `i(B)` is the
        number of internally active elements of `B`, and `e(B)` is the number
        of externally active elements of `B`.

        INPUT:

        - ``x`` -- (optional) a variable or numerical argument.
        - ``y`` -- (optional) a variable or numerical argument.

        OUTPUT:

        The Tutte-polynomial `T(x, y)`, where `x` and `y` are substituted with
        any values provided as input.

        .. TODO::

            Make implementation more efficient, e.g. generalizing the
            approach from :trac:`1314` from graphs to matroids.

        EXAMPLES::

            sage: M = matroids.named_matroids.Fano()
            sage: M.tutte_polynomial()
            y^4 + x^3 + 3*y^3 + 4*x^2 + 7*x*y + 6*y^2 + 3*x + 3*y
            sage: M.tutte_polynomial(1, 1) == M.bases_count()
            True

        ALGORITHM:

        Enumerate the bases and compute the internal and external activities
        for each `B`.
        """
        a = x
        b = y
        R = ZZ['x, y']
        x, y = R._first_ngens(2)
        T = R(0)
        for B in self.bases():
            T += x ** len(self._internal(B)) * y ** len(self._external(B))
        if a is not None and b is not None:
            T = T(a, b)
        return T

    cpdef flat_cover(self):
        """
        Return a minimum-size cover of the nonbases by non-spanning flats.

        A *nonbasis* is a subset that has the size of a basis, yet is
        dependent. A *flat* is a closed set.

        .. SEEALSO::

            :meth:`M.nonbases() <sage.matroids.matroid.Matroid.nonbases>`,
            :meth:`M.flats() <sage.matroids.matroid.Matroid.flats>`

        EXAMPLES::

            sage: from sage.matroids.advanced import setprint
            sage: M = matroids.named_matroids.Fano()
            sage: setprint(M.flat_cover())
            [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'},
             {'b', 'c', 'd'}, {'b', 'e', 'g'}, {'c', 'f', 'g'},
             {'d', 'e', 'f'}]

        """
        NB = self.nonbases()
        if len(NB) == 0:
            return []
        FF = []
        for r in range(self.full_rank()):
            FF.extend(self.flats(r))

        MIP = MixedIntegerLinearProgram(maximization=False)
        f = MIP.new_variable(binary=True)
        MIP.set_objective(sum([f[F] for F in FF]))
        for N in NB:
            MIP.add_constraint(sum([f[F] for F in FF if len(F.intersection(N)) > self.rank(F)]), min=1)
        opt = MIP.solve()

        fsol = MIP.get_values(f)
        eps = 0.00000001

        return [F for F in FF if fsol[F] > 1 - eps]
        
    cpdef plot(self, B=None, lineorders=None, pos_method=None,pos_dict=None,save_pos=False):
        """
        Return geometric representation as a sage graphics object.

        INPUT:

        - ``B`` -- (optional) a list containing a basis.
          If internal point placement is used, these elements will be placed as vertices of a triangle.
        - ``lineorders`` -- (optional) A list of lists where each of the inner lists
          specify ground set elements in a certain order which will be used to draw the
          corresponding line in geometric representation (if it exists).
        - ``pos_method`` -- An integer specifying positioning method.
            - ``0``: default positioning
            - ``1``: use pos_dict if it is not ``None``
            - ``2``: Force directed (Not yet implemented).
            
        - ``pos_dict``: A dictionary mapping ground set elements to their (x,y) positions.
        - ``save_pos``: A boolean indicating that point placements (either internal or user provided) and
          line orders (if provided) will be cached in the matroid (``M._cached_info``) and can be used for
          reproducing the geometric representation during the same session

        OUTPUT:

        A sage graphics object of type <class 'sage.plot.graphics.Graphics'> that
        corresponds to the geometric representation of the matroid

        EXAMPLES::

            sage: M=matroids.named_matroids.Fano()
            sage: G=M.plot()
            sage: type(G)
            <class 'sage.plot.graphics.Graphics'>
            sage: G.show()

        """
        import matroids_plot_helpers
        if pos_method == 1  and pos_dict != None:
        # check sanity of pos_dict and add it to cached info if sane
            if matroids_plot_helpers.posdict_is_sane(self, pos_dict) == True: 
                self._cached_info={'plot_positions':pos_dict, 'lineorders':lineorders}
        # placeholder for aditional placement methods. Only need to compute positions and update self._cached_info
        elif pos_method == 2:
            raise NotImplementedError

        if self._cached_info == None:
            self._cached_info={'plot_positions':None,'lineorders': None}
        if 'plot_positions' not in self._cached_info.keys():
            self._cached_info['plot_positions'] = None
        if 'lineorders'  not in self._cached_info.keys():
            self._cached_info['lineorders'] = None

        if self.rank() > 3:
            raise NotImplementedError
        elif B == None:
            B = list(self.basis())
        elif B != None and self.is_basis(B)==False:
            return
        lineorders2=matroids_plot_helpers.lineorders_union(self._cached_info['lineorders'],lineorders)
        return matroids_plot_helpers.geomrep(self,B,lineorders2,pd=pos_dict, sp=save_pos)

    cpdef show(self,B=None,lineorders=None,pos_method=None,pos_dict=None,save_pos=False,lims=None):
        """
        Show the geometric representation of the matroid.

        INPUT:

        - ``B`` -- (optional) a list containing elements of the groundset not in any particular order.
          If internal point placement is used, these elements will be placed as vertices of a triangle.
        - ``lineorders`` -- (optional) A list of lists where each of the inner lists
          specify ground set elements in a certain order which will be used to draw the
          corresponding line in geometric representation (if it exists).
        - ``pos_method`` -- An integer specifying positioning method
            - ``0``: default positioning
            - ``1``: use pos_dict if it is not ``None``
            - ``2``: Force directed (Not yet implemented).
            
        - ``pos_dict`` -- A dictionary mapping ground set elements to their (x,y) positions.
        - ``save_pos`` -- A boolean indicating that point placements (either internal or user provided) and
          line orders (if provided) will be cached in the matroid (``M._cached_info``) and can be used for
          reproducing the geometric representation during the same session
        - ``lims`` -- A list of 4 elements ``[xmin,xmax,ymin,ymax]``

        EXAMPLES::

            sage: M=matroids.named_matroids.TernaryDowling3()
            sage: M.show(B=['a','b','c'])
            sage: M.show(B=['a','b','c'],lineorders=[['f','e','i']])
            sage: pos = {'a':(0,0), 'b': (0,1), 'c':(1,0), 'd':(1,1), 'e':(1,-1), 'f':(-1,1), 'g':(-1,-1),'h':(2,0), 'i':(0,2)}
            sage: M.show(pos_method=1, pos_dict=pos,lims=[-3,3,-3,3])
        """
        if self.rank() > 3:
            raise NotImplementedError
        elif B == None:
            B = list(self.basis())
        elif B != None and self.is_basis(B)==False:
            return
        B1=B
        lineorders1=lineorders
        pm=pos_method
        pd=pos_dict
        sp=save_pos
        G=self.plot(B1,lineorders1,pm,pd,sp)
        if lims == None:
            G.show()
        else:
            G.show(xmin=lims[0], xmax=lims[1], ymin=lims[2], ymax=lims[3])
        return

    cpdef _fix_positions(self,pos_dict=None,lineorders=None):
        """
        Cache point positions and line orders without actually plotting

        INPUT:

        - ``pos_dict`` -- (optional) A dictionary mapping ground set elements to their (x,y) positions.
        - ``lineorders`` -- (optional) A list of lists where each of the inner lists
          specify ground set elements in a certain order which will be used to draw the
          corresponding line in geometric representation (if it exists).

        EXAMPLES::

            sage: M=matroids.named_matroids.BetsyRoss()
            sage: pos={}
            sage: s="abcde"
            sage: t="fghij"
            sage: x=1.61
            sage: y=1/1.61
            sage: for i in range(5):
            ....:         pos[s[i]]=(RR(x*sin(2*pi*i/5)), RR(x*cos(2*pi*i/5)))
            ....:         pos[t[i]]=(RR(y*sin(2*pi*(i+1/2)/5)), RR(y*cos(2*pi*(i+1/2)/5)))
            ....:
            sage: pos['k']=(0,0)
            sage: M._fix_positions(pos_dict=pos)
            sage: M._cached_info['lineorders'] is None
            True
            sage: M._cached_info['plot_positions']['k']
            (0, 0)
        """
        if self.rank() > 3:
            raise NotImplementedError
        # check sanity of pos_dict and add it to cached info if sane
        if(pos_dict!=None):
            import matroids_plot_helpers
            if matroids_plot_helpers.posdict_is_sane(self,pos_dict) ==True:
                self._cached_info={'plot_positions':pos_dict,'lineorders':lineorders}
        return

    def broken_circuit_complex(self, ordering=None):
        r"""
        Return the broken circuit complex of ``self``.

        The broken circuit complex of a matroid with a total ordering `<`
        on the ground set is obtained from the
        :meth:`NBC sets <no_broken_circuits_sets>` under subset inclusion.

        INPUT:

        - ``ordering`` -- a total ordering of the groundset given as a list

        OUTPUT:

        A simplicial complex of the NBC sets under inclusion.

        EXAMPLES::

            sage: M = Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]])
            sage: M.broken_circuit_complex()
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 4), (1, 3, 5), (1, 2, 5), (1, 2, 4)}
            sage: M.broken_circuit_complex([5,4,3,2,1])
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 4, 5), (2, 3, 5), (1, 3, 5), (2, 4, 5)}
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        return SimplicialComplex(self.no_broken_circuits_sets(ordering))

